/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
** Copyright (c) 2005, Monash Cluster Computing 
** All rights reserved.
** Redistribution and use in source and binary forms, with or without modification,
** are permitted provided that the following conditions are met:
**
** 		* Redistributions of source code must retain the above copyright notice, 
** 			this list of conditions and the following disclaimer.
** 		* Redistributions in binary form must reproduce the above copyright 
**			notice, this list of conditions and the following disclaimer in the 
**			documentation and/or other materials provided with the distribution.
** 		* Neither the name of the Monash University nor the names of its contributors 
**			may be used to endorse or promote products derived from this software 
**			without specific prior written permission.
**
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
** THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR 
** PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS 
** BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
** CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
** SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) 
** HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT 
** LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT 
** OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
**
**
** Contact:
*%		Louis Moresi - Louis.Moresi@sci.monash.edu.au
*%
** Contributors:
*+		Robert Turnbull
*+		Vincent Lemiale
*+		Louis Moresi
*+		David May
*+		David Stegman
*+		Mirko Velic
*+		Patrick Sunter
*+		Julian Giordani
*+
** $Id: ConvectionData.c 182 2006-05-01 12:32:01Z RobertTurnbull $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include <Underworld/Underworld.h>
#include <Underworld/BoundaryLayers/BoundaryLayers.h>
#include <StgFEM/FrequentOutput/FrequentOutput.h>
#include <assert.h>
#include <string.h>

typedef struct {
	__Codelet
	OperatorFeVariable*        velocitySquaredField;
	Underworld_BoundaryLayers* boundaryLayersPlugin;
	char*                      rheologyName;
	double                     stressExponent;
	double                     eta0;
	double                     Ra;
	double                     vrms;
	double                     diffusivity;
        double                     horizontalDimOfConvecCell;
} Underworld_ConvectionData;

double Underworld_ConvectionData_XZPlaneVrms( UnderworldContext* context, double yCoord_Of_XZPlane ); 

const Type Underworld_ConvectionData_Type = "Underworld_ConvectionData";
Stream* dataStream;

void Underworld_ConvectionData_Setup( void* _context ) {
	UnderworldContext*          context       = (UnderworldContext*) _context;
	Underworld_ConvectionData*  self;
	Arrhenius*                  arrhenius;
	FrankKamenetskii*           frankKamenetskii;
	Rheology*                   rheology;
	NonNewtonian*               nonNewtonian;
	char*  filename;

	dataStream = Journal_Register( Info_Type, "ConvectionData Info Stream" );
	Stream_SetPrintingRank( dataStream, 0 ); /** Only prints to main proccessor */
	
	self = (Underworld_ConvectionData*)LiveComponentRegister_Get(
					context->CF->LCRegister,
					Underworld_ConvectionData_Type );
	
	rheology = (Rheology*)Stg_ComponentFactory_ConstructByName( context->CF, "temperatureDependence", Rheology, False, 0 /* dummy */ );	
	nonNewtonian = (NonNewtonian*)LiveComponentRegister_Get( context->CF->LCRegister, NonNewtonian_Type );
		//(NonNewtonian*)Stg_ComponentFactory_ConstructByName( context->CF, "nonNewtonian", NonNewtonian, True, 0 /* dummy */ );	
	(nonNewtonian == NULL) ?  (self->stressExponent = 1) :
		       	(self->stressExponent = nonNewtonian->stressExponent) ;

	if( !strcmp( rheology->type, "FrankKamenetskii" ) ) {
		frankKamenetskii = (FrankKamenetskii*)rheology; 
		arrhenius=NULL;
	} else {
		arrhenius = (Arrhenius*)rheology;
		frankKamenetskii=NULL;
	}
	
	self->diffusivity = Stg_ComponentFactory_GetRootDictDouble( context->CF, "diffusivity", 1.0 );
	self->horizontalDimOfConvecCell = Stg_ComponentFactory_GetRootDictDouble( context->CF, "horizontalDimOfConvecCell", 1.0 );
	self->Ra = Stg_ComponentFactory_GetRootDictDouble( context->CF, "Ra", 1.0 );
	
	self->eta0 = ( arrhenius != NULL ? arrhenius->eta0 : frankKamenetskii->eta0 );
//	self->stressExponent = nonNewtonian->stressExponent;
	self->rheologyName = rheology->type;//( arrhenius != NULL ? StG_Strdup(arrhenius->type) : StG_Strdup(frankKamenetskii->type) );
	self->boundaryLayersPlugin = (Underworld_BoundaryLayers*)LiveComponentRegister_Get(
					context->CF->LCRegister,
					Underworld_BoundaryLayers_Type );
	Journal_Firewall( self->boundaryLayersPlugin != NULL, Underworld_Error, "Error in %s. Cannot find the BoundaryLayers Plugin. Make sure <param>Underworld_BoundaryLayers</param> is in your plugins list\n");

	Journal_Firewall( 
			context->gaussSwarm != NULL, 
			Underworld_Error,
			"Cannot find gauss swarm. Cannot use %s.\n", CURR_MODULE_NAME );
	Journal_Firewall( 
			context->velocityField != NULL, 
			Underworld_Error,
			"Cannot find velocityField. Cannot use %s.\n", CURR_MODULE_NAME );

	/* Create new Field Variable */
	self->velocitySquaredField = OperatorFeVariable_NewUnary( 
			"VelocitySquaredField", 
			context->velocityField, 
			"VectorSquare" );

	Stg_asprintf( &filename, "ConvectionData.%dof%d.dat", context->rank, context->nproc );
	Stream_RedirectFile_WithPrependedPath( dataStream, context->outputPath, filename );
	Stream_SetAutoFlush( dataStream, True );
	// Print Header
	Journal_Printf( dataStream, "#Rheology | etaContrast | Ra | UpperVrms | LowerVrms | UTBL | LTBL | SurfaceMobility\n");
	Memory_Free( filename );

}

void Underworld_ConvectionData_Dump( void* _context ) {
	UnderworldContext*         context = (UnderworldContext*) _context;
	Underworld_ConvectionData* self;
	BlockGeometry*             geometry;      
	double                     topVrms;
	double                     bottomVrms;
	double                     Upper_tbl_Thinckness;
	double                     Lower_tbl_Thinckness;
	double                     surfaceMobility;
	double                     deltaViscosity;

	self = (Underworld_ConvectionData*)LiveComponentRegister_Get(
					context->CF->LCRegister,
					Underworld_ConvectionData_Type );

	geometry = Stg_CheckType( self->velocitySquaredField->feMesh->layout->elementLayout->geometry, BlockGeometry );

	assert( self->stressExponent != 0 );
	deltaViscosity = pow(self->eta0, (-1/self->stressExponent) );
	Journal_Printf( dataStream, "%s %g %g ", self->rheologyName, deltaViscosity, self->Ra ); 
	
	// Prints out Surface Vrms 
	topVrms = Underworld_ConvectionData_XZPlaneVrms( context, geometry->max[ J_AXIS ] );
	Journal_Printf( dataStream, " %g", topVrms );

	bottomVrms = Underworld_ConvectionData_XZPlaneVrms( context, geometry->min[ J_AXIS ] );
	Journal_Printf( dataStream, " %g", bottomVrms );
	
	Upper_tbl_Thinckness = self->boundaryLayersPlugin->hotLayerThickness;
	Journal_Printf( dataStream, " %g", Upper_tbl_Thinckness );

	Lower_tbl_Thinckness = self->boundaryLayersPlugin->coldLayerThickness;
	Journal_Printf( dataStream, " %g", Lower_tbl_Thinckness );

	surfaceMobility = ( Upper_tbl_Thinckness * topVrms ) / ( self->diffusivity * self->horizontalDimOfConvecCell );
	Journal_Printf( dataStream, " %g", surfaceMobility );
/*
*/
	Journal_Printf( dataStream, "\n" );
}
	

double Underworld_ConvectionData_XZPlaneVrms( UnderworldContext* context, double yCoord_Of_XZPlane ) {
	BlockGeometry*                       geometry;
	double                               integral;
	double                               samplingSpace = 0.0;
	Dimension_Index                      dim           = context->dim;

	Underworld_ConvectionData* self;

	/*
	TODO:PAT help:
	   Here I would like to get the VelocitySquareField instead of getting the whole plugin
	 */
	self = (Underworld_ConvectionData*)LiveComponentRegister_Get(
					context->CF->LCRegister,
					Underworld_ConvectionData_Type );
	
	geometry = Stg_CheckType( self->velocitySquaredField->feMesh->layout->elementLayout->geometry, BlockGeometry );

	/* Sum integral */
	integral = FeVariable_IntegratePlane( self->velocitySquaredField, J_AXIS, yCoord_Of_XZPlane );

	/* Get Volume of Mesh - TODO Make general for irregular meshes */
	samplingSpace = ( geometry->max[ I_AXIS ] - geometry->min[ I_AXIS ] ); 
		
	if ( dim == 3 ) 
		samplingSpace *= geometry->max[ K_AXIS ] - geometry->min[ K_AXIS ];

	/* Calculate ConvectionData 
	 * V_{rms} = \sqrt{ \frac{ \int_\Omega \mathbf{u . u} d\Omega }{\Omega} } */

	return ( sqrt( integral / samplingSpace ) ) ;

}
	
void _Underworld_ConvectionData_Construct( void* component, Stg_ComponentFactory* cf, void* data ) {
	UnderworldContext*  context;

	context = Stg_ComponentFactory_ConstructByName( cf, "context", UnderworldContext, True, data );

	ContextEP_Append( context, AbstractContext_EP_ConstructExtensions, Underworld_ConvectionData_Setup );
	ContextEP_Append( context, AbstractContext_EP_FrequentOutput, Underworld_ConvectionData_Dump );
}

void _Underworld_ConvectionData_Delete( void* component ) {
	Underworld_ConvectionData* self = (Underworld_ConvectionData*) component;

	if( self->rheologyName )
		Memory_Free( self->rheologyName );
	_Codelet_Delete( self );
}
void* _Underworld_ConvectionData_DefaultNew( Name name ) {
	return _Codelet_New(
		sizeof(Underworld_ConvectionData),
		Underworld_ConvectionData_Type,
		_Underworld_ConvectionData_Delete,
		//_Codelet_Delete,
		_Codelet_Print,
		_Codelet_Copy,
		_Underworld_ConvectionData_DefaultNew,
		_Underworld_ConvectionData_Construct,
		_Codelet_Build,
		_Codelet_Initialise,
		_Codelet_Execute,
		_Codelet_Destroy,
		name );
}

Index Underworld_ConvectionData_Register( PluginsManager* pluginsManager ) {
	return PluginsManager_Submit( pluginsManager, Underworld_ConvectionData_Type, "0", _Underworld_ConvectionData_DefaultNew );
}
