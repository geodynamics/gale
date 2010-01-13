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
** $Id:  $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include <Underworld/Underworld.h>
#include <assert.h>


const Type Underworld_DensityChange_Type = "Underworld_DensityChange";

typedef struct { 
	__Codelet
	BuoyancyForceTerm_MaterialExt* bftExt;
	IntegrationPointsSwarm*     swarm;
	Material*  material;
	double     height;
	double     newDensity;
} Underworld_DensityChange;

void Underworld_DensityChange_Check( UnderworldContext* context ) {
	/* Function runs each timestep to check if centroid is at certain height */

	double volume;
	double centroid[3];
	static int densityChangeIsDone = 0;

	/* test if this has already happened, if so leave function */
	if ( densityChangeIsDone )
		return;

	/* Get self (the plugin) */
	Underworld_DensityChange* self = (Underworld_DensityChange*)LiveComponentRegister_Get( context->CF->LCRegister, (Name)Underworld_DensityChange_Type  );

	/* get centroid coordinate */
	volume = Material_Volume( self->material, (IntegrationPointsSwarm*)self->swarm, centroid );

	/* test if centroid height (y-coord) is >= input height
	 * 	Note following code is only run once */
	if( centroid[1] >= self->height ) {
		self->bftExt->density = self->newDensity;
		densityChangeIsDone = 1;
		Journal_RPrintf( Journal_Register( Info_Type, (Name)"DensityChange"  ), 
				"************** Density Change for material %s, new density is %g**************\n",
				self->material->name, self->bftExt->density );
	}

}
void Underworld_DensityChange_Setup( UnderworldContext* context ) {
		/* Function pulls and checks user input from the xml file */
	BuoyancyForceTerm*  bft = NULL;
	BuoyancyForceTerm_MaterialExt* materialExt = NULL;;
	Stream* stream = Journal_Register( Info_Type, (Name)"cows" );
	Name   materialName = NULL;
	int materialIndex;

	/* Get self (the plugin) */
	Underworld_DensityChange* self = (Underworld_DensityChange* )LiveComponentRegister_Get( context->CF->LCRegister, (Name)Underworld_DensityChange_Type  ); 

	/* Initialise plugin data */
	self->bftExt = NULL;
	self->swarm = NULL;
	self->material = NULL;
	self->height = 0;
	self->newDensity = 0;

	/* Need buoyancy force term, to get density infomation from bft extension + the integration swarm */
	bft = Stg_ComponentFactory_ConstructByName( context->CF, (Name)"buoyancyForceTerm", BuoyancyForceTerm, True, NULL  );

	/* Read in input from root xml dictionary */
	self->height = Dictionary_GetDouble_WithDefault( context->dictionary, (Dictionary_Entry_Key)"materialDensityChangeHeight", 0.5  );
	self->newDensity = Dictionary_GetDouble_WithDefault( context->dictionary, (Dictionary_Entry_Key)"materialDensityNewDensity", 0.3 );
	self->swarm = (IntegrationPointsSwarm* )bft->integrationSwarm;
	materialName = Dictionary_GetString( context->dictionary, (Dictionary_Entry_Key)"materialDensityToChange"  );
	self->material = Materials_Register_GetByName( self->swarm->materials_Register, materialName );

	/* check if material index exists */
	if( self->material==NULL ) {
		printf("Error\nCounld find the material with index %d\n", materialIndex ); exit(0);
	}
	materialExt = ExtensionManager_Get( self->material->extensionMgr, self->material, bft->materialExtHandle );
	Journal_RPrintf( stream, "Will change %s's density at height %g from %g to %g\n", 
			self->material->name, self->height, materialExt->density, self->newDensity );

	self->bftExt = materialExt;
}

void _Underworld_DensityChange_AssignFromXML( void* component, Stg_ComponentFactory* cf, void* data ) {
	UnderworldContext* context;

	context = (UnderworldContext*)Stg_ComponentFactory_ConstructByName( cf, (Name)"context", UnderworldContext, True, data  ); 

	/* Add functions to entry points */
	ContextEP_Append( context, AbstractContext_EP_Initialise, Underworld_DensityChange_Setup );
	ContextEP_Append( context, AbstractContext_EP_FrequentOutput, Underworld_DensityChange_Check );
	
}

void* _Underworld_DensityChange_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                              _sizeOfSelf = sizeof( Underworld_DensityChange );
	Type                                                      type = Underworld_DensityChange_Type;
	Stg_Class_DeleteFunction*                              _delete = _Codelet_Delete;
	Stg_Class_PrintFunction*                                _print = _Codelet_Print;
	Stg_Class_CopyFunction*                                  _copy = _Codelet_Copy;
	Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _Underworld_DensityChange_DefaultNew;
	Stg_Component_ConstructFunction*                    _construct = _Underworld_DensityChange_AssignFromXML;
	Stg_Component_BuildFunction*                            _build = _Codelet_Build;
	Stg_Component_InitialiseFunction*                  _initialise = _Codelet_Initialise;
	Stg_Component_ExecuteFunction*                        _execute = _Codelet_Execute;
	Stg_Component_DestroyFunction*                        _destroy = _Codelet_Destroy;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

	return _Codelet_New(  CODELET_PASSARGS  );
}

Index Underworld_DensityChangeAtDepth_Register( PluginsManager* pluginsManager ) {
	return PluginsManager_Submit( pluginsManager, Underworld_DensityChange_Type, (Name)"0", _Underworld_DensityChange_DefaultNew  );
}


