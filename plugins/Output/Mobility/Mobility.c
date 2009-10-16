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
#include <assert.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include <Underworld/Underworld.h>
#include <StgFEM/FrequentOutput/FrequentOutput.h>

#include "Mobility.h"

const Type Underworld_Mobility_Type = "Underworld_Mobility";

void _Underworld_Mobility_Construct( void* component, Stg_ComponentFactory* cf, void* data ) {
	UnderworldContext*  context;

	context = Stg_ComponentFactory_ConstructByName( cf, "context", UnderworldContext, True, data ); 

	Underworld_Mobility_PrintHeaderToFile( context );
	ContextEP_Append( context, AbstractContext_EP_ConstructExtensions, Underworld_Mobility_Setup );
	ContextEP_Append( context, AbstractContext_EP_FrequentOutput     , Underworld_Mobility_Dump );
}

void _Underworld_Mobility_Build( void* component, void* data ) {
	Underworld_Mobility*	self = (Underworld_Mobility*)component;

	assert( self );

	Stg_Component_Build( self->velocitySquaredField, data, False );
}

void _Underworld_Mobility_Initialise( void* component, void* data ) {
	Underworld_Mobility*	self = (Underworld_Mobility*)component;

	assert( self );

	Stg_Component_Initialise( self->velocitySquaredField, data, False );
}

void* _Underworld_Mobility_DefaultNew( Name name ) {
	return _Codelet_New(
		sizeof(Underworld_Mobility),
		Underworld_Mobility_Type,
		_Codelet_Delete,
		_Codelet_Print,
		_Codelet_Copy,
		_Underworld_Mobility_DefaultNew,
		_Underworld_Mobility_Construct,
		_Underworld_Mobility_Build,
		_Underworld_Mobility_Initialise,
		_Codelet_Execute,
		_Codelet_Destroy,
		name );
}

Index Underworld_Mobility_Register( PluginsManager* pluginsManager ) {
	return PluginsManager_Submit( pluginsManager, Underworld_Mobility_Type, "0", _Underworld_Mobility_DefaultNew );
}

void Underworld_Mobility_Setup( void* _context ) {
	UnderworldContext*                context       = (UnderworldContext*) _context;
	Swarm* gaussSwarm = (Swarm*)LiveComponentRegister_Get( context->CF->LCRegister, "gaussSwarm" );
	FeVariable* velocityField = (FeVariable*)LiveComponentRegister_Get( context->CF->LCRegister, "VelocityField" );

	Underworld_Mobility* self;

	self = (Underworld_Mobility*)LiveComponentRegister_Get(
					context->CF->LCRegister,
					Underworld_Mobility_Type );

	Journal_Firewall( 
			gaussSwarm != NULL, 
			Underworld_Error,
			"Cannot find gauss swarm. Cannot use %s.\n", CURR_MODULE_NAME );
	Journal_Firewall( 
			velocityField != NULL, 
			Underworld_Error,
			"Cannot find velocityField. Cannot use %s.\n", CURR_MODULE_NAME );

	/* Create new Field Variable */
	self->velocitySquaredField = OperatorFeVariable_NewUnary( "VelocitySquaredField", velocityField, "VectorSquare" );
}

/* Integrate Every Step and dump to file */
void Underworld_Mobility_Dump( void* _context ) {
	UnderworldContext*                   context       = (UnderworldContext*) _context;
	Swarm* gaussSwarm = (Swarm*)LiveComponentRegister_Get( context->CF->LCRegister, "gaussSwarm" );
	Mesh*			   	     mesh;
	double		    	  	     maxCrd[3], minCrd[3];
	double                               integral;
	double                               topLayerAverage;
	double                               mobility;
	double                               vrms;
	double                               volume        = 0.0;
	Dimension_Index                      dim           = context->dim;

	Underworld_Mobility* self;

	self = (Underworld_Mobility*)LiveComponentRegister_Get(
					context->CF->LCRegister,
					Underworld_Mobility_Type );

	mesh = (Mesh*)self->velocitySquaredField->feMesh;
	Mesh_GetGlobalCoordRange( mesh, minCrd, maxCrd );
	
	/* Sum integral */
	integral        = FeVariable_Integrate( self->velocitySquaredField, gaussSwarm );
	topLayerAverage = FeVariable_AverageTopLayer( self->velocitySquaredField, gaussSwarm, J_AXIS );
	/* Get Volume of Mesh - TODO Make general for irregular meshes */
	volume = ( maxCrd[ I_AXIS ] - minCrd[ I_AXIS ] ) * 
		( maxCrd[ J_AXIS ] - minCrd[ J_AXIS ] );
	if ( dim == 3 ) 
		volume *= maxCrd[ K_AXIS ] - minCrd[ K_AXIS ];

	/* Calculate Vrms 
	 * V_{rms} = \sqrt{ \frac{ \int_\Omega \mathbf{u . u} d\Omega }{\Omega} } */
	vrms = sqrt( integral / volume );

	/* Calculate Mobility */	
	mobility = sqrt(topLayerAverage)/vrms;

	/* Print data to file */
	StgFEM_FrequentOutput_PrintValue( context, mobility );

}

void Underworld_Mobility_PrintHeaderToFile( void* context ) {
	StgFEM_FrequentOutput_PrintString( context, "Mobility" );
}

