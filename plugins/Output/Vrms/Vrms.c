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
** $Id: Vrms.c 610 2007-10-11 08:09:29Z SteveQuenette $
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

#include "Vrms.h"

const Type Underworld_Vrms_Type = "Underworld_Vrms";

void _Underworld_Vrms_AssignFromXML( void* component, Stg_ComponentFactory* cf, void* data ) {
	Underworld_Vrms*	self 		= (Underworld_Vrms*)component;

	self->context       = (AbstractContext*)Stg_ComponentFactory_PluginConstructByKey( cf, self, "Context", UnderworldContext, True, data );
	self->gaussSwarm    = Stg_ComponentFactory_PluginConstructByKey( cf, self, "GaussSwarm", Swarm, True, data );
	self->velocityField = Stg_ComponentFactory_PluginConstructByKey( cf, self, "VelocityField", FeVariable, True, data );

	/* Create new Field Variable */
	self->velocitySquaredField = OperatorFeVariable_NewUnary( "VelocitySquaredField", (DomainContext*)self->context, self->velocityField, "VectorSquare" );
	self->velocitySquaredField->context = (DomainContext*)self->context;

	Underworld_Vrms_PrintHeaderToFile( self->context );
	ContextEP_Append( self->context, AbstractContext_EP_FrequentOutput, Underworld_Vrms_Dump );
}

void _Underworld_Vrms_Build( void* component, void* data ) {
	Underworld_Vrms*	self = (Underworld_Vrms*)component;

	assert( self );

   Stg_Component_Build( self->gaussSwarm, data, False );
   Stg_Component_Build( self->velocityField, data, False );
	Stg_Component_Build( self->velocitySquaredField, data, False );
   
   _Codelet_Build( self, data );
}

void _Underworld_Vrms_Initialise( void* component, void* data ) {
	Underworld_Vrms*	self = (Underworld_Vrms*)component;

	assert( self );

   Stg_Component_Initialise( self->gaussSwarm, data, False );
   Stg_Component_Initialise( self->velocityField, data, False );
	Stg_Component_Initialise( self->velocitySquaredField, data, False );
   
   _Codelet_Initialise( self, data );

}

void _Underworld_Vrms_Destroy( void* component, void* data ) {
	Underworld_Vrms*	self = (Underworld_Vrms*)component;

	assert( self );

   _Codelet_Destroy( self, data );
   
   Stg_Component_Destroy( self->gaussSwarm, data, False );
   Stg_Component_Destroy( self->velocityField, data, False );
	Stg_Component_Destroy( self->velocitySquaredField, data, False );

}

void* _Underworld_Vrms_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                              _sizeOfSelf = sizeof(Underworld_Vrms);
	Type                                                      type = Underworld_Vrms_Type;
	Stg_Class_DeleteFunction*                              _delete = _Codelet_Delete;
	Stg_Class_PrintFunction*                                _print = _Codelet_Print;
	Stg_Class_CopyFunction*                                  _copy = _Codelet_Copy;
	Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _Underworld_Vrms_DefaultNew;
	Stg_Component_ConstructFunction*                    _construct = _Underworld_Vrms_AssignFromXML;
	Stg_Component_BuildFunction*                            _build = _Underworld_Vrms_Build;
	Stg_Component_InitialiseFunction*                  _initialise = _Underworld_Vrms_Initialise;
	Stg_Component_ExecuteFunction*                        _execute = _Codelet_Execute;
	Stg_Component_DestroyFunction*                        _destroy = _Underworld_Vrms_Destroy;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = ZERO;

	return _Codelet_New(  CODELET_PASSARGS  );
}

Index Underworld_Vrms_Register( PluginsManager* pluginsManager ) {
	return PluginsManager_Submit( pluginsManager, Underworld_Vrms_Type, "0", _Underworld_Vrms_DefaultNew );
}

/* Integrate Every Step and dump to file */
void Underworld_Vrms_Dump( void* _context ) {
	UnderworldContext*                   context       = (UnderworldContext*) _context;
	Mesh*			   	     mesh;
	double		    	  	     maxCrd[3], minCrd[3];
	double                               integral;
	double                               vrms;
	double                               volume        = 0.0;
	Dimension_Index                      dim           = context->dim;

	Underworld_Vrms* self;

	self = (Underworld_Vrms*)LiveComponentRegister_Get( context->CF->LCRegister, Underworld_Vrms_Type );

	mesh = (Mesh*)self->velocitySquaredField->feMesh;
	Mesh_GetGlobalCoordRange( mesh, minCrd, maxCrd );
	
	/* Sum integral */
	integral = FeVariable_Integrate( self->velocitySquaredField, self->gaussSwarm );

	/* Get Volume of Mesh - TODO Make general for irregular meshes */
	volume = ( maxCrd[ I_AXIS ] - minCrd[ I_AXIS ] ) * 
		( maxCrd[ J_AXIS ] - minCrd[ J_AXIS ] );
	if ( dim == 3 ) 
		volume *= maxCrd[ K_AXIS ] - minCrd[ K_AXIS ];

	/* Calculate Vrms 
	 * V_{rms} = \sqrt{ \frac{ \int_\Omega \mathbf{u . u} d\Omega }{\Omega} } */
	vrms = sqrt( integral / volume );

	/* Print data to file */
	StgFEM_FrequentOutput_PrintValue( context, vrms );

	/* Put Value onto context */
	self->vrms = vrms;
}

void Underworld_Vrms_PrintHeaderToFile( void* context ) {
	StgFEM_FrequentOutput_PrintString( context, "Vrms" );
}



