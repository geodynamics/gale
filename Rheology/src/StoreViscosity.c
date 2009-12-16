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
** $Id: StoreViscosity.c 788 2008-08-15 04:20:57Z LukeHodkinson $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>

#include "types.h"
#include "RheologyClass.h"
#include "StoreViscosity.h"
#include "ConstitutiveMatrix.h"

#include <assert.h>

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type StoreVisc_Type = "StoreVisc";

/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
StoreVisc* _StoreVisc_New(  STOREVISC_DEFARGS  ) 
{
	StoreVisc*					self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( _sizeOfSelf >= sizeof(StoreVisc) );
	self = (StoreVisc*) _Rheology_New(  RHEOLOGY_PASSARGS  );
	
	return self;
}

void _StoreVisc_Init(
		StoreVisc*                                         self,
		MaterialPointsSwarm*                               materialPointsSwarm )
{
	StandardParticle           particle;
	StoreVisc_ParticleExt*     particleExt;
	SwarmVariable*             swarmVariable;
		
	/* Assign Pointers */
	self->materialPointsSwarm           = materialPointsSwarm;

	self->particleExtHandle = 
			ExtensionManager_Add( materialPointsSwarm->particleExtensionMgr, self->type, sizeof( StoreVisc_ParticleExt ) );
	
	/* Add SwarmVariables for plotting */
	particleExt = ExtensionManager_Get( materialPointsSwarm->particleExtensionMgr, &particle, self->particleExtHandle );
	
	self->swarmVariable = Swarm_NewScalarVariable(
		materialPointsSwarm,
		"Viscosity",
		(ArithPointer) &particleExt->effVisc - (ArithPointer) &particle,
		Variable_DataType_Double );
}

void* _StoreVisc_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                                     _sizeOfSelf = sizeof(StoreVisc);
	Type                                                             type = StoreVisc_Type;
	Stg_Class_DeleteFunction*                                     _delete = _Rheology_Delete;
	Stg_Class_PrintFunction*                                       _print = _Rheology_Print;
	Stg_Class_CopyFunction*                                         _copy = _Rheology_Copy;
	Stg_Component_DefaultConstructorFunction*         _defaultConstructor = _StoreVisc_DefaultNew;
	Stg_Component_ConstructFunction*                           _construct = _StoreVisc_AssignFromXML;
	Stg_Component_BuildFunction*                                   _build = _Rheology_Build;
	Stg_Component_InitialiseFunction*                         _initialise = _Rheology_Initialise;
	Stg_Component_ExecuteFunction*                               _execute = _Rheology_Execute;
	Stg_Component_DestroyFunction*                               _destroy = _Rheology_Destroy;
	Rheology_ModifyConstitutiveMatrixFunction*  _modifyConstitutiveMatrix = _StoreVisc_ModifyConstitutiveMatrix;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

	return (void*) _StoreVisc_New(  STOREVISC_PASSARGS  );
}

void _StoreVisc_AssignFromXML( void* rheology, Stg_ComponentFactory* cf, void* data ){
	StoreVisc*              self              = (StoreVisc*)rheology;
	MaterialPointsSwarm* materialPointsSwarm;

	/* Construct Parent */
	_Rheology_AssignFromXML( self, cf, data );

	materialPointsSwarm = Stg_ComponentFactory_ConstructByKey( 
		cf, 
		self->name, 
		"MaterialPointsSwarm", 
		MaterialPointsSwarm, 
		True,
		data );

	_StoreVisc_Init( self, materialPointsSwarm );
}
	
void _StoreVisc_Build( void* _self, void* data ) {
	StoreVisc*  self               = (StoreVisc*) _self;

	/* Build parent */
	_Rheology_Build( self, data );

	Stg_Component_Build(	self->swarmVariable, data, False);	
	Stg_Component_Build(	self->materialPointsSwarm, data, False);

}

void _StoreVisc_Initialise( void* _self, void* data ) {
	StoreVisc*  self               = (StoreVisc*) _self;

	/* Initialise parent */
	_Rheology_Initialise( self, data );

	Stg_Component_Initialise(	self->swarmVariable, data, False);	
	Stg_Component_Initialise(	self->materialPointsSwarm, data, False);

}

void _StoreVisc_Destroy( void* _self, void* data ) {
	StoreVisc*  self               = (StoreVisc*) _self;

	Stg_Component_Destroy(	self->swarmVariable, data, False);	
	Stg_Component_Destroy(	self->materialPointsSwarm, data, False);

	/* Destroy parent */
	_Rheology_Destroy( self, data );

}

void _StoreVisc_ModifyConstitutiveMatrix( 
		void*                                              rheology, 
		ConstitutiveMatrix*                                constitutiveMatrix,
		MaterialPointsSwarm*                               materialPointsSwarm,
		Element_LocalIndex                                 lElement_I,
		MaterialPoint*                                     materialPoint,
		Coord                                              xi )
{
	StoreVisc*	                      self              = (StoreVisc*) rheology;
	StoreVisc_ParticleExt*            particleExt;

#if 0
	/*
	** HAXOR: Throwing this flag in here to try and prevent REP from
	** overwriting the previously calculated viscosity value. */
	if ( !constitutiveMatrix->previousSolutionExists )
	  return;
#endif
                          
	/* Get Parameters From Material Extension */
	particleExt          = ExtensionManager_Get( materialPointsSwarm->particleExtensionMgr, materialPoint, self->particleExtHandle );

	particleExt->effVisc =  ConstitutiveMatrix_GetIsotropicViscosity(constitutiveMatrix);
	
}



