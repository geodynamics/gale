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
** $Id: StoreViscosity.c 747 2008-07-04 01:36:54Z LukeHodkinson $
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
StoreVisc* _StoreVisc_New( 
		SizeT                                              sizeOfSelf,
		Type                                               type,
		Stg_Class_DeleteFunction*                          _delete,
		Stg_Class_PrintFunction*                           _print,
		Stg_Class_CopyFunction*                            _copy, 
		Stg_Component_DefaultConstructorFunction*          _defaultConstructor,
		Stg_Component_ConstructFunction*                   _construct,
		Stg_Component_BuildFunction*                       _build,
		Stg_Component_InitialiseFunction*                  _initialise,
		Stg_Component_ExecuteFunction*                     _execute,
		Stg_Component_DestroyFunction*                     _destroy,
		Rheology_ModifyConstitutiveMatrixFunction*         _modifyConstitutiveMatrix,
		Name                                               name ) 
{
	StoreVisc*					self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( sizeOfSelf >= sizeof(StoreVisc) );
	self = (StoreVisc*) _Rheology_New( 
			sizeOfSelf,
			type, 
			_delete,
			_print,
			_copy,
			_defaultConstructor,
			_construct,
			_build,
			_initialise,
			_execute,
			_destroy,
			_modifyConstitutiveMatrix,
			name );
	
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
	
	swarmVariable = Swarm_NewScalarVariable(
		materialPointsSwarm,
		"Viscosity",
		(ArithPointer) &particleExt->effVisc - (ArithPointer) &particle,
		Variable_DataType_Double );
}

void* _StoreVisc_DefaultNew( Name name ) {
	return (void*) _StoreVisc_New(
			sizeof(StoreVisc),
		 StoreVisc_Type,
		_Rheology_Delete,
		_Rheology_Print,
		_Rheology_Copy,
		_StoreVisc_DefaultNew,
		_StoreVisc_Construct,
		_Rheology_Build,
		_Rheology_Initialise,
		_Rheology_Execute,
		_Rheology_Destroy,
		_StoreVisc_ModifyConstitutiveMatrix,
		name );
}

void _StoreVisc_Construct( void* rheology, Stg_ComponentFactory* cf, void* data ){
	StoreVisc*              self              = (StoreVisc*)rheology;
	MaterialPointsSwarm* materialPointsSwarm;

	/* Construct Parent */
	_Rheology_Construct( self, cf, data );

	materialPointsSwarm = Stg_ComponentFactory_ConstructByKey( 
		cf, 
		self->name, 
		"MaterialPointsSwarm", 
		MaterialPointsSwarm, 
		True,
		data );

	_StoreVisc_Init( self, materialPointsSwarm );
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

	/*
	** HAXOR: Throwing this flag in here to try and prevent REP from
	** overwriting the previously calculated viscosity value. */
	if ( !constitutiveMatrix->previousSolutionExists )
	  return;
                          
	/* Get Parameters From Material Extension */
	particleExt          = ExtensionManager_Get( materialPointsSwarm->particleExtensionMgr, materialPoint, self->particleExtHandle );

	particleExt->effVisc =  ConstitutiveMatrix_GetIsotropicViscosity(constitutiveMatrix);
	
}

