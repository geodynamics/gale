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
** $Id: StoreStress.c 747 2008-07-04 01:36:54Z LukeHodkinson $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>

#include "types.h"
#include "RheologyClass.h"
#include "StoreStress.h"
#include "ConstitutiveMatrix.h"

#include <assert.h>

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type StoreStress_Type = "StoreStress";

/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
StoreStress* _StoreStress_New(  STORESTRESS_DEFARGS  ) 
{
	StoreStress*					self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( _sizeOfSelf >= sizeof(StoreStress) );
	self = (StoreStress*) _Rheology_New(  RHEOLOGY_PASSARGS  );
	
	return self;
}

void _StoreStress_Init(
	StoreStress*			self,
	MaterialPointsSwarm*	materialPointsSwarm,
	FeVariable*				strainRateField )
{
	Name                       variableName[6];
	StandardParticle           particle;
	StoreStress_ParticleExt*   particleExt;
	Index                      variable_I;
	Dimension_Index            dim = materialPointsSwarm->dim;
		
	/* Assign Pointers */
	self->materialPointsSwarm = materialPointsSwarm;
	self->strainRateField = strainRateField;

	self->particleExtHandle = ExtensionManager_Add( materialPointsSwarm->particleExtensionMgr, (Name)self->type, sizeof( StoreStress_ParticleExt )  );
	
	/* Add SwarmVariables for plotting */
	particleExt = ExtensionManager_Get( materialPointsSwarm->particleExtensionMgr, &particle, self->particleExtHandle );

	if ( dim == 2 ) {
		variableName[0] = StG_Strdup( "tau_xx" );
		variableName[1] = StG_Strdup( "tau_yy" );
		variableName[2] = StG_Strdup( "tau_xy" );
	}
	else {
		variableName[0] = StG_Strdup( "tau_xx" );
		variableName[1] = StG_Strdup( "tau_yy" );
		variableName[2] = StG_Strdup( "tau_zz" );
		variableName[3] = StG_Strdup( "tau_xy" );
		variableName[4] = StG_Strdup( "tau_xz" );
		variableName[5] = StG_Strdup( "tau_yz" );
	}
	
	self->materialPointsSwarmVariable = Swarm_NewVectorVariable( materialPointsSwarm, (Name)"StressTensor", (ArithPointer) &particleExt->stress - (ArithPointer) &particle, 
		Variable_DataType_Double, 
		StGermain_nSymmetricTensorVectorComponents(dim), 
		variableName[0],variableName[1],variableName[2],
		variableName[3],variableName[4],variableName[5]);
	
	for( variable_I = 0; variable_I < StGermain_nSymmetricTensorVectorComponents(dim) ; variable_I++ )
		Memory_Free( variableName[ variable_I ] );

}

void* _StoreStress_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                                     _sizeOfSelf = sizeof(StoreStress);
	Type                                                             type = StoreStress_Type;
	Stg_Class_DeleteFunction*                                     _delete = _Rheology_Delete;
	Stg_Class_PrintFunction*                                       _print = _Rheology_Print;
	Stg_Class_CopyFunction*                                         _copy = _Rheology_Copy;
	Stg_Component_DefaultConstructorFunction*         _defaultConstructor = _StoreStress_DefaultNew;
	Stg_Component_ConstructFunction*                           _construct = _StoreStress_AssignFromXML;
	Stg_Component_BuildFunction*                                   _build = _StoreStress_Build;
	Stg_Component_InitialiseFunction*                         _initialise = _StoreStress_Initialise;
	Stg_Component_ExecuteFunction*                               _execute = _Rheology_Execute;
	Stg_Component_DestroyFunction*                               _destroy = _StoreStress_Destroy;
	Rheology_ModifyConstitutiveMatrixFunction*  _modifyConstitutiveMatrix = _StoreStress_ModifyConstitutiveMatrix;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

	return (void*) _StoreStress_New(  STORESTRESS_PASSARGS  );
}

void _StoreStress_AssignFromXML( void* rheology, Stg_ComponentFactory* cf, void* data ){
	StoreStress*              self              = (StoreStress*)rheology;
	MaterialPointsSwarm*      materialPointsSwarm;
	FeVariable*               strainRateField;

	/* Construct Parent */
	_Rheology_AssignFromXML( self, cf, data );

	materialPointsSwarm = Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"MaterialPointsSwarm", MaterialPointsSwarm, True, data  );
	/* TODO : 'KeyFallback' soon to be deprecated/updated */
	strainRateField = Stg_ComponentFactory_ConstructByNameWithKeyFallback( cf, self->name, (Name)"StrainRateField", (Dictionary_Entry_Key)"StrainRateField", FeVariable, True, data  );
	/*
	strainRateField = Stg_ComponentFactory_ConstructByKey( cf, self->name, 
			"StrainRateField", FeVariable, True );
	*/

	_StoreStress_Init( self, materialPointsSwarm, strainRateField );
}

void _StoreStress_Build( void* _self, void* data ) {
	StoreStress*  self               = (StoreStress*) _self;

	/* Build parent */
	_Rheology_Build( self, data );

	Stg_Component_Build(	self->strainRateField, data, False);	
	Stg_Component_Build(	self->materialPointsSwarm, data, False);
	Stg_Component_Build(	self->materialPointsSwarmVariable, data, False);

}

void _StoreStress_Initialise( void* _self, void* data ) {
	StoreStress*  self               = (StoreStress*) _self;

	/* Initialise parent */
	_Rheology_Initialise( self, data );

	Stg_Component_Initialise(	self->strainRateField, data, False);	
	Stg_Component_Initialise(	self->materialPointsSwarm, data, False);
	Stg_Component_Initialise(	self->materialPointsSwarmVariable, data, False);

}

void _StoreStress_Destroy( void* _self, void* data ) {
	StoreStress*  self               = (StoreStress*) _self;

	Stg_Component_Destroy(	self->strainRateField, data, False);	
	Stg_Component_Destroy(	self->materialPointsSwarm, data, False);
	Stg_Component_Destroy(	self->materialPointsSwarmVariable, data, False);

	/* Destroy parent */
	_Rheology_Destroy( self, data );

}

void _StoreStress_ModifyConstitutiveMatrix( 
		void*                                              rheology, 
		ConstitutiveMatrix*                                constitutiveMatrix,
		MaterialPointsSwarm*                               materialPointsSwarm,
		Element_LocalIndex                                 lElement_I,
		MaterialPoint*                                     materialPoint,
		Coord                                              xi )
{
	StoreStress*	                  self              = (StoreStress*) rheology;
	StoreStress_ParticleExt*          particleExt;
	SymmetricTensor                   strainRate;

	/*
	** HAXOR: Throwing this flag in here to try and prevent REP from
	** overwriting the previously calculated viscosity value. */
	if ( !constitutiveMatrix->previousSolutionExists )
	  return;

	particleExt      = ExtensionManager_Get( materialPointsSwarm->particleExtensionMgr, materialPoint, self->particleExtHandle );

	FeVariable_InterpolateWithinElement( self->strainRateField, lElement_I, xi, strainRate );
	ConstitutiveMatrix_CalculateStress( constitutiveMatrix, strainRate, particleExt->stress );
}


