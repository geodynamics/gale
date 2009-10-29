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
StoreStress* _StoreStress_New( 
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
	StoreStress*					self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( sizeOfSelf >= sizeof(StoreStress) );
	self = (StoreStress*) _Rheology_New( 
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

void _StoreStress_Init(
		StoreStress*                                       self,
		MaterialPointsSwarm*                               materialPointsSwarm,
		FeVariable*                                        strainRateField )
{
	Name                       variableName[6];
	SwarmVariable*             materialPointsSwarmVariable;
	StandardParticle           particle;
	StoreStress_ParticleExt*   particleExt;
	Index                      variable_I;
	Dimension_Index            dim            = materialPointsSwarm->dim;
		
	/* Assign Pointers */
	self->materialPointsSwarm           = materialPointsSwarm;
	self->strainRateField = strainRateField;

	self->particleExtHandle = 
		ExtensionManager_Add( materialPointsSwarm->particleExtensionMgr, self->type, sizeof( StoreStress_ParticleExt ) );
	
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
	
	materialPointsSwarmVariable = Swarm_NewVectorVariable(
		materialPointsSwarm,
		"StressTensor",
		(ArithPointer) &particleExt->stress - (ArithPointer) &particle, 
		Variable_DataType_Double, 
		StGermain_nSymmetricTensorVectorComponents(dim), 
		variableName[0],variableName[1],variableName[2],
		variableName[3],variableName[4],variableName[5]);
	
	for( variable_I = 0; variable_I < StGermain_nSymmetricTensorVectorComponents(dim) ; variable_I++ )
		Memory_Free( variableName[ variable_I ] );

}

void* _StoreStress_DefaultNew( Name name ) {
	return (void*) _StoreStress_New(
		sizeof(StoreStress),
		StoreStress_Type,
		_Rheology_Delete,
		_Rheology_Print,
		_Rheology_Copy,
		_StoreStress_DefaultNew,
		_StoreStress_AssignFromXML,
		_Rheology_Build,
		_Rheology_Initialise,
		_Rheology_Execute,
		_Rheology_Destroy,
		_StoreStress_ModifyConstitutiveMatrix,
		name );
}

void _StoreStress_AssignFromXML( void* rheology, Stg_ComponentFactory* cf, void* data ){
	StoreStress*              self              = (StoreStress*)rheology;
	MaterialPointsSwarm*      materialPointsSwarm;
	FeVariable*               strainRateField;

	/* Construct Parent */
	_Rheology_AssignFromXML( self, cf, data );

	materialPointsSwarm = Stg_ComponentFactory_ConstructByKey( 
		cf, 
		self->name, 
		"MaterialPointsSwarm", 
		MaterialPointsSwarm, 
		True,
		data );
	/* TODO : 'KeyFallback' soon to be deprecated/updated */
	strainRateField = Stg_ComponentFactory_ConstructByNameWithKeyFallback( 
		cf, 
		self->name,
                "StrainRateField", 
		"StrainRateField", 
		FeVariable, 
		True,
		data );
	/*
	strainRateField = Stg_ComponentFactory_ConstructByKey( cf, self->name, 
			"StrainRateField", FeVariable, True );
	*/

	_StoreStress_Init( self, materialPointsSwarm, strainRateField );
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
