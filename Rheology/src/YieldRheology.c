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
** $Id: YieldRheology.c 779 2008-08-06 15:50:41Z LukeHodkinson $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>

#include "types.h"
#include "ConstitutiveMatrix.h"
#include "RheologyClass.h"
#include "StrainWeakening.h"
#include "YieldRheology.h"

#include <assert.h>

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type YieldRheology_Type = "YieldRheology";

/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
YieldRheology* _YieldRheology_New( 
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
		YieldRheology_GetYieldCriterionFunction*           _getYieldCriterion,
		YieldRheology_GetYieldIndicatorFunction*           _getYieldIndicator,
		YieldRheology_HasYieldedFunction*                  _hasYielded,
		Name                                               name ) 
{
	YieldRheology*					self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( sizeOfSelf >= sizeof(YieldRheology) );
	self = (YieldRheology*) _Rheology_New( 
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
	
	/* Function pointers for this class that are not on the parent class should be set here */
	self->_getYieldCriterion = _getYieldCriterion;
	self->_getYieldIndicator = _getYieldIndicator;
	self->_hasYielded        = _hasYielded;
	
	return self;
}

void _YieldRheology_Init( YieldRheology* self, StrainWeakening* strainWeakening, MaterialPointsSwarm* materialPointsSwarm )
{
	ExtensionInfo_Index     handle             = (ExtensionInfo_Index) -1;
	self->strainWeakening = strainWeakening;

	Rheology_SetToNonLinear( self );
	
	if ( materialPointsSwarm ) {
		ArithPointer offset; 
		
		/* See if the YieldRheology Type has already added an extension to the particle 
		 * If handle is given a value of '-1' - that means that no extension has been added
		 * with the YieldRheology_Type
		 * We should then add the extension */
		
		handle = ExtensionManager_GetHandle( materialPointsSwarm->particleExtensionMgr, YieldRheology_Type );
		
		if ( handle == (ExtensionInfo_Index) -1 ) {
			handle = ExtensionManager_Add(materialPointsSwarm->particleExtensionMgr, YieldRheology_Type, sizeof(Particle_Bool));

			/* Adding variable for plotting purpose */
			offset = (ArithPointer) 
				ExtensionManager_Get( materialPointsSwarm->particleExtensionMgr, NULL, handle );

			self->hasYieldedVariable = Swarm_NewScalarVariable(
				materialPointsSwarm,
				"HasYielded",
				offset,
				Variable_DataType_Char );
		}
		else {
			/* if the variable has already been created - then just store the pointer */
			Name hasYieldedVariableName = Stg_Object_AppendSuffix( materialPointsSwarm, "HasYielded" );
			self->hasYieldedVariable = SwarmVariable_Register_GetByName( materialPointsSwarm->swarmVariable_Register, hasYieldedVariableName );
			Memory_Free( hasYieldedVariableName );
		}
	}
	
	/* Store value of particle extension handle
	 * if there are no material points - this will be '-1'
	 * if there are material points - all YieldRheology objects will refer to the same handle */
	self->hasYieldedParticleExtHandle = handle;

        self->yieldCriterion = 0.0;
        self->minVisc = 0.0;
}


void _YieldRheology_Delete( void* rheology ) {
	YieldRheology*					self = (YieldRheology*)rheology;
	_Stg_Component_Delete( self );
}
void _YieldRheology_Print( void* rheology, Stream* stream ) {}

void* _YieldRheology_Copy( void* rheology, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	YieldRheology*	self = (YieldRheology*)rheology;

	/* TODO */ abort();
	return (void*) self;
}

void _YieldRheology_AssignFromXML( void* rheology, Stg_ComponentFactory* cf, void* data ){
	YieldRheology*           self                 = (YieldRheology*)rheology;
	StrainWeakening*         strainWeakening;
	MaterialPointsSwarm*  materialPoints;
	
	_Rheology_AssignFromXML( self, cf, data );

	strainWeakening =  Stg_ComponentFactory_ConstructByKey( 
		cf, 
		self->name,  
		"StrainWeakening", 
		StrainWeakening,  
		False,
		data  ) ;

	materialPoints =  Stg_ComponentFactory_ConstructByKey(  
		cf,  
		self->name,  
		"MaterialPointsSwarm", 
		MaterialPointsSwarm,  
		False,
		data  ) ;

	_YieldRheology_Init( self, strainWeakening, materialPoints );

	self->minVisc = Stg_ComponentFactory_GetDouble( cf, self->name, "minimumViscosity", 0.0 );
}

void _YieldRheology_Build( void* rheology, void* data ) {
	YieldRheology*                   self          = (YieldRheology*) rheology;

	_Rheology_Build( rheology, data );

	/* This variable only needs to be built if there are material points (hasYieldedVariable is created
	 * in _YieldRheology_Init only in that case) */
	
	if ( self->hasYieldedVariable ) {
		Stg_Component_Build( self->hasYieldedVariable, data, False );
	}
}
void _YieldRheology_Initialise( void* rheology, void* data ) {
	YieldRheology*                   self          = (YieldRheology*) rheology;
	
	_Rheology_Initialise( rheology, data );

	/* This variable only needs to be initialised if there are material points (hasYieldedVariable is created
	 * in _YieldRheology_Init only in that case) */
	if ( self->hasYieldedVariable ) {
		Stg_Component_Initialise( self->hasYieldedVariable, data, False );
	}
}
void _YieldRheology_Execute( void* rheology, void* data ) {}
void _YieldRheology_Destroy( void* rheology, void* data ) {}

void _YieldRheology_ModifyConstitutiveMatrix( 
		void*                                              yieldRheology, 
		ConstitutiveMatrix*                                constitutiveMatrix,
		MaterialPointsSwarm*                               materialPointsSwarm,
		Element_LocalIndex                                 lElement_I,
		MaterialPoint*                                     materialPoint,
		Coord                                              xi )
{
	YieldRheology*                   self          = (YieldRheology*) yieldRheology;
	double                           yieldCriterion;
	double                           yieldIndicator; /* A particle will yield if yieldCriterion < yieldIndicator */
	
	/* Don't want to yield on the first ever solve */
	if ( !constitutiveMatrix->previousSolutionExists )
		return;
		
	yieldIndicator = YieldRheology_GetYieldIndicator( self, constitutiveMatrix, materialPointsSwarm, lElement_I, materialPoint, xi );
	yieldCriterion = YieldRheology_GetYieldCriterion( self, constitutiveMatrix, materialPointsSwarm, lElement_I, materialPoint, xi );
	
	/* This function is called before YieldRheology_HasYielded  because :
	 * If a particle has not yielded, we still need to update the strain weakening for it (it will heal
	 * in that case rather than weaken)
	 * If a particle has yielded, the initial viscosity is needed to update the strain weakening. But since this 
	 * viscosity can be modified in the constitutive matrix by YieldRheology_HasYielded function, we need to do the update 
	 * before */
	
	if ( self->strainWeakening ) {
		StrainWeakening_AssignIncrement( self->strainWeakening, constitutiveMatrix, materialPointsSwarm,
			lElement_I, materialPoint, yieldCriterion, yieldIndicator );
	}

	/* Set a bool to TRUE or FLAG depending on whether a particle has failed or not */
	YieldRheology_SetParticleFlag( self, materialPointsSwarm, materialPoint, (yieldCriterion < yieldIndicator) );

	/* Test to see if it has yielded */
	if ( yieldCriterion < yieldIndicator ) {
		YieldRheology_HasYielded( self, constitutiveMatrix, materialPointsSwarm, lElement_I, materialPoint, yieldCriterion, yieldIndicator );
	}
}

