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
** $Id: VonMises.c 803 2008-09-11 05:22:20Z LukeHodkinson $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>

#include "types.h"
#include "RheologyClass.h"
#include "StrainWeakening.h"
#include "YieldRheology.h"
#include "VonMises.h"
#include "ConstitutiveMatrix.h"

#include <assert.h>

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type VonMises_Type = "VonMises";

/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
VonMises* _VonMises_New( 
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
	VonMises*					self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( sizeOfSelf >= sizeof(VonMises) );
	self = (VonMises*) _YieldRheology_New( 
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
			_getYieldCriterion,
			_getYieldIndicator,
			_hasYielded,
			name );
	
	/* Function pointers for this class that are not on the parent class should be set here */
	
	return self;
}

void _VonMises_Init( 
		VonMises*   self, 
		FeVariable* strainRateField, 
		double      cohesion, 
		double      cohesionAfterSoftening,
		Bool        strainRateSoftening)
{
	self->strainRateField        = strainRateField;
	self->cohesion               = cohesion;
	
	/* Strain softening of Cohesion - (linear weakening is assumed) */
	/* needs a softening factor between +0 and 1 and a reference strain > 0 */
	self->cohesionAfterSoftening = cohesionAfterSoftening;
	self->strainRateSoftening    = strainRateSoftening;
}

void* _VonMises_DefaultNew( Name name ) {
	return (void*) _VonMises_New(
		sizeof(VonMises),
		VonMises_Type,
		_YieldRheology_Delete,
		_YieldRheology_Print,
		_YieldRheology_Copy,
		_VonMises_DefaultNew,
		_VonMises_Construct,
		_YieldRheology_Build,
		_YieldRheology_Initialise,
		_YieldRheology_Execute,
		_YieldRheology_Destroy,
		_YieldRheology_ModifyConstitutiveMatrix,
		_VonMises_GetYieldCriterion,
		_VonMises_GetYieldIndicator,
		_VonMises_HasYielded,
		name );
}

void _VonMises_Construct( void* rheology, Stg_ComponentFactory* cf, void* data ){
	VonMises*          self           = (VonMises*)rheology;
	FeVariable*        strainRateField;

	/* Construct Parent */
	_YieldRheology_Construct( self, cf, data );
	
	strainRateField = Stg_ComponentFactory_ConstructByKey(  
			cf, self->name, "StrainRateField", FeVariable, False, data );
	self->swarmStrainRate = Stg_ComponentFactory_ConstructByKey(
	   cf, self->name, "swarmStrainRate", SwarmVariable, False, data );
	assert( strainRateField || self->swarmStrainRate );
	
	_VonMises_Init( 
			self, 
			strainRateField,
			Stg_ComponentFactory_GetDouble( cf, self->name, "cohesion", 0.0 ),
			Stg_ComponentFactory_GetDouble( cf, self->name, "cohesionAfterSoftening", 0.0 ),
			Stg_ComponentFactory_GetBool(   cf, self->name, "strainRateSoftening", False ) );
}

double _VonMises_GetYieldCriterion( 
		void*                            rheology,
		ConstitutiveMatrix*              constitutiveMatrix,
		MaterialPointsSwarm*             materialPointsSwarm,
		Element_LocalIndex               lElement_I,
		MaterialPoint*                   materialPoint,
		Coord                            xi )
{
	VonMises*                         self             = (VonMises*) rheology;
	double                            strainWeakeningRatio;
	double                            effectiveCohesion;
	
	/* Strain softening of yield stress */
	strainWeakeningRatio = StrainWeakening_CalcRatio( self->strainWeakening, materialPoint );
		
	effectiveCohesion =  self->cohesion * (1.0 - strainWeakeningRatio) + 
			self->cohesionAfterSoftening * strainWeakeningRatio;

	return effectiveCohesion;
}

double _VonMises_GetYieldIndicator( 
		void*                            rheology,
		ConstitutiveMatrix*              constitutiveMatrix,
		MaterialPointsSwarm*             materialPointsSwarm,
		Element_LocalIndex               lElement_I,
		MaterialPoint*                   materialPoint,
		Coord                            xi )
{
	VonMises*                         self             = (VonMises*) rheology;
	SymmetricTensor                   stress;
	SymmetricTensor                   strainRate;
	double                            stressInv;
	
	/* Get Strain Rate */
	if( self->strainRateField ) {
	   FeVariable_InterpolateWithinElement(
	      self->strainRateField, lElement_I, xi, strainRate );
	}
	else {
	   SwarmVariable_ValueAt( self->swarmStrainRate,
				  constitutiveMatrix->currentParticleIndex,
				  strainRate );
	}

	/* Get Stress */
	ConstitutiveMatrix_CalculateStress( constitutiveMatrix, strainRate, stress );
	stressInv = SymmetricTensor_2ndInvariant( stress, constitutiveMatrix->dim );

	return stressInv;
}

void _VonMises_HasYielded( 
		void*                            rheology,
		ConstitutiveMatrix*              constitutiveMatrix,
		MaterialPointsSwarm*             materialPointsSwarm,
		Element_LocalIndex               lElement_I,
		MaterialPoint*                   materialPoint,
		double                           yieldCriterion,
		double                           yieldIndicator )
{
	VonMises*                 self             = (VonMises*) rheology;
	double                    viscosity        = ConstitutiveMatrix_GetIsotropicViscosity( constitutiveMatrix );

	if ( self->strainRateSoftening ) {
		/* TODO - Put this function in a plugin */
		viscosity = 2.0 * yieldCriterion * yieldCriterion * viscosity / 
				( yieldCriterion * yieldCriterion + yieldIndicator * yieldIndicator );
		ConstitutiveMatrix_SetIsotropicViscosity( constitutiveMatrix, viscosity );
	}	
	else {
		double beta = 1.0 - yieldCriterion/yieldIndicator;
		double corr = -viscosity * beta;

		if( (viscosity + corr) < self->minVisc )
		   corr = self->minVisc - viscosity;
		ConstitutiveMatrix_IsotropicCorrection( constitutiveMatrix, corr );
	}
}

