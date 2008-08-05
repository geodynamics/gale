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
** $Id: BuiterStrainWeakening.c 98 2005-12-15 12:56:32Z VincentLemiale $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/* Modified 2006 Walter Landry to implement Buiter Strain Weakening */



#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>

#include "types.h"
#include "RheologyClass.h"
#include "YieldRheology.h"
#include "VonMises.h"
#include "StrainWeakening.h"
#include "BuiterStrainWeakening.h"
#include "ConstitutiveMatrix.h"

#include <assert.h>

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type BuiterStrainWeakening_Type = "BuiterStrainWeakening";

/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
BuiterStrainWeakening* _BuiterStrainWeakening_New( 
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
		TimeIntegratee_CalculateTimeDerivFunction*         _calculateTimeDeriv,
		TimeIntegratee_IntermediateFunction*               _intermediate,
		StrainWeakening_CalcIncrementFunction*             _calcIncrement,
		Name                                               name ) 
{
	BuiterStrainWeakening*					self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( sizeOfSelf >= sizeof(BuiterStrainWeakening) );
	self = (BuiterStrainWeakening*) _StrainWeakening_New( 
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
			_calculateTimeDeriv,
			_intermediate,
			_calcIncrement, 
			name );
	
	/* Function pointers for this class that are not on the parent class should be set here */
	
	return self;
}

void _BuiterStrainWeakening_Init(
		BuiterStrainWeakening*                             self,
		MaterialPointsSwarm*                               swarm,
		double                                             healingRate,
		double                                             softeningStrain,
		double                                             initialDamageFraction,
		double                                             initialDamageWavenumber,
		double                                             initialDamageWavenumberSinI,
		double                                             initialDamageWavenumberCosI,
		double                                             initialDamageWavenumberSinJ,
		double                                             initialDamageWavenumberCosJ,
		double                                             initialDamageFactor,
		long int                                           randomSeed,
		Stg_Shape*                                         initialStrainShape  )
{
	_StrainWeakening_Init( (StrainWeakening*)self, 
				swarm, healingRate, 
				softeningStrain, 
				initialDamageFraction,
				initialDamageWavenumber,
				initialDamageWavenumberSinI,
				initialDamageWavenumberCosI,
				initialDamageWavenumberSinJ,
				initialDamageWavenumberCosJ,
				initialDamageFactor, randomSeed,
				initialStrainShape );
}

void* _BuiterStrainWeakening_DefaultNew( Name name ) {
	return (void*) _BuiterStrainWeakening_New(
		sizeof(BuiterStrainWeakening),
		BuiterStrainWeakening_Type,
		_TimeIntegratee_Delete,
		_TimeIntegratee_Print,
		_TimeIntegratee_Copy,
		_BuiterStrainWeakening_DefaultNew,
		_BuiterStrainWeakening_Construct,
		_BuiterStrainWeakening_Build,
		_BuiterStrainWeakening_Initialise,
		_TimeIntegratee_Execute,
		_TimeIntegratee_Destroy,
		_StrainWeakening_TimeDerivative,
		_TimeIntegratee_Intermediate,
		_StrainWeakening_CalcIncrementIsotropic,
		name );
}

void _BuiterStrainWeakening_Construct( void* strainWeakening, Stg_ComponentFactory* cf, void* data ){
	BuiterStrainWeakening*        self           = (BuiterStrainWeakening*) strainWeakening;

	/* Construct Parent */
	_StrainWeakening_Construct( self, cf, data );
}

void _BuiterStrainWeakening_Build( void* strainWeakening, void* data ) {
	BuiterStrainWeakening*                       self               = (BuiterStrainWeakening*) strainWeakening;

	/* Build parent */
	_StrainWeakening_Build( self, data );
}

void _BuiterStrainWeakening_Initialise( void* strainWeakening, void* data ) {
	BuiterStrainWeakening*                       self               = (BuiterStrainWeakening*) strainWeakening;
	
	/* Initialise Parent */
	_StrainWeakening_Initialise( self, data );
}
