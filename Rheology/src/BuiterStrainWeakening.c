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

/* Public Constructor */
BuiterStrainWeakening* BuiterStrainWeakening_New(
      Name                                               name,
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
   BuiterStrainWeakening* self = (BuiterStrainWeakening*) _BuiterStrainWeakening_DefaultNew( name );

   _BuiterStrainWeakening_Init(
	       self,
	       swarm,
	       healingRate,
	       softeningStrain,
	       initialDamageFraction,
	       initialDamageWavenumber,
	       initialDamageWavenumberSinI,
	       initialDamageWavenumberCosI,
	       initialDamageWavenumberSinJ,
	       initialDamageWavenumberCosJ,
	       initialDamageFactor,
	       randomSeed,
	       initialStrainShape  );
   self->isConstructed = True;
   return self;
}

/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
BuiterStrainWeakening* _BuiterStrainWeakening_New(  BUITERSTRAINWEAKENING_DEFARGS  ) 
{
	BuiterStrainWeakening*					self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( _sizeOfSelf >= sizeof(BuiterStrainWeakening) );
	self = (BuiterStrainWeakening*) _StrainWeakening_New(  STRAINWEAKENING_PASSARGS  );
	
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
	/* Variables set in this function */
	SizeT                                               _sizeOfSelf = sizeof(BuiterStrainWeakening);
	Type                                                       type = BuiterStrainWeakening_Type;
	Stg_Class_DeleteFunction*                               _delete = _TimeIntegratee_Delete;
	Stg_Class_PrintFunction*                                 _print = _TimeIntegratee_Print;
	Stg_Class_CopyFunction*                                   _copy = _TimeIntegratee_Copy;
	Stg_Component_DefaultConstructorFunction*   _defaultConstructor = _BuiterStrainWeakening_DefaultNew;
	Stg_Component_ConstructFunction*                     _construct = _StrainWeakening_AssignFromXML;
	Stg_Component_BuildFunction*                             _build = _StrainWeakening_Build;
	Stg_Component_InitialiseFunction*                   _initialise = _StrainWeakening_Initialise;
	Stg_Component_ExecuteFunction*                         _execute = _TimeIntegratee_Execute;
	Stg_Component_DestroyFunction*                         _destroy = _TimeIntegratee_Destroy;
	TimeIntegratee_CalculateTimeDerivFunction*  _calculateTimeDeriv = _StrainWeakening_TimeDerivative;
	TimeIntegratee_IntermediateFunction*              _intermediate = _TimeIntegratee_Intermediate;
	StrainWeakening_CalcIncrementFunction*           _calcIncrement = _StrainWeakening_CalcIncrementIsotropic;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = ZERO;

	return (void*) _BuiterStrainWeakening_New(  BUITERSTRAINWEAKENING_PASSARGS  );
}



