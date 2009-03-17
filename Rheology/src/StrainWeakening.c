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
** $Id: StrainWeakening.c 776 2008-08-05 02:39:26Z LouisMoresi $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#include <stdlib.h>
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
#include "ConstitutiveMatrix.h"

#include <assert.h>

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type StrainWeakening_Type = "StrainWeakening";

/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
StrainWeakening* _StrainWeakening_New( 
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
	StrainWeakening*					self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( sizeOfSelf >= sizeof(StrainWeakening) );
	self = (StrainWeakening*) _TimeIntegratee_New( 
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
			name );
	
	/* Function pointers for this class that are not on the parent class should be set here */
	self->_calcIncrement = _calcIncrement;
	
	return self;
}

void _StrainWeakening_Init(
		StrainWeakening*                                   self,
		MaterialPointsSwarm*                               swarm,
		double                                             healingRate,
		double                                             softeningStrain,
		double                                             initialDamageFraction,
		double                                             initialDamageWavenumber,
		double                                             initialDamageWavenumberSinI,
		double                                             initialDamageWavenumberCosI,
		double                                             initialDamageWavenumberSinK,
		double                                             initialDamageWavenumberCosK,
		double                                             initialDamageFactor,
		long int                                           randomSeed,
		Stg_Shape*                                         initialStrainShape )
{
	/* Assign Values */
	self->swarm                    = swarm;
	self->healingRate              = healingRate;
	self->softeningStrain          = softeningStrain;
	self->initialDamageFraction    = initialDamageFraction;
	self->initialDamageWavenumber  = initialDamageWavenumber;
	self->initialDamageWavenumberSinI = initialDamageWavenumberSinI;
	self->initialDamageWavenumberCosI = initialDamageWavenumberCosI;
	self->initialDamageWavenumberSinK = initialDamageWavenumberSinK;
	self->initialDamageWavenumberCosK = initialDamageWavenumberCosK;
	self->initialDamageFactor      = initialDamageFactor;
	self->randomSeed               = randomSeed;
	self->initialStrainShape       = initialStrainShape;
	
	/****** Setup Variables *****/

	/* First check to see if a particle extension has already been created for this swarm */
	self->particleExtHandle = ExtensionManager_GetHandle( swarm->particleExtensionMgr, StrainWeakening_Type );

	/* If there isn't one then create the particle extension - otherwise just use the one already there*/
	if ( self->particleExtHandle == (ExtensionInfo_Index) -1 ) {
		StandardParticle             particle;
		StrainWeakening_ParticleExt* particleExt;

		/* Add particle extension */
		self->particleExtHandle = 
			ExtensionManager_Add( swarm->particleExtensionMgr, StrainWeakening_Type, sizeof(StrainWeakening_ParticleExt) );	

		particleExt = ExtensionManager_Get( swarm->particleExtensionMgr, &particle, self->particleExtHandle );

		/*Add variables for vizualization / analysis purposes */
		
		
		
		self->postFailureWeakening = Swarm_NewScalarVariable(
			swarm,
			"PostFailureWeakening",
			(ArithPointer) &particleExt->postFailureWeakening - (ArithPointer) &particle,
			Variable_DataType_Double );
		
		self->postFailureWeakeningIncrement = Swarm_NewScalarVariable(
			swarm,
			"PostFailureWeakeningIncrement",
			(ArithPointer) &particleExt->postFailureWeakeningIncrement - (ArithPointer) &particle,
			Variable_DataType_Double );
	}
	else {
		Name variableName;

		/* Get Variables already created */
		variableName = Stg_Object_AppendSuffix( swarm, "PostFailureWeakening" );
		self->postFailureWeakening = SwarmVariable_Register_GetByName( swarm->swarmVariable_Register, variableName );
		assert( self->postFailureWeakening );
		Memory_Free( variableName );
		
		variableName = Stg_Object_AppendSuffix( swarm, "PostFailureWeakeningIncrement" );
		self->postFailureWeakeningIncrement = SwarmVariable_Register_GetByName( swarm->swarmVariable_Register, variableName );
		assert( self->postFailureWeakeningIncrement );
		Memory_Free( variableName );
	}
	
	/* The strain weakening class inherits from the TimeIntegratee class - this class needs a 'Variable' to 
	 * integrate through time. For the StrainWeakening component this variable is the 'PostFailureWeakening'
	 * we need to set this explicitly here */
	self->variable = self->postFailureWeakening->variable;
	
	/* Add function to entry point which is called at the end of the time integration steps - 
	 * this will make all the negative values zero */
	TimeIntegrator_AppendFinishEP( self->timeIntegrator,
		"StrainWeakening_MakeValuesPositive", _StrainWeakening_MakeValuesPositive, self->name, self );

}

void* _StrainWeakening_DefaultNew( Name name ) {
	return (void*) _StrainWeakening_New(
			sizeof(StrainWeakening),
		StrainWeakening_Type,
		_TimeIntegratee_Delete,
		_TimeIntegratee_Print,
		_TimeIntegratee_Copy,
		_StrainWeakening_DefaultNew,
		_StrainWeakening_Construct,
		_StrainWeakening_Build,
		_StrainWeakening_Initialise,
		_TimeIntegratee_Execute,
		_TimeIntegratee_Destroy,
		_StrainWeakening_TimeDerivative,
		_TimeIntegratee_Intermediate,
		_StrainWeakening_CalcIncrementIsotropic,
		name );
}

void _StrainWeakening_Construct( void* strainWeakening, Stg_ComponentFactory* cf, void* data ){
	StrainWeakening*        self           = (StrainWeakening*) strainWeakening;
	MaterialPointsSwarm*    materialPointsSwarm;
	double                  healingRate;
	double                  softeningStrain;
	double                  initialDamageFraction;
	double                  initialDamageWavenumber;
	double                  initialDamageWavenumberSinI;
	double                  initialDamageWavenumberCosI;
	double                  initialDamageWavenumberSinK;
	double                  initialDamageWavenumberCosK;
	double                  initialDamageFactor;
	long int                randomSeed;
	Stg_Shape*              initialStrainShape;

	/* Construct Parent */
	_TimeIntegratee_Construct( self, cf, data );
	
	materialPointsSwarm     = Stg_ComponentFactory_ConstructByKey(
		cf, 
		self->name, 
		"MaterialPointsSwarm", 
		MaterialPointsSwarm, 
		True,
		data );

	healingRate              = Stg_ComponentFactory_GetDouble( cf, self->name, "healingRate",              0.0 );
	softeningStrain          = Stg_ComponentFactory_GetDouble( cf, self->name, "softeningStrain",          HUGE_VAL );
	initialDamageFraction    = Stg_ComponentFactory_GetDouble( cf, self->name, "initialDamageFraction",    0.0 );
	initialDamageWavenumber  = Stg_ComponentFactory_GetDouble( cf, self->name, "initialDamageWavenumber",  -1.0 );
	initialDamageWavenumberSinI = Stg_ComponentFactory_GetDouble( cf, self->name, "initialDamageWavenumberSinI", -1.0 );
	initialDamageWavenumberCosI = Stg_ComponentFactory_GetDouble( cf, self->name, "initialDamageWavenumberCosI", -1.0 );
	initialDamageWavenumberSinK = Stg_ComponentFactory_GetDouble( cf, self->name, "initialDamageWavenumberSinK", -1.0 );
	initialDamageWavenumberCosK = Stg_ComponentFactory_GetDouble( cf, self->name, "initialDamageWavenumberCosK", -1.0 );
	initialDamageFactor      = Stg_ComponentFactory_GetDouble( cf, self->name, "initialDamageFactor",       1.0 );
	randomSeed               = (long int) Stg_ComponentFactory_GetInt( cf, self->name, "randomSeed",        0 );
	initialStrainShape       = Stg_ComponentFactory_ConstructByKey( cf, self->name, "initialStrainShape", Stg_Shape, False, data );

	_StrainWeakening_Init(
			self, 
			materialPointsSwarm, 
			healingRate,
			softeningStrain,
			initialDamageFraction,
			initialDamageWavenumber, 
			initialDamageWavenumberSinI, 
			initialDamageWavenumberCosI, 
			initialDamageWavenumberSinK, 
			initialDamageWavenumberCosK, 
			initialDamageFactor,
			randomSeed,
			initialStrainShape );
}

void _StrainWeakening_Build( void* strainWeakening, void* data ) {
	StrainWeakening*                       self               = (StrainWeakening*) strainWeakening;

	/* Build parent */
	_TimeIntegratee_Build( self, data );

	Stg_Component_Build( self->postFailureWeakeningIncrement, data, False );
	/* The postFailureWeakening doesn't need to be built here because it has already been
	 * built in the TimeIntegratee class
	 * (self->variable = self->postFailureWeakening->variable in _StrainWeakening_Init function) */
}

void _StrainWeakening_Initialise( void* strainWeakening, void* data ) {
	StrainWeakening*                       self  = (StrainWeakening*) strainWeakening;
	Particle_Index                         lParticle_I;
	Particle_Index                         particleLocalCount;
	Variable*                              positionVariable   = self->swarm->particleCoordVariable->variable;
	double                                 postFailureWeakening;
	double*                                coord;
	AbstractContext*                       context = (AbstractContext*)data;

	int myrank;
	
	/* Initialise Parent */
	_TimeIntegratee_Initialise( self, data );

	Stg_Component_Initialise( self->swarm, data, False );
	Stg_Component_Initialise( self->postFailureWeakeningIncrement, data, False );

	/* Update variables */
	Variable_Update( positionVariable );
	Variable_Update( self->variable );
	Variable_Update( self->postFailureWeakeningIncrement->variable );

	particleLocalCount = self->variable->arraySize;

	/* We should only set initial conditions if in regular non-restart mode. If in restart mode, then
	the particle-based variables will be set correcty when we re-load the Swarm. */
	
	if ( !(context && (True == context->loadFromCheckPoint)) ) {
		/* Initialise random number generator */
		if(self->randomSeed==0) {
                    MPI_Comm_rank( context->communicator, &myrank);
		    srand( self->randomSeed );
		}
		else {
			srand( self->randomSeed );
		}
	
		for ( lParticle_I = 0 ; lParticle_I < particleLocalCount ; lParticle_I++ ) {
			/* Initialise Increment to Zero */
			Variable_SetValueDouble( self->postFailureWeakeningIncrement->variable, lParticle_I, 0.0 );
			Variable_SetValueDouble( self->variable, lParticle_I, 0.0 );

			/* Set initial damage parameter
			 * There is a certain fraction of the number of particles which are given initial strain */
			
			postFailureWeakening = 0.0;
			
			if ( rand() < RAND_MAX*self->initialDamageFraction ) {

				coord = Variable_GetPtrDouble( positionVariable, lParticle_I );
				
				if ( self->initialStrainShape && !Stg_Shape_IsCoordInside( self->initialStrainShape, coord ) ) {
					Variable_SetValueDouble( self->variable, lParticle_I, postFailureWeakening );
						continue;
				}
				
				postFailureWeakening = self->initialDamageFactor * rand() * self->softeningStrain/RAND_MAX;

				/* Modulate the initial weakening by a harmonic-squared function with wavenumber(s) specified by
					the user. */
				
				/* Use old definition if new one is not set */
				
				if ( self->initialDamageWavenumber > 0.0 && self->initialDamageWavenumberCosI == -1.0 ) {				
						postFailureWeakening *= 
							pow(0.5+0.5*cos(M_PI * coord[ I_AXIS ] * self->initialDamageWavenumber),2.0);
					}
					
				/* Alternate phase is appropriate for different bc's and choice of origin */	
					
				if ( self->initialDamageWavenumberCosI > 0.0 ) {				
						postFailureWeakening *= 
							pow(0.5+0.5*cos(M_PI * coord[ I_AXIS ] * self->initialDamageWavenumberCosI),2.0);
				}
				
				if ( self->initialDamageWavenumberSinI > 0.0 ) {				
					postFailureWeakening *= 
						pow(0.5+0.5*sin(M_PI * coord[ I_AXIS ] * self->initialDamageWavenumberSinI),2.0);
				}
				
				if ( self->initialDamageWavenumberCosK > 0.0 ) {				
					postFailureWeakening *= 
						pow(0.5+0.5*cos(M_PI * coord[ K_AXIS ] * self->initialDamageWavenumberCosK),2.0);
				}
				if ( self->initialDamageWavenumberSinK > 0.0 ) {				
					postFailureWeakening *= 
						pow(0.5+0.5*sin(M_PI * coord[ K_AXIS ] * self->initialDamageWavenumberSinK),2.0);
				}
			}
		
			Variable_SetValueDouble( self->variable, lParticle_I, postFailureWeakening );
		}
	}	
}

Bool _StrainWeakening_TimeDerivative( void* strainWeakening, Index lParticle_I, double* timeDeriv ) {
	StrainWeakening*  self   = (StrainWeakening*) strainWeakening;

	Variable_Update( self->postFailureWeakeningIncrement->variable );
	*timeDeriv = Variable_GetValueDouble( self->postFailureWeakeningIncrement->variable, lParticle_I );

	return True;
}

/* This function is called after the time integration steps - 
 * here we just want to make sure that each value is positive */
void _StrainWeakening_MakeValuesPositive( void* timeIntegrator, void* strainWeakening ) {
	StrainWeakening*   self   = (StrainWeakening*) strainWeakening;
	double*            value;
	Particle_Index     lParticle_I;
	Particle_Index     particleCount;

	Variable_Update( self->variable );
	
	particleCount = self->variable->arraySize;
	
	for ( lParticle_I = 0 ; lParticle_I < particleCount ; lParticle_I++ ) {
		value = Variable_GetPtrDouble( self->variable, lParticle_I );
		if (*value < 0.0)
			*value = 0.0;
	}
}

double _StrainWeakening_CalcIncrementIsotropic( 
		void*                            strainWeakening,
		ConstitutiveMatrix*              constitutiveMatrix,
		MaterialPointsSwarm*             swarm,
		Element_LocalIndex               lElement_I,
		MaterialPoint*                   particle,
		double                           yieldCriterion,
		double                           yieldIndicator )
{
	StrainWeakening*               self             = (StrainWeakening*) strainWeakening;
	StrainWeakening_ParticleExt*   particleExt;
	double                         beta             = 1.0 - yieldCriterion / yieldIndicator;
	double                         healingRate      = self->healingRate;
	double                         viscosity        = ConstitutiveMatrix_GetIsotropicViscosity( constitutiveMatrix );
	double 						   postFailureWeakening;
	
	
	particleExt = ExtensionManager_Get( self->swarm->particleExtensionMgr, particle, self->particleExtHandle );
	postFailureWeakening = particleExt->postFailureWeakening;
	
	if (beta < 0.0)
		beta = 0.0;

	return(    (yieldCriterion/viscosity) * ( beta/(1.0-beta)    /* growth term depends on plastic strain rate */
	 		 - healingRate * postFailureWeakening)               /* decay term */
			);
}

double StrainWeakening_CalcRatio( void* strainWeakening, void* particle ) {
	StrainWeakening*             self   = (StrainWeakening*) strainWeakening;
	StrainWeakening_ParticleExt* particleExt;
	double                       strainWeakeningRatio;

	if ( self == NULL ) 
		return 0.0;

	particleExt = ExtensionManager_Get( self->swarm->particleExtensionMgr, particle, self->particleExtHandle );

	if ( particleExt->postFailureWeakening < 0.0 ) 
		particleExt->postFailureWeakening = 0.0;

	strainWeakeningRatio = particleExt->postFailureWeakening / self->softeningStrain;

	if (strainWeakeningRatio > 1.0)
		strainWeakeningRatio = 1.0;
	
	return strainWeakeningRatio;
}

void StrainWeakening_AssignIncrement( 
		void*                            strainWeakening,
		ConstitutiveMatrix*              constitutiveMatrix,
		MaterialPointsSwarm*             swarm,
		Element_LocalIndex               lElement_I,
		MaterialPoint*                   particle,
		double                           yieldCriterion,
		double                           yieldIndicator ) 
{
	StrainWeakening*             self   = (StrainWeakening*) strainWeakening;
	StrainWeakening_ParticleExt* particleExt;

	particleExt = ExtensionManager_Get( swarm->particleExtensionMgr, particle, self->particleExtHandle );
	
	particleExt->postFailureWeakeningIncrement = 
		self->_calcIncrement( self, constitutiveMatrix, swarm, lElement_I, particle, yieldCriterion, yieldIndicator );
}

double StrainWeakening_GetPostFailureWeakening( void* strainWeakening, void* particle ) {
	StrainWeakening*             self   = (StrainWeakening*) strainWeakening;
	StrainWeakening_ParticleExt* particleExt;

	particleExt = ExtensionManager_Get( self->swarm->particleExtensionMgr, particle, self->particleExtHandle );

	return particleExt->postFailureWeakening;
}

double StrainWeakening_GetInitialDamageFraction( void* strainWeakening, void* particle ) {
	StrainWeakening*             self   = (StrainWeakening*) strainWeakening;

	return self->initialDamageFraction;
}

