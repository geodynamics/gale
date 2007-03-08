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
** $Id: StrainWeakening.h 354 2006-10-12 08:19:27Z SteveQuenette $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#ifndef __Underworld_Rheology_StrainWeakening_h__
#define __Underworld_Rheology_StrainWeakening_h__

	/** Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
	extern const Type StrainWeakening_Type;

	typedef struct {
		double postFailureWeakening;
		double postFailureWeakeningIncrement;
	}  StrainWeakening_ParticleExt; 

	/* virtual function interface */
	typedef double (StrainWeakening_CalcIncrementFunction) ( 		
		void*                            strainWeakening,
		ConstitutiveMatrix*              constitutiveMatrix,
		MaterialPointsSwarm*             swarm,
		Element_LocalIndex               lElement_I,
		MaterialPoint*                   particle,
		double                           yieldCriterion,
		double                           yieldIndicator );
	
	/** Rheology class contents - this is defined as a macro so that sub-classes of this class can use this macro at the start of the definition of their struct */
	#define __StrainWeakening \
		/* Parent info */ \
 		__TimeIntegratee \
		/* Virtual functions go here */ \
		StrainWeakening_CalcIncrementFunction* _calcIncrement;                    \
		/* General Info */\
		SwarmVariable*                         postFailureWeakening;              \
		SwarmVariable*                         postFailureWeakeningIncrement;     \
		ExtensionInfo_Index                    particleExtHandle;                 \
		/* Param passed in */ \
		MaterialPointsSwarm*                   swarm;                             \
		double                                 healingRate;                       \
		double                                 softeningStrain;                   \
		double                                 initialDamageFraction;             \
		double                                 initialDamageWavenumber;           \
		double                                 initialDamageFactor;               \
		Stg_Shape*                             initialStrainShape;                \
		long int                               randomSeed;
				
	struct StrainWeakening { __StrainWeakening };

	/** Private Constructor: This will accept all the virtual functions for this class as arguments. */
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
		Name                                               name ) ;

	void _StrainWeakening_Init(
		StrainWeakening*                                   self,
		MaterialPointsSwarm*                               swarm,
		double                                             healingRate,
		double                                             softeningStrain,
		double                                             initialDamageFraction,
		double                                             initialDamageWavenumber,
		double                                             initialDamageFactor,
		long int                                           randomSeed,
		Stg_Shape*                                         initialStrainShape  );
	
	/* 'Stg_Component' implementations */
	void* _StrainWeakening_DefaultNew( Name name ) ;
	void _StrainWeakening_Construct( void* rheology, Stg_ComponentFactory* cf, void* data );
	void _StrainWeakening_Build( void* strainWeakening, void* data ) ;
	void _StrainWeakening_Initialise( void* strainWeakening, void* data ) ;
	
	Bool _StrainWeakening_TimeDerivative( void* _strainWeakening, Index lParticle_I, double* timeDeriv );

	/* Private Functions */
	void _StrainWeakening_MakeValuesPositive( void* timeIntegrator, void* strainWeakening ) ;

	double _StrainWeakening_CalcIncrementIsotropic( 
		void*                            strainWeakening,
		ConstitutiveMatrix*              constitutiveMatrix,
		MaterialPointsSwarm*             swarm,
		Element_LocalIndex               lElement_I,
		MaterialPoint*                   particle,
		double                           yieldCriterion,
		double                           yieldIndicator );

	/* Public Functions */
	double StrainWeakening_CalcRatio( void* strainWeakening, void* particle ) ;

	void StrainWeakening_AssignIncrement( 
		void*                            strainWeakening,
		ConstitutiveMatrix*              constitutiveMatrix,
		MaterialPointsSwarm*             swarm,
		Element_LocalIndex               lElement_I,
		MaterialPoint*                   particle,
		double                           yieldCriterion,
		double                           yieldIndicator ) ;

	double StrainWeakening_GetPostFailureWeakening( void* strainWeakening, void* particle ) ;
	double StrainWeakening_GetInitialDamageFraction( void* strainWeakening, void* particle );

#endif
