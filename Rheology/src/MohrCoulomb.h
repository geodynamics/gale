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
** $Id: MohrCoulomb.h 11823 2008-04-17 17:17:38Z walter $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#ifndef __Underworld_Rheology_MohrCoulomb_h__
#define __Underworld_Rheology_MohrCoulomb_h__

	/** Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
	extern const Type MohrCoulomb_Type;
	
	/** Rheology class contents - this is defined as a macro so that sub-classes of this class can use this macro at the start of the definition of their struct */
	#define __MohrCoulomb \
		/* Parent info */ \
 		__YieldRheology \
		/* Virtual functions go here */ \
		/* Material Parameters */\
		ExtensionInfo_Index                                 particleExtHandle;                     \
		/* Param passed in */ \
		MaterialPointsSwarm*                                materialPointsSwarm;                                 \
		FeVariable*                                         pressureField;                         \
		FeVariable*                                         strainRateField;\
		SwarmVariable*                                      swarmStrainRate;\
		/* Director component is used to update the normal */\
		Director*                                           director;                              \
		double                                              cohesion;                              \
		double                                              cohesionAfterSoftening;                \
		double                                              frictionCoefficient;                   \
		double                                              frictionCoefficientAfterSoftening;     \
		double                                              minimumYieldStress;                    \
		/* Stored values that are calculated once for each particle */ \
		Eigenvector                                         currentEigenvectorList[3];             \
		SymmetricTensor                                     currentStrainRate; \
		SymmetricTensor                                     currentStress;                         \
		double                                              currentPressure;                       \
		double                                              strainRateSecondInvariant;

	
	struct MohrCoulomb { __MohrCoulomb };
	
	/** Private Constructor: This will accept all the virtual functions for this class as arguments. */
	MohrCoulomb* _MohrCoulomb_New( 
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
		Name                                               name ) ;
	
	/* 'Stg_Component' implementations */
	void* _MohrCoulomb_DefaultNew( Name name ) ;
	void _MohrCoulomb_AssignFromXML( void* rheology, Stg_ComponentFactory* cf,
                                     void *data );
   void _MohrCoulomb_Init(    
		MohrCoulomb*          self,
		FeVariable*           pressureField,
		FeVariable*           strainRateField,
		SwarmVariable*        swarmStrainRate,
		MaterialPointsSwarm*  materialPointsSwarm,
		double                cohesion,
		double                cohesionAfterSoftening,
		double                frictionCoefficient,
		double                frictionCoefficientAfterSoftening,
		double                minimumYieldStress);
	void _MohrCoulomb_Build( void* rheology, void* data );
	void _MohrCoulomb_Initialise( void* rheology, void* data );
   void _MohrCoulomb_Destroy( void* rheology, void* data );
	/* 'YieldRheology' implementations */
	void _MohrCoulomb_ModifyConstitutiveMatrix( 
		void*                                              rheology, 
		ConstitutiveMatrix*                                constitutiveMatrix,
		MaterialPointsSwarm*                               materialPointsSwarm,
		Element_LocalIndex                                 lElement_I,
		MaterialPoint*                                     materialPoint,
		Coord                                              xi );
		
	double _MohrCoulomb_GetYieldCriterion( 
		void*                                              rheology,
		ConstitutiveMatrix*                                constitutiveMatrix,
		MaterialPointsSwarm*                               materialPointsSwarm,
		Element_LocalIndex                                 lElement_I,
		MaterialPoint*                                     materialPoint,
		Coord                                              xi );
	
	double _MohrCoulomb_GetYieldIndicator( 
		void*                                              rheology,
		ConstitutiveMatrix*                                constitutiveMatrix,
		MaterialPointsSwarm*                               materialPointsSwarm,
		Element_LocalIndex                                 lElement_I,
		MaterialPoint*                                     materialPoint,
		Coord                                              xi );
	
	void _MohrCoulomb_HasYielded( 
		void*                                              rheology,
		ConstitutiveMatrix*                                constitutiveMatrix,
		MaterialPointsSwarm*                               materialPointsSwarm,
		Element_LocalIndex                                 lElement_I,
		MaterialPoint*                                     materialPoint,
		double                                             yieldCriterion,
		double                                             yieldIndicator );
	
	double _MohrCoulomb_EffectiveCohesion( void* rheology, void* materialPoint ) ;
	double _MohrCoulomb_EffectiveFrictionCoefficient( void* rheology, void* materialPoint ) ;
	double _MohrCoulomb_Sigma_nn( void* rheology, void* materialPoint, Dimension_Index dim ) ;

	void _MohrCoulomb_StoreCurrentParameters( 
		void*                                              rheology,
		ConstitutiveMatrix*                                constitutiveMatrix, 
		MaterialPointsSwarm*                               materialPointsSwarm,
		Element_LocalIndex                                 lElement_I, 
		MaterialPoint*                                     materialPoint,
		Coord                                              xi ) ;
	
#endif
