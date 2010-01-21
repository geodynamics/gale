/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
** Copyright (c) 2005, Monash Cluster Computing
** All rights reserved.
** Redistribution and use in source and binary forms, with or without modification,
** are permitted provided that the following conditions are met:
**
**       * Redistributions of source code must retain the above copyright notice,
**          this list of conditions and the following disclaimer.
**       * Redistributions in binary form must reproduce the above copyright
**       notice, this list of conditions and the following disclaimer in the
**       documentation and/or other materials provided with the distribution.
**       * Neither the name of the Monash University nor the names of its contributors
**       may be used to endorse or promote products derived from this software
**       without specific prior written permission.
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
*%    Louis Moresi - Louis.Moresi@sci.monash.edu.au
*%
** Contributors:
*+    Robert Turnbull
*+    Vincent Lemiale
*+    Louis Moresi
*+    David May
*+    David Stegman
*+    Mirko Velic
*+    Patrick Sunter
*+    Julian Giordani
*+
** $Id: DruckerPrager.h 743 2008-06-23 01:49:43Z JulianGiordani $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#ifndef __Underworld_Rheology_DruckerPrager_h__
#define __Underworld_Rheology_DruckerPrager_h__

   /** Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
   extern const Type DruckerPrager_Type;

   typedef struct {
      float          brightness;
      float          opacity;
      float          diameter;
      Particle_Bool  tensileFailure;
   }  DruckerPrager_Particle;

   /** Rheology class contents - this is defined as a macro so that sub-classes of this class can use this macro at the start of the definition of their struct */
   #define __DruckerPrager \
      /* Parent info */ \
      __VonMises \
      /* Virtual functions go here */ \
      /* General Info */\
      SwarmVariable*                                      brightness;                            \
      SwarmVariable*                                      opacity;                               \
      SwarmVariable*                                      diameter;                              \
      SwarmVariable*                                      tensileFailure;                        \
      ExtensionInfo_Index                                 particleExtHandle;                     \
      /* Param passed in */\
      FeVariable*                                         pressureField;                         \
      SwarmVariable*                                      swarmPressure;                         \
      double                                              minimumYieldStress;                    \
      double                                              frictionCoefficient;                   \
      double                                              frictionCoefficientAfterSoftening;     \
		double                                              boundaryCohesion;                              \
		double                                              boundaryCohesionAfterSoftening;                \
		double                                              boundaryFrictionCoefficient;                   \
		double                                              boundaryFrictionCoefficientAfterSoftening;     \
                Bool                                                boundaryTop; \
                Bool                                                boundaryBottom; \
                Bool                                                boundaryLeft; \
                Bool                                                boundaryRight; \
                Bool                                                boundaryFront; \
                Bool                                                boundaryBack; \
		/* Stored values that are calculated once for each particle */ \
                double                                              trace; \
		TensorArray                                         currentVelocityGradient;               \
		double                                              currentPressure;                       \
		double                                              strainRateSecondInvariant;                   \
                HydrostaticTerm*                                    hydrostaticTerm; \
      double curFrictionCoef;


   struct DruckerPrager { __DruckerPrager };

   /* Public Constructor */
   DruckerPrager* DruckerPrager_New(
      Name                  name,
      AbstractContext*      context,
      StrainWeakening*      strainWeakening,
      MaterialPointsSwarm*  materialPointsSwarm,
      double                minVisc,
      FeVariable*           strainRateField,
      SwarmVariable*        swarmStrainRate,
      double                cohesion,
      double                cohesionAfterSoftening,
      Bool                  strainRateSoftening,
      FeVariable*           pressureField,
      SwarmVariable*        swarmPressure,
      double                minimumYieldStress,
      double                frictionCoefficient,
      double                frictionCoefficientAfterSoftening,
      double                                              boundaryCohesion,
      double                                              boundaryCohesionAfterSoftening,
      double                                              boundaryFrictionCoefficient,
      double                                              boundaryFrictionCoefficientAfterSoftening,
      Bool                                                boundaryTop,
      Bool                                                boundaryBottom,
      Bool                                                boundaryLeft,
      Bool                                                boundaryRight,
      Bool                                                boundaryFront,
      Bool                                                boundaryBack,
      double                                              trace,
      TensorArray                                         currentVelocityGradient,
      double                                              currentPressure,
      double                                              strainRateSecondInvariant,
      HydrostaticTerm*                                    hydrostaticTerm);

   /** Private Constructor: This will accept all the virtual functions for this class as arguments. */
	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define DRUCKERPRAGER_DEFARGS \
                VONMISES_DEFARGS

	#define DRUCKERPRAGER_PASSARGS \
                VONMISES_PASSARGS

   DruckerPrager* _DruckerPrager_New(  DRUCKERPRAGER_DEFARGS  ) ;

   /* 'Stg_Component' implementations */
   void* _DruckerPrager_DefaultNew( Name name ) ;
   void _DruckerPrager_AssignFromXML( void* rheology, Stg_ComponentFactory* cf, void* data );

   void _DruckerPrager_Build( void* rheology, void* data );
   void _DruckerPrager_Initialise( void* rheology, void* data ) ;
   void _DruckerPrager_Destroy( void* rheology, void* data ) ;
   void _DruckerPrager_Init(
         DruckerPrager*                                     self,
         FeVariable*                                        pressureField,
         SwarmVariable*                                     swarmPressure,
         MaterialPointsSwarm*                               materialPointsSwarm,
         double                                             minimumYieldStress,
         double                                             frictionCoefficient,
         double                                             frictionCoefficientAfterSoftening );


   /* 'YieldRheology' implementations */
   double _DruckerPrager_GetYieldCriterion(
         void*                            druckerPrager,
         ConstitutiveMatrix*              constitutiveMatrix,
         MaterialPointsSwarm*             materialPointsSwarm,
         Element_LocalIndex               lElement_I,
         MaterialPoint*                   materialPoint,
         Coord                            xi );

   void _DruckerPrager_HasYielded(
         void*                            rheology,
         ConstitutiveMatrix*              constitutiveMatrix,
         MaterialPointsSwarm*             materialPointsSwarm,
         Element_LocalIndex               lElement_I,
         MaterialPoint*                   materialPoint,
         double                           yieldCriterion,
         double                           yieldIndicator );

   void _DruckerPrager_UpdateDrawParameters( void* rheology ) ;

	double _DruckerPrager_EffectiveCohesion( void* rheology, void* materialPoint, Bool inside_boundary ) ;
	double _DruckerPrager_EffectiveFrictionCoefficient( void* rheology, void* materialPoint, Bool inside_boundary ) ;
	
#endif

