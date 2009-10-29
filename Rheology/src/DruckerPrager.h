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
		double                                              frictionCoefficientAfterSoftening; \
                                                                        \
                double curFrictionCoef;

				
	struct DruckerPrager { __DruckerPrager };
	
	/** Private Constructor: This will accept all the virtual functions for this class as arguments. */
	DruckerPrager* _DruckerPrager_New( 
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
	void* _DruckerPrager_DefaultNew( Name name ) ;
	void _DruckerPrager_AssignFromXML( void* rheology, Stg_ComponentFactory* cf, void* data );

	void _DruckerPrager_Build( void* rheology, void* data );
	void _DruckerPrager_Initialise( void* rheology, void* data ) ;
	
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
	
#endif
