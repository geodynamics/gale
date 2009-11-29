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
** $Id: BuiterStrainWeakening.h 98 2005-12-15 12:56:32Z VincentLemiale $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/* Modified 2006 Walter Landry to implement Buiter Strain Weakening */


#ifndef __Underworld_Rheology_BuiterStrainWeakening_h__
#define __Underworld_Rheology_BuiterStrainWeakening_h__

	/** Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
	extern const Type BuiterStrainWeakening_Type;

	/** Rheology class contents - this is defined as a macro so that sub-classes of this class can use this macro at the start of the definition of their struct */
	#define __BuiterStrainWeakening \
		/* Parent info */ \
 		__StrainWeakening \
		/* Virtual functions go here */ \
						\
		/* General Info */ \
				
	struct BuiterStrainWeakening { __BuiterStrainWeakening };

	/** Public Constructor */
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
		Stg_Shape*                                         initialStrainShape  );

	/** Private Constructor: This will accept all the virtual functions for this class as arguments. */
	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define BUITERSTRAINWEAKENING_DEFARGS \
                STRAINWEAKENING_DEFARGS

	#define BUITERSTRAINWEAKENING_PASSARGS \
                STRAINWEAKENING_PASSARGS

	BuiterStrainWeakening* _BuiterStrainWeakening_New(  BUITERSTRAINWEAKENING_DEFARGS  ) ;
	
	/* 'Stg_Component' implementations */
	void* _BuiterStrainWeakening_DefaultNew( Name name ) ;
	void _BuiterStrainWeakening_AssignFromXML( void* rheology, Stg_ComponentFactory* cf, void* data );
	void _BuiterStrainWeakening_Build( void* strainWeakening, void* data ) ;
	void _BuiterStrainWeakening_Initialise( void* strainWeakening, void* data ) ;
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
		Stg_Shape*                                         initialStrainShape  );

#endif

