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
** $Id: YieldRheology.h 779 2008-08-06 15:50:41Z LukeHodkinson $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#ifndef __Underworld_Rheology_YieldRheology_h__
#define __Underworld_Rheology_YieldRheology_h__

	/** Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
	extern const Type YieldRheology_Type;
	
	/* virtual function interface */
	typedef double (YieldRheology_GetYieldCriterionFunction) ( 		
		void*                            yieldRheology,
		ConstitutiveMatrix*              constitutiveMatrix,
		MaterialPointsSwarm*             materialPointsSwarm,
		Element_LocalIndex               lElement_I,
		MaterialPoint*                   materialPoint,
		Coord                            xi );

	typedef double (YieldRheology_GetYieldIndicatorFunction) (
		void*                            yieldRheology,
		ConstitutiveMatrix*              constitutiveMatrix,
		MaterialPointsSwarm*             materialPointsSwarm,
		Element_LocalIndex               lElement_I,
		MaterialPoint*                   materialPoint,
		Coord                            xi );

	typedef void   (YieldRheology_HasYieldedFunction)        ( 
		void*                            yieldRheology,
		ConstitutiveMatrix*              constitutiveMatrix,
		MaterialPointsSwarm*             materialPointsSwarm,
		Element_LocalIndex               lElement_I,
		MaterialPoint*                   materialPoint,
		double                           yieldCriterion,
		double                           yieldIndicator );
		
	/** YieldRheology class contents - this is defined as a macro so that sub-classes of this class can use this macro at the start of the definition of their struct */
	#define __YieldRheology \
		/* Macro defining parent goes here - This means you can cast this class as its parent */ \
		__Rheology \
		/* Virtual functions go here */ \
		YieldRheology_GetYieldCriterionFunction*           _getYieldCriterion;          \
		YieldRheology_GetYieldIndicatorFunction*           _getYieldIndicator;          \
		YieldRheology_HasYieldedFunction*                  _hasYielded;                 \
		/* Other info */ \
		StrainWeakening*                                   strainWeakening;             \
		ExtensionInfo_Index                                hasYieldedParticleExtHandle; \
		SwarmVariable*                                     hasYieldedVariable;          \
                double                                             yieldCriterion;              \
                double                                             yieldIndicator;              \
		double                                             minVisc;

	struct YieldRheology { __YieldRheology };

	/** Private Constructor: This will accept all the virtual functions for this class as arguments. */
	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define YIELDRHEOLOGY_DEFARGS \
                RHEOLOGY_DEFARGS, \
                YieldRheology_GetYieldCriterionFunction*  _getYieldCriterion, \
                YieldRheology_GetYieldIndicatorFunction*  _getYieldIndicator, \
                YieldRheology_HasYieldedFunction*                _hasYielded

	#define YIELDRHEOLOGY_PASSARGS \
                RHEOLOGY_PASSARGS, \
	        _getYieldCriterion, \
	        _getYieldIndicator, \
	        _hasYielded       

	YieldRheology* _YieldRheology_New(  YIELDRHEOLOGY_DEFARGS  );

	/* 'Stg_Class' implementations */
	void _YieldRheology_Delete( void* rheology );
	void _YieldRheology_Print( void* rheology, Stream* stream );
	#define YieldRheology_Copy( self ) \
		(YieldRheology*)Stg_Class_Copy( self, NULL, False, NULL, NULL )
	#define YieldRheology_DeepCopy( self ) \
		(YieldRheology*)Stg_Class_Copy( self, NULL, True, NULL, NULL )
	void* _YieldRheology_Copy( void* rheology, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );
	
	/* 'Stg_Component' implementations */
	void* _YieldRheology_DefaultNew( Name name ) ;
	void _YieldRheology_AssignFromXML( void* rheology, Stg_ComponentFactory* cf, void* data );
	void _YieldRheology_Build( void* rheology, void* data );
	void _YieldRheology_Initialise( void* rheology, void* data );
	void _YieldRheology_Execute( void* rheology, void* data );
	void _YieldRheology_Destroy( void* rheology, void* data );
	void _YieldRheology_Init( YieldRheology* self, StrainWeakening* strainWeakening, MaterialPointsSwarm* materialPointsSwarm, double minVisc );
	
	void _YieldRheology_ModifyConstitutiveMatrix( 
		void*                                              rheology, 
		ConstitutiveMatrix*                                constitutiveMatrix,
		MaterialPointsSwarm*                               materialPointsSwarm,
		Element_LocalIndex                                 lElement_I,
		MaterialPoint*                                     materialPoint,
		Coord                                              xi );

	/** Public Functions - Defining these as macros - for speed (compared to a wrapper function)  
	 *  The reason why these functions are used rather than using self->_name_of_the_virtual_function
	 *  directly is for consistency wrt the rest of the code. It's done that way everywhere. */

	#define YieldRheology_GetYieldCriterion( yieldRheology, constitutiveMatrix, materialPointsSwarm, lElement_I, materialPoint, xi ) \
		( ( (YieldRheology*)yieldRheology )->_getYieldCriterion( \
			yieldRheology, constitutiveMatrix, materialPointsSwarm, lElement_I, materialPoint, xi ) )

	#define YieldRheology_GetYieldIndicator( yieldRheology, constitutiveMatrix, materialPointsSwarm, lElement_I, materialPoint, xi ) \
		( ( (YieldRheology*)yieldRheology )->_getYieldIndicator( \
			yieldRheology, constitutiveMatrix, materialPointsSwarm, lElement_I, materialPoint, xi ) )

	#define YieldRheology_HasYielded( yieldRheology, constitutiveMatrix, materialPointsSwarm, lElement_I, materialPoint, yieldCriterion, yieldIndicator ) \
		( ( (YieldRheology*)yieldRheology )->_hasYielded( \
			yieldRheology, constitutiveMatrix, materialPointsSwarm, lElement_I, materialPoint, yieldCriterion, yieldIndicator ) )

        /* This macro sets a bool to TRUE or FLAG depending on whether a materialPoint has yielded or not
	 * If hasYieldedParticleExtHandle has a value of '-1', then it means that there are no material points
	 * Otherwise the boolean pointed by the returned pointer of ExtensionManager_Get function is set to 'flag' value*/ 
	#define YieldRheology_SetParticleFlag( self, materialPointsSwarm, materialPoint, flag )\
		if ( (self)->hasYieldedParticleExtHandle != (ExtensionInfo_Index) -1 ) \
		* (Particle_Bool*) ExtensionManager_Get( (materialPointsSwarm)->particleExtensionMgr, (materialPoint), (self)->hasYieldedParticleExtHandle ) \
			= (Particle_Bool) (flag)
	
#endif

