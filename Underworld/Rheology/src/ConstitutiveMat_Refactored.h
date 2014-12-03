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
** $Id: ConstitutiveMat_Refactored.h 803 2008-09-11 05:22:20Z LukeHodkinson $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#ifndef __Underworld_Rheology_ConstitutiveMat_Refactored_h__
#define __Underworld_Rheology_ConstitutiveMat_Refactored_h__
	
	/* typedefs for virtual functions: */
	typedef double (ConstitutiveMat_Refactored_GetValueFunc) ( void* constitutiveMatrix );
	typedef void (ConstitutiveMat_Refactored_SetValueFunc) ( void* constitutiveMatrix, double value );
	typedef void (ConstitutiveMat_Refactored_SetSecondViscosityFunc) ( void* constitutiveMatrix, double value, XYZ vector );
	typedef void (ConstitutiveMat_Refactored_Assemble_D_B_Func) ( void* constitutiveMatrix, double** GNx, Node_Index node_I, double** D_B );
	typedef void (ConstitutiveMat_Refactored_CalculateStressFunc) ( void* constitutiveMatrix, SymmetricTensor strainRate, SymmetricTensor stress );
	
	extern const Type ConstitutiveMat_Refactored_Type;
	
	/* ConstitutiveMat_Refactored information */
	#define __ConstitutiveMat_Refactored  \
		/* Parent info */ \
		__Stg_Component \
		PICelleratorContext*											context; \
		\
		/* Virtual functions go here */ \
		ConstitutiveMat_Refactored_SetValueFunc*				_setValue; \
		ConstitutiveMat_Refactored_GetValueFunc*				_getViscosity; \
		ConstitutiveMat_Refactored_SetValueFunc*				_isotropicCorrection; \
		ConstitutiveMat_Refactored_SetSecondViscosityFunc*	_setSecondViscosity; \
		ConstitutiveMat_Refactored_Assemble_D_B_Func*		_assemble_D_B; \
		ConstitutiveMat_Refactored_CalculateStressFunc*		_calculateStress;\
		\
		/* ConstitutiveMat_Refactored info */ \
		Stream*															debug; \
		double**															matrixData; \
		double															derivs[9]; \
		Dimension_Index												dim; \
		Materials_Register*											materials_Register; \
		Bool																isDiagonal; \
		Index																columnSize; \
		Index																rowSize; \
		Bool																previousSolutionExists; \
		SystemLinearEquations*										sle; \
		Iteration_Index												sleNonLinearIteration_I;
		
	struct ConstitutiveMat_Refactored { __ConstitutiveMat_Refactored };



	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define CONSTITUTIVEMAT_REFACTORED_DEFARGS \
                STG_COMPONENT_DEFARGS, \
                ConstitutiveMat_Refactored_SetValueFunc*                       _setValue, \
                ConstitutiveMat_Refactored_GetValueFunc*                   _getViscosity, \
                ConstitutiveMat_Refactored_SetValueFunc*            _isotropicCorrection, \
                ConstitutiveMat_Refactored_SetSecondViscosityFunc*   _setSecondViscosity, \
                ConstitutiveMat_Refactored_Assemble_D_B_Func*              _assemble_D_B, \
                ConstitutiveMat_Refactored_CalculateStressFunc*         _calculateStress

	#define CONSTITUTIVEMAT_REFACTORED_PASSARGS \
                STG_COMPONENT_PASSARGS, \
	        _setValue,            \
	        _getViscosity,        \
	        _isotropicCorrection, \
	        _setSecondViscosity,  \
	        _assemble_D_B,        \
	        _calculateStress    

	ConstitutiveMat_Refactored* _ConstitutiveMat_Refactored_New(  CONSTITUTIVEMAT_REFACTORED_DEFARGS  );

	void _ConstitutiveMat_Refactored_Init(                                                                                
		void*                constitutiveMatrix,                                                                           
		Dimension_Index      dim,                                                                                          
		PICelleratorContext* context,                                                                                      
		Materials_Register*  materials_Register );
	
	/* 'Stg_Class' Virtual Functions */
	void _ConstitutiveMat_Refactored_Delete( void* constitutiveMatrix );

	void _ConstitutiveMat_Refactored_Print( void* constitutiveMatrix, Stream* stream );

	#define ConstitutiveMat_Refactored_Copy( self ) \
		(ConstitutiveMat_Refactored*)Stg_Class_Copy( self, NULL, False, NULL, NULL )
	#define ConstitutiveMat_Refactored_DeepCopy( self ) \
		(ConstitutiveMat_Refactored*)Stg_Class_Copy( self, NULL, True, NULL, NULL )

	void* _ConstitutiveMat_Refactored_Copy( const void* constitutiveMatrix, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );
	
	/* 'Stg_Component' Virtual Functions */
	void _ConstitutiveMat_Refactored_AssignFromXML( void* constitutiveMatrix, Stg_ComponentFactory* cf, void* data );

	void _ConstitutiveMat_Refactored_Build( void* constitutiveMatrix, void* data );

	void _ConstitutiveMat_Refactored_Initialise( void* constitutiveMatrix, void* data );

	void _ConstitutiveMat_Refactored_Execute( void* constitutiveMatrix, void* data );

	void _ConstitutiveMat_Refactored_Destroy( void* constitutiveMatrix, void* data );

	/* Wrapper macros to virtual functions - These must be macros for the sake of speed */
	#define ConstitutiveMat_Refactored_SetValueInAllEntries( constitutiveMatrix, value ) \
		(((ConstitutiveMat_Refactored*) constitutiveMatrix)->_setValue( constitutiveMatrix, value ))

	#define ConstitutiveMat_Refactored_SetSecondViscosity( constitutiveMatrix, deltaViscosity, director ) \
		(((ConstitutiveMat_Refactored*) constitutiveMatrix)->_setSecondViscosity( constitutiveMatrix, deltaViscosity, director ))

	#define ConstitutiveMat_Refactored_GetIsotropicViscosity( constitutiveMatrix ) \
		(((ConstitutiveMat_Refactored*) constitutiveMatrix)->_getViscosity( constitutiveMatrix ))

	#define ConstitutiveMat_Refactored_IsotropicCorrection( constitutiveMatrix, isotropicCorrection ) \
		(((ConstitutiveMat_Refactored*) constitutiveMatrix)->_isotropicCorrection( constitutiveMatrix, isotropicCorrection ))

	#define ConstitutiveMat_Refactored_Assemble_D_B( constitutiveMatrix, GNx, node_I, D_B ) \
		(((ConstitutiveMat_Refactored*) constitutiveMatrix)->_assemble_D_B( constitutiveMatrix, GNx, node_I, D_B ))

	#define ConstitutiveMat_Refactored_CalculateStress( constitutiveMatrix, strainRate, stress ) \
		(((ConstitutiveMat_Refactored*) constitutiveMatrix)->_calculateStress( constitutiveMatrix, strainRate, stress ) )

	/* +++ Public Functions +++ */
	void ConstitutiveMat_Refactored_MultiplyByValue( void* constitutiveMatrix, double factor );

	void ConstitutiveMat_Refactored_PrintContents( void* constitutiveMatrix, Stream* stream );
	
	void ConstitutiveMat_Refactored_ZeroMatrix( void* constitutiveMatrix );

	void ConstitutiveMat_Refactored_SetIsotropicViscosity( void* constitutiveMatrix, double viscosity );

	void ConstitutiveMat_Refactored_Assemble( 
		void*							constitutiveMatrix,
		Element_LocalIndex		lElement_I,
		IntegrationPointsSwarm*	swarm,
		IntegrationPoint*			particle );

	#define ConstitutiveMat_Refactored_SetToNonLinear( constitutiveMatrix ) \
		StiffnessMatrix_SetToNonLinear( constitutiveMatrix->stiffnessMatrix )

	#define ConstitutiveMat_Refactored_GetMesh( constitutiveMatrix ) \
		( (constitutiveMatrix)->stiffnessMatrix->rowVariable->feMesh )

#endif /* __Underworld_Rheology_ConstitutiveMat_Refactored_h__ */

