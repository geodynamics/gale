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
** $Id: ConstitutiveMatCartesian_Refactored.h 430 2007-02-07 00:10:36Z PatrickSunter $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __Underworld_Rheology_ConstitutiveMatCartesian_Refactored_h__
#define __Underworld_Rheology_ConstitutiveMatCartesian_Refactored_h__

	/** Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
	extern const Type ConstitutiveMatCartesian_Refactored_Type;

	/** ConstitutiveMatCartesian_Refactored class contents - this is defined as a macro so that sub-classes of this class can use this macro at the start of the definition of their struct */
	#define __ConstitutiveMatCartesian_Refactored \
		/* Macro defining parent goes here - This means you can cast this class as its parent */ \
		__ConstitutiveMat_Refactored \
		\
		/* Virtual functions go here */ \
		\
		/* ConstitutiveMatCartesian_Refactored info */ \
                int maxNumNodes;                       \
                double** GNx;                          \
                double* Ni;

	struct ConstitutiveMatCartesian_Refactored { __ConstitutiveMatCartesian_Refactored };

	ConstitutiveMatCartesian_Refactored* ConstitutiveMatCartesian_Refactored_New( 
		Name                                                name,
		Dimension_Index                                     dim,
		FiniteElementContext*                               context,	
		Materials_Register*                                 materials_Register );

	ConstitutiveMatCartesian_Refactored* _ConstitutiveMatCartesian_Refactored_New( 
		SizeT                                               sizeOfSelf,  
		Type                                                type,
		Stg_Class_DeleteFunction*                           _delete,
		Stg_Class_PrintFunction*                            _print,
		Stg_Class_CopyFunction*                             _copy, 
		Stg_Component_DefaultConstructorFunction*           _defaultConstructor,
		Stg_Component_ConstructFunction*                    _construct,
		Stg_Component_BuildFunction*                        _build,
		Stg_Component_InitialiseFunction*                   _initialise,
		Stg_Component_ExecuteFunction*                      _execute,
		Stg_Component_DestroyFunction*                      _destroy,
		ConstitutiveMat_Refactored_SetValueFunc*            _setValue,
		ConstitutiveMat_Refactored_GetValueFunc*            _getViscosity,
		ConstitutiveMat_Refactored_SetValueFunc*            _isotropicCorrection,
		ConstitutiveMat_Refactored_SetSecondViscosityFunc*  _setSecondViscosity,
		ConstitutiveMat_Refactored_Assemble_D_B_Func*       _assemble_D_B,
		ConstitutiveMat_Refactored_CalculateStressFunc*     _calculateStress,
		Name                                                name );
	
	void ConstitutiveMatCartesian_Refactored_InitAll( 
		void*                                               constitutiveMatrix,
		Dimension_Index                                     dim,
		FiniteElementContext*                               context,
		Materials_Register*                                 materials_Register );

	void _ConstitutiveMatCartesian_Refactored_Delete( void* constitutiveMatrix );
	void _ConstitutiveMatCartesian_Refactored_Print( void* constitutiveMatrix, Stream* stream );

	void* _ConstitutiveMatCartesian_Refactored_DefaultNew( Name name ) ;
	void _ConstitutiveMatCartesian_Refactored_Construct( void* constitutiveMatrix, Stg_ComponentFactory* cf, void* data ) ;
	void _ConstitutiveMatCartesian_Refactored_Build( void* constitutiveMatrix, void* data ) ;
	void _ConstitutiveMatCartesian_Refactored_Initialise( void* constitutiveMatrix, void* data ) ;
	void _ConstitutiveMatCartesian_Refactored_Execute( void* constitutiveMatrix, void* data ) ;
	void _ConstitutiveMatCartesian_Refactored_Destroy( void* constitutiveMatrix, void* data ) ;

	void _ConstitutiveMatCartesian_Refactored_SetValueInAllEntries( void* constitutiveMatrix, double value ) ;
	void _ConstitutiveMatCartesian_Refactored2D_SetValueInAllEntries( void* constitutiveMatrix, double value ) ;
	void _ConstitutiveMatCartesian_Refactored3D_SetValueInAllEntries( void* constitutiveMatrix, double value ) ;

	double _ConstitutiveMatCartesian_Refactored_GetIsotropicViscosity( void* constitutiveMatrix ) ;
	double _ConstitutiveMatCartesian_Refactored2D_GetIsotropicViscosity( void* constitutiveMatrix ) ;
	double _ConstitutiveMatCartesian_Refactored3D_GetIsotropicViscosity( void* constitutiveMatrix ) ;

	void _ConstitutiveMatCartesian_Refactored_IsotropicCorrection( void* constitutiveMatrix, double isotropicCorrection ) ;
	void _ConstitutiveMatCartesian_Refactored2D_IsotropicCorrection( void* constitutiveMatrix, double isotropicCorrection ) ;
	void _ConstitutiveMatCartesian_Refactored3D_IsotropicCorrection( void* constitutiveMatrix, double isotropicCorrection ) ;

	void _ConstitutiveMatCartesian_Refactored_SetSecondViscosity( void* constitutiveMatrix, double deltaViscosity, XYZ director ) ;
	void _ConstitutiveMatCartesian_Refactored2D_SetSecondViscosity( void* constitutiveMatrix, double deltaViscosity, XYZ director );
	void _ConstitutiveMatCartesian_Refactored3D_SetSecondViscosity( void* constitutiveMatrix, double deltaViscosity, XYZ director );

	void _ConstitutiveMatCartesian_Refactored_Assemble_D_B( void* constitutiveMatrix, double** GNx, Node_Index node_I, double** D_B ) ;
	void _ConstitutiveMatCartesian_Refactored2D_Assemble_D_B( void* constitutiveMatrix, double** GNx, Node_Index node_I, double** D_B );
	void _ConstitutiveMatCartesian_Refactored3D_Assemble_D_B( void* constitutiveMatrix, double** GNx, Node_Index node_I, double** D_B );

	void _ConstitutiveMatCartesian_Refactored_CalculateStress( void* constitutiveMatrix, SymmetricTensor strainRate, SymmetricTensor stress ) ;
	void _ConstitutiveMatCartesian_Refactored2D_CalculateStress( void* constitutiveMatrix, SymmetricTensor strainRate, SymmetricTensor stress ) ;
	void _ConstitutiveMatCartesian_Refactored3D_CalculateStress( void* constitutiveMatrix, SymmetricTensor strainRate, SymmetricTensor stress ) ;

#endif
