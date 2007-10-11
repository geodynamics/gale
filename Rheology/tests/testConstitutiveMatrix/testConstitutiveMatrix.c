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
** $Id: testConstitutiveMatrix.c 610 2007-10-11 08:09:29Z SteveQuenette $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>

#include <Underworld/Rheology/Rheology.h>

const Type Underworld_testConstitutiveMatrix_Type = "Underworld_testConstitutiveMatrix";
typedef struct {
	__Codelet
} Underworld_testConstitutiveMatrix;

void SetMatrixWithSecondViscosity2( ConstitutiveMatrix* constitutiveMatrix, Element_LocalIndex lElement_I, Particle_InCellIndex cParticle_I, void* data ) {
	XYZ director = { 1.0 , 2.0 , 3.0 };

	StGermain_VectorNormalise( director, constitutiveMatrix->dim );

	ConstitutiveMatrix_SetSecondViscosity( constitutiveMatrix, 2.0, director );
}

void testConstitutiveMatrix( FiniteElementContext* context ) {
	Stream*                     stream     = Journal_Register( Info_Type, "testConstitutiveMatrix" );
	SymmetricTensor             stress;
	SymmetricTensor             strainRate = {1, 2, 3, 4, 5, 6};
	ConstitutiveMatrix*         constitutiveMatrix;

	constitutiveMatrix = (ConstitutiveMatrix*) LiveComponentRegister_Get( context->CF->LCRegister, "constitutiveMatrix" );

	/* Create Constitutive Matrix */
	Stream_RedirectFile_WithPrependedPath( stream, context->outputPath, "output.dat" );
	
	Journal_Printf( stream, "Constitutive Matrix after initialisation.\n" );
	Stg_Class_Print( constitutiveMatrix, stream );
	
	Journal_Printf( stream, "Constitutive Matrix after func 'SetMatrixWithThrees'.\n" );
	ConstitutiveMatrix_SetValueInAllEntries( constitutiveMatrix, 3.0 );
	ConstitutiveMatrix_PrintContents( constitutiveMatrix, stream );
	Journal_PrintBool( stream, constitutiveMatrix->isDiagonal );
	
	ConstitutiveMatrix_CalculateStress( constitutiveMatrix, strainRate, stress );
	Journal_PrintSymmetricTensor( stream, strainRate, constitutiveMatrix->dim );
	Journal_PrintSymmetricTensor( stream, stress, constitutiveMatrix->dim );
	
	Journal_Printf( stream, "Constitutive Matrix after func 'SetMatrixWithViscosity7'.\n" );
	ConstitutiveMatrix_SetIsotropicViscosity( constitutiveMatrix, 7.0 );
	ConstitutiveMatrix_PrintContents( constitutiveMatrix, stream );
	Journal_PrintBool( stream, constitutiveMatrix->isDiagonal );

	Journal_Printf( stream, "Viscosity = %4g\n", ConstitutiveMatrix_GetIsotropicViscosity( constitutiveMatrix ) );
	ConstitutiveMatrix_CalculateStress( constitutiveMatrix, strainRate, stress );
	Journal_PrintSymmetricTensor( stream, strainRate, constitutiveMatrix->dim );
	Journal_PrintSymmetricTensor( stream, stress, constitutiveMatrix->dim );
	constitutiveMatrix->isDiagonal = False;
	ConstitutiveMatrix_CalculateStress( constitutiveMatrix, strainRate, stress );
	Journal_PrintSymmetricTensor( stream, stress, constitutiveMatrix->dim );
	
	Journal_Printf( stream, "Constitutive Matrix after func 'SetMatrixWithSecondViscosity2'.\n" );
	SetMatrixWithSecondViscosity2( constitutiveMatrix, 0, 0, NULL );
	ConstitutiveMatrix_PrintContents( constitutiveMatrix, stream );
	Journal_PrintBool( stream, constitutiveMatrix->isDiagonal );
	
	ConstitutiveMatrix_CalculateStress( constitutiveMatrix, strainRate, stress );
	Journal_PrintSymmetricTensor( stream, strainRate, constitutiveMatrix->dim );
	Journal_PrintSymmetricTensor( stream, stress, constitutiveMatrix->dim );
}

void _Underworld_testConstitutiveMatrix_Construct( void* component, Stg_ComponentFactory* cf, void* data ) {
	FiniteElementContext* context;

	context = (FiniteElementContext*)Stg_ComponentFactory_ConstructByName( cf, "context", FiniteElementContext, True, data ) ;

	ContextEP_ReplaceAll( context, AbstractContext_EP_Execute, testConstitutiveMatrix );
}

void* _Underworld_testConstitutiveMatrix_DefaultNew( Name name ) {
	return Codelet_New(
			Underworld_testConstitutiveMatrix_Type,
			_Underworld_testConstitutiveMatrix_DefaultNew,
			_Underworld_testConstitutiveMatrix_Construct,
			_Codelet_Build,
			_Codelet_Initialise,
			_Codelet_Execute,
			_Codelet_Destroy,
			name );
}

Index Underworld_testConstitutiveMatrix_Register( PluginsManager* pluginsManager ) {
	return PluginsManager_Submit( pluginsManager, Underworld_testConstitutiveMatrix_Type, "0",
	_Underworld_testConstitutiveMatrix_DefaultNew );
}
