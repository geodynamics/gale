#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pcu/pcu.h"
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include "Underworld/Underworld.h"

/* silly stgermain, I must define this */
#define CURR_MODULE_NAME "UnderworldContext.c"

typedef struct {
	UnderworldContext* context;
} ConstitutiveMatrixSuiteData;

void ConstitutiveMatrixSuite_Setup( ConstitutiveMatrixSuiteData* data ) { 
}

void ConstitutiveMatrixSuite_Teardown( ConstitutiveMatrixSuiteData* data ) {
	UnderworldContext* context = data->context;
	char rFile[PCU_PATH_MAX];
	int err;
	/* Now clean output path */
	sprintf(rFile, "%s/output.dat", context->outputPath );
	err = remove( rFile );
	if( err == -1 ) printf("Error in %s, can't delete the input.xml\n", __func__);
}

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


void ConstitutiveMatrixSuite_CartesianMatrix2D( ConstitutiveMatrixSuiteData* data ) {
	StiffnessMatrix*			stiffnessMatrix;
	SystemLinearEquations*	sle;
	Dictionary*					dictionary;
	UnderworldContext*		context;
	Stg_ComponentFactory*	cf;
	PetscViewer					bviewer;
	PetscReal					matrixNorm, errorNorm, test;
	Mat							expected;
	char							expected_file[PCU_PATH_MAX], output_file[PCU_PATH_MAX];
	char							*filename;
	double						tolerance;
	char							xml_input[PCU_PATH_MAX];

	pcu_docstring( "This test creates a 2D cartesian constitutive matrix data structure and tests its functionality against an expected file" );

	/* read in the xml input file */
	pcu_filename_input( "testConstitutiveMatrix_Cartesian2D.xml", xml_input );
	cf = stgMainInitFromXML( xml_input, MPI_COMM_WORLD, NULL );
	context = (UnderworldContext*)LiveComponentRegister_Get( cf->LCRegister, "context" );
	data->context = context;
	dictionary = context->dictionary;

	/* replace the Execute EP with the function to test the ConstitutiveMatrix */
	ContextEP_ReplaceAll( context, AbstractContext_EP_Execute, testConstitutiveMatrix );

	/* build and initialise all components */
	stgMainBuildAndInitialise( cf );

	/* now run the all EPs */
	stgMainLoop( cf );

	/* create the expected file string */
	pcu_filename_expected( "testConstitutiveMatrix_Cartesian2D.expected", expected_file );

	/* create the string for the output.dat */
	strcpy( output_file, context->outputPath );
	strcat( output_file, "/output.dat" );

	/* check the output file against the expected output file */
	pcu_check_fileEq( output_file, expected_file );

	/* end the run */
	stgMainDestroy( cf );
}

void ConstitutiveMatrixSuite( pcu_suite_t* suite ) {
   pcu_suite_setData( suite, ConstitutiveMatrixSuiteData );
   pcu_suite_setFixtures( suite, ConstitutiveMatrixSuite_Setup, ConstitutiveMatrixSuite_Teardown );
   pcu_suite_addTest( suite, ConstitutiveMatrixSuite_CartesianMatrix2D );
}
