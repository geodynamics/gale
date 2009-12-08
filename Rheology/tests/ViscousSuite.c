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
} ViscousSuiteData;

void ViscousSuite_Setup( ViscousSuiteData* data ) { 
	data->context=NULL;
}

void ViscousSuite_Teardown( ViscousSuiteData* data ) {
}

void ViscousSuite_ArrheniusStiffnessMatrix2D( ViscousSuiteData* data ) {
	StiffnessMatrix*			stiffnessMatrix;
	SystemLinearEquations*	sle;
	Dictionary*					dictionary;
	UnderworldContext*		context;
	Stg_ComponentFactory*	cf;
	PetscViewer					expViewer;
	PetscReal					matrixNorm, errorNorm, test;
	Mat							expected;
	char							expected_file[PCU_PATH_MAX];
	char							*filename;
	double						tolerance;
	char							xml_input[PCU_PATH_MAX];
	char							rFile[PCU_PATH_MAX];
	int							err;


	pcu_docstring( "This test compares a Stiffness matrix against a previously generated stiffness matrix"
		"The stiffness matrix is generated from a 2D FEM model for a flow with an Arrhenius rheology."
		"See testArrhenius2D.xml for the actual xml used"	);

	/* read in the xml input file */
	pcu_filename_input( "testArrhenius2D.xml", xml_input );
	cf = stgMainInitFromXML( xml_input, MPI_COMM_WORLD, NULL );
	context = (UnderworldContext*)LiveComponentRegister_Get( cf->LCRegister, "context" );
	data->context = context;
	dictionary = context->dictionary;

	stgMainBuildAndInitialise( cf );

	/* Assemble */
	stiffnessMatrix = (StiffnessMatrix*) LiveComponentRegister_Get( context->CF->LCRegister, "k_matrix" );
	sle = (SystemLinearEquations*)       LiveComponentRegister_Get( context->CF->LCRegister, "stokesEqn" );
	assert( stiffnessMatrix );
	StiffnessMatrix_Assemble( stiffnessMatrix, False, sle, context );
	stiffnessMatrix = stiffnessMatrix;

	/* Test is to check the relative error between an
		 1) expected stiffness matrix, (made years ago)
		 2) the current stiffness matrix.

		 both matricies are built using only an Arrhenius rheology 
	 */

	filename = Dictionary_GetString( dictionary, "StiffnessMatrixCompareFilename" );
	tolerance = Dictionary_GetDouble_WithDefault( dictionary, "StiffnessMatrixCompareTolerance", 1e-4 );

	pcu_filename_expected( filename, expected_file );
	PetscViewerBinaryOpen( context->communicator, expected_file, FILE_MODE_READ, &expViewer );

	MatLoad( expViewer, MATAIJ, &expected );

	/* 
		 To view the expected and computed matricies uncomment this
	PetscViewerASCIIOpen(context->communicator, "numerical.dat",&currViewer);
	PetscViewerASCIIOpen(context->communicator, "expected.dat",&parallelViewer);
	MatView( stiffnessMatrix->matrix, currViewer ); //PETSC_VIEWER_STDOUT_WORLD );
	MatView( expected, parallelViewer ); //PETSC_VIEWER_STDOUT_WORLD );
	PetscViewerDestroy(currViewer);
	PetscViewerDestroy(parallelViewer);
	*/

	MatNorm( expected, NORM_FROBENIUS, &matrixNorm );
	assert( matrixNorm != 0 );

	MatAXPY( expected, -1, (stiffnessMatrix->matrix) , DIFFERENT_NONZERO_PATTERN );
	MatNorm( expected, NORM_FROBENIUS, &errorNorm );
	test = errorNorm / matrixNorm;

	pcu_check_lt( test, tolerance );

	if( data->context->rank == 0 ) {
		/* Now clean output path */
		sprintf(rFile, "%s/input.xml", data->context->outputPath );
		err = remove( rFile );
		if( err == -1 ) printf("Error in %s, can't delete the input.xml\n", __func__);
	}

   PetscViewerDestroy( expViewer );
   MatDestroy( expected);
	stgMainDestroy( cf );
}

void ViscousSuite_FrankKamenetskiiStiffnessMatrix2D( ViscousSuiteData* data ) {
	StiffnessMatrix*			stiffnessMatrix;
	SystemLinearEquations*	sle;
	Dictionary*					dictionary;
	UnderworldContext*		context;
	Stg_ComponentFactory*	cf;
	PetscViewer					expViewer;
	PetscReal					matrixNorm, errorNorm, test;
	Mat							expected;
	char							expected_file[PCU_PATH_MAX];
	char							*filename;
	double						tolerance;
	char							xml_input[PCU_PATH_MAX];
	char							rFile[PCU_PATH_MAX];
	int							err;

	pcu_docstring( "This test compares a Stiffness matrix against a previously generated stiffness matrix"
		"The stiffness matrix is generated from a 2D FEM model for a flow with an FrankKamenetskii rheology."
		"See testFrankKamenetskii2D.xml for the actual xml used"	);

	/* read in the xml input file */
	pcu_filename_input( "testFrankKamenetskii2D.xml", xml_input );
	cf = stgMainInitFromXML( xml_input, MPI_COMM_WORLD, NULL );
	context = (UnderworldContext*)LiveComponentRegister_Get( cf->LCRegister, "context" );
	data->context = context;
	dictionary = context->dictionary;

	stgMainBuildAndInitialise( cf );

	/* Assemble */
	stiffnessMatrix = (StiffnessMatrix*) LiveComponentRegister_Get( context->CF->LCRegister, "k_matrix" );
	sle = (SystemLinearEquations*)       LiveComponentRegister_Get( context->CF->LCRegister, "stokesEqn" );
	assert( stiffnessMatrix );
	StiffnessMatrix_Assemble( stiffnessMatrix, False, sle, context );
	stiffnessMatrix = stiffnessMatrix;

	/* Test is to check the relative error between an
		 1) expected stiffness matrix, (made years ago)
		 2) the current stiffness matrix.

		 both matricies are built using only an Arrhenius rheology 
	 */

	filename = Dictionary_GetString( dictionary, "StiffnessMatrixCompareFilename" );
	tolerance = Dictionary_GetDouble_WithDefault( dictionary, "StiffnessMatrixCompareTolerance", 1e-4 );

	pcu_filename_expected( filename, expected_file );
	PetscViewerBinaryOpen( context->communicator, expected_file, FILE_MODE_READ, &expViewer );

	MatLoad( expViewer, MATAIJ, &expected );

	/* 
		 To view the expected and computed matricies uncomment this
	PetscViewerASCIIOpen(context->communicator, "numerical.dat",&currViewer);
	PetscViewerASCIIOpen(context->communicator, "expected.dat",&parallelViewer);
	MatView( stiffnessMatrix->matrix, currViewer ); //PETSC_VIEWER_STDOUT_WORLD );
	MatView( expected, parallelViewer ); //PETSC_VIEWER_STDOUT_WORLD );
	PetscViewerDestroy(currViewer);
	PetscViewerDestroy(parallelViewer);
	*/

	MatNorm( expected, NORM_FROBENIUS, &matrixNorm );
	assert( matrixNorm != 0 );

	MatAXPY( expected, -1, (stiffnessMatrix->matrix) , DIFFERENT_NONZERO_PATTERN );
	MatNorm( expected, NORM_FROBENIUS, &errorNorm );
	test = errorNorm / matrixNorm;

	pcu_check_lt( test, tolerance );

	if( data->context->rank == 0 ) {
		/* Now clean output path */
		sprintf(rFile, "%s/input.xml", data->context->outputPath );
		err = remove( rFile );
		if( err == -1 ) printf("Error in %s, can't delete the input.xml\n", __func__);
	}

   PetscViewerDestroy( expViewer );
   MatDestroy( expected);
	stgMainDestroy( cf );
}

void ViscousSuite_MaterialViscosityStiffnessMatrix2D( ViscousSuiteData* data ) {
	StiffnessMatrix*			stiffnessMatrix;
	SystemLinearEquations*	sle;
	Dictionary*					dictionary;
	UnderworldContext*		context;
	Stg_ComponentFactory*	cf;
	PetscViewer					parallelViewer, expViewer, currViewer;
	PetscReal					matrixNorm, errorNorm, test;
	Mat							expected;
	char							expected_file[PCU_PATH_MAX];
	char							*filename;
	double						tolerance;
	char							xml_input[PCU_PATH_MAX];
	char							rFile[PCU_PATH_MAX];
	int							err;

	pcu_docstring( "This test compares a Stiffness matrix against a previously generated stiffness matrix"
		"The stiffness matrix is generated from a 2D FEM model for a flow with an constant isoviscous rheology."
		"See testMaterialViscosity2D.xml for the actual xml used"	);

	/* read in the xml input file */
	pcu_filename_input( "testMaterialViscosity2D.xml", xml_input );
	cf = stgMainInitFromXML( xml_input, MPI_COMM_WORLD, NULL );
	context = (UnderworldContext*)LiveComponentRegister_Get( cf->LCRegister, "context" );
	data->context = context;
	dictionary = context->dictionary;

	stgMainBuildAndInitialise( cf );

	/* Assemble */
	stiffnessMatrix = (StiffnessMatrix*) LiveComponentRegister_Get( context->CF->LCRegister, "k_matrix" );
	sle = (SystemLinearEquations*)       LiveComponentRegister_Get( context->CF->LCRegister, "stokesEqn" );
	assert( stiffnessMatrix );
	StiffnessMatrix_Assemble( stiffnessMatrix, False, sle, context );
	stiffnessMatrix = stiffnessMatrix;

/* Test is to check the relative error between an
	 1) expected stiffness matrix, (made years ago)
	 2) the current stiffness matrix.

	 both matricies are built using only an Arrhenius rheology 
 */

	filename = Dictionary_GetString( dictionary, "StiffnessMatrixCompareFilename" );
	tolerance = Dictionary_GetDouble_WithDefault( dictionary, "StiffnessMatrixCompareTolerance", 1e-4 );

	pcu_filename_expected( filename, expected_file );
	PetscViewerBinaryOpen( context->communicator, expected_file, FILE_MODE_READ, &expViewer );

	MatLoad( expViewer, MATAIJ, &expected );
	/* 
		 To view the expected and computed matricies uncomment this
	PetscViewerASCIIOpen(context->communicator, "numerical.dat",&currViewer);
	PetscViewerASCIIOpen(context->communicator, "expected.dat",&parallelViewer);
	MatView( stiffnessMatrix->matrix, currViewer ); //PETSC_VIEWER_STDOUT_WORLD );
	MatView( expected, parallelViewer ); //PETSC_VIEWER_STDOUT_WORLD );
	PetscViewerDestroy(currViewer);
	PetscViewerDestroy(parallelViewer);
	*/

	MatNorm( expected, NORM_FROBENIUS, &matrixNorm );
	assert( matrixNorm != 0 );

	MatAXPY( expected, -1, (stiffnessMatrix->matrix) , DIFFERENT_NONZERO_PATTERN );
	MatNorm( expected, NORM_FROBENIUS, &errorNorm );
	test = errorNorm / matrixNorm;

	pcu_check_lt( test, tolerance );

	if( data->context->rank == 0 ) {
		/* Now clean output path */
		sprintf(rFile, "%s/input.xml", data->context->outputPath );
		err = remove( rFile );
		if( err == -1 ) printf("Error in %s, can't delete the input.xml\n", __func__);
	}

   PetscViewerDestroy( expViewer );
   MatDestroy( expected);
	stgMainDestroy( cf );
}

void ViscousSuite_ArrheniusStiffnessMatrix2D_DualMesh( ViscousSuiteData* data ) {
	StiffnessMatrix*			stiffnessMatrix;
	SystemLinearEquations*	sle;
	Dictionary*					dictionary;
	UnderworldContext*		context;
	Stg_ComponentFactory*	cf;
	PetscViewer					expViewer;
	PetscReal					matrixNorm, errorNorm, test;
	Mat							expected;
	char							expected_file[PCU_PATH_MAX];
	char							*filename;
	double						tolerance;
	char							xml_input[PCU_PATH_MAX];
	char							rFile[PCU_PATH_MAX];
	int							err;

	pcu_docstring( "This test compares a Stiffness matrix against a previously generated stiffness matrix"
		"The stiffness matrix is generated from a 2D FEM model for a flow with an FrankKamenetskii rheology."
		"See testFrankKamenetskii2D.xml for the actual xml used"	);

	/* read in the xml input file */
	pcu_filename_input( "testArrhenius2D-DualMesh.xml", xml_input );
	cf = stgMainInitFromXML( xml_input, MPI_COMM_WORLD, NULL );
	context = (UnderworldContext*)LiveComponentRegister_Get( cf->LCRegister, "context" );
	data->context = context;
	dictionary = context->dictionary;

	stgMainBuildAndInitialise( cf );

	/* Assemble */
	stiffnessMatrix = (StiffnessMatrix*) LiveComponentRegister_Get( context->CF->LCRegister, "k_matrix" );
	sle = (SystemLinearEquations*)       LiveComponentRegister_Get( context->CF->LCRegister, "stokesEqn" );
	assert( stiffnessMatrix );
	StiffnessMatrix_Assemble( stiffnessMatrix, False, sle, context );
	stiffnessMatrix = stiffnessMatrix;

	/* Test is to check the relative error between an
		 1) expected stiffness matrix, (made years ago)
		 2) the current stiffness matrix.

		 both matricies are built using only an Arrhenius rheology 
	 */

	filename = Dictionary_GetString( dictionary, "StiffnessMatrixCompareFilename" );
	tolerance = Dictionary_GetDouble_WithDefault( dictionary, "StiffnessMatrixCompareTolerance", 1e-4 );

	pcu_filename_expected( filename, expected_file );
	PetscViewerBinaryOpen( context->communicator, expected_file, FILE_MODE_READ, &expViewer );

	MatLoad( expViewer, MATAIJ, &expected );

	/* 
		 To view the expected and computed matricies uncomment this
	PetscViewerASCIIOpen(context->communicator, "numerical.dat",&currViewer);
	PetscViewerASCIIOpen(context->communicator, "expected.dat",&parallelViewer);
	MatView( stiffnessMatrix->matrix, currViewer ); //PETSC_VIEWER_STDOUT_WORLD );
	MatView( expected, parallelViewer ); //PETSC_VIEWER_STDOUT_WORLD );
	PetscViewerDestroy(currViewer);
	PetscViewerDestroy(parallelViewer);
	*/

	MatNorm( expected, NORM_FROBENIUS, &matrixNorm );
	assert( matrixNorm != 0 );

	MatAXPY( expected, -1, (stiffnessMatrix->matrix) , DIFFERENT_NONZERO_PATTERN );
	MatNorm( expected, NORM_FROBENIUS, &errorNorm );
	test = errorNorm / matrixNorm;

	pcu_check_lt( test, tolerance );

	if( data->context->rank == 0 ) {
		/* Now clean output path */
		sprintf(rFile, "%s/input.xml", data->context->outputPath );
		err = remove( rFile );
		if( err == -1 ) printf("Error in %s, can't delete the input.xml\n", __func__);
	}

   PetscViewerDestroy( expViewer );
   MatDestroy( expected);
	stgMainDestroy( cf );
}

void ViscousSuite_FrankKamenetskiiStiffnessMatrix2D_DualMesh( ViscousSuiteData* data ) {
	StiffnessMatrix*			stiffnessMatrix;
	SystemLinearEquations*	sle;
	Dictionary*					dictionary;
	UnderworldContext*		context;
	Stg_ComponentFactory*	cf;
	PetscViewer					expViewer;
	PetscReal					matrixNorm, errorNorm, test;
	Mat							expected;
	char							expected_file[PCU_PATH_MAX];
	char							*filename;
	double						tolerance;
	char							xml_input[PCU_PATH_MAX];
	char							rFile[PCU_PATH_MAX];
	int							err;

	pcu_docstring( "This test compares a Stiffness matrix against a previously generated stiffness matrix"
		"The stiffness matrix is generated from a 2D FEM model for a flow with an FrankKamenetskii rheology."
		"See testFrankKamenetskii2D.xml for the actual xml used"	);

	/* read in the xml input file */
	pcu_filename_input( "testFrankKamenetskii2D-DualMesh.xml", xml_input );
	cf = stgMainInitFromXML( xml_input, MPI_COMM_WORLD, NULL );
	context = (UnderworldContext*)LiveComponentRegister_Get( cf->LCRegister, "context" );
	data->context = context;
	dictionary = context->dictionary;

	stgMainBuildAndInitialise( cf );

	/* Assemble */
	stiffnessMatrix = (StiffnessMatrix*) LiveComponentRegister_Get( context->CF->LCRegister, "k_matrix" );
	sle = (SystemLinearEquations*)       LiveComponentRegister_Get( context->CF->LCRegister, "stokesEqn" );
	assert( stiffnessMatrix );
	StiffnessMatrix_Assemble( stiffnessMatrix, False, sle, context );
	stiffnessMatrix = stiffnessMatrix;

	/* Test is to check the relative error between an
		 1) expected stiffness matrix, (made years ago)
		 2) the current stiffness matrix.

		 both matricies are built using only an Arrhenius rheology 
	 */

	filename = Dictionary_GetString( dictionary, "StiffnessMatrixCompareFilename" );
	tolerance = Dictionary_GetDouble_WithDefault( dictionary, "StiffnessMatrixCompareTolerance", 1e-4 );

	pcu_filename_expected( filename, expected_file );
	PetscViewerBinaryOpen( context->communicator, expected_file, FILE_MODE_READ, &expViewer );

	MatLoad( expViewer, MATAIJ, &expected );

	/* 
		 To view the expected and computed matricies uncomment this
	PetscViewerASCIIOpen(context->communicator, "numerical.dat",&currViewer);
	PetscViewerASCIIOpen(context->communicator, "expected.dat",&parallelViewer);
	MatView( stiffnessMatrix->matrix, currViewer ); //PETSC_VIEWER_STDOUT_WORLD );
	MatView( expected, parallelViewer ); //PETSC_VIEWER_STDOUT_WORLD );
	PetscViewerDestroy(currViewer);
	PetscViewerDestroy(parallelViewer);
	*/

	MatNorm( expected, NORM_FROBENIUS, &matrixNorm );
	assert( matrixNorm != 0 );

	MatAXPY( expected, -1, (stiffnessMatrix->matrix) , DIFFERENT_NONZERO_PATTERN );
	MatNorm( expected, NORM_FROBENIUS, &errorNorm );
	test = errorNorm / matrixNorm;

	pcu_check_lt( test, tolerance );

	if( data->context->rank == 0 ) {
		/* Now clean output path */
		sprintf(rFile, "%s/input.xml", data->context->outputPath );
		err = remove( rFile );
		if( err == -1 ) printf("Error in %s, can't delete the input.xml\n", __func__);
	}

   PetscViewerDestroy( expViewer );
   MatDestroy( expected);
	stgMainDestroy( cf );
}

void ViscousSuite( pcu_suite_t* suite ) {
   pcu_suite_setData( suite, ViscousSuiteData );
   pcu_suite_setFixtures( suite, ViscousSuite_Setup, ViscousSuite_Teardown );
   pcu_suite_addTest( suite, ViscousSuite_ArrheniusStiffnessMatrix2D );
   pcu_suite_addTest( suite, ViscousSuite_FrankKamenetskiiStiffnessMatrix2D );
   pcu_suite_addTest( suite, ViscousSuite_MaterialViscosityStiffnessMatrix2D );
   pcu_suite_addTest( suite, ViscousSuite_ArrheniusStiffnessMatrix2D_DualMesh );
   pcu_suite_addTest( suite, ViscousSuite_FrankKamenetskiiStiffnessMatrix2D_DualMesh );
}


