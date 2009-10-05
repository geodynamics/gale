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
} ArrheniusSuiteData;

void ArrheniusSuite_Setup( ArrheniusSuiteData* data ) { 
}

void ArrheniusSuite_Teardown( ArrheniusSuiteData* data ) {
}

void ArrheniusSuite_TestStiffnessMatrix2D( ArrheniusSuiteData* data ) {
	StiffnessMatrix*			stiffnessMatrix;
	SystemLinearEquations*	sle;
	Dictionary*					dictionary;
	UnderworldContext*		context;
	Stg_ComponentFactory*	cf;
	PetscViewer					bviewer;
	PetscReal					matrixNorm, errorNorm, test;
	Mat							expected;
	char							expected_file[PCU_PATH_MAX];
	char							*filename;
	double						tolerance;
	char							xml_input[PCU_PATH_MAX];

	/* read in the xml input file */
	pcu_filename_input( "testArrhenius2D.xml", xml_input );
	cf = stgMainInitFromXML( xml_input, MPI_COMM_WORLD, NULL );
	context = (UnderworldContext*)LiveComponentRegister_Get( cf->LCRegister, "context" );
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

	dictionary = context->dictionary;
	filename = Dictionary_GetString( dictionary, "StiffnessMatrixCompareFilename" );
	tolerance = Dictionary_GetDouble_WithDefault( dictionary, "StiffnessMatrixCompareTolerance", 1e-4 );

	pcu_filename_expected( filename, expected_file );
	PetscViewerBinaryOpen( MPI_COMM_WORLD, expected_file, FILE_MODE_READ, &bviewer );
	MatLoad( bviewer, MATAIJ, &expected );

	MatNorm( expected, NORM_FROBENIUS, &matrixNorm );
	MatAXPY( expected, -1, (stiffnessMatrix->matrix) , SAME_NONZERO_PATTERN );
	MatNorm( expected, NORM_FROBENIUS, &errorNorm );

	assert( matrixNorm != 0 );
	test = errorNorm / matrixNorm;

	pcu_check_lt( test, tolerance );
	stgMainDestroy( cf );

}

void ArrheniusSuite( pcu_suite_t* suite ) {
   pcu_suite_setData( suite, ArrheniusSuiteData );
   pcu_suite_setFixtures( suite, ArrheniusSuite_Setup, ArrheniusSuite_Teardown );
   pcu_suite_addTest( suite, ArrheniusSuite_TestStiffnessMatrix2D );
}
