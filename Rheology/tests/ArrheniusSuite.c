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
	Stg_ComponentFactory* cf;
	StiffnessMatrix* stiffnessMatrix;
} ArrheniusSuiteData;

void ArrheniusSuite_Setup( ArrheniusSuiteData* data ) { 
	StiffnessMatrix*            stiffnessMatrix;
	SystemLinearEquations*      sle;
	Dictionary*      dictionary;
	char xml_input[PCU_PATH_MAX];

	/* read in the xml input file */
	pcu_filename_input( "testArrhenius2D.xml", xml_input );
	data->cf = stgMainInitFromXML( xml_input, MPI_COMM_WORLD, NULL );
	data->context = (UnderworldContext*)LiveComponentRegister_Get( data->cf->LCRegister, "context" );
	dictionary = data->context->dictionary;

	/* Assemble */
	stiffnessMatrix = (StiffnessMatrix*) LiveComponentRegister_Get( data->context->CF->LCRegister, "k_matrix" );
	sle = (SystemLinearEquations*)       LiveComponentRegister_Get( data->context->CF->LCRegister, "stokesEqn" );
	assert( stiffnessMatrix );
	StiffnessMatrix_Assemble( stiffnessMatrix, False, sle, data->context );
	data->stiffnessMatrix = stiffnessMatrix;
}

void ArrheniusSuite_Teardown( ArrheniusSuiteData* data ) {
	Stg_Component_Destroy( data->context, 0 /* dummy */, False );
	Stg_Class_Delete( data->context );
}

void ArrheniusSuite_TestStiffnessMatrix2D( ArrheniusSuiteData* data ) {
	/* Test is to check the relative error between an
		 1) expected stiffness matrix, (made years ago)
		 2) the current stiffness matrix.

		 both matricies are built using only an Arrhenius rheology 
	 */
	StiffnessMatrix *stiffnessMatrix = data->stiffnessMatrix;
	Dictionary* dictionary = data->context->dictionary;
	PetscViewer bviewer;
	Mat expected;
	char expected_file[PCU_PATH_MAX];
	char *filename = Dictionary_GetString( dictionary, "StiffnessMatrixCompareFilename" );
	double tolerance = Dictionary_GetDouble_WithDefault( dictionary, "StiffnessMatrixCompareTolerance", 1e-4 );
	PetscReal matrixNorm, errorNorm, test;

	pcu_filename_expected( filename, expected_file );
	PetscViewerBinaryOpen(MPI_COMM_WORLD, expected_file, FILE_MODE_READ, &bviewer);
	MatLoad(bviewer, MATAIJ, &expected);

	MatNorm( expected, NORM_FROBENIUS, &matrixNorm );
	MatAXPY(expected, -1, (stiffnessMatrix->matrix) , SAME_NONZERO_PATTERN);
	MatNorm( expected, NORM_FROBENIUS, &errorNorm );

	assert( matrixNorm != 0 );
	test = errorNorm / matrixNorm;

	pcu_check_lt( test, tolerance );
}

void ArrheniusSuite( pcu_suite_t* suite ) {
   pcu_suite_setData( suite, ArrheniusSuiteData );
   pcu_suite_setFixtures( suite, ArrheniusSuite_Setup, ArrheniusSuite_Teardown );
   pcu_suite_addTest( suite, ArrheniusSuite_TestStiffnessMatrix2D );
}
