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

const Type Underworld_testYieldCriterion_Type = "Underworld_testYieldCriterion";

/* Define plugin structure */
typedef struct {
	__Codelet
	YieldRheology_HasYieldedFunction* realHasYieldedFunction;
	FeMesh* 		          mesh;
	XYZ                   min;
	XYZ                   max;
	Bool                  hasYielded;
} Underworld_testYieldCriterion;

/* define global vars .... that's very ugly, but it's only a pcu test */
Underworld_testYieldCriterion globalSelf;

typedef struct {
} YieldCriterionSuiteData;

void YieldCriterionSuite_Setup( YieldCriterionSuiteData* data ) { 
}

void YieldCriterionSuite_Teardown( YieldCriterionSuiteData* data ) {
}

double Underworld_testYieldCriterion_dt( FiniteElementContext* context ) {
	if ( context->currentTime >= 0.65 ) {
		return 0.01;
	}
	return HUGE_VAL;
}

void testYieldCriterion_HasYielded( 
		void*                            yieldRheology,
		ConstitutiveMatrix*              constitutiveMatrix,
		MaterialPointsSwarm*             materialPointsSwarm,
		Element_LocalIndex               lElement_I,
		MaterialPoint*                   materialPoint,
		double                           yieldCriterion,
		double                           yieldIndicator ) 
{
	Dimension_Index dim_I;

	/* Call real 'HasYielded' function */
	globalSelf.realHasYieldedFunction( 
			yieldRheology, constitutiveMatrix, materialPointsSwarm, lElement_I, materialPoint, yieldCriterion, yieldIndicator );

	/* Don't output information if this is the first non-linear iteration */
	if ( constitutiveMatrix->sleNonLinearIteration_I == 0 ) {
		return;
	}

	/* Store information */
	globalSelf.hasYielded = True;
	for ( dim_I = 0 ; dim_I < constitutiveMatrix->dim ; dim_I++ ) {
		if ( materialPoint->coord[ dim_I ] < globalSelf.min[ dim_I ] )
			globalSelf.min[ dim_I ] = materialPoint->coord[ dim_I ];
		if ( materialPoint->coord[ dim_I ] > globalSelf.max[ dim_I ] )
			globalSelf.max[ dim_I ] = materialPoint->coord[ dim_I ];
	}
}


void Underworld_testYieldCriterion_Check( FiniteElementContext* context ) {
	Stream* stream = Journal_Register( Dump_Type, Underworld_testYieldCriterion_Type );

	/* Don't do anything if nothing has yielded yet */
	if ( !globalSelf.hasYielded ) {
		return;
	}

	/* Get Calculation to stop */
	context->maxTimeSteps = context->timeStep;

	/* Set the stream to point to our output file (so we can do a diff on it later) */
	Stream_Enable( stream, True );
	Stream_RedirectFile_WithPrependedPath( stream, context->outputPath, "testYieldCriterion.dat" );

	Journal_Printf( stream, "Material yielded at time %4g (step %u) within:\n", context->currentTime, context->timeStep ); 

	/* Output information */
	Journal_Printf( stream, "\tx: %12.4g - %12.4g\n", globalSelf.min[ I_AXIS ], globalSelf.max[ I_AXIS ] );
	Journal_Printf( stream, "\ty: %12.4g - %12.4g\n", globalSelf.min[ J_AXIS ], globalSelf.max[ J_AXIS ] );
	if ( context->dim == 3 ) {
		Journal_Printf( stream, "\tz: %12.4g - %12.4g\n", globalSelf.min[ K_AXIS ], globalSelf.max[ K_AXIS ] );
	}
}

void YieldCriterionSuite_Byerlee2D( YieldCriterionSuiteData* data ) {
	UnderworldContext* context;
	Dictionary*					dictionary;
	YieldRheology*          yieldRheology;
	Stg_ComponentFactory*	cf;
	char							expected_file[PCU_PATH_MAX];
	char							*filename;
	double						tolerance;
	char							xml_input[PCU_PATH_MAX];

	/* read in the xml input file */
	pcu_filename_input( "testByerleeYieldCriterion.xml", xml_input );
	cf = stgMainInitFromXML( xml_input, MPI_COMM_WORLD, NULL );
	context = (UnderworldContext*)LiveComponentRegister_Get( cf->LCRegister, "context" );
	dictionary = context->dictionary;

	stgMainBuildAndInitialise( cf );


	/* get pointer to the mesh */
	globalSelf.mesh = Stg_ComponentFactory_ConstructByName( cf, "linearMesh", FeMesh, True, NULL ); 
	
	/* Get a pointer the yield rheology that we are trying to test */
	yieldRheology = (YieldRheology*) LiveComponentRegister_Get( context->CF->LCRegister, "yieldRheology" );
	
	/* Store the pointer to the original 'HasYielded' function */
	globalSelf.realHasYieldedFunction = yieldRheology->_hasYielded;

	/* Reset this function pointer with our own */
	yieldRheology->_hasYielded = testYieldCriterion_HasYielded;

	ContextEP_Append( context, AbstractContext_EP_Step, Underworld_testYieldCriterion_Check );
	EP_AppendClassHook( Context_GetEntryPoint( context, FiniteElementContext_EP_CalcDt ), 
			Underworld_testYieldCriterion_dt, context );

	globalSelf.min[ I_AXIS ] = HUGE_VAL;
	globalSelf.min[ J_AXIS ] = HUGE_VAL;
	globalSelf.min[ K_AXIS ] = HUGE_VAL;
	globalSelf.max[ I_AXIS ] = -HUGE_VAL;
	globalSelf.max[ J_AXIS ] = -HUGE_VAL;
	globalSelf.max[ K_AXIS ] = -HUGE_VAL;

	stgMainLoop( cf );
	/* 
		 Somewhere here we do a test
	pcu_filename_expected( filename, expected_file );

	pcu_check_lt( test, tolerance );
	*/
	stgMainDestroy( cf );

}

void YieldCriterionSuite( pcu_suite_t* suite ) {
   pcu_suite_setData( suite, YieldCriterionSuiteData );
   pcu_suite_setFixtures( suite, YieldCriterionSuite_Setup, YieldCriterionSuite_Teardown );
   pcu_suite_addTest( suite, YieldCriterionSuite_Byerlee2D );
}
