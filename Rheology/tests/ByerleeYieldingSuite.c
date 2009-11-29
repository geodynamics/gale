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

const Type ByerleeYieldingSuite_testByerlee2D_Type = "Underworld_testYieldCriterion";

/* Define plugin structure */
typedef struct {
   __Codelet
   YieldRheology_HasYieldedFunction* realHasYieldedFunction;
   FeMesh*               mesh;
   XYZ                   min;
   XYZ                   max;
   Bool                  hasYielded;
} ByerleeYielding2D_testYieldCriterion;

/* define global vars .... that's very ugly, but it's only a pcu test */
ByerleeYielding2D_testYieldCriterion globalSelf;

typedef struct {
} ByerleeYieldingSuiteData;

void ByerleeYieldingSuite_Setup( ByerleeYieldingSuiteData* data ) { 
}

void ByerleeYieldingSuite_Teardown( ByerleeYieldingSuiteData* data ) {
}

double ByerleeYieldingSuite_dt( FiniteElementContext* context ) {
	if ( context->currentTime >= 0.65 ) {
		return 0.01;
	}
	return HUGE_VAL;
}

void testByerlee2D_HasYielded( 
		void*						yieldRheology,
		ConstitutiveMatrix*	constitutiveMatrix,
		MaterialPointsSwarm*	materialPointsSwarm,
		Element_LocalIndex	lElement_I,
		MaterialPoint*			materialPoint,
		double					yieldCriterion,
		double					yieldIndicator ) 
{
	Dimension_Index dim_I;

	/* Call real 'HasYielded' function */
	globalSelf.realHasYieldedFunction( yieldRheology, constitutiveMatrix, materialPointsSwarm, lElement_I, materialPoint, yieldCriterion, yieldIndicator );

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


void ByerleeYieldingSuite_Check( FiniteElementContext* context ) {
	Stream* stream = Journal_Register( Dump_Type, ByerleeYieldingSuite_testByerlee2D_Type );

   if(context->rank == 0){
      /* Don't do anything if nothing has yielded yet */
      if ( !globalSelf.hasYielded ) {
         return;
      }
   
      /* Set the stream to point to our output file (so we can do a diff on it later) */
      Stream_Enable( stream, True );
      Stream_RedirectFile_WithPrependedPath( stream, context->outputPath, "testByerlee2D.dat" );
   
      Journal_Printf( stream, "Material yielded at time %12.4g (step %u) within:\n", context->currentTime, context->timeStep ); 
   
      /* Output information */
      Journal_Printf( stream, "\tx: %12.4g - %12.4g\n", globalSelf.min[ I_AXIS ], globalSelf.max[ I_AXIS ] );
      Journal_Printf( stream, "\ty: %12.4g - %12.4g\n", globalSelf.min[ J_AXIS ], globalSelf.max[ J_AXIS ] );
      if ( context->dim == 3 ) {
         Journal_Printf( stream, "\tz: %12.4g - %12.4g\n", globalSelf.min[ K_AXIS ], globalSelf.max[ K_AXIS ] );
      }
      Stream_CloseAndFreeFile( stream );
   }
}

void ByerleeYieldingSuite_Check_Sync( FiniteElementContext* context ) {
   Coord mins;
   Coord maxs;
   
   /* get true global max/min */

   MPI_Allreduce( &globalSelf.min[ I_AXIS ], &mins[I_AXIS], 1, MPI_DOUBLE, MPI_MIN, context->communicator );
   MPI_Allreduce( &globalSelf.min[ J_AXIS ], &mins[J_AXIS], 1, MPI_DOUBLE, MPI_MIN, context->communicator );
   if ( context->dim == 3 )
      MPI_Allreduce( &globalSelf.min[ K_AXIS ], &mins[K_AXIS], 1, MPI_DOUBLE, MPI_MIN, context->communicator );

   MPI_Allreduce( &globalSelf.max[ I_AXIS ], &maxs[I_AXIS], 1, MPI_DOUBLE, MPI_MAX, context->communicator );
   MPI_Allreduce( &globalSelf.max[ J_AXIS ], &maxs[J_AXIS], 1, MPI_DOUBLE, MPI_MAX, context->communicator );
   if ( context->dim == 3 )
      MPI_Allreduce( &globalSelf.max[ K_AXIS ], &maxs[K_AXIS], 1, MPI_DOUBLE, MPI_MAX, context->communicator );
   
   /* now set to these values */
   globalSelf.min[ I_AXIS ] = mins[ I_AXIS ];
   globalSelf.min[ J_AXIS ] = mins[ J_AXIS ];
   if ( context->dim == 3 )
      globalSelf.min[ K_AXIS ] = mins[ K_AXIS ];
   globalSelf.max[ I_AXIS ] = maxs[ I_AXIS ];
   globalSelf.max[ J_AXIS ] = maxs[ J_AXIS ];
   if ( context->dim == 3 )
      globalSelf.max[ K_AXIS ] = maxs[ K_AXIS ];
}

void ByerleeYieldingSuite_TestByerlee2D( ByerleeYieldingSuiteData* data ) {
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
	yieldRheology->_hasYielded = testByerlee2D_HasYielded;

	/* Get Calculation after one step */
	context->maxTimeSteps = 1;

	ContextEP_Append( context, AbstractContext_EP_Step, ByerleeYieldingSuite_Check_Sync );
	ContextEP_Append( context, AbstractContext_EP_Step, ByerleeYieldingSuite_Check );
	EP_AppendClassHook( Context_GetEntryPoint( context, FiniteElementContext_EP_CalcDt ), ByerleeYieldingSuite_dt, context );

	globalSelf.min[ I_AXIS ] = HUGE_VAL;
	globalSelf.min[ J_AXIS ] = HUGE_VAL;
	globalSelf.min[ K_AXIS ] = HUGE_VAL;
	globalSelf.max[ I_AXIS ] = -HUGE_VAL;
	globalSelf.max[ J_AXIS ] = -HUGE_VAL;
	globalSelf.max[ K_AXIS ] = -HUGE_VAL;

	//stgMainLoop( cf );
	Stg_Component_Execute( context, 0, True );

	MPI_Barrier( context->communicator );
	pcu_filename_expected( "testByerleeYieldCriterion.expected", expected_file );
	pcu_check_fileEq( "output/testByerlee2D.dat", expected_file );

   if(context->rank == 0){
      remove("output/testByerlee2D.dat");
      remove("output");
   }

	stgMainDestroy( cf );
}

void ByerleeYieldingSuite( pcu_suite_t* suite ) {
   pcu_suite_setData( suite, ByerleeYieldingSuiteData );
   pcu_suite_setFixtures( suite, ByerleeYieldingSuite_Setup, ByerleeYieldingSuite_Teardown );
   pcu_suite_addTest( suite, ByerleeYieldingSuite_TestByerlee2D );
}


