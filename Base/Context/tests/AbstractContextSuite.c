/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street, Melbourne, 3053, Australia.
**
** Authors:
**   Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**   Patrick D. Sunter, Software Engineer, VPAC. (pds@vpac.org)
**   Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**   Siew-Ching Tan, Software Engineer, VPAC. (siew@vpac.org)
**   Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
**   Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
**
**  This library is free software; you can redistribute it and/or
**  modify it under the terms of the GNU Lesser General Public
**  License as published by the Free Software Foundation; either
**  version 2.1 of the License, or (at your option) any later version.
**
**  This library is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
**  Lesser General Public License for more details.
**
**  You should have received a copy of the GNU Lesser General Public
**  License along with this library; if not, write to the Free Software
**  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
**
**
** $Id: testJournal-Dictionary.c 2745 2005-03-05 08:12:18Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "pcu/pcu.h"
#include "StGermain/Base/Foundation/Foundation.h"
#include "StGermain/Base/IO/IO.h"
#include "StGermain/Base/Container/Container.h"
#include "StGermain/Base/Automation/Automation.h"
#include "StGermain/Base/Extensibility/Extensibility.h"
#include "StGermain/Base/Context/Context.h"
#include "AbstractContextSuite.h"

/* Temporarily needed until EP shortcuts are fixed */
#define  CURR_MODULE_NAME "AbstractContextSuite"

/* need to allocate memory for this stream */
Stream* stream;

#define __TestContext \
	__AbstractContext \
	unsigned int	buildHookCalled; \
	unsigned int	icHookCalled; \
	unsigned int	dtHookCalled; \
	unsigned int	solveHookCalled; \
	unsigned int	solve2HookCalled; \
	unsigned int	syncHookCalled; \
	unsigned int	outputHookCalled; \
	unsigned int	dumpHookCalled; \
	unsigned int	checkpointHookCalled; \
	double			computedValue; 
struct TestContext { __TestContext };
typedef struct TestContext TestContext;

double dt = 2.0f;
#define MAX_TIME_STEPS 1000
double GLOBAL_COMP_VALUE[MAX_TIME_STEPS];

void TestSetDt( void* context, double _dt ) {
   dt = _dt;
}


typedef struct {
   TestContext*	ctx;
   Dictionary*		dict;
} AbstractContextSuiteData;


TestContext* TestContext_New(
	Name			name,
	double		startTime,
	double		stopTime,
	MPI_Comm		communicator,
	Dictionary*	dictionary ) 
{
   TestContext* ctx;

   ctx = (TestContext*)_AbstractContext_New( 
      sizeof(TestContext), 
      "TestContext", 
      _AbstractContext_Delete, 
      _AbstractContext_Print, 
      NULL,
      NULL, 
      NULL, 
      _AbstractContext_Build, 
      _AbstractContext_Initialise, 
      _AbstractContext_Execute, 
      _AbstractContext_Destroy, 
      name, 
      NON_GLOBAL, 
      TestSetDt, 
      startTime, 
      stopTime, 
      communicator, 
      dictionary );

   ctx->buildHookCalled = 0;  
   ctx->icHookCalled = 0;
   ctx->dtHookCalled = 0;
   ctx->solveHookCalled = 0;
   ctx->solve2HookCalled = 0;
   ctx->syncHookCalled = 0;
   ctx->outputHookCalled = 0;
   ctx->dumpHookCalled = 0;
   ctx->checkpointHookCalled = 0;
   ctx->computedValue = 0;
   
   return ctx;
}


void TestBuild( void* context ) {
   TestContext* self = (TestContext*)context;
   self->buildHookCalled++;
}

void TestInitialConditions( void* context ) {
   TestContext* self = (TestContext*)context;
   self->icHookCalled++;
   /* Since the current convention for loading from checkpoint is that there's no special entry point and 
    * it should be done by the init() (possibly in data objects themselves), follow that here */   
   if (self->loadFromCheckPoint) {
      self->computedValue = GLOBAL_COMP_VALUE[self->restartTimestep];
   }
}


double TestDt( void* context ) {
   TestContext* self = (TestContext*)context;
   self->dtHookCalled++;
   return dt;
}

void TestSolve( void* context ) {
   TestContext* self = (TestContext*)context;
   self->solveHookCalled++;
   self->computedValue = pow( 1.1, self->timeStep );
}

void TestSolve2( void* context ) {
   TestContext* self = (TestContext*)context;
   self->solve2HookCalled++;
}

void TestSync( void* context ) {
   TestContext* self = (TestContext*)context;
   self->syncHookCalled++;
}

void TestOutput( void* context ) {
   TestContext* self = (TestContext*)context;
   self->outputHookCalled++;
}

void TestCheckpoint( void* context ) {
   TestContext* self = (TestContext*)context;
   self->checkpointHookCalled++;
   GLOBAL_COMP_VALUE[self->timeStep] = self->computedValue;
}

void TestDump( void* context ) {
   TestContext* self = (TestContext*)context;
   self->dumpHookCalled++;
}

void AbstractContextSuite_Setup( AbstractContextSuiteData* data ) {
   Stg_ComponentFactory* cf;
   MPI_Comm       CommWorld;
   Index          ii;

   data->dict = Dictionary_New();

   Dictionary_Add( data->dict, "outputPath", Dictionary_Entry_Value_FromString( "output" ) );
   Dictionary_Add( data->dict, "checkpointEvery", Dictionary_Entry_Value_FromUnsignedInt( 5 ) );
   Dictionary_Add( data->dict, "dumpEvery", Dictionary_Entry_Value_FromUnsignedInt( 2 ) );
   Dictionary_Add( data->dict, "maxTimeSteps", Dictionary_Entry_Value_FromUnsignedInt( 10 ) );
   
   cf = Stg_ComponentFactory_New( data->dict, NULL );

   for (ii=0; ii < MAX_TIME_STEPS; ii++) {
      GLOBAL_COMP_VALUE[ii] = 0.0;
   }

   MPI_Comm_dup( MPI_COMM_WORLD, &CommWorld );
   /* Build the context */
   data->ctx = TestContext_New( 
      "context", 
      0, 
      0, 
      CommWorld, 
      data->dict );

	_AbstractContext_Init( data->ctx );
   _AbstractContext_AssignFromXML( data->ctx, cf, NULL );

   Stream_Enable( data->ctx->info, False );
}

void AbstractContextSuite_Teardown( AbstractContextSuiteData* data ) {
   Stg_Class_Delete( data->dict );
}

  
void AbstractContextSuite_TestDefaultEPs( AbstractContextSuiteData* data ) {
   ContextEntryPoint*      contextEP=NULL;

   /* Assert that default EPs are set up correctly, eg for saving */
   contextEP = (ContextEntryPoint*)AbstractContext_GetEntryPoint( data->ctx, "Context_Construct" );
   pcu_check_true( contextEP->hooks->count == 1 );
   pcu_check_streq( ((Hook*)contextEP->hooks->data[0])->name, "default" );
   contextEP = (ContextEntryPoint*)AbstractContext_GetEntryPoint( data->ctx, "Context_ConstructExtensions" );
   pcu_check_true( contextEP->hooks->count == 0 );
   contextEP = (ContextEntryPoint*)AbstractContext_GetEntryPoint( data->ctx, "Context_Build" );
   pcu_check_true( contextEP->hooks->count == 0 );
   //pcu_check_streq( ((Hook*)contextEP->hooks->data[0])->name, "BuildAllLiveComponents" );
   contextEP = (ContextEntryPoint*)AbstractContext_GetEntryPoint( data->ctx, "Context_Initialise" );
   pcu_check_true( contextEP->hooks->count == 0 );
   //pcu_check_streq( ((Hook*)contextEP->hooks->data[0])->name, "InitialiseAllLiveComponents" );
   contextEP = (ContextEntryPoint*)AbstractContext_GetEntryPoint( data->ctx, "Context_Execute" );
   pcu_check_true( contextEP->hooks->count == 1 );
   pcu_check_streq( ((Hook*)contextEP->hooks->data[0])->name, "default" );
   contextEP = (ContextEntryPoint*)AbstractContext_GetEntryPoint( data->ctx, "Context_Destroy" );
   pcu_check_true( contextEP->hooks->count == 0 );
   contextEP = (ContextEntryPoint*)AbstractContext_GetEntryPoint( data->ctx, "Context_DestroyExtensions" );
   pcu_check_true( contextEP->hooks->count == 0 );
   contextEP = (ContextEntryPoint*)AbstractContext_GetEntryPoint( data->ctx, "Context_Dt" );
   pcu_check_true( contextEP->hooks->count == 0 );
   contextEP = (ContextEntryPoint*)AbstractContext_GetEntryPoint( data->ctx, "Context_Step" );
   pcu_check_true( contextEP->hooks->count == 1 );
   pcu_check_streq( ((Hook*)contextEP->hooks->data[0])->name, "default" );
   contextEP = (ContextEntryPoint*)AbstractContext_GetEntryPoint( data->ctx, "Context_Solve" );
   pcu_check_true( contextEP->hooks->count == 0 );
   contextEP = (ContextEntryPoint*)AbstractContext_GetEntryPoint( data->ctx, "Context_UpdateClass" );
   pcu_check_true( contextEP->hooks->count == 0 );
   contextEP = (ContextEntryPoint*)AbstractContext_GetEntryPoint( data->ctx, "Context_Sync" );
   pcu_check_true( contextEP->hooks->count == 0 );
   contextEP = (ContextEntryPoint*)AbstractContext_GetEntryPoint( data->ctx, "Context_FrequentOutput" );
   pcu_check_true( contextEP->hooks->count == 0 );
   contextEP = (ContextEntryPoint*)AbstractContext_GetEntryPoint( data->ctx, "Context_Dump" );
   pcu_check_true( contextEP->hooks->count == 0 );
   contextEP = (ContextEntryPoint*)AbstractContext_GetEntryPoint( data->ctx, "Context_DumpClass" );
   pcu_check_true( contextEP->hooks->count == 0 );
   contextEP = (ContextEntryPoint*)AbstractContext_GetEntryPoint( data->ctx, "Context_Save" );
   pcu_check_true( contextEP->hooks->count == 2 );
   pcu_check_streq( ((Hook*)contextEP->hooks->data[0])->name, "CreateCheckpointDirectory" );
   pcu_check_streq( ((Hook*)contextEP->hooks->data[1])->name, "SaveTimeInfo" );
   contextEP = (ContextEntryPoint*)AbstractContext_GetEntryPoint( data->ctx, "Context_SaveClass" );
   pcu_check_true( contextEP->hooks->count == 0 );
}


void AbstractContextSuite_TestRunBasic( AbstractContextSuiteData* data ) {
   ContextEP_ReplaceAll( data->ctx, AbstractContext_EP_Build, TestBuild );
   ContextEP_ReplaceAll( data->ctx, AbstractContext_EP_Initialise, TestInitialConditions );
   ContextEP_ReplaceAll( data->ctx, AbstractContext_EP_Solve, TestSolve ); 
   ContextEP_ReplaceAll( data->ctx, AbstractContext_EP_Sync, TestSync ); 
   ContextEP_ReplaceAll( data->ctx, AbstractContext_EP_Dt, TestDt ); 
   ContextEP_ReplaceAll( data->ctx, AbstractContext_EP_Dump, TestDump ); 
   ContextEP_ReplaceAll( data->ctx, AbstractContext_EP_Save, TestCheckpoint ); 
   ContextEP_ReplaceAll( data->ctx, AbstractContext_EP_FrequentOutput, TestOutput ); 
   ContextEP_Append( data->ctx, AbstractContext_EP_Solve, TestSolve2 );

   Stg_Component_Build( data->ctx, 0 /* dummy */, False );
   Stg_Component_Initialise( data->ctx, 0 /* dummy */, False );
   Stg_Component_Execute( data->ctx, 0 /* dummy */, False );

   pcu_check_true( data->ctx->buildHookCalled == 1 ); 
   pcu_check_true( data->ctx->icHookCalled == 1 );
   pcu_check_true( data->ctx->dtHookCalled == 10 );
   pcu_check_true( data->ctx->solveHookCalled == 10 );
   pcu_check_true( data->ctx->solve2HookCalled == 10 );
   pcu_check_true( data->ctx->syncHookCalled == 10 );
   pcu_check_true( data->ctx->outputHookCalled == 10 );
   pcu_check_true( data->ctx->dumpHookCalled == 10/Dictionary_GetUnsignedInt(data->dict, "dumpEvery" ) );
   pcu_check_true( data->ctx->checkpointHookCalled == 10/Dictionary_GetUnsignedInt(data->dict, "checkpointEvery" ) );

   Stg_Component_Destroy( data->ctx, 0 /* dummy */, False );
}


void AbstractContextSuite_TestRunNoDtDefined( AbstractContextSuiteData* data ) {
   ContextEP_ReplaceAll( data->ctx, AbstractContext_EP_Build, TestBuild );
   ContextEP_ReplaceAll( data->ctx, AbstractContext_EP_Initialise, TestInitialConditions );
   ContextEP_ReplaceAll( data->ctx, AbstractContext_EP_Solve, TestSolve ); 

   ContextEP_Append( data->ctx, AbstractContext_EP_Solve, TestSolve2 );

   Stg_Component_Build( data->ctx, 0 /* dummy */, False );
   Stg_Component_Initialise( data->ctx, 0 /* dummy */, False );

   stJournal->enable = False;
   pcu_check_assert( Stg_Component_Execute( data->ctx, 0 /* dummy */, False ) );
   stJournal->enable = True;

   Stg_Component_Destroy( data->ctx, 0 /* dummy */, False );
}


void AbstractContextSuite_TestRestartFromCheckpoint( AbstractContextSuiteData* data ) {
   Stg_ComponentFactory* cf;
   MPI_Comm CommWorld;

   ContextEP_ReplaceAll( data->ctx, AbstractContext_EP_Build, TestBuild );
   ContextEP_ReplaceAll( data->ctx, AbstractContext_EP_Initialise, TestInitialConditions );
   ContextEP_ReplaceAll( data->ctx, AbstractContext_EP_Solve, TestSolve ); 
   ContextEP_ReplaceAll( data->ctx, AbstractContext_EP_Dt, TestDt ); 

   Stg_Component_Build( data->ctx, 0 /* dummy */, False );
   Stg_Component_Initialise( data->ctx, 0 /* dummy */, False );
   Stg_Component_Execute( data->ctx, 0 /* dummy */, False );
   Stg_Component_Destroy( data->ctx, 0 /* dummy */, False );

   /* ReBuild the context */
   Dictionary_Set( data->dict, "maxTimeSteps", Dictionary_Entry_Value_FromUnsignedInt( 20 ) );
   Dictionary_Set( data->dict, "restartTimestep", Dictionary_Entry_Value_FromUnsignedInt( 5 ) );
   MPI_Comm_dup( MPI_COMM_WORLD, &CommWorld );
   cf = Stg_ComponentFactory_New( data->dict, NULL );

   data->ctx = TestContext_New( 
      "context", 
      0, 
      0, 
      CommWorld, 
      data->dict );

	 _AbstractContext_Init( data->ctx );
   _AbstractContext_AssignFromXML( data->ctx, cf, NULL );

   Stream_Enable( data->ctx->info, False );

   /* add hooks to existing entry points */
   ContextEP_ReplaceAll( data->ctx, AbstractContext_EP_Build, TestBuild );
   ContextEP_ReplaceAll( data->ctx, AbstractContext_EP_Initialise, TestInitialConditions );
   ContextEP_ReplaceAll( data->ctx, AbstractContext_EP_Solve, TestSolve );
   ContextEP_ReplaceAll( data->ctx, AbstractContext_EP_Dt, TestDt );

   /* Run the context for the second time */
   Stg_Component_Build( data->ctx, 0 /* dummy */, False );
   Stg_Component_Initialise( data->ctx, 0 /* dummy */, False );
   Stg_Component_Execute( data->ctx, 0 /* dummy */, False );

   /* As a fairly simple test of basic CP infrastructure, we know that the timesteps should equal
    *  run1_ts + run2_ts, and computed value should equal 1.1 to power (run1_ts + run2_ts) */  
   pcu_check_true( data->ctx->timeStep == (5 + 20) );
   pcu_check_true( abs(data->ctx->computedValue - pow( 1.1, (5 + 20) )) < 1e-8 );

   Stg_Component_Destroy( data->ctx, 0 /* dummy */, False );
}


void AbstractContextSuite( pcu_suite_t* suite ) {
   pcu_suite_setData( suite, AbstractContextSuiteData );
   pcu_suite_setFixtures( suite, AbstractContextSuite_Setup, AbstractContextSuite_Teardown );
   pcu_suite_addTest( suite, AbstractContextSuite_TestDefaultEPs );
   pcu_suite_addTest( suite, AbstractContextSuite_TestRunBasic );
   pcu_suite_addTest( suite, AbstractContextSuite_TestRunNoDtDefined );
   pcu_suite_addTest( suite, AbstractContextSuite_TestRestartFromCheckpoint );
}
