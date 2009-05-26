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

#include "pcu/pcu.h"
#include "StGermain/Base/Foundation/Foundation.h"
#include "StGermain/Base/IO/IO.h"
#include "StGermain/Base/Container/Container.h"
#include "StGermain/Base/Automation/Automation.h"
#include "StGermain/Base/Extensibility/Extensibility.h"
#include "EntryPointSuite.h"

Stream* stream;


#define NUM_TEST_FUNCS 10

typedef struct {
   EntryPoint* ep;
   Bool        testFuncsRan[NUM_TEST_FUNCS];
} EntryPointSuiteData;

void TestHook0( EntryPointSuiteData* data ) {
   data->testFuncsRan[0] = True;
}

void TestHook1( EntryPointSuiteData* data ) {
   data->testFuncsRan[1] = True;
}

void TestHook2( EntryPointSuiteData* data ) {
   data->testFuncsRan[2] = True;
}

void TestHook3( EntryPointSuiteData* data ) {
   data->testFuncsRan[3] = True;
}

void TestHook4( EntryPointSuiteData* data ) {
   data->testFuncsRan[4] = True;
}

void TestHook5( EntryPointSuiteData* data ) {
   data->testFuncsRan[5] = True;
}

void TestHook6( EntryPointSuiteData* data ) {
   data->testFuncsRan[6] = True;
}

void TestHook7( EntryPointSuiteData* data ) {
   data->testFuncsRan[7] = True;
}

void TestHook8( EntryPointSuiteData* data ) {
   data->testFuncsRan[8] = True;
}

void TestHook9( EntryPointSuiteData* data ) {
   data->testFuncsRan[9] = True;
}


void EntryPointSuite_Setup( EntryPointSuiteData* data ) {
   Index ii;

   data->ep = NULL;
   for (ii=0; ii < NUM_TEST_FUNCS; ii++ ) {
      data->testFuncsRan[ii] = False;
   }
}


void EntryPointSuite_Teardown( EntryPointSuiteData* data ) {
   Stg_Class_Delete( data->ep );
}


void EntryPointSuite_TestRunEmpty( EntryPointSuiteData* data ) {
   data->ep = EntryPoint_New( "test", EntryPoint_VoidPtr_CastType );
   pcu_check_true( data->ep->hooks->count == 0 );
   ((EntryPoint_VoidPtr_CallCast*) data->ep->run)( data->ep, NULL );
}


void EntryPointSuite_TestAppendPrepend( EntryPointSuiteData* data ) {
   Index    ii=0;

   data->ep = EntryPoint_New( "testEP", EntryPoint_0_CastType );
   EntryPoint_Append( data->ep, "TestHook0", (void*)TestHook0, "testCode" );
   /* TestHook0 */
   EntryPoint_Prepend( data->ep, "TestHook1", (void*)TestHook1, "testCode" );
   /* TestHook1, TestHook0 */

   pcu_check_true( data->ep->hooks->count == 2 );
   pcu_check_true( 0 == strcmp( data->ep->hooks->data[0]->name, "TestHook1" ) );
   pcu_check_true( ((Hook*)data->ep->hooks->data[0])->funcPtr == TestHook1 );
   pcu_check_true( 0 == strcmp( data->ep->hooks->data[1]->name, "TestHook0" ) );
   pcu_check_true( ((Hook*)data->ep->hooks->data[1])->funcPtr == TestHook0 );
   for (ii=0; ii < data->ep->hooks->count; ii++ ) {
      pcu_check_true( 0 == strcmp( ((Hook*)data->ep->hooks->data[ii])->addedBy, "testCode" ) );
   }
}


void EntryPointSuite_TestInsertBeforeAfterReplace( EntryPointSuiteData* data ) {
   data->ep = EntryPoint_New( "testEP", EntryPoint_0_CastType );
   EntryPoint_Append( data->ep, "TestHook2", (void*)TestHook2, "testCode" );
   EntryPoint_Append( data->ep, "TestHook3", (void*)TestHook3, "testCode" );
   /* TestHook2, TestHook3 */
   EntryPoint_Prepend( data->ep, "TestHook4", (void*)TestHook4, "testCode" );
   /* TestHook4, TestHook2, TestHook3 */
   EntryPoint_InsertBefore( data->ep, "TestHook3", "TestHook5", (void*)TestHook5, "testCode" );
   /* TestHook4, TestHook2, TestHook5, TestHook3 */
   EntryPoint_InsertAfter( data->ep, "TestHook4", "TestHook6", (void*)TestHook6, "testCode" );
   /* TestHook4, TestHook6, TestHook2, TestHook5, TestHook3 */
   EntryPoint_Replace( data->ep, "TestHook5", "TestHook7", (void*)TestHook7, "testCode" );
   /* TestHook4, TestHook6, TestHook2, TestHook7, TestHook3 */

   pcu_check_true( data->ep->hooks->count == 5 );
   pcu_check_true( 0 == strcmp( data->ep->hooks->data[0]->name, "TestHook4" ) );
   pcu_check_true( ((Hook*)data->ep->hooks->data[0])->funcPtr == TestHook4 );
   pcu_check_true( 0 == strcmp( data->ep->hooks->data[1]->name, "TestHook6" ) );
   pcu_check_true( ((Hook*)data->ep->hooks->data[1])->funcPtr == TestHook6 );
   pcu_check_true( 0 == strcmp( data->ep->hooks->data[2]->name, "TestHook2" ) );
   pcu_check_true( ((Hook*)data->ep->hooks->data[2])->funcPtr == TestHook2 );
   pcu_check_true( 0 == strcmp( data->ep->hooks->data[3]->name, "TestHook7" ) );
   pcu_check_true( ((Hook*)data->ep->hooks->data[3])->funcPtr == TestHook7 );
   pcu_check_true( 0 == strcmp( data->ep->hooks->data[4]->name, "TestHook3" ) );
   pcu_check_true( ((Hook*)data->ep->hooks->data[4])->funcPtr == TestHook3 );
}


void EntryPointSuite_TestAlwaysFirstLast( EntryPointSuiteData* data ) {
   data->ep = EntryPoint_New( "testEP", EntryPoint_0_CastType );
   EntryPoint_Append_AlwaysLast( data->ep, "TestHook8", (void*)TestHook8, "testCode" );
   /* - TestHook8 */
   EntryPoint_Append( data->ep, "TestHook9", (void*)TestHook9, "testCode" );
   /* TestHook9 - TestHook8 */
   EntryPoint_Prepend_AlwaysFirst( data->ep, "TestHook0", (void*)TestHook0, "testCode" );
   /* TestHook0 - TestHook9 - TestHook8 */
   EntryPoint_Prepend( data->ep, "TestHook1", (void*)TestHook1, "testCode" );
   /* TestHook0 - TestHook1, TestHook9 - TestHook8 */

   pcu_check_true( data->ep->hooks->count == 4 );
   pcu_check_true( 0 == strcmp( data->ep->hooks->data[0]->name, "TestHook0" ) );
   pcu_check_true( ((Hook*)data->ep->hooks->data[0])->funcPtr == TestHook0 );
   pcu_check_true( 0 == strcmp( data->ep->hooks->data[1]->name, "TestHook1" ) );
   pcu_check_true( ((Hook*)data->ep->hooks->data[1])->funcPtr == TestHook1 );
   pcu_check_true( 0 == strcmp( data->ep->hooks->data[2]->name, "TestHook9" ) );
   pcu_check_true( ((Hook*)data->ep->hooks->data[2])->funcPtr == TestHook9 );
   pcu_check_true( 0 == strcmp( data->ep->hooks->data[3]->name, "TestHook8" ) );
   pcu_check_true( ((Hook*)data->ep->hooks->data[3])->funcPtr == TestHook8 );
   pcu_check_true( 0 == strcmp( data->ep->alwaysFirstHook->name, "TestHook0" ) );
   pcu_check_true( ((Hook*)data->ep->alwaysFirstHook)->funcPtr == TestHook0 );
   pcu_check_true( 0 == strcmp( data->ep->alwaysLastHook->name, "TestHook8" ) );
   pcu_check_true( ((Hook*)data->ep->alwaysLastHook)->funcPtr == TestHook8 );
}


void EntryPointSuite_TestReplaceAll( EntryPointSuiteData* data ) {
   data->ep = EntryPoint_New( "testEP", EntryPoint_0_CastType );
   EntryPoint_Append( data->ep, "TestHook0", (void*)TestHook0, "testCode" );
   EntryPoint_Append( data->ep, "TestHook1", (void*)TestHook0, "testCode" );
   /* TestHook0, TestHook1 */
   EntryPoint_ReplaceAll( data->ep, "TestHook2", (void*)TestHook2, "testCode" );
   /* TestHook2 */
   pcu_check_true( data->ep->hooks->count == 1 );
   pcu_check_true( 0 == strcmp( data->ep->hooks->data[0]->name, "TestHook2" ) );
   pcu_check_true( ((Hook*)data->ep->hooks->data[0])->funcPtr == TestHook2 );
}


void EntryPointSuite_TestPurge( EntryPointSuiteData* data ) {
   data->ep = EntryPoint_New( "testEP", EntryPoint_0_CastType );
   EntryPoint_Append( data->ep, "TestHook2", (void*)TestHook2, "testCode" );
   EntryPoint_Append( data->ep, "TestHook3", (void*)TestHook3, "testCode" );
   /* TestHook2, TestHook3 */
   EntryPoint_Purge( data->ep );
   /* */
   pcu_check_true( data->ep->hooks->count == 0 );
}


void EntryPointSuite_TestRun( EntryPointSuiteData* data ) {
   Hook_Index hookIndex;

   data->ep = EntryPoint_New( "testEP", EntryPoint_VoidPtr_CastType );
   EntryPoint_Append( data->ep, "TestHook0", (void*)TestHook0, "testCode" );
   EntryPoint_Append( data->ep, "TestHook1", (void*)TestHook1, "testCode" );
   EntryPoint_Append( data->ep, "TestHook2", (void*)TestHook2, "testCode" );
   EntryPoint_Append( data->ep, "TestHook3", (void*)TestHook3, "testCode" );
   EntryPoint_Append( data->ep, "TestHook4", (void*)TestHook4, "testCode" );
   
   pcu_check_true( data->ep->hooks->count == 5 );

   ((EntryPoint_VoidPtr_CallCast*) data->ep->run)( data->ep, data );

   for (hookIndex = 0; hookIndex < data->ep->hooks->count; hookIndex++ ) {
      pcu_check_true( data->testFuncsRan[hookIndex] == True );
   }
   for (hookIndex = data->ep->hooks->count; hookIndex < NUM_TEST_FUNCS; hookIndex++ ) {
      pcu_check_true( data->testFuncsRan[hookIndex] == False );
   }
}

/* testEP-ClassHook-test with some class hooks involved */

void EntryPointSuite( pcu_suite_t* suite ) {
   pcu_suite_setData( suite, EntryPointSuiteData );
   pcu_suite_setFixtures( suite, EntryPointSuite_Setup, EntryPointSuite_Teardown );
   pcu_suite_addTest( suite, EntryPointSuite_TestRunEmpty );
   pcu_suite_addTest( suite, EntryPointSuite_TestAppendPrepend );
   pcu_suite_addTest( suite, EntryPointSuite_TestInsertBeforeAfterReplace );
   pcu_suite_addTest( suite, EntryPointSuite_TestAlwaysFirstLast );
   pcu_suite_addTest( suite, EntryPointSuite_TestReplaceAll );
   pcu_suite_addTest( suite, EntryPointSuite_TestPurge );
   pcu_suite_addTest( suite, EntryPointSuite_TestRun );
}
