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
#include "CallGraphSuite.h"

typedef struct {
} CallGraphSuiteData;

const Type TestFunc1_Type = "TestFunc1";
void TestFunc1( void ) {}

const Type TestFunc2_Type = "TestFunc2";
void TestFunc2( void ) {}


void CallGraphSuite_Setup( CallGraphSuiteData* data ) {
}

void CallGraphSuite_Teardown( CallGraphSuiteData* data ) {
}
   

void CallGraphSuite_TestPushPopBasics( CallGraphSuiteData* data ) {
   Stg_CallGraph*	cg0;
   Stg_CallGraph 	cg1;
   
   cg0 = Stg_CallGraph_New();
   Stg_CallGraph_Init( &cg1 );
   
   /* Test 2: Can we push the first call on the stack? */
   Stg_CallGraph_Push( cg0, TestFunc1, TestFunc1_Type );
   Stg_CallGraph_Push( &cg1, TestFunc1, TestFunc1_Type );
   pcu_check_true( 
      cg0->_stack && 
      !cg0->_stack->pop && 
      cg0->_stack->entry_I == 0 && 
      cg0->_stack->functionPtr == TestFunc1 &&
      cg0->_tableCount == 1 && 
      cg0->table[0].name == TestFunc1_Type &&
      cg0->table[0].functionPtr == TestFunc1 &&
      cg0->table[0].parentPtr == 0 &&
      cg0->table[0].returned == 0 &&
      cg0->table[0].called == 1 &&
      cg1._stack && 
      !cg1._stack->pop && 
      cg1._stack->entry_I == 0 && 
      cg1._stack->functionPtr == TestFunc1 &&
      cg1._tableCount == 1 && 
      cg1.table[0].name == TestFunc1_Type &&
      cg1.table[0].functionPtr == TestFunc1 &&
      cg1.table[0].parentPtr == 0 &&
      cg1.table[0].returned == 0 &&
      cg1.table[0].called == 1 );
   
   
   /* Test 3: Can we pop the first call on the stack? */
   Stg_CallGraph_Pop( cg0 );
   Stg_CallGraph_Pop( &cg1 );

   pcu_check_true( 
      !cg0->_stack && 
      cg0->table[0].name == TestFunc1_Type &&
      cg0->table[0].functionPtr == TestFunc1 &&
      cg0->table[0].parentPtr == 0 &&
      cg0->table[0].returned == 1 &&
      cg0->table[0].called == 1 &&
      !cg1._stack && 
      cg1.table[0].name == TestFunc1_Type &&
      cg1.table[0].functionPtr == TestFunc1 &&
      cg1.table[0].parentPtr == 0 &&
      cg1.table[0].returned == 1 &&
      cg1.table[0].called == 1 );

   Stg_Class_Delete( &cg1 );
   Stg_Class_Delete(  cg0 );
}


void CallGraphSuite_TestFuncParentNameHandling( CallGraphSuiteData* data ) {
   Stg_CallGraph*	cg2;
   cg2 = Stg_CallGraph_New();
   
   /* Test 4: Ensure that each function pointer - parent pointer - name is a unique entry */
   Stg_CallGraph_Push( cg2, TestFunc1, TestFunc1_Type ); /* p0: f:1 n:1  Should add */
   Stg_CallGraph_Push( cg2, TestFunc2, TestFunc2_Type ); /* p1: f:2 n:2  Should add */
   Stg_CallGraph_Pop( cg2 );
   Stg_CallGraph_Push( cg2, TestFunc2, TestFunc2_Type ); /* p1: f:2 n:2  Should increment */
   Stg_CallGraph_Pop( cg2 );
   Stg_CallGraph_Push( cg2, TestFunc1, TestFunc2_Type ); /* p1: f:1 n:2  Should add (recursive case) */
   Stg_CallGraph_Pop( cg2 );
   Stg_CallGraph_Push( cg2, TestFunc2, TestFunc1_Type ); /* p1: f:2 n:1  Should add */
   Stg_CallGraph_Push( cg2, TestFunc1, TestFunc1_Type ); /* p2: f:1 n:1  Should add */
   Stg_CallGraph_Pop( cg2 );
   Stg_CallGraph_Pop( cg2 );
   Stg_CallGraph_Pop( cg2 );
   pcu_check_true( 
      !cg2->_stack && 
      cg2->_tableCount == 5 && 
      cg2->table[0].returned == 1 &&
      cg2->table[0].called == 1 &&
      cg2->table[1].returned == 2 &&
      cg2->table[1].called == 2 &&
      cg2->table[2].returned == 1 &&
      cg2->table[2].called == 1 &&
      cg2->table[3].returned == 1 &&
      cg2->table[3].called == 1 &&
      cg2->table[4].returned == 1 &&
      cg2->table[4].called == 1 );

   Stg_Class_Delete(  cg2 );
}
   

void CallGraphSuite_TestTableRealloc( CallGraphSuiteData* data ) {
   Stg_CallGraph*	cg3;
   int  		count=0;
   Index		size=0;
   int  		ii=0;

   cg3 = Stg_CallGraph_New();
   /* Test 5: Force a realloc of the table */
   count = (int)((double)1.5 * cg3->_tableSize);
   size = cg3->_tableSize;
   for( ii = 0; ii < count; ii++ ) {
      /* Use "i" as a unique string (unique pointer value)... don't try to print! */
     Stg_CallGraph_Push( cg3, TestFunc1, (Name)(NULL + ii) );
   }
   pcu_check_true( 
      cg3->_tableCount == count && 
      cg3->_tableSize == (size * 2) );
   
   Stg_Class_Delete(  cg3 );
}


void CallGraphSuite_TestCopy( CallGraphSuiteData* data ) {
   Stg_CallGraph*	cg3;
   Stg_CallGraph*	cg3deep;
   int		        count=0;
   int		        ii=0;

   cg3 = Stg_CallGraph_New();
   count = (int)((double)1.5 * cg3->_tableSize);
   for( ii = 0; ii < count; ii++ ) {
      /* Use "i" as a unique string (unique pointer value)... don't try to print! */
     Stg_CallGraph_Push( cg3, TestFunc1, (Name)(NULL + ii) );
   }

   /* Shallow copying not yet implemented */
   cg3deep = Stg_Class_Copy( cg3, 0, True, 0, 0 );
   pcu_check_true(
      cg3->_tableCount == cg3deep->_tableCount &&
      cg3->_tableSize == cg3deep->_tableSize &&
      /* TODO: check not just the table, but the stack too */
      memcmp( cg3->table, cg3deep->table, sizeof(_Stg_CallGraph_Entry) * cg3->_tableCount ) == 0 );
   
   Stg_Class_Delete(  cg3 );
   Stg_Class_Delete(  cg3deep );
}


void CallGraphSuite( pcu_suite_t* suite ) {
   pcu_suite_setData( suite, CallGraphSuiteData );
   pcu_suite_setFixtures( suite, CallGraphSuite_Setup, CallGraphSuite_Teardown );
   pcu_suite_addTest( suite, CallGraphSuite_TestPushPopBasics );
   pcu_suite_addTest( suite, CallGraphSuite_TestFuncParentNameHandling );
   pcu_suite_addTest( suite, CallGraphSuite_TestTableRealloc );
   pcu_suite_addTest( suite, CallGraphSuite_TestCopy );
}
