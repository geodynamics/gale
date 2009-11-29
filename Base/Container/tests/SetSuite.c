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
** $Id: testSet.c 2136 2004-09-30 02:47:13Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#include "pcu/pcu.h"
#include "StGermain/Base/Foundation/Foundation.h"
#include "StGermain/Base/IO/IO.h"
#include "StGermain/Base/Container/Container.h"
#include "SetSuite.h"

#define NUM_ITEMS 100

typedef struct {
   Set*    setA;
   Set*    setB;
   Bool    inSet[NUM_ITEMS];
} SetSuiteData;


static int compareFunc( void* left, void* right ) {
   if( *(int*)left < *(int*)right ) {
      return -1;
   }
   else if( *(int*)left > *(int*)right ) {
      return 1;
   }
   else {
      return 0;
   }
}


static void copyFunc( void** newData, void* data, SizeT size ) {
   *newData = Memory_Alloc_Bytes_Unnamed( size, "char" );
   assert( *newData );

   *(int*)(*newData) = *(int*)data;
}


static void deleteFunc( void* data ) {
   assert( data );
   Memory_Free( data );
}


void SetSuite_Setup( SetSuiteData* data ) {
   Index         idx;

   data->setA = Set_New( NULL, int, compareFunc, copyFunc, deleteFunc );
   data->setB = Set_New( NULL, int, compareFunc, copyFunc, deleteFunc );

   for( idx=0; idx < NUM_ITEMS; idx++ ) {
      data->inSet[idx] = False;
   }
}


void SetSuite_Teardown( SetSuiteData* data ) {
   Stg_Class_Delete( data->setA );
   Stg_Class_Delete( data->setB );
}


static void markArray( void* setItem, void* args ) {
   SetSuiteData*  data = (SetSuiteData*)args;
   assert( data );

   data->inSet[*(int*)setItem] = True;
}


void SetSuite_TestInsertTraverse( SetSuiteData* data ) {
   unsigned   int_I;

   pcu_docstring( "This test inserts a set of ints into a Set, then tests that a simple function "
      "when passed to Set_Traverse() is executed properly." );

   for( int_I = 0; int_I < NUM_ITEMS; int_I++ ) {
      Set_Insert( data->setA, &int_I );
   }

   Set_Traverse( data->setA, markArray, data );
   for( int_I = 0; int_I < NUM_ITEMS; int_I++ ) {
      pcu_check_true( data->inSet[int_I] == True );
   }
}


void SetSuite_TestUnion( SetSuiteData* data ) {
   Set*        setC=NULL;
   unsigned    int_I;

   pcu_docstring( "Checks that the Union of two overlapping sets results in a new set that "
      "contains correct Union of set elements." );

   for( int_I = NUM_ITEMS*1/4; int_I < NUM_ITEMS*5/8; int_I++ ) {
      Set_Insert( data->setA, &int_I );
   }
   for( int_I = NUM_ITEMS*3/8; int_I < NUM_ITEMS*3/4; int_I++ ) {
      Set_Insert( data->setB, &int_I );
   }

   setC = Set_Union( data->setA, data->setB );
   Set_Traverse( setC, markArray, data );
   
   for( int_I = 0; int_I < NUM_ITEMS*1/4; int_I++ ) {
      pcu_check_true( data->inSet[int_I] == False );
   }
   for( int_I = NUM_ITEMS*1/4; int_I < NUM_ITEMS*3/4; int_I++ ) {
      pcu_check_true( data->inSet[int_I] == True );
   }
   for( int_I = NUM_ITEMS*3/4; int_I < NUM_ITEMS; int_I++ ) {
      pcu_check_true( data->inSet[int_I] == False );
   }
}


void SetSuite_TestIntersection( SetSuiteData* data ) {
   Set*        setC=NULL;
   unsigned    int_I;

   for( int_I = NUM_ITEMS*1/4; int_I < NUM_ITEMS*5/8; int_I++ ) {
      Set_Insert( data->setA, &int_I );
   }
   for( int_I = NUM_ITEMS*3/8; int_I < NUM_ITEMS*3/4; int_I++ ) {
      Set_Insert( data->setB, &int_I );
   }

   setC = Set_Intersection( data->setA, data->setB );
   Set_Traverse( setC, markArray, data );
   
   for( int_I = 0; int_I < NUM_ITEMS*3/8; int_I++ ) {
      pcu_check_true( data->inSet[int_I] == False );
   }
   for( int_I = NUM_ITEMS*3/8; int_I < NUM_ITEMS*5/8; int_I++ ) {
      pcu_check_true( data->inSet[int_I] == True );
   }
   for( int_I = NUM_ITEMS*5/8; int_I < NUM_ITEMS; int_I++ ) {
      pcu_check_true( data->inSet[int_I] == False );
   }
}


void SetSuite_TestSubtraction( SetSuiteData* data ) {
   Set*        setC=NULL;
   unsigned    int_I;

   for( int_I = NUM_ITEMS*1/4; int_I < NUM_ITEMS*5/8; int_I++ ) {
      Set_Insert( data->setA, &int_I );
   }
   for( int_I = NUM_ITEMS*3/8; int_I < NUM_ITEMS*3/4; int_I++ ) {
      Set_Insert( data->setB, &int_I );
   }

   setC = Set_Subtraction( data->setA, data->setB );
   Set_Traverse( setC, markArray, data );
   
   for( int_I = 0; int_I < NUM_ITEMS*1/4; int_I++ ) {
      pcu_check_true( data->inSet[int_I] == False );
   }
   for( int_I = NUM_ITEMS*1/4; int_I < NUM_ITEMS*3/8; int_I++ ) {
      pcu_check_true( data->inSet[int_I] == True );
   }
   for( int_I = NUM_ITEMS*3/8; int_I < NUM_ITEMS; int_I++ ) {
      pcu_check_true( data->inSet[int_I] == False );
   }
}


void SetSuite_TestPerformance( SetSuiteData* data ) {
   unsigned    int_I;
   double      startTime=0;
   double      stopTime=0;
   double      timeSpent=0;

   /* Raq: I expect the set to be able to insert 100,000 items in the worst case
      scenario in a reasonable amount of time. */
   startTime = MPI_Wtime();
   for( int_I = 0; int_I < 100000; int_I++ ) {
      Set_Insert( data->setA, &int_I );
   }
   stopTime = MPI_Wtime();
   timeSpent = stopTime - startTime;

   /* 5 seconds is arbitrary, take into account slower systems */
   pcu_check_true( timeSpent < 5.0 );
}


void SetSuite( pcu_suite_t* suite ) {
   pcu_suite_setData( suite, SetSuiteData );
   pcu_suite_setFixtures( suite, SetSuite_Setup, SetSuite_Teardown );
   pcu_suite_addTest( suite, SetSuite_TestInsertTraverse );
   pcu_suite_addTest( suite, SetSuite_TestUnion );
   pcu_suite_addTest( suite, SetSuite_TestIntersection );
   pcu_suite_addTest( suite, SetSuite_TestSubtraction );
   pcu_suite_addTest( suite, SetSuite_TestPerformance );
}


