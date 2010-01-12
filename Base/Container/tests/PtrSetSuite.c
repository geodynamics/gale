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
#include "PtrSetSuite.h"

#define NUM_ITEMS 100

typedef struct {
   PtrSet*     setA;
   PtrSet*     setB;
   unsigned    dataArray[NUM_ITEMS];
   Bool        inSet[NUM_ITEMS];
} PtrSetSuiteData;


void PtrSetSuite_Setup( PtrSetSuiteData* data ) {
   Index         idx;

   data->setA = PtrSet_New( NULL );
   data->setB = PtrSet_New( NULL );

   for( idx=0; idx < NUM_ITEMS; idx++ ) {
      data->inSet[idx] = False;
      /* We deliberately want the actual _numbers_ to be random, as the ptrset should do all its
       * comparisons based on ptr values, not the numbers themselves */
      data->dataArray[idx] = (unsigned)rand();
   }
}


void PtrSetSuite_Teardown( PtrSetSuiteData* data ) {
   Stg_Class_Delete( data->setA );
   Stg_Class_Delete( data->setB );
}


static void markArray( void* setItem, void* args ) {
   PtrSetSuiteData*  data = (PtrSetSuiteData*)args;
   unsigned int      ptrIndex = 0;

   assert( data );

   ptrIndex = ((ArithPointer)setItem - (ArithPointer)data->dataArray) / sizeof(unsigned);
   data->inSet[ptrIndex] = True;
}


void PtrSetSuite_TestInsertTraverse( PtrSetSuiteData* data ) {
   unsigned    ptr_I;

   for( ptr_I = 0; ptr_I < NUM_ITEMS; ptr_I++ ) {
      Set_Insert( data->setA, &data->dataArray[ptr_I] );
   }

   Set_Traverse( data->setA, markArray, data );
   for( ptr_I = 0; ptr_I < NUM_ITEMS; ptr_I++ ) {
      pcu_check_true( data->inSet[ptr_I] == True );
   }
}


void PtrSetSuite_TestUnion( PtrSetSuiteData* data ) {
   Set*        setC=NULL;
   unsigned    ptr_I;

   for( ptr_I = NUM_ITEMS*1/4; ptr_I < NUM_ITEMS*5/8; ptr_I++ ) {
      Set_Insert( data->setA, &data->dataArray[ptr_I] );
   }
   for( ptr_I = NUM_ITEMS*3/8; ptr_I < NUM_ITEMS*3/4; ptr_I++ ) {
      Set_Insert( data->setB, &data->dataArray[ptr_I] );
   }

   setC = Set_Union( data->setA, data->setB );
   Set_Traverse( setC, markArray, data );
   
   for( ptr_I = 0; ptr_I < NUM_ITEMS*1/4; ptr_I++ ) {
      pcu_check_true( data->inSet[ptr_I] == False );
   }
   for( ptr_I = NUM_ITEMS*1/4; ptr_I < NUM_ITEMS*3/4; ptr_I++ ) {
      pcu_check_true( data->inSet[ptr_I] == True );
   }
   for( ptr_I = NUM_ITEMS*3/4; ptr_I < NUM_ITEMS; ptr_I++ ) {
      pcu_check_true( data->inSet[ptr_I] == False );
   }
}


void PtrSetSuite_TestIntersection( PtrSetSuiteData* data ) {
   Set*        setC=NULL;
   unsigned    ptr_I;

   for( ptr_I = NUM_ITEMS*1/4; ptr_I < NUM_ITEMS*5/8; ptr_I++ ) {
      Set_Insert( data->setA, &data->dataArray[ptr_I] );
   }
   for( ptr_I = NUM_ITEMS*3/8; ptr_I < NUM_ITEMS*3/4; ptr_I++ ) {
      Set_Insert( data->setB, &data->dataArray[ptr_I] );
   }

   setC = Set_Intersection( data->setA, data->setB );
   Set_Traverse( setC, markArray, data );
   
   for( ptr_I = 0; ptr_I < NUM_ITEMS*3/8; ptr_I++ ) {
      pcu_check_true( data->inSet[ptr_I] == False );
   }
   for( ptr_I = NUM_ITEMS*3/8; ptr_I < NUM_ITEMS*5/8; ptr_I++ ) {
      pcu_check_true( data->inSet[ptr_I] == True );
   }
   for( ptr_I = NUM_ITEMS*5/8; ptr_I < NUM_ITEMS; ptr_I++ ) {
      pcu_check_true( data->inSet[ptr_I] == False );
   }
}


void PtrSetSuite_TestSubtraction( PtrSetSuiteData* data ) {
   Set*        setC=NULL;
   unsigned    ptr_I;

   for( ptr_I = NUM_ITEMS*1/4; ptr_I < NUM_ITEMS*5/8; ptr_I++ ) {
      Set_Insert( data->setA, &data->dataArray[ptr_I] );
   }
   for( ptr_I = NUM_ITEMS*3/8; ptr_I < NUM_ITEMS*3/4; ptr_I++ ) {
      Set_Insert( data->setB, &data->dataArray[ptr_I] );
   }

   setC = Set_Subtraction( data->setA, data->setB );
   Set_Traverse( setC, markArray, data );
   
   for( ptr_I = 0; ptr_I < NUM_ITEMS*1/4; ptr_I++ ) {
      pcu_check_true( data->inSet[ptr_I] == False );
   }
   for( ptr_I = NUM_ITEMS*1/4; ptr_I < NUM_ITEMS*3/8; ptr_I++ ) {
      pcu_check_true( data->inSet[ptr_I] == True );
   }
   for( ptr_I = NUM_ITEMS*3/8; ptr_I < NUM_ITEMS; ptr_I++ ) {
      pcu_check_true( data->inSet[ptr_I] == False );
   }
}


void PtrSetSuite_TestPerformance( PtrSetSuiteData* data ) {
   unsigned    ptr_I;
   double      startTime=0;
   double      stopTime=0;
   double      timeSpent=0;

   /* Raq: I expect the set to be able to insert 100,000 items in the worst case
      scenario in a reasonable amount of time. */
   startTime = MPI_Wtime();
   for( ptr_I = 0; ptr_I < 100000; ptr_I++ ) {
      Set_Insert( data->setA, &data->dataArray[ptr_I] );
   }
   stopTime = MPI_Wtime();
   timeSpent = stopTime - startTime;

   /* 5 seconds is arbitrary, take into account slower systems */
   pcu_check_true( timeSpent < 5.0 );
}


void PtrSetSuite( pcu_suite_t* suite ) {
   pcu_suite_setData( suite, PtrSetSuiteData );
   pcu_suite_setFixtures( suite, PtrSetSuite_Setup, PtrSetSuite_Teardown );
   pcu_suite_addTest( suite, PtrSetSuite_TestInsertTraverse );
   pcu_suite_addTest( suite, PtrSetSuite_TestUnion );
   pcu_suite_addTest( suite, PtrSetSuite_TestIntersection );
   pcu_suite_addTest( suite, PtrSetSuite_TestSubtraction );
   pcu_suite_addTest( suite, PtrSetSuite_TestPerformance );
}


