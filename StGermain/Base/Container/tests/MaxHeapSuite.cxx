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
** $Id: testMaxHeap.c 3462 2006-02-19 06:53:24Z WalterLandry $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "pcu/pcu.h"
#include "StGermain/Base/Foundation/Foundation.h"
#include "StGermain/Base/IO/IO.h"
#include "StGermain/Base/Container/Container.h"
#include "MaxHeapSuite.h"

typedef struct {
   MaxHeap*       heap;
   int*           dataArray;
   int**          keys;
} MaxHeapSuiteData;

#define NUM_DATA 10
#define NUM_INITIAL_DATA 5

int compareFunction(void *data1, void *data2)
{
   int *d1 = NULL, *d2 = NULL;

   d1 = (int*)data1;
   d2 = (int*)data2;

   if (d1 == NULL || d2 == NULL){
      return 0;   
   }
   
   if (*d1 > *d2){
      return  1;
   }
   else if (*d1 == *d2){
      return 0;
   }
   else{
      return -1;
   }
}

void** extendArray( int newCount, void ***array )
{
   assert( array );
   /* TODO Check if this is executed by creating memory */
   (*(int***)array) = (int**)Memory_Realloc_Array((*(int***)array), int**, newCount );
   if( ((*(int***)array) == NULL) ){
      Journal_Firewall( 0, Journal_Register( ErrorStream_Type, (Name)"testMaxHeap" ), "Memory allocation failed in '%s'!!\n Aborting..!!\n", __func__ );
   
   }
   else{
      return *array;
   }

   return NULL;
}

void keySwap( void **a, void **b )
{
   int *temp;

   temp = (*((int**)a));

   (*((int**)a)) = (*((int**)b));
   (*((int**)b)) = temp;
}


void MaxHeapSuite_Setup( MaxHeapSuiteData* data ) {
   Index    ii=0;

   data->dataArray = Memory_Alloc_Array_Unnamed( int, NUM_DATA );
   data->keys = Memory_Alloc_Array_Unnamed( int*, NUM_INITIAL_DATA );
   
   for(ii=0; ii<NUM_INITIAL_DATA; ii++){
      data->keys[ii] = &(data->dataArray[ii]);
   }
   for(ii=0; ii<NUM_DATA; ii++){
      data->dataArray[ii] = ii;
   }
   
   data->heap = MaxHeap_New(
            (void**)(data->keys), sizeof(int),
            NUM_INITIAL_DATA,
            keySwap,
            compareFunction,
            extendArray );
}


void MaxHeapSuite_Teardown( MaxHeapSuiteData* data ) {
   Stg_Class_Delete( data->heap );
   Memory_Free( data->dataArray );
   /* Note: _Heap_Delete() (Heap.c:144) already frees the keys array. Not sure this is entirely logical - needs to
    *  be well doco'd at least */
   /*Memory_Free( data->keys );*/
}


void MaxHeapSuite_TestCreationExtraction( MaxHeapSuiteData* data ) {
   Index    ii=0;

   /* These initial totals due to set up above */
   pcu_check_true( data->heap->numHeapElements == NUM_INITIAL_DATA );
   pcu_check_true( data->heap->numArrayElements == NUM_INITIAL_DATA );

   for( ii=0; ii<NUM_INITIAL_DATA; ii++ ){
      /* Since we are always extracting the max, expect the order to be reversed */
      pcu_check_true( *(int*)MaxHeap_Extract( data->heap ) == data->dataArray[(NUM_INITIAL_DATA-1)-ii] );
   }

   pcu_check_true( data->heap->numHeapElements == 0 );
   pcu_check_true( data->heap->numArrayElements == NUM_INITIAL_DATA );
}


void MaxHeapSuite_TestInsertionExtraction( MaxHeapSuiteData* data ) {
   Index    ii=0;

   /*Inserting more entries into the Heap*/
   for( ii=NUM_INITIAL_DATA; ii<NUM_DATA; ii++ ){
      MaxHeap_Insert( data->heap, &(data->dataArray[ii]) );
   }
   
   pcu_check_true( data->heap->numHeapElements == NUM_DATA );
   pcu_check_true( data->heap->numArrayElements == NUM_DATA );

   for( ii=0; ii<NUM_DATA; ii++ ){
      /* Since we are always extracting the max, expect the order to be reversed */
      pcu_check_true( *(int*)MaxHeap_Extract( data->heap ) == data->dataArray[(NUM_DATA-1)-ii] );
   }

   pcu_check_true( data->heap->numHeapElements == 0 );
   pcu_check_true( data->heap->numArrayElements == NUM_DATA );
}


void MaxHeapSuite( pcu_suite_t* suite ) {
   pcu_suite_setData( suite, MaxHeapSuiteData );
   pcu_suite_setFixtures( suite, MaxHeapSuite_Setup, MaxHeapSuite_Teardown );
   pcu_suite_addTest( suite, MaxHeapSuite_TestCreationExtraction );
   pcu_suite_addTest( suite, MaxHeapSuite_TestInsertionExtraction );
}


