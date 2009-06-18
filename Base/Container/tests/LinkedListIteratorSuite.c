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
** $Id: testLinkedList.c 3462 2006-02-19 06:53:24Z WalterLandry $
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
#include "LinkedListSuite.h"
#include "LinkedListIteratorSuite.h"

#define NUM_DATA 100

typedef struct {
   LinkedList*       numList;
   int*              array[NUM_DATA];
} LinkedListIteratorSuiteData;


void LinkedListIteratorSuite_Setup( LinkedListIteratorSuiteData* data ) {
   Index          ii = 0;

   data->numList = LinkedList_New(
            LinkedListSuite_CompareFunction,
            LinkedListSuite_DataCopyFunction,
            LinkedListSuite_DataPrintFunction,
            NULL,
            LINKEDLIST_UNSORTED);

   for(ii=0; ii<NUM_DATA; ii++){
      data->array[ii] = Memory_Alloc(int, "TestLinkedList_ArrayEntry");
      *data->array[ii] = ii;
   }
}

void LinkedListIteratorSuite_Teardown( LinkedListIteratorSuiteData* data ) {
   Index          ii = 0;

   Stg_Class_Delete( data->numList );
   for(ii=0; ii < NUM_DATA; ii++){
      Memory_Free( data->array[ii] );
   }
}


void LinkedListIteratorSuite_TestRetreive( LinkedListIteratorSuiteData* data ) {
	LinkedListIterator *iterator = NULL;
	LinkedListIterator *iterator1 = NULL;
	void *result = NULL, *result1 = NULL;
	int ii = 0;
	int jj = 0;

   iterator = LinkedListIterator_New( data->numList );
   iterator1 = LinkedListIterator_New( data->numList );

   for(ii=0; ii<NUM_DATA; ii++){
      LinkedList_InsertNode( data->numList, data->array[ii], sizeof(int));
   }

   ii = 0;
   for( result = LinkedListIterator_First( iterator ); result; result = LinkedListIterator_Next( iterator ) ){
      pcu_check_true( *((int*)result) == (NUM_DATA-1 - *data->array[ii]) );
      jj=0;
      for( result1 = LinkedListIterator_First( iterator1 ); result1; result1 = LinkedListIterator_Next( iterator1 ) ){
         pcu_check_true( *((int*)result1) == (NUM_DATA-1 - *data->array[jj]) );
         jj++;
      }
      pcu_check_true( jj == NUM_DATA );
      result = LinkedListIterator_Next( iterator1 ); 
      pcu_check_true( NULL == result );
      ii++;
   }
   pcu_check_true( ii == NUM_DATA );
   result = LinkedListIterator_Next( iterator );
   pcu_check_true( NULL == result );
}


void LinkedListIteratorSuite( pcu_suite_t* suite ) {
   pcu_suite_setData( suite, LinkedListIteratorSuiteData );
   pcu_suite_setFixtures( suite, LinkedListIteratorSuite_Setup, LinkedListIteratorSuite_Teardown );
   pcu_suite_addTest( suite, LinkedListIteratorSuite_TestRetreive );
}
