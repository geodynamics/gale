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
** $Id: testIArray.c 2136 2004-09-30 02:47:13Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pcu/pcu.h"
#include "StGermain/Base/Foundation/Foundation.h"
#include "StGermain/Base/IO/IO.h"
#include "StGermain/Base/Container/Container.h"
#include "IArraySuite.h"

#define NUM_ITEMS 5

typedef struct {
   IArray*  iArray;
   int      arrayData[NUM_ITEMS];
} IArraySuiteData;


void IArraySuite_Setup( IArraySuiteData* data ) {
   Index         idx;

   data->iArray = IArray_New();
   for( idx = 0; idx < NUM_ITEMS; idx++ ) {
      data->arrayData[idx] = idx;
   }
}


void IArraySuite_Teardown( IArraySuiteData* data ) {
   NewClass_Delete( data->iArray );
}


void IArraySuite_TestSet( IArraySuiteData* data ) {
   const int* ptr;
   int i_i;

   IArray_Set( data->iArray, NUM_ITEMS, data->arrayData );
   pcu_check_true( IArray_GetSize( data->iArray ) == NUM_ITEMS );
   pcu_check_true( (ptr = IArray_GetPtr( data->iArray )) != 0 );
   for( i_i = 0; i_i < NUM_ITEMS; i_i++ ) {
      pcu_check_true( ptr[i_i] == i_i );
   }
}


void IArraySuite_TestAdd( IArraySuiteData* data ) {
   const int* ptr;
   int i_i;

   IArray_Set( data->iArray, 3, data->arrayData );
   IArray_Add( data->iArray, 2, data->arrayData + 3 );
   pcu_check_true( IArray_GetSize( data->iArray ) == NUM_ITEMS );
   pcu_check_true( (ptr = IArray_GetPtr( data->iArray )) != 0 );
   for( i_i = 0; i_i < NUM_ITEMS; i_i++ ) {
      pcu_check_true( ptr[i_i] == i_i );
   }
}


void IArraySuite_TestRemove( IArraySuiteData* data ) {
   const int*  ptr;
   IMap*       map = NULL;

   IArray_Set( data->iArray, NUM_ITEMS, data->arrayData );
   data->arrayData[0] = 1; data->arrayData[1] = 3;
   map = IMap_New();
   IArray_Remove( data->iArray, 2, data->arrayData, map );
   pcu_check_true( IArray_GetSize( data->iArray ) == 3 );
   pcu_check_true( (ptr = IArray_GetPtr( data->iArray )) != 0 );
   pcu_check_true( ptr[0] == 0 && ptr[1] == 4 && ptr[2] == 2 );
   pcu_check_true( IMap_GetSize( map ) == 1 );
   pcu_check_true( IMap_Has( map, 4 ) );
   pcu_check_true( IMap_Map( map, 4 ) == 1 );

   NewClass_Delete( map );
}


void IArraySuite( pcu_suite_t* suite ) {
   pcu_suite_setData( suite, IArraySuiteData );
   pcu_suite_setFixtures( suite, IArraySuite_Setup, IArraySuite_Teardown );
   pcu_suite_addTest( suite, IArraySuite_TestSet );
   pcu_suite_addTest( suite, IArraySuite_TestAdd );
   pcu_suite_addTest( suite, IArraySuite_TestRemove );
}


