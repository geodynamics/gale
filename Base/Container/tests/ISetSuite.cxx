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
** $Id: testISet.c 2136 2004-09-30 02:47:13Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pcu/pcu.h"
#include "StGermain/Base/Foundation/Foundation.h"
#include "StGermain/Base/IO/IO.h"
#include "StGermain/Base/Container/Container.h"
#include "ISetSuite.h"

typedef struct {
   ISet*    iSet0;
   ISet*    iSet1;
} ISetSuiteData;


void ISetSuite_Setup( ISetSuiteData* data ) {
   data->iSet0 = ISet_New();
   data->iSet1 = ISet_New();
}


void ISetSuite_Teardown( ISetSuiteData* data ) {
   NewClass_Delete( data->iSet0 );
   NewClass_Delete( data->iSet1 );
}


void ISetSuite_TestConstruct( ISetSuiteData* data ) {
   pcu_check_true( data->iSet0 );
   pcu_check_true( data->iSet0->maxSize == 0 );
   pcu_check_true( data->iSet0->curSize == 0 );
   /* Same situation as IMapSuite_Construct - PatrickSunter, 3 June 2009 */
#if 0
   pcu_check_true( data->iSet0->tblSize == 0 );
   pcu_check_true( data->iSet0->tbl == NULL );
   pcu_check_true( data->iSet0->used == NULL );
#endif
}


void ISetSuite_TestInsert( ISetSuiteData* data ) {
   int i_i;

   ISet_SetMaxSize( data->iSet0, 20 );
   for( i_i = 0; i_i < 20; i_i += 2 ) {
      ISet_Insert( data->iSet0, i_i );
   }
   for( i_i = 1; i_i < 20; i_i += 2 ) {
      ISet_Insert( data->iSet0, i_i );
   }
   pcu_check_true( ISet_GetSize( data->iSet0 ) == 20 );
   pcu_check_assert( ISet_Insert( data->iSet0, 0 ) );
   ISet_TryInsert( data->iSet0, 0 );
   for( i_i = 0; i_i < 20; i_i++ ) {
      pcu_check_true( ISet_Has( data->iSet0, i_i ) );
   }
   pcu_check_true( !ISet_Has( data->iSet0, 21 ) );
}


void ISetSuite_TestUseArray( ISetSuiteData* data ) {
   int array[10] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
   int i_i;

   ISet_UseArray( data->iSet0, 10, array );
   pcu_check_true( ISet_GetSize( data->iSet0 ) == 10 );
   for( i_i = 0; i_i < 10; i_i++ ) {
      pcu_check_true( ISet_Has( data->iSet0, i_i ) );
   }
   pcu_check_true( !ISet_Has( data->iSet0, 10 ) );
}


void ISetSuite_TestUnion( ISetSuiteData* data ) {
   int array0[5] = {0, 1, 2, 3, 4};
   int array1[5] = {3, 4, 5, 6, 7};
   int i_i;

   ISet_UseArray( data->iSet0, 5, array0 );
   ISet_UseArray( data->iSet1, 5, array1 );
   ISet_Union( data->iSet0, data->iSet1 );
   pcu_check_true( ISet_GetSize( data->iSet0 ) == 8 );
   for( i_i = 0; i_i < 8; i_i++ ) {
      pcu_check_true( ISet_Has( data->iSet0, i_i ) );
   }
   pcu_check_true( !ISet_Has( data->iSet0, 10 ) );
}


void ISetSuite_TestIsect( ISetSuiteData* data ) {
   int array0[5] = {0, 1, 2, 3, 4};
   int array1[5] = {3, 4, 5, 6, 7};
   int i_i;

   ISet_UseArray( data->iSet0, 5, array0 );
   ISet_UseArray( data->iSet1, 5, array1 );
   ISet_Isect( data->iSet0, data->iSet1 );
   pcu_check_true( ISet_GetSize( data->iSet0 ) == 2 );
   for( i_i = 3; i_i < 5; i_i++ ) {
      pcu_check_true( ISet_Has( data->iSet0, i_i ) );
   }
   pcu_check_true( !ISet_Has( data->iSet0, 1 ) );
}


void ISetSuite_TestSubtr( ISetSuiteData* data ) {
   int array0[5] = {0, 1, 2, 3, 4};
   int array1[5] = {3, 4, 5, 6, 7};
   int i_i;

   ISet_UseArray( data->iSet0, 5, array0 );
   ISet_UseArray( data->iSet1, 5, array1 );
   ISet_Subtr( data->iSet0, data->iSet1 );
   pcu_check_true( ISet_GetSize( data->iSet0 ) == 3 );
   for( i_i = 0; i_i < 3; i_i++ ) {
      pcu_check_true( ISet_Has( data->iSet0, i_i ) );
   }
   pcu_check_true( !ISet_Has( data->iSet0, 3 ) );
}


void ISetSuite( pcu_suite_t* suite ) {
   pcu_suite_setData( suite, ISetSuiteData );
   pcu_suite_setFixtures( suite, ISetSuite_Setup, ISetSuite_Teardown );
   pcu_suite_addTest( suite, ISetSuite_TestConstruct );
   pcu_suite_addTest( suite, ISetSuite_TestInsert );
   pcu_suite_addTest( suite, ISetSuite_TestUseArray );
   pcu_suite_addTest( suite, ISetSuite_TestUnion );
   pcu_suite_addTest( suite, ISetSuite_TestIsect );
   pcu_suite_addTest( suite, ISetSuite_TestSubtr );
}


