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
** $Id: testUIntMap.c 3462 2006-02-19 06:53:24Z WalterLandry $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pcu/pcu.h"
#include "StGermain/Base/Foundation/Foundation.h"
#include "StGermain/Base/IO/IO.h"
#include "StGermain/Base/Container/Container.h"
#include "UIntMapSuite.h"

typedef struct {
   UIntMap*    map;
   unsigned	   size;
} UIntMapSuiteData;


void UIntMapSuite_Setup( UIntMapSuiteData* data ) {
   data->map = UIntMap_New();
   data->size = 5;
}


void UIntMapSuite_Teardown( UIntMapSuiteData* data ) {
   FreeObject( data->map );
}


void UIntMapSuite_FillMap( UIntMapSuiteData* data ) {
	unsigned	i;

	for( i = 0; i < data->size; i++ ) {
		UIntMap_Insert( data->map, i, data->size + i );
   }
}


void UIntMapSuite_TestInsert( UIntMapSuiteData* data ) {
   UIntMapSuite_FillMap( data );
   pcu_check_true( data->map->size == data->size );
}


void UIntMapSuite_TestMap( UIntMapSuiteData* data ) {
   unsigned	val;
   unsigned	i;

   UIntMapSuite_FillMap( data );

   for( i = 0; i < data->size; i++ ) {
      pcu_check_true( UIntMap_Map( data->map, i, &val ) );
      pcu_check_true( val == data->size + i );
   }
}


void UIntMapSuite_TestMemoryMap( UIntMapSuiteData* data ) {
   unsigned	nReps = 10;
   unsigned	val;
   unsigned	r_i;

   for( r_i = 0; r_i < nReps; r_i++ ) {
      unsigned	i;

      UIntMapSuite_FillMap( data );
      for( i = 0; i < data->size; i++ ) {
         pcu_check_true( UIntMap_Map( data->map, i, &val ) );
         pcu_check_true( val == data->size + i );
      }
      FreeObject( data->map );
      data->map = UIntMap_New();
   }
}


void UIntMapSuite( pcu_suite_t* suite ) {
   pcu_suite_setData( suite, UIntMapSuiteData );
   pcu_suite_setFixtures( suite, UIntMapSuite_Setup, UIntMapSuite_Teardown );
   pcu_suite_addTest( suite, UIntMapSuite_TestInsert );
   pcu_suite_addTest( suite, UIntMapSuite_TestMap );
   pcu_suite_addTest( suite, UIntMapSuite_TestMemoryMap );
}
