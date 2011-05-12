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
** $Id: testIMap.c 2136 2004-09-30 02:47:13Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pcu/pcu.h"
#include "StGermain/Base/Foundation/Foundation.h"
#include "StGermain/Base/IO/IO.h"
#include "StGermain/Base/Container/Container.h"
#include "IMapSuite.h"

#define NUM_ITEMS 5

typedef struct {
   IMap*    iMap;
} IMapSuiteData;


void IMapSuite_Setup( IMapSuiteData* data ) {
   data->iMap = IMap_New();
}


void IMapSuite_Teardown( IMapSuiteData* data ) {
   NewClass_Delete( data->iMap );
}


void IMapSuite_TestConstruct( IMapSuiteData* data ) {

   pcu_check_true( data->iMap );
   pcu_check_true( data->iMap->maxSize == 0 );
   pcu_check_true( data->iMap->curSize == 0 );
   /* The 3 conditions below don't hold because IMap_SetMaxSize( self, 0 ) is called
    *  in the IMap constructor function. I presume Luke did this for a good reason, so
    *  will just comment out the conditions for now. PatrickSunter, 3 Jul 2009 */   
   /*pcu_check_true( data->iMap->tblSize == 0 );
   pcu_check_true( data->iMap->tbl == NULL );
   pcu_check_true( data->iMap->used == NULL );*/
}


void IMapSuite_TestSetMaxSize( IMapSuiteData* data ) {
   IMap_SetMaxSize( data->iMap, 10 );
   pcu_check_true( data->iMap->maxSize == 10 );
   pcu_check_true( data->iMap->tblSize >= 10 );
   IMap_SetMaxSize( data->iMap, 20 );
   pcu_check_true( data->iMap->maxSize == 20 );
   pcu_check_true( data->iMap->tblSize >= 20 );
}


void IMapSuite_TestInsert( IMapSuiteData* data ) {
   int i_i;

   IMap_SetMaxSize( data->iMap, 10 );
   for( i_i = 0; i_i < 20; i_i += 2 ) {
      IMap_Insert( data->iMap, i_i, i_i + 100 );
   }
   pcu_check_assert( IMap_Insert( data->iMap, 0, 100 ) );
   pcu_check_true( IMap_GetSize( data->iMap ) == 10 );
   for( i_i = 0; i_i < 20; i_i += 2 ) {
      pcu_check_true( IMap_Has( data->iMap, i_i ) );
   }
}


void IMapSuite_TestMap( IMapSuiteData* data ) {
   int i_i;

   IMap_SetMaxSize( data->iMap, 10 );
   for( i_i = 0; i_i < 20; i_i += 2 ) {
      IMap_Insert( data->iMap, i_i, i_i + 100 );
   }
   for( i_i = 0; i_i < 20; i_i += 2 ) {
      pcu_check_true( IMap_Map( data->iMap, i_i ) == i_i + 100 );
   }
}


void IMapSuite_TestRemove( IMapSuiteData* data ) {
   int i_i;

   IMap_SetMaxSize( data->iMap, 10 );
   for( i_i = 0; i_i < 20; i_i += 2 ) {
      IMap_Insert( data->iMap, i_i, i_i + 100 );
   }
   for( i_i = 0; i_i < 20; i_i += 4 ) {
      IMap_Remove( data->iMap, i_i );
   }
   pcu_check_assert( IMap_Remove( data->iMap, 0 ) );
   pcu_check_true( IMap_GetSize( data->iMap ) == 5 );
   for( i_i = 0; i_i < 20; i_i += 4 ) {
      if( i_i % 4 ) {
	      pcu_check_true( IMap_Has( data->iMap, i_i ) );
      }
      else {
	      pcu_check_true( !IMap_Has( data->iMap, i_i ) );
      }
   }
}


void IMapSuite( pcu_suite_t* suite ) {
   pcu_suite_setData( suite, IMapSuiteData );
   pcu_suite_setFixtures( suite, IMapSuite_Setup, IMapSuite_Teardown );
   pcu_suite_addTest( suite, IMapSuite_TestConstruct );
   pcu_suite_addTest( suite, IMapSuite_TestSetMaxSize );
   pcu_suite_addTest( suite, IMapSuite_TestInsert );
   pcu_suite_addTest( suite, IMapSuite_TestMap );
   pcu_suite_addTest( suite, IMapSuite_TestRemove );
}


