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
** $Id: testSTree.c 2136 2004-09-30 02:47:13Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pcu/pcu.h"
#include "StGermain/Base/Foundation/Foundation.h"
#include "StGermain/Base/IO/IO.h"
#include "StGermain/Base/Container/Container.h"
#include "STreeMapSuite.h"

typedef struct {
   STreeMap*      sTreeMap;
} STreeMapSuiteData;


void STreeMapSuite_Setup( STreeMapSuiteData* data ) {
   data->sTreeMap = STreeMap_New();
   STreeMap_SetItemSize( data->sTreeMap, sizeof(int), sizeof(int) );
   STree_SetIntCallbacks( data->sTreeMap );
}


void STreeMapSuite_Teardown( STreeMapSuiteData* data ) {
   NewClass_Delete( data->sTreeMap );
}


void STreeMapSuite_TestConstruct( STreeMapSuiteData* data ) {
   pcu_check_true( data->sTreeMap );
}


void STreeMapSuite_TestInsert( STreeMapSuiteData* data ) {
   int c_i;

   for ( c_i = 0; c_i < 10; c_i++ ) {
      int tmp = 10 * c_i;
      STreeMap_Insert( data->sTreeMap, &c_i, &tmp );
   }
   pcu_check_true( STree_GetSize( data->sTreeMap ) == 10 );

}


void STreeMapSuite_TestMap( STreeMapSuiteData* data ) {
   int c_i;

   for ( c_i = 0; c_i < 10; c_i++ ) {
      int tmp = 10 * c_i;
      STreeMap_Insert( data->sTreeMap, &c_i, &tmp );
   }
   for( c_i = 0; c_i < 10; c_i++ ) {
      pcu_check_true( *((int*)STreeMap_Map( data->sTreeMap, &c_i )) == 10 * c_i );
   }

}


void STreeMapSuite( pcu_suite_t* suite ) {
   pcu_suite_setData( suite, STreeMapSuiteData );
   pcu_suite_setFixtures( suite, STreeMapSuite_Setup, STreeMapSuite_Teardown );
   pcu_suite_addTest( suite, STreeMapSuite_TestConstruct );
   pcu_suite_addTest( suite, STreeMapSuite_TestInsert );
   pcu_suite_addTest( suite, STreeMapSuite_TestMap );
}
