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
** $Id: testIndexMap.c 3462 2006-02-19 06:53:24Z WalterLandry $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pcu/pcu.h"
#include "StGermain/Base/Foundation/Foundation.h"
#include "StGermain/Base/IO/IO.h"
#include "StGermain/Base/Container/Container.h"
#include "IndexMapSuite.h"

typedef struct {
   IndexMap*      map;
} IndexMapSuiteData;


void IndexMapSuite_Setup( IndexMapSuiteData* data ) {
   data->map = IndexMap_New();
}

void IndexMapSuite_Teardown( IndexMapSuiteData* data ) {
   Stg_Class_Delete( data->map );
}


void IndexMapSuite_TestAppendFind( IndexMapSuiteData* data ) {
   Index          idx;

   pcu_docstring( "This test inserts a set of indices to an IndexSet in reverse order, "
      " then checks they can be found at the correct indices." );
   
   for( idx = 0; idx < 100; idx++ ) {
      IndexMap_Append( data->map, idx + 1, 100 - idx );
   }
   
   for( idx = 0; idx < 100; idx++ ) {
      pcu_check_true( IndexMap_Find( data->map, idx + 1 ) == (100 - idx) );
   }
}


void IndexMapSuite( pcu_suite_t* suite ) {
   pcu_suite_setData( suite, IndexMapSuiteData );
   pcu_suite_setFixtures( suite, IndexMapSuite_Setup, IndexMapSuite_Teardown );
   pcu_suite_addTest( suite, IndexMapSuite_TestAppendFind );
}
