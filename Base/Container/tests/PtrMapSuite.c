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
** $Id: testPtrMap.c 3462 2006-02-19 06:53:24Z WalterLandry $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pcu/pcu.h"
#include "StGermain/Base/Foundation/Foundation.h"
#include "StGermain/Base/IO/IO.h"
#include "StGermain/Base/Container/Container.h"
#include "PtrMapSuite.h"

typedef struct {
   PtrMap*    map;
} PtrMapSuiteData;


void PtrMapSuite_Setup( PtrMapSuiteData* data ) {
   data->map = PtrMap_New( 10 );
}


void PtrMapSuite_Teardown( PtrMapSuiteData* data ) {
   Stg_Class_Delete( data->map );
}


void PtrMapSuite_TestAppendFind( PtrMapSuiteData* data ) {
   ArithPointer		idx;
   
   for( idx = 0; idx < 100; idx++ ) {
      PtrMap_Append( data->map, (void*)(idx + 1), (void*)(100 - idx) );
   }
   
   for( idx = 0; idx < 100; idx++ ) {
      pcu_check_true( (ArithPointer)PtrMap_Find( data->map, (void*)(idx + 1) )
         == (100 - idx) );
   }
}


void PtrMapSuite( pcu_suite_t* suite ) {
   pcu_suite_setData( suite, PtrMapSuiteData );
   pcu_suite_setFixtures( suite, PtrMapSuite_Setup, PtrMapSuite_Teardown );
   pcu_suite_addTest( suite, PtrMapSuite_TestAppendFind );
}
