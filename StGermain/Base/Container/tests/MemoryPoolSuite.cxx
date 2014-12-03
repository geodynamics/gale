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
** $Id: testMemoryPool.c 3462 2006-02-19 06:53:24Z WalterLandry $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pcu/pcu.h"
#include "StGermain/Base/Foundation/Foundation.h"
#include "StGermain/Base/IO/IO.h"
#include "StGermain/Base/Container/Container.h"
#include "MemoryPoolSuite.h"

typedef struct Plane_t{
	double normal[3];
	double k;
} Plane;

#define CACHE_SIZE 200

typedef struct {
   MemoryPool*    pool;
   Plane*         planeRefs[CACHE_SIZE*500];
} MemoryPoolSuiteData;



void MemoryPoolSuite_Setup( MemoryPoolSuiteData* data ) {
   Index    ii=0;

   data->pool = MemoryPool_New( Plane, CACHE_SIZE, CACHE_SIZE/2 );
   for( ii=0; ii<CACHE_SIZE; ii++ ){
      data->planeRefs[ii] = NULL;
   }
}


void MemoryPoolSuite_Teardown( MemoryPoolSuiteData* data ) {
   Stg_Class_Delete( data->pool );
}


void MemoryPoolSuite_TestAllocation( MemoryPoolSuiteData* data ) {
   Plane*      p = NULL;
   int         i = 0;
   Bool        passed = False;
   int         objCounter = 0;

   passed = True;
   /* Testing memory allocation from the Memory Pool..  */
   for( i=0; i<CACHE_SIZE; i++ ){
      p = NULL;
      p = MemoryPool_NewObject( Plane, data->pool );
      if( !p ){
         passed = False;
         break;
      }
      else{
         objCounter++;
         data->planeRefs[i] = p;
         pcu_check_true( data->pool->numElementsFree == CACHE_SIZE-(i+1) );
      }
   }
   pcu_check_true( passed );

   pcu_check_true( data->pool->numInitialElements == CACHE_SIZE );
   pcu_check_true( data->pool->numElements == CACHE_SIZE );
   pcu_check_true( data->pool->numElementsFree == 0 );
   pcu_check_true( data->pool->numMemChunks == 1 );

}
   

void MemoryPoolSuite_TestOverAllocation( MemoryPoolSuiteData* data ) {
   Plane*      p = NULL;
   int         i = 0;
   Bool        passed = False;
   int         objCounter = 0;

   passed = True;
   for( i=0; i<CACHE_SIZE*100; i++ ){
      p = MemoryPool_NewObject( Plane, data->pool );
      data->planeRefs[objCounter++] = p;
      if( !p ){
         passed = False;
         break;
      }
   }
   pcu_check_true( passed );
   pcu_check_true( data->pool->numInitialElements == CACHE_SIZE );
   pcu_check_true( data->pool->numElements == CACHE_SIZE*100 );
   pcu_check_true( data->pool->numElementsFree == 0 );
   pcu_check_true( data->pool->numMemChunks == CACHE_SIZE*100/(CACHE_SIZE/2)-1 );
}


void MemoryPoolSuite_TestDeallocation( MemoryPoolSuiteData* data ) {
   Plane*      p = NULL;
   int         i = 0;
   Bool        passed = False;
   int         objCounter = 0;

   passed = True;

   for( i=0; i<CACHE_SIZE; i++ ){
      p = NULL;
      p = MemoryPool_NewObject( Plane, data->pool );
      objCounter++;
      data->planeRefs[i] = p;
   }

   for( i=0; i<CACHE_SIZE; i++ ){
      if(False == MemoryPool_DeleteObject( data->pool, data->planeRefs[i] )){
         passed = False;
         break;
      }
   }
   pcu_check_true( passed );

   pcu_check_true( data->pool->numInitialElements == CACHE_SIZE );
   pcu_check_true( data->pool->numElements == 0 );
   pcu_check_true( data->pool->numElementsFree == 0 );
   pcu_check_true( data->pool->numMemChunks == 0 );
}
   

void MemoryPoolSuite_TestIllegalDeallocation( MemoryPoolSuiteData* data ) {
   Plane*      p = NULL;
   int         i = 0;
   int         objCounter = 0;
   int*        junkRefs[CACHE_SIZE];
   int         testData[CACHE_SIZE];

   for( i=0; i<CACHE_SIZE; i++ ){
      p = NULL;
      p = MemoryPool_NewObject( Plane, data->pool );
      objCounter++;
      data->planeRefs[i] = p;
   }

   for( i=0; i<CACHE_SIZE; i++ ){
      junkRefs[i] = &testData[i];
   }

   for( i=0; i<CACHE_SIZE/4; i++ ){
      junkRefs[i] = (int*)(junkRefs+i+1);
      pcu_check_true(False == MemoryPool_DeleteObject( data->pool, junkRefs[i] ));
   }

   pcu_check_true( data->pool->numInitialElements == CACHE_SIZE );
   pcu_check_true( data->pool->numElements == CACHE_SIZE );
   pcu_check_true( data->pool->numElementsFree == 0 );
   pcu_check_true( data->pool->numMemChunks == 1 );
}


void MemoryPoolSuite( pcu_suite_t* suite ) {
   pcu_suite_setData( suite, MemoryPoolSuiteData );
   pcu_suite_setFixtures( suite, MemoryPoolSuite_Setup, MemoryPoolSuite_Teardown );
   pcu_suite_addTest( suite, MemoryPoolSuite_TestAllocation );
   pcu_suite_addTest( suite, MemoryPoolSuite_TestOverAllocation );
   pcu_suite_addTest( suite, MemoryPoolSuite_TestDeallocation );
   pcu_suite_addTest( suite, MemoryPoolSuite_TestIllegalDeallocation );
}


