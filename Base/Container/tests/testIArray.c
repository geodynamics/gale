/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003, Victorian Partnership for Advanced Computing (VPAC) Ltd, 
110 Victoria Street, Melbourne, 3053, Australia.
**
** Authors:
**	Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**	Patrick D. Sunter, Software Engineer, VPAC. (pds@vpac.org)
**	Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**	Siew-Ching Tan, Software Engineer, VPAC. (siew@vpac.org)
**	Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
**	Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
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
** $Id: IMap.c 3952 2007-01-09 06:24:06Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <mpi.h>
#include "Base/Foundation/Foundation.h"
#include "Base/IO/IO.h"
#include "Base/Container/Container.h"

#include "Base/Foundation/TestBegin.h"


void testSetup( int* argc, char** argv[] ) {
   BaseFoundation_Init( argc, argv );
   BaseIO_Init( argc, argv );
   BaseContainer_Init( argc, argv );
}

void testTeardown() {
   BaseContainer_Finalise();
   BaseIO_Finalise();
   BaseFoundation_Finalise();
}

TestBegin( Construct ) {
   IArray* arr;

   TestNoAssert( arr = IArray_New() );
   TestTrue( arr );
   TestTrue( arr->size == 0 );
   TestTrue( arr->ptr == NULL );

  done:
   NewClass_Delete( arr );
}
TestEnd

TestBegin( Set ) {
   IArray* arr;
   int itms[5] = {0, 1, 2, 3, 4};
   const int* ptr;
   int i_i;

   arr = IArray_New();
   TestNoAssert( IArray_Set( arr, 5, itms ) );
   TestTrue( IArray_GetSize( arr ) == 5 );
   TestTrue( (ptr = IArray_GetPtr( arr )) != 0 );
   for( i_i = 0; i_i < 5; i_i++ ) {
      TestTrue( ptr[i_i] == i_i );
   }

  done:
   NewClass_Delete( arr );
}
TestEnd

TestBegin( Add ) {
   IArray* arr;
   int itms[5] = {0, 1, 2, 3, 4};
   const int* ptr;
   int i_i;

   arr = IArray_New();
   IArray_Set( arr, 3, itms );
   TestNoAssert( IArray_Add( arr, 2, itms + 3 ) );
   TestTrue( IArray_GetSize( arr ) == 5 );
   TestTrue( (ptr = IArray_GetPtr( arr )) != 0 );
   for( i_i = 0; i_i < 5; i_i++ ) {
      TestTrue( ptr[i_i] == i_i );
   }

  done:
   NewClass_Delete( arr );
}
TestEnd

TestBegin( Remove ) {
   IArray* arr;
   int itms[5] = {0, 1, 2, 3, 4};
   const int* ptr;
   IMap mapObj, *map = &mapObj;

   arr = IArray_New();
   IArray_Set( arr, 5, itms );
   itms[0] = 1; itms[1] = 3;
   IMap_Init( map );
   TestNoAssert( IArray_Remove( arr, 2, itms, map ) );
   TestTrue( IArray_GetSize( arr ) == 3 );
   TestTrue( (ptr = IArray_GetPtr( arr )) != 0 );
   TestTrue( ptr[0] == 0 && ptr[1] == 4 && ptr[2] == 2 );
   TestTrue( IMap_GetSize( map ) == 1 );
   TestTrue( IMap_Has( map, 4 ) );
   TestTrue( IMap_Map( map, 4 ) == 1 );

  done:
   NewClass_Delete( arr );
   IMap_Destruct( map );
}
TestEnd


#define nTests 4
TestSuite_Test tests[nTests] = {{"construct", testConstruct}, 
				{"set", testSet}, 
				{"add", testAdd}, 
				{"remove", testRemove}};


#include "Base/Foundation/TestEnd.h"
