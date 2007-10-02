/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003, Victorian Partnership for Advanced Computing (VPAC) Ltd, 
** 110 Victoria Street, Melbourne, 3053, Australia.
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

TestBegin( Insert ) {
   STreeMap* map;
   int c_i;

   map = STreeMap_New();
   STreeMap_SetItemSize( map, sizeof(int), sizeof(int) );
   STree_SetIntCallbacks( map );

   for ( c_i = 0; c_i < 10; c_i++ ) {
      int tmp = 10 * c_i;
      STreeMap_Insert( map, &c_i, &tmp );
   }
   TestTrue( STree_GetSize( map ) == 10 );

  done:
   NewClass_Delete( map );
}
TestEnd

TestBegin( Map ) {
   STreeMap* map;
   int c_i;

   map = STreeMap_New();
   STreeMap_SetItemSize( map, sizeof(int), sizeof(int) );
   STree_SetIntCallbacks( map );

   for ( c_i = 0; c_i < 10; c_i++ ) {
      int tmp = 10 * c_i;
      STreeMap_Insert( map, &c_i, &tmp );
   }
   for( c_i = 0; c_i < 10; c_i++ ) {
      TestTrue( *((int*)STreeMap_Map( map, &c_i )) == 10 * c_i );
   }

  done:
   NewClass_Delete( map );
}
TestEnd


#define nTests 2
TestSuite_Test tests[nTests] = {{"construct", testInsert},
				{"has", testMap}};


#include "Base/Foundation/TestEnd.h"
