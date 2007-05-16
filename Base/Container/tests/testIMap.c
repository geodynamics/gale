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
   IMap* map;

   TestNoAssert( map = IMap_New() );
   TestTrue( map );
   TestTrue( map->maxItms == 0 );
   TestTrue( map->nItms == 0 );
   TestTrue( map->tblSize == 0 );
   TestTrue( map->tbl == NULL );
   TestTrue( map->used == NULL );

  done:
   NewClass_Delete( map );
}
TestEnd

TestBegin( SetMaxItms ) {
   IMap* map;

   map = IMap_New();
   TestNoAssert( IMap_SetMaxItems( map, 10 ) );
   TestTrue( map->maxItms == 10 );
   TestTrue( map->tblSize >= 10 );
   TestNoAssert( IMap_SetMaxItems( map, 20 ) );
   TestTrue( map->maxItms == 20 );
   TestTrue( map->tblSize >= 20 );

  done:
   NewClass_Delete( map );
}
TestEnd

TestBegin( Insert ) {
   IMap* map;
   int i_i;

   map = IMap_New();
   IMap_SetMaxItems( map, 10 );
   for( i_i = 0; i_i < 20; i_i += 2 ) {
      TestNoAssert( IMap_Insert( map, i_i, i_i + 100 ) );
   }
   TestAssert( IMap_Insert( map, 0, 100 ) );
   TestTrue( IMap_GetNumItems( map ) == 10 );
   for( i_i = 0; i_i < 20; i_i += 2 ) {
      TestTrue( IMap_Has( map, i_i ) );
   }

  done:
   NewClass_Delete( map );
}
TestEnd

TestBegin( Map ) {
   IMap* map;
   int i_i;

   map = IMap_New();
   IMap_SetMaxItems( map, 10 );
   for( i_i = 0; i_i < 20; i_i += 2 ) {
      IMap_Insert( map, i_i, i_i + 100 );
   }
   for( i_i = 0; i_i < 20; i_i += 2 ) {
      TestTrue( IMap_Map( map, i_i ) == i_i + 100 );
   }

  done:
   NewClass_Delete( map );
}
TestEnd

TestBegin( Remove ) {
   IMap* map;
   int i_i;

   map = IMap_New();
   IMap_SetMaxItems( map, 10 );
   for( i_i = 0; i_i < 20; i_i += 2 ) {
      IMap_Insert( map, i_i, i_i + 100 );
   }
   for( i_i = 0; i_i < 20; i_i += 4 ) {
      TestNoAssert( IMap_Remove( map, i_i ) );
   }
   TestAssert( IMap_Remove( map, 0 ) );
   TestTrue( IMap_GetNumItems( map ) == 5 );
   for( i_i = 0; i_i < 20; i_i += 4 ) {
      if( i_i % 4 ) {
	 TestTrue( IMap_Has( map, i_i ) );
      }
      else {
	 TestTrue( !IMap_Has( map, i_i ) );
      }
   }

  done:
   NewClass_Delete( map );
}
TestEnd


#define nTests 5
TestSuite_Test tests[nTests] = {{"construct", testConstruct}, 
				{"set maximum items", testSetMaxItms}, 
				{"insert", testInsert}, 
				{"map", testMap}, 
				{"remove", testRemove}};


#include "Base/Foundation/TestEnd.h"
