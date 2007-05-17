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
   ISet* set;

   TestNoAssert( set = ISet_New() );
   TestTrue( set );
   TestTrue( set->maxSize == 0 );
   TestTrue( set->curSize == 0 );
   TestTrue( set->tblSize == 0 );
   TestTrue( set->tbl == NULL );
   TestTrue( set->used == NULL );

  done:
   NewClass_Delete( set );
}
TestEnd

TestBegin( Insert ) {
   ISet* set;
   int i_i;

   set = ISet_New();
   ISet_SetMaxSize( set, 20 );
   for( i_i = 0; i_i < 20; i_i += 2 ) {
      TestNoAssert( ISet_Insert( set, i_i ) );
   }
   for( i_i = 1; i_i < 20; i_i += 2 ) {
      TestNoAssert( ISet_Insert( set, i_i ) );
   }
   TestTrue( ISet_GetSize( set ) == 20 );
   TestAssert( ISet_Insert( set, 0 ) );
   TestNoAssert( ISet_TryInsert( set, 0 ) );
   for( i_i = 0; i_i < 20; i_i++ ) {
      TestTrue( ISet_Has( set, i_i ) );
   }
   TestTrue( !ISet_Has( set, 21 ) );

  done:
   NewClass_Delete( set );
}
TestEnd

TestBegin( UseArray ) {
   int array[10] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
   ISet* set;
   int i_i;

   set = ISet_New();
   TestNoAssert( ISet_UseArray( set, 10, array ) );
   TestTrue( ISet_GetSize( set ) == 10 );
   for( i_i = 0; i_i < 10; i_i++ ) {
      TestTrue( ISet_Has( set, i_i ) );
   }
   TestTrue( !ISet_Has( set, 10 ) );

  done:
   NewClass_Delete( set );
}
TestEnd

TestBegin( Union ) {
   int array0[5] = {0, 1, 2, 3, 4};
   int array1[5] = {3, 4, 5, 6, 7};
   ISet *set0, *set1;
   int i_i;

   set0 = ISet_New();
   ISet_UseArray( set0, 5, array0 );
   set1 = ISet_New();
   ISet_UseArray( set1, 5, array1 );
   ISet_Union( set0, set1 );
   TestTrue( ISet_GetSize( set0 ) == 8 );
   for( i_i = 0; i_i < 8; i_i++ ) {
      TestTrue( ISet_Has( set0, i_i ) );
   }
   TestTrue( !ISet_Has( set0, 10 ) );

  done:
   NewClass_Delete( set0 );
   NewClass_Delete( set1 );
}
TestEnd

TestBegin( Isect ) {
   int array0[5] = {0, 1, 2, 3, 4};
   int array1[5] = {3, 4, 5, 6, 7};
   ISet *set0, *set1;
   int i_i;

   set0 = ISet_New();
   ISet_UseArray( set0, 5, array0 );
   set1 = ISet_New();
   ISet_UseArray( set1, 5, array1 );
   ISet_Isect( set0, set1 );
   TestTrue( ISet_GetSize( set0 ) == 2 );
   for( i_i = 3; i_i < 5; i_i++ ) {
      TestTrue( ISet_Has( set0, i_i ) );
   }
   TestTrue( !ISet_Has( set0, 1 ) );

  done:
   NewClass_Delete( set0 );
   NewClass_Delete( set1 );
}
TestEnd

TestBegin( Subtr ) {
   int array0[5] = {0, 1, 2, 3, 4};
   int array1[5] = {3, 4, 5, 6, 7};
   ISet *set0, *set1;
   int i_i;

   set0 = ISet_New();
   ISet_UseArray( set0, 5, array0 );
   set1 = ISet_New();
   ISet_UseArray( set1, 5, array1 );
   ISet_Subtr( set0, set1 );
   TestTrue( ISet_GetSize( set0 ) == 3 );
   for( i_i = 0; i_i < 3; i_i++ ) {
      TestTrue( ISet_Has( set0, i_i ) );
   }
   TestTrue( !ISet_Has( set0, 3 ) );

  done:
   NewClass_Delete( set0 );
   NewClass_Delete( set1 );
}
TestEnd


#define nTests 6
TestSuite_Test tests[nTests] = {{"construct", testConstruct}, 
				{"insert", testInsert}, 
				{"use array", testUseArray}, 
				{"union", testUnion}, 
				{"intersection", testIsect}, 
				{"subtraction", testSubtr}};


#include "Base/Foundation/TestEnd.h"
