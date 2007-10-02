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


void calcDepth( const STreeNode* node, int* curDepth, int* depth ) {
   if( !node )
      return;
   if( ++(*curDepth) > *depth )
      *depth = *curDepth;
   if( node->left )
      calcDepth( node->left, curDepth, depth );
   if( node->right )
      calcDepth( node->right, curDepth, depth );
   (*curDepth)--;
}


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
   STree* st;

   TestNoAssert( st = STree_New() );
   TestTrue( STree_GetSize( st ) == 0 );

  done:
   NewClass_Delete( st );
}
TestEnd

TestBegin( Insert ) {
   STree* st;
   int depth = 0, curDepth = 0;
   int i_i;

   st = STree_New();
   STree_SetItemSize( st, sizeof(int) );
   STree_SetIntCallbacks( st );
   STree_SetAlpha( st, 0.5 );
   for( i_i = 0; i_i < 15; i_i++ )
      STree_Insert( st, &i_i );
   TestTrue( STree_GetSize( st ) == 15 );
   TestTrue( (calcDepth( STree_GetRoot( st ), &curDepth, &depth), depth) == 4 );

  done:
   NewClass_Delete( st );
}
TestEnd

TestBegin( Remove ) {
   STree* st;
   int depth = 0, curDepth = 0;
   int i_i;

   st = STree_New();
   STree_SetItemSize( st, sizeof(int) );
   STree_SetIntCallbacks( st );
   STree_SetAlpha( st, 0.5 );
   for( i_i = 0; i_i < 30; i_i++ )
      STree_Insert( st, &i_i );
   for( i_i = 0; i_i < 30; i_i += 2 )
      STree_Remove( st, &i_i );
   TestTrue( STree_GetSize( st ) == 15 );
   TestTrue( (calcDepth( STree_GetRoot( st ), &curDepth, &depth), depth) == 4 );

  done:
   NewClass_Delete( st );
}
TestEnd

TestBegin( Has ) {
   STree* st;
   int i_i;

   st = STree_New();
   STree_SetItemSize( st, sizeof(int) );
   STree_SetIntCallbacks( st );
   STree_SetAlpha( st, 0.5 );
   for( i_i = 0; i_i < 15; i_i++ )
      STree_Insert( st, &i_i );
   for( i_i = 0; i_i < 15; i_i++ ) {
      TestTrue( STree_Has( st, &i_i ) );
   }
   for( i_i = 15; i_i < 30; i_i++ ) {
      TestTrue( !STree_Has( st, &i_i ) );
   }

  done:
   NewClass_Delete( st );
}
TestEnd


#define nTests 4
TestSuite_Test tests[nTests] = {{"construct", testConstruct}, 
				{"insert", testInsert}, 
				{"remove", testRemove},
				{"has", testHas}};


#include "Base/Foundation/TestEnd.h"
