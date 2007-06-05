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
** $Id: testDecomp.c 3952 2007-01-09 06:24:06Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include "StGermain/Base/Base.h"
#include "StGermain/Discretisation/Mesh/Mesh.h"

#include "StGermain/Base/Foundation/TestBegin.h"


void testSetup( int* argc, char** argv[] ) {
   Base_Init( argc, argv );
   DiscretisationMesh_Init( argc, argv );
}

void testTeardown() {
   DiscretisationMesh_Finalise();
   Base_Finalise();
}

TestBegin( Construct ) {
   Decomp* decomp;
   int nLocs;
   const int *locs;

   TestNoAssert( decomp = Decomp_New() );
   TestTrue( decomp );
   TestTrue( Decomp_GetMPIComm( decomp ) == MPI_COMM_WORLD );
   TestTrue( Decomp_GetNumGlobals( decomp ) == 0 );
   TestTrue( Decomp_GetNumLocals( decomp ) == 0 );
   TestNoAssert( Decomp_GetLocals( decomp, &nLocs, &locs ) );
   TestTrue( nLocs == 0 && locs == NULL );

  done:
   NewClass_Delete( decomp );
}
TestEnd

TestBegin( Owners ) {
   Decomp* decomp;
   int rank, nRanks;
   int nLocs, *locs, *ranks;
   int l_i, g_i;

   insist( MPI_Comm_rank( MPI_COMM_WORLD, &rank ), == MPI_SUCCESS );
   insist( MPI_Comm_size( MPI_COMM_WORLD, &nRanks ), == MPI_SUCCESS );

   nLocs = 10;
   locs = MemArray( int, nLocs, "testDecomp" );

   decomp = Decomp_New();
   for( l_i = 0; l_i < nLocs; l_i++ )
      locs[l_i] = rank * nLocs + l_i;
   TestNoAssert( Decomp_SetLocals( decomp, nLocs, locs ) );
   for( g_i = 0; g_i < nRanks * nLocs; g_i++ ) {
      if( g_i >= rank * nLocs && g_i < (rank + 1) * nLocs ) {
	 TestTrue( IMap_Map( decomp->owners, g_i ) == rank );
      }
      else {
	 TestTrue( !IMap_Has( decomp->owners, g_i ) );
      }
   }

   for( l_i = 0; l_i < nLocs; l_i++ ) {
      locs[l_i] = (rank * nLocs + nLocs / 2 + l_i) % 
	 (nRanks * nLocs);
   }
   TestNoAssert( Decomp_SetLocals( decomp, nLocs, locs ) );
   for( g_i = 0; g_i < nRanks * nLocs; g_i++ ) {
      if( g_i >= rank * nLocs && g_i < (rank + 1) * nLocs ) {
	 if( g_i < rank * nLocs + nLocs / 2 ) {
	    if( rank > 0 ) {
	       TestTrue( IMap_Map( decomp->owners, g_i ) == rank - 1 );
	    }
	    else {
	       TestTrue( IMap_Map( decomp->owners, g_i ) == nRanks - 1 );
	    }
	 }
	 else {
	    TestTrue( IMap_Map( decomp->owners, g_i ) == rank );
	 }
      }
      else {
	 TestTrue( !IMap_Has( decomp->owners, g_i ) );
      }
   }

   locs = MemRearray( locs, int, nRanks * nLocs, "testDecomp" );
   ranks = MemArray( int, nRanks * nLocs, "testDecomp" );
   for( g_i = 0; g_i < nRanks * nLocs; g_i++ )
      locs[g_i] = g_i;
   TestNoAssert( Decomp_FindOwners( decomp, nRanks * nLocs, locs, ranks ) );

  done:
   NewClass_Delete( decomp );
   MemFree( locs );
   MemFree( ranks );
}
TestEnd


#define nTests 2
TestSuite_Test tests[nTests] = {{"construct", testConstruct}, 
				{"owners", testOwners }};


#include "Base/Foundation/TestEnd.h"
