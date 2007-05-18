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
** $Id: testMeshTopology.c 3952 2007-01-09 06:24:06Z LukeHodkinson $
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
   /*DiscretisationMesh_Init( argc, argv );*/
}

void testTeardown() {
   /*DiscretisationMesh_Finalise();*/
   Base_Finalise();
}

TestBegin( Construct ) {
   MeshTopology* topo;

   TestNoAssert( topo = MeshTopology_New() );
   TestTrue( topo );

  done:
   NewClass_Delete( topo );
}
TestEnd

TestBegin( SetDims ) {
   MeshTopology* topo;
   int d_i, d_j;

   topo = MeshTopology_New();
   TestNoAssert( MeshTopology_SetNumDims( topo, 3 ) );
   TestTrue( MeshTopology_GetNumDims( topo ) == 3 && topo->nTDims == 4 );
   TestTrue( topo->locals && topo->remotes );
   TestTrue( topo->nIncEls && topo->incEls );
   for( d_i = 0; d_i < 4; d_i++ ) {
      TestTrue( topo->locals[d_i] && topo->remotes[d_i] );
      TestTrue( Sync_GetDecomp( topo->remotes[d_i] ) == topo->locals[d_i] );
      TestTrue( topo->nIncEls[d_i] && topo->incEls[d_i] );
      TestTrue( topo->incEls[d_i] );
      for( d_j = 0; d_j < 4; d_j++ ) {
	 TestTrue( !topo->nIncEls[d_i][d_j] && !topo->incEls[d_i][d_j] );
      }
   }

   TestNoAssert( MeshTopology_SetNumDims( topo, 0 ) );
   TestTrue( !MeshTopology_GetNumDims( topo ) && !topo->nTDims );
   TestTrue( !topo->locals && !topo->remotes );
   TestTrue( !topo->nIncEls && !topo->incEls );

  done:
   NewClass_Delete( topo );
}
TestEnd

TestBegin( SetComm ) {
   MeshTopology* topo;
   Comm* comm;
   int d_i;

   comm = Comm_New();
   topo = MeshTopology_New();
   MeshTopology_SetNumDims( topo, 3 );

   TestNoAssert( MeshTopology_SetComm( topo, comm ) );
   TestTrue( MeshTopology_GetComm( topo ) == comm );
   for( d_i = 0; d_i < 4; d_i++ ) {
      TestTrue( Sync_GetComm( topo->remotes[d_i] ) == comm );
   }

  done:
   NewClass_Delete( topo );
}
TestEnd


#define nTests 3
TestSuite_Test tests[nTests] = {{"construct", testConstruct}, 
				{"set number of dimensions", testSetDims}, 
				{"set communicator", testSetComm}};


#include "Base/Foundation/TestEnd.h"
