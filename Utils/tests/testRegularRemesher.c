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
** $Id: testRegularRemesher.c 3952 2007-01-09 06:24:06Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/Geometry/Geometry.h>
#include <StgDomain/Shape/Shape.h>
#include <StgDomain/Mesh/Mesh.h>
#include <StgDomain/Utils/Utils.h>

#include "StGermain/Base/Foundation/TestBegin.h"


Mesh* buildMesh() {
   CartesianGenerator* gen;
   int nRanks;
   unsigned sizes[3];
   double minCrd[3];
   double maxCrd[3];
   Mesh* mesh;

   insist( MPI_Comm_size( MPI_COMM_WORLD, &nRanks ), == MPI_SUCCESS );
   sizes[0] = sizes[1] = sizes[2] = 2 * nRanks;
   minCrd[0] = minCrd[1] = minCrd[2] = 0.0;
   maxCrd[0] = maxCrd[1] = maxCrd[2] = (double)nRanks;

   gen = CartesianGenerator_New( "" );
   MeshGenerator_SetDimSize( gen, 3 );
   CartesianGenerator_SetShadowDepth( gen, 1 );
   CartesianGenerator_SetTopologyParams( gen, sizes, 0, NULL, NULL );
   CartesianGenerator_SetGeometryParams( gen, minCrd, maxCrd );

   mesh = Mesh_New( "" );
   CartesianGenerator_Generate( gen, mesh );
   FreeObject( gen );
   Stg_Component_Build( mesh, NULL, False );
   Stg_Component_Initialise( mesh, NULL, False );

   return mesh;
}


void testSetup( int* argc, char** argv[] ) {
   StGermain_Init( argc, argv );
   StgDomainMesh_Init( argc, argv );
   StgDomainUtils_Init( argc, argv );
}

void testTeardown() {
   StgDomainUtils_Finalise();
   StgDomainMesh_Finalise();
   StGermain_Finalise();
}

TestBegin( Construct ) {
   RegularRemesher* remesh;

   TestNoAssert( remesh = RegularRemesher_New() );
   TestTrue( remesh );

  done:
   NewClass_Delete( remesh );
}
TestEnd

TestBegin( Build ) {
   RegularRemesher* remesh;
   Mesh* mesh;

   mesh = buildMesh();
   remesh = RegularRemesher_New();
   NewRemesher_SetMesh( remesh, mesh );
   RegularRemesher_Build( remesh );

   RegularRemesher_SetRemeshState( remesh, 0, True );
   RegularRemesher_SetRemeshState( remesh, 1, True );
   RegularRemesher_SetRemeshState( remesh, 2, True );
   RegularRemesher_Build( remesh );

   RegularRemesher_SetStaticWall( remesh, 0, 0, True );
   RegularRemesher_SetStaticWall( remesh, 0, 1, True );
   RegularRemesher_SetStaticWall( remesh, 2, 0, True );
   RegularRemesher_SetStaticWall( remesh, 2, 1, True );
   RegularRemesher_Build( remesh );

  done:
   NewClass_Delete( remesh );
   FreeObject( mesh );
}
TestEnd

TestBegin( Remesh ) {
   RegularRemesher* remesh;
   Mesh* mesh;
   Grid* vGrid;
   int nVerts;
   double** oldVerts;
   int ind, inds[3];
   int v_i, d_i;

   mesh = buildMesh();
   remesh = RegularRemesher_New();
   NewRemesher_SetMesh( remesh, mesh );

   RegularRemesher_SetRemeshState( remesh, 0, True );
   RegularRemesher_SetRemeshState( remesh, 1, True );
   RegularRemesher_SetRemeshState( remesh, 2, True );
   RegularRemesher_SetStaticWall( remesh, 0, 0, True );
   RegularRemesher_SetStaticWall( remesh, 0, 1, True );
   RegularRemesher_SetStaticWall( remesh, 2, 0, True );
   RegularRemesher_SetStaticWall( remesh, 2, 1, True );
   RegularRemesher_Build( remesh );

   nVerts = Mesh_GetDomainSize( mesh, 0 );
   oldVerts = MemArray2D( double, nVerts, 3, "testRegularRemesher" );

   memcpy( oldVerts[0], mesh->verts[0], nVerts * 3 * sizeof(double) );
   RegularRemesher_Remesh( remesh );
   for( v_i = 0; v_i < nVerts * 3; v_i++ ) {
      if( !Num_Approx( oldVerts[0][v_i], mesh->verts[0][v_i] ) )
	 break;
   }
   TestTrue( v_i == nVerts * 3 );

   vGrid = *Mesh_GetExtension( mesh, Grid**, "vertexGrid" );
   for( v_i = 0; v_i < nVerts; v_i++ ) {
      ind = Mesh_DomainToGlobal( mesh, 0, v_i );
      Grid_Lift( vGrid, ind, (unsigned*)inds );
      for( d_i = 0; d_i < 3; d_i++ ) {
	 if( inds[d_i] == 0 || inds[d_i] == vGrid->sizes[d_i] - 1 )
	    continue;
	 mesh->verts[v_i][d_i] = -1.0;
      }
   }
   RegularRemesher_Remesh( remesh );
   for( v_i = 0; v_i < nVerts * 3; v_i++ ) {
      if( !Num_Approx( oldVerts[0][v_i], mesh->verts[0][v_i] ) )
	 break;
   }
   TestTrue( v_i == nVerts * 3 );

  done:
   NewClass_Delete( remesh );
   FreeObject( mesh );
   MemFree( oldVerts );
}
TestEnd


#define nTests 3
TestSuite_Test tests[nTests] = {{"construct", testConstruct}, 
				{"build", testBuild}, 
				{"remesh", testRemesh}};


#include "StGermain/Base/Foundation/TestEnd.h"
