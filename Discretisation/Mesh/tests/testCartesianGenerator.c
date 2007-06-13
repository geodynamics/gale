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
** $Id: testSync.c 3952 2007-01-09 06:24:06Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include "StGermain/Base/Base.h"
#include "StGermain/Discretisation/Mesh/Mesh.h"

#include "StGermain/Base/Foundation/TestBegin.h"


Mesh* buildMesh() {
   CartesianGenerator* gen;
   int nRanks;
   unsigned sizes[3];
   double minCrd[3];
   double maxCrd[3];
   Mesh* mesh;
   int rank;

   insist( MPI_Comm_size( MPI_COMM_WORLD, &nRanks ), == MPI_SUCCESS );
   sizes[0] = sizes[1] = sizes[2] = nRanks * 4;
   minCrd[0] = minCrd[1] = minCrd[2] = 0.0;
   maxCrd[0] = minCrd[1] = minCrd[2] = (double)nRanks;

   gen = CartesianGenerator_New( "" );
   MeshGenerator_SetDimSize( gen, 3 );
   CartesianGenerator_SetShadowDepth( gen, 1 );
   CartesianGenerator_SetTopologyParams( gen, sizes, 0, NULL, NULL );
   CartesianGenerator_SetGeometryParams( gen, minCrd, maxCrd );

   mesh = Mesh_New( "" );
   CartesianGenerator_Generate( gen, mesh );
   FreeObject( gen );

   return mesh;
}

void testSetup( int* argc, char** argv[] ) {
   Base_Init( argc, argv );
   DiscretisationMesh_Init( argc, argv );
   Stream_Enable( Journal_GetTypedStream( Debug_Type ), False );
   Stream_Enable( Journal_GetTypedStream( Info_Type ), False );
   Stream_Enable( Journal_GetTypedStream( Dump_Type ), False );
   Stream_Enable( Journal_GetTypedStream( Error_Type ), False );
}

void testTeardown() {
   DiscretisationMesh_Finalise();
   Base_Finalise();
}

TestBegin( Construct ) {
   CartesianGenerator* gen;

   TestNoAssert( gen = CartesianGenerator_New( "" ) );
   TestTrue( gen );

  done:
   FreeObject( gen );
}
TestEnd

TestBegin( Gen ) {
   CartesianGenerator* gen;
   int nRanks;
   unsigned sizes[3];
   double minCrd[3];
   double maxCrd[3];
   Mesh* mesh;
   int rank;
   unsigned long netMem, mem, syncMem;

   insist( MPI_Comm_size( MPI_COMM_WORLD, &nRanks ), == MPI_SUCCESS );
   sizes[0] = sizes[1] = sizes[2] = nRanks * 4;
   minCrd[0] = minCrd[1] = minCrd[2] = 0.0;
   maxCrd[0] = minCrd[1] = minCrd[2] = (double)nRanks;

   gen = CartesianGenerator_New( "" );
   MeshGenerator_SetDimSize( gen, 3 );
   CartesianGenerator_SetShadowDepth( gen, 1 );
   CartesianGenerator_SetTopologyParams( gen, sizes, 0, NULL, NULL );
   CartesianGenerator_SetGeometryParams( gen, minCrd, maxCrd );

   mesh = Mesh_New( "" );
   TestNoAssert( CartesianGenerator_Generate( gen, mesh ) );

  done:
   FreeObject( gen );
   FreeObject( mesh );
}
TestEnd

TestBegin( Inc ) {
   Mesh* mesh;
   MeshTopology* topo;
   const Sync *elSync, *vertSync;
   Grid *elGrid, *vertGrid;
   int nEls, nDims;
   int* elParam;
   int nIncEls;
   const int* incEls;
   int elGlobal, incGlobal, vertGlobal;
   int e_i;

   mesh = buildMesh();
   topo = mesh->topo;
   nDims = MeshTopology_GetNumDims( topo );
   elSync = MeshTopology_GetDomain( topo, nDims );
   elGrid = *Mesh_GetExtension( mesh, Grid**, "elementGrid" );
   vertSync = MeshTopology_GetDomain( topo, 0 );
   vertGrid = *Mesh_GetExtension( mesh, Grid**, "vertexGrid" );
   nEls = Sync_GetNumDomains( elSync );
   elParam = MemArray( int, nDims, "testInc" );
   for( e_i = 0; e_i < nEls; e_i++ ) {
      elGlobal = Sync_DomainToGlobal( elSync, e_i );
      Grid_Lift( elGrid, elGlobal, elParam );
      MeshTopology_GetIncidence( topo, nDims, e_i, 0, &nIncEls, &incEls );
      if( nIncEls != 8 )
	 break;

      vertGlobal = Grid_Project( vertGrid, elParam );
      incGlobal = Sync_DomainToGlobal( vertSync, incEls[0] );
      if( incGlobal != vertGlobal )
	 break;

      elParam[0]++;
      vertGlobal = Grid_Project( vertGrid, elParam );
      incGlobal = Sync_DomainToGlobal( vertSync, incEls[1] );
      if( incGlobal != vertGlobal )
	 break;
      elParam[0]--;

      elParam[1]++;
      vertGlobal = Grid_Project( vertGrid, elParam );
      incGlobal = Sync_DomainToGlobal( vertSync, incEls[2] );
      if( incGlobal != vertGlobal )
	 break;

      elParam[0]++;
      vertGlobal = Grid_Project( vertGrid, elParam );
      incGlobal = Sync_DomainToGlobal( vertSync, incEls[3] );
      if( incGlobal != vertGlobal )
	 break;
      elParam[0]--;
      elParam[1]--;

      elParam[2]++;
      vertGlobal = Grid_Project( vertGrid, elParam );
      incGlobal = Sync_DomainToGlobal( vertSync, incEls[4] );
      if( incGlobal != vertGlobal )
	 break;

      elParam[0]++;
      vertGlobal = Grid_Project( vertGrid, elParam );
      incGlobal = Sync_DomainToGlobal( vertSync, incEls[5] );
      if( incGlobal != vertGlobal )
	 break;
      elParam[0]--;

      elParam[1]++;
      vertGlobal = Grid_Project( vertGrid, elParam );
      incGlobal = Sync_DomainToGlobal( vertSync, incEls[6] );
      if( incGlobal != vertGlobal )
	 break;

      elParam[0]++;
      vertGlobal = Grid_Project( vertGrid, elParam );
      incGlobal = Sync_DomainToGlobal( vertSync, incEls[7] );
      if( incGlobal != vertGlobal )
	 break;
      elParam[0]--;
      elParam[1]--;
      elParam[2]--;
   }
   TestTrue( e_i == nEls );

  done:
   MemFree( elParam );
   FreeObject( mesh );
}
TestEnd


#define nTests 3
TestSuite_Test tests[nTests] = {{"construct", testConstruct}, 
				{"generate", testGen}, 
				{"incidence", testInc}};


#include "Base/Foundation/TestEnd.h"
