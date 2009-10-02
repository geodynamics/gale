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
** $Id: SpatialTree.c 3952 2007-01-09 06:24:06Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/Mesh/Mesh.h>

#include "StGermain/Base/Foundation/TestBegin.h"


void testSetup( int* argc, char** argv[] ) {
   StGermain_Init( argc, argv );
   StgDomainMesh_Init( argc, argv );
   Stream_Enable( Journal_GetTypedStream( Debug_Type ), False );
   Stream_Enable( Journal_GetTypedStream( Info_Type ), False );
   Stream_Enable( Journal_GetTypedStream( Dump_Type ), False );
   Stream_Enable( Journal_GetTypedStream( Error_Type ), False );
}

void testTeardown() {
   StGermain_Finalise();
}

TestBegin( All ) {
   SpatialTree* tree;
   Mesh* mesh;
   int nRanks;
   int sizes[3];
   double minCrd[3], maxCrd[3];
   CartesianGenerator* gen;
   double pnt[3];
   int nDims;
   int nEls, *els;

   MPI_Comm_size( MPI_COMM_WORLD, &nRanks );
   sizes[0] = sizes[1] = sizes[2] = 4 * nRanks;
   minCrd[0] = minCrd[1] = minCrd[2] = 0.0;
   maxCrd[0] = maxCrd[1] = maxCrd[2] = (double)nRanks;

   nDims = 3;
   gen = CartesianGenerator_New( "" );
   MeshGenerator_SetDimSize( gen, nDims );
   CartesianGenerator_SetShadowDepth( gen, 1 );
   CartesianGenerator_SetTopologyParams( gen, (unsigned*)sizes, 0, NULL, NULL );
   CartesianGenerator_SetGeometryParams( gen, minCrd, maxCrd );
   mesh = Mesh_New( "" );
   Mesh_SetGenerator( mesh, gen );
   Stg_Component_Build( mesh, NULL, False );

   tree = SpatialTree_New();
   SpatialTree_SetMesh( tree, mesh );
   TestTrue( tree );

   SpatialTree_Rebuild( tree );

   pnt[0] = pnt[1] = pnt[2] = 0.0;
   TestTrue( SpatialTree_Search( tree, pnt, &nEls, &els ) );

  done:
   NewClass_Delete( tree );
   Stg_Class_Delete( mesh );
}
TestEnd


#define nTests 1
TestSuite_Test tests[nTests] = {{"all", testAll}};

#include "StGermain/Base/Foundation/TestEnd.h"
