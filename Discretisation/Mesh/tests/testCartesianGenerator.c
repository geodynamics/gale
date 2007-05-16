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


void testSetup( int* argc, char** argv[] ) {
   Base_Init( argc, argv );
   DiscretisationMesh_Init( argc, argv );
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
   sizes[0] = sizes[1] = sizes[2] = nRanks;
   minCrd[0] = minCrd[1] = minCrd[2] = 0.0;
   maxCrd[0] = minCrd[1] = minCrd[2] = (double)nRanks;

   gen = CartesianGenerator_New( "" );
   MeshGenerator_SetDimSize( gen, 3 );
   CartesianGenerator_SetShadowDepth( gen, 0 );
   CartesianGenerator_SetTopologyParams( gen, sizes, 0, NULL, NULL );
   CartesianGenerator_SetGeometryParams( gen, minCrd, maxCrd );

   mesh = Mesh_New( "" );
   TestNoAssert( CartesianGenerator_Generate( gen, mesh ) );

  done:
   FreeObject( gen );
   FreeObject( mesh );
}
TestEnd


#define nTests 2
TestSuite_Test tests[nTests] = {{"construct", testConstruct}, 
				{"generate", testGen}};


#include "Base/Foundation/TestEnd.h"
