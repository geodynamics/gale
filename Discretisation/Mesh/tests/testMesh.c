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

TestBegin( NearVert ) {
   CartesianGenerator* gen;
   Mesh* mesh;
   int sizes[3];
   double minCrd[3];
   double maxCrd[3];
   int nRanks;
   int nInc, *inc;
   double* vert;
   int e_i, inc_i;

   insist( MPI_Comm_size( MPI_COMM_WORLD, &nRanks ), == MPI_SUCCESS );
   sizes[0] = sizes[1] = sizes[2] = nRanks;
   minCrd[0] = minCrd[1] = minCrd[2] = 0.0;
   maxCrd[0] = maxCrd[1] = maxCrd[2] = (double)nRanks;

   gen = CartesianGenerator_New( "" );
   MeshGenerator_SetDimSize( gen, 3 );
   CartesianGenerator_SetShadowDepth( gen, 1 );
   CartesianGenerator_SetTopologyParams( gen, sizes, 0, NULL, NULL );
   CartesianGenerator_SetGeometryParams( gen, minCrd, maxCrd );
   mesh = Mesh_New( "" );
   Mesh_SetGenerator( mesh, gen );
   Stg_Component_Build( mesh, NULL, False );

   for( e_i = 0; e_i < Mesh_GetDomainSize( mesh, MT_VOLUME ); e_i++ ) {
      Mesh_GetIncidence( mesh, MT_VOLUME, e_i, MT_VERTEX, &nInc, &inc );
      for( inc_i = 0; inc_i < nInc; inc_i++ ) {
	 vert = Mesh_GetVertex( mesh, inc[inc_i] );
	 TestTrue( Mesh_NearestVertex( mesh, vert ) == inc[inc_i] );
      }
   }

done:
   FreeObject( gen );
   FreeObject( mesh );
}
TestEnd

TestBegin( ElSearch3D ) {
   CartesianGenerator* gen;
   Mesh* mesh;
   int sizes[3];
   double minCrd[3];
   double maxCrd[3];
   int nRanks;
   int nInc, *inc;
   double* vert;
   double point[3];
   unsigned elDim, elInd;
   int nEdgeInc, *edgeInc;
   int nFaceInc, *faceInc;
   int e_i, inc_i, inc_j;

   insist( MPI_Comm_size( MPI_COMM_WORLD, &nRanks ), == MPI_SUCCESS );
   sizes[0] = sizes[1] = sizes[2] = nRanks;
   minCrd[0] = minCrd[1] = minCrd[2] = 0.0;
   maxCrd[0] = maxCrd[1] = maxCrd[2] = (double)nRanks;

   gen = CartesianGenerator_New( "" );
   MeshGenerator_SetDimSize( gen, 3 );
   CartesianGenerator_SetShadowDepth( gen, 1 );
   CartesianGenerator_SetTopologyParams( gen, sizes, 0, NULL, NULL );
   CartesianGenerator_SetGeometryParams( gen, minCrd, maxCrd );
   mesh = Mesh_New( "" );
   Mesh_SetGenerator( mesh, gen );
   Stg_Component_Build( mesh, NULL, False );

   for( e_i = 0; e_i < Mesh_GetDomainSize( mesh, MT_VOLUME ); e_i++ ) {
      Mesh_GetIncidence( mesh, MT_VOLUME, e_i, MT_VERTEX, &nInc, &inc );
      for( inc_i = 0; inc_i < nInc; inc_i++ ) {
	 vert = Mesh_GetVertex( mesh, inc[inc_i] );
	 TestTrue( Mesh_Search( mesh, vert, &elDim, &elInd ) );
	 TestTrue( elDim == MT_VERTEX );
	 TestTrue( elInd == inc[inc_i] );
      }

      Mesh_GetIncidence( mesh, MT_VOLUME, e_i, MT_EDGE, &nInc, &inc );
      for( inc_i = 0; inc_i < nInc; inc_i++ ) {
	 Mesh_GetIncidence( mesh, MT_EDGE, inc[inc_i], MT_VERTEX, &nEdgeInc, &edgeInc );
	 point[0] = point[1] = point[2] = 0.0;
	 for( inc_j = 0; inc_j < nEdgeInc; inc_j++ ) {
	    point[0] += Mesh_GetVertex( mesh, edgeInc[inc_j] )[0];
	    point[1] += Mesh_GetVertex( mesh, edgeInc[inc_j] )[1];
	    point[2] += Mesh_GetVertex( mesh, edgeInc[inc_j] )[2];
	 }
	 point[0] /= (double)nEdgeInc;
	 point[1] /= (double)nEdgeInc;
	 point[2] /= (double)nEdgeInc;
	 TestTrue( Mesh_Search( mesh, point, &elDim, &elInd ) );
	 TestTrue( elDim == MT_EDGE );
	 TestTrue( elInd == inc[inc_i] );

	 Mesh_GetIncidence( mesh, MT_VOLUME, e_i, MT_FACE, &nInc, &inc );
	 for( inc_i = 0; inc_i < nInc; inc_i++ ) {
	    Mesh_GetIncidence( mesh, MT_FACE, inc[inc_i], MT_VERTEX, &nFaceInc, &faceInc );
	    point[0] = point[1] = point[2] = 0.0;
	    for( inc_j = 0; inc_j < nFaceInc; inc_j++ ) {
	       point[0] += Mesh_GetVertex( mesh, faceInc[inc_j] )[0];
	       point[1] += Mesh_GetVertex( mesh, faceInc[inc_j] )[1];
	       point[2] += Mesh_GetVertex( mesh, faceInc[inc_j] )[2];
	    }
	    point[0] /= (double)nFaceInc;
	    point[1] /= (double)nFaceInc;
	    point[2] /= (double)nFaceInc;
	    TestTrue( Mesh_Search( mesh, point, &elDim, &elInd ) );
	    TestTrue( elDim == MT_FACE );
	    TestTrue( elInd == inc[inc_i] );

	    Mesh_GetIncidence( mesh, MT_VOLUME, e_i, MT_VERTEX, &nInc, &inc );
	    point[0] = point[1] = point[2] = 0.0;
	    for( inc_i = 0; inc_i < nInc; inc_i++ ) {
	       point[0] += Mesh_GetVertex( mesh, inc[inc_i] )[0];
	       point[1] += Mesh_GetVertex( mesh, inc[inc_i] )[1];
	       point[2] += Mesh_GetVertex( mesh, inc[inc_i] )[2];
	    }
	    point[0] /= (double)nInc;
	    point[1] /= (double)nInc;
	    point[2] /= (double)nInc;
	    TestTrue( Mesh_Search( mesh, point, &elDim, &elInd ) );
	    TestTrue( elDim == MT_VOLUME );
	    TestTrue( elInd == e_i );
	 }
      }
   }

done:
   FreeObject( gen );
   FreeObject( mesh );
}
TestEnd


#define nTests 2
TestSuite_Test tests[nTests] = {{"nearest vertex", testNearVert}, 
				{"element search 3D", testElSearch3D}};


#include "Base/Foundation/TestEnd.h"
