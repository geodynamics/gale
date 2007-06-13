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
   Stream_Enable( Journal_GetTypedStream( Debug_Type ), False );
   Stream_Enable( Journal_GetTypedStream( Info_Type ), False );
   Stream_Enable( Journal_GetTypedStream( Dump_Type ), False );
   Stream_Enable( Journal_GetTypedStream( Error_Type ), False );
}

void testTeardown() {
   DiscretisationMesh_Finalise();
   Base_Finalise();
}

TestBegin( NearVert ) {
   CartesianGenerator* gen;
   Mesh* mesh;
   int nDims;
   int sizes[3];
   double minCrd[3];
   double maxCrd[3];
   int nRanks;
   int nInc, *inc;
   double* vert;
   int e_i, inc_i;

   insist( MPI_Comm_size( MPI_COMM_WORLD, &nRanks ), == MPI_SUCCESS );
   sizes[0] = sizes[1] = sizes[2] = 4 * nRanks;
   minCrd[0] = minCrd[1] = minCrd[2] = 0.0;
   maxCrd[0] = maxCrd[1] = maxCrd[2] = (double)nRanks;

   nDims = 1;
   gen = CartesianGenerator_New( "" );
   MeshGenerator_SetDimSize( gen, nDims );
   CartesianGenerator_SetShadowDepth( gen, 1 );
   CartesianGenerator_SetTopologyParams( gen, sizes, 0, NULL, NULL );
   CartesianGenerator_SetGeometryParams( gen, minCrd, maxCrd );
   mesh = Mesh_New( "" );
   Mesh_SetGenerator( mesh, gen );
   Stg_Component_Build( mesh, NULL, False );

   for( e_i = 0; e_i < Mesh_GetDomainSize( mesh, nDims ); e_i++ ) {
      Mesh_GetIncidence( mesh, nDims, e_i, MT_VERTEX, &nInc, &inc );
      for( inc_i = 0; inc_i < nInc; inc_i++ ) {
	 vert = Mesh_GetVertex( mesh, inc[inc_i] );
	 if( !Mesh_NearestVertex( mesh, vert ) == inc[inc_i] ) break;
      }
   }
   TestTrue( e_i == Mesh_GetDomainSize( mesh, nDims ) );

   nDims = 2;
   MeshGenerator_SetDimSize( gen, nDims );
   CartesianGenerator_SetShadowDepth( gen, 1 );
   CartesianGenerator_SetTopologyParams( gen, sizes, 0, NULL, NULL );
   CartesianGenerator_SetGeometryParams( gen, minCrd, maxCrd );
   Mesh_SetGenerator( mesh, gen );
   Stg_Component_Build( mesh, NULL, True );

   for( e_i = 0; e_i < Mesh_GetDomainSize( mesh, nDims ); e_i++ ) {
      Mesh_GetIncidence( mesh, nDims, e_i, MT_VERTEX, &nInc, &inc );
      for( inc_i = 0; inc_i < nInc; inc_i++ ) {
	 vert = Mesh_GetVertex( mesh, inc[inc_i] );
	 if( !Mesh_NearestVertex( mesh, vert ) == inc[inc_i] ) break;
      }
   }
   TestTrue( e_i == Mesh_GetDomainSize( mesh, nDims ) );

   nDims = 3;
   MeshGenerator_SetDimSize( gen, nDims );
   CartesianGenerator_SetShadowDepth( gen, 1 );
   CartesianGenerator_SetTopologyParams( gen, sizes, 0, NULL, NULL );
   CartesianGenerator_SetGeometryParams( gen, minCrd, maxCrd );
   Mesh_SetGenerator( mesh, gen );
   Stg_Component_Build( mesh, NULL, True );

   for( e_i = 0; e_i < Mesh_GetDomainSize( mesh, nDims ); e_i++ ) {
      Mesh_GetIncidence( mesh, nDims, e_i, MT_VERTEX, &nInc, &inc );
      for( inc_i = 0; inc_i < nInc; inc_i++ ) {
	 vert = Mesh_GetVertex( mesh, inc[inc_i] );
	 if( !Mesh_NearestVertex( mesh, vert ) == inc[inc_i] ) break;
      }
   }
   TestTrue( e_i == Mesh_GetDomainSize( mesh, nDims ) );

done:
   FreeObject( gen );
   FreeObject( mesh );
}
TestEnd

TestBegin( Search ) {
   CartesianGenerator* gen;
   Mesh* mesh;
   int nDims;
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
   int e_i, inc_i, inc_j, d_i;

   insist( MPI_Comm_size( MPI_COMM_WORLD, &nRanks ), == MPI_SUCCESS );
   sizes[0] = sizes[1] = sizes[2] = 4 * nRanks;
   minCrd[0] = minCrd[1] = minCrd[2] = 0.0;
   maxCrd[0] = maxCrd[1] = maxCrd[2] = (double)nRanks;

   nDims = 1;
   gen = CartesianGenerator_New( "" );
   MeshGenerator_SetDimSize( gen, nDims );
   MeshGenerator_SetFullIncidence( gen );
   CartesianGenerator_SetShadowDepth( gen, 1 );
   CartesianGenerator_SetTopologyParams( gen, sizes, 0, NULL, NULL );
   CartesianGenerator_SetGeometryParams( gen, minCrd, maxCrd );
   mesh = Mesh_New( "" );
   Mesh_SetGenerator( mesh, gen );
   Stg_Component_Build( mesh, NULL, False );

   for( e_i = 0; e_i < Mesh_GetDomainSize( mesh, nDims ); e_i++ ) {
      Mesh_GetIncidence( mesh, nDims, e_i, MT_VERTEX, &nInc, &inc );
      for( inc_i = 0; inc_i < nInc; inc_i++ ) {
	 vert = Mesh_GetVertex( mesh, inc[inc_i] );
	 if( !Mesh_Search( mesh, vert, &elDim, &elInd ) ) break;
	 if( elDim != MT_VERTEX ) break;
	 if( elInd != inc[inc_i] ) break;
      }

      Mesh_GetIncidence( mesh, nDims, e_i, MT_VERTEX, &nInc, &inc );
      for( d_i = 0; d_i < nDims; d_i++ )
	 point[d_i] = 0.0;
      for( inc_i = 0; inc_i < nInc; inc_i++ ) {
	 for( d_i = 0; d_i < nDims; d_i++ )
	    point[d_i] += Mesh_GetVertex( mesh, inc[inc_i] )[d_i];
      }
      for( d_i = 0; d_i < nDims; d_i++ )
	 point[d_i] /= (double)nInc;
      if( !Mesh_Search( mesh, point, &elDim, &elInd ) ) break;
      if( elDim != nDims ) break;
      if( elInd != e_i ) break;
   }
   TestTrue( e_i == Mesh_GetDomainSize( mesh, nDims ) );

   nDims = 2;
   MeshGenerator_SetDimSize( gen, nDims );
   MeshGenerator_SetFullIncidence( gen );
   CartesianGenerator_SetShadowDepth( gen, 1 );
   CartesianGenerator_SetTopologyParams( gen, sizes, 0, NULL, NULL );
   CartesianGenerator_SetGeometryParams( gen, minCrd, maxCrd );
   Mesh_SetGenerator( mesh, gen );
   Stg_Component_Build( mesh, NULL, True );

   for( e_i = 0; e_i < Mesh_GetDomainSize( mesh, nDims ); e_i++ ) {
      Mesh_GetIncidence( mesh, nDims, e_i, MT_VERTEX, &nInc, &inc );
      for( inc_i = 0; inc_i < nInc; inc_i++ ) {
	 vert = Mesh_GetVertex( mesh, inc[inc_i] );
	 if( !Mesh_Search( mesh, vert, &elDim, &elInd ) ) break;
	 if( elDim != MT_VERTEX ) break;
	 if( elInd != inc[inc_i] ) break;
      }

      Mesh_GetIncidence( mesh, nDims, e_i, MT_EDGE, &nInc, &inc );
      for( inc_i = 0; inc_i < nInc; inc_i++ ) {
	 Mesh_GetIncidence( mesh, MT_EDGE, inc[inc_i], MT_VERTEX, &nEdgeInc, &edgeInc );
	 for( d_i = 0; d_i < nDims; d_i++ )
	    point[d_i] = 0.0;
	 for( inc_j = 0; inc_j < nEdgeInc; inc_j++ ) {
	    for( d_i = 0; d_i < nDims; d_i++ )
	       point[d_i] += Mesh_GetVertex( mesh, edgeInc[inc_j] )[d_i];
	 }
	 for( d_i = 0; d_i < nDims; d_i++ )
	    point[d_i] /= (double)nEdgeInc;
	 if( !Mesh_Search( mesh, point, &elDim, &elInd ) ) break;
	 if( elDim != MT_EDGE ) break;
	 if( elInd != inc[inc_i] ) break;
      }

      Mesh_GetIncidence( mesh, nDims, e_i, MT_VERTEX, &nInc, &inc );
      for( d_i = 0; d_i < nDims; d_i++ )
	 point[d_i] = 0.0;
      for( inc_i = 0; inc_i < nInc; inc_i++ ) {
	 for( d_i = 0; d_i < nDims; d_i++ )
	    point[d_i] += Mesh_GetVertex( mesh, inc[inc_i] )[d_i];
      }
      for( d_i = 0; d_i < nDims; d_i++ )
	 point[d_i] /= (double)nInc;
      if( !Mesh_Search( mesh, point, &elDim, &elInd ) ) break;
      if( elDim != nDims ) break;
      if( elInd != e_i ) break;
   }
   TestTrue( e_i == Mesh_GetDomainSize( mesh, nDims ) );

   nDims = 3;
   MeshGenerator_SetDimSize( gen, nDims );
   MeshGenerator_SetFullIncidence( gen );
   CartesianGenerator_SetShadowDepth( gen, 1 );
   CartesianGenerator_SetTopologyParams( gen, sizes, 0, NULL, NULL );
   CartesianGenerator_SetGeometryParams( gen, minCrd, maxCrd );
   Mesh_SetGenerator( mesh, gen );
   Stg_Component_Build( mesh, NULL, True );

   for( e_i = 0; e_i < Mesh_GetDomainSize( mesh, nDims ); e_i++ ) {
      Mesh_GetIncidence( mesh, nDims, e_i, MT_VERTEX, &nInc, &inc );
      for( inc_i = 0; inc_i < nInc; inc_i++ ) {
	 vert = Mesh_GetVertex( mesh, inc[inc_i] );
	 if( !Mesh_Search( mesh, vert, &elDim, &elInd ) ) break;
	 if( elDim != MT_VERTEX ) break;
	 if( elInd != inc[inc_i] ) break;
      }

      Mesh_GetIncidence( mesh, nDims, e_i, MT_EDGE, &nInc, &inc );
      for( inc_i = 0; inc_i < nInc; inc_i++ ) {
	 Mesh_GetIncidence( mesh, MT_EDGE, inc[inc_i], MT_VERTEX, &nEdgeInc, &edgeInc );
	 for( d_i = 0; d_i < nDims; d_i++ )
	    point[d_i] = 0.0;
	 for( inc_j = 0; inc_j < nEdgeInc; inc_j++ ) {
	    for( d_i = 0; d_i < nDims; d_i++ )
	       point[d_i] += Mesh_GetVertex( mesh, edgeInc[inc_j] )[d_i];
	 }
	 for( d_i = 0; d_i < nDims; d_i++ )
	    point[d_i] /= (double)nEdgeInc;
	 if( !Mesh_Search( mesh, point, &elDim, &elInd ) ) break;
	 if( elDim != MT_EDGE ) break;
	 if( elInd != inc[inc_i] ) break;
      }

      Mesh_GetIncidence( mesh, nDims, e_i, MT_FACE, &nInc, &inc );
      for( inc_i = 0; inc_i < nInc; inc_i++ ) {
	 Mesh_GetIncidence( mesh, MT_FACE, inc[inc_i], MT_VERTEX, &nFaceInc, &faceInc );
	 for( d_i = 0; d_i < nDims; d_i++ )
	    point[d_i] = 0.0;
	 for( inc_j = 0; inc_j < nFaceInc; inc_j++ ) {
	    for( d_i = 0; d_i < nDims; d_i++ )
	       point[d_i] += Mesh_GetVertex( mesh, faceInc[inc_j] )[d_i];
	 }
	 for( d_i = 0; d_i < nDims; d_i++ )
	    point[d_i] /= (double)nFaceInc;
	 if( !Mesh_Search( mesh, point, &elDim, &elInd ) ) break;
	 if( elDim != MT_FACE ) break;
	 if( elInd != inc[inc_i] ) break;
      }

      Mesh_GetIncidence( mesh, nDims, e_i, MT_VERTEX, &nInc, &inc );
      for( d_i = 0; d_i < nDims; d_i++ )
	 point[d_i] = 0.0;
      for( inc_i = 0; inc_i < nInc; inc_i++ ) {
	 for( d_i = 0; d_i < nDims; d_i++ )
	    point[d_i] += Mesh_GetVertex( mesh, inc[inc_i] )[d_i];
      }
      for( d_i = 0; d_i < nDims; d_i++ )
	 point[d_i] /= (double)nInc;
      if( !Mesh_Search( mesh, point, &elDim, &elInd ) ) break;
      if( elDim != nDims ) break;
      if( elInd != e_i ) break;
   }
   TestTrue( e_i == Mesh_GetDomainSize( mesh, nDims ) );

done:
   FreeObject( gen );
   FreeObject( mesh );
}
TestEnd


#define nTests 2
TestSuite_Test tests[nTests] = {{"nearest vertex", testNearVert}, 
				{"search", testSearch}};


#include "Base/Foundation/TestEnd.h"
