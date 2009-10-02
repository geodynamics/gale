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
** $Id: testTrilinearElementType.c 896 2007-07-03 04:46:51Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include "StgFEM/Discretisation/Discretisation.h"

#include "StGermain/Base/Foundation/TestBegin.h"


FeMesh* buildMesh() {
   CartesianGenerator* gen;
   int nRanks;
   unsigned sizes[3];
   double minCrd[3];
   double maxCrd[3];
   FeMesh* mesh;

   insist( MPI_Comm_size( MPI_COMM_WORLD, &nRanks ), == MPI_SUCCESS );
   sizes[0] = sizes[1] = sizes[2] = nRanks * 4;
   minCrd[0] = minCrd[1] = minCrd[2] = 0.0;
   maxCrd[0] = minCrd[1] = minCrd[2] = (double)nRanks;

   gen = CartesianGenerator_New( "" );
   MeshGenerator_SetDimSize( gen, 3 );
   CartesianGenerator_SetShadowDepth( gen, 1 );
   CartesianGenerator_SetTopologyParams( gen, sizes, 0, NULL, NULL );
   CartesianGenerator_SetGeometryParams( gen, minCrd, maxCrd );

   mesh = FeMesh_New( "" );
   Mesh_SetGenerator( mesh, gen );
   FeMesh_SetElementFamily( mesh, "linear" );
   Stg_Component_Build( mesh, NULL, False );

   return mesh;
}


void testSetup( int* argc, char** argv[] ) {
   StGermain_Init( argc, argv );
   StgDomain_Init( argc, argv );
   StgFEM_Discretisation_Init( argc, argv );
   Stream_Enable( Journal_GetTypedStream( Debug_Type ), False );
   Stream_Enable( Journal_GetTypedStream( Info_Type ), False );
   Stream_Enable( Journal_GetTypedStream( Dump_Type ), False );
   Stream_Enable( Journal_GetTypedStream( Error_Type ), False );
}

void testTeardown() {
   StgFEM_Discretisation_Finalise();
   StgDomain_Finalise();
   StGermain_Finalise();
}

TestBegin( Shape ) {
   FeMesh* mesh;
   int nEls, nVerts, nDims;
   const int *verts;
   double* vert;
   double lCrd[3], basis[8];
   IArray* inc;
   int e_i, v_i, v_j;

   mesh = buildMesh();
   TestTrue( mesh );

   nDims = Mesh_GetDimSize( mesh );
   nEls = Mesh_GetDomainSize( mesh, nDims );
   inc = IArray_New();
   for( e_i = 0; e_i < nEls; e_i++ ) {
      Mesh_GetIncidence( mesh, nDims, e_i, 0, inc );
      nVerts = IArray_GetSize( inc );
      verts = IArray_GetPtr( inc );
      for( v_i = 0; v_i < nVerts; v_i++ ) {
	 vert = Mesh_GetVertex( mesh, verts[v_i] );
	 FeMesh_CoordGlobalToLocal( mesh, e_i, vert, lCrd );
	 FeMesh_EvalBasis( mesh, e_i, lCrd, basis );
	 for( v_j = 0; v_j < nVerts; v_j++ ) {
	    if( (v_i == v_j && !Num_Approx( basis[v_j], 1.0 )) || 
		(v_i != v_j && !Num_Approx( basis[v_j], 0.0 )) )
	    {
	       break;
	    }
	 }
	 if( v_j < nVerts )
	    break;
      }
      if( v_i < nVerts )
	 break;
   }
   TestTrue( e_i == nEls );

   NewClass_Delete( inc );

  done:
   FreeObject( mesh );
}
TestEnd


#define nTests 1
TestSuite_Test tests[nTests] = {{"shape functions", testShape}};


#include "StGermain/Base/Foundation/TestEnd.h"
