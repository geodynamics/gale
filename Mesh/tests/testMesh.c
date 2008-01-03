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
#include <StGermain/StGermain.h>
#include <StgDomain/Mesh/Mesh.h>

#include "StGermain/Base/Foundation/TestBegin.h"


int findOwner( Mesh* mesh, int vert ) {
   IArray* inc;
   int lowest, cur;
   int nDims;
   int ii;

   inc = IArray_New();

   nDims = Mesh_GetDimSize( mesh );
   Mesh_GetIncidence( mesh, 0, vert, nDims, inc );
   lowest = Mesh_DomainToGlobal( mesh, nDims, IArray_GetPtr( inc )[0] );
   for( ii = 1; ii < IArray_GetSize( inc ); ii++ ) {
      cur = Mesh_DomainToGlobal( mesh, nDims, IArray_GetPtr( inc )[ii] );
      if( cur < lowest )
	 lowest = cur;
   }

   NewClass_Delete( inc );

   insist( Mesh_GlobalToDomain( mesh, nDims, lowest, &lowest ), == True );
   return lowest;
}


void testSetup( int* argc, char** argv[] ) {
   StGermain_Init( argc, argv );
   StgDomainMesh_Init( argc, argv );
   Stream_Enable( Journal_GetTypedStream( Debug_Type ), False );
   Stream_Enable( Journal_GetTypedStream( Info_Type ), False );
   Stream_Enable( Journal_GetTypedStream( Dump_Type ), False );
   Stream_Enable( Journal_GetTypedStream( Error_Type ), False );
}

void testTeardown() {
   StgDomainMesh_Finalise();
   StGermain_Finalise();
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
   IArray* incArray;
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
   incArray = IArray_New();

   for( e_i = 0; e_i < Mesh_GetDomainSize( mesh, nDims ); e_i++ ) {
      Mesh_GetIncidence( mesh, nDims, e_i, MT_VERTEX, incArray );
      nInc = IArray_GetSize( incArray );
      inc = IArray_GetPtr( incArray );
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
      Mesh_GetIncidence( mesh, nDims, e_i, MT_VERTEX, incArray );
      nInc = IArray_GetSize( incArray );
      inc = IArray_GetPtr( incArray );
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
      Mesh_GetIncidence( mesh, nDims, e_i, MT_VERTEX, incArray );
      nInc = IArray_GetSize( incArray );
      inc = IArray_GetPtr( incArray );
      for( inc_i = 0; inc_i < nInc; inc_i++ ) {
	 vert = Mesh_GetVertex( mesh, inc[inc_i] );
	 if( !Mesh_NearestVertex( mesh, vert ) == inc[inc_i] ) break;
      }
   }
   TestTrue( e_i == Mesh_GetDomainSize( mesh, nDims ) );

   NewClass_Delete( incArray );

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
   int el;
   int ii;

   insist( MPI_Comm_size( MPI_COMM_WORLD, &nRanks ), == MPI_SUCCESS );
   sizes[0] = sizes[1] = sizes[2] = 2 * nRanks;
   minCrd[0] = minCrd[1] = minCrd[2] = 0.0;
   maxCrd[0] = maxCrd[1] = maxCrd[2] = (double)nRanks;

   nDims = 3;
   gen = CartesianGenerator_New( "" );
   MeshGenerator_SetDimSize( gen, nDims );
   CartesianGenerator_SetShadowDepth( gen, 1 );
   CartesianGenerator_SetTopologyParams( gen, sizes, 0, NULL, NULL );
   CartesianGenerator_SetGeometryParams( gen, minCrd, maxCrd );
   mesh = Mesh_New( "" );
   Mesh_SetGenerator( mesh, gen );
   Stg_Component_Build( mesh, NULL, False );

   for( ii = 0; ii < Mesh_GetLocalSize( mesh, 0 ); ii++ ) {
      if( !Mesh_SearchElements( mesh, Mesh_GetVertex( mesh, ii ), &el ) )
	 break;
      if( el != findOwner( mesh, ii ) )
	 break;
   }
   TestTrue( ii == Mesh_GetLocalSize( mesh, 0 ) );


done:
   FreeObject( gen );
   FreeObject( mesh );
}
TestEnd


#define nTests 2
TestSuite_Test tests[nTests] = {{"nearest vertex", testNearVert}, 
				{"search", testSearch}};


#include "StGermain/Base/Foundation/TestEnd.h"
