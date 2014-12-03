/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street, Melbourne, 3053, Australia.
**
** Authors:
**   Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**   Patrick D. Sunter, Software Engineer, VPAC. (pds@vpac.org)
**   Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**   Siew-Ching Tan, Software Engineer, VPAC. (siew@vpac.org)
**   Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
**   Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
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
** Role:
**   Tests the TrilinearElementTypeSuite
**
** $Id: testTrilinearElementType.c 3462 2006-02-19 06:53:24Z WalterLandry $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>

#include "pcu/pcu.h"
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/Discretisation/Discretisation.h>
#include "TrilinearElementTypeSuite.h"

typedef struct {
} TrilinearElementTypeSuiteData;

FeMesh* buildMesh() {
   CartesianGenerator*	gen;
   int						nRanks;
   unsigned					sizes[3];
   double					minCrd[3];
   double					maxCrd[3];
   FeMesh*					mesh;

   insist( MPI_Comm_size( MPI_COMM_WORLD, &nRanks ), == MPI_SUCCESS );
   sizes[0] = sizes[1] = sizes[2] = nRanks * 4;
   minCrd[0] = minCrd[1] = minCrd[2] = 0.0;
   maxCrd[0] = minCrd[1] = minCrd[2] = (double)nRanks;

   gen = CartesianGenerator_New( "", NULL );
   MeshGenerator_SetDimSize( gen, 3 );
   CartesianGenerator_SetShadowDepth( gen, 1 );
   CartesianGenerator_SetTopologyParams( gen, sizes, 0, NULL, NULL );
   CartesianGenerator_SetGeometryParams( gen, minCrd, maxCrd );

   mesh = FeMesh_New( "", NULL );
   Mesh_SetGenerator( mesh, gen );
   FeMesh_SetElementFamily( mesh, "linear" );
   Stg_Component_Build( mesh, NULL, False );

   return mesh;
}

void TrilinearElementTypeSuite_Setup( TrilinearElementTypeSuiteData* data ) {
	Journal_Enable_AllTypedStream( False );

	//Stream_RedirectAllToFile( "TrilinearElementTypeSuite" );
}

void TrilinearElementTypeSuite_Teardown( TrilinearElementTypeSuiteData* data ) {
	//Stream_PurgeAllRedirectedFiles();
}

void TrilinearElementTypeSuite_TestShape( TrilinearElementTypeSuiteData* data ) {
   FeMesh*		mesh = NULL;
   int			nEls, nVerts, nDims;
   const int	*verts;
   double*		vert = NULL;
   double		lCrd[3] = { 0.0, 0.0, 0.0 };
	double		basis[8] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
   IArray*		inc;
   int			e_i, v_i, v_j;

   mesh = buildMesh();
   pcu_check_true( mesh );
   Stg_Component_Initialise( mesh, data, True );

   nDims = Mesh_GetDimSize( mesh );
   nEls = Mesh_GetDomainSize( mesh, (MeshTopology_Dim)nDims );
   inc = IArray_New();

   for( e_i = 0; e_i < nEls; e_i++ ) {
     Mesh_GetIncidence( mesh, (MeshTopology_Dim)nDims, e_i, (MeshTopology_Dim)0, inc );
      nVerts = IArray_GetSize( inc );
      verts = IArray_GetPtr( inc );

      for( v_i = 0; v_i < nVerts; v_i++ ) {
			vert = Mesh_GetVertex( mesh, verts[v_i] );
			FeMesh_CoordGlobalToLocal( mesh, e_i, vert, lCrd );
			FeMesh_EvalBasis( mesh, e_i, lCrd, basis );

			for( v_j = 0; v_j < nVerts; v_j++ ) {
				if( (v_i == v_j && !Num_Approx( basis[v_j], 1.0 )) || (v_i != v_j && !Num_Approx( basis[v_j], 0.0 )) ) {
					break;
				}
			}
			if( v_j < nVerts )
				break;
      }
      if( v_i < nVerts )
			break;
   }
   pcu_check_true( e_i == nEls );

   NewClass_Delete( inc );

   _Stg_Component_Delete( mesh );
}

void TrilinearElementTypeSuite( pcu_suite_t* suite ) {
   pcu_suite_setData( suite, TrilinearElementTypeSuiteData );
   pcu_suite_setFixtures( suite, TrilinearElementTypeSuite_Setup, TrilinearElementTypeSuite_Teardown );
   pcu_suite_addTest( suite, TrilinearElementTypeSuite_TestShape );
}


