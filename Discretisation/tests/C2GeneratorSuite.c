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
**   Tests the C2GeneratorSuite
**
** $Id: testTemplate.c 3462 2006-02-19 06:53:24Z WalterLandry $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>

#include "pcu/pcu.h"
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/Discretisation/Discretisation.h>
#include "C2GeneratorSuite.h"

typedef struct {
	Mesh*	mesh;
} C2GeneratorSuiteData;

void C2GeneratorSuite_Setup( C2GeneratorSuiteData* data ) {
	C2Generator* gen;
	int nRanks;
	unsigned sizes[3];
	double minCrd[3];
	double maxCrd[3];
	int rank;

	Journal_Enable_AllTypedStream( False );

	insist( MPI_Comm_size( MPI_COMM_WORLD, &nRanks ), == MPI_SUCCESS );
	sizes[0] = sizes[1] = sizes[2] = nRanks * 4;
	minCrd[0] = minCrd[1] = minCrd[2] = 0.0;
	maxCrd[0] = minCrd[1] = minCrd[2] = (double)nRanks;

	gen = C2Generator_New( "", NULL );
	MeshGenerator_SetDimSize( gen, 3 );
	CartesianGenerator_SetShadowDepth( gen, 1 );
	C2Generator_SetTopologyParams( gen, sizes, 0, NULL, NULL );
	CartesianGenerator_SetGeometryParams( gen, minCrd, maxCrd );
   Stg_Component_Build( gen, NULL, False );
   Stg_Component_Initialise( gen, NULL, False );

	data->mesh = Mesh_New( "", NULL );
	CartesianGenerator_Generate( gen, data->mesh, NULL );
   /* need 2 lines below to setup mesh */
   Stg_Component_Build( data->mesh, NULL, False );
   Stg_Component_Initialise( data->mesh, NULL, False );
}

void C2GeneratorSuite_Teardown( C2GeneratorSuiteData* data ) {
	Stg_Component_Destroy( data->mesh, NULL, True );
}

void C2GeneratorSuite_TestElementVertexInc( C2GeneratorSuiteData* data ) {
	unsigned	dim		= Mesh_GetDimSize( data->mesh );
	Sync*		elSync		= (Sync*)IGraph_GetDomain( (IGraph*)data->mesh->topo, dim );
	Sync*		vertSync	= (Sync*)IGraph_GetDomain( (IGraph*)data->mesh->topo, MT_VERTEX );
	Grid*		elGrid		= *(Grid**)Mesh_GetExtension( data->mesh, Grid**, "elementGrid" );
	Grid*		vertGrid	= *(Grid**)Mesh_GetExtension( data->mesh, Grid**, "vertexGrid" );
	IArray*		inc		= IArray_New();
	unsigned	vert;
	unsigned*	incVerts;
	unsigned	nIncVerts;
	unsigned	el_i;
	unsigned	gEl;
	unsigned	dimInds[3];
	unsigned	gNode0, gNode1, gNode2;
	int		checkNodes;

	for( el_i = 0; el_i < Mesh_GetLocalSize( data->mesh, dim ); el_i++ ) {
		gEl = Sync_DomainToGlobal( elSync, el_i );
		Grid_Lift( elGrid, gEl, dimInds );
		dimInds[0] *= 2; dimInds[1] *= 2; dimInds[2] *= 2;

		MeshTopology_GetIncidence( (IGraph*)data->mesh->topo, dim, el_i, MT_VERTEX, inc );
		nIncVerts = IArray_GetSize( inc );
		incVerts = IArray_GetPtr( inc );

		pcu_check_true( nIncVerts == 27 );

		gNode0 = Grid_Project( vertGrid, dimInds );
		gNode1 = Grid_Project( vertGrid, dimInds );
		gNode2 = Sync_DomainToGlobal( vertSync, incVerts[0] );
		checkNodes = (gNode0 == gNode1) && (gNode1 == gNode2);
		pcu_check_true( checkNodes );

		gNode0++;
		dimInds[0]++;
		gNode1 = Grid_Project( vertGrid, dimInds );
		gNode2 = Sync_DomainToGlobal( vertSync, incVerts[1] );
		checkNodes = (gNode0 == gNode1) && (gNode1 == gNode2);
		pcu_check_true( checkNodes );

		gNode0++;
		dimInds[0]++;
		gNode1 = Grid_Project( vertGrid, dimInds );
		gNode2 = Sync_DomainToGlobal( vertSync, incVerts[2] );
		checkNodes = (gNode0 == gNode1) && (gNode1 == gNode2);
		pcu_check_true( checkNodes );

		gNode0 -= 2; gNode0 += vertGrid->sizes[0];
		dimInds[0] -= 2; dimInds[1]++;
		gNode1 = Grid_Project( vertGrid, dimInds );
		gNode2 = Sync_DomainToGlobal( vertSync, incVerts[3] );
		checkNodes = (gNode0 == gNode1) && (gNode1 == gNode2);
		pcu_check_true( checkNodes );

		gNode0++;
		dimInds[0]++;
		gNode1 = Grid_Project( vertGrid, dimInds );
		gNode2 = Sync_DomainToGlobal( vertSync, incVerts[4] );
		checkNodes = (gNode0 == gNode1) && (gNode1 == gNode2);
		pcu_check_true( checkNodes );

		gNode0++;
		dimInds[0]++;
		gNode1 = Grid_Project( vertGrid, dimInds );
		gNode2 = Sync_DomainToGlobal( vertSync, incVerts[5] );
		checkNodes = (gNode0 == gNode1) && (gNode1 == gNode2);
		pcu_check_true( checkNodes );

		gNode0 -= 2; gNode0 += vertGrid->sizes[0];
		dimInds[0] -= 2; dimInds[1]++;
		gNode1 = Grid_Project( vertGrid, dimInds );
		gNode2 = Sync_DomainToGlobal( vertSync, incVerts[6] );
		checkNodes = (gNode0 == gNode1) && (gNode1 == gNode2);
		pcu_check_true( checkNodes );

		gNode0++;
		dimInds[0]++;
		gNode1 = Grid_Project( vertGrid, dimInds );
		gNode2 = Sync_DomainToGlobal( vertSync, incVerts[7] );
		checkNodes = (gNode0 == gNode1) && (gNode1 == gNode2);
		pcu_check_true( checkNodes );

		gNode0++;
		dimInds[0]++;
		gNode1 = Grid_Project( vertGrid, dimInds );
		gNode2 = Sync_DomainToGlobal( vertSync, incVerts[8] );
		checkNodes = (gNode0 == gNode1) && (gNode1 == gNode2);
		pcu_check_true( checkNodes );

		gNode0 -= 2; gNode0 -= 2 * vertGrid->sizes[0]; gNode0 += vertGrid->sizes[0] * vertGrid->sizes[1];
		dimInds[0] -= 2; dimInds[1] -= 2; dimInds[2]++;
		gNode1 = Grid_Project( vertGrid, dimInds );
		gNode2 = Sync_DomainToGlobal( vertSync, incVerts[9] );
		checkNodes = (gNode0 == gNode1) && (gNode1 == gNode2);
		pcu_check_true( checkNodes );

		gNode0++;
		dimInds[0]++;
		gNode1 = Grid_Project( vertGrid, dimInds );
		gNode2 = Sync_DomainToGlobal( vertSync, incVerts[10] );
		checkNodes = (gNode0 == gNode1) && (gNode1 == gNode2);
		pcu_check_true( checkNodes );

		gNode0++;
		dimInds[0]++;
		gNode1 = Grid_Project( vertGrid, dimInds );
		gNode2 = Sync_DomainToGlobal( vertSync, incVerts[11] );
		checkNodes = (gNode0 == gNode1) && (gNode1 == gNode2);
		pcu_check_true( checkNodes );

		gNode0 -= 2; gNode0 += vertGrid->sizes[0];
		dimInds[0] -= 2; dimInds[1]++;
		gNode1 = Grid_Project( vertGrid, dimInds );
		gNode2 = Sync_DomainToGlobal( vertSync, incVerts[12] );
		checkNodes = (gNode0 == gNode1) && (gNode1 == gNode2);
		pcu_check_true( checkNodes );

		gNode0++;
		dimInds[0]++;
		gNode1 = Grid_Project( vertGrid, dimInds );
		gNode2 = Sync_DomainToGlobal( vertSync, incVerts[13] );
		checkNodes = (gNode0 == gNode1) && (gNode1 == gNode2);
		pcu_check_true( checkNodes );

		gNode0++;
		dimInds[0]++;
		gNode1 = Grid_Project( vertGrid, dimInds );
		gNode2 = Sync_DomainToGlobal( vertSync, incVerts[14] );
		checkNodes = (gNode0 == gNode1) && (gNode1 == gNode2);
		pcu_check_true( checkNodes );

		gNode0 -= 2; gNode0 += vertGrid->sizes[0];
		dimInds[0] -= 2; dimInds[1]++;
		gNode1 = Grid_Project( vertGrid, dimInds );
		gNode2 = Sync_DomainToGlobal( vertSync, incVerts[15] );
		checkNodes = (gNode0 == gNode1) && (gNode1 == gNode2);
		pcu_check_true( checkNodes );

		gNode0++;
		dimInds[0]++;
		gNode1 = Grid_Project( vertGrid, dimInds );
		gNode2 = Sync_DomainToGlobal( vertSync, incVerts[16] );
		checkNodes = (gNode0 == gNode1) && (gNode1 == gNode2);
		pcu_check_true( checkNodes );

		gNode0++;
		dimInds[0]++;
		gNode1 = Grid_Project( vertGrid, dimInds );
		gNode2 = Sync_DomainToGlobal( vertSync, incVerts[17] );
		checkNodes = (gNode0 == gNode1) && (gNode1 == gNode2);
		pcu_check_true( checkNodes );

		gNode0 -= 2; gNode0 -= 2 * vertGrid->sizes[0]; gNode0 += vertGrid->sizes[0] * vertGrid->sizes[1];
		dimInds[0] -= 2; dimInds[1] -= 2; dimInds[2]++;
		gNode1 = Grid_Project( vertGrid, dimInds );
		gNode2 = Sync_DomainToGlobal( vertSync, incVerts[18] );
		checkNodes = (gNode0 == gNode1) && (gNode1 == gNode2);
		pcu_check_true( checkNodes );

		gNode0++;
		dimInds[0]++;
		gNode1 = Grid_Project( vertGrid, dimInds );
		gNode2 = Sync_DomainToGlobal( vertSync, incVerts[19] );
		checkNodes = (gNode0 == gNode1) && (gNode1 == gNode2);
		pcu_check_true( checkNodes );

		gNode0++;
		dimInds[0]++;
		gNode1 = Grid_Project( vertGrid, dimInds );
		gNode2 = Sync_DomainToGlobal( vertSync, incVerts[20] );
		checkNodes = (gNode0 == gNode1) && (gNode1 == gNode2);
		pcu_check_true( checkNodes );

		gNode0 -= 2; gNode0 += vertGrid->sizes[0];
		dimInds[0] -= 2; dimInds[1]++;
		gNode1 = Grid_Project( vertGrid, dimInds );
		gNode2 = Sync_DomainToGlobal( vertSync, incVerts[21] );
		checkNodes = (gNode0 == gNode1) && (gNode1 == gNode2);
		pcu_check_true( checkNodes );

		gNode0++;
		dimInds[0]++;
		gNode1 = Grid_Project( vertGrid, dimInds );
		gNode2 = Sync_DomainToGlobal( vertSync, incVerts[22] );
		checkNodes = (gNode0 == gNode1) && (gNode1 == gNode2);
		pcu_check_true( checkNodes );

		gNode0++;
		dimInds[0]++;
		gNode1 = Grid_Project( vertGrid, dimInds );
		gNode2 = Sync_DomainToGlobal( vertSync, incVerts[23] );
		checkNodes = (gNode0 == gNode1) && (gNode1 == gNode2);
		pcu_check_true( checkNodes );

		gNode0 -= 2; gNode0 += vertGrid->sizes[0];
		dimInds[0] -= 2; dimInds[1]++;
		gNode1 = Grid_Project( vertGrid, dimInds );
		gNode2 = Sync_DomainToGlobal( vertSync, incVerts[24] );
		checkNodes = (gNode0 == gNode1) && (gNode1 == gNode2);
		pcu_check_true( checkNodes );

		gNode0++;
		dimInds[0]++;
		gNode1 = Grid_Project( vertGrid, dimInds );
		gNode2 = Sync_DomainToGlobal( vertSync, incVerts[25] );
		checkNodes = (gNode0 == gNode1) && (gNode1 == gNode2);
		pcu_check_true( checkNodes );

		gNode0++;
		dimInds[0]++;
		gNode1 = Grid_Project( vertGrid, dimInds );
		gNode2 = Sync_DomainToGlobal( vertSync, incVerts[26] );
		checkNodes = (gNode0 == gNode1) && (gNode1 == gNode2);
		pcu_check_true( checkNodes );
	}

	NewClass_Delete( inc );
}

void C2GeneratorSuite_TestEdgeVertexInc( C2GeneratorSuiteData* data ) {
	Grid*		elGrid		= *(Grid**)Mesh_GetExtension( data->mesh, Grid**, "elementGrid" );
	Grid*		vertGrid	= *(Grid**)Mesh_GetExtension( data->mesh, Grid**, "vertexGrid" );
	Grid*		edgeGrid_0	= Grid_New();
	Grid*		edgeGrid_1	= Grid_New();
	Grid*		edgeGrid_2	= Grid_New();
	Sync*		vertSync	= (Sync*)IGraph_GetDomain( (IGraph*)data->mesh->topo, MT_VERTEX );
	Sync*		edgeSync	= (Sync*)IGraph_GetDomain( (IGraph*)data->mesh->topo, MT_EDGE );
	IArray*		inc		= IArray_New();
	unsigned*	incVerts;
	unsigned	nIncVerts;
	unsigned	dim		= ((IGraph*)data->mesh->topo)->nDims;
	unsigned	sizes[3];
	unsigned	edge_i;
	unsigned	gEdge;
	unsigned	dimInds[3];
	unsigned	gNode0, gNode1, gNode2;
	int		checkNodes;

	sizes[0] = elGrid->sizes[0];
	sizes[1] = elGrid->sizes[1] + 1;
	sizes[2] = elGrid->sizes[2] + 1;
	Grid_SetNumDims( edgeGrid_0, dim );
	Grid_SetSizes( edgeGrid_0, sizes );
	sizes[0] = elGrid->sizes[0] + 1;
	sizes[1] = elGrid->sizes[1];
	sizes[2] = elGrid->sizes[2] + 1;
	Grid_SetNumDims( edgeGrid_1, dim );
	Grid_SetSizes( edgeGrid_1, sizes );
	sizes[0] = elGrid->sizes[0] + 1;
	sizes[1] = elGrid->sizes[1] + 1;
	sizes[2] = elGrid->sizes[2];
	Grid_SetNumDims( edgeGrid_2, dim );
	Grid_SetSizes( edgeGrid_2, sizes );

	for( edge_i = 0; edge_i < Sync_GetNumDomains( edgeSync ); edge_i++ ) {
		gEdge = Sync_DomainToGlobal( edgeSync, edge_i );

		MeshTopology_GetIncidence( (IGraph*)data->mesh->topo, MT_EDGE, edge_i, MT_VERTEX, inc );
		nIncVerts = IArray_GetSize( inc );
		incVerts = IArray_GetPtr( inc );

		pcu_check_true( nIncVerts == 3 );

		if( gEdge < edgeGrid_0->nPoints ) {
			Grid_Lift( edgeGrid_0, gEdge, dimInds );
			dimInds[0] *= 2; dimInds[1] *= 2; dimInds[2] *= 2;

			gNode0 = Grid_Project( vertGrid, dimInds );
			gNode1 = Grid_Project( vertGrid, dimInds );
			gNode2 = Sync_DomainToGlobal( vertSync, incVerts[0] );
			checkNodes = (gNode0 == gNode1) && (gNode1 == gNode2);
			pcu_check_true( checkNodes );

			gNode0++;
			dimInds[0]++;
			gNode1 = Grid_Project( vertGrid, dimInds );
			gNode2 = Sync_DomainToGlobal( vertSync, incVerts[1] );
			checkNodes = (gNode0 == gNode1) && (gNode1 == gNode2);
			pcu_check_true( checkNodes );

			gNode0++;
			dimInds[0]++;
			gNode1 = Grid_Project( vertGrid, dimInds );
			gNode2 = Sync_DomainToGlobal( vertSync, incVerts[2] );
			checkNodes = (gNode0 == gNode1) && (gNode1 == gNode2);
			pcu_check_true( checkNodes );
		}
		else if( gEdge < edgeGrid_0->nPoints + edgeGrid_1->nPoints ) {
			Grid_Lift( edgeGrid_1, gEdge - edgeGrid_0->nPoints, dimInds );
			dimInds[0] *= 2; dimInds[1] *= 2; dimInds[2] *= 2;

			gNode0 = Grid_Project( vertGrid, dimInds );
			gNode1 = Grid_Project( vertGrid, dimInds );
			gNode2 = Sync_DomainToGlobal( vertSync, incVerts[1] );
			checkNodes = (gNode0 == gNode1) && (gNode1 == gNode2);
			pcu_check_true( checkNodes );

			gNode0 += vertGrid->sizes[0];
			dimInds[1]++;
			gNode1 = Grid_Project( vertGrid, dimInds );
			gNode2 = Sync_DomainToGlobal( vertSync, incVerts[1] );
			checkNodes = (gNode0 == gNode1) && (gNode1 == gNode2);
			pcu_check_true( checkNodes );

			gNode0 += vertGrid->sizes[0];
			dimInds[1]++;
			gNode1 = Grid_Project( vertGrid, dimInds );
			gNode2 = Sync_DomainToGlobal( vertSync, incVerts[2] );
			checkNodes = (gNode0 == gNode1) && (gNode1 == gNode2);
			pcu_check_true( checkNodes );
		}
		else if( gEdge < edgeGrid_0->nPoints + edgeGrid_1->nPoints + edgeGrid_2->nPoints ) {
			Grid_Lift( edgeGrid_2, gEdge - edgeGrid_0->nPoints - edgeGrid_1->nPoints, dimInds );
			dimInds[0] *= 2; dimInds[1] *= 2; dimInds[2] *= 2;

			gNode0 = Grid_Project( vertGrid, dimInds );
			gNode1 = Grid_Project( vertGrid, dimInds );
			gNode2 = Sync_DomainToGlobal( vertSync, incVerts[0] );
			checkNodes = (gNode0 == gNode1) && (gNode1 == gNode2);
			pcu_check_true( checkNodes );

			gNode0 += vertGrid->sizes[0] * vertGrid->sizes[1];
			dimInds[2]++;
			gNode1 = Grid_Project( vertGrid, dimInds );
			gNode2 = Sync_DomainToGlobal( vertSync, incVerts[1] );
			checkNodes = (gNode0 == gNode1) && (gNode1 == gNode2);
			pcu_check_true( checkNodes );

			gNode0 += vertGrid->sizes[0] * vertGrid->sizes[1];
			dimInds[2]++;
			gNode1 = Grid_Project( vertGrid, dimInds );
			gNode2 = Sync_DomainToGlobal( vertSync, incVerts[2] );
			checkNodes = (gNode0 == gNode1) && (gNode1 == gNode2);
			pcu_check_true( checkNodes );
		}
	}

	FreeObject( edgeGrid_0 );
	FreeObject( edgeGrid_1 );
	FreeObject( edgeGrid_2 );

	NewClass_Delete( inc );
}

void C2GeneratorSuite_TestFaceVertexInc( C2GeneratorSuiteData* data ) {
	Grid*		elGrid		= *(Grid**)Mesh_GetExtension( data->mesh, Grid**, "elementGrid" );
	Grid*		vertGrid	= *(Grid**)Mesh_GetExtension( data->mesh, Grid**, "vertexGrid" );
	Grid*		faceGrid_0	= Grid_New();
	Grid*		faceGrid_1	= Grid_New();
	Grid*		faceGrid_2	= Grid_New();
	Sync*		vertSync	= (Sync*)IGraph_GetDomain( (IGraph*)data->mesh->topo, MT_VERTEX );
	Sync*		faceSync	= (Sync*)IGraph_GetDomain( (IGraph*)data->mesh->topo, MT_FACE );
	IArray*		inc		= IArray_New();
	unsigned*	incVerts;
	unsigned	nIncVerts;
	unsigned	dim		= ((IGraph*)data->mesh->topo)->nDims;
	unsigned	sizes[3];
	unsigned	face_i;
	unsigned	gFace;
	unsigned	dimInds[3];
	unsigned	gNode0, gNode1, gNode2;
	int		checkNodes;

	sizes[0] = elGrid->sizes[0];
	sizes[1] = elGrid->sizes[1];
	sizes[2] = elGrid->sizes[2] + 1;
	Grid_SetNumDims( faceGrid_0, dim );
	Grid_SetSizes( faceGrid_0, sizes );
	sizes[0] = elGrid->sizes[0];
	sizes[1] = elGrid->sizes[1] + 1;
	sizes[2] = elGrid->sizes[2];
	Grid_SetNumDims( faceGrid_1, dim );
	Grid_SetSizes( faceGrid_1, sizes );
	sizes[0] = elGrid->sizes[0] + 1;
	sizes[1] = elGrid->sizes[1];
	sizes[2] = elGrid->sizes[2];
	Grid_SetNumDims( faceGrid_2, dim );
	Grid_SetSizes( faceGrid_2, sizes );

	for( face_i = 0; face_i < ((IGraph*)data->mesh->topo)->remotes[MT_FACE]->nDomains; face_i++ ) {
		gFace = Sync_DomainToGlobal( faceSync, face_i );

		MeshTopology_GetIncidence( (IGraph*)data->mesh->topo, MT_FACE, face_i, MT_VERTEX, inc );
		nIncVerts = IArray_GetSize( inc );
		incVerts = IArray_GetPtr( inc );

		pcu_check_true( nIncVerts == 9 );

		if( gFace < faceGrid_0->nPoints ) {
			Grid_Lift( faceGrid_0, gFace, dimInds );
			dimInds[0] *= 2; dimInds[1] *= 2; dimInds[2] *= 2;

			gNode0 = Grid_Project( vertGrid, dimInds );
			gNode1 = Grid_Project( vertGrid, dimInds );
			gNode2 = Sync_DomainToGlobal( vertSync, incVerts[0] );
			checkNodes = (gNode0 == gNode1) && (gNode1 == gNode2);
			pcu_check_true( checkNodes );

			gNode0++;
			dimInds[0]++;
			gNode1 = Grid_Project( vertGrid, dimInds );
			gNode2 = Sync_DomainToGlobal( vertSync, incVerts[1] );
			checkNodes = (gNode0 == gNode1) && (gNode1 == gNode2);
			pcu_check_true( checkNodes );

			gNode0++;
			dimInds[0]++;
			gNode1 = Grid_Project( vertGrid, dimInds );
			gNode2 = Sync_DomainToGlobal( vertSync, incVerts[2] );
			checkNodes = (gNode0 == gNode1) && (gNode1 == gNode2);
			pcu_check_true( checkNodes );

			gNode0 -= 2; gNode0 += vertGrid->sizes[0];
			dimInds[0] -= 2; dimInds[1]++;
			gNode1 = Grid_Project( vertGrid, dimInds );
			gNode2 = Sync_DomainToGlobal( vertSync, incVerts[3] );
			checkNodes = (gNode0 == gNode1) && (gNode1 == gNode2);
			pcu_check_true( checkNodes );

			gNode0++;
			dimInds[0]++;
			gNode1 = Grid_Project( vertGrid, dimInds );
			gNode2 = Sync_DomainToGlobal( vertSync, incVerts[4] );
			checkNodes = (gNode0 == gNode1) && (gNode1 == gNode2);
			pcu_check_true( checkNodes );

			gNode0++;
			dimInds[0]++;
			gNode1 = Grid_Project( vertGrid, dimInds );
			gNode2 = Sync_DomainToGlobal( vertSync, incVerts[5] );
			checkNodes = (gNode0 == gNode1) && (gNode1 == gNode2);
			pcu_check_true( checkNodes );

			gNode0 -= 2; gNode0 += vertGrid->sizes[0];
			dimInds[0] -= 2; dimInds[1]++;
			gNode1 = Grid_Project( vertGrid, dimInds );
			gNode2 = Sync_DomainToGlobal( vertSync, incVerts[6] );
			checkNodes = (gNode0 == gNode1) && (gNode1 == gNode2);
			pcu_check_true( checkNodes );

			gNode0++;
			dimInds[0]++;
			gNode1 = Grid_Project( vertGrid, dimInds );
			gNode2 = Sync_DomainToGlobal( vertSync, incVerts[7] );
			checkNodes = (gNode0 == gNode1) && (gNode1 == gNode2);
			pcu_check_true( checkNodes );
			
			gNode0++;
			dimInds[0]++;
			gNode1 = Grid_Project( vertGrid, dimInds );
			gNode2 = Sync_DomainToGlobal( vertSync, incVerts[8] );
			checkNodes = (gNode0 == gNode1) && (gNode1 == gNode2);
			pcu_check_true( checkNodes );
		}
		else if( gFace < faceGrid_0->nPoints + faceGrid_1->nPoints ) {
			Grid_Lift( faceGrid_1, gFace - faceGrid_0->nPoints, dimInds );
			dimInds[0] *= 2; dimInds[1] *= 2; dimInds[2] *= 2;

			gNode0 = Grid_Project( vertGrid, dimInds );
			gNode1 = Grid_Project( vertGrid, dimInds );
			gNode2 = Sync_DomainToGlobal( vertSync, incVerts[0] );
			checkNodes = (gNode0 == gNode1) && (gNode1 == gNode2);
			pcu_check_true( checkNodes );

			gNode0++;
			dimInds[0]++;
			gNode1 = Grid_Project( vertGrid, dimInds );
			gNode2 = Sync_DomainToGlobal( vertSync, incVerts[1] );
			checkNodes = (gNode0 == gNode1) && (gNode1 == gNode2);
			pcu_check_true( checkNodes );

			gNode0++;
			dimInds[0]++;
			gNode1 = Grid_Project( vertGrid, dimInds );
			gNode2 = Sync_DomainToGlobal( vertSync, incVerts[2] );
			checkNodes = (gNode0 == gNode1) && (gNode1 == gNode2);
			pcu_check_true( checkNodes );

			gNode0 -= 2; gNode0 += vertGrid->sizes[0] * vertGrid->sizes[1];
			dimInds[0] -= 2; dimInds[2]++;
			gNode1 = Grid_Project( vertGrid, dimInds );
			gNode2 = Sync_DomainToGlobal( vertSync, incVerts[3] );
			checkNodes = (gNode0 == gNode1) && (gNode1 == gNode2);
			pcu_check_true( checkNodes );

			gNode0++;
			dimInds[0]++;
			gNode1 = Grid_Project( vertGrid, dimInds );
			gNode2 = Sync_DomainToGlobal( vertSync, incVerts[4] );
			checkNodes = (gNode0 == gNode1) && (gNode1 == gNode2);
			pcu_check_true( checkNodes );

			gNode0++;
			dimInds[0]++;
			gNode1 = Grid_Project( vertGrid, dimInds );
			gNode2 = Sync_DomainToGlobal( vertSync, incVerts[5] );
			checkNodes = (gNode0 == gNode1) && (gNode1 == gNode2);
			pcu_check_true( checkNodes );

			gNode0 -= 2; gNode0 += vertGrid->sizes[0] * vertGrid->sizes[1];
			dimInds[0] -= 2; dimInds[2]++;
			gNode1 = Grid_Project( vertGrid, dimInds );
			gNode2 = Sync_DomainToGlobal( vertSync, incVerts[6] );
			checkNodes = (gNode0 == gNode1) && (gNode1 == gNode2);
			pcu_check_true( checkNodes );

			gNode0++;
			dimInds[0]++;
			gNode1 = Grid_Project( vertGrid, dimInds );
			gNode2 = Sync_DomainToGlobal( vertSync, incVerts[7] );
			checkNodes = (gNode0 == gNode1) && (gNode1 == gNode2);
			pcu_check_true( checkNodes );

			gNode0++;
			dimInds[0]++;
			gNode1 = Grid_Project( vertGrid, dimInds );
			gNode2 = Sync_DomainToGlobal( vertSync, incVerts[8] );
			checkNodes = (gNode0 == gNode1) && (gNode1 == gNode2);
			pcu_check_true( checkNodes );
		}
		else if( gFace < faceGrid_0->nPoints + faceGrid_1->nPoints + faceGrid_2->nPoints ) {
			Grid_Lift( faceGrid_2, gFace - faceGrid_0->nPoints - faceGrid_1->nPoints, dimInds );
			dimInds[0] *= 2; dimInds[1] *= 2; dimInds[2] *= 2;

			gNode0 = Grid_Project( vertGrid, dimInds );
			gNode1 = Grid_Project( vertGrid, dimInds );
			gNode2 = Sync_DomainToGlobal( vertSync, incVerts[0] );
			checkNodes = (gNode0 == gNode1) && (gNode1 == gNode2);
			pcu_check_true( checkNodes );

			gNode0 += vertGrid->sizes[0];
			dimInds[1]++;
			gNode1 = Grid_Project( vertGrid, dimInds );
			gNode2 = Sync_DomainToGlobal( vertSync, incVerts[1] );
			checkNodes = (gNode0 == gNode1) && (gNode1 == gNode2);
			pcu_check_true( checkNodes );

			gNode0 += vertGrid->sizes[0];
			dimInds[1]++;
			gNode1 = Grid_Project( vertGrid, dimInds );
			gNode2 = Sync_DomainToGlobal( vertSync, incVerts[2] );
			checkNodes = (gNode0 == gNode1) && (gNode1 == gNode2);
			pcu_check_true( checkNodes );

			gNode0 -= 2 * vertGrid->sizes[0]; gNode0 += vertGrid->sizes[0] * vertGrid->sizes[1];
			dimInds[1] -= 2; dimInds[2]++;
			gNode1 = Grid_Project( vertGrid, dimInds );
			gNode2 = Sync_DomainToGlobal( vertSync, incVerts[3] );
			checkNodes = (gNode0 == gNode1) && (gNode1 == gNode2);
			pcu_check_true( checkNodes );

			gNode0 += vertGrid->sizes[0];
			dimInds[1]++;
			gNode1 = Grid_Project( vertGrid, dimInds );
			gNode2 = Sync_DomainToGlobal( vertSync, incVerts[4] );
			checkNodes = (gNode0 == gNode1) && (gNode1 == gNode2);
			pcu_check_true( checkNodes );

			gNode0 += vertGrid->sizes[0];
			dimInds[1]++;
			gNode1 = Grid_Project( vertGrid, dimInds );
			gNode2 = Sync_DomainToGlobal( vertSync, incVerts[5] );
			checkNodes = (gNode0 == gNode1) && (gNode1 == gNode2);
			pcu_check_true( checkNodes );

			gNode0 -= 2 * vertGrid->sizes[0]; gNode0 += vertGrid->sizes[0] * vertGrid->sizes[1];
			dimInds[1] -= 2; dimInds[2]++;
			gNode1 = Grid_Project( vertGrid, dimInds );
			gNode2 = Sync_DomainToGlobal( vertSync, incVerts[6] );
			checkNodes = (gNode0 == gNode1) && (gNode1 == gNode2);
			pcu_check_true( checkNodes );

			gNode0 += vertGrid->sizes[0];
			dimInds[1]++;
			gNode1 = Grid_Project( vertGrid, dimInds );
			gNode2 = Sync_DomainToGlobal( vertSync, incVerts[7] );
			checkNodes = (gNode0 == gNode1) && (gNode1 == gNode2);
			pcu_check_true( checkNodes );

			gNode0 += vertGrid->sizes[0];
			dimInds[1]++;
			gNode1 = Grid_Project( vertGrid, dimInds );
			gNode2 = Sync_DomainToGlobal( vertSync, incVerts[8] );
			checkNodes = (gNode0 == gNode1) && (gNode1 == gNode2);
			pcu_check_true( checkNodes );
		}
	}

	FreeObject( faceGrid_0 );
	FreeObject( faceGrid_1 );
	FreeObject( faceGrid_2 );

	NewClass_Delete( inc );
}
void C2GeneratorSuite( pcu_suite_t* suite ) {
   pcu_suite_setData( suite, C2GeneratorSuiteData );
   pcu_suite_setFixtures( suite, C2GeneratorSuite_Setup, C2GeneratorSuite_Teardown );
   pcu_suite_addTest( suite, C2GeneratorSuite_TestElementVertexInc );
   pcu_suite_addTest( suite, C2GeneratorSuite_TestEdgeVertexInc );
   pcu_suite_addTest( suite, C2GeneratorSuite_TestFaceVertexInc );
}


