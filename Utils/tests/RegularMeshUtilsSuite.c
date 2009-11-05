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
**   Tests the RegularMeshUtilsSuite
**
** $Id: testTemplate.c 3462 2006-02-19 06:53:24Z WalterLandry $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>

#include "pcu/pcu.h"
#include <StGermain/StGermain.h> 
#include "StgDomain/Geometry/Geometry.h"
#include "StgDomain/Shape/Shape.h"
#include "StgDomain/Mesh/Mesh.h" 
#include "StgDomain/Utils/Utils.h"
#include "StgDomain/Swarm/Swarm.h"

#include "RegularMeshUtilsSuite.h"

typedef struct {
	MPI_Comm	comm;
	unsigned	rank;
	unsigned	nProcs;
} RegularMeshUtilsSuiteData;

Mesh* RegularMeshUtilsSuite_buildMesh( unsigned nDims, unsigned* size, double* minCrds, double* maxCrds, ExtensionManager_Register* emReg ) {
	CartesianGenerator*	gen;
	Mesh*						mesh;
	unsigned					maxDecomp[3] = {0, 1, 1};

	gen = CartesianGenerator_New( "" );
	CartesianGenerator_SetDimSize( gen, nDims );
	CartesianGenerator_SetTopologyParams( gen, size, 0, NULL, maxDecomp );
	CartesianGenerator_SetGeometryParams( gen, minCrds, maxCrds );

	mesh = Mesh_New( "" );
	Mesh_SetExtensionManagerRegister( mesh, emReg );
	Mesh_SetGenerator( mesh, gen );

	Stg_Component_Build( mesh, NULL, False );
	Stg_Component_Initialise( mesh, NULL, False );

	FreeObject( mesh->generator );

	return mesh;
}


void RegularMeshUtilsSuite_Setup( RegularMeshUtilsSuiteData* data ) {
	/* MPI Initializations */
	data->comm = MPI_COMM_WORLD;
	MPI_Comm_rank( data->comm, &data->rank );
	MPI_Comm_size( data->comm, &data->nProcs );
}

void RegularMeshUtilsSuite_Teardown( RegularMeshUtilsSuiteData* data ) {
}

void RegularMeshUtilsSuite_TestMeshUtils( RegularMeshUtilsSuiteData* data ) {
	ExtensionManager_Register*	extensionMgr_Register;
	unsigned							nDims = 3;
	unsigned							meshSize[3] = {6, 6, 6};
	double							minCrds[3] = {0.0, 0.0, 0.0};
	double							maxCrds[3] = {1.0, 1.0, 1.0};
	Mesh*								mesh;
	char								expected_file[PCU_PATH_MAX];
	int								procToWatch;
	Stream*							stream; 	

	Journal_Enable_NamedStream( Info_Type, CartesianGenerator_Type, False );
	stream = Journal_Register( Info_Type, "RegularMeshUtilsStream" );

	procToWatch = data->nProcs >=2 ? 1 : 0;

	extensionMgr_Register = ExtensionManager_Register_New();
	mesh = RegularMeshUtilsSuite_buildMesh( nDims, meshSize, minCrds, maxCrds, extensionMgr_Register );

	if( data->rank == procToWatch ) {
		unsigned		currElementNodesCount=0;
		unsigned*	currElementNodes = NULL;
		unsigned		element_dI = 0;
		unsigned		refNode_eI = 0;
		unsigned		node_Diagonal = 0;
		unsigned		node_Diagonal_gI = 0;
		IArray*		inc;

		inc = IArray_New();
		
		Stream_RedirectFile( stream, "regularMeshUtils.dat" );

		for (element_dI=0; element_dI < Mesh_GetDomainSize( mesh, nDims ); element_dI++) {
			Mesh_GetIncidence( mesh, nDims, element_dI, MT_VERTEX, inc );
			currElementNodesCount = IArray_GetSize( inc );
			currElementNodes = IArray_GetPtr( inc );

			for (refNode_eI = 0; refNode_eI < currElementNodesCount; refNode_eI++ ) {
				node_Diagonal = RegularMeshUtils_GetDiagOppositeAcrossElementNodeIndex(mesh, element_dI, currElementNodes[refNode_eI]) ;
				node_Diagonal_gI = Mesh_DomainToGlobal( mesh, MT_VERTEX, node_Diagonal );
				/*print message stating: Element #, curr node #, diag opp node #*/
				Journal_Printf( stream, "Element #: %d, Current Node #: %d, Diagonal Node #: %d, (%d) \n", element_dI, currElementNodes[refNode_eI], node_Diagonal, node_Diagonal_gI );
			}
		}
		pcu_filename_expected( "testRegularMeshUtilsOutput.expected", expected_file );
		pcu_check_fileEq( "regularMeshUtils.dat", expected_file );
		NewClass_Delete( inc );
		remove( "regularMeshUtils.dat" );
	}
}
	
void RegularMeshUtilsSuite( pcu_suite_t* suite ) {
	pcu_suite_setData( suite, RegularMeshUtilsSuiteData );
	pcu_suite_setFixtures( suite, RegularMeshUtilsSuite_Setup, RegularMeshUtilsSuite_Teardown );
	pcu_suite_addTest( suite, RegularMeshUtilsSuite_TestMeshUtils );
}
