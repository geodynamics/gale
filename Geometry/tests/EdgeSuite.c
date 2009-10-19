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
**   Tests the EdgeSuite
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

#include "EdgeSuite.h"

typedef struct {
	MPI_Comm comm;
	unsigned rank;
	unsigned nProcs;
} EdgeSuiteData;

void EdgeSuite_Setup( EdgeSuiteData* data ) {
	/* MPI Initializations */
	data->comm = MPI_COMM_WORLD;
	MPI_Comm_rank( data->comm, &data->rank );
	MPI_Comm_size( data->comm, &data->nProcs );
}

void EdgeSuite_Teardown( EdgeSuiteData* data ) {
}

void EdgeSuite_TestEdge( EdgeSuiteData* data ) {
	Stream*	stream;
	char		expected_file[PCU_PATH_MAX];
	unsigned	procToWatch;

	procToWatch = data->nProcs >=2 ? 1 : 0;
	
   if( data->rank == procToWatch ) {
		Triangle			tri[100];
      Triangle_Index	triCount = 100;
		Edge_Index		edgeCount;
		Edge_List		edge;
		EdgeFaces_List	edgeFace;
		Triangle_Index	tri_I;
		Edge_Index		edge_I;

		stream = Journal_Register(Info_Type, "EdgeStream");
		Stream_RedirectFile( stream, "testEdge.dat" );
 
		tri[0][0] = 0;
		tri[0][1] = 1;
		tri[0][2] = 2;
		for( tri_I = 1; tri_I < triCount; tri_I++ ) {
			tri[tri_I][0] = tri[tri_I - 1][0] + 1;
			tri[tri_I][1] = tri[tri_I - 1][2] + 1;
			tri[tri_I][2] = tri[tri_I - 1][1] + 1;
		}

		edgeCount = Edge_BuildList_FromTriangles( tri, triCount, &edge, &edgeFace );

		for( edge_I = 0; edge_I < edgeCount; edge_I++ ) {
			Journal_Printf( stream, "Edge: %u, points: { %u, %u }, ", edge_I, edge[edge_I][0], edge[edge_I][1] );
			Journal_Printf( stream, "faces: { %u, %u }\n", edgeFace[edge_I][0], edgeFace[edge_I][1] );
		}

		pcu_filename_expected( "testEdge.expected", expected_file );
		pcu_check_fileEq( "testEdge.dat", expected_file );
		remove( "testEdge.dat" );

		if( edge ) Memory_Free( edge );
		if( edgeFace ) Memory_Free( edgeFace );
	}
}

void EdgeSuite( pcu_suite_t* suite ) {
	pcu_suite_setData( suite, EdgeSuiteData );
   pcu_suite_setFixtures( suite, EdgeSuite_Setup, EdgeSuite_Teardown );
   pcu_suite_addTest( suite, EdgeSuite_TestEdge );
}
