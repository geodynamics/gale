/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street, Melbourne, 3053, Australia.
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
** $Id: testSemiRegDeform.c 3124 2005-07-25 04:52:06Z RobertTurnbull $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <mpi.h>

#include <StGermain/StGermain.h>
#include <StgDomain/Geometry/Geometry.h>
#include <StgDomain/Shape/Shape.h>
#include <StgDomain/Mesh/Mesh.h>

#include "StgDomain/Utils/types.h"
#include "StgDomain/Utils/SemiRegDeform.h"


struct _Node {
	double tmp;
};

struct _Element {
	double tmp;
};


Mesh* buildMesh( unsigned nDims, unsigned* size, 
		     double* minCrds, double* maxCrds, 
		     ExtensionManager_Register* emReg )
{
	CartesianGenerator*	gen;
	Mesh*			mesh;

	gen = CartesianGenerator_New( "" );
	CartesianGenerator_SetDimSize( gen, nDims );
	CartesianGenerator_SetTopologyParams( gen, size, 0, NULL, NULL );
	CartesianGenerator_SetGeometryParams( gen, minCrds, maxCrds );

	mesh = Mesh_New( "" );
	Mesh_SetExtensionManagerRegister( mesh, emReg );
	Mesh_SetGenerator( mesh, gen );

	Build( mesh, NULL, False );
	Initialise( mesh, NULL, False );

	KillObject( mesh->generator );

	return mesh;
}


int main( int argc, char* argv[] ) {
	MPI_Comm			commWorld;
	int				rank;
	int				nProcs;
	int				procToWatch;
	ExtensionManager_Register*	extensionMgr_Register;
	Stream*				stream;
	SemiRegDeform*			srd;

	unsigned	nDims = 3;
	unsigned	meshSize[3] = {4, 4, 4};
	double		minCrds[3] = {0.0, 0.0, 0.0};
	double		maxCrds[3] = {1.0, 1.0, 1.0};
	Mesh*		mesh;
	Grid*		vertGrid;

	/*
	** Initialise MPI, StGermain Base, get world info.
	*/

	MPI_Init( &argc, &argv );
	MPI_Comm_dup( MPI_COMM_WORLD, &commWorld );
	MPI_Comm_size( commWorld, &nProcs );
	MPI_Comm_rank( commWorld, &rank );

	StGermain_Init( &argc, &argv );
	StgDomainGeometry_Init( &argc, &argv );
	StgDomainShape_Init( &argc, &argv );
	StgDomainMesh_Init( &argc, &argv );
	MPI_Barrier( commWorld ); /* Ensures copyright info always come first in output */

	stream = Journal_Register( Info_Type, "myStream" );
	procToWatch = argc >= 2 ? atoi(argv[1]) : 0;

	/*
	** Create the mesh.
	*/

	mesh = buildMesh( nDims, meshSize, minCrds, maxCrds, extensionMgr_Register );
	vertGrid = *(Grid**)ExtensionManager_Get( mesh->info, mesh, 
						  ExtensionManager_GetHandle( mesh->info, "vertexGrid" ) );

	/*
	** Create the deformation.
	*/

	{
		IJK		ijk;
		unsigned	lower, upper;

		srd = SemiRegDeform_New( "SemiRegDeform" );
		SemiRegDeform_SetMesh( srd, mesh );

		/* Set up strips to remesh in the y direction. */
		for( ijk[2] = 0; ijk[2] < meshSize[2]; ijk[2]++ ) {
			for( ijk[0] = 0; ijk[0] < meshSize[0]; ijk[0]++ ) {
				ijk[1] = 0;
				lower = Grid_Project( vertGrid, ijk );

				ijk[1] = meshSize[1] - 1;
				upper = Grid_Project( vertGrid, ijk );

				SemiRegDeform_AddStrip( srd, lower, upper );
			}
		}

		/* Build and initialise. */
		Build( srd, 0, False );
		Initialise( srd, 0, False );

		/* Execute the deformation. */
		SemiRegDeform_Deform( srd );

		/* Check the deformation. */
		if (rank == procToWatch) {
		}

		/* Kill it. */
		Stg_Class_Delete( srd );
	}

	{
		IJK		ijk;
		unsigned	lower, upper;

		srd = SemiRegDeform_New( "SemiRegDeform" );
		SemiRegDeform_SetMesh( srd, mesh );

		/* Set up strips to remesh in the y direction. */
		for( ijk[2] = 0; ijk[2] < meshSize[2]; ijk[2]++ ) {
			for( ijk[0] = 0; ijk[0] < meshSize[0]; ijk[0]++ ) {
				ijk[1] = 0;
				lower = Grid_Project( vertGrid, ijk );
				ijk[1] = meshSize[1] - 2;
				upper = Grid_Project( vertGrid, ijk );
				SemiRegDeform_AddStrip( srd, lower, upper );

				lower = upper;
				ijk[1] = meshSize[1] - 1;
				upper = Grid_Project( vertGrid, ijk );
				SemiRegDeform_AddStrip( srd, lower, upper );
			}
		}

		/* Build and initialise. */
		Build( srd, 0, False );
		Initialise( srd, 0, False );

		/* Check the deformation. */
		if (rank == procToWatch) {
		}

		/* Kill it. */
		Stg_Class_Delete( srd );
	}


	/*
	** Cleanup.
	*/
	
	Stg_Class_Delete( mesh );

	StgDomainMesh_Finalise();
	StgDomainShape_Finalise();
	StgDomainGeometry_Finalise();
	StGermain_Finalise();

	/* Close off MPI */
	MPI_Finalize();

	return EXIT_SUCCESS;
}
