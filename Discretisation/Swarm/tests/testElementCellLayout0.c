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
** Role:
**	Test that the ElementCellLayout has the same layout and geometry as the mesh's element layout.
**
** Assumptions:
**	None as yet.
**
** Comments:
**	None as yet.
**
** $Id: testElementCellLayout0.c 4081 2007-04-27 06:20:07Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include "Base/Base.h"

#include "Discretisation/Geometry/Geometry.h"
#include "Discretisation/Shape/Shape.h"
#include "Discretisation/Mesh/Mesh.h"
#include "Discretisation/Utils/Utils.h"
#include "Discretisation/Swarm/Swarm.h"

#include <stdio.h>
#include <stdlib.h>

#define ConvertNode(node) node


struct _Node {
	double temp;
};

struct _Element {
	double temp;
};

Mesh* buildMesh( unsigned nDims, unsigned* size, 
		     double* minCrds, double* maxCrds, 
		     ExtensionManager_Register* emReg )
{
	CartesianGenerator*	gen;
	Mesh*			mesh;
	unsigned		maxDecomp[3] = {1, 0, 1};

	gen = CartesianGenerator_New( "" );
	CartesianGenerator_SetDimSize( gen, nDims );
	CartesianGenerator_SetTopologyParams( gen, size, 0, NULL, maxDecomp );
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
	MPI_Comm			CommWorld;
	int				rank;
	int				numProcessors;
	int				procToWatch;
	unsigned	nDims = 3;
	unsigned	meshSize[3] = {2, 3, 2};
	double		minCrds[3] = {0.0, 0.0, 0.0};
	double		maxCrds[3] = {300.0, 12.0, 300.0};
	ExtensionManager_Register*		extensionMgr_Register;
	Mesh*				mesh;
	ElementCellLayout*		elementCellLayout;

	/* Initialise MPI, get world info */
	MPI_Init( &argc, &argv );
	MPI_Comm_dup( MPI_COMM_WORLD, &CommWorld );
	MPI_Comm_size( CommWorld, &numProcessors );
	MPI_Comm_rank( CommWorld, &rank );
	
	Base_Init( &argc, &argv );
	
	DiscretisationGeometry_Init( &argc, &argv );
	DiscretisationShape_Init( &argc, &argv );
	DiscretisationMesh_Init( &argc, &argv );
	DiscretisationUtils_Init( &argc, &argv );
	DiscretisationSwarm_Init( &argc, &argv );
	MPI_Barrier( CommWorld ); /* Ensures copyright info always come first in output */

	if( argc >= 2 ) {
		procToWatch = atoi( argv[1] );
	}
	else {
		procToWatch = 0;
	}
	if( rank == procToWatch ) printf( "Watching rank: %i\n", rank );
	
	/* Init mesh */
	extensionMgr_Register = ExtensionManager_Register_New();
	mesh = buildMesh( nDims, meshSize, minCrds, maxCrds, extensionMgr_Register );
	
	/* Configure the element-cell-layout */
	elementCellLayout = ElementCellLayout_New( "elementCellLayout", mesh );
	
	/* Build the mesh */
	Build( mesh, 0, False );
	Initialise( mesh, 0, False );
	
	if( rank == procToWatch ) {
		Cell_Index cell;
		Element_DomainIndex	element;
		GlobalParticle          testParticle;
		
		for( element = 0; element < Mesh_GetLocalSize( mesh, nDims ); element++ ) {
			Cell_PointIndex			point;
			Cell_PointIndex			count;
			double***			cellPoints;
			Bool				result;

			cell = CellLayout_MapElementIdToCellId( elementCellLayout, element );

			if ( cell != element ) { 
				printf( "Wrong result: CellLayout_MapElementIdToCellId returned %d, when %d "
					"expected.\n", cell, element );
				exit(0);
			}

			count = elementCellLayout->_pointCount( elementCellLayout, cell );
			printf( "cellPointTbl  [%2u][0-%u]:\n", cell, count );
			cellPoints = Memory_Alloc_Array( double**, count, "cellPoints" );
			elementCellLayout->_initialisePoints( elementCellLayout, cell, count, cellPoints );
			for( point = 0; point < count; point++ ) {
				printf( "\t{%.3g %.3g %.3g}\n", (*cellPoints[ConvertNode(point)])[0], (*cellPoints[ConvertNode(point)])[1], 
					(*cellPoints[ConvertNode(point)])[2] );
			}
			printf( "\n" );

			testParticle.coord[0] = ( (*cellPoints[0])[0] + (*cellPoints[1])[0] ) / 2;
			testParticle.coord[1] = ( (*cellPoints[0])[1] + (*cellPoints[2])[1] ) / 2;
			testParticle.coord[2] = ( (*cellPoints[0])[2] + (*cellPoints[4])[2] ) / 2;
			printf( "Testing if test particle at (%f,%f,%f) is in the cell: ",
				testParticle.coord[0], testParticle.coord[1], testParticle.coord[2] );
			result = CellLayout_IsInCell( elementCellLayout, cell, &testParticle );
			printf( "%d\n\n", result );

			testParticle.coord[0] = (*cellPoints[count-2])[0] + 1;
			testParticle.coord[1] = (*cellPoints[count-2])[1] + 1;
			testParticle.coord[2] = (*cellPoints[count-2])[2] + 1;
			printf( "Testing if test particle at (%f,%f,%f) is in the cell: ",
				testParticle.coord[0], testParticle.coord[1], testParticle.coord[2] );
			result = CellLayout_IsInCell( elementCellLayout, cell, &testParticle );
			printf( "%d\n\n", result );

			Memory_Free( cellPoints );
		}
	}
	
	/* Destroy stuff */
	Stg_Class_Delete( elementCellLayout );
	Stg_Class_Delete( mesh );
	Stg_Class_Delete( extensionMgr_Register );
	
	DiscretisationSwarm_Finalise();
	DiscretisationUtils_Finalise();
	DiscretisationMesh_Finalise();
	DiscretisationShape_Finalise();
	DiscretisationGeometry_Finalise();
	
	Base_Finalise();
	
	/* Close off MPI */
	MPI_Finalize();
	
	return 0; /* success */
}
