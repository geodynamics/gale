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
** $Id: testGaussLayout.c 4175 2007-08-16 03:39:26Z DavidLee $
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

struct _Node {
	Coord				coord;
};

struct _Element {
	Coord				coord;
};

struct _Particle {
	__IntegrationPoint
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

	Stg_Component_Build( mesh, NULL, False );
	Stg_Component_Initialise( mesh, NULL, False );

	KillObject( mesh->generator );

	return mesh;
}

int main( int argc, char* argv[] ) {
	MPI_Comm			CommWorld;
	int				rank;
	int				numProcessors;
	int				procToWatch;
	Dictionary*	dictionary;
	unsigned	nDims = 3;
	unsigned	meshSize[3] = {2, 3, 2};
	double		minCrds[3] = {0.0, 0.0, 0.0};
	double		maxCrds[3] = {300.0, 12.0, 300.0};
	ExtensionManager_Register*		extensionMgr_Register;
	Mesh*				mesh;
	GaussParticleLayout*		gaussParticleLayout;
	ElementCellLayout*		elementCellLayout;
	Swarm*				swarm;
	Stream*				stream;
	Dimension_Index     dim;
	
	Cell_PointIndex			count;
	double x,y,z;
	int p;
	LocalParticle* particle;
	Coord minCell;
	Coord maxCell;
	Particle_InCellIndex particlesPerDim[3];
	
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

	stream = Journal_Register (Info_Type, "myStream");

	if( argc >= 2 ) {
		procToWatch = atoi( argv[1] );
	}
	else {
		procToWatch = 0;
	}
	if( rank == procToWatch ) printf( "Watching rank: %i\n", rank );
	
	/* Read input */
	dictionary = Dictionary_New();
	Dictionary_Add( dictionary, "rank", Dictionary_Entry_Value_FromUnsignedInt( rank ) );
	Dictionary_Add( dictionary, "numProcessors", Dictionary_Entry_Value_FromUnsignedInt( numProcessors ) );
	Dictionary_Add( dictionary, "gaussParticlesX", Dictionary_Entry_Value_FromUnsignedInt( 2 ) );
	Dictionary_Add( dictionary, "gaussParticlesY", Dictionary_Entry_Value_FromUnsignedInt( 2 ) );
	Dictionary_Add( dictionary, "gaussParticlesZ", Dictionary_Entry_Value_FromUnsignedInt( 2 ) );
	Dictionary_Add( dictionary, "dim", Dictionary_Entry_Value_FromUnsignedInt( 3 ) );
	
	/* Init mesh */
	extensionMgr_Register = ExtensionManager_Register_New();
	mesh = buildMesh( nDims, meshSize, minCrds, maxCrds, extensionMgr_Register );
	
	/* Configure the element-cell-layout */
	elementCellLayout = ElementCellLayout_New( "elementCellLayout", mesh );
	
	/* Build the mesh */
	Stg_Component_Build( mesh, 0, False );
	Stg_Component_Initialise( mesh, 0, False );
	
	/* Configure the gauss-particle-layout */
	dim = Dictionary_GetUnsignedInt( dictionary, "dim" );
	particlesPerDim[0] = Dictionary_GetUnsignedInt( dictionary, "gaussParticlesX" );
	particlesPerDim[1] = Dictionary_GetUnsignedInt( dictionary, "gaussParticlesY" );
	particlesPerDim[2] = Dictionary_GetUnsignedInt( dictionary, "gaussParticlesZ" );
	gaussParticleLayout = GaussParticleLayout_New( "gaussParticleLayout", dim, particlesPerDim );
	
	/* Configure the swarm */
	swarm = Swarm_New( "testGaussSwarm", elementCellLayout, gaussParticleLayout, dim, sizeof(Particle),
		extensionMgr_Register, NULL, CommWorld, NULL );
	
	
	/* Build the swarm */
	Stg_Component_Build( swarm, 0, False );
	Stg_Component_Initialise( swarm, 0, False );
	
	if( rank == procToWatch ) {
		Stg_Class_Print( gaussParticleLayout, stream );
	}
	
	
	/* Print out the particles on any cell (they should all be the same), lets try 4 */
	count = swarm->cellParticleCountTbl[4];
	printf("count = %d \n",count );
	
	Swarm_GetCellMinMaxCoords( swarm, 4, minCell, maxCell );
	
	for( p = 0; p < count; p++ ) {
		particle = (LocalParticle*)Swarm_ParticleInCellAt( swarm, 4, p  );

		/* convert to global coords */
		x = 0.5 * ( maxCell[0] - minCell[0] ) * ( particle->xi[0] + 1.0 ) + minCell[0];
		y = 0.5 * ( maxCell[1] - minCell[1] ) * ( particle->xi[1] + 1.0 ) + minCell[1];
		z = 0.5 * ( maxCell[2] - minCell[2] ) * ( particle->xi[2] + 1.0 ) + minCell[2];
		
		printf( "pId=%d : coords = { %f, %f, %f }, xi = { %f, %f, %f }\n", 
				p, x, y, z, particle->xi[0], particle->xi[1], particle->xi[2] );
	}	
	
	/* Destroy stuff */
	Stg_Class_Delete( gaussParticleLayout );
	Stg_Class_Delete( elementCellLayout );
	Stg_Class_Delete( swarm );
	Stg_Class_Delete( mesh );
	Stg_Class_Delete( extensionMgr_Register );
	Stg_Class_Delete( dictionary );
	
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
