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
** $Id: testSwarm.c 4118 2007-05-22 02:01:25Z RaquibulHassan $
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
	__GlobalParticle
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

void TestParticleSearchFunc( Swarm* swarm, Coord coord, Stream* stream ) {
	double              distance;
	Particle_Index      closestParticle_I;
	Particle_Index      lParticle_I;
	GlobalParticle*     particle;

	Journal_Printf( stream, "Testing coord %g %g %g\n", coord[ I_AXIS ], coord[ J_AXIS ], coord[ K_AXIS ] );

	closestParticle_I = Swarm_FindClosestParticle( swarm, 3, coord, &distance );

	if ( closestParticle_I < swarm->particleLocalCount ) {
		Journal_Printf( stream, "Closest Particle is %u with distance %g\n", closestParticle_I, distance );

		for ( lParticle_I = 0 ; lParticle_I < swarm->particleLocalCount ; lParticle_I++ ) {
			/* don't bother with testing same particle */
			if (lParticle_I == closestParticle_I)
				continue;

			particle = (GlobalParticle*)Swarm_ParticleAt( swarm , lParticle_I );

			/* Test if particle is closer */
			if (distance > StGermain_DistanceBetweenPoints( coord, particle->coord, 3 ) )
				Journal_Printf( stream, "Particle %u is closer to coord and has distance %g.\n" , 
						lParticle_I, StGermain_DistanceBetweenPoints( coord, particle->coord, 3 ) );
		}
	}
	else
		Journal_Printf( stream, "Coord not found on this processor.\n" );
}

int main( int argc, char* argv[] ) {
	MPI_Comm			CommWorld;
	int				rank;
	int				numProcessors;
	int				procToWatch;
	Dictionary*			dictionary;
	unsigned	nDims = 3;
	unsigned	meshSize[3] = {2, 3, 2};
	double		minCrds[3] = {0.0, 0.0, 0.0};
	double		maxCrds[3] = {300.0, 12.0, 300.0};
	ExtensionManager_Register*		extensionMgr_Register;
	Mesh*				mesh;
	ElementCellLayout*		elementCellLayout;
	RandomParticleLayout*		randomParticleLayout;
	Swarm*				swarm;
	Stream*				stream;
	Coord                coord = { 120, 2, 210 };

	
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
	
	/* Init mesh */
	extensionMgr_Register = ExtensionManager_Register_New();
	mesh = buildMesh( nDims, meshSize, minCrds, maxCrds, extensionMgr_Register );
	
	/* Build the mesh */
	Build( mesh, 0, False );
	Initialise( mesh, 0, False );
	
	/* Configure the element-cell-layout */
	elementCellLayout = ElementCellLayout_New( "elementCellLayout", mesh );
	
	/* Configure the random-particle-layout */
	randomParticleLayout = RandomParticleLayout_New( "randomParticleCellLayout", 4, 13 );
	
	/* Configure the swarm */
	swarm = Swarm_New( "testSwarm", elementCellLayout, randomParticleLayout, 3, sizeof(Particle),
		extensionMgr_Register, NULL, CommWorld );
	
	/* Build the swarm */
	Build( swarm, 0, False );
	Initialise( swarm, 0, False );
	
	if( rank == procToWatch ) {
		Stg_Class_Print( swarm, stream );

		Journal_Printf( stream, "\n" );
		Swarm_PrintParticleCoords( swarm, stream );
		Journal_Printf( stream, "\n" );
		Swarm_PrintParticleCoords_ByCell( swarm, stream );
	
		/* Test Particle Searching Function */
		if (numProcessors == 1) {
			coord[ I_AXIS ] = 0.6 * 300; coord[ J_AXIS ] = 0.2*12; coord[ K_AXIS ] = 0.820*300;
			TestParticleSearchFunc( swarm, coord, stream );

			/* Test New Coord */
			coord[ I_AXIS ] = 0.2*300; coord[ J_AXIS ] = 0.9*12; coord[ K_AXIS ] = 0.120*300;
			TestParticleSearchFunc( swarm, coord, stream );
		}
	}
	
	/* Destroy stuff */
	Stg_Class_Delete( swarm );
	Stg_Class_Delete( randomParticleLayout );
	Stg_Class_Delete( elementCellLayout );
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
