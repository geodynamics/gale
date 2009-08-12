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
**	Tests that particles can be saved to file, then re-loaded onto a new context with exactly
**	the same positions and values.
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

#include "SwarmDumpAndLoadSuite.h"

struct _Particle {
	__GlobalParticle
	Coord				xi;
	unsigned int	testValue;
};

typedef struct {
	MPI_Comm							comm;
	unsigned int					rank;
	unsigned int					nProcs;
} SwarmDumpAndLoadSuiteData;

Mesh* SwarmDumpAndLoadSuite_BuildMesh( unsigned nDims, unsigned* size, double* minCrds, double* maxCrds, ExtensionManager_Register* emReg ) {
	CartesianGenerator*	gen;
	Mesh*						mesh;

	gen = CartesianGenerator_New( "" );
	CartesianGenerator_SetDimSize( gen, nDims );
	MeshGenerator_SetDimState( gen, 2, False );
	MeshGenerator_SetDimState( gen, 1, True );
	MeshGenerator_ClearIncidenceStates( gen );
	MeshGenerator_SetIncidenceState( gen, 3, 0, True );
	/*	MeshGenerator_SetIncidenceState( gen, 1, 0, True );*/
	MeshGenerator_SetIncidenceState( gen, 0, 3, True );	
	/*	MeshGenerator_SetIncidenceState( gen, 0, 1, True );*/
	MeshGenerator_SetIncidenceState( gen, 0, 0, True );
	CartesianGenerator_SetTopologyParams( gen, size, 0, NULL, NULL );
	CartesianGenerator_SetGeometryParams( gen, minCrds, maxCrds );

	mesh = Mesh_New( "" );
	Mesh_SetExtensionManagerRegister( mesh, emReg );
	Mesh_SetGenerator( mesh, gen );

	Stg_Component_Build( mesh, NULL, False );
	Stg_Component_Initialise( mesh, NULL, False );

	KillObject( mesh->generator );

	return mesh;
}


void SwarmDumpAndLoadSuite_UpdateParticlePositionsTowardsAttractor( Swarm* swarm, Coord attractorPoint, Processor_Index rank, Processor_Index procToWatch ) {
	Cell_LocalIndex		lCell_I;
	Particle_InCellIndex	cParticle_I;
	Particle*				currParticle;
	Index						dim_I;

	for ( lCell_I=0; lCell_I < swarm->cellLocalCount; lCell_I++ ) {
		if( rank == procToWatch ) {
			//printf("\tUpdating Particles positions in local cell %d:\n", lCell_I );
		}	
		for ( cParticle_I=0; cParticle_I < swarm->cellParticleCountTbl[lCell_I]; cParticle_I++ ) {
			Coord movementVector = {0,0,0};
			Coord newParticleCoord = {0,0,0};
			Coord* oldCoord;

			currParticle = (Particle*)Swarm_ParticleInCellAt( swarm, lCell_I, cParticle_I );
			oldCoord = &currParticle->coord;
			if( rank == procToWatch ) {
				//printf("\t\tUpdating particleInCell %d:\n", cParticle_I );
			}	

			for ( dim_I=0; dim_I < 3; dim_I++ ) {
				movementVector[dim_I] = ( attractorPoint[dim_I] - (*oldCoord)[dim_I] ) / 3;
				newParticleCoord[dim_I] = (*oldCoord)[dim_I] + movementVector[dim_I];
			}

			if( rank == procToWatch ) {
				//printf("\t\tChanging its coords from (%f,%f,%f) to (%f,%f,%f):\n", (*oldCoord)[0], (*oldCoord)[1], (*oldCoord)[2], newParticleCoord[0], newParticleCoord[1], newParticleCoord[2] );
			}		

			for ( dim_I=0; dim_I < 3; dim_I++ ) {
				currParticle->coord[dim_I] = newParticleCoord[dim_I];
			}
		}
	}
}
void SwarmDumpAndLoadSuite_Setup( SwarmDumpAndLoadSuiteData* data ) {
	
	/* MPI Initializations */
	data->comm = MPI_COMM_WORLD;  
	MPI_Comm_rank( data->comm, &data->rank );
	MPI_Comm_size( data->comm, &data->nProcs );
   
}

void SwarmDumpAndLoadSuite_Teardown( SwarmDumpAndLoadSuiteData* data ) {
	/* Destroy stuff */
}

void SwarmDumpAndLoadSuite_TestSwarmDumpAndLoad( SwarmDumpAndLoadSuiteData* data ) {
	unsigned							nDims;
	unsigned							meshSize[3];
	double							minCrds[3];
	double							maxCrds[3];
	ExtensionManager_Register*	extensionMgr_Register;
	SwarmVariable_Register*		swarmVariable_Register;
	Mesh*								mesh;
	Swarm*							swarm;
	ElementCellLayout*			elementCellLayout;
	RandomParticleLayout*		randomParticleLayout;
	ParticleMovementHandler		*movementHandler;
	Dictionary*						dictionary;
	int								procToWatch;
	Stream*							stream;
	Index								dim_I;
	Index								timeStep;
	Coord								attractorPoint;
	Particle_Index					lParticle_I = 0;
	AbstractContext*				context = NULL;
	SwarmDump*						swarmDumper = NULL;
	Swarm*							newSwarm = NULL;
	Swarm*							swarmList[1];
	FileParticleLayout*			fileParticleLayout = NULL;
	Index 							errorCount = 0;
	char								output_file[PCU_PATH_MAX];

	if( data->nProcs >= 2 ) {
		procToWatch = 1;
	}
	else {
		procToWatch = 0;
	}
	if( data->rank == procToWatch ) printf( "Watching rank: %i\n", data->rank );
	
	nDims = 3;
	meshSize[0] = 32;	meshSize[1] = 32;	meshSize[2] = 32;
	minCrds[0] = 0.0; minCrds[1] = 0.0; minCrds[2] = 0.0;
	maxCrds[0] = 1.0; maxCrds[1] = 1.0; maxCrds[2] = 1.0;

	/* Dictionary */
	dictionary = Dictionary_New();
	Dictionary_Add( dictionary, "particlesPerCell", Dictionary_Entry_Value_FromUnsignedInt( 1 ) );
	Dictionary_Add( dictionary, "seed", Dictionary_Entry_Value_FromUnsignedInt( 13 ) );
	Dictionary_Add( dictionary, "shadowDepth", Dictionary_Entry_Value_FromUnsignedInt( 1 ) );
	
	/* Init mesh */
	extensionMgr_Register = ExtensionManager_Register_New();
	swarmVariable_Register = SwarmVariable_Register_New( NULL );
	mesh = SwarmDumpAndLoadSuite_BuildMesh( nDims, meshSize, minCrds, maxCrds, extensionMgr_Register );
	
	/* Configure the element-cell-layout */
	elementCellLayout = ElementCellLayout_New( "elementCellLayout", mesh );
	
	/* Configure the random-particle-layout */
	randomParticleLayout = RandomParticleLayout_New( "randomParticleLayout", 1, 13 );

	/* Configure the swarm */	
	swarm = Swarm_New( "testSwarm", elementCellLayout, randomParticleLayout, 3, sizeof(Particle), extensionMgr_Register, NULL, data->comm, NULL );
	
	movementHandler = ParticleMovementHandler_New( "movementHandler", True );
	Swarm_AddCommHandler( swarm, movementHandler );
	Stream_Enable( Journal_Register( Info_Type, ParticleMovementHandler_Type ), False );

	/* Build the mesh, swarm */
	Stg_Component_Build( mesh, 0, False );
	Stg_Component_Build( swarm, 0, False );
	
	/* Initialize the mesh, swarm */
	Stg_Component_Initialise( mesh, 0, False );
	Stg_Component_Initialise( swarm, 0, False );

	/* Stream */
	stream = Journal_Register( Info_Type, "SwarmADumpAndLoadStream" );


	
	for ( dim_I=0; dim_I < 3; dim_I++ ) {
		attractorPoint[dim_I] = ( maxCrds[dim_I] - minCrds[dim_I] ) / 3;
	}
	if( data->rank == procToWatch ) {
		printf("Calculated attractor point is at (%f,%f,%f):\n", attractorPoint[0], attractorPoint[1], attractorPoint[2] );
	}
	for ( lParticle_I = 0; lParticle_I < swarm->particleLocalCount; lParticle_I++ ) {
		swarm->particles[lParticle_I].testValue = rand() % 1000;
	}
	/* Start a sample app, where each timestep we move the particles towards the attractor point */
	for ( timeStep=1; timeStep <= 2; timeStep++ ) {
		if( data->rank == procToWatch ) {
			printf("\nStarting timestep %d:\n", timeStep );
		}	

		SwarmDumpAndLoadSuite_UpdateParticlePositionsTowardsAttractor( swarm, attractorPoint, data->rank, procToWatch );
		Swarm_UpdateAllParticleOwners( swarm );
	}

	Dictionary_Add( dictionary, "outputPath", Dictionary_Entry_Value_FromString( "./output" ) );

	/* Now we dump the swarm values, then create a new swarm and load the dumped values onto it,
		and check to see that they're the same */
	context = _AbstractContext_New( 
		sizeof(AbstractContext),
		AbstractContext_Type,
		_AbstractContext_Delete,
		_AbstractContext_Print,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		NULL,
		"testContext",
		True,
		NULL,
		0,
		0,
		data->comm,
		dictionary );

   Stg_ObjectList_ClassAppend( context->CF->registerRegister, (void*)extensionMgr_Register, "ExtensionManager_Register" );
   Stg_ObjectList_ClassAppend( context->CF->registerRegister, (void*)swarmVariable_Register, "SwarmVariable_Register" );
	
	swarmList[0] = swarm;
	swarmDumper = SwarmDump_New( "swarmDumper", context, swarmList, 1, True );
	SwarmDump_Execute( swarmDumper, context );

	sprintf( output_file, "%s/%s.%05d.dat", context->outputPath, swarm->name, context->timeStep ); 
	fileParticleLayout = FileParticleLayout_New( "fileParticleLayout", output_file );
	newSwarm = Swarm_New( "testSwarm2", elementCellLayout, fileParticleLayout, 3, sizeof(Particle), extensionMgr_Register, NULL, data->comm, NULL );
	Stg_Component_Build( newSwarm, 0, False );
	Stg_Component_Initialise( newSwarm, 0, False );

	assert( newSwarm->particleLocalCount == swarm->particleLocalCount );

	if( data->rank == procToWatch ) {
		printf( "\nComparing the %d local Particles between old and new swarms:\n", swarm->particleLocalCount );
		for ( lParticle_I = 0; lParticle_I < swarm->particleLocalCount; lParticle_I++ ) {
			if ( ( swarm->particles[lParticle_I].coord[I_AXIS] != newSwarm->particles[lParticle_I].coord[I_AXIS] )
				|| ( swarm->particles[lParticle_I].coord[J_AXIS] != newSwarm->particles[lParticle_I].coord[J_AXIS] )
				|| ( swarm->particles[lParticle_I].coord[K_AXIS] != newSwarm->particles[lParticle_I].coord[K_AXIS] ) ) {
				printf( "Error: Co-ords at particle %d don't match between old and new swarms.\n", lParticle_I );
				errorCount++;
			}

			if ( swarm->particles[lParticle_I].owningCell != newSwarm->particles[lParticle_I].owningCell ) {
				printf( "Error: owningCell at particle %d doesn't match between old and new swarms.\n", lParticle_I );
				errorCount++;
         }

			if ( ( swarm->particles[lParticle_I].xi[I_AXIS] != newSwarm->particles[lParticle_I].xi[I_AXIS] )
				|| ( swarm->particles[lParticle_I].xi[J_AXIS] != newSwarm->particles[lParticle_I].xi[J_AXIS] )
				|| ( swarm->particles[lParticle_I].xi[K_AXIS] != newSwarm->particles[lParticle_I].xi[K_AXIS] ) ) {
				printf( "Error: Xi values at particle %d don't match between old and new swarms.\n", lParticle_I );
				errorCount++;
			}

			if ( swarm->particles[lParticle_I].testValue != newSwarm->particles[lParticle_I].testValue ) {
				printf( "Error: testValue at particle %d doesn't match between old and new swarms.\n", lParticle_I );
				errorCount++;
			}
		}

		if ( 0 == errorCount ) {
			printf( "\tPassed: swarms are identical.\n" );
		}
		else {
			printf( "\tFailed: %d differences detected.\n", errorCount );
      }
   }
	/* Destroy stuff */
   Stg_Class_Delete( swarm );
   Stg_Class_Delete( randomParticleLayout );
   Stg_Class_Delete( elementCellLayout );
   Stg_Class_Delete( mesh );
   Stg_Class_Delete( extensionMgr_Register );
   Stg_Class_Delete( dictionary );

}

void SwarmDumpAndLoadSuite( pcu_suite_t* suite ) {
	pcu_suite_setData( suite, SwarmDumpAndLoadSuiteData );
	pcu_suite_setFixtures( suite, SwarmDumpAndLoadSuite_Setup, SwarmDumpAndLoadSuite_Teardown );
	pcu_suite_addTest( suite, SwarmDumpAndLoadSuite_TestSwarmDumpAndLoad );
}
