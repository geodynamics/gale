/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003-2006, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street,
**	Melbourne, 3053, Australia.
** Copyright (c) 2005-2006, Monash Cluster Computing, Building 28, Monash University Clayton Campus,
**	Victoria, 3800, Australia
**
** Primary Contributing Organisations:
**	Victorian Partnership for Advanced Computing Ltd, Computational Software Development - http://csd.vpac.org
**	Australian Computational Earth Systems Simulator - http://www.access.edu.au
**	Monash Cluster Computing - http://www.mcc.monash.edu.au
**
** Contributors:
**	Robert Turnbull, Research Assistant, Monash University. (robert.turnbull@sci.monash.edu.au)
**	Patrick D. Sunter, Software Engineer, VPAC. (patrick@vpac.org)
**	Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
**	Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**	David May, PhD Student, Monash University (david.may@sci.monash.edu.au)
**	Vincent Lemiale, Postdoctoral Fellow, Monash University. (vincent.lemiale@sci.monash.edu.au)
**	Julian Giordani, Research Assistant, Monash University. (julian.giordani@sci.monash.edu.au)
**	Louis Moresi, Associate Professor, Monash University. (louis.moresi@sci.monash.edu.au)
**	Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**	Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
**	David Stegman, Postdoctoral Fellow, Monash University. (david.stegman@sci.monash.edu.au)
**	Wendy Sharples, PhD Student, Monash University (wendy.sharples@sci.monash.edu.au)
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
** $Id: testPeriodicBCs.c 503 2007-08-02 08:41:48Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgFEM/StgFEM.h>

#include <PICellerator/Voronoi/Voronoi.h>
#include <PICellerator/PopulationControl/PopulationControl.h>
#include <PICellerator/Weights/Weights.h>

#include "MaterialPoints/MaterialPoints.h"

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

void UpdateParticlePositionsToLeft( 
		Swarm* swarm,
		Processor_Index rank,
		Processor_Index procToWatch );

FeMesh* configureFeMesh( double minCrd[3], double maxCrd[3], unsigned int sizes[3], unsigned nDims ) {
	CartesianGenerator*	gen;
	FeMesh*			feMesh;
	unsigned		maxDecomp[3] = {0, 1, 1};

	sizes[0] = sizes[1] = 6;
	sizes[2] = 1;
	minCrd[0] = minCrd[1] = minCrd[2] = 0.0;
	maxCrd[0] = maxCrd[1] = maxCrd[2] = 1.2;

	gen = CartesianGenerator_New( "" );
	CartesianGenerator_SetDimSize( gen, nDims );
	CartesianGenerator_SetTopologyParams( gen, sizes, 0, NULL, maxDecomp );
	CartesianGenerator_SetGeometryParams( gen, minCrd, maxCrd );
	CartesianGenerator_SetShadowDepth( gen, 0 );

	feMesh = FeMesh_New( "" );
	Mesh_SetGenerator( feMesh, gen );
	FeMesh_SetElementFamily( feMesh, "linear" );

	return feMesh;
}


int main( int argc, char* argv[] ) {
	MPI_Comm			CommWorld;
	int				rank;
	int				numProcessors;
	int				procToWatch;
	Dictionary*			dictionary;
	ExtensionManager_Register*	extensionMgr_Register;
	Mesh*				mesh;
	ElementCellLayout*		elementCellLayout;
	RandomParticleLayout*		randomParticleLayout;
	Swarm*				swarm;
	Stream*				stream;
	Index				timeStep;
	EntryPoint_Register*		epRegister;
	PeriodicBoundariesManager*	perBCsManager;
	Index                           decompDims;
	char*                           directory;	
	unsigned                        numElements[3];
	double                          minCoords[3];
	double                          maxCoords[3];
	unsigned int                    dim = 2;

	/* Initialise MPI, get world info */
	MPI_Init( &argc, &argv );
	MPI_Comm_dup( MPI_COMM_WORLD, &CommWorld );
	MPI_Comm_size( CommWorld, &numProcessors );
	MPI_Comm_rank( CommWorld, &rank );
	
	StGermain_Init( &argc, &argv );
	StgFEM_Init( &argc, &argv );
        PICellerator_Voronoi_Init( &argc, &argv );
        PICellerator_PopulationControl_Init( &argc, &argv );
        PICellerator_Weights_Init( &argc, &argv );
        PICellerator_MaterialPoints_Init( &argc, &argv );	

	MPI_Barrier( CommWorld ); /* Ensures copyright info always come first in output */
	/* Add the PICellerator path to the global xml path dictionary */
	directory = Memory_Alloc_Array( char, 200, "xmlDirectory" ) ;
	sprintf(directory, "%s%s", LIB_DIR, "/StGermain" );
	XML_IO_Handler_AddDirectory( "PICellerator", directory );
	Memory_Free(directory);
						
	/* Add the plugin path to the global plugin list */
	ModulesManager_AddDirectory( "PICellerator", LIB_DIR );

	stream = Journal_Register (Info_Type, "myStream");

	/* *** Journal stuff *** */
	Journal_Enable_TypedStream( DebugStream_Type, False );
	Stream_EnableBranch( Journal_Register( Info_Type, ParticleCommHandler_Type ), False );
	/*
	perBCsManagerStream = Journal_Register( Debug_Type, PeriodicBoundariesManager_Type );
	Stream_EnableBranch( perBCsManagerStream, True );
	Stream_SetLevelBranch( perBCsManagerStream, 3 );
	perBCsManagerStream = Journal_Register( Debug_Type, Swarm_Type );
	Stream_EnableBranch( perBCsManagerStream, True );
	Stream_SetLevelBranch( perBCsManagerStream, 3 );
	perBCsManagerStream = Journal_Register( Debug_Type, Mesh_Type );
	Stream_EnableBranch( perBCsManagerStream, True );
	Stream_SetLevelBranch( perBCsManagerStream, 3 );
	*/

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
	numElements[0] = 12;
	numElements[1] = 12;
	numElements[2] = 1;
	minCoords[0] = 0.0f;
	minCoords[1] = 0.0f;
	minCoords[2] = 0.0f;
	maxCoords[0] = 1.0f;
	maxCoords[1] = 1.0f;
	maxCoords[2] = 1.0f;
	Dictionary_Add( dictionary, "particlesPerCell", Dictionary_Entry_Value_FromUnsignedInt( 1 ) );
	Dictionary_Add( dictionary, "seed", Dictionary_Entry_Value_FromUnsignedInt( 13 ) );
	/*  TODO: a 2nd test with the periodic shadowing enabled. Its handy to keep the orig one */
	/* 	without it though. */
	/* Dictionary_Add( dictionary, "isPeriodicI", Dictionary_Entry_Value_FromBool( True ) ); */
	/* Dictionary_Add( dictionary, "isPeriodicJ", Dictionary_Entry_Value_FromBool( True ) ); */
	decompDims = 1;

	/* Run the mesher */
	mesh = (Mesh*)configureFeMesh( minCoords, maxCoords, numElements, dim );

	/* Init mesh */
	extensionMgr_Register = ExtensionManager_Register_New();
	
	/* Configure the element-cell-layout */
	elementCellLayout = ElementCellLayout_New( "elementCellLayout", mesh );
	
	/* Configure the random-particle-layout */
	randomParticleLayout = RandomParticleLayout_New( "randomParticleLayout", 1, 13 );
	
	/* Configure the swarm */
	swarm = Swarm_New( "testSwarm", elementCellLayout, randomParticleLayout, 3, sizeof(Particle),
		extensionMgr_Register, NULL, CommWorld );
	
	epRegister = EntryPoint_Register_New();
	
	/* Configure the perBCs manager */
	perBCsManager = PeriodicBoundariesManager_New( "perBCsManager", mesh, swarm,
		dictionary );
	
	/* +++ BUILD PHASE +++ */
	
	/* Build the mesh */
	Stg_Component_Build( mesh, 0, False );
	/* Build the swarm */
	Stg_Component_Build( swarm, 0, False );

	Stg_Component_Build( perBCsManager, 0, False );
	PeriodicBoundariesManager_AddPeriodicBoundary( perBCsManager, I_AXIS );
	/* +++ INITIALISE PHASE +++ */

	Stg_Component_Initialise( mesh, 0, False );
	Stg_Component_Initialise( swarm, 0, False );
	Stg_Component_Initialise( perBCsManager, 0, False );
	
	if( rank == procToWatch ) {
		Stg_Class_Print( swarm, stream );
	}	

	/* +++ RUN PHASE +++ */

	/* Start a sample app, where each timestep we move the particles towards the attractor point */
	for ( timeStep=1; timeStep <= 5; timeStep++ ) {
		Particle_Index lParticle_I;
		
		if( rank == procToWatch ) {
			printf("\nStarting timestep %d:\n", timeStep );
		}	

		UpdateParticlePositionsToLeft( swarm, rank, procToWatch );
	
		if( rank == procToWatch ) {
			printf("\nUpdating periodic BCs\n" );
		}	
		for ( lParticle_I = 0; lParticle_I < swarm->particleLocalCount ; lParticle_I++ ) {
			PeriodicBoundariesManager_UpdateParticle( perBCsManager, lParticle_I );
		}
		
		Swarm_UpdateAllParticleOwners( swarm );
	}
	
	/* Destroy stuff */
	Stg_Class_Delete( swarm );
	Stg_Class_Delete( randomParticleLayout );
	Stg_Class_Delete( elementCellLayout );
	Stg_Class_Delete( mesh );
	Stg_Class_Delete( extensionMgr_Register );
	Stg_Class_Delete( dictionary );
	
        PICellerator_MaterialPoints_Finalise();
        PICellerator_Weights_Finalise();
        PICellerator_PopulationControl_Finalise();
        PICellerator_Voronoi_Finalise();
	StgFEM_Finalise();
	StGermain_Finalise();

	/* Close off MPI */
	MPI_Finalize();
	
	return 0; /* success */
}


void UpdateParticlePositionsToLeft( 
		Swarm* swarm,
		Processor_Index rank,
		Processor_Index procToWatch )
{
	Cell_LocalIndex			lCell_I;
	Particle_InCellIndex		cParticle_I;
	GlobalParticle* 		currParticle;
	Index				dim_I;

	for ( lCell_I=0; lCell_I < swarm->cellLocalCount; lCell_I++ ) {
		if( rank == procToWatch ) {
			printf("\tUpdating Particles positions in local cell %d:\n", lCell_I );
		}	
		for ( cParticle_I=0; cParticle_I < swarm->cellParticleCountTbl[lCell_I]; cParticle_I++ ) {
			Coord movementVector = {0,0,0};
			Coord newParticleCoord = {0,0,0};
			Coord* oldCoord;

			currParticle = (GlobalParticle*)Swarm_ParticleInCellAt( swarm, lCell_I, cParticle_I );
			oldCoord = &currParticle->coord;
			if( rank == procToWatch ) {
				printf("\t\tUpdating particleInCell %d:\n", cParticle_I );
			}	

			movementVector[I_AXIS] = -0.1;

			for ( dim_I=0; dim_I < 3; dim_I++ ) {
				newParticleCoord[dim_I] = (*oldCoord)[dim_I] + movementVector[dim_I];
			}

			if( rank == procToWatch ) {
				printf("\t\tChanging its coords from (%f,%f,%f) to (%f,%f,%f):\n",
					(*oldCoord)[0], (*oldCoord)[1], (*oldCoord)[2],
					newParticleCoord[0], newParticleCoord[1], newParticleCoord[2] );
			}		

			for ( dim_I=0; dim_I < 3; dim_I++ ) {
				currParticle->coord[dim_I] = newParticleCoord[dim_I];
			}
		}
	}
}
