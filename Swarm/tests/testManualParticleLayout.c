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
**	Test that the ManualParticleLayout initialises particle positions properly.
**
** Assumptions:
**	None as yet.
**
** Comments:
**	None as yet.
**
** $Id: testManualParticleLayout.c 4175 2007-08-16 03:39:26Z DavidLee $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>

#include <StgDomain/Geometry/Geometry.h>
#include <StgDomain/Shape/Shape.h>
#include <StgDomain/Mesh/Mesh.h>
#include <StgDomain/Utils/Utils.h>
#include <StgDomain/Swarm/Swarm.h>

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

	gen = CartesianGenerator_New( "" );
	CartesianGenerator_SetDimSize( gen, nDims );
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

int main( int argc, char* argv[] ) {
	MPI_Comm                    CommWorld;
	int                         rank;
	int                         numProcessors;
	int                         procToWatch;
	Dictionary*                 dictionary;
	unsigned	nDims = 3;
	unsigned	meshSize[3] = {4, 2, 1};
	double		minCrds[3] = {0.0, 0.0, 0.0};
	double		maxCrds[3] = {1.0, 1.0, 1.0};
	ExtensionManager_Register*  extensionMgr_Register;
	Mesh*                       mesh;
	ManualParticleLayout*       particleLayout;
	ElementCellLayout*          elementCellLayout;
	Swarm*                      swarm;
	Stream*                     stream;
	Dimension_Index             dim;
	Dictionary_Entry_Value*     particlePositionsList;
	Dictionary_Entry_Value*     particlePositionEntry;
	
	/* Initialise MPI, get world info */
	MPI_Init( &argc, &argv );
	MPI_Comm_dup( MPI_COMM_WORLD, &CommWorld );
	MPI_Comm_size( CommWorld, &numProcessors );
	MPI_Comm_rank( CommWorld, &rank );
	
	StGermain_Init( &argc, &argv );
	
	StgDomainGeometry_Init( &argc, &argv );
	StgDomainShape_Init( &argc, &argv );
	StgDomainMesh_Init( &argc, &argv );
	StgDomainUtils_Init( &argc, &argv );
	StgDomainSwarm_Init( &argc, &argv );
	MPI_Barrier( CommWorld ); /* Ensures copyright info always come first in output */

	stream = Journal_Register (Info_Type, "myStream");

	if( argc >= 2 ) {
		procToWatch = atoi( argv[1] );
	}
	else {
		procToWatch = 0;
	}
/*	if( rank == procToWatch ) printf( "Watching rank: %i\n", rank );*/
	
	/* Read input */
	dictionary = Dictionary_New();
	particlePositionsList = Dictionary_Entry_Value_NewList();
	Dictionary_Add( dictionary, "manualParticlePositions", particlePositionsList );
	particlePositionEntry = Dictionary_Entry_Value_NewStruct();
	Dictionary_Entry_Value_AddElement( particlePositionsList, particlePositionEntry );
	Dictionary_Entry_Value_AddMember( particlePositionEntry, "x", Dictionary_Entry_Value_FromDouble( 0.4 ) );
	Dictionary_Entry_Value_AddMember( particlePositionEntry, "y", Dictionary_Entry_Value_FromDouble( 0.3 ) );
	Dictionary_Entry_Value_AddMember( particlePositionEntry, "z", Dictionary_Entry_Value_FromDouble( 0.2 ) );
	particlePositionEntry = Dictionary_Entry_Value_NewStruct();
	Dictionary_Entry_Value_AddElement( particlePositionsList, particlePositionEntry );
	Dictionary_Entry_Value_AddMember( particlePositionEntry, "x", Dictionary_Entry_Value_FromDouble( 0.7 ) );
	Dictionary_Entry_Value_AddMember( particlePositionEntry, "y", Dictionary_Entry_Value_FromDouble( 0.6 ) );
	Dictionary_Entry_Value_AddMember( particlePositionEntry, "z", Dictionary_Entry_Value_FromDouble( 0.5 ) );
	particlePositionEntry = Dictionary_Entry_Value_NewStruct();
	Dictionary_Entry_Value_AddElement( particlePositionsList, particlePositionEntry );
	Dictionary_Entry_Value_AddMember( particlePositionEntry, "x", Dictionary_Entry_Value_FromDouble( 0.8 ) );
	Dictionary_Entry_Value_AddMember( particlePositionEntry, "y", Dictionary_Entry_Value_FromDouble( 0.1 ) );
	Dictionary_Entry_Value_AddMember( particlePositionEntry, "z", Dictionary_Entry_Value_FromDouble( 0.3 ) );
	particlePositionEntry = Dictionary_Entry_Value_NewStruct();
	Dictionary_Entry_Value_AddElement( particlePositionsList, particlePositionEntry );
	Dictionary_Entry_Value_AddMember( particlePositionEntry, "x", Dictionary_Entry_Value_FromDouble( 0.9 ) );
	Dictionary_Entry_Value_AddMember( particlePositionEntry, "y", Dictionary_Entry_Value_FromDouble( 0.4 ) );
	Dictionary_Entry_Value_AddMember( particlePositionEntry, "z", Dictionary_Entry_Value_FromDouble( 0.1 ) );
	
	/* Init mesh */
	extensionMgr_Register = ExtensionManager_Register_New();
	mesh = buildMesh( nDims, meshSize, minCrds, maxCrds, extensionMgr_Register );
	
	/* Configure the element-cell-layout */
	elementCellLayout = ElementCellLayout_New( "elementCellLayout", mesh );
	
	/* Build the mesh */
	Stg_Component_Build( mesh, 0, False );
	Stg_Component_Initialise( mesh, 0, False );
	
	/* Configure the gauss-particle-layout */
	particleLayout = ManualParticleLayout_New( "manualParticleLayout", dictionary );
	
	/* Configure the swarm */
	swarm = Swarm_New(  "testSwarm", elementCellLayout, particleLayout, dim, sizeof(Particle),
		extensionMgr_Register, NULL, CommWorld, NULL );
	
	/* Build the swarm */
	Stg_Component_Build( swarm, 0, False );
	Stg_Component_Initialise( swarm, 0, False );
	
	if( rank == procToWatch ) {
		Stg_Class_Print( particleLayout, stream );
		/* Print out the particles on all cells */
		Swarm_PrintParticleCoords_ByCell( swarm, stream );
	}
	

	/* Destroy stuff */
	Stg_Class_Delete( particleLayout );
	Stg_Class_Delete( elementCellLayout );
	Stg_Class_Delete( swarm );
	Stg_Class_Delete( mesh );
	Stg_Class_Delete( extensionMgr_Register );
	Stg_Class_Delete( dictionary );
	
	StgDomainSwarm_Finalise();
	StgDomainUtils_Finalise();
	StgDomainMesh_Finalise();
	StgDomainShape_Finalise();
	StgDomainGeometry_Finalise();
	
	StGermain_Finalise();
	
	/* Close off MPI */
	MPI_Finalize();
	
	return 0; /* success */
}
