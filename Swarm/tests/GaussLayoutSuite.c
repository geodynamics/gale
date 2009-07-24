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
**   Tests the GaussLayoutSuite
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

#include "GaussLayoutSuite.h"

struct _Particle {
	__IntegrationPoint
};

typedef struct {
	unsigned							nDims;
	unsigned							meshSize[3];
	double							minCrds[3];
	double							maxCrds[3];
	ExtensionManager_Register*	extensionMgr_Register;
	Mesh*								mesh;
	GaussParticleLayout*			gaussParticleLayout;
	ElementCellLayout*			elementCellLayout;
	Swarm*							swarm;
	MPI_Comm       				comm;
   unsigned int   				rank;
   unsigned int   				nProcs;
} GaussLayoutSuiteData;

Mesh* GaussLayoutSuite_BuildMesh( unsigned nDims, unsigned* size, double* minCrds, double* maxCrds, ExtensionManager_Register* emReg ) {
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

void GaussLayoutSuite_Setup( GaussLayoutSuiteData* data ) {
	/* MPI Initializations */
	data->comm = MPI_COMM_WORLD;  
   MPI_Comm_rank( data->comm, &data->rank );
   MPI_Comm_size( data->comm, &data->nProcs );
   
	data->nDims = 3;
	data->meshSize[0] = 2;	data->meshSize[1] = 3;	data->meshSize[2] = 2;
	data->minCrds[0] = 0.0; data->minCrds[1] = 0.0; data->minCrds[2] = 0.0;
	data->maxCrds[0] = 300.0; data->maxCrds[1] = 12.0; data->maxCrds[2] = 300.0;
	
	/* Init mesh */
	data->extensionMgr_Register = ExtensionManager_Register_New();
	data->mesh = GaussLayoutSuite_BuildMesh( data->nDims, data->meshSize, data->minCrds, data->maxCrds, data->extensionMgr_Register );
	
	/* Configure the element-cell-layout */
	data->elementCellLayout = ElementCellLayout_New( "elementCellLayout", data->mesh );
}

void GaussLayoutSuite_Teardown( GaussLayoutSuiteData* data ) {
	/* Destroy stuff */
	Stg_Class_Delete( data->gaussParticleLayout );
	Stg_Class_Delete( data->elementCellLayout );
	Stg_Class_Delete( data->swarm );
	Stg_Class_Delete( data->mesh );
	Stg_Class_Delete( data->extensionMgr_Register );
	remove( "1ParticlePerDim_3D.dat" );
   remove( "2ParticlesPerDim_3D.dat" );
   remove( "3ParticlesPerDim_3D.dat" );
}

void GaussLayoutSuite_Test1ParticlePerDim_3D( GaussLayoutSuiteData* data ) {
	Cell_PointIndex			count;
	double 						x,y,z;
	unsigned int 				p, i, len;
	int							procToWatch;
	LocalParticle* 			particle;
	Coord 						minCell;
	Coord 						maxCell;
	Particle_InCellIndex 	particlesPerDim[3] = {1, 1, 1};
	Stream*						stream;
	char 							expected_file[PCU_PATH_MAX];
	
	if( data->nProcs >= 2 ) {
		procToWatch = 1;
	}
	else {
		procToWatch = 0;
	}
	if( data->rank == procToWatch ) printf( "Watching rank: %i\n", data->rank );

	stream = Journal_Register( Info_Type, "1ParticlePerDim_3D" );
	Stream_RedirectFile( stream, "1ParticlePerDim_3D.dat" );
	
	data->gaussParticleLayout = GaussParticleLayout_New( "gaussParticleLayout", data->nDims, particlesPerDim );
	/* Configure the swarm */
	data->swarm = Swarm_New( "testGaussSwarm", data->elementCellLayout, data->gaussParticleLayout, 
				data->nDims, sizeof(Particle), data->extensionMgr_Register, NULL, MPI_COMM_WORLD, NULL );
	/* Build the swarm */
	Stg_Component_Build( data->swarm, 0, False );
	Stg_Component_Initialise( data->swarm, 0, False );
	
	len = (int) sizeof( data->swarm->cellParticleCountTbl );
	count = 0;
	
	/* Checks that the particule count on each cell are the same. */
	for( i = 0; i < len; i++ ) {
			count = data->swarm->cellParticleCountTbl[i];
			pcu_check_true( count == 1 );
	}
	Swarm_GetCellMinMaxCoords( data->swarm, 4, minCell, maxCell );
	Journal_Printf( stream, "Particle per dim: 1 1 1\n");	
	
	for( p = 0; p < count; p++ ) {
		particle = (LocalParticle*)Swarm_ParticleInCellAt( data->swarm, 4, p  );

		/* convert to global coords */
		x = 0.5 * ( maxCell[0] - minCell[0] ) * ( particle->xi[0] + 1.0 ) + minCell[0];
		y = 0.5 * ( maxCell[1] - minCell[1] ) * ( particle->xi[1] + 1.0 ) + minCell[1];
		z = 0.5 * ( maxCell[2] - minCell[2] ) * ( particle->xi[2] + 1.0 ) + minCell[2];
		
		printf( "pId=%d : coords = { %.12f, %.12f, %.12f }, xi = { %.12f, %.12f, %.12f }\n", 
			p, x, y, z, particle->xi[0], particle->xi[1], particle->xi[2] );
		Journal_Printf( stream, "pId=%d : coords = { %.12f, %.12f, %.12f }, xi = { %.12f, %.12f, %.12f }\n", 
			p, x, y, z, particle->xi[0], particle->xi[1], particle->xi[2] );	
	}
	pcu_filename_expected( "testGaussLayoutOutput1ParticlePerDim.expected", expected_file );
	pcu_check_fileEq( "1ParticlePerDim_3D.dat", expected_file );
}

void GaussLayoutSuite_Test2ParticlesPerDim_3D( GaussLayoutSuiteData* data ) {
	Cell_PointIndex			count;
	double 						x,y,z;
	unsigned int				p, i, len;
	int							procToWatch;
	LocalParticle* 			particle;
	Coord 						minCell;
	Coord 						maxCell;
	Particle_InCellIndex 	particlesPerDim[3] = {2, 2, 2};
	double						xi[8][3];
	double						coord[8][3];
	Stream*						stream;
	char 							expected_file[PCU_PATH_MAX];
	
	if( data->nProcs >= 2 ) {
		procToWatch = 1;
	}
	else {
		procToWatch = 0;
	}
	if( data->rank == procToWatch ) printf( "Watching rank: %i\n", data->rank );

	stream = Journal_Register( Info_Type, "2ParticlesPerDim_3D" );
	Stream_RedirectFile( stream, "2ParticlesPerDim_3D.dat" );
	
	data->gaussParticleLayout = GaussParticleLayout_New( "gaussParticleLayout", data->nDims, particlesPerDim );
	/* Configure the swarm */
	data->swarm = Swarm_New( "testGaussSwarm", data->elementCellLayout, data->gaussParticleLayout, 
				data->nDims, sizeof(Particle), data->extensionMgr_Register, NULL, MPI_COMM_WORLD, NULL );
	/* Build the swarm */
	Stg_Component_Build( data->swarm, 0, False );
	Stg_Component_Initialise( data->swarm, 0, False );

	len = (int) sizeof( data->swarm->cellParticleCountTbl );
	count = 0;
	
	/* Checks that the particule count on each cell are the same. */
	for( i = 0; i < len; i++ ) {
		count = data->swarm->cellParticleCountTbl[i];
		pcu_check_true( count == 8 );
	}
	Swarm_GetCellMinMaxCoords( data->swarm, 4, minCell, maxCell );
	Journal_Printf( stream, "Particle per dim: 2 2 2\n");	

	for( p = 0; p < count; p++ ) {
		particle = (LocalParticle*)Swarm_ParticleInCellAt( data->swarm, 4, p  );

		/* convert to global coords */
		x = 0.5 * ( maxCell[0] - minCell[0] ) * ( particle->xi[0] + 1.0 ) + minCell[0];
		y = 0.5 * ( maxCell[1] - minCell[1] ) * ( particle->xi[1] + 1.0 ) + minCell[1];
		z = 0.5 * ( maxCell[2] - minCell[2] ) * ( particle->xi[2] + 1.0 ) + minCell[2];
		
		printf( "pId=%d : coords = { %.12f, %.12f, %.12f }, xi = { %.12f, %.12f, %.12f }\n", 
			p, x, y, z, particle->xi[0], particle->xi[1], particle->xi[2] );
		Journal_Printf( stream, "pId=%d : coords = { %.12f, %.12f, %.12f }, xi = { %.12f, %.12f, %.12f }\n", 
			p, x, y, z, particle->xi[0], particle->xi[1], particle->xi[2] );	
	}
	pcu_filename_expected( "testGaussLayoutOutput2ParticlesPerDim.expected", expected_file );
	pcu_check_fileEq( "2ParticlesPerDim_3D.dat", expected_file );
}


void GaussLayoutSuite_Test3ParticlesPerDim_3D( GaussLayoutSuiteData* data ) {
	Cell_PointIndex			count;
	double		 				x,y,z;
	unsigned int				p, i, len;
	int							procToWatch;
	LocalParticle* 			particle;
	Coord 						minCell;
	Coord 						maxCell;
	Particle_InCellIndex 	particlesPerDim[3] = {3, 3, 3};
	double						xi[27][3];
	double						coord[27][3];
	Stream*						stream;
	char 							expected_file[PCU_PATH_MAX];
	
	if( data->nProcs >= 2 ) {
		procToWatch = 1;
	}
	else {
		procToWatch = 0;
	}
	if( data->rank == procToWatch ) printf( "Watching rank: %i\n", data->rank );
	
	stream = Journal_Register( Info_Type, "3ParticlesPerDim_3D" );
	Stream_RedirectFile( stream, "3ParticlesPerDim_3D.dat" );
	
	data->gaussParticleLayout = GaussParticleLayout_New( "gaussParticleLayout", data->nDims, particlesPerDim );
	/* Configure the swarm */
	data->swarm = Swarm_New( "testGaussSwarm", data->elementCellLayout, data->gaussParticleLayout, 
				data->nDims, sizeof(Particle), data->extensionMgr_Register, NULL, MPI_COMM_WORLD, NULL );
	/* Build the swarm */
	Stg_Component_Build( data->swarm, 0, False );
	Stg_Component_Initialise( data->swarm, 0, False );
		
	len = (int) sizeof( data->swarm->cellParticleCountTbl );
	count = 0;
	
	/* Checks that the particule count on each cell are the same. */
	for( i = 0; i < len; i++ ) {
		count = data->swarm->cellParticleCountTbl[i];
			pcu_check_true( count == 27 );
	}	
	Swarm_GetCellMinMaxCoords( data->swarm, 4, minCell, maxCell );
	Journal_Printf( stream, "Particle per dim: 3 3 3\n");	
			
	for( p = 0; p < count; p++ ) {
		particle = (LocalParticle*)Swarm_ParticleInCellAt( data->swarm, 4, p  );

		/* convert to global coords */
		x = 0.5 * ( maxCell[0] - minCell[0] ) * ( particle->xi[0] + 1.0 ) + minCell[0];
		y = 0.5 * ( maxCell[1] - minCell[1] ) * ( particle->xi[1] + 1.0 ) + minCell[1];
		z = 0.5 * ( maxCell[2] - minCell[2] ) * ( particle->xi[2] + 1.0 ) + minCell[2];
		
		printf( "pId=%d : coords = { %.12f, %.12f, %.12f }, xi = { %.12f, %.12f, %.12f }\n", 
			p, x, y, z, particle->xi[0], particle->xi[1], particle->xi[2] );
		Journal_Printf( stream, "pId=%d : coords = { %.12f, %.12f, %.12f }, xi = { %.12f, %.12f, %.12f }\n", 
			p, x, y, z, particle->xi[0], particle->xi[1], particle->xi[2] );	
	}
	pcu_filename_expected( "testGaussLayoutOutput3ParticlesPerDim.expected", expected_file );
	pcu_check_fileEq( "3ParticlesPerDim_3D.dat", expected_file );
}

void GaussLayoutSuite( pcu_suite_t* suite ) {
   pcu_suite_setData( suite, GaussLayoutSuiteData );
   pcu_suite_setFixtures( suite, GaussLayoutSuite_Setup, GaussLayoutSuite_Teardown );
   pcu_suite_addTest( suite, GaussLayoutSuite_Test1ParticlePerDim_3D );
   pcu_suite_addTest( suite, GaussLayoutSuite_Test2ParticlesPerDim_3D );
   pcu_suite_addTest( suite, GaussLayoutSuite_Test3ParticlesPerDim_3D );
}
