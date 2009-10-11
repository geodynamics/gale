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
**   Tests the SpaceFillerParticleLayoutSuite
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

#include "SpaceFillerParticleLayoutSuite.h"

struct _Particle {
	__IntegrationPoint
};

typedef struct {
	MPI_Comm comm;
	unsigned rank;
	unsigned nProcs;
} SpaceFillerParticleLayoutSuiteData;

Mesh* SpaceFillerParticleLayoutSuite_BuildMesh( unsigned nDims, unsigned* size, double* minCrds, double* maxCrds, ExtensionManager_Register* emReg ) {
	CartesianGenerator* gen;
	Mesh*						mesh;

	gen = CartesianGenerator_New( "" );
	CartesianGenerator_SetDimSize( gen, nDims );
	CartesianGenerator_SetTopologyParams( gen, size, 0, NULL, NULL );
	CartesianGenerator_SetGeometryParams( gen, minCrds, maxCrds );

	mesh = Mesh_New( "" );
	Mesh_SetExtensionManagerRegister( mesh, emReg );
	Mesh_SetGenerator( mesh, gen );
	Mesh_SetAlgorithms( mesh, Mesh_RegularAlgorithms_New( "" ) );

	Stg_Component_Build( mesh, NULL, False );
	Stg_Component_Initialise( mesh, NULL, False );

	KillObject( mesh->generator );

	return mesh;
}

void SpaceFillerParticleLayoutSuite_Setup( SpaceFillerParticleLayoutSuiteData* data ) {
	/* MPI Initializations */
	data->comm = MPI_COMM_WORLD;  
	MPI_Comm_rank( data->comm, &data->rank );
	MPI_Comm_size( data->comm, &data->nProcs );
   
}

void SpaceFillerParticleLayoutSuite_Teardown( SpaceFillerParticleLayoutSuiteData* data ) {
}

void SpaceFillerParticleLayoutSuite_TestSpaceFillerParticle( SpaceFillerParticleLayoutSuiteData* data ) {
	ExtensionManager_Register*	extensionMgr_Register;
	SpaceFillerParticleLayout*	particleLayout;
	ElementCellLayout*			elementCellLayout;
	Dimension_Index				dim;
	Mesh*								mesh;
	Swarm*							swarm;
	Stream*							stream;
	unsigned							nDims;
	unsigned							meshSize[3];
	double							minCrds[3];
	double							maxCrds[3];
	int								procToWatch = data->nProcs > 1 ? 1 : 0;
	char								expected_file[PCU_PATH_MAX];

	if( data->rank == procToWatch ) {
		stream = Journal_Register( Info_Type, "TestSpaceFillerParticle" );
		Stream_RedirectFile( stream, "spaceFillerParticle.dat" );

		nDims = 3;
		meshSize[0] = 4;	meshSize[1] = 2;	meshSize[2] = 1;
		minCrds[0] = 0.0; minCrds[1] = 0.0; minCrds[2] = 0.0;
		maxCrds[0] = 400.0; maxCrds[1] = 200.0; maxCrds[2] = 100.0;

		extensionMgr_Register = ExtensionManager_Register_New();
		mesh = SpaceFillerParticleLayoutSuite_BuildMesh( nDims, meshSize, minCrds, maxCrds, extensionMgr_Register );
		
		elementCellLayout = ElementCellLayout_New( "spaceFillerParticlElementCellLayout", mesh );
		particleLayout = SpaceFillerParticleLayout_New( "spaceFillerParticleLayout", nDims, SpaceFillerParticleLayout_Invalid, 20 );
	
		swarm = Swarm_New( "testSpaceFIllerParticle", elementCellLayout, particleLayout, dim, sizeof(Particle),
			extensionMgr_Register, NULL, data->comm, NULL );
 
		Stg_Component_Build( swarm, 0, False );
		Stg_Component_Initialise( swarm, 0, False );

		Swarm_PrintParticleCoords_ByCell( swarm, stream );

		pcu_filename_expected( "testSpaceFillerParticleLayoutOutput.expected", expected_file );
		pcu_check_fileEq( "spaceFillerParticle.dat", expected_file );
		remove( "spaceFillerParticle.dat" );

		Stg_Class_Delete( extensionMgr_Register );
		Stg_Component_Destroy( mesh, NULL, True );
		Stg_Component_Destroy( elementCellLayout, NULL, True );
		Stg_Component_Destroy( particleLayout, NULL, True );
		Stg_Component_Destroy( swarm, NULL, True );
	}
}

void SpaceFillerParticleLayoutSuite( pcu_suite_t* suite ) {
	pcu_suite_setData( suite, SpaceFillerParticleLayoutSuiteData );
	pcu_suite_setFixtures( suite, SpaceFillerParticleLayoutSuite_Setup, SpaceFillerParticleLayoutSuite_Teardown );
	pcu_suite_addTest( suite, SpaceFillerParticleLayoutSuite_TestSpaceFillerParticle );
}
