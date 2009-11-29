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
**   Tests the WithinShapeParticleLayoutSuite
**
** $Id: testTemplate.c 3462 2006-02-19 06:53:24Z WalterLandry $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include "pcu/pcu.h"
#include <StGermain/StGermain.h>
#include "StgDomain/Geometry/Geometry.h"
#include "StgDomain/Shape/Shape.h"
#include "StgDomain/Mesh/Mesh.h"
#include "StgDomain/Utils/Utils.h"
#include "StgDomain/Swarm/Swarm.h"

#include "WithinShapeParticleLayoutSuite.h"

struct _Particle {
	__GlobalParticle
};

typedef struct {
	MPI_Comm comm;
	unsigned rank;
	unsigned nProcs;
} WithinShapeParticleLayoutSuiteData;

Mesh* WithinShapeParticleLayoutSuite_buildMesh( unsigned nDims, unsigned* size, double* minCrds, double* maxCrds, ExtensionManager_Register* emReg ) {
	CartesianGenerator*	gen;
	Mesh*						mesh;

	gen = CartesianGenerator_New( "", NULL );
	CartesianGenerator_SetDimSize( gen, nDims );
	CartesianGenerator_SetTopologyParams( gen, size, 0, NULL, NULL );
	CartesianGenerator_SetGeometryParams( gen, minCrds, maxCrds );

	mesh = Mesh_New( "", NULL );
	Mesh_SetExtensionManagerRegister( mesh, emReg );
	Mesh_SetGenerator( mesh, gen );

	Stg_Component_Build( mesh, NULL, False );
	Stg_Component_Initialise( mesh, NULL, False );

	FreeObject( mesh->generator );

	return mesh;
}


void WithinShapeParticleLayoutSuite_Setup( WithinShapeParticleLayoutSuiteData* data ) {
	/* MPI Initializations */	
	data->comm = MPI_COMM_WORLD;  
	MPI_Comm_rank( data->comm, &data->rank );
	MPI_Comm_size( data->comm, &data->nProcs );

}

void WithinShapeParticleLayoutSuite_Teardown( WithinShapeParticleLayoutSuiteData* data ) {
}

void WithinShapeParticleLayoutSuite_TestWithinShapeSphere( WithinShapeParticleLayoutSuiteData* data ) {
	ExtensionManager_Register*	extensionMgr_Register;
	WithinShapeParticleLayout*	particleLayout;
	ElementCellLayout*			elementCellLayout;
	Mesh*								mesh;
	Swarm*							swarm;
	unsigned							nDims;
	unsigned							meshSize[3];
	double							minCrds[3];
	double							maxCrds[3];
	XYZ								centre;
	Stg_Shape*						shape;
	int								procToWatch;
	int								particleCount = 10;
	Index								i;

	procToWatch = data->nProcs > 1 ? 1 : 0;	

	if( data->rank == procToWatch ) {
		nDims = 3;
		meshSize[0] = 5; meshSize[1] = 3; meshSize[2] = 2;
		minCrds[0] = 0; minCrds[1] = 0; minCrds[2] = 0;
		maxCrds[0] = 1; maxCrds[1] = 1; maxCrds[2] = 1;	

		/* Init mesh */
		extensionMgr_Register = ExtensionManager_Register_New();
		mesh = WithinShapeParticleLayoutSuite_buildMesh( nDims, meshSize, minCrds, maxCrds, extensionMgr_Register );

		/* Configure the element-cell-layout */
		elementCellLayout = ElementCellLayout_New( "elementCellLayout", NULL, mesh );

		/* Build the mesh */
		Stg_Component_Build( mesh, 0, False );
		Stg_Component_Initialise( mesh, 0, False );	

		centre[I_AXIS] = centre[J_AXIS] = centre[K_AXIS] = 0.5;
		shape = (Stg_Shape*)Sphere_New( "testSphere", nDims, centre, 0, 0, 0, 0.05 );

		/* Configure the gauss-particle-layout */
		particleLayout = WithinShapeParticleLayout_New( "withinShapeParticleLayoutSphere", NULL, 
			GlobalCoordSystem, False, particleCount, 0.0, nDims, shape );

		/* Configure the swarm */
		swarm = Swarm_New( "testSwarm", NULL, elementCellLayout, particleLayout, nDims, sizeof(Particle),
		extensionMgr_Register, NULL, data->comm, NULL );
	
		/* Build the swarm */
		Stg_Component_Build( swarm, 0, False );
		Stg_Component_Initialise( swarm, 0, False );

		pcu_check_true( swarm->cellDomainCount == 30 );
		pcu_check_true( swarm->cellLocalCount == 30 );
		pcu_check_true( swarm->particleLocalCount == 10 );
		pcu_check_true( swarm->particles[0].owningCell );
		pcu_check_true( swarm->particles[1].owningCell );
		pcu_check_true( swarm->particles[3].owningCell );
		pcu_check_true( swarm->particles[4].owningCell );
		pcu_check_true( swarm->particles[7].owningCell );
		pcu_check_true( swarm->particles[8].owningCell );
		pcu_check_true( swarm->particles[2].owningCell );
		pcu_check_true( swarm->particles[5].owningCell );
		pcu_check_true( swarm->particles[6].owningCell );
		pcu_check_true( swarm->particles[9].owningCell );

		Stg_Class_Delete( extensionMgr_Register );
		Stg_Component_Destroy( elementCellLayout, NULL, True );
		Stg_Component_Destroy( particleLayout, NULL, True );
		/*Stg_Component_Destroy( mesh, NULL, True );*/
		Stg_Component_Destroy( swarm, NULL, True );	
	}
}

void WithinShapeParticleLayoutSuite( pcu_suite_t* suite ) {
	pcu_suite_setData( suite, WithinShapeParticleLayoutSuiteData );
	pcu_suite_setFixtures( suite, WithinShapeParticleLayoutSuite_Setup, WithinShapeParticleLayoutSuite_Teardown );
	pcu_suite_addTest( suite, WithinShapeParticleLayoutSuite_TestWithinShapeSphere );
}


