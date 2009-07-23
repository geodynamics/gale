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
	unsigned			nDims;
	unsigned			meshSize[3];
	double				minCrds[3];
	double				maxCrds[3];
	ExtensionManager_Register*	extensionMgr_Register;
	Mesh*				mesh;
	GaussParticleLayout*		gaussParticleLayout;
	ElementCellLayout*		elementCellLayout;
	Swarm*				swarm;
	//Dimension_Index     		dim;
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
}

void GaussLayoutSuite_1ParticlePerDim_3D( GaussLayoutSuiteData* data ) {
	Cell_PointIndex		count;
	double 					x,y,z;
	unsigned int 			p, i, len;
	LocalParticle* 		particle;
	Coord 			minCell;
	Coord 			maxCell;
	Particle_InCellIndex 	particlesPerDim[3] = {1, 1, 1};
	unsigned		dim = 3;

	data->gaussParticleLayout = GaussParticleLayout_New( "gaussParticleLayout", dim, particlesPerDim );
	/* Configure the swarm */
	data->swarm = Swarm_New( "testGaussSwarm", data->elementCellLayout, data->gaussParticleLayout, 
				dim, sizeof(Particle), data->extensionMgr_Register, NULL, MPI_COMM_WORLD, NULL );
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
	
	for( p = 0; p < count; p++ ) {
		particle = (LocalParticle*)Swarm_ParticleInCellAt( data->swarm, 4, p  );

		/* convert to global coords */
		x = 0.5 * ( maxCell[0] - minCell[0] ) * ( particle->xi[0] + 1.0 ) + minCell[0];
		y = 0.5 * ( maxCell[1] - minCell[1] ) * ( particle->xi[1] + 1.0 ) + minCell[1];
		z = 0.5 * ( maxCell[2] - minCell[2] ) * ( particle->xi[2] + 1.0 ) + minCell[2];
		
		printf( "pId=%d : coords = { %.12f, %.12f, %.12f }, xi = { %.12f, %.12f, %.12f }\n", 
			p, x, y, z, particle->xi[0], particle->xi[1], particle->xi[2] );

		pcu_check_lt( ( particle->xi[0] - 0.0 ), 1.0E-6);
		pcu_check_lt( ( particle->xi[1] - 0.0 ), 1.0E-6);
		pcu_check_lt( ( particle->xi[2] - 0.0 ), 1.0E-6);
				
		pcu_check_lt( ( x - 75.0 ), 1.0E-6 );
		pcu_check_lt( ( y - 10.0 ), 1.0E-6 );
		pcu_check_lt( ( z - 75.0 ), 1.0E-6 );
	}
}

void GaussLayoutSuite_2ParticlePerDim_3D( GaussLayoutSuiteData* data ) {
	Cell_PointIndex		count;
	double 					x,y,z;
	unsigned int			p, i, len;
	LocalParticle* 		particle;
	Coord 					minCell;
	Coord 					maxCell;
	Particle_InCellIndex 	particlesPerDim[3] = {2, 2, 2};
	unsigned					dim = 3;
	double					xi[8][3];
	double					coord[8][3];

	data->gaussParticleLayout = GaussParticleLayout_New( "gaussParticleLayout", dim, particlesPerDim );
	/* Configure the swarm */
	data->swarm = Swarm_New( "testGaussSwarm", data->elementCellLayout, data->gaussParticleLayout, 
				dim, sizeof(Particle), data->extensionMgr_Register, NULL, MPI_COMM_WORLD, NULL );
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
	
	xi[0][0] = -0.577350269190;	xi[0][1] = -0.577350269190;	xi[0][2] = -0.577350269190;
	xi[1][0] = +0.577350269190;	xi[1][1] = -0.577350269190;	xi[1][2] = -0.577350269190;
	xi[2][0] = -0.577350269190;	xi[2][1] = +0.577350269190;	xi[2][2] = -0.577350269190;
	xi[3][0] = +0.577350269190;	xi[3][1] = +0.577350269190;	xi[3][2] = -0.577350269190;
	xi[4][0] = -0.577350269190;	xi[4][1] = -0.577350269190;	xi[4][2] = +0.577350269190;
	xi[5][0] = +0.577350269190;	xi[5][1] = -0.577350269190;	xi[5][2] = +0.577350269190;
	xi[6][0] = -0.577350269190;	xi[6][1] = +0.577350269190;	xi[6][2] = +0.577350269190;
	xi[7][0] = +0.577350269190;	xi[7][1] = +0.577350269190;	xi[7][2] = +0.577350269190;

	coord[0][0] = 31.698729810778;	coord[0][1] = 8.845299461621;		coord[0][2] = 31.698729810778;
	coord[1][0] = 118.301270189222;	coord[1][1] = 8.845299461621;		coord[1][2] = 31.698729810778;
	coord[2][0] = 31.698729810778;	coord[2][1] = 11.154700538379;	coord[2][2] = 31.698729810778;
	coord[3][0] = 118.301270189222;	coord[3][1] = 11.154700538379;	coord[3][2] = 31.698729810778;
	coord[4][0] = 31.698729810778;	coord[4][1] = 8.845299461621;		coord[4][2] = 118.301270189222;
	coord[5][0] = 118.301270189222;	coord[5][1] = 8.845299461621;		coord[5][2] = 118.301270189222;
	coord[6][0] = 31.698729810778;	coord[6][1] = 11.154700538379;	coord[6][2] = 118.301270189222;
	coord[7][0] = 118.301270189222;	coord[7][1] = 11.154700538379;	coord[7][2] = 118.301270189222;

	for( p = 0; p < count; p++ ) {
		particle = (LocalParticle*)Swarm_ParticleInCellAt( data->swarm, 4, p  );

		/* convert to global coords */
		x = 0.5 * ( maxCell[0] - minCell[0] ) * ( particle->xi[0] + 1.0 ) + minCell[0];
		y = 0.5 * ( maxCell[1] - minCell[1] ) * ( particle->xi[1] + 1.0 ) + minCell[1];
		z = 0.5 * ( maxCell[2] - minCell[2] ) * ( particle->xi[2] + 1.0 ) + minCell[2];
		
		printf( "pId=%d : coords = { %.12f, %.12f, %.12f }, xi = { %.12f, %.12f, %.12f }\n", 
			p, x, y, z, particle->xi[0], particle->xi[1], particle->xi[2] );
				
		pcu_check_lt( ( particle->xi[0] - xi[p][0] ), 1.0E-6 );
		pcu_check_lt( ( particle->xi[1] - xi[p][1] ), 1.0E-6 );
		pcu_check_lt( ( particle->xi[2] - xi[p][2] ), 1.0E-6 );
		
		pcu_check_lt( ( x - coord[p][0] ), 1.0E-6 );
		pcu_check_lt( ( y - coord[p][1] ), 1.0E-6 );
		pcu_check_lt( ( z - coord[p][2] ), 1.0E-6 );
	}
}


void GaussLayoutSuite_3ParticlePerDim_3D( GaussLayoutSuiteData* data ) {
	Cell_PointIndex		count;
	double		 			x,y,z;
	unsigned int			p, i, len;
	LocalParticle* 		particle;
	Coord 					minCell;
	Coord 					maxCell;
	Particle_InCellIndex 	particlesPerDim[3] = {3, 3, 3};
	unsigned					dim = 3;
	double					xi[27][3];
	double					coord[27][3];

	data->gaussParticleLayout = GaussParticleLayout_New( "gaussParticleLayout", dim, particlesPerDim );
	/* Configure the swarm */
	data->swarm = Swarm_New( "testGaussSwarm", data->elementCellLayout, data->gaussParticleLayout, 
				dim, sizeof(Particle), data->extensionMgr_Register, NULL, MPI_COMM_WORLD, NULL );
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
	
	xi[0][0] = -0.774596669241;	xi[0][1] = -0.774596669241;	xi[0][2] = -0.774596669241;
	xi[1][0] = 0.000000000000;		xi[1][1] = -0.774596669241;	xi[1][2] = -0.774596669241;
	xi[2][0] = 0.774596669241;		xi[2][1] = -0.774596669241;	xi[2][2] = -0.774596669241;
	xi[3][0] = -0.774596669241;	xi[3][1] = 0.000000000000;		xi[3][2] = -0.774596669241;
	xi[4][0] = 0.000000000000;		xi[4][1] = 0.000000000000;		xi[4][2] = -0.774596669241;
	xi[5][0] = 0.774596669241;		xi[5][1] = 0.000000000000;		xi[5][2] = -0.774596669241;
	xi[6][0] = -0.774596669241;	xi[6][1] = 0.774596669241;		xi[6][2] = -0.774596669241;
	xi[7][0] = 0.000000000000;		xi[7][1] = 0.774596669241;		xi[7][2] = -0.774596669241;
	xi[8][0] = 0.774596669241;		xi[8][1] = 0.774596669241;		xi[8][2] = -0.774596669241;
	xi[9][0] = -0.774596669241;	xi[9][1] = -0.774596669241;	xi[9][2] = 0.000000000000;
	xi[10][0] = 0.000000000000;	xi[10][1] = -0.774596669241;	xi[10][2] = 0.000000000000;
	xi[11][0] = 0.774596669241;	xi[11][1] = -0.774596669241;	xi[11][2] = 0.000000000000;
	xi[12][0] = -0.774596669241;	xi[12][1] = 0.000000000000;	xi[12][2] = 0.000000000000;
	xi[13][0] = 0.000000000000;	xi[13][1] = 0.000000000000;	xi[13][2] = 0.000000000000;
	xi[14][0] = 0.774596669241;	xi[14][1] = 0.000000000000;	xi[14][2] = 0.000000000000;
	xi[15][0] = -0.774596669241;	xi[15][1] = 0.774596669241;	xi[15][2] = 0.000000000000;
	xi[16][0] = 0.000000000000;	xi[16][1] = 0.774596669241;	xi[16][2] = 0.000000000000;
	xi[17][0] = 0.774596669241;	xi[17][1] = 0.774596669241;	xi[17][2] = 0.000000000000;
	xi[18][0] = -0.774596669241;	xi[18][1] = -0.774596669241;	xi[18][2] = 0.774596669241;
	xi[19][0] = 0.000000000000;	xi[19][1] = -0.774596669241;	xi[19][2] = 0.774596669241;
	xi[20][0] = 0.774596669241;	xi[20][1] = -0.774596669241;	xi[20][2] = 0.774596669241;
	xi[21][0] = -0.774596669241;	xi[21][1] = 0.000000000000;	xi[21][2] = 0.774596669241;
	xi[22][0] = 0.000000000000;	xi[22][1] = 0.000000000000;	xi[22][2] = 0.774596669241;
	xi[23][0] = 0.774596669241;	xi[23][1] = 0.000000000000;	xi[23][2] = 0.774596669241;
	xi[24][0] = -0.774596669241;	xi[24][1] = 0.774596669241;	xi[24][2] = 0.774596669241;
	xi[25][0] = 0.000000000000;	xi[25][1] = 0.774596669241;	xi[25][2] = 0.774596669241;
	xi[26][0] = 0.774596669241;	xi[26][1] = 0.774596669241;	xi[26][2] = 0.774596669241;

	coord[0][0] = 16.905249806889;	coord[0][1] = 8.450806661517;		coord[0][2] = 16.905249806889;
	coord[1][0] = 75.000000000000;	coord[1][1] = 8.450806661517;		coord[1][2] = 16.905249806889;
	coord[2][0] = 133.094750193111;	coord[2][1] = 8.450806661517;		coord[2][2] = 16.905249806889;
	coord[3][0] = 16.905249806889;	coord[3][1] = 10.000000000000;	coord[3][2] = 16.905249806889;
	coord[4][0] = 75.000000000000;	coord[4][1] = 10.000000000000;	coord[4][2] = 16.905249806889;
	coord[5][0] = 133.094750193111;	coord[5][1] = 10.000000000000;	coord[5][2] = 16.905249806889;
	coord[6][0] = 16.905249806889;	coord[6][1] = 11.549193338483;	coord[6][2] = 16.905249806889;
	coord[7][0] = 75.000000000000;	coord[7][1] = 11.549193338483;	coord[7][2] = 16.905249806889;
	coord[8][0] = 133.094750193111;	coord[8][1] = 11.549193338483;	coord[8][2] = 16.905249806889;
	coord[9][0] = 16.905249806889;	coord[9][1] = 8.450806661517;		coord[9][2] = 75.000000000000;
	coord[10][0] = 75.000000000000;	coord[10][1] = 8.450806661517;	coord[10][2] = 75.000000000000;
	coord[11][0] = 133.094750193111;	coord[11][1] = 8.450806661517;	coord[11][2] = 75.000000000000;
	coord[12][0] = 16.905249806889;	coord[12][1] = 10.000000000000;	coord[12][2] = 75.000000000000;
	coord[13][0] = 75.000000000000;	coord[13][1] = 10.000000000000;	coord[13][2] = 75.000000000000;
	coord[14][0] = 133.094750193111;	coord[14][1] = 10.000000000000;	coord[14][2] = 75.000000000000;
	coord[15][0] = 16.905249806889;	coord[15][1] = 11.549193338483;	coord[15][2] = 75.000000000000;
	coord[16][0] = 75.000000000000;	coord[16][1] = 11.549193338483;	coord[16][2] = 75.000000000000;
	coord[17][0] = 133.094750193111;	coord[17][1] = 11.549193338483;	coord[17][2] = 75.000000000000;
	coord[18][0] = 16.905249806889;	coord[18][1] = 8.450806661517;	coord[18][2] = 133.094750193111;
	coord[19][0] = 75.000000000000;	coord[19][1] = 8.450806661517;	coord[19][2] = 133.094750193111;
	coord[20][0] = 133.094750193111;	coord[20][1] = 8.845299461621;	coord[20][2] = 133.094750193111;
	coord[21][0] = 16.905249806889;	coord[21][1] = 10.000000000000;	coord[21][2] = 133.094750193111;
	coord[22][0] = 75.000000000000;	coord[22][1] = 10.000000000000;	coord[22][2] = 133.094750193111;
	coord[23][0] = 133.094750193111;	coord[23][1] = 10.000000000000;	coord[23][2] = 133.094750193111;
	coord[24][0] = 16.905249806889;	coord[24][1] = 11.549193338483;	coord[24][2] = 133.094750193111;
	coord[25][0] = 75.000000000000;	coord[25][1] = 11.549193338483;	coord[25][2] = 133.094750193111;
	coord[26][0] = 133.094750193111;	coord[26][1] = 11.549193338483;	coord[26][2] = 133.094750193111;

	for( p = 0; p < count; p++ ) {
		particle = (LocalParticle*)Swarm_ParticleInCellAt( data->swarm, 4, p  );

		/* convert to global coords */
		x = 0.5 * ( maxCell[0] - minCell[0] ) * ( particle->xi[0] + 1.0 ) + minCell[0];
		y = 0.5 * ( maxCell[1] - minCell[1] ) * ( particle->xi[1] + 1.0 ) + minCell[1];
		z = 0.5 * ( maxCell[2] - minCell[2] ) * ( particle->xi[2] + 1.0 ) + minCell[2];
		
		printf( "pId=%d : coords = { %.12f, %.12f, %.12f }, xi = { %.12f, %.12f, %.12f }\n", 
			p, x, y, z, particle->xi[0], particle->xi[1], particle->xi[2] );
		
		pcu_check_lt( ( particle->xi[0] - xi[p][0] ), 1.0E-6 );
		pcu_check_lt( ( particle->xi[1] - xi[p][1] ), 1.0E-6 );
		pcu_check_lt( ( particle->xi[2] - xi[p][2] ), 1.0E-6 );
		
		pcu_check_lt( ( x - coord[p][0] ), 1.0E-6);
		pcu_check_lt( ( y - coord[p][1] ), 1.0E-6);
		pcu_check_lt( ( z - coord[p][2] ), 1.0E-6); 
	}
}

void GaussLayoutSuite( pcu_suite_t* suite ) {
   pcu_suite_setData( suite, GaussLayoutSuiteData );
   pcu_suite_setFixtures( suite, GaussLayoutSuite_Setup, GaussLayoutSuite_Teardown );
   pcu_suite_addTest( suite, GaussLayoutSuite_1ParticlePerDim_3D );
   pcu_suite_addTest( suite, GaussLayoutSuite_2ParticlePerDim_3D );
   pcu_suite_addTest( suite, GaussLayoutSuite_3ParticlePerDim_3D );
}
