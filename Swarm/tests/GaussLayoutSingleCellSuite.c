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
**   Tests the GaussLayoutSingleCellSuite
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

#include "GaussLayoutSingleCellSuite.h"

struct _Particle {
	__IntegrationPoint
};

typedef struct {
	unsigned							nDims;
	unsigned							gaussParticles[3];
	double							minCrds[3];
	double							maxCrds[3];
	ExtensionManager_Register*	extensionMgr_Register;
	GaussParticleLayout*			gaussParticleLayout;
	SingleCellLayout*				singleCellLayout;
	Swarm*							swarm;
	MPI_Comm							comm;
	unsigned int					rank;
	unsigned int					nProcs;
} GaussLayoutSingleCellSuiteData;

void GaussLayoutSingleCellSuite_Setup( GaussLayoutSingleCellSuiteData* data ) {
	/* MPI Initializations */
	data->comm = MPI_COMM_WORLD;  
	MPI_Comm_rank( data->comm, &data->rank );
	MPI_Comm_size( data->comm, &data->nProcs );
   
	data->nDims = 3;
	data->gaussParticles[0] = 2; data->gaussParticles[1] = 1; data->gaussParticles[2] = 3;
	data->minCrds[0] = 0.0; data->minCrds[1] = 0.0; data->minCrds[2] = 0.0;
	data->maxCrds[0] = 300.0; data->maxCrds[1] = 12.0; data->maxCrds[2] = 300.0;
	
	data->extensionMgr_Register = ExtensionManager_Register_New();
}

void GaussLayoutSingleCellSuite_Teardown( GaussLayoutSingleCellSuiteData* data ) {
	/* Destroy stuff */
	Stg_Class_Delete( data->gaussParticleLayout );
	Stg_Class_Delete( data->singleCellLayout );
	Stg_Class_Delete( data->swarm );
	Stg_Class_Delete( data->extensionMgr_Register );
	/* Clean generated output files */
	remove( "1ParticlePerDim_3D.dat" );
	remove( "2ParticlesPerDim_3D.dat" );
	remove( "3ParticlesPerDim_3D.dat" );
}

void GaussLayoutSingleCellSuite_Test1ParticlePerDim_3D( GaussLayoutSingleCellSuiteData* data ) {
	Cell_PointIndex		count;
	double					x,y,z,w;
	unsigned int			p, i, len;
	LocalParticle*			particle;
	Coord						minCell;
	Coord						maxCell;
	Stream*					stream;
	Particle_InCellIndex	particlesPerDim[3] = {1, 1, 1};
	Bool						dimExists[] = { True, True, True };	
	char						expected_file[PCU_PATH_MAX];
	
	if( data->rank == 0 ) {	
		/* Configure the element-cell-layout */
		data->singleCellLayout = SingleCellLayout_New( "singleCellLayout", dimExists, NULL, NULL );
	
		/* Configure the gauss-particle-layout */
		data->gaussParticleLayout = GaussParticleLayout_New( "gaussParticleLayout", data->nDims , particlesPerDim );
	
		/* Configure the swarm */
		data->swarm = Swarm_New( "testGaussSwarmSingleCell", data->singleCellLayout, data->gaussParticleLayout, data->nDims,
			sizeof(Particle), data->extensionMgr_Register, NULL, data->comm, NULL );
		
		/* Build the swarm */
		Stg_Component_Build( data->swarm, 0, False );
		Stg_Component_Initialise( data->swarm, 0, False );

		count = data->swarm->cellParticleCountTbl[0];
	 
		stream = Journal_Register( Info_Type, "1ParticlePerDim_3D" );
		Stream_RedirectFile( stream, "1ParticlePerDim_3D.dat" );
	
		for( p = 0; p < count; p++ ) {
			x = ((IntegrationPoint*)Swarm_ParticleInCellAt( data->swarm, 0, p ))->xi[0]; 
			y = ((IntegrationPoint*)Swarm_ParticleInCellAt( data->swarm, 0, p ))->xi[1]; 
			z = ((IntegrationPoint*)Swarm_ParticleInCellAt( data->swarm, 0, p ))->xi[2]; 	
			w = ((IntegrationPoint*)Swarm_ParticleInCellAt( data->swarm, 0, p ))->weight;
			Journal_Printf( stream, "pId=%d : xi = { %f, %f, %f } weight = %f\n",p,x,y,z,w );	
		}	
		pcu_filename_expected( "testGaussLayoutSingleCell1ParticlePerDimOutput.expected", expected_file );
		pcu_check_fileEq( "1ParticlePerDim_3D.dat", expected_file );
	}
}


void GaussLayoutSingleCellSuite_Test2ParticlesPerDim_3D( GaussLayoutSingleCellSuiteData* data ) {
	Cell_PointIndex		count;
	double					x,y,z,w;
	unsigned int			p, i, len;
	LocalParticle*			particle;
	Coord						minCell;
	Coord						maxCell;
	Stream*					stream;
	Particle_InCellIndex	particlesPerDim[3] = {2, 2, 2};
	Bool						dimExists[] = { True, True, True };	
	char 						expected_file[PCU_PATH_MAX];
	
	if( data->rank == 0 ) {	
		/* Configure the element-cell-layout */
		data->singleCellLayout = SingleCellLayout_New( "singleCellLayout", dimExists, NULL, NULL );
	
		/* Configure the gauss-particle-layout */
		data->gaussParticleLayout = GaussParticleLayout_New( "gaussParticleLayout", data->nDims , particlesPerDim );
		
		/* Configure the swarm */
		data->swarm = Swarm_New( "testGaussSwarmSingleCell", data->singleCellLayout, data->gaussParticleLayout, data->nDims,
			sizeof(Particle), data->extensionMgr_Register, NULL, data->comm, NULL );
		
		/* Build the swarm */
		Stg_Component_Build( data->swarm, 0, False );
		Stg_Component_Initialise( data->swarm, 0, False );

		count = data->swarm->cellParticleCountTbl[0];
		 
		stream = Journal_Register( Info_Type, "2ParticlesPerDim_3D" );
		Stream_RedirectFile( stream, "2ParticlesPerDim_3D.dat" );
	
		for( p = 0; p < count; p++ ) {
			x = ((IntegrationPoint*)Swarm_ParticleInCellAt( data->swarm, 0, p ))->xi[0]; 
			y = ((IntegrationPoint*)Swarm_ParticleInCellAt( data->swarm, 0, p ))->xi[1]; 
			z = ((IntegrationPoint*)Swarm_ParticleInCellAt( data->swarm, 0, p ))->xi[2]; 	
			w = ((IntegrationPoint*)Swarm_ParticleInCellAt( data->swarm, 0, p ))->weight;
			Journal_Printf( stream, "pId=%d : xi = { %f, %f, %f } weight = %f\n",p,x,y,z,w );	
		}	
		pcu_filename_expected( "testGaussLayoutSingleCell2ParticlesPerDimOutput.expected", expected_file );
		pcu_check_fileEq( "2ParticlesPerDim_3D.dat", expected_file );
	}
}


void GaussLayoutSingleCellSuite_Test3ParticlesPerDim_3D( GaussLayoutSingleCellSuiteData* data ) {
	Cell_PointIndex		count;
	double					x,y,z,w;
	unsigned int			p, i, len;
	LocalParticle*			particle;
	Coord						minCell;
	Coord						maxCell;
	Stream*					stream;
	Particle_InCellIndex	particlesPerDim[3] = {3, 3, 3};
	Bool						dimExists[] = { True, True, True };	
	char						expected_file[PCU_PATH_MAX];

	if( data->rank == 0 ) {	
		/* Configure the element-cell-layout */
		data->singleCellLayout = SingleCellLayout_New( "singleCellLayout", dimExists, NULL, NULL );
	
		/* Configure the gauss-particle-layout */
		data->gaussParticleLayout = GaussParticleLayout_New( "gaussParticleLayout", data->nDims , particlesPerDim );
	
		/* Configure the swarm */
		data->swarm = Swarm_New( "testGaussSwarmSingleCell", data->singleCellLayout, data->gaussParticleLayout, data->nDims,
			sizeof(Particle), data->extensionMgr_Register, NULL, data->comm, NULL );
		
		/* Build the swarm */
		Stg_Component_Build( data->swarm, 0, False );
		Stg_Component_Initialise( data->swarm, 0, False );
	
		count = data->swarm->cellParticleCountTbl[0];
	 
		stream = Journal_Register( Info_Type, "3ParticlesPerDim_3D" );
		Stream_RedirectFile( stream, "3ParticlesPerDim_3D.dat" );
	
		for( p = 0; p < count; p++ ) {
			x = ((IntegrationPoint*)Swarm_ParticleInCellAt( data->swarm, 0, p ))->xi[0]; 
			y = ((IntegrationPoint*)Swarm_ParticleInCellAt( data->swarm, 0, p ))->xi[1]; 
			z = ((IntegrationPoint*)Swarm_ParticleInCellAt( data->swarm, 0, p ))->xi[2]; 	
			w = ((IntegrationPoint*)Swarm_ParticleInCellAt( data->swarm, 0, p ))->weight;
			Journal_Printf( stream, "pId=%d : xi = { %f, %f, %f } weight = %f\n",p,x,y,z,w );	
		}	
		pcu_filename_expected( "testGaussLayoutSingleCell3ParticlesPerDimOutput.expected", expected_file );
		pcu_check_fileEq( "3ParticlesPerDim_3D.dat", expected_file );
	}
}


void GaussLayoutSingleCellSuite( pcu_suite_t* suite ) {
	pcu_suite_setData( suite, GaussLayoutSingleCellSuiteData );
	pcu_suite_setFixtures( suite, GaussLayoutSingleCellSuite_Setup, GaussLayoutSingleCellSuite_Teardown );
	pcu_suite_addTest( suite, GaussLayoutSingleCellSuite_Test1ParticlePerDim_3D );
	pcu_suite_addTest( suite, GaussLayoutSingleCellSuite_Test2ParticlesPerDim_3D );
	pcu_suite_addTest( suite, GaussLayoutSingleCellSuite_Test3ParticlesPerDim_3D );
}

