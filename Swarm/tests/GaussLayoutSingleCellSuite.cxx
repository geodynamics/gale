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
	MPI_Comm comm;
	int rank;
	int nProcs;
} GaussLayoutSingleCellSuiteData;

void GaussLayoutSingleCellSuite_Setup( GaussLayoutSingleCellSuiteData* data ) {
	/* MPI Initializations */
	data->comm = MPI_COMM_WORLD;  
	MPI_Comm_rank( data->comm, &data->rank );
	MPI_Comm_size( data->comm, &data->nProcs );
}

void GaussLayoutSingleCellSuite_Teardown( GaussLayoutSingleCellSuiteData* data ) {
}

void GaussLayoutSingleCellSuite_Test1ParticlePerDim_3D( GaussLayoutSingleCellSuiteData* data ) {
	unsigned							nDims;
	ExtensionManager_Register*	extensionMgr_Register;
	GaussParticleLayout*			gaussParticleLayout;
	SingleCellLayout*				singleCellLayout;
	Swarm*							swarm;
	int								procToWatch = data->nProcs > 1 ? 1 : 0;
	Cell_PointIndex				count;
	double							x,y,z,w;
	unsigned int					p;
	Stream*							stream;
	Particle_InCellIndex			particlesPerDim[3] = {1, 1, 1};
	Bool								dimExists[] = { True, True, True };	
	char								expected_file[PCU_PATH_MAX];
	
	if( data->rank == procToWatch ) {	
		stream = Journal_Register( Info_Type, (Name)"1ParticlePerDim_3D"  );
		Stream_RedirectFile( stream, "1ParticlePerDim_3D.dat" );

		nDims = 3;
	
		extensionMgr_Register = ExtensionManager_Register_New();

		/* Configure the element-cell-layout */
		singleCellLayout = SingleCellLayout_New( "singleCellLayout", NULL, dimExists, NULL, NULL );
	
		/* Configure the gauss-particle-layout */
		gaussParticleLayout = GaussParticleLayout_New( "gaussParticleLayout", NULL,
           LocalCoordSystem, True, nDims, particlesPerDim );

		swarm = Swarm_New( "testGaussSwarmSingleCell", NULL, singleCellLayout, gaussParticleLayout, nDims,
			sizeof(Particle), extensionMgr_Register, NULL, data->comm, NULL );
		
		/* Build the swarm */
		Stg_Component_Build( swarm, 0, False );
		Stg_Component_Initialise( swarm, 0, False );

		count = swarm->cellParticleCountTbl[0];
	 
		for( p = 0; p < count; p++ ) {
			x = ((IntegrationPoint*)Swarm_ParticleInCellAt( swarm, 0, p ))->xi[0]; 
			y = ((IntegrationPoint*)Swarm_ParticleInCellAt( swarm, 0, p ))->xi[1]; 
			z = ((IntegrationPoint*)Swarm_ParticleInCellAt( swarm, 0, p ))->xi[2]; 	
			w = ((IntegrationPoint*)Swarm_ParticleInCellAt( swarm, 0, p ))->weight;
			Journal_Printf( stream, "pId=%d : xi = { %f, %f, %f } weight = %f\n",p,x,y,z,w );	
		}	
		pcu_filename_expected( "testGaussLayoutSingleCell1ParticlePerDimOutput.expected", expected_file );
		pcu_check_fileEq( "1ParticlePerDim_3D.dat", expected_file );
		remove( "1ParticlePerDim_3D.dat" );

		Stg_Class_Delete( extensionMgr_Register );
		Stg_Component_Destroy( singleCellLayout, NULL, True );
		Stg_Component_Destroy( gaussParticleLayout, NULL, True );
		Stg_Component_Destroy( swarm, NULL, True );
	}
}


void GaussLayoutSingleCellSuite_Test2ParticlesPerDim_3D( GaussLayoutSingleCellSuiteData* data ) {
	unsigned							nDims;
	ExtensionManager_Register*	extensionMgr_Register;
	GaussParticleLayout*			gaussParticleLayout;
	SingleCellLayout*				singleCellLayout;
	Swarm*							swarm;
	int								procToWatch = data->nProcs > 1 ? 1 : 0;
	Cell_PointIndex				count;
	double							x,y,z,w;
	unsigned int					p;
	Stream*							stream;
	Particle_InCellIndex			particlesPerDim[3] = {2, 2, 2};
	Bool								dimExists[] = { True, True, True };	
	char								expected_file[PCU_PATH_MAX];
	
	if( data->rank == procToWatch ) {	
		stream = Journal_Register( Info_Type, (Name)"2ParticlesPerDim_3D"  );
		Stream_RedirectFile( stream, "2ParticlesPerDim_3D.dat" );

		nDims = 3;
		extensionMgr_Register = ExtensionManager_Register_New();

		/* Configure the element-cell-layout */
		singleCellLayout = SingleCellLayout_New( "singleCellLayout", NULL, dimExists, NULL, NULL );
	
		/* Configure the gauss-particle-layout */
		gaussParticleLayout = GaussParticleLayout_New( "gaussParticleLayout", NULL,
           LocalCoordSystem, True, nDims, particlesPerDim );
		
		/* Configure the swarm */
		swarm = Swarm_New( "testGaussSwarmSingleCell", NULL, singleCellLayout, gaussParticleLayout, nDims,
			sizeof(Particle), extensionMgr_Register, NULL, data->comm, NULL );
		
		/* Build the swarm */
		Stg_Component_Build( swarm, 0, False );
		Stg_Component_Initialise( swarm, 0, False );

		count = swarm->cellParticleCountTbl[0];
		 
		for( p = 0; p < count; p++ ) {
			x = ((IntegrationPoint*)Swarm_ParticleInCellAt( swarm, 0, p ))->xi[0]; 
			y = ((IntegrationPoint*)Swarm_ParticleInCellAt( swarm, 0, p ))->xi[1]; 
			z = ((IntegrationPoint*)Swarm_ParticleInCellAt( swarm, 0, p ))->xi[2]; 	
			w = ((IntegrationPoint*)Swarm_ParticleInCellAt( swarm, 0, p ))->weight;
			Journal_Printf( stream, "pId=%d : xi = { %f, %f, %f } weight = %f\n",p,x,y,z,w );	
		}	
		pcu_filename_expected( "testGaussLayoutSingleCell2ParticlesPerDimOutput.expected", expected_file );
		pcu_check_fileEq( "2ParticlesPerDim_3D.dat", expected_file );
		remove( "2ParticlesPerDim_3D.dat" );

		Stg_Class_Delete( extensionMgr_Register );
		Stg_Component_Destroy( singleCellLayout, NULL, True );
		Stg_Component_Destroy( gaussParticleLayout, NULL, True );
		Stg_Component_Destroy( swarm, NULL, True );
	}
}

void GaussLayoutSingleCellSuite_Test3ParticlesPerDim_3D( GaussLayoutSingleCellSuiteData* data ) {
	unsigned							nDims;
	ExtensionManager_Register*	extensionMgr_Register;
	GaussParticleLayout*			gaussParticleLayout;
	SingleCellLayout*				singleCellLayout;
	Swarm*							swarm;
	int								procToWatch = data->nProcs > 1 ? 1 : 0;
	Cell_PointIndex				count;
	double							x,y,z,w;
	unsigned int					p;
	Stream*							stream;
	Particle_InCellIndex			particlesPerDim[3] = {3, 3, 3};
	Bool								dimExists[] = { True, True, True };	
	char								expected_file[PCU_PATH_MAX];

	if( data->rank == procToWatch ) {	
		stream = Journal_Register( Info_Type, (Name)"3ParticlesPerDim_3D"  );
		Stream_RedirectFile( stream, "3ParticlesPerDim_3D.dat" );

		nDims = 3;
		extensionMgr_Register = ExtensionManager_Register_New();

		/* Configure the element-cell-layout */
		singleCellLayout = SingleCellLayout_New( "singleCellLayout", NULL, dimExists, NULL, NULL );
	
		/* Configure the gauss-particle-layout */
		gaussParticleLayout = GaussParticleLayout_New( "gaussParticleLayout", NULL,
           LocalCoordSystem, True, nDims, particlesPerDim );
	
		/* Configure the swarm */
		swarm = Swarm_New( "testGaussSwarmSingleCell", NULL, singleCellLayout, gaussParticleLayout, nDims,
			sizeof(Particle), extensionMgr_Register, NULL, data->comm, NULL );
		
		/* Build the swarm */
		Stg_Component_Build( swarm, 0, False );
		Stg_Component_Initialise( swarm, 0, False );
	
		count = swarm->cellParticleCountTbl[0];
	 
		for( p = 0; p < count; p++ ) {
			x = ((IntegrationPoint*)Swarm_ParticleInCellAt( swarm, 0, p ))->xi[0]; 
			y = ((IntegrationPoint*)Swarm_ParticleInCellAt( swarm, 0, p ))->xi[1]; 
			z = ((IntegrationPoint*)Swarm_ParticleInCellAt( swarm, 0, p ))->xi[2]; 	
			w = ((IntegrationPoint*)Swarm_ParticleInCellAt( swarm, 0, p ))->weight;
			Journal_Printf( stream, "pId=%d : xi = { %f, %f, %f } weight = %f\n",p,x,y,z,w );	
		}	
		pcu_filename_expected( "testGaussLayoutSingleCell3ParticlesPerDimOutput.expected", expected_file );
		pcu_check_fileEq( "3ParticlesPerDim_3D.dat", expected_file );
		remove( "3ParticlesPerDim_3D.dat" );

		Stg_Class_Delete( extensionMgr_Register );
		Stg_Component_Destroy( singleCellLayout, NULL, True );
		Stg_Component_Destroy( gaussParticleLayout, NULL, True );
		Stg_Component_Destroy( swarm, NULL, True );
	}
}


void GaussLayoutSingleCellSuite( pcu_suite_t* suite ) {
	pcu_suite_setData( suite, GaussLayoutSingleCellSuiteData );
	pcu_suite_setFixtures( suite, GaussLayoutSingleCellSuite_Setup, GaussLayoutSingleCellSuite_Teardown );
	pcu_suite_addTest( suite, GaussLayoutSingleCellSuite_Test1ParticlePerDim_3D );
	pcu_suite_addTest( suite, GaussLayoutSingleCellSuite_Test2ParticlesPerDim_3D );
	pcu_suite_addTest( suite, GaussLayoutSingleCellSuite_Test3ParticlesPerDim_3D );
}



