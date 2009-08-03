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
**   Tests the ParticleCoordsSuite
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

#include "ParticleCoordsSuite.h"

struct _Particle {
	__IntegrationPoint
};

typedef struct {
	ExtensionManager_Register*	extensionMgr_Register;
	Variable_Register*		variable_Register;
	Swarm*				swarm;
	MPI_Comm			comm;
	unsigned int			rank;
	unsigned int			nProcs;
} ParticleCoordsSuiteData;

void ParticleCoordsSuite_Setup( ParticleCoordsSuiteData* data ) {
	char input_file[PCU_PATH_MAX];

	/* MPI Initializations */
	data->comm = MPI_COMM_WORLD;  
	MPI_Comm_rank( data->comm, &data->rank );
	MPI_Comm_size( data->comm, &data->nProcs );

	data->extensionMgr_Register = ExtensionManager_Register_New();   
	data->variable_Register = Variable_Register_New();
}

void ParticleCoordsSuite_Teardown( ParticleCoordsSuiteData* data ) {
	/* Destroy stuff */
	Stg_Class_Delete( data->swarm );
	Stg_Class_Delete( data->extensionMgr_Register );
	Stg_Class_Delete( data->variable_Register );
}

void ParticleCoordsSuite_TestLineParticle( ParticleCoordsSuiteData* data ) {
	int		procToWatch;
	Stream*		stream;
	Dictionary*	dictionary;
	DomainContext*	context;
	char		input_file[PCU_PATH_MAX];
	char		expected_file[PCU_PATH_MAX];
	
	if( data->nProcs >= 2 ) {
		procToWatch = 1;
	}
	else {
		procToWatch = 0;
	}
	if( data->rank == procToWatch ) printf( "Watching rank: %i\n", data->rank );

	/* read in the xml input file */
	pcu_filename_input( "testLineParticleLayout.xml", input_file );
	context = (DomainContext*)stgMainInitFromXML( input_file, data->comm );
	dictionary = context->dictionary;

	stream = Journal_Register( Info_Type, "ParticleCoords" );
	procToWatch = Dictionary_GetUnsignedInt( dictionary, "procToWatch" );
	
	if( data->rank == procToWatch ) {
		Stream_RedirectFile_WithPrependedPath( stream, Dictionary_GetString( dictionary, "outputPath" ), "output.dat" );

		data->swarm = (Swarm*) LiveComponentRegister_Get( context->CF->LCRegister, "swarm" );
		pcu_check_true( data->swarm );

		Swarm_PrintParticleCoords( data->swarm, stream );

		pcu_filename_expected( "testLineParticleLayout.0of1.output.dat.expected", expected_file );
		pcu_check_fileEq( "output/output.dat", expected_file );
	}
}

void ParticleCoordsSuite( pcu_suite_t* suite ) {
	pcu_suite_setData( suite, ParticleCoordsSuiteData );
	pcu_suite_setFixtures( suite, ParticleCoordsSuite_Setup, ParticleCoordsSuite_Teardown );
	pcu_suite_addTest( suite, ParticleCoordsSuite_TestLineParticle );
}
