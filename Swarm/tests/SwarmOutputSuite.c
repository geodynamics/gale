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

#include "SwarmOutputSuite.h"

#define CURR_MODULE_NAME "SwarmOutputSuite"

typedef struct {
	MPI_Comm comm;
	int		rank;
	int		nProcs;
} SwarmOutputSuiteData;

double SwarmOutputSuite_Dt( void* context ) {
	return 2.0;
}

void _SwarmOutputSuite_SetDt( void* context, double dt ) {
}

void SwarmOutputSuite_MoveParticles( DomainContext* context ) {
	Swarm*				swarm = (Swarm*) LiveComponentRegister_Get( context->CF->LCRegister, "swarm" );
	Particle_Index		lParticle_I;
	GlobalParticle*	particle;
	double				x,y;

	for ( lParticle_I = 0 ; lParticle_I < swarm->particleLocalCount; lParticle_I++ ) {
		particle = (GlobalParticle*)Swarm_ParticleAt( swarm, lParticle_I );

		x = particle->coord[ I_AXIS ];
		y = particle->coord[ J_AXIS ];
		particle->coord[ I_AXIS ] = 1.0 - y;
		particle->coord[ J_AXIS ] = x;
	}

	Swarm_UpdateAllParticleOwners( swarm );
}

void SwarmOutputSuite_Setup( SwarmOutputSuiteData* data ) {
	Journal_Enable_AllTypedStream( False );

	/* MPI Initializations */
	data->comm = MPI_COMM_WORLD;  
	MPI_Comm_rank( data->comm, &data->rank );
	MPI_Comm_size( data->comm, &data->nProcs );
}

void SwarmOutputSuite_Teardown( SwarmOutputSuiteData* data ) {
	Journal_Enable_AllTypedStream( True );
}

void SwarmOutputSuite_TestSwarmOutput( SwarmOutputSuiteData* data ) {
	int procToWatch = data->nProcs > 1 ? 1 : 0;

	if( data->rank == procToWatch ) {
		Dictionary*					dictionary;
		Dictionary*					componentDict;
		Stg_ComponentFactory*	cf;
		DomainContext*				context;
		char							input_file[PCU_PATH_MAX];
		Swarm*						swarm;
		SwarmOutput*				swarmOutput;
		SpaceFillerParticleLayout* particleLayout;

		pcu_filename_input( "testSwarmOutput.xml", input_file );
		cf = stgMainInitFromXML( input_file, data->comm, NULL );
		context = (DomainContext*) LiveComponentRegister_Get( cf->LCRegister, "context" );
		dictionary = context->dictionary;

		Journal_ReadFromDictionary( dictionary );
		stgMainBuildAndInitialise( cf );

		ContextEP_Append( context, AbstractContext_EP_Dt, SwarmOutputSuite_Dt );
		ContextEP_Append( context, AbstractContext_EP_Step, SwarmOutputSuite_MoveParticles );

		componentDict = Dictionary_GetDictionary( dictionary, "components" );
		assert( componentDict );

		AbstractContext_Dump( context );
		Stg_Component_Execute( context, 0, False );

		particleLayout = (SpaceFillerParticleLayout*) LiveComponentRegister_Get( context->CF->LCRegister, "particleLayout" );
		swarmOutput = (SwarmOutput*) LiveComponentRegister_Get( context->CF->LCRegister, "swarmOutput" );
		swarm = (Swarm*) LiveComponentRegister_Get( context->CF->LCRegister, "swarm" );

		pcu_check_true( particleLayout->isConstructed && particleLayout->isBuilt && particleLayout->isInitialised );
		pcu_check_true( swarmOutput->isConstructed && swarmOutput->isBuilt && swarmOutput->isInitialised );
		pcu_check_true( swarm->isConstructed && swarm->isBuilt && swarm->isInitialised );

		pcu_check_streq( swarm->name, swarmOutput->swarm->name );

		_Stg_Component_Delete( context );
		stgMainDestroy( cf );
	}
}

void SwarmOutputSuite( pcu_suite_t* suite ) {
	pcu_suite_setData( suite, SwarmOutputSuiteData );
	pcu_suite_setFixtures( suite, SwarmOutputSuite_Setup, SwarmOutputSuite_Teardown );
	pcu_suite_addTest( suite, SwarmOutputSuite_TestSwarmOutput );
}


