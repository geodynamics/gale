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
**   Tests the AdvectionSuite
**
** $Id: testAdvection.c 3462 2006-02-19 06:53:24Z WalterLandry $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>

#include "pcu/pcu.h"
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include "AdvectionSuite.h"

#define CURR_MODULE_NAME "AdvectionSuite"

ExtensionInfo_Index handle;

typedef struct {
} AdvectionSuiteData;


void AdvectionSuite_Setup( AdvectionSuiteData* data ) {
}

void AdvectionSuite_Teardown( AdvectionSuiteData* data ) {
}

double dt( void* class, PICelleratorContext* context ) {
	return Dictionary_GetDouble_WithDefault( context->dictionary, "dt", 0.01 );
}

void check( PICelleratorContext* context ) {
	MaterialPointsSwarm*   materialPointsSwarm;
	GlobalParticle*   particle;
	double*           originalCoord;
	double*           coord;
	Particle_Index    lParticle_I;
	double            maxRadiusError       = 0.0;
	double            maxThetaError        = 0.0;
	double            maxDepthError        = 0.0;
	double            maxRadiusErrorGlobal;
	double            maxThetaErrorGlobal;
	double            maxDepthErrorGlobal;
	double            currentRadius;
	double            originalRadius;
	double            originalTheta;
	Coord             analyticCoord;
	double            time                 = context->currentTime + context->dt;
	Dictionary*       dictionary           = context->dictionary;
	Stream*           stream               = Journal_Register( Info_Type, CURR_MODULE_NAME );
	double		  depthErrorTolerance  = Dictionary_GetDouble( dictionary, "depthErrorTolerance" );
	double		  radiusErrorTolerance = Dictionary_GetDouble( dictionary, "radiusErrorTolerance" );
	double		  thetaErrorTolerance  = Dictionary_GetDouble( dictionary, "thetaErrorTolerance" );
	double		  theta;

	/* Add original pos to particle */
	materialPointsSwarm = (MaterialPointsSwarm*) LiveComponentRegister_Get( context->CF->LCRegister, "materialPointsSwarm" );

	for ( lParticle_I = 0 ; lParticle_I < materialPointsSwarm->particleLocalCount ; lParticle_I++ ) {
		particle      = (GlobalParticle*)Swarm_ParticleAt( materialPointsSwarm, lParticle_I );
		coord         = particle->coord;
		originalCoord = ExtensionManager_Get( materialPointsSwarm->particleExtensionMgr, particle, handle );

		currentRadius  = StGermain_VectorMagnitude( coord, 2 );
		originalRadius = StGermain_VectorMagnitude( originalCoord, 2 );
		
		originalTheta = acos( originalCoord[ I_AXIS ]/originalRadius );
		if ( originalCoord[ J_AXIS ] < 0.0 )
			originalTheta = 2 * M_PI - originalTheta;

		analyticCoord[ I_AXIS ] = originalRadius*cos( time + originalTheta );
		analyticCoord[ J_AXIS ] = originalRadius*sin( time + originalTheta );
		analyticCoord[ K_AXIS ] = originalCoord[ K_AXIS ];

		maxDepthError  = MAX( maxDepthError,  fabs( originalCoord[ K_AXIS ] - coord[ K_AXIS ] ) );
		maxRadiusError = MAX( maxRadiusError, fabs( currentRadius - originalRadius ) );

		theta = StGermain_AngleBetweenVectors( analyticCoord, coord, 2 );
		if( theta > maxThetaError )
			maxThetaError = theta;
		//maxThetaError  = MAX( maxThetaError, theta );
	}

	MPI_Allreduce( &maxDepthError,  &maxDepthErrorGlobal,  1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
	MPI_Allreduce( &maxRadiusError, &maxRadiusErrorGlobal, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
	MPI_Allreduce( &maxThetaError,  &maxThetaErrorGlobal,  1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );

	pcu_check_true( maxDepthErrorGlobal < depthErrorTolerance );
	pcu_check_true( maxRadiusErrorGlobal < radiusErrorTolerance );
	pcu_check_true( maxThetaErrorGlobal < thetaErrorTolerance );
}

void AdvectionSuite_Euler( AdvectionSuiteData* data ) {
	Stg_ComponentFactory* 	cf;
	PICelleratorContext* 	context;
	char			filename[PCU_PATH_MAX];
	MaterialPointsSwarm*	matSwarm;
	unsigned		particle_i;
	GlobalParticle*		particle;
	double*			origCoord;
	
	pcu_filename_input( "testEuler.xml", filename );
	context = (PICelleratorContext*)_PICelleratorContext_DefaultNew( "context" );
	cf = stgMainInitFromXML( filename, MPI_COMM_WORLD, context );

	matSwarm = (MaterialPointsSwarm*)LiveComponentRegister_Get( context->CF->LCRegister, "materialPointsSwarm" );
	handle = ExtensionManager_Add( matSwarm->particleExtensionMgr, CURR_MODULE_NAME, sizeof(Coord) );

	Stg_ComponentFactory_BuildComponents( cf, NULL );
	Stg_ComponentFactory_InitialiseComponents( cf, NULL );

	for( particle_i = 0; particle_i < matSwarm->particleLocalCount; particle_i++ ) {
		particle = (GlobalParticle*)Swarm_ParticleAt( matSwarm, particle_i );
		origCoord = ExtensionManager_Get( matSwarm->particleExtensionMgr, particle, handle );
		memcpy( origCoord, particle->coord, sizeof(Coord) );
	}

	ContextEP_Prepend( context, AbstractContext_EP_Step, check );
	EP_AppendClassHook( Context_GetEntryPoint( context, FiniteElementContext_EP_CalcDt ), dt, context );

	stgMainLoop( cf );

	Stg_ComponentFactory_DestroyComponents( cf, data );
}

void AdvectionSuite_RK2( AdvectionSuiteData* data ) {
	Stg_ComponentFactory* 	cf;
	PICelleratorContext* 	context;
	char			filename[PCU_PATH_MAX];
	MaterialPointsSwarm*	matSwarm;
	unsigned		particle_i;
	GlobalParticle*		particle;
	double*			origCoord;
	
	pcu_filename_input( "testRK2.xml", filename );
	context = (PICelleratorContext*)_PICelleratorContext_DefaultNew( "context" );
	cf = stgMainInitFromXML( filename, MPI_COMM_WORLD, context );

	matSwarm = (MaterialPointsSwarm*)LiveComponentRegister_Get( context->CF->LCRegister, "materialPointsSwarm" );
	handle = ExtensionManager_Add( matSwarm->particleExtensionMgr, CURR_MODULE_NAME, sizeof(Coord) );

	Stg_ComponentFactory_BuildComponents( cf, NULL );
	Stg_ComponentFactory_InitialiseComponents( cf, NULL );

	for( particle_i = 0; particle_i < matSwarm->particleLocalCount; particle_i++ ) {
		particle = (GlobalParticle*)Swarm_ParticleAt( matSwarm, particle_i );
		origCoord = ExtensionManager_Get( matSwarm->particleExtensionMgr, particle, handle );
		memcpy( origCoord, particle->coord, sizeof(Coord) );
	}

	ContextEP_Prepend( context, AbstractContext_EP_Step, check );
	EP_AppendClassHook( Context_GetEntryPoint( context, FiniteElementContext_EP_CalcDt ), dt, context );

	stgMainLoop( cf );

	Stg_ComponentFactory_DestroyComponents( cf, data );
}

void AdvectionSuite_RK4( AdvectionSuiteData* data ) {
	Stg_ComponentFactory* 	cf;
	PICelleratorContext* 	context;
	char			filename[PCU_PATH_MAX];
	MaterialPointsSwarm*	matSwarm;
	unsigned		particle_i;
	GlobalParticle*		particle;
	double*			origCoord;
	
	pcu_filename_input( "testRK4.xml", filename );
	context = (PICelleratorContext*)_PICelleratorContext_DefaultNew( "context" );
	cf = stgMainInitFromXML( filename, MPI_COMM_WORLD, context );

	matSwarm = (MaterialPointsSwarm*)LiveComponentRegister_Get( context->CF->LCRegister, "materialPointsSwarm" );
	handle = ExtensionManager_Add( matSwarm->particleExtensionMgr, CURR_MODULE_NAME, sizeof(Coord) );

	Stg_ComponentFactory_BuildComponents( cf, NULL );
	Stg_ComponentFactory_InitialiseComponents( cf, NULL );

	for( particle_i = 0; particle_i < matSwarm->particleLocalCount; particle_i++ ) {
		particle = (GlobalParticle*)Swarm_ParticleAt( matSwarm, particle_i );
		origCoord = ExtensionManager_Get( matSwarm->particleExtensionMgr, particle, handle );
		memcpy( origCoord, particle->coord, sizeof(Coord) );
	}

	ContextEP_Prepend( context, AbstractContext_EP_Step, check );
	EP_AppendClassHook( Context_GetEntryPoint( context, FiniteElementContext_EP_CalcDt ), dt, context );

	stgMainLoop( cf );

	Stg_ComponentFactory_DestroyComponents( cf, data );
}

void AdvectionSuite( pcu_suite_t* suite ) {
   pcu_suite_setData( suite, AdvectionSuiteData );
   pcu_suite_setFixtures( suite, AdvectionSuite_Setup, AdvectionSuite_Teardown );
   pcu_suite_addTest( suite, AdvectionSuite_Euler );
   pcu_suite_addTest( suite, AdvectionSuite_RK2 );
   pcu_suite_addTest( suite, AdvectionSuite_RK4 );
}
