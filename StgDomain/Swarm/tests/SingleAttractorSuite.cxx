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
**   Tests the SingleAttractorSuite
**
** $Id: testSwarm.c 3462 2006-02-19 06:53:24Z WalterLandry $
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

#include "SingleAttractorSuite.h"

/* silly stgermain, I must define this */
#define CURR_MODULE_NAME "DomainContext.c"

struct _Particle {
	__GlobalParticle
	double velocity[3];
	double randomColour;
};

void SingleAttractorSuite_SaveSwarms( void* context );
void SingleAttractorSuite_SingleAttractor( DomainContext* context );
double SingleAttractorSuite_Dt( void* context ) {
	return 2.0;
}

typedef struct {
	MPI_Comm	comm;
	int		rank;
	int		nProcs;
} SingleAttractorSuiteData;

void SingleAttractorSuite_Setup( SingleAttractorSuiteData* data ) {
	Journal_Enable_AllTypedStream( False );

	data->comm = MPI_COMM_WORLD;  
	MPI_Comm_rank( data->comm, &data->rank );
	MPI_Comm_size( data->comm, &data->nProcs );
}

void SingleAttractorSuite_Teardown( SingleAttractorSuiteData* data ) {
	Journal_Enable_AllTypedStream( True );
}

void SingleAttractorSuite_TestSingleAttractor( SingleAttractorSuiteData* data ) {
	Dictionary*					dictionary;
	Stg_ComponentFactory*	cf;
	Swarm*						swarm = NULL;
	Particle						particle;
	Particle*					currParticle = NULL;
	Particle_Index				lParticle_I = 0;
	Dimension_Index			dim_I = 0;
	DomainContext*				context;
	char							input_file[PCU_PATH_MAX];
	char							expected_file[PCU_PATH_MAX];

	pcu_filename_input( "testSwarmParticleAdvectionSingleAttractor.xml", input_file );

	cf = stgMainInitFromXML( input_file, data->comm, NULL );
	context = (DomainContext*)LiveComponentRegister_Get( cf->LCRegister, (Name)"context"  );
	Stream_Enable( cf->infoStream, False );
	Stream_Enable( context->info, False );
	Stream_Enable( context->verbose, False );
	Stream_Enable( context->debug, False );

	dictionary = context->dictionary;
	Journal_ReadFromDictionary( dictionary );
	assert( Dictionary_GetDictionary( dictionary, "components" ) );

	KeyCall( context, context->constructExtensionsK, EntryPoint_VoidPtr_CallCast* )( KeyHandle(context,context->constructExtensionsK), context );

	swarm = (Swarm*) LiveComponentRegister_Get( context->CF->LCRegister, (Name)"swarm"  );
	ExtensionManager_Add( swarm->particleExtensionMgr, (Name)"ParticleVelocity", sizeof(double[3])  );
	ExtensionManager_Add( swarm->particleExtensionMgr, (Name)"ParticleColour", sizeof(double)  );

	Swarm_NewVectorVariable(
		swarm,
		"Velocity",
		(ArithPointer) &particle.velocity - (ArithPointer) &particle,
		Variable_DataType_Double,
		swarm->dim,
		"VelocityX",
		"VelocityY",
		"VelocityZ" );

	Swarm_NewScalarVariable(
		swarm,
		"RandomColour",
		(ArithPointer) &particle.randomColour - (ArithPointer) &particle,
		Variable_DataType_Double );

	stgMainBuildAndInitialise( cf );

	/* for each particle, set a random colour */
	for ( lParticle_I=0; lParticle_I < swarm->particleLocalCount; lParticle_I++ ) {
		currParticle = (Particle*)Swarm_ParticleAt( swarm, lParticle_I );
		for ( dim_I=0; dim_I < 3; dim_I++ ) {
			currParticle->velocity[dim_I] = 0;
		}	
		currParticle->randomColour = ( (double)  rand() ) / RAND_MAX;
	}
	
	ContextEP_ReplaceAll( context, AbstractContext_EP_Dt, SingleAttractorSuite_Dt );
	ContextEP_Append( context, AbstractContext_EP_Save, SingleAttractorSuite_SaveSwarms );
	ContextEP_ReplaceAll( context, AbstractContext_EP_Solve, SingleAttractorSuite_SingleAttractor );

	Stg_Component_Execute( context, 0, False );
	Stg_Component_Destroy( context, 0, False );

	pcu_filename_expected( "testSwarmParticleAdvectionSingleAttractor.expected", expected_file );
	pcu_check_fileEq( "testSingleAttractor.dat", expected_file );
	remove( "testSingleAttractor.dat" );

   stgMainDestroy( cf );
}

void SingleAttractorSuite( pcu_suite_t* suite ) {
	pcu_suite_setData( suite, SingleAttractorSuiteData );
	pcu_suite_setFixtures( suite, SingleAttractorSuite_Setup, SingleAttractorSuite_Teardown );
	pcu_suite_addTest( suite, SingleAttractorSuite_TestSingleAttractor );
}

void SingleAttractorSuite_SaveSwarms( void* context ) {
   Swarm_Register_SaveAllRegisteredSwarms( Swarm_Register_GetSwarm_Register(), context );
}

void SingleAttractorSuite_SingleAttractor( DomainContext* context ) {
	Cell_LocalIndex		lCell_I;
	Particle_InCellIndex	cParticle_I;
	Particle*				currParticle;
	Index						dim_I;
	Swarm*					swarm = (Swarm*) LiveComponentRegister_Get( context->CF->LCRegister, (Name)"swarm"  );
	Coord						attractorPoint;
	Mesh*						mesh;
	Stream*					stream = Journal_Register( Info_Type, (Name)"particleUpdate"  );
	unsigned int			movementSpeedDivisor = 0;
	int						movementSign = 1;
	unsigned int			explosionPeriod = 20;
	double					minCrd[3], maxCrd[3];

	Stream_RedirectFile( stream, "testSingleAttractor.dat" );

	Stream_SetPrintingRank( stream, Dictionary_GetUnsignedInt_WithDefault( context->dictionary, "procToWatch", 0 ) );
	movementSpeedDivisor = Dictionary_GetDouble_WithDefault( context->dictionary, (Dictionary_Entry_Key)"movementSpeedDivisor", 10 );
	
	mesh = (Mesh* )LiveComponentRegister_Get( context->CF->LCRegister, (Name)"mesh-linear"  );
	Mesh_GetGlobalCoordRange( mesh, minCrd, maxCrd );
	for ( dim_I=0; dim_I < 3; dim_I++ ) {
		attractorPoint[dim_I] = (maxCrd[dim_I] - minCrd[dim_I]) / 3;
	}
	Journal_Printf( stream, "Calculated attractor point is at (%f,%f,%f):\n", attractorPoint[0], attractorPoint[1], attractorPoint[2] );
	
	/* Now decide if we are attracting or repelling */
	if ( ( ( context->timeStep / explosionPeriod ) % 2 ) == 0 ) {
		Journal_Printf( stream, "Timestep %d - Implosive mode\n", context->timeStep );
		movementSign = 1;
	}
	else {
		Journal_Printf( stream, "Timestep %d - Explosive mode\n", context->timeStep );
		movementSign = -1;
	}	

	for ( lCell_I=0; lCell_I < swarm->cellLocalCount; lCell_I++ ) {
		Journal_Printf( stream, "\tUpdating Particles positions in local cell %d:\n", lCell_I );
		for ( cParticle_I=0; cParticle_I < swarm->cellParticleCountTbl[lCell_I]; cParticle_I++ ) {
			Coord movementVector = {0,0,0};
			Coord newParticleCoord = {0,0,0};
			Coord* oldCoord;

			currParticle = (Particle*)Swarm_ParticleInCellAt( swarm, lCell_I, cParticle_I );
			oldCoord = &currParticle->coord;
			Journal_Printf( stream, "\t\tUpdating particleInCell %d:\n", cParticle_I );

			for ( dim_I=0; dim_I < 3; dim_I++ ) {
				movementVector[dim_I] = ( attractorPoint[dim_I] - (*oldCoord)[dim_I] ) / movementSpeedDivisor;
				movementVector[dim_I] *= movementSign;	
				if ( movementSign == -1 ) {
					movementVector[dim_I] *= (float)movementSpeedDivisor / (movementSpeedDivisor-1); 
				}
				newParticleCoord[dim_I] = (*oldCoord)[dim_I] + movementVector[dim_I];
			}
			memcpy( currParticle->velocity, movementVector, 3*sizeof(double) ); 

			Journal_Printf( stream, "\t\tChanging its coords from (%f,%f,%f) to (%f,%f,%f):\n",
				(*oldCoord)[0], (*oldCoord)[1], (*oldCoord)[2],
				newParticleCoord[0], newParticleCoord[1], newParticleCoord[2] );

			for ( dim_I=0; dim_I < 3; dim_I++ ) {
				currParticle->coord[dim_I] = newParticleCoord[dim_I];
			}
		}
	}
	Swarm_UpdateAllParticleOwners( swarm );
}


