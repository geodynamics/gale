/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003-2006, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street,
**	Melbourne, 3053, Australia.
** Copyright (c) 2005-2006, Monash Cluster Computing, Building 28, Monash University Clayton Campus,
**	Victoria, 3800, Australia
**
** Primary Contributing Organisations:
**	Victorian Partnership for Advanced Computing Ltd, Computational Software Development - http://csd.vpac.org
**	Australian Computational Earth Systems Simulator - http://www.access.edu.au
**	Monash Cluster Computing - http://www.mcc.monash.edu.au
**
** Contributors:
**	Robert Turnbull, Research Assistant, Monash University. (robert.turnbull@sci.monash.edu.au)
**	Patrick D. Sunter, Software Engineer, VPAC. (patrick@vpac.org)
**	Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
**	Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**	David May, PhD Student, Monash University (david.may@sci.monash.edu.au)
**	Vincent Lemiale, Postdoctoral Fellow, Monash University. (vincent.lemiale@sci.monash.edu.au)
**	Julian Giordani, Research Assistant, Monash University. (julian.giordani@sci.monash.edu.au)
**	Louis Moresi, Associate Professor, Monash University. (louis.moresi@sci.monash.edu.au)
**	Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**	Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
**	David Stegman, Postdoctoral Fellow, Monash University. (david.stegman@sci.monash.edu.au)
**	Wendy Sharples, PhD Student, Monash University (wendy.sharples@sci.monash.edu.au)
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
** $Id: testAdvection.c 456 2007-04-27 06:21:01Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/Voronoi/Voronoi.h>
#include <PICellerator/PopulationControl/PopulationControl.h>
#include <PICellerator/Weights/Weights.h>
#include <PICellerator/MaterialPoints/MaterialPoints.h>

#include <math.h>
#include <string.h>
#include <assert.h>

ExtensionInfo_Index handle;

const Type TestAdvection_Type = "TestAdvection";

double dt( void* class, PICelleratorContext* context ) {
	return Dictionary_GetDouble_WithDefault( context->dictionary, "dt", 0.01 );
}

void construct( PICelleratorContext* context ) {
	MaterialPointsSwarm* materialPointsSwarm;

	/* Add original pos to particle */
	materialPointsSwarm = (MaterialPointsSwarm*) LiveComponentRegister_Get( context->CF->LCRegister, "materialPointsSwarm" );
	handle = ExtensionManager_Add( materialPointsSwarm->particleExtensionMgr, CURR_MODULE_NAME, sizeof( Coord ) );
}

void storeOriginalPos( PICelleratorContext* context ) {
	MaterialPointsSwarm*   materialPointsSwarm;
	GlobalParticle*   particle;
	double*           originalCoord;
	Particle_Index    lParticle_I;

	materialPointsSwarm = (MaterialPointsSwarm*) LiveComponentRegister_Get( context->CF->LCRegister, "materialPointsSwarm" );

	for ( lParticle_I = 0 ; lParticle_I < materialPointsSwarm->particleLocalCount ; lParticle_I++ ) {
		particle = (GlobalParticle*)Swarm_ParticleAt( materialPointsSwarm, lParticle_I );
		originalCoord = ExtensionManager_Get( materialPointsSwarm->particleExtensionMgr, particle, handle );

		memcpy( originalCoord, particle->coord, sizeof(Coord) );
	}
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

	/* Add original pos to particle */
	materialPointsSwarm = (MaterialPointsSwarm*) LiveComponentRegister_Get( context->CF->LCRegister, "materialPointsSwarm" );

	for ( lParticle_I = 0 ; lParticle_I < materialPointsSwarm->particleLocalCount ; lParticle_I++ ) {
		particle      = (GlobalParticle*)Swarm_ParticleAt( materialPointsSwarm, lParticle_I );
		coord         = particle->coord;
		originalCoord = ExtensionManager_Get( materialPointsSwarm->particleExtensionMgr, particle, handle );

		currentRadius  = StGermain_VectorMagnitude( coord, 2 );
		originalRadius = StGermain_VectorMagnitude( originalCoord, 2 );
		
		/* if ( originalRadius >= 1.0 || currentRadius >= 1.0 ) */
		/* 	continue; */

		originalTheta = acos( originalCoord[ I_AXIS ]/originalRadius );
		if ( originalCoord[ J_AXIS ] < 0.0 )
			originalTheta = 2 * M_PI - originalTheta;

		analyticCoord[ I_AXIS ] = originalRadius*cos( time + originalTheta );
		analyticCoord[ J_AXIS ] = originalRadius*sin( time + originalTheta );
		analyticCoord[ K_AXIS ] = originalCoord[ K_AXIS ];

		maxDepthError  = MAX( maxDepthError,  fabs( originalCoord[ K_AXIS ] - coord[ K_AXIS ] ) );
		maxRadiusError = MAX( maxRadiusError, fabs( currentRadius - originalRadius ) );
		maxThetaError  = MAX( maxThetaError,  StGermain_AngleBetweenVectors( analyticCoord, coord, 2 ) );
	}

	MPI_Allreduce( &maxDepthError,  &maxDepthErrorGlobal,  1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
	MPI_Allreduce( &maxRadiusError, &maxRadiusErrorGlobal, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );
	MPI_Allreduce( &maxThetaError,  &maxThetaErrorGlobal,  1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );

	Journal_Printf( stream, "Timestep = %u\n", context->timeStep );
	Journal_Printf( stream, "%s Depth Test\n",
			maxDepthErrorGlobal < Dictionary_GetDouble( dictionary, "depthErrorTolerance" ) ? "Passed" : "Failed" );
	Journal_Printf( stream, "%s Radius Test\n",
			maxRadiusErrorGlobal < Dictionary_GetDouble( dictionary, "radiusErrorTolerance" ) ? "Passed" : "Failed" );
	Journal_Printf( stream, "%s Theta Test\n",
			maxThetaErrorGlobal < Dictionary_GetDouble( dictionary, "thetaErrorTolerance" ) ? "Passed" : "Failed" );
}


void _testAdvection_Construct( void* component, Stg_ComponentFactory* cf, void* data ) {
	DiscretisationContext* context;
	Stream*           stream               = Journal_Register( Info_Type, CURR_MODULE_NAME );

	context = (DiscretisationContext*)Stg_ComponentFactory_ConstructByName( cf, "context", DiscretisationContext, True, data ); 

	ContextEP_Append( context, AbstractContext_EP_ConstructExtensions, construct );
	ContextEP_Append( context, AbstractContext_EP_Initialise, storeOriginalPos );
	ContextEP_Prepend( context, AbstractContext_EP_Step, check );
	EP_AppendClassHook( Context_GetEntryPoint( context, FiniteElementContext_EP_CalcDt ), dt, context );

	Stream_RedirectFile_WithPrependedPath( stream, context->outputPath, "output.dat" );
	Stream_SetPrintingRank( stream, 0 );
}

void* _testAdvection_DefaultNew( Name name ) {
	return _Codelet_New(
			sizeof( Codelet ),
			TestAdvection_Type,
			_Codelet_Delete,
			_Codelet_Print,
			_Codelet_Copy,
			_testAdvection_DefaultNew,
			_testAdvection_Construct,
			_Codelet_Build,
			_Codelet_Initialise,
			_Codelet_Execute,
			_Codelet_Destroy,
			name );
}


Index testAdvection_Register( PluginsManager* pluginsManager ) {
	Index result;

	result = PluginsManager_Submit( pluginsManager, TestAdvection_Type, "0",
		_testAdvection_DefaultNew );

	return result;
}
