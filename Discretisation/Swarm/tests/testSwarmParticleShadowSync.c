/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street, Melbourne, 3053, Australia.
**
** Authors:
**	Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**	Patrick D. Sunter, Software Engineer, VPAC. (pds@vpac.org)
**	Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**	Siew-Ching Tan, Software Engineer, VPAC. (siew@vpac.org)
**	Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
**	Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
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
**	Tests that particles can be successfully moved between cells. The problem is set up with a 
**	"gravitational attractor" in the exact middle of the domain - all particles are sucked in
**	towards it each timestep.
**
** Assumptions:
**	None as yet.
**
** Comments:
**	None as yet.
**
** $Id: testSwarmParticleAdvection.c 4081 2007-04-27 06:20:07Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include "Base/Base.h"

#include "Discretisation/Geometry/Geometry.h"
#include "Discretisation/Shape/Shape.h"
#include "Discretisation/Mesh/Mesh.h"
#include "Discretisation/Utils/Utils.h"
#include "Discretisation/Swarm/Swarm.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>


struct _Node {
	Coord				coord;
};

struct _Element {
	Coord				coord;
};

struct _Particle {
	__GlobalParticle
	double  velocity[3];
	double  randomColour;
};

double Dt( void* context ) {
	return 2.0;
}

void ValidateShadowing( void *dummy, void* context );

/** Global so other funcs can use */
Index procToWatch = 0;

int main( int argc, char* argv[] ) {
	DiscretisationContext*          context;
	MPI_Comm			CommWorld;
	int				rank;
	int				numProcessors;
	Dictionary*			dictionary;
	Dictionary*                     componentDict;
	Stg_ComponentFactory*           cf;
	XML_IO_Handler*                 ioHandler;
	ExtensionManager_Register*      extensionMgr_Register;
	SwarmVariable_Register*         swarmVariable_Register;
	Stream*                         stream;
	Swarm*                          swarm = NULL;
	Particle                        particle;
	Particle*                       currParticle = NULL;
	Particle_Index                  lParticle_I = 0;
	Dimension_Index                 dim_I = 0;
	char*                           directory;
	
	/* Initialise MPI, get world info */
	MPI_Init( &argc, &argv );
	MPI_Comm_dup( MPI_COMM_WORLD, &CommWorld );
	MPI_Comm_size( CommWorld, &numProcessors );
	MPI_Comm_rank( CommWorld, &rank );
	
	Base_Init( &argc, &argv );
	
	DiscretisationGeometry_Init( &argc, &argv );
	DiscretisationShape_Init( &argc, &argv );
	DiscretisationMesh_Init( &argc, &argv );
	DiscretisationUtils_Init( &argc, &argv );
	DiscretisationSwarm_Init( &argc, &argv );
	MPI_Barrier( CommWorld ); /* Ensures copyright info always come first in output */

	/* Add the StGermain path to the global xml path dictionary */
	directory = Memory_Alloc_Array( char, 200, "xmlDirectory" ) ;
	sprintf(directory, "%s%s", LIB_DIR, "/StGermain" );
	XML_IO_Handler_AddDirectory( "StGermain", directory  );
	Memory_Free(directory);
	/* Add the plugin path to the global plugin list */
	PluginsManager_AddDirectory( "StGermain", LIB_DIR );

	stream = Journal_Register (Info_Type, "myStream");

	dictionary = Dictionary_New();
	ioHandler = XML_IO_Handler_New();
	IO_Handler_ReadAllFromCommandLine( ioHandler, argc, argv, dictionary );

	/* TODO: temporary hack until Al gets the journal read from file going again */
	if ( False == Dictionary_GetBool_WithDefault( dictionary, "particleCommInfo", True ) ) {
		Stream_Enable( Journal_Register( Info_Type, ParticleCommHandler_Type ), False );
	}

	Journal_ReadFromDictionary( dictionary );

	/* *** Journal stuff *** */
	Journal_Enable_TypedStream( DebugStream_Type, False );
	Stream_EnableBranch( Swarm_Debug, True );
	Stream_SetLevelBranch( Swarm_Debug, 3 );

	if( argc >= 2 ) {
		procToWatch = atoi( argv[1] );
	}
	else {
		procToWatch = 0;
	}
	if( rank == procToWatch ) printf( "Watching rank: %i\n", rank );
	Stream_SetPrintingRank( Journal_Register( InfoStream_Type, "Context" ), procToWatch);
	/* For plugins to read */
	Dictionary_Add( dictionary, "procToWatch", Dictionary_Entry_Value_FromUnsignedInt( procToWatch ) );
	
/* Construction phase -------------------------------------------------------------------------------------------*/
	
	/* Create the Context */
	context = DiscretisationContext_New(
			"context",
			0,
			0,
			MPI_COMM_WORLD,
			dictionary );

	componentDict = Dictionary_GetDictionary( dictionary, "components" );
	
	assert( componentDict );

	cf = context->CF = Stg_ComponentFactory_New( dictionary, componentDict, context->register_Register );

	LiveComponentRegister_Add( cf->LCRegister, (Stg_Component*) context );
	PluginsManager_Load( context->plugins, context, dictionary );

	extensionMgr_Register = ExtensionManager_Register_New();
	swarmVariable_Register = SwarmVariable_Register_New( NULL );
	Stg_ObjectList_ClassAppend( cf->registerRegister, (void*)extensionMgr_Register, "ExtensionManager_Register" );
	Stg_ObjectList_ClassAppend( cf->registerRegister, (void*)swarmVariable_Register, "SwarmVariable_Register" );

	Stg_ComponentFactory_CreateComponents( cf );
	Stg_ComponentFactory_ConstructComponents( cf, 0 /* dummy */ );
	PluginsManager_ConstructPlugins( context->plugins, context->CF, 0 /* dummy */ );

	KeyCall( context, context->constructExtensionsK, EntryPoint_VoidPtr_CallCast* )( KeyHandle(context,context->constructExtensionsK), context );

	swarm = (Swarm*) LiveComponentRegister_Get( context->CF->LCRegister, "swarm" );
	ExtensionManager_Add( swarm->particleExtensionMgr, "ParticleVelocity", sizeof(double[3]) );
	ExtensionManager_Add( swarm->particleExtensionMgr, "ParticleColour", sizeof(double) );

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

	LiveComponentRegister_BuildAll( cf->LCRegister, context );
	LiveComponentRegister_InitialiseAll( cf->LCRegister, context );

	/* for each particle, set a random colour */
	for ( lParticle_I=0; lParticle_I < swarm->particleLocalCount; lParticle_I++ ) {
		currParticle = (Particle*)Swarm_ParticleAt( swarm, lParticle_I );
		for ( dim_I=0; dim_I < 3; dim_I++ ) {
			currParticle->velocity[dim_I] = 0;
		}	
		currParticle->randomColour = ( (double)  rand() ) / RAND_MAX;
	}
	
	if( rank == procToWatch ) {
          /* Print( swarm, stream ); */
	}	

	Stg_Component_Build( context, 0 /* dummy */, False );
	Stg_Component_Initialise( context, 0 /* dummy */, False );
	
	/* +++ RUN PHASE +++ */
	AbstractContext_Dump( context );

	ContextEP_ReplaceAll( context, AbstractContext_EP_Dt, Dt );

	ContextEP_Append( context, AbstractContext_EP_Sync, ValidateShadowing );
	Stg_Component_Execute( context, 0 /* dummy */, False );
	Stg_Component_Destroy( context, 0 /* dummy */, False );

	/* Delete stuff */
	/* Deleting the component factory automatically deletes all components in it */
	/* TODO: should the component factory be renamed a comp. manager? Since it deletes */
	/* 	components as well? */

	if( rank == procToWatch ) {
		Journal_Printf( Journal_Register( Info_Type, "success" ), "Shadow particle validation: passed\n" );
	}

	Stg_Class_Delete( cf );
	/* Remaining registers etc that don't live on the context or anything */
	Stg_Class_Delete( extensionMgr_Register );
	Stg_Class_Delete( swarmVariable_Register );
	/* Input/Output stuff */
	Stg_Class_Delete( dictionary );
	Stg_Class_Delete( ioHandler );
	
	DiscretisationSwarm_Finalise();
	DiscretisationUtils_Finalise();
	DiscretisationMesh_Finalise();
	DiscretisationShape_Finalise();
	DiscretisationGeometry_Finalise();
	
	Base_Finalise();
	
	/* Close off MPI */
	MPI_Finalize();

	return 0; /* success */
}

int listCompareFunction( void *a, void *b ){
	if( a>b )
		return 1;
	else if( a<b )
		return -1;
	else
		return 0;
}

void listDeleteFunction( void *a ){

}

void ValidateShadowing( void *dummy, void* context ) {
	DiscretisationContext *self = (DiscretisationContext*)context;
	Swarm*                          swarm = (Swarm*) LiveComponentRegister_Get( self->CF->LCRegister, "swarm" );

	Swarm_UpdateAllParticleOwners( swarm );

	return 0;
	if(swarm->nProc > 1)
	{
		int ii = 0, jj = 0;
		ShadowInfo*			cellShadowInfo = CellLayout_GetShadowInfo( swarm->cellLayout );
		ProcNbrInfo*		procNbrInfo = cellShadowInfo->procNbrInfo;
		
		{
			MemoryPool *requestPool = MemoryPool_New( MPI_Request, 100, 10 );
			MemoryPool *particlePool = MemoryPool_NewFunc( swarm->particleExtensionMgr->finalSize, 100, 10 );
			LinkedList *list = LinkedList_New( listCompareFunction, NULL, NULL, listDeleteFunction, LINKEDLIST_UNSORTED );
			LinkedList *particleList= LinkedList_New( listCompareFunction, NULL, NULL, listDeleteFunction, LINKEDLIST_UNSORTED );
			LinkedListIterator *iter = LinkedListIterator_New( list );
			void *data = NULL;

			{
				Neighbour_Index		        nbr_I;
				int i = 0, j = 0;
				int shadowCell = 0;

				for ( nbr_I=0; nbr_I < procNbrInfo->procNbrCnt; nbr_I++ ) {
					for( i=0; i<cellShadowInfo->procShadowCnt[nbr_I]; i++ ){
						MPI_Request *req;
						char *particle = NULL;
						int proc = 0;
					
						proc = procNbrInfo->procNbrTbl[nbr_I];

						shadowCell = cellShadowInfo->procShadowTbl[nbr_I][i];
						for( j=0; j<swarm->shadowCellParticleCountTbl[shadowCell]; j++ ){
							
							req = MemoryPool_NewObject( MPI_Request, requestPool );
							particle = MemoryPool_NewObjectFunc( swarm->particleExtensionMgr->finalSize, particlePool );
							LinkedList_InsertNode( list, req, sizeof(void*) );
							LinkedList_InsertNode( particleList, particle, sizeof(void*) );
						
							MPI_Irecv( particle, swarm->particleExtensionMgr->finalSize, MPI_BYTE, proc, 2001, MPI_COMM_WORLD, req );
						}
					}
				}
			
			}

			{
				Neighbour_Index		        nbr_I;
				int i = 0, j = 0;
				int shadowedCell = 0;
				char *array = malloc( swarm->particleExtensionMgr->finalSize );

				for ( nbr_I=0; nbr_I < procNbrInfo->procNbrCnt; nbr_I++ ) {
					for( i=0; i<cellShadowInfo->procShadowedCnt[nbr_I]; i++ ){
						MPI_Request *req;
						int proc = 0;

						proc = procNbrInfo->procNbrTbl[nbr_I];

						shadowedCell = cellShadowInfo->procShadowedTbl[nbr_I][i];
						for( j=0; j<swarm->cellParticleCountTbl[shadowedCell]; j++ ){
							int pIndex = 0;

							pIndex = swarm->cellParticleTbl[shadowedCell][j];
							Swarm_CopyParticleOffSwarm( swarm,
								array, 0,
								pIndex );

							req = MemoryPool_NewObject( MPI_Request, requestPool );
							LinkedList_InsertNode( list, req, sizeof(void*) );
							MPI_Isend( array, swarm->particleExtensionMgr->finalSize, MPI_BYTE, proc, 2001, MPI_COMM_WORLD, req );
						}
					}
				}
				free( array );
			}

			for( data=LinkedListIterator_First(iter); data!=NULL; data=LinkedListIterator_Next(iter) ){
				MPI_Status s;
				MPI_Wait( ((MPI_Request*)data), &s );
			}

			{
				LinkedListIterator *pIter = LinkedListIterator_New( particleList );
				for( data=LinkedListIterator_First(pIter); data!=NULL; data=LinkedListIterator_Next(pIter) ){
					Particle *p = (Particle*)data;
					int found = 0;
					double epsilon=1e-5;

					for( ii=0; ii<procNbrInfo->procNbrCnt; ii++ ){
						for( jj=0; jj<cellShadowInfo->procShadowCnt[ii]; jj++ ){
							int shadowCell = 0;
							int kk = 0;

							shadowCell = cellShadowInfo->procShadowTbl[ii][jj];
							for( kk=0; kk<swarm->shadowCellParticleCountTbl[shadowCell]; kk++ ){
								int shadowParticleIdx = 0;
								GlobalParticle *gp = NULL;

								shadowParticleIdx = swarm->shadowCellParticleTbl[shadowCell][kk];
								gp = ((GlobalParticle*)Swarm_ShadowParticleAt( swarm, shadowParticleIdx));
								if( (fabs(p->coord[0]-gp->coord[0])<epsilon) && (fabs(p->coord[1]-gp->coord[1])<epsilon) && (fabs(p->coord[2]-gp->coord[2])<epsilon) ){
									found = 1;
								}
							}
						}
					}
					Journal_Firewall( found, swarm->debug, "Shadow particle validation failed" );
				}
			}
			
			Stg_Class_Delete( list );
			Stg_Class_Delete( iter );
			Stg_Class_Delete( requestPool );
			Stg_Class_Delete( particlePool );
			Stg_Class_Delete( particleList );
		}
	}
}

