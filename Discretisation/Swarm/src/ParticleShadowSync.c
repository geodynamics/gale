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
** $Id: ParticleShadowSync.c 4001 2007-02-09 01:26:29Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include "Base/Base.h"

#include "Discretisation/Geometry/Geometry.h"
#include "Discretisation/Shape/Shape.h"
#include "Discretisation/Mesh/Mesh.h"
#include "Discretisation/Utils/Utils.h"

#include "ShadowInfo.h"
#include "types.h"
#include "shortcuts.h"
#include "ParticleCommHandler.h"
#include "ParticleShadowSync.h"

#include "SwarmClass.h"
#include "CellLayout.h"
#include "ElementCellLayout.h"
#include "StandardParticle.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

const Type ParticleShadowSync_Type = "ParticleShadowSync";

/* MPI tags */
static const int SHADOW_PARTICLE_COUNTS_PER_CELL = 10;
static const int SHADOW_PARTICLES = 20;

void* ParticleShadowSync_DefaultNew( Name name )
{
	return _ParticleShadowSync_New( sizeof(ParticleShadowSync), ParticleShadowSync_Type,
		_ParticleShadowSync_Delete, _ParticleShadowSync_Print, _ParticleShadowSync_CopyFunc,
		(Stg_Component_DefaultConstructorFunction*)ParticleShadowSync_DefaultNew,
		_ParticleShadowSync_Construct, _ParticleShadowSync_Build, _ParticleShadowSync_Initialise,
		_ParticleShadowSync_Execute, _ParticleShadowSync_Destroy, name, False, NULL );
}


ParticleShadowSync* ParticleShadowSync_New( 
		Name name,
		void* swarm
		)
{
		return _ParticleShadowSync_New( sizeof(ParticleShadowSync), ParticleShadowSync_Type,
		_ParticleShadowSync_Delete, _ParticleShadowSync_Print, _ParticleShadowSync_CopyFunc,
		(Stg_Component_DefaultConstructorFunction*)ParticleShadowSync_DefaultNew,
		_ParticleShadowSync_Construct, _ParticleShadowSync_Build, _ParticleShadowSync_Initialise,
		_ParticleShadowSync_Execute, _ParticleShadowSync_Destroy, name, True, swarm );
}


ParticleShadowSync* _ParticleShadowSync_New( 
		SizeT                                                           _sizeOfSelf,
		Type                                                            type,
		Stg_Class_DeleteFunction*                                       _delete,
		Stg_Class_PrintFunction*                                        _print,
		Stg_Class_CopyFunction*                                         _copy, 
		Stg_Component_DefaultConstructorFunction*                       _defaultConstructor,
		Stg_Component_ConstructFunction*                                _construct,
		Stg_Component_BuildFunction*                                    _build,
		Stg_Component_InitialiseFunction*                               _initialise,
		Stg_Component_ExecuteFunction*                                  _execute,
		Stg_Component_DestroyFunction*                                  _destroy,
		Name                                                            name,
		Bool                                                            initFlag,
		void*                                                           swarm
		)
{
	ParticleShadowSync* self;
	
	/* Allocate memory */
	
	self = (ParticleShadowSync*)_ParticleCommHandler_New( sizeof(ParticleShadowSync), ParticleShadowSync_Type,
		_ParticleShadowSync_Delete, _ParticleShadowSync_Print, _ParticleShadowSync_CopyFunc,
		(Stg_Component_DefaultConstructorFunction*)ParticleShadowSync_DefaultNew,
		_ParticleShadowSync_Construct, _ParticleShadowSync_Build, _ParticleShadowSync_Initialise,
		_ParticleShadowSync_Execute, _ParticleShadowSync_Destroy, name, initFlag,
		_ParticleCommHandler_AllocateOutgoingCountArrays,
		NULL,
		_ParticleCommHandler_FreeOutgoingArrays,
		_ParticleCommHandler_AllocateIncomingCountArrays,
		NULL,
		_ParticleCommHandler_FreeIncomingArrays,
		_ParticleCommHandler_BeginReceiveOfIncomingParticleCounts,
		_ParticleShadowSync_FinishReceiveOfIncomingParticleCounts,
		_ParticleShadowSync_BeginReceiveOfIncomingParticles,
		_ParticleShadowSync_FinishReceiveOfIncomingParticles,
		_ParticleShadowSync_SendParticleTotalsInShadowCellsToNbrs,
		_ParticleShadowSync_SendShadowParticles,
		_ParticleCommHandler_ConfirmOutgoingSendsCompleted,
		ParticleShadowSync_HandleParticleMovementBetweenProcs,
		swarm );
	
	/* General info */
	/* Virtual info */
	
	/* ParticleShadowSync info */
	if( initFlag ){
		_ParticleShadowSync_Init( self, swarm );
	}
	
	return self;
}


void _ParticleShadowSync_Init(
		ParticleShadowSync*     self,
		void*                   swarm
		)
{
	_ParticleCommHandler_Init( (ParticleCommHandler*)self, swarm );
	self->particlesOutsideDomainIndices = NULL;
	_ParticleCommHandler_ZeroShadowCommStrategyCounters( (ParticleCommHandler*)self );
}


void _ParticleShadowSync_Delete(void* pCommsHandler ) {

	_ParticleCommHandler_Delete( pCommsHandler );
}


void _ParticleShadowSync_Print( void* pCommsHandler, Stream* stream ) {
	ParticleShadowSync*	self = (ParticleShadowSync*)pCommsHandler;
	
	/* General info */
	Journal_Printf( stream, "ParticleShadowSync (ptr): %p\n", self );
	
	/* Parent class info */
	_ParticleCommHandler_Print( self, stream );

	/* Virtual info */
}


void* _ParticleShadowSync_CopyFunc( void* particleMovementHandler, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	
	return _ParticleCommHandler_Copy( particleMovementHandler, dest, deep,
		   nameExt, ptrMap );
}

void _ParticleShadowSync_Construct( void* pCommsHandler, Stg_ComponentFactory* cf, void* data ){
	Swarm *swarm = NULL;
	ParticleShadowSync *self = (ParticleShadowSync*)pCommsHandler;

	swarm =  Stg_ComponentFactory_ConstructByKey(  cf,  self->name,  Swarm_Type, Swarm,  True, data );
	assert( swarm );

	self->isConstructed = True;
	_ParticleShadowSync_Init( self, swarm );
}
	
void _ParticleShadowSync_Build( void* pCommsHandler, void *data ){
	ParticleShadowSync		*self = (ParticleShadowSync*)pCommsHandler;

	Stg_Component_Build( self->swarm, NULL, True );
	_Swarm_BuildShadowParticles( self->swarm );
}
	
void _ParticleShadowSync_Initialise( void* pCommsHandler, void *data ){
	
}
	
void _ParticleShadowSync_Execute( void* pCommsHandler, void *data ){
	ParticleCommHandler*	self = (ParticleCommHandler*)pCommsHandler;
	
	Stg_Component_Build( self, NULL, False );
	self->_commFunction( self );	
}

void _ParticleShadowSync_Destroy( void* pCommsHandler, void *data ){
	
}

void _ParticleShadowSync_SendParticleTotalsInShadowCellsToNbrs( ParticleCommHandler* self )
{	
	Processor_Index			proc_I;
	ShadowInfo*			cellShadowInfo = CellLayout_GetShadowInfo( self->swarm->cellLayout );
	ProcNbrInfo*			procNbrInfo = cellShadowInfo->procNbrInfo;
	Neighbour_Index		nbr_I, nbrCount;
	Cell_Index              lCellID, cellParticleCount, shadowedCell_I;	

	Journal_DPrintfL( self->debug, 1, "In %s():\n", __func__ );
	Stream_IndentBranch( Swarm_Debug );

	self->shadowParticlesLeavingMeTotalCount = 0;

	nbrCount = procNbrInfo->procNbrCnt;
	for( nbr_I = 0 ; nbr_I < nbrCount ; nbr_I++ ) {
		proc_I = procNbrInfo->procNbrTbl[nbr_I];
		
		self->shadowParticlesLeavingMeTotalCounts[nbr_I] = 0;
		for( shadowedCell_I = 0 ; shadowedCell_I < cellShadowInfo->procShadowedCnt[ nbr_I ]; shadowedCell_I++ ) {
			lCellID = cellShadowInfo->procShadowedTbl[nbr_I][shadowedCell_I];
			cellParticleCount = self->swarm->cellParticleCountTbl[ lCellID ];
			self->shadowParticlesLeavingMeCountsPerCell[nbr_I][shadowedCell_I] = cellParticleCount;

			self->shadowParticlesLeavingMeTotalCounts[nbr_I] += cellParticleCount;
		}
		MPI_Ssend( self->shadowParticlesLeavingMeCountsPerCell[nbr_I], cellShadowInfo->procShadowCnt[nbr_I], MPI_UNSIGNED,
			proc_I, SHADOW_PARTICLE_COUNTS_PER_CELL, self->swarm->comm );
	}

	Stream_UnIndentBranch( Swarm_Debug );
}

void _ParticleShadowSync_FinishReceiveOfIncomingParticleCounts( ParticleCommHandler* self ) {
	MPI_Status			status;
	Processor_Index			proc_I;
	ShadowInfo*		        cellShadowInfo = CellLayout_GetShadowInfo( self->swarm->cellLayout );
	ProcNbrInfo*		        procNbrInfo = cellShadowInfo->procNbrInfo;
	Neighbour_Index		        nbr_I;
	int i = 0;

	self->swarm->shadowParticleCount = 0;
	/* TODO: may be worth converting the below into an MPI_Test loop */
	for ( nbr_I=0; nbr_I < procNbrInfo->procNbrCnt; nbr_I++ ) {
		proc_I = procNbrInfo->procNbrTbl[nbr_I];
		
		MPI_Wait( self->particlesArrivingFromNbrShadowCellCountsHandles[nbr_I], &status );

		self->particlesArrivingFromNbrShadowCellsTotalCounts[nbr_I] = 0;
		for( i=0; i<cellShadowInfo->procShadowCnt[nbr_I]; i++ ){
			Index shadowCell = 0;

			shadowCell = cellShadowInfo->procShadowTbl[nbr_I][i];
			self->swarm->shadowCellParticleCountTbl[shadowCell] = 0;

			self->swarm->shadowCellParticleTbl[shadowCell] = Memory_Realloc_Array( self->swarm->shadowCellParticleTbl[shadowCell], Particle_Index,
																	(self->particlesArrivingFromNbrShadowCellCounts[nbr_I][i]) );

			self->particlesArrivingFromNbrShadowCellsTotalCounts[nbr_I] += self->particlesArrivingFromNbrShadowCellCounts[nbr_I][i];
		}
	
		self->swarm->shadowParticleCount+= self->particlesArrivingFromNbrShadowCellsTotalCounts[nbr_I];
	}
}

void _ParticleShadowSync_BeginReceiveOfIncomingParticles( ParticleCommHandler* pCommHandler ) {
	ParticleShadowSync		*self = (ParticleShadowSync*)pCommHandler;
	ShadowInfo*		        cellShadowInfo = CellLayout_GetShadowInfo( self->swarm->cellLayout );
	ProcNbrInfo*		        procNbrInfo = cellShadowInfo->procNbrInfo;
	long                           incomingViaShadowArrayBytes = 0;
	Neighbour_Index		        nbr_I;
	Processor_Index			proc_I;
	int i = 0;
	char* recvLocation  = NULL;

	self->particlesArrivingFromNbrShadowCellsHandles = Memory_Alloc_Array_Unnamed( MPI_Request*, procNbrInfo->procNbrCnt );
	for( i=0; i<procNbrInfo->procNbrCnt; i++ ){
		self->particlesArrivingFromNbrShadowCellsHandles[i] = Memory_Alloc_Array_Unnamed( MPI_Request, 1 );
	}

	self->swarm->shadowParticles = Memory_Realloc( self->swarm->shadowParticles,
			self->swarm->particleExtensionMgr->finalSize*(self->swarm->shadowParticleCount) );
	
	recvLocation = (char*)self->swarm->shadowParticles;
	for ( nbr_I=0; nbr_I < procNbrInfo->procNbrCnt; nbr_I++ ) {
		
		if ( self->particlesArrivingFromNbrShadowCellsTotalCounts[nbr_I] != 0 ) {
			proc_I = procNbrInfo->procNbrTbl[nbr_I];

			/* start non-blocking recv of particles */
			incomingViaShadowArrayBytes = self->swarm->particleExtensionMgr->finalSize * 
				self->particlesArrivingFromNbrShadowCellsTotalCounts[nbr_I];
			
			/*printf( "receiving %ld bytes\n", incomingViaShadowArrayBytes );*/
			MPI_Irecv( recvLocation, incomingViaShadowArrayBytes, MPI_BYTE,
				proc_I, SHADOW_PARTICLES, self->swarm->comm,
				self->particlesArrivingFromNbrShadowCellsHandles[nbr_I] );
			
			recvLocation += incomingViaShadowArrayBytes;
		}
	}
}

void _ParticleShadowSync_FinishReceiveOfIncomingParticles( ParticleCommHandler* pCommHandler ) {
	ParticleShadowSync *self = (ParticleShadowSync*)pCommHandler;
	MPI_Status			status;
	Processor_Index			proc_I;
	ShadowInfo*		        cellShadowInfo = CellLayout_GetShadowInfo( self->swarm->cellLayout );
	ProcNbrInfo*		        procNbrInfo = cellShadowInfo->procNbrInfo;
	Neighbour_Index		        nbr_I;
	int i = 0, j = 0;
	int shadowParticleCounter;
	int shadowCell = 0;

	shadowParticleCounter = 0;
	/* TODO: may be worth converting the below into an MPI_Test loop */
	for ( nbr_I=0; nbr_I < procNbrInfo->procNbrCnt; nbr_I++ ) {
		proc_I = procNbrInfo->procNbrTbl[nbr_I];
		
		if( self->particlesArrivingFromNbrShadowCellsTotalCounts[nbr_I] > 0 ){
			MPI_Wait( self->particlesArrivingFromNbrShadowCellsHandles[nbr_I], &status );
		}
	}

	for ( nbr_I=0; nbr_I < procNbrInfo->procNbrCnt; nbr_I++ ) {
		for( i=0; i<cellShadowInfo->procShadowCnt[nbr_I]; i++ ){
			
			shadowCell = cellShadowInfo->procShadowTbl[nbr_I][i];
			for( j=0; j<self->particlesArrivingFromNbrShadowCellCounts[nbr_I][i]; j++ ){

				Swarm_AddShadowParticleToShadowCell( self->swarm, shadowCell, shadowParticleCounter );
				shadowParticleCounter++;
			}
		}
	}
}

void _ParticleShadowSync_SendShadowParticles( ParticleCommHandler *self )
{
	ShadowInfo*		        cellShadowInfo = CellLayout_GetShadowInfo( self->swarm->cellLayout );
	ProcNbrInfo*		        procNbrInfo = cellShadowInfo->procNbrInfo;
	Processor_Index			proc_I;
	int i = 0, j = 0, k = 0, cell = 0;
	unsigned int	arrayIndex = 0;
	long			arraySize = 0;
	unsigned int	pIndex = 0;

	self->shadowParticlesLeavingMeHandles = Memory_Alloc_Array_Unnamed( MPI_Request*, procNbrInfo->procNbrCnt );
	self->shadowParticlesLeavingMe = Memory_Alloc_Array_Unnamed( Particle*, procNbrInfo->procNbrCnt );

	for( i=0; i<procNbrInfo->procNbrCnt; i++ ){
		proc_I = procNbrInfo->procNbrTbl[i];

		if( self->shadowParticlesLeavingMeTotalCounts[i] != 0 ){

			self->shadowParticlesLeavingMeHandles[i] = Memory_Alloc_Array_Unnamed( MPI_Request, 1 );

			arraySize =  self->swarm->particleExtensionMgr->finalSize * self->shadowParticlesLeavingMeTotalCounts[i];
			self->shadowParticlesLeavingMe[i] = Memory_Alloc_Bytes( arraySize, "Particle", "pCommHandler->outgoingPArray" );
			memset( self->shadowParticlesLeavingMe[i], 0, arraySize );

			arrayIndex = 0;
			for( j=0; j<cellShadowInfo->procShadowedCnt[i]; j++ ){
				cell = cellShadowInfo->procShadowedTbl[i][j];

				for( k=0; k<self->swarm->cellParticleCountTbl[cell]; k++ ){
					pIndex = self->swarm->cellParticleTbl[cell][k];
				
					Swarm_CopyParticleOffSwarm( self->swarm,
							self->shadowParticlesLeavingMe[i], arrayIndex++,
							pIndex );
				}
			}
			
			/*printf( "sending %ld bytes\n", self->shadowParticlesLeavingMeTotalCounts[i] * self->swarm->particleExtensionMgr->finalSize );*/
			MPI_Issend( self->shadowParticlesLeavingMe[i],
			self->shadowParticlesLeavingMeTotalCounts[i] * self->swarm->particleExtensionMgr->finalSize,
			MPI_BYTE, proc_I, SHADOW_PARTICLES, self->swarm->comm,
			self->shadowParticlesLeavingMeHandles[i] );
		}
	}
}

void ParticleShadowSync_HandleParticleMovementBetweenProcs( ParticleCommHandler* pCommsHandler ) {
	ParticleShadowSync*	self = (ParticleShadowSync*)pCommsHandler;

	Journal_DPrintfL( self->debug, 1, "In %s(), for swarm \"%s\":\n", __func__, self->swarm->name );
	if ( 1 == self->swarm->nProc ) {
		Journal_DPrintfL( self->debug, 1, "Serial run -> nothing to communicate in %s, returning.\n", __func__ );
		Stream_UnIndentBranch( Swarm_Debug );
		return;
	}

	Stream_IndentBranch( Swarm_Debug );
	
	if ( self->swarm->cellShadowCount > 0 ) {
		
		self->allocateIncomingCountArrays( (ParticleCommHandler*)self );
		self->allocateOutgoingCountArrays( (ParticleCommHandler*)self );



		self->beginReceiveOfIncomingParticleCounts( (ParticleCommHandler*)self );


		self->sendOutgoingParticleCounts( (ParticleCommHandler*)self );

		
		self->finishReceiveOfIncomingParticleCounts( (ParticleCommHandler*)self );


		self->beginReceiveOfIncomingParticles( (ParticleCommHandler*)self );
		self->beginSendingParticles( (ParticleCommHandler*)self );


		self->confirmOutgoingSendsCompleted( (ParticleCommHandler*)self );
		self->finishReceiveOfIncomingParticlesAndUpdateIndices( (ParticleCommHandler*)self );


		self->freeIncomingArrays( (ParticleCommHandler*)self );
		self->freeOutgoingArrays( (ParticleCommHandler*)self );
	}
	
	MPI_Barrier( self->swarm->comm );

	_ParticleCommHandler_ZeroShadowCommStrategyCounters( (ParticleCommHandler*)self );
	
	Stream_UnIndentBranch( Swarm_Debug );
}


