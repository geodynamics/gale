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
** $Id: ParticleCommHandler.c 4119 2007-05-22 07:35:46Z RaquibulHassan $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include "Base/Base.h"

#include "Discretisation/Geometry/Geometry.h"
#include "Discretisation/Shape/Shape.h"
#include "Discretisation/Mesh/Mesh.h"
#include "Discretisation/Utils/Utils.h"

#include "types.h"
#include "shortcuts.h"
#include "ParticleCommHandler.h"

#include "SwarmClass.h"
#include "ShadowInfo.h"
#include "CellLayout.h"
#include "ElementCellLayout.h"
#include "StandardParticle.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

const Type ParticleCommHandler_Type = "ParticleCommHandler";

/* MPI tags */
static const int SHADOW_PARTICLE_COUNTS_PER_CELL = 10;
static const int SHADOW_PARTICLES = 20;

ParticleCommHandler* _ParticleCommHandler_New( 
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
		ParticleCommHandler_AllocateOutgoingCountArrays						*_allocateOutgoingCountArrays,
		ParticleCommHandler_AllocateOutgoingParticleArrays						*_allocateOutgoingParticleArrays,
		ParticleCommHandler_FreeOutgoingArrays									*_freeOutgoingArrays,
		ParticleCommHandler_AllocateIncomingCountArrays						*_allocateIncomingCountArrays,
		ParticleCommHandler_AllocateIncomingParticleArrays						*_allocateIncomingParticleArrays,
		ParticleCommHandler_FreeIncomingArrays									*_freeIncomingArrays,
		ParticleCommHandler_BeginReceiveOfIncomingParticleCounts				*_beginReceiveOfIncomingParticleCounts,
		ParticleCommHandler_FinishReceiveOfIncomingParticleCounts				*_finishReceiveOfIncomingParticleCounts,
		ParticleCommHandler_BeginReceiveOfIncomingParticles					*_beginReceiveOfIncomingParticles,
		ParticleCommHandler_FinishReceiveOfIncomingParticlesAndUpdateIndices	*_finishReceiveOfIncomingParticlesAndUpdateIndices,
		ParticleCommHandler_SendOutgoingParticleCounts							*_sendOutgoingParticleCounts,
		ParticleCommHandler_BeginSendingParticles								*_beginSendingParticles,
		ParticleCommHandler_ConfirmOutgoingSendsCompleted						*_confirmOutgoingSendsCompleted,
		ParticleCommHandler_CommFunction										*_commFunction )
{
	ParticleCommHandler* self;
	
	/* Allocate memory */
	self = (ParticleCommHandler*)_Stg_Component_New( _sizeOfSelf, type, _delete, _print, _copy, _defaultConstructor,
		   _construct, _build, _initialise, _execute, _destroy, name, NON_GLOBAL );
	
	/* General info */
	/* Virtual info */
	self->allocateOutgoingCountArrays = _allocateOutgoingCountArrays;
	self->allocateOutgoingParticleArrays	= _allocateOutgoingParticleArrays,
	self->freeOutgoingArrays = _freeOutgoingArrays,
	self->allocateIncomingCountArrays = _allocateIncomingCountArrays,
	self->allocateIncomingParticleArrays = _allocateIncomingParticleArrays,
	self->freeIncomingArrays = _freeIncomingArrays,
	self->beginReceiveOfIncomingParticleCounts = _beginReceiveOfIncomingParticleCounts,
	self->finishReceiveOfIncomingParticleCounts = _finishReceiveOfIncomingParticleCounts,
	self->beginReceiveOfIncomingParticles = _beginReceiveOfIncomingParticles,
	self->finishReceiveOfIncomingParticlesAndUpdateIndices = _finishReceiveOfIncomingParticlesAndUpdateIndices,
	self->sendOutgoingParticleCounts = _sendOutgoingParticleCounts,
	self->beginSendingParticles = _beginSendingParticles,
	self->confirmOutgoingSendsCompleted = _confirmOutgoingSendsCompleted,
	self->_commFunction = _commFunction;
	
	/* ParticleCommHandler info */
	if( initFlag ){
		_ParticleCommHandler_Init( self );
	}
	
	return self;
}


void _ParticleCommHandler_Init(
		ParticleCommHandler*     self )
{
	self->isConstructed = True;
	self->swarm = NULL;
	self->debug = Stream_RegisterChild( Swarm_Debug, self->type );

	_ParticleCommHandler_ZeroShadowCommStrategyCounters( self );

	self->shadowParticlesLeavingMeIndices = NULL;
	self->shadowParticlesLeavingMeCountsPerCell = NULL;
	self->shadowParticlesLeavingMeTotalCounts = NULL;
	self->shadowParticlesLeavingMe = NULL;
	self->shadowParticlesLeavingMeHandles = NULL;

	self->particlesArrivingFromNbrShadowCellCounts = NULL;
	self->particlesArrivingFromNbrShadowCellCountsHandles = NULL;
	self->particlesArrivingFromNbrShadowCellsTotalCounts = NULL;
	self->particlesArrivingFromNbrShadowCells = NULL;
	self->particlesArrivingFromNbrShadowCellsHandles = NULL;
}


void _ParticleCommHandler_Delete(void* pCommsHandler ) {
	/* Nothing to do, as the temporary arrays are deleted each step as soon as they are finished with. */
	_Stg_Class_Delete( pCommsHandler );
}


void _ParticleCommHandler_Print( void* pCommsHandler, Stream* stream ) {
	ParticleCommHandler*	self = (ParticleCommHandler*)pCommsHandler;
	
	/* General info */
	Journal_Printf( stream, "ParticleCommHandler (ptr): %p\n", self );
	
	/* Parent class info */
	_Stg_Class_Print( self, stream );

	/* Virtual info */
	Journal_Printf( stream, "self->_commFunction(ptr): %p\n", 
		self->_commFunction );

	/* class info */	
	Journal_Printf( stream, "self->swarm (ptr): %p\n", self->swarm ); 
	Journal_Printf( stream, "self->shadowParticlesLeavingMeTotalCount: %d\n", self->shadowParticlesLeavingMeTotalCount ); 
	Journal_Printf( stream, "self->shadowParticlesLeavingMeUnfilledCount: %d\n", self->shadowParticlesLeavingMeUnfilledCount ); 
	Journal_Printf( stream, "self->currShadowParticleLeavingMeIndex: %d\n", self->currShadowParticleLeavingMeIndex ); 
	Journal_Printf( stream, "self->shadowParticlesLeavingMeIndices (ptr): %p\n", self->shadowParticlesLeavingMeIndices ); 
	Journal_Printf( stream, "self->currParticleLeavingMyDomainIndex: %d\n", self->currParticleLeavingMyDomainIndex ); 
}


void* _ParticleCommHandler_Copy( void* particleCommHandler, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {

	// TODO: put back in.
	assert(0);

}

void _ParticleCommHandler_Construct( void* pCommsHandler, Stg_ComponentFactory* cf, void* data ){
	// TODO : will be necessary when we allow this to be constructed as a component
	assert(0);
}
	
void _ParticleCommHandler_Build( void* pCommsHandler, void *data ){
	
}
	
void _ParticleCommHandler_Initialise( void* pCommsHandler, void *data ){
	
}
	
void _ParticleCommHandler_Execute( void* pCommsHandler, void *data ){
}

void _ParticleCommHandler_Destroy( void* pCommsHandler, void *data ){
	
}

void _ParticleCommHandler_AllocateOutgoingParticlesArrays( ParticleCommHandler* self ) {
	ShadowInfo*		        cellShadowInfo = CellLayout_GetShadowInfo( self->swarm->cellLayout );
	ProcNbrInfo*		        procNbrInfo = cellShadowInfo->procNbrInfo;
	Neighbour_Index		        nbrCount = procNbrInfo->procNbrCnt;
	Neighbour_Index		        nbr_I;
	SizeT				particlesArrayBytes;

	/* Allocate the list of shadow particle indices to all nbrs */
	self->shadowParticlesLeavingMeIndices = Memory_Alloc_Array( Particle_Index, self->shadowParticlesLeavingMeTotalCount,
		"self->shadowParticlesLeavingMeIndices" );
	self->shadowParticlesLeavingMeHandles = Memory_Alloc_Array(
		MPI_Request*, procNbrInfo->procNbrCnt,
		"ParticleCommHandler->shadowParticlesLeavingMeHandles");
	self->shadowParticlesLeavingMe = Memory_Alloc_Array(
		Particle*, procNbrInfo->procNbrCnt,
		"ParticleCommHandler->shadowParticlesLeavingMe" );

	for ( nbr_I=0; nbr_I < nbrCount; nbr_I++ ) {

		if ( self->shadowParticlesLeavingMeTotalCounts[nbr_I] == 0 ) {
			/* If we're not sending any particles to this proc, skip to next */
			self->shadowParticlesLeavingMeHandles[nbr_I] = NULL;
			self->shadowParticlesLeavingMe[nbr_I] = NULL;
		}
		else {
			self->shadowParticlesLeavingMeHandles[nbr_I] = Memory_Alloc( MPI_Request,
				"ParticleCommHandler->shadowParticlesLeavingMeHandles[]" );
			particlesArrayBytes = self->swarm->particleExtensionMgr->finalSize * 
				self->shadowParticlesLeavingMeTotalCounts[nbr_I];
			self->shadowParticlesLeavingMe[nbr_I] = Memory_Alloc_Bytes( particlesArrayBytes,
				"Particle", "ParticleCommHandler->shadowParticlesLeavingMe[]" );
		}
	}
}


void _ParticleCommHandler_BeginReceiveOfIncomingParticleCounts( ParticleCommHandler* self ) {
	ShadowInfo*		        cellShadowInfo = CellLayout_GetShadowInfo( self->swarm->cellLayout );
	ProcNbrInfo*		        procNbrInfo = cellShadowInfo->procNbrInfo;
	Neighbour_Index		        nbr_I;
	Processor_Index		        proc_I;

	for ( nbr_I=0; nbr_I < procNbrInfo->procNbrCnt; nbr_I++ ) {
		proc_I = procNbrInfo->procNbrTbl[nbr_I];
		MPI_Irecv( self->particlesArrivingFromNbrShadowCellCounts[nbr_I], cellShadowInfo->procShadowedCnt[nbr_I],
			MPI_UNSIGNED, proc_I, SHADOW_PARTICLE_COUNTS_PER_CELL, self->swarm->comm,
			self->particlesArrivingFromNbrShadowCellCountsHandles[nbr_I] );
	}
}


void _ParticleCommHandler_FinishReceiveOfIncomingParticleCounts( ParticleCommHandler* self ) {
	Cell_ShadowTransferIndex	stCell_I;
	MPI_Status			status;
	Processor_Index			proc_I;
	Particle_Index			incomingCellParticleCount = 0;
	ShadowInfo*		        cellShadowInfo = CellLayout_GetShadowInfo( self->swarm->cellLayout );
	ProcNbrInfo*		        procNbrInfo = cellShadowInfo->procNbrInfo;
	Neighbour_Index		        nbr_I;

	/* TODO: may be worth converting the below into an MPI_Test loop */
	for ( nbr_I=0; nbr_I < procNbrInfo->procNbrCnt; nbr_I++ ) {
		proc_I = procNbrInfo->procNbrTbl[nbr_I];
		
		MPI_Wait( self->particlesArrivingFromNbrShadowCellCountsHandles[nbr_I], &status );

		Journal_DPrintfL( self->debug, 1, "proc %d: recv shadow counts from nbr %d (rank %d):\n",
			self->swarm->myRank, nbr_I, proc_I );
		Journal_DPrintfL( self->debug, 2, "\tare [" );
		self->particlesArrivingFromNbrShadowCellsTotalCounts[nbr_I] = 0;

		for ( stCell_I=0; stCell_I < cellShadowInfo->procShadowedCnt[nbr_I]; stCell_I++ ) {
			incomingCellParticleCount = self->particlesArrivingFromNbrShadowCellCounts[nbr_I][stCell_I];
			Journal_DPrintfL( self->debug, 2, "%d, ", incomingCellParticleCount );
			self->particlesArrivingFromNbrShadowCellsTotalCounts[nbr_I] += incomingCellParticleCount;
		}
		Journal_DPrintfL( self->debug, 2, "]\n" );
		Journal_DPrintfL( self->debug, 1, "(Proc %d):....totalled to %d\n", self->swarm->myRank,
			self->particlesArrivingFromNbrShadowCellsTotalCounts[nbr_I] );
	}
}


void _ParticleCommHandler_AllocateIncomingParticlesArrays( ParticleCommHandler* self ) {
	ShadowInfo*		        cellShadowInfo = CellLayout_GetShadowInfo( self->swarm->cellLayout );
	ProcNbrInfo*		        procNbrInfo = cellShadowInfo->procNbrInfo;
	Neighbour_Index		        nbr_I;
	SizeT                           incomingViaShadowArrayBytes = 0;

	self->particlesArrivingFromNbrShadowCellsHandles = Memory_Alloc_Array(
		MPI_Request*, procNbrInfo->procNbrCnt,
		"ParticleCommHandler->particlesArrivingFromNbrShadowCellsHandles");
	self->particlesArrivingFromNbrShadowCells = Memory_Alloc_Array(
		Particle*, procNbrInfo->procNbrCnt,
		"ParticleCommHandler->particlesArrivingFromNbrShadowCells" );

	for ( nbr_I=0; nbr_I < procNbrInfo->procNbrCnt; nbr_I++ ) {

		if ( self->particlesArrivingFromNbrShadowCellsTotalCounts[nbr_I] == 0 ) {
			/* No particles to receive from this proc -> just clear recv ptr */
			self->particlesArrivingFromNbrShadowCellsHandles[nbr_I] = NULL;
			self->particlesArrivingFromNbrShadowCells[nbr_I] = NULL;
		}	
		else { 
			self->particlesArrivingFromNbrShadowCellsHandles[nbr_I] =
				Memory_Alloc(
					MPI_Request,
					"ParticleCommHandler->particlesArrivingFromNbrShadowCellsHandles[]" );

			/* allocate particles recv array to right size */
			incomingViaShadowArrayBytes = self->swarm->particleExtensionMgr->finalSize * 
				self->particlesArrivingFromNbrShadowCellsTotalCounts[nbr_I];
			self->particlesArrivingFromNbrShadowCells[nbr_I] = Memory_Alloc_Bytes( incomingViaShadowArrayBytes,
				"Particle", "particleCommHandler->particlesArrivingFromNbrShadowCells[]" );
		}
	}
}


void _ParticleCommHandler_BeginReceiveOfIncomingParticles( ParticleCommHandler* self ) {
	ShadowInfo*		        cellShadowInfo = CellLayout_GetShadowInfo( self->swarm->cellLayout );
	ProcNbrInfo*		        procNbrInfo = cellShadowInfo->procNbrInfo;
	SizeT                           incomingViaShadowArrayBytes = 0;
	Neighbour_Index		        nbr_I;
	Processor_Index			proc_I;

	for ( nbr_I=0; nbr_I < procNbrInfo->procNbrCnt; nbr_I++ ) {
		if ( self->particlesArrivingFromNbrShadowCellsTotalCounts[nbr_I] != 0 ) {
			proc_I = procNbrInfo->procNbrTbl[nbr_I];

			/* start non-blocking recv of particles */
			incomingViaShadowArrayBytes = self->swarm->particleExtensionMgr->finalSize * 
				self->particlesArrivingFromNbrShadowCellsTotalCounts[nbr_I];
			MPI_Irecv( self->particlesArrivingFromNbrShadowCells[nbr_I], incomingViaShadowArrayBytes, MPI_BYTE,
				proc_I, SHADOW_PARTICLES, self->swarm->comm,
				self->particlesArrivingFromNbrShadowCellsHandles[nbr_I] );
		}
	}
}

void _ParticleCommHandler_AllocateOutgoingCountArrays( ParticleCommHandler* self ) {
	ShadowInfo*		        cellShadowInfo = CellLayout_GetShadowInfo( self->swarm->cellLayout );
	ProcNbrInfo*		        procNbrInfo = cellShadowInfo->procNbrInfo;

	Neighbour_Index		nbrCount = procNbrInfo->procNbrCnt;

	self->shadowParticlesLeavingMeCountsPerCell = Memory_Alloc_2DComplex(
		Particle_Index, nbrCount, cellShadowInfo->procShadowCnt,
		"ParticleCommHandler->shadowParticlesLeavingMeCountsPerCell" );
	self->shadowParticlesLeavingMeTotalCounts = Memory_Alloc_Array( Particle_Index, nbrCount,
		"ParticleCommHandler->shadowParticlesLeavingMeTotalCounts" );
}


void _ParticleCommHandler_ConfirmOutgoingSendsCompleted( ParticleCommHandler* self ) {
	ShadowInfo*		cellShadowInfo = CellLayout_GetShadowInfo( self->swarm->cellLayout );
	ProcNbrInfo*		procNbrInfo = cellShadowInfo->procNbrInfo;
	Neighbour_Index		nbr_I;
	MPI_Status		status;

	/* MPI_Wait for all shadow sends to complete */
	for ( nbr_I=0; nbr_I < procNbrInfo->procNbrCnt; nbr_I++ ) {
		if ( self->shadowParticlesLeavingMeTotalCounts[nbr_I] > 0 ) {
			MPI_Wait( self->shadowParticlesLeavingMeHandles[nbr_I], &status );
		}	
	}
}


void _ParticleCommHandler_AllocateIncomingCountArrays( ParticleCommHandler* self )
{
	ShadowInfo*		cellShadowInfo = CellLayout_GetShadowInfo( self->swarm->cellLayout );
	ProcNbrInfo*		procNbrInfo = cellShadowInfo->procNbrInfo;
	Neighbour_Index		nbr_I;
	Neighbour_Index         nbrCount = procNbrInfo->procNbrCnt;

	self->particlesArrivingFromNbrShadowCellCounts = Memory_Alloc_2DComplex(
		Particle_Index, nbrCount, cellShadowInfo->procShadowedCnt,
		"ParticleCommHandler->particlesArrivingFromNbrShadowCellCounts" );
	self->particlesArrivingFromNbrShadowCellCountsHandles = Memory_Alloc_Array( MPI_Request*, nbrCount,
		"ParticleCommHandler->particlesArrivingFromNbrShadowCellCountsHandles" );
	for ( nbr_I=0; nbr_I < procNbrInfo->procNbrCnt; nbr_I++ ) {
		self->particlesArrivingFromNbrShadowCellCountsHandles[nbr_I] = Memory_Alloc( MPI_Request,
			"ParticleCommHandler->particlesArrivingFromNbrShadowCellCountsHandles[]" );
	}
	self->particlesArrivingFromNbrShadowCellsTotalCounts = Memory_Alloc_Array( Particle_Index, nbrCount,
		"ParticleCommHandler->particlesArrivingFromNbrShadowCellsTotalCounts" );
}


void _ParticleCommHandler_ZeroShadowCommStrategyCounters( ParticleCommHandler* self ) {
	self->currShadowParticleLeavingMeIndex = 0;
	self->shadowParticlesLeavingMeTotalCount = 0;
	self->shadowParticlesLeavingMeUnfilledCount = 0;
}


void _ParticleCommHandler_FreeIncomingArrays( ParticleCommHandler* self ) {
	ShadowInfo*		cellShadowInfo = CellLayout_GetShadowInfo( self->swarm->cellLayout );
	ProcNbrInfo*		procNbrInfo = cellShadowInfo->procNbrInfo;
	Neighbour_Index		nbr_I;

	for ( nbr_I=0; nbr_I < procNbrInfo->procNbrCnt; nbr_I++ ) {
		if ( self->particlesArrivingFromNbrShadowCellsTotalCounts[nbr_I] > 0 ) {

			if( self->particlesArrivingFromNbrShadowCells ){
				Memory_Free( self->particlesArrivingFromNbrShadowCells[nbr_I] );
			}
			Memory_Free( self->particlesArrivingFromNbrShadowCellsHandles[nbr_I] );
		}
	}
	
	if( self->particlesArrivingFromNbrShadowCells ){
		Memory_Free( self->particlesArrivingFromNbrShadowCells );
	}

	Memory_Free( self->particlesArrivingFromNbrShadowCellsHandles );

	for ( nbr_I=0; nbr_I < procNbrInfo->procNbrCnt; nbr_I++ ) {
		Memory_Free( self->particlesArrivingFromNbrShadowCellCountsHandles[nbr_I] );
	}
	Memory_Free( self->particlesArrivingFromNbrShadowCellCountsHandles );

	Memory_Free( self->particlesArrivingFromNbrShadowCellCounts );
	Memory_Free( self->particlesArrivingFromNbrShadowCellsTotalCounts );
}


void _ParticleCommHandler_FreeOutgoingArrays( ParticleCommHandler* self ) {
	ShadowInfo*		cellShadowInfo = CellLayout_GetShadowInfo( self->swarm->cellLayout );
	ProcNbrInfo*		procNbrInfo = cellShadowInfo->procNbrInfo;
	Neighbour_Index		nbr_I;

	Memory_Free( self->shadowParticlesLeavingMeCountsPerCell );
	Memory_Free( self->shadowParticlesLeavingMeTotalCounts );

	if( self->shadowParticlesLeavingMeIndices ){
		Memory_Free( self->shadowParticlesLeavingMeIndices );
		self->shadowParticlesLeavingMeIndices = NULL;
	}

	for ( nbr_I=0; nbr_I < procNbrInfo->procNbrCnt; nbr_I++ ) {
		if ( self->shadowParticlesLeavingMeTotalCounts[nbr_I] > 0 ) {
			Memory_Free( self->shadowParticlesLeavingMe[nbr_I] );
			Memory_Free( self->shadowParticlesLeavingMeHandles[nbr_I] );
		}	
	}
	Memory_Free( self->shadowParticlesLeavingMe );
	Memory_Free( self->shadowParticlesLeavingMeHandles );
}

void _ParticleCommHandler_SendParticleTotalsInShadowCellsToNbrs( ParticleCommHandler* self )
{	
	Cell_ShadowTransferIndex	stCell_I;
	Cell_DomainIndex		dCell_I;
	Index				nbr_I;
	Processor_Index			proc_I;
	Cell_ShadowTransferIndex	shadowCellsToProcCount;
	ShadowInfo*			cellShadowInfo = CellLayout_GetShadowInfo( self->swarm->cellLayout );
	ProcNbrInfo*			procNbrInfo = cellShadowInfo->procNbrInfo;
	Cell_PointIndex			currCellParticleCount;

	Journal_DPrintfL( self->debug, 1, "In %s():\n", __func__ );
	Stream_IndentBranch( Swarm_Debug );

	self->shadowParticlesLeavingMeTotalCount = 0;

	for ( nbr_I = 0; nbr_I < procNbrInfo->procNbrCnt; nbr_I++ ) {
		proc_I = procNbrInfo->procNbrTbl[nbr_I];
	
		shadowCellsToProcCount = cellShadowInfo->procShadowCnt[nbr_I];	
		Journal_DPrintfL( self->debug, 3, "Saving particle count in %d shadow cells going to nbr %d (proc %d):\n\t",
			shadowCellsToProcCount, nbr_I, proc_I );

		self->shadowParticlesLeavingMeTotalCounts[nbr_I] = 0;

		for ( stCell_I=0; stCell_I < shadowCellsToProcCount; stCell_I++ ) {
			dCell_I = cellShadowInfo->procShadowTbl[nbr_I][stCell_I];
			
			currCellParticleCount =  self->swarm->cellParticleCountTbl[dCell_I];
			Journal_DPrintfL( self->debug, 3, "(stCell_I=%d, dCell_I=%d, cnt=%d), ",
				stCell_I, dCell_I, currCellParticleCount );
			self->shadowParticlesLeavingMeCountsPerCell[nbr_I][stCell_I] = currCellParticleCount;
			self->shadowParticlesLeavingMeTotalCounts[nbr_I] += currCellParticleCount;
			self->shadowParticlesLeavingMeTotalCount += currCellParticleCount;
		}	
		Journal_DPrintfL( self->debug, 3, "\n" );

		MPI_Ssend( self->shadowParticlesLeavingMeCountsPerCell[nbr_I], shadowCellsToProcCount, MPI_UNSIGNED,
			proc_I, SHADOW_PARTICLE_COUNTS_PER_CELL, self->swarm->comm );
	}	

	self->shadowParticlesLeavingMeUnfilledCount = self->shadowParticlesLeavingMeTotalCount;
	Stream_UnIndentBranch( Swarm_Debug );
}


void _ParticleCommHandler_BeginSendingParticlesInShadowCellsToNbrs( ParticleCommHandler* self ) {	
	Cell_ShadowTransferIndex	stCell_I;
	Cell_DomainIndex		dCell_I;
	Neighbour_Index			nbr_I;
	Processor_Index			proc_I;
	Cell_ShadowTransferIndex	shadowCellsToProcCount;
	ShadowInfo*			cellShadowInfo = CellLayout_GetShadowInfo( self->swarm->cellLayout );
	ProcNbrInfo*			procNbrInfo = cellShadowInfo->procNbrInfo;
	Neighbour_Index			nbrCount = procNbrInfo->procNbrCnt;
	Particle_InCellIndex		currCellParticleCount;
	Particle_InCellIndex		cParticle_I;
	Particle_Index			lParticle_I;
	Particle_Index			tParticle_I=0; /*Index into the particle transfer array */
	Index				leavingParticle_I=0; /*Index into the array of all leaving particle indices */
	#if DEBUG
	GlobalParticle*                 currParticle;
	#endif
	#if CAUTIOUS
	Bool*                           cellsClearedForTransfer = NULL;
	Neighbour_Index*                cellsClearedForTransferDests = NULL;
	Stream*                         errorStream = Journal_Register( Error_Type, self->type );
	#endif

	Journal_DPrintfL( self->debug, 1, "In %s():\n", __func__ );
	Stream_Indent( self->debug );

	#if CAUTIOUS
	/* These guys are used to check that we never try and send the same cell of particles twice */
	cellsClearedForTransfer = Memory_Alloc_Array( Bool, self->swarm->cellDomainCount,
		"cellsClearedForTransfer" );
	cellsClearedForTransferDests = Memory_Alloc_Array( Neighbour_Index, self->swarm->cellDomainCount,
		"cellsClearedForTransferDests" );
	for ( dCell_I = 0; dCell_I < self->swarm->cellDomainCount; dCell_I++ ) {
		cellsClearedForTransfer[dCell_I] = False;
		cellsClearedForTransferDests[dCell_I] = 0;
	}
	#endif

	Journal_DPrintfL( self->debug, 1, "Sending the shadow particles going to my %d neighbours:\n", nbrCount );

	for ( nbr_I=0; nbr_I < nbrCount; nbr_I++ ) {
		tParticle_I=0; /* Reset index for new neighbour processor */

		if ( self->shadowParticlesLeavingMeTotalCounts[nbr_I] != 0 ) {
			#ifdef CAUTIOUS
			Neighbour_Index      prevNbr_I;
			#endif

			shadowCellsToProcCount = cellShadowInfo->procShadowCnt[nbr_I];	
			proc_I = procNbrInfo->procNbrTbl[nbr_I];
			Journal_DPrintfL( self->debug, 3, "nbr %d (proc %d) - %d shadow cells going to it:\n",
				nbr_I, proc_I, shadowCellsToProcCount ); 
			#ifdef CAUTIOUS
			for ( prevNbr_I=0; prevNbr_I < nbr_I; prevNbr_I++ ) {
				Journal_Firewall( proc_I != procNbrInfo->procNbrTbl[prevNbr_I], errorStream,
					"Error - in %s(), on proc %u: found in our Swarm's ProcNbrInfo "
					"that our nbr %u is proc %u, but we already sent particles to "
					"that processor as nbr %u! Can't send same particles to same "
					"proc twice.\n", __func__, self->swarm->myRank, nbr_I, proc_I,
					prevNbr_I );
			}
			#endif

			Stream_Indent( self->debug );


			for ( stCell_I=0; stCell_I < shadowCellsToProcCount; stCell_I++ ) {
				currCellParticleCount = self->shadowParticlesLeavingMeCountsPerCell[nbr_I][stCell_I];
				dCell_I = cellShadowInfo->procShadowTbl[nbr_I][stCell_I];
				#ifdef CAUTIOUS
				Journal_Firewall( cellsClearedForTransfer[dCell_I] == False, errorStream,
					"Error - in %s(), on proc %u: while trying to send shadow particles to "
					"nbr %u (proc %u), tried to copy particles from domain cell %u, but "
					"this cell has already had all its particles cleared for send to "
					"nbr %u (proc %u).\n", __func__, self->swarm->myRank, nbr_I, proc_I,
					dCell_I, cellsClearedForTransferDests[dCell_I],
					procNbrInfo->procNbrTbl[cellsClearedForTransferDests[dCell_I]] );
				#endif

				Journal_DPrintfL( self->debug, 3, "Processing Cell %d (%d particles):\n", dCell_I,
					currCellParticleCount );
			
				Stream_Indent( self->debug );

				for ( cParticle_I=0; cParticle_I < currCellParticleCount; cParticle_I++ ) {

					lParticle_I = self->swarm->cellParticleTbl[dCell_I][cParticle_I];
					#if DEBUG
					currParticle = (GlobalParticle*)Swarm_ParticleAt( self->swarm, lParticle_I );
					Journal_DPrintfL( self->debug, 3, "Copying PIC %d, particle %d at "
						"(%.2f,%.2f,%.2g) to shadowParticlesLeavingMe[%d][%d]\n",
						cParticle_I, lParticle_I,
						currParticle->coord[0], currParticle->coord[1], currParticle->coord[2],
						nbr_I, tParticle_I );
					#endif	
					Swarm_CopyParticleOffSwarm( self->swarm,
						self->shadowParticlesLeavingMe[nbr_I], tParticle_I++,
						lParticle_I );
					/* Note: we have no guarantee that the local particle index of where these
					shadow cells are leaving from is monotonically increasing: thus do an insertion
					at the right place. */
					{
						Index		insertionIndex = 0;
						Particle_Index*	currInsertionPtr;
						for ( ;insertionIndex < leavingParticle_I; insertionIndex++ ) {
							currInsertionPtr = &self->shadowParticlesLeavingMeIndices[insertionIndex];
							if (lParticle_I < (*currInsertionPtr) )
							{
								memmove(
									(Pointer)(((ArithPointer)currInsertionPtr) + sizeof(Particle_Index)), 
									currInsertionPtr,
									(leavingParticle_I - insertionIndex) * sizeof(Particle_Index) );

								(*currInsertionPtr) = lParticle_I;
								break;
							}
						}
						if ( insertionIndex == leavingParticle_I) {
							self->shadowParticlesLeavingMeIndices[leavingParticle_I] = lParticle_I;
						}
						leavingParticle_I++;
					}

				}
				Stream_UnIndent( self->debug );

				#ifdef CAUTIOUS
				cellsClearedForTransfer[dCell_I] = True;
				cellsClearedForTransferDests[dCell_I] = nbr_I;
				#endif
				/* Remember to clear the entries for that cell now. */
				self->swarm->cellParticleCountTbl[dCell_I] = 0;
				self->swarm->cellParticleSizeTbl[dCell_I] = 0;
				if ( self->swarm->cellParticleTbl[dCell_I] ) {
					Memory_Free( self->swarm->cellParticleTbl[dCell_I] );
				}
				self->swarm->cellParticleTbl[dCell_I] = NULL;

			}
			Stream_UnIndent( self->debug );

			/* non blocking send out particles */
			MPI_Issend( self->shadowParticlesLeavingMe[nbr_I],
				self->shadowParticlesLeavingMeTotalCounts[nbr_I] * self->swarm->particleExtensionMgr->finalSize,
				MPI_BYTE, proc_I, SHADOW_PARTICLES, self->swarm->comm,
				self->shadowParticlesLeavingMeHandles[nbr_I] );
		}
	}
	#if CAUTIOUS
	Memory_Free( cellsClearedForTransfer );
	Memory_Free( cellsClearedForTransferDests );
	#endif
	Stream_UnIndent( self->debug );
}

/* +++ Statistics Printing +++ */


void _ParticleCommHandler_PrintCommunicationVolumeStats( ParticleCommHandler* self, double startTime, Stream* stream )
{
	double                  myProcTime = 0;
	double*                 procTimes = NULL;
	double                  maxProcTime = 0;
	Particle_Index          totalParticlesRecvdViaShadowFromNbrs = 0;
	ShadowInfo*             cellShadowInfo = CellLayout_GetShadowInfo( self->swarm->cellLayout );
	Processor_Index         proc_I = 0;
	ProcNbrInfo*            procNbrInfo = cellShadowInfo->procNbrInfo;
	Neighbour_Index         nbrCount = procNbrInfo->procNbrCnt;

	myProcTime = MPI_Wtime() - startTime;


	totalParticlesRecvdViaShadowFromNbrs = 0;

	if ( self->swarm->cellShadowCount > 0 ) {
		ShadowInfo*	        cellShadowInfo = CellLayout_GetShadowInfo( self->swarm->cellLayout );
		ProcNbrInfo*		procNbrInfo = cellShadowInfo->procNbrInfo;
		Neighbour_Index		nbr_I;

		for ( nbr_I=0; nbr_I < procNbrInfo->procNbrCnt; nbr_I++ ) {
			totalParticlesRecvdViaShadowFromNbrs += 
				self->particlesArrivingFromNbrShadowCellsTotalCounts[nbr_I];
		}	
	}

	for ( proc_I = 0; proc_I < self->swarm->nProc; proc_I++ ) {
		MPI_Barrier( self->swarm->comm );

		if ( self->swarm->myRank == proc_I ) {			
			Journal_Printf( stream, "...proc %d finished particle communication:\n", self->swarm->myRank );
			Stream_Indent( stream );
			Journal_Printf( stream, "- Particle comm totals via shadow cells (%d nbr procs):"
				" sent %d, recvd %d\n",
				nbrCount, self->shadowParticlesLeavingMeTotalCount,
				totalParticlesRecvdViaShadowFromNbrs );
			Journal_Printf( stream, "- time taken = %.2f (secs)\n", myProcTime );
			Stream_UnIndent( stream );
		}
	}

	procTimes = Memory_Alloc_Array( double, self->swarm->nProc, "procTimes" );
	MPI_Gather( &myProcTime, 1, MPI_DOUBLE, procTimes, 1, MPI_DOUBLE, 0, self->swarm->comm );
	if (self->swarm->myRank == 0 ) {
		for ( proc_I = 0; proc_I < self->swarm->nProc; proc_I++ ) {
			if ( procTimes[proc_I] > maxProcTime ) {
				maxProcTime = procTimes[proc_I];
			}
		}
		// TODO: print some stats on max particles sent/recvd and total sent/recvd
		Journal_Printf( stream, "...Max Communication time by any proc was %.2f (secs)\n", maxProcTime );
	}
	Memory_Free( procTimes );
}

