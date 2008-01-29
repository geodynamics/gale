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
** $Id: Decomp_Sync.c 3584 2006-05-16 11:11:07Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <mpi.h>

#include <StGermain/StGermain.h>
#include "Mesh.h"


/* Textual name of this class */
const Type Decomp_Sync_Type = "Decomp_Sync";


/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/

Decomp_Sync* Decomp_Sync_New( Name name ) {
	return _Decomp_Sync_New( sizeof(Decomp_Sync), 
				 Decomp_Sync_Type, 
				 _Decomp_Sync_Delete, 
				 _Decomp_Sync_Print, 
				 NULL );
}

Decomp_Sync* _Decomp_Sync_New( DECOMP_SYNC_DEFARGS ) {
	Decomp_Sync* self;
	
	/* Allocate memory */
	assert( sizeOfSelf >= sizeof(Decomp_Sync) );
	self = (Decomp_Sync*)_Stg_Class_New( STG_CLASS_PASSARGS );

	/* Virtual info */

	/* Decomp_Sync info */
	_Decomp_Sync_Init( self );

	return self;
}

void _Decomp_Sync_Init( Decomp_Sync* self ) {
	self->decomp = NULL;
	self->commTopo = CommTopology_New();
	Stg_Class_AddRef( self->commTopo );

	self->nDomains = 0;
	self->nRemotes = 0;
	self->remotes = NULL;
	self->nShared = 0;
	self->shared = NULL;
	self->nSharers = NULL;
	self->sharers = NULL;
	self->owners = NULL;

	self->grMap = UIntMap_New();
	self->dsMap = UIntMap_New();

	self->netSrcs = 0;
	self->nSrcs = NULL;
	self->srcs = NULL;
	self->netSnks = 0;
	self->nSnks = NULL;
	self->snks = NULL;

	self->arrays = List_New();
	List_SetItemSize( self->arrays, sizeof(Decomp_Sync_Array*) );
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _Decomp_Sync_Delete( void* sync ) {
	Decomp_Sync*	self = (Decomp_Sync*)sync;

	Decomp_Sync_Destruct( self );
	FreeObject( self->arrays );

	/* Delete the parent. */
	_Stg_Class_Delete( self );
}

void _Decomp_Sync_Print( void* sync, Stream* stream ) {
	Decomp_Sync*	self = (Decomp_Sync*)sync;
	
	/* Set the Journal for printing informations */
	Stream* syncStream;
	syncStream = Journal_Register( InfoStream_Type, "Decomp_SyncStream" );

	/* Print parent */
	Journal_Printf( stream, "Decomp_Sync (ptr): (%p)\n", self );
	_Stg_Class_Print( self, stream );
}


/*--------------------------------------------------------------------------------------------------------------------------
** Public Functions
*/

void Decomp_Sync_SetDecomp( void* sync, Decomp* decomp ) {
	Decomp_Sync*	self = (Decomp_Sync*)sync;

	assert( self );
	assert( Decomp_Sync_ValidateComms( self ) );

	Decomp_Sync_DestructDecomp( self );
	self->decomp = decomp;
	assert( Decomp_Sync_ValidateComms( self ) );

	if( decomp ) {
		self->nDomains = Decomp_GetLocalSize( decomp );
		List_Append( Decomp_GetSyncList( decomp ), &self );
	}
	Decomp_Sync_InitArrays( self );
}

void Decomp_Sync_SetCommTopology( void* sync, CommTopology* commTopo ) {
	Decomp_Sync*	self = (Decomp_Sync*)sync;

	assert( self );

	Decomp_Sync_DestructComm( self );
	if( commTopo ) {
		self->commTopo = commTopo;
		Stg_Class_AddRef( commTopo );
	}
	else {
		self->commTopo = CommTopology_New();
		Stg_Class_AddRef( commTopo );
		if( self->decomp )
			CommTopology_SetComm( self->commTopo, Decomp_GetComm( self->decomp ) );
	}
	assert( Decomp_Sync_ValidateComms( self ) );

	Decomp_Sync_InitArrays( self );
}

void Decomp_Sync_AddRemoteRanks( void* sync, unsigned nRanks, unsigned* ranks ) {
	Decomp_Sync*	self = (Decomp_Sync*)sync;

	assert( self );
	assert( !nRanks || ranks );

	CommTopology_AddIncidence( self->commTopo, nRanks, ranks );
	Decomp_Sync_ExpandArrays( self, nRanks );
}

void Decomp_Sync_SetRemotes( void* sync, unsigned nRemotes, unsigned* remotes ) {
	Decomp_Sync*	self = (Decomp_Sync*)sync;

	assert( self );

	Decomp_Sync_DestructRemotes( self );
	Decomp_Sync_AddRemotes( self, nRemotes, remotes );
}

void Decomp_Sync_AddRemotes( void* sync, unsigned nRemotes, unsigned* remotes ) {
	Decomp_Sync*	self = (Decomp_Sync*)sync;
	unsigned	nIncRanks, *incRanks, *ranks;
	unsigned	*nSrcRemotes, **srcRemotes;
	unsigned	*nSnkRemotes, **snkRemotes;
	unsigned	p_i;

	assert( self );
	assert( self->decomp );
	assert( !nRemotes || remotes );
	assert( Decomp_Sync_ValidateRemotes( self, nRemotes, remotes ) );

	/* If we're not connected to anything, just exit. */
	CommTopology_GetIncidence( self->commTopo, &nIncRanks, &incRanks );
	if( !nIncRanks )
		return;

	/* Prepare source and sink arrays. */
	nSrcRemotes = AllocArray( unsigned, nIncRanks );
	nSnkRemotes = AllocArray( unsigned, nIncRanks );
	srcRemotes = AllocArray( unsigned*, nIncRanks );
	snkRemotes = AllocArray( unsigned*, nIncRanks );

	/* Split these remotes for each processor. */
	Decomp_Sync_SplitRemotes( self, nRemotes, remotes, 
				  nSrcRemotes, srcRemotes, 
				  nSnkRemotes, snkRemotes );

	/* Add them. */
	ranks = AllocArray( unsigned, nIncRanks );
	for( p_i = 0; p_i < nIncRanks; p_i++ )
		ranks[p_i] = p_i;
	Decomp_Sync_AddSources( self, nIncRanks, ranks, nSrcRemotes, srcRemotes );
	Decomp_Sync_AddSinks( self, nIncRanks, ranks, nSnkRemotes, snkRemotes );
	Decomp_Sync_BuildShared( self );
	FreeArray( ranks );

	/* Free the space. */
	FreeArray( nSrcRemotes );
	FreeArray( nSnkRemotes );
	FreeArray2D( nIncRanks, srcRemotes );
	FreeArray2D( nIncRanks, snkRemotes );
}

void Decomp_Sync_SetSources( void* sync, unsigned nRanks, unsigned* ranks, 
			     unsigned* nSources, unsigned** sources )
{
	Decomp_Sync*	self = (Decomp_Sync*)sync;
	unsigned	p_i, s_i;

	assert( nRanks <= CommTopology_GetIncidenceSize( self->commTopo ) );
	assert( !nRanks || ranks );

	for( p_i = 0; p_i < nRanks; p_i++ ) {
		assert( ranks[p_i] < CommTopology_GetIncidenceSize( self->commTopo ) );
		for( s_i = 0; s_i < self->nSrcs[ranks[p_i]]; s_i++ )
			UIntMap_Remove( self->grMap, Decomp_Sync_DomainToGlobal( self, self->srcs[ranks[p_i]][s_i] ) );
		KillArray( self->srcs[ranks[p_i]] );
		self->netSrcs -= self->nSrcs[ranks[p_i]];
		self->nDomains -= self->nSrcs[ranks[p_i]];
		self->nRemotes -= self->nSrcs[ranks[p_i]];
		self->nSrcs[ranks[p_i]] = 0;
	}

	Decomp_Sync_AddSources( self, nRanks, ranks, nSources, sources );
}

void Decomp_Sync_SetSinks( void* sync, unsigned nRanks, unsigned* ranks, 
			   unsigned* nSinks, unsigned** sinks )
{
	Decomp_Sync*	self = (Decomp_Sync*)sync;
	unsigned	p_i;

	assert( nRanks <= CommTopology_GetIncidenceSize( self->commTopo ) );
	assert( !nRanks || ranks );

	for( p_i = 0; p_i < nRanks; p_i++ ) {
		assert( ranks[p_i] < CommTopology_GetIncidenceSize( self->commTopo ) );
		KillArray( self->snks[ranks[p_i]] );
		self->netSnks -= self->nSnks[ranks[p_i]];
		self->nSnks[ranks[p_i]] = 0;
	}

	Decomp_Sync_AddSinks( self, nRanks, ranks, nSinks, sinks );
}

void Decomp_Sync_AddSources( void* sync, unsigned nRanks, unsigned* ranks, 
			     unsigned* nSources, unsigned** sources )
{
	Decomp_Sync*	self = (Decomp_Sync*)sync;
	unsigned	nNewRemotes;
	unsigned	domainInd;
	unsigned	offs;
	unsigned	p_i;

	assert( self );
	assert( self->decomp );
	assert( nRanks <= CommTopology_GetIncidenceSize( self->commTopo ) );
	assert( !nRanks || ranks );
	assert( nSources && sources );

	nNewRemotes = 0;
	for( p_i = 0; p_i < nRanks; p_i++ ) {
		assert( ranks[p_i] < CommTopology_GetIncidenceSize( self->commTopo ) );
		assert( Decomp_Sync_ValidateRemotes( self, nSources[p_i], sources[p_i] ) );
		nNewRemotes += nSources[p_i];
	}

	if( nNewRemotes ) {
		self->remotes = ReallocNamedArray( self->remotes, unsigned, self->nRemotes + nNewRemotes, 
						   "Decomp_Sync::remotes" );
		self->owners = ReallocNamedArray( self->owners, unsigned, self->nRemotes + nNewRemotes, 
						  "Decomp_Sync::owners" );
		offs = self->nRemotes;
		for( p_i = 0; p_i < nRanks; p_i++ ) {
			unsigned	r_i;

			for( r_i = 0; r_i < nSources[p_i]; r_i++ ) {
				self->remotes[offs] = sources[p_i][r_i];
				UIntMap_Insert( self->grMap, sources[p_i][r_i], offs );
				self->owners[offs++] = ranks[p_i];
			}
		}
		self->nRemotes += nNewRemotes;
		self->nDomains += nNewRemotes;
	}

	for( p_i = 0; p_i < nRanks; p_i++ ) {
		unsigned	rank = ranks[p_i];
		unsigned	s_i;

		self->srcs[rank] = ReallocArray( self->srcs[rank], unsigned, self->nSrcs[rank] + nSources[p_i] );
		for( s_i = 0; s_i < nSources[p_i]; s_i++ ) {
			insist( Decomp_Sync_GlobalToDomain( self, sources[p_i][s_i], &domainInd ), == True );
			self->srcs[rank][self->nSrcs[rank] + s_i] = domainInd;
		}
		self->nSrcs[rank] += nSources[p_i];
		self->netSrcs += nSources[p_i];
	}
}

void Decomp_Sync_AddSinks( void* sync, unsigned nRanks, unsigned* ranks, 
			   unsigned* nSinks, unsigned** sinks )
{
	Decomp_Sync*	self = (Decomp_Sync*)sync;
	unsigned	domainInd;
	unsigned	p_i;

	assert( self );
	assert( self->decomp );
	assert( nRanks <= CommTopology_GetIncidenceSize( self->commTopo ) );
	assert( !nRanks || ranks );
	assert( nSinks && sinks );

	for( p_i = 0; p_i < nRanks; p_i++ ) {
		unsigned	rank = ranks[p_i];
		unsigned	s_i;

		assert( rank < CommTopology_GetIncidenceSize( self->commTopo ) );
		assert( Decomp_Sync_ValidateSinks( self, nSinks[p_i], sinks[p_i] ) );
		self->snks[rank] = ReallocArray( self->snks[rank], unsigned, self->nSnks[rank] + nSinks[p_i] );
		for( s_i = 0; s_i < nSinks[p_i]; s_i++ ) {
			insist( Decomp_Sync_GlobalToDomain( self, sinks[p_i][s_i], &domainInd ), == True );
			self->snks[rank][self->nSnks[rank] + s_i] = domainInd;
		}
		self->nSnks[rank] += nSinks[p_i];
		self->netSnks += nSinks[p_i];
	}
}

void Decomp_Sync_SetRequired( void* sync, unsigned nRequired, unsigned* required ) {
	Decomp_Sync*	self = (Decomp_Sync*)sync;

	assert( self );

	Decomp_Sync_DestructRemotes( self );
	Decomp_Sync_AddRequired( self, nRequired, required );
}

void Decomp_Sync_AddRequired( void* sync, unsigned nRequired, unsigned* required ) {
	Decomp_Sync*	self = (Decomp_Sync*)sync;
	MPI_Comm	comm;
	unsigned	nProcs, rank;
	RangeSet	*reqSet, *remReqSet, *lSet;
	unsigned	nLocals, *locals;
	unsigned	nOldIncRanks, nIncRanks;
	unsigned	nBytes;
	Stg_Byte*	bytes;
	unsigned	tag = 8849;
	unsigned	p_i;

	assert( self && Stg_CheckType( self, Decomp_Sync ) );
	assert( Decomp_Sync_ValidateRemotes( self, nRequired, required ) );

	comm = CommTopology_GetComm( self->commTopo );
	MPI_Comm_rank( comm, (int*)&rank );
	MPI_Comm_size( comm, (int*)&nProcs );

	reqSet = RangeSet_New();
	RangeSet_SetIndices( reqSet, nRequired, required );

	Decomp_GetLocals( self->decomp, &nLocals, &locals );
	lSet = RangeSet_New();
	RangeSet_SetIndices( lSet, nLocals, locals );

	nOldIncRanks = nIncRanks;
	remReqSet = RangeSet_New();
	for( p_i = 0; p_i < nProcs; p_i++ ) {
		Stg_Byte	state;
		unsigned	nInds, *inds;
		unsigned	localRank;

		if( rank == p_i )
			RangeSet_Pickle( reqSet, &nBytes, &bytes );
		MPI_Bcast( &nBytes, 1, MPI_UNSIGNED, p_i, comm );
		if( rank != p_i )
			bytes = AllocArray( Stg_Byte, nBytes );
		MPI_Bcast( bytes, nBytes, MPI_BYTE, p_i, comm );

		if( rank != p_i ) {
			RangeSet_Unpickle( remReqSet, nBytes, bytes );
			FreeArray( bytes );

			RangeSet_Intersection( remReqSet, lSet );
			if( RangeSet_GetSize( remReqSet ) ) {
				state = 1;
				MPI_Gather( &state, 1, MPI_BYTE, NULL, 1, MPI_BYTE, p_i, comm );

				RangeSet_Pickle( remReqSet, &nBytes, &bytes );
				MPI_Send( &nBytes, 1, MPI_UNSIGNED, p_i, tag, comm );
				MPI_Send( bytes, nBytes, MPI_BYTE, p_i, tag, comm );
				FreeArray( bytes );

				inds = NULL;
				RangeSet_GetIndices( remReqSet, &nInds, &inds );
				if( !CommTopology_GlobalToLocal( self->commTopo, p_i, &localRank ) ) {
					Decomp_Sync_AddRemoteRanks( self, 1, &p_i );
					insist( CommTopology_GlobalToLocal( self->commTopo, p_i, &localRank ), == True );
					nIncRanks = CommTopology_GetIncidenceSize( self->commTopo );
				}
				Decomp_Sync_AddSinks( sync, 1, &localRank, &nInds, &inds );
				FreeArray( inds );
			}
			else {
				state = 0;
				MPI_Gather( &state, 1, MPI_BYTE, NULL, 1, MPI_BYTE, p_i, comm );
			}
		}
		else {
			Stg_Byte*	states;
			unsigned	nActive, *active;
			MPI_Status	status;
			unsigned	nNewRanks, *newRanks;
			unsigned	p_j;

			FreeArray( bytes );

			state = 0;
			states = AllocArray( Stg_Byte, nProcs );
			MPI_Gather( &state, 1, MPI_BYTE, states, 1, MPI_BYTE, p_i, comm );
			nActive = 0;
			nNewRanks = 0;
			for( p_j = 0; p_j < nProcs; p_j++ ) {
				if( states[p_j] ) {
					nActive++;
					if( !CommTopology_GlobalToLocal( self->commTopo, p_j, &localRank ) )
						nNewRanks++;
				}
			}
			active = AllocArray( unsigned, nActive );
			newRanks = AllocArray( unsigned, nNewRanks );
			nActive = 0;
			nNewRanks = 0;
			for( p_j = 0; p_j < nProcs; p_j++ ) {
				if( states[p_j] ) {
					active[nActive++] = p_j;
					if( !CommTopology_GlobalToLocal( self->commTopo, p_j, &localRank ) )
						newRanks[nNewRanks++] = p_j;
				}
			}
			FreeArray( states );

			Decomp_Sync_AddRemoteRanks( sync, nNewRanks, newRanks );
			if( nNewRanks )
				nIncRanks = CommTopology_GetIncidenceSize( self->commTopo );
			FreeArray( newRanks );

			for( p_j = 0; p_j < nActive; p_j++ ) {
				MPI_Recv( &nBytes, 1, MPI_UNSIGNED, active[p_j], tag, comm, &status );
				bytes = AllocArray( Stg_Byte, nBytes );
				MPI_Recv( bytes, nBytes, MPI_BYTE, active[p_j], tag, comm, &status );
				RangeSet_Unpickle( reqSet, nBytes, bytes );
				FreeArray( bytes );
				inds = NULL;
				RangeSet_GetIndices( reqSet, &nInds, &inds );
				insist( CommTopology_GlobalToLocal( self->commTopo, active[p_j], &localRank ), == True );
				Decomp_Sync_AddSources( sync, 1, &localRank, &nInds, &inds );
				FreeArray( inds );
			}

			FreeArray( active );
		}
	}

	FreeObject( reqSet );
	FreeObject( remReqSet );
	FreeObject( lSet );

	Decomp_Sync_BuildShared( self );
}

void Decomp_Sync_BuildShared( void* sync ) {
	Decomp_Sync*	self = (Decomp_Sync*)sync;
	unsigned	nIncRanks;
	RangeSet	*dSet, **rSets, *allSet;
	unsigned	nLocals, *locals;
	unsigned	nSendBytes, *nRecvBytes;
	Stg_Byte	*sendBytes, **recvBytes;
	unsigned	**sharers;
	unsigned	nInds, *inds;
	unsigned	r_i, s_i, ind_i;

	assert( self && Stg_CheckType( self, Decomp_Sync ) );

	/* Free shared info. */
	FreeArray( self->shared );
	FreeArray( self->nSharers );
	FreeArray( self->sharers );

	nIncRanks = CommTopology_GetIncidenceSize( self->commTopo );

	/* Communicate a range set of my domain indices to all neighbours. */
	dSet = RangeSet_New();
	RangeSet_SetIndices( dSet, self->nRemotes, self->remotes );
	Decomp_GetLocals( self->decomp, &nLocals, &locals );
	RangeSet_AddIndices( dSet, nLocals, locals );
	RangeSet_Pickle( dSet, &nSendBytes, &sendBytes );
	CommTopology_Allgather( self->commTopo, 
				nSendBytes, sendBytes, 
				&nRecvBytes, &recvBytes, 
				sizeof(Stg_Byte) );
	FreeArray( sendBytes );
	rSets = AllocArray( RangeSet*, nIncRanks );
	allSet = RangeSet_New();
	for( r_i = 0; r_i < nIncRanks; r_i++ ) {
		rSets[r_i] = RangeSet_New();
		RangeSet_Unpickle( rSets[r_i], nRecvBytes[r_i], recvBytes[r_i] );
		RangeSet_Intersection( rSets[r_i], dSet );
		RangeSet_Union( allSet, rSets[r_i] );
		FreeArray( recvBytes[r_i] );
	}
	FreeArray( recvBytes );
	FreeArray( nRecvBytes );

	/* Use the 'allSet' to generate the shared indices. */
	self->shared = NULL;
	RangeSet_GetIndices( allSet, &self->nShared, &self->shared );
	FreeObject( allSet );
	for( s_i = 0; s_i < self->nShared; s_i++ ) {
		insist( Decomp_Sync_GlobalToDomain( self, self->shared[s_i], self->shared + s_i ), == True );
		UIntMap_Insert( self->dsMap, self->shared[s_i], s_i );
	}

	/* Now that we have each intersection, convert to shared arrays. */
	self->nSharers = AllocArray( unsigned, self->nShared );
	memset( self->nSharers, 0, self->nShared * sizeof(unsigned) );
	sharers = AllocArray2D( unsigned, self->nShared, nIncRanks );
	for( r_i = 0; r_i < nIncRanks; r_i++ ) {
		inds = NULL;
		RangeSet_GetIndices( rSets[r_i], &nInds, &inds );
		FreeObject( rSets[r_i] );
		for( ind_i = 0; ind_i < nInds; ind_i++ ) {
			insist( Decomp_Sync_GlobalToDomain( self, inds[ind_i], inds + ind_i ), == True );
			insist( Decomp_Sync_DomainToShared( self, inds[ind_i], inds + ind_i ), == True );
			sharers[inds[ind_i]][self->nSharers[inds[ind_i]]++] = r_i;
		}
		FreeArray( inds );
	}
	FreeArray( rSets );

	/* Store final array. */
	self->sharers = AllocComplex2D( unsigned, self->nShared, self->nSharers );
	for( s_i = 0; s_i < self->nShared; s_i++ )
		memcpy( self->sharers[s_i], sharers[s_i], self->nSharers[s_i] * sizeof(unsigned) );
	FreeArray( sharers );

#if 0
	/* Create a range set of all sinks and sources. */
	nIncRanks = CommTopology_GetIncidenceSize( self->commTopo );
	shareSet = RangeSet_New();
	tmpSet = RangeSet_New();
	for( p_i = 0; p_i < nIncRanks; p_i++ ) {
		RangeSet_SetIndices( tmpSet, self->nSnks[p_i], self->snks[p_i] );
		RangeSet_Union( shareSet, tmpSet );
		RangeSet_SetIndices( tmpSet, self->nSrcs[p_i], self->srcs[p_i] );
		RangeSet_Union( shareSet, tmpSet );
	}
	FreeObject( tmpSet );

	/* These indices are the shared elements. */
	self->shared = NULL;
	RangeSet_GetIndices( shareSet, &self->nShared, &self->shared );
	FreeObject( shareSet );

	/* If there are no shared indices, exit now. */
	if( !self->nShared )
		return;

	/* Build the shared mapping. */
	for( s_i = 0; s_i < self->nShared; s_i++ )
		UIntMap_Insert( self->dsMap, self->shared[s_i], s_i );

	/* Create temporary storage for the sharers. */
	self->nSharers = AllocNamedArray( unsigned, self->nShared, "Decomp_Sync::nSharers" );
	sharers = Memory_Alloc_2DArray_Unnamed( unsigned, self->nShared, nIncRanks );
	memset( self->nSharers, 0, self->nShared * sizeof(unsigned) );

	/* Build the sharer lists. */
	for( p_i = 0; p_i < nIncRanks; p_i++ ) {
		for( s_i = 0; s_i < self->nSnks[p_i]; s_i++ ) {
			unsigned	sharedInd;
			unsigned	curSharer;

			insist( UIntMap_Map( self->dsMap, self->snks[p_i][s_i], &sharedInd ), == True );
			curSharer = self->nSharers[sharedInd]++;
			sharers[sharedInd][curSharer] = p_i;
		}
	}

	/* Transfer to final storage. */
	self->sharers = Memory_Alloc_2DComplex( unsigned, self->nShared, self->nSharers, "Decomp_Sync::sharers" );
	for( s_i = 0; s_i < self->nShared; s_i++ )
		memcpy( self->sharers[s_i], sharers[s_i], self->nSharers[s_i] * sizeof(unsigned) );

	/* Free the old sharers array. */
	FreeArray( sharers );
#endif
}

unsigned Decomp_Sync_GetGlobalSize( void* sync ) {
	Decomp_Sync*	self = (Decomp_Sync*)sync;

	assert( self );
	assert( self->decomp );

	return Decomp_GetGlobalSize( self->decomp );
}

unsigned Decomp_Sync_GetLocalSize( void* sync ) {
	Decomp_Sync*	self = (Decomp_Sync*)sync;

	assert( self );
	assert( self->decomp );

	return Decomp_GetLocalSize( self->decomp );
}

unsigned Decomp_Sync_GetRemoteSize( void* sync ) {
	Decomp_Sync*	self = (Decomp_Sync*)sync;

	assert( self );

	return self->nRemotes;
}

unsigned Decomp_Sync_GetDomainSize( void* sync ) {
	Decomp_Sync*	self = (Decomp_Sync*)sync;

	assert( self );

	return self->nDomains;
}

unsigned Decomp_Sync_GetSharedSize( void* sync ) {
	Decomp_Sync*	self = (Decomp_Sync*)sync;

	assert( self );

	return self->nShared;
}

void Decomp_Sync_GetRemotes( void* sync, unsigned* nRemotes, unsigned** remotes ) {
	Decomp_Sync*	self = (Decomp_Sync*)sync;

	assert( self );
	assert( nRemotes );
	assert( remotes );

	*nRemotes = self->nRemotes;
	*remotes = self->remotes;
}

Bool Decomp_Sync_GlobalToDomain( void* sync, unsigned global, unsigned* domain ) {
	Decomp_Sync*	self = (Decomp_Sync*)sync;

	assert( self );
	assert( self->decomp );
	assert( global < self->decomp->nGlobals );
	assert( self->decomp->glMap );
	assert( self->grMap );
	assert( domain );

	if( Decomp_GlobalToLocal( self->decomp, global, domain ) )
		return True;
	else if( UIntMap_Map( self->grMap, global, domain ) ) {
		*domain += Decomp_GetLocalSize( self->decomp );
		return True;
	}

	return False;
}

unsigned Decomp_Sync_DomainToGlobal( void* sync, unsigned domain ) {
	Decomp_Sync*	self = (Decomp_Sync*)sync;

	assert( self );
	assert( self->decomp );
	assert( domain < self->decomp->nLocals + self->nRemotes );

	if( domain < self->decomp->nLocals )
		return Decomp_LocalToGlobal( self->decomp, domain );
	else
		return self->remotes[domain - Decomp_GetLocalSize( self->decomp )];
}

Bool Decomp_Sync_DomainToShared( void* sync, unsigned domain, unsigned* shared ) {
	Decomp_Sync*	self = (Decomp_Sync*)sync;

	assert( self );
	assert( self->decomp );
	assert( domain < self->nDomains );
	assert( self->dsMap );
	assert( shared );

	return UIntMap_Map( self->dsMap, domain, shared );
}

unsigned Decomp_Sync_SharedToDomain( void* sync, unsigned shared ) {
	Decomp_Sync*	self = (Decomp_Sync*)sync;

	assert( self );
	assert( shared < self->nShared );
	assert( self->shared );

	return self->shared[shared];
}

unsigned Decomp_Sync_GetOwner( void* sync, unsigned remote ) {
	Decomp_Sync*	self = (Decomp_Sync*)sync;

	assert( self );
	assert( remote < self->nRemotes );
	assert( self->owners );

	return self->owners[remote];
}

void Decomp_Sync_GetSharers( void* sync, unsigned shared, 
			     unsigned* nSharers, unsigned** sharers )
{
	Decomp_Sync*	self = (Decomp_Sync*)sync;

	assert( self );
	assert( shared < self->nShared );
	assert( self->nSharers );
	assert( self->sharers );
	assert( nSharers );
	assert( sharers );

	*nSharers = self->nSharers[shared];
	*sharers = self->sharers[shared];
}

void Decomp_Sync_GetSources( void* sync, unsigned rank, unsigned* nSources, unsigned** sources ) {
	Decomp_Sync*	self = (Decomp_Sync*)sync;

	assert( self );
	assert( rank < CommTopology_GetIncidenceSize( self->commTopo ) );
	assert( nSources && sources );
	assert( self->nSrcs && self->srcs );

	*nSources = self->nSrcs[rank];
	*sources = self->srcs[rank];
}

void Decomp_Sync_GetSinks( void* sync, unsigned rank, unsigned* nSinks, unsigned** sinks ) {
	Decomp_Sync*	self = (Decomp_Sync*)sync;

	assert( self );
	assert( rank < CommTopology_GetIncidenceSize( self->commTopo ) );
	assert( nSinks && sinks );
	assert( self->nSnks && self->snks );

	*nSinks = self->nSnks[rank];
	*sinks = self->snks[rank];
}

Decomp* Decomp_Sync_GetDecomp( void* sync ) {
	Decomp_Sync*	self = (Decomp_Sync*)sync;

	assert( self );

	return self->decomp;
}

CommTopology* Decomp_Sync_GetCommTopology( void* sync ) {
	Decomp_Sync*	self = (Decomp_Sync*)sync;

	assert( self );

	return self->commTopo;
}

void Decomp_Sync_SyncArray( void* sync, Decomp_Sync_Array* array ) {
	Decomp_Sync*	self = (Decomp_Sync*)sync;
	unsigned	nInc;
	unsigned*	inc;
	Stg_Byte*	snkArray;
	Stg_Byte*	srcArray;
	unsigned	p_i;

	/* Sanity checks. */
	assert( self );
	assert( self->decomp );
	assert( array );

	/* Pack from locals to a contiguous array. */
	if( self->netSnks > 0 ) {
		unsigned	snk_i;

		snkArray = AllocArray( Stg_Byte, self->netSnks * array->itemSize );
		for( snk_i = 0; snk_i < self->netSnks; snk_i++ ) {
			memcpy( snkArray + snk_i * array->itemSize, 
				(Stg_Byte*)array->snkArray + array->snkOffs[snk_i], 
				array->itemSize );
		}
	}
	else
		snkArray = NULL;

	/* Allocate for sources. */
	srcArray = AllocArray( Stg_Byte, self->netSrcs * array->itemSize );

	/* Get incidence. */
	CommTopology_GetIncidence( self->commTopo, &nInc, &inc );

	/* Transfer. */
	for( p_i = 0; p_i < nInc; p_i++ ) {
		int		snkSize = array->snkSizes[p_i];
		int		snkDisp = array->snkDisps[p_i];
		int		srcSize = array->srcSizes[p_i];
		int		srcDisp = array->srcDisps[p_i];
		MPI_Status	status;
		unsigned	tag = 6669;

		MPI_Sendrecv( snkArray + snkDisp, snkSize, MPI_BYTE, inc[p_i], tag, 
			      srcArray + srcDisp, srcSize, MPI_BYTE, inc[p_i], tag, 
			      CommTopology_GetComm( self->commTopo ), &status );
	}

	/* Free the sink array. */
	FreeArray( snkArray );

	/* Unpack sources. */
	if( self->netSrcs > 0 ) {
		unsigned	src_i;

		for( src_i = 0; src_i < self->netSrcs; src_i++ ) {
			memcpy( (Stg_Byte*)array->srcArray + array->srcOffs[src_i], 
				srcArray + src_i * array->itemSize, 
				array->itemSize );
		}
	}

	/* Free source array. */
	FreeArray( srcArray );
}

void Decomp_Sync_Update( void* sync ) {
	Decomp_Sync*	self = (Decomp_Sync*)sync;

	assert( self );

	Decomp_Sync_DestructRemotes( self );
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Private Functions
*/

void Decomp_Sync_InitArrays( Decomp_Sync* self ) {
	unsigned	nIncRanks;

	assert( self );

	if( self->commTopo ) {
		nIncRanks = CommTopology_GetIncidenceSize( self->commTopo );
		if( nIncRanks ) {
			self->nSrcs = AllocNamedArray( unsigned, nIncRanks, "Decomp_Sync::nSrcs" );
			self->srcs = AllocNamedArray( unsigned*, nIncRanks, "Decomp_Sync::srcs" );
			self->nSnks = AllocNamedArray( unsigned, nIncRanks, "Decomp_Sync::nSnks" );
			self->snks = AllocNamedArray( unsigned*, nIncRanks, "Decomp_Sync::snks" );
			memset( self->nSrcs, 0, nIncRanks * sizeof(unsigned) );
			memset( self->srcs, 0, nIncRanks * sizeof(unsigned*) );
			memset( self->nSnks, 0, nIncRanks * sizeof(unsigned) );
			memset( self->snks, 0, nIncRanks * sizeof(unsigned*) );
		}
	}
}

void Decomp_Sync_ExpandArrays( Decomp_Sync* self, unsigned nNewRanks ) {
	unsigned	nIncRanks, nOldIncRanks;

	assert( self );
	assert( !self->nShared || (self->nSharers && self->sharers) );
	assert( !self->nRemotes || self->owners );

	nIncRanks = CommTopology_GetIncidenceSize( self->commTopo );
	nOldIncRanks = nIncRanks - nNewRanks;
	self->nSrcs = ReallocNamedArray( self->nSrcs, unsigned, nIncRanks, "Decomp_Sync::nSrcs" );
	self->srcs = ReallocNamedArray( self->srcs, unsigned*, nIncRanks, "Decomp_Sync::srcs" );
	self->nSnks = ReallocNamedArray( self->nSnks, unsigned, nIncRanks, "Decomp_Sync::nSnks" );
	self->snks = ReallocNamedArray( self->snks, unsigned*, nIncRanks, "Decomp_Sync::snks" );
	memset( self->nSrcs + nOldIncRanks, 0, nNewRanks * sizeof(unsigned) );
	memset( self->srcs + nOldIncRanks, 0, nNewRanks * sizeof(unsigned*) );
	memset( self->nSnks + nOldIncRanks, 0, nNewRanks * sizeof(unsigned) );
	memset( self->snks + nOldIncRanks, 0, nNewRanks * sizeof(unsigned*) );
}

void Decomp_Sync_SplitRemotes( Decomp_Sync* self, unsigned nRemotes, unsigned* remotes, 
			       unsigned* nSrcs, unsigned** srcs, 
			       unsigned* nSnks, unsigned** snks )
{
	unsigned	nLocals, *locals;
	unsigned	nIncRanks;
	RangeSet	*remSet, *locSet;
	unsigned	nBytes, *nRemBytes, *nFndBytes;
	Stg_Byte	*bytes, **remBytes, **fndBytes;
	unsigned	p_i;

	assert( self );
	assert( self->decomp );

	/* Make a range set of our remotes. */
	remSet = RangeSet_New();
	RangeSet_SetIndices( remSet, nRemotes, remotes );
	RangeSet_Pickle( remSet, &nBytes, &bytes );

	/* Collect neighbouring remotes. */
	CommTopology_Allgather( self->commTopo, 
				nBytes, bytes, 
				&nRemBytes, (void***)&remBytes, 
				sizeof(Stg_Byte) );

	/* Done with the bytes. */
	FreeArray( bytes );

	/* Build a range set of our locals. */
	locSet = RangeSet_New();
	Decomp_GetLocals( self->decomp, &nLocals, &locals );
	RangeSet_SetIndices( locSet, nLocals, locals );

	/* Intersect our locals and our neighbours remotes to build our sink sets. */
	nIncRanks = CommTopology_GetIncidenceSize( self->commTopo );
	remSet = RangeSet_New();
	for( p_i = 0; p_i < nIncRanks; p_i++ ) {
		RangeSet_Unpickle( remSet, nRemBytes[p_i], remBytes[p_i] );

		/* Done with the remote bytes. */
		FreeArray( remBytes[p_i] );

		RangeSet_Intersection( remSet, locSet );
		snks[p_i] = NULL;
		RangeSet_GetIndices( remSet, nSnks + p_i, snks + p_i );
		RangeSet_Pickle( remSet, nRemBytes + p_i, remBytes + p_i );
	}

	/* Local set is no longer needed. */
	FreeObject( locSet );

	/* Return what we've found. */
	CommTopology_Alltoall( self->commTopo, 
			       nRemBytes, (void**)remBytes, 
			       &nFndBytes, (void***)&fndBytes, 
			       sizeof(Stg_Byte) );

	/* Done with the remote bytes arrays. */
	FreeArray( nRemBytes );
	FreeArray2D( nIncRanks, remBytes );

	/* Construct our source arrays with what was found. */
	for( p_i = 0; p_i < nIncRanks; p_i++ ) {
		RangeSet_Unpickle( remSet, nFndBytes[p_i], fndBytes[p_i] );

		/* Empty the internal found arrays. */
		FreeArray( fndBytes[p_i] );

		srcs[p_i] = NULL;
		RangeSet_GetIndices( remSet, nSrcs + p_i, srcs + p_i );
	}

	/* Free everything else. */
	FreeArray( nFndBytes );
	FreeArray( fndBytes );
	FreeObject( remSet );
}

void Decomp_Sync_Destruct( Decomp_Sync* self ) {
	assert( self );

	Decomp_Sync_DestructDecomp( self );
	Decomp_Sync_DestructComm( self );
}

void Decomp_Sync_DestructDecomp( Decomp_Sync* self ) {
	assert( self );

	Decomp_Sync_DestructRemotes( self );
	if( self->decomp ) {
		List_Remove( Decomp_GetSyncList( self->decomp ), &self );
		self->decomp = NULL;
	}
}

void Decomp_Sync_DestructComm( Decomp_Sync* self ) {
	assert( self );

	Decomp_Sync_DestructRemotes( self );
	Stg_Class_RemoveRef( self->commTopo );
	self->commTopo = NULL;
	KillArray( self->nSrcs );
	KillArray( self->srcs );
	KillArray( self->nSnks );
	KillArray( self->snks );
}

void Decomp_Sync_DestructRemotes( Decomp_Sync* self ) {
	assert( self );

	Decomp_Sync_DestructArrays( self );
	Decomp_Sync_DestructSources( self );
	Decomp_Sync_DestructSinks( self );
}

void Decomp_Sync_DestructSources( Decomp_Sync* self ) {
	unsigned	nIncRanks;
	unsigned	p_i;

	assert( self );

	if( self->commTopo )
		nIncRanks = CommTopology_GetIncidenceSize( self->commTopo );

	if( self->decomp )
		self->nDomains = Decomp_GetLocalSize( self->decomp );

	KillArray( self->remotes );
	self->nRemotes = 0;
	UIntMap_Clear( self->grMap );
	self->netSrcs = 0;
	for( p_i = 0; p_i < nIncRanks; p_i++ )
		KillArray( self->srcs[p_i] );
}

void Decomp_Sync_DestructSinks( Decomp_Sync* self ) {
	unsigned	nIncRanks;
	unsigned	p_i;

	if( self->commTopo )
		nIncRanks = CommTopology_GetIncidenceSize( self->commTopo );

	KillArray( self->shared );
	self->nShared = 0;
	KillArray( self->nSharers );
	KillArray( self->sharers );
	UIntMap_Clear( self->dsMap );
	self->netSnks = 0;
	for( p_i = 0; p_i < nIncRanks; p_i++ )
		KillArray( self->snks[p_i] );
}

void Decomp_Sync_DestructArrays( Decomp_Sync* self ) {
	unsigned	a_i;

	for( a_i = 0; a_i < List_GetSize( self->arrays ); a_i++ ) {
		Decomp_Sync_Array*	array = *(Decomp_Sync_Array**)List_GetItem( self->arrays, a_i );

		Decomp_Sync_Array_SetSync( array, self );
	}

	List_Clear( self->arrays );
}

#ifndef NDEBUG
Bool Decomp_Sync_ValidateRemotes( Decomp_Sync* self, unsigned nRemotes, unsigned* remotes ) {
	unsigned	domainInd;
	unsigned	r_i;

	for( r_i = 0; r_i < nRemotes; r_i++ ) {
		if( Decomp_Sync_GlobalToDomain( self, remotes[r_i], &domainInd ) )
			return False;
	}

	return True;
}

Bool Decomp_Sync_ValidateSinks( Decomp_Sync* self, unsigned nSinks, unsigned* sinks ) {
	RangeSet	*lSet, *sSet;
	unsigned	nLocals, *locals;
	unsigned	nInds;

	Decomp_GetLocals( self->decomp, &nLocals, &locals );

	lSet = RangeSet_New();
	sSet = RangeSet_New();
	RangeSet_SetIndices( lSet, nLocals, locals );
	RangeSet_SetIndices( sSet, nSinks, sinks );
	RangeSet_Subtraction( sSet, lSet );
	nInds = RangeSet_GetSize( sSet );
	FreeObject( lSet );
	FreeObject( sSet );

	return !nInds;
}

Bool Decomp_Sync_ValidateComms( Decomp_Sync* self ) {
	assert( self );

	if( self->decomp && self->commTopo )
		return Decomp_GetComm( self->decomp ) == CommTopology_GetComm( self->commTopo );
	else
		return True;
}
#endif
