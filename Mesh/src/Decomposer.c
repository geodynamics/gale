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
** $Id: Decomposer.c 3584 2006-05-16 11:11:07Z PatrickSunter $
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
const Type Decomposer_Type = "Decomposer";


/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/

Decomposer* Decomposer_New( Name name ) {
	return _Decomposer_New( sizeof(Decomposer), 
				Decomposer_Type, 
				_Decomposer_Delete, 
				_Decomposer_Print, 
				NULL, 
				_Decomposer_Decompose );
}

Decomposer* _Decomposer_New( DECOMPOSER_DEFARGS ) {
	Decomposer* self;
	
	/* Allocate memory */
	assert( sizeOfSelf >= sizeof(Decomposer) );
	self = (Decomposer*)_Stg_Class_New( STG_CLASS_PASSARGS );

	/* Virtual info */
	self->decomposeFunc = decomposeFunc;

	/* Decomposer info */
	_Decomposer_Init( self );

	return self;
}

void _Decomposer_Init( Decomposer* self ) {
	self->comm = MPI_COMM_WORLD;
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _Decomposer_Delete( void* decomposer ) {
	Decomposer*	self = (Decomposer*)decomposer;

	/* Delete the parent. */
	_Stg_Class_Delete( self );
}

void _Decomposer_Print( void* decomposer, Stream* stream ) {
	Decomposer*	self = (Decomposer*)decomposer;
	
	/* Set the Journal for printing informations */
	Stream* decomposerStream;
	decomposerStream = Journal_Register( InfoStream_Type, "DecomposerStream" );

	/* Print parent */
	Journal_Printf( stream, "Decomposer (ptr): (%p)\n", self );
	_Stg_Class_Print( self, stream );
}

void _Decomposer_Decompose( void* decomposer, unsigned nDomains, unsigned* domains, 
			    CommTopology** commTopo, Decomp** decomp, Decomp_Sync** sync )
{
	Decomposer*	self = (Decomposer*)decomposer;
	CommTopology*	tmpCommTopo;
	Decomp*		tmpDecomp;
	RangeSet**	isects;

	assert( self );
	assert( sync );

	if( commTopo && *commTopo ) {
		Decomposer_BuildLocalIntersections( self, nDomains, domains, 
						    *commTopo, &isects );
		tmpCommTopo = *commTopo;
	}
	else {
		Decomposer_BuildCommTopology( self, nDomains, domains, 
					      &tmpCommTopo, &isects );

		/* Keep communication topology? */
		if( commTopo )
			*commTopo = tmpCommTopo;
	}

	Decomposer_Claim( self, tmpCommTopo, isects, nDomains, domains, 
			  &tmpDecomp, sync );

	/* Interested in keeping decomp? */
	if( decomp )
		*decomp = tmpDecomp;
}


/*--------------------------------------------------------------------------------------------------------------------------
** Public Functions
*/

void Decomposer_SetComm( void* decomposer, MPI_Comm comm ) {
	Decomposer*	self = (Decomposer*)decomposer;

	assert( self );

	self->comm = comm;
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Private Functions
*/

void Decomposer_BuildCommTopology( Decomposer* self, unsigned nDomains, unsigned* domains, 
				   CommTopology** commTopo, RangeSet*** isects )
{
	unsigned	rank, nProcs;
	RangeSet*	lSet;
	MPI_Group	worldGroup;
	unsigned*	worldRanks;
	unsigned*	subRanks;
	unsigned	nCommInc;
	unsigned*	commInc;
	RangeSet**	iSets;
	unsigned	nInc;
	unsigned*	inc;
	unsigned	p_i;

	assert( self );
	assert( !nDomains || domains );
	assert( commTopo );
	assert( isects );

	/* Get basic MPI info. */
	MPI_Comm_rank( self->comm, (int*)&rank );
	MPI_Comm_size( self->comm, (int*)&nProcs );

	/* We'll need to modify the world group. */
	MPI_Comm_group( self->comm, &worldGroup );
	worldRanks = Memory_Alloc_Array_Unnamed( unsigned, nProcs );
	subRanks = Memory_Alloc_Array_Unnamed( unsigned, nProcs );

	/* We need space to store index intersections. */
	iSets = Memory_Alloc_Array_Unnamed( RangeSet*, nProcs );
	memset( iSets, 0, nProcs * sizeof(RangeSet*) );

	/* Create a local set of required indices. */
	lSet = RangeSet_New();
	RangeSet_SetIndices( lSet, nDomains, domains );

	/* Tackle each processor one at a time. */
	for( p_i = 0; p_i < nProcs - 1; p_i++ ) {
		int		groupRange[3];
		MPI_Group	subGroup;
		MPI_Comm	subComm;
		unsigned	p_j;

		/* Set the processor range. */
		groupRange[0] = p_i;
		groupRange[1] = nProcs - 1;
		groupRange[2] = 1;

		/* We'll need a new group, as we only want to communicate using a triangular scheme. */
		MPI_Group_range_incl( worldGroup, 1, &groupRange, &subGroup );
		MPI_Comm_create( self->comm, subGroup, &subComm );

		/* Only continue if we're part of the sub-communicator. */
		if( rank >= p_i ) {
			unsigned	nBytes;
			Stg_Byte*	bytes;
			unsigned*	nFounds;
			Stg_Byte**	founds;

			/* Create a mapping between ranks. */
			for( p_j = p_i; p_j < nProcs; p_j++ )
				subRanks[p_j] = p_j - p_i;
			MPI_Group_translate_ranks( subGroup, nProcs - p_i, (int*)(subRanks + p_i), 
						   worldGroup, (int*)worldRanks );

			if( p_i == rank )
				RangeSet_Pickle( lSet, &nBytes, &bytes );

			MPIArray_Bcast( &nBytes, (void**)&bytes, sizeof(Stg_Byte), subRanks[p_i], subComm );

			if( p_i != rank ) {
				/* Create the intersection. */
				iSets[p_i] = RangeSet_New();
				RangeSet_Unpickle( iSets[p_i], nBytes, bytes );
				RangeSet_Intersection( iSets[p_i], lSet );

				/* Pickle the intersection to send back. */
				FreeArray( bytes );
				RangeSet_Pickle( iSets[p_i], &nBytes, &bytes );
			}
			else {
				KillArray( bytes );
				nBytes = 0;
			}

			/* Retrieve the results and unpickle each of them. */
			MPIArray_Gather( nBytes, bytes, &nFounds, (void***)&founds, sizeof(Stg_Byte), 
					 subRanks[p_i], subComm );
			if( p_i == rank ) {
				for( p_j = 0; p_j < nProcs - p_i; p_j++ ) {
					if( !nFounds[p_j] )
						continue;

					iSets[worldRanks[p_j]] = RangeSet_New();
					RangeSet_Unpickle( iSets[worldRanks[p_j]], nFounds[p_j], founds[p_j] );
				}

				/* Free the found arrays. */
				FreeArray( nFounds );
				FreeArray( founds );
			}
			else {
				/* Free pickled range set. */
				FreeArray( bytes );
			}

			/* Destroy the sub-communicator. */
			MPI_Comm_free( &subComm );
		}

		/* Destroy the sub-group. */
		MPI_Group_free( &subGroup );
	}

	/* Free rank translation arrays and local range set. */
	FreeArray( worldRanks );
	FreeArray( subRanks );

	/* Build a set of communication incidence. */
	nCommInc = 0;
	commInc = Memory_Alloc_Array_Unnamed( unsigned, nProcs );
	for( p_i = 0; p_i < nProcs; p_i++ ) {
		if( iSets[p_i] && iSets[p_i]->nInds )
			commInc[nCommInc++] = p_i;
	}

	/* Create the communication topology, unless one has already been specified. */
	*commTopo = CommTopology_New();
	CommTopology_SetComm( *commTopo, self->comm );
	CommTopology_SetIncidence( *commTopo, nCommInc, commInc );
	FreeArray( commInc );

	/* Build final intersections. */
	CommTopology_GetIncidence( *commTopo, &nInc, &inc );
	if( nInc ) {
		unsigned	inc_i;

		*isects = Memory_Alloc_Array_Unnamed( RangeSet*, nInc );
		for( inc_i = 0; inc_i < nInc; inc_i++ )
			(*isects)[inc_i] = iSets[inc[inc_i]];
	}
	else
		*isects = NULL;

	/* Free intersection array. */
	FreeArray( iSets );
}

void Decomposer_BuildLocalIntersections( Decomposer* self, unsigned nDomains, unsigned* domains, 
					 CommTopology* commTopo, RangeSet*** isects )
{
	unsigned	nIncRanks;
	RangeSet*	lSet;
	unsigned	nBytes, *nRemBytes, *nFndBytes;
	Stg_Byte	*bytes, **remBytes, **fndBytes;
	unsigned	p_i;

	assert( self );
	assert( !nDomains || domains );
	assert( commTopo );
	assert( isects );

	/* Get some crap. */
	nIncRanks = CommTopology_GetIncidenceSize( commTopo );

	/* We need space to store index intersections. */
	*isects = Memory_Alloc_Array_Unnamed( RangeSet*, nIncRanks );
	memset( *isects, 0, nIncRanks * sizeof(RangeSet*) );

	/* Create a local set of required indices, send to everyone. */
	lSet = RangeSet_New();
	RangeSet_SetIndices( lSet, nDomains, domains );
	RangeSet_Pickle( lSet, &nBytes, &bytes );
	CommTopology_Allgather( commTopo, 
				nBytes, bytes, 
				&nRemBytes, &remBytes, 
				sizeof(Stg_Byte) );

	/* Done with bytes. */
	FreeArray( bytes );

	/* Build intersections. */
	for( p_i = 0; p_i < nIncRanks; p_i++ ) {
		(*isects)[p_i] = RangeSet_New();
		RangeSet_Unpickle( (*isects)[p_i], nRemBytes[p_i], remBytes[p_i] );
		RangeSet_Intersection( (*isects)[p_i], lSet );
		FreeArray( remBytes[p_i] );
		RangeSet_Pickle( (*isects)[p_i], nRemBytes + p_i, remBytes + p_i );
	}

	/* Send back. */
	CommTopology_Alltoall( commTopo, 
			       nRemBytes, remBytes, 
			       &nFndBytes, &fndBytes, 
			       sizeof(Stg_Byte) );

	/* Free unused arrays. */
	FreeArray( nRemBytes );
	FreeArray2D( nIncRanks, remBytes );

	/* Extract our intersections. */
	for( p_i = 0; p_i < nIncRanks; p_i++ )
		RangeSet_Unpickle( (*isects)[p_i], nFndBytes[p_i], fndBytes[p_i] );

	/* Free remaining arrays. */
	FreeArray( nFndBytes );
	FreeArray2D( nIncRanks, fndBytes );
}

void Decomposer_Claim( Decomposer* self, CommTopology* commTopo, RangeSet** isects, 
		       unsigned nDomains, unsigned* domains, 
		       Decomp** decomp, Decomp_Sync** sync )
{
	unsigned	rank, nProcs;
	MPI_Comm	comm;
	unsigned	nInc;
	unsigned*	inc;
	RangeSet*	lSet;
	unsigned	nBytes;
	Stg_Byte*	bytes;
	RangeSet*	tmpClaimed;
	unsigned	tag = 6669;
	unsigned	p_i, p_j;

	/* Get basic MPI info. */
	comm = CommTopology_GetComm( commTopo );
	MPI_Comm_rank( comm, (int*)&rank );
	MPI_Comm_size( comm, (int*)&nProcs );

	/* Build a set of domains. */
	lSet = RangeSet_New();
	RangeSet_SetIndices( lSet, nDomains, domains );

	/* Extract our neighbouring processors. */
	CommTopology_GetIncidence( commTopo, &nInc, &inc );

	/* Figure out where info is coming from and going to. Note that the incidence
	   is always ordered from lowest to highest rank. */
	tmpClaimed = RangeSet_New();
	for( p_i = 0; p_i < nInc; p_i++ ) {
		MPI_Status	status;

		if( inc[p_i] > rank )
			break;

		/* Receive from neighbour which indices it has taken. */
		MPI_Recv( &nBytes, 1, MPI_UNSIGNED, inc[p_i], tag, comm, &status );
		bytes = Memory_Alloc_Array_Unnamed( Stg_Byte, nBytes );
		MPI_Recv( bytes, nBytes, MPI_BYTE, inc[p_i], tag, comm, &status );
		RangeSet_Unpickle( tmpClaimed, nBytes, bytes );
		FreeArray( bytes );

		/* Subtract from our claimed set. */
		RangeSet_Subtraction( lSet, tmpClaimed );
	}
	FreeObject( tmpClaimed );

	/* Update remaining neighbours as to which indices we've taken. */
	for( p_j = p_i; p_j < nInc; p_j++ ) {
		RangeSet*	intersect;

		intersect = RangeSet_DeepCopy( isects[p_j] );
		RangeSet_Intersection( intersect, lSet );
		RangeSet_Pickle( intersect, &nBytes, &bytes );
		MPI_Send( &nBytes, 1, MPI_UNSIGNED, inc[p_j], tag, comm );
		MPI_Send( bytes, nBytes, MPI_BYTE, inc[p_j], tag, comm );
		FreeArray( bytes );
		FreeObject( intersect );
	}

	/* Extract our ownership. */
	Decomposer_BuildIndices( self, nDomains, domains, lSet, 
				 commTopo, decomp, sync );

	/* Free local range set. */
	FreeObject( lSet );
}

void Decomposer_BuildIndices( Decomposer* self, unsigned nDomains, unsigned* domains, RangeSet* claimed, 
			      CommTopology* commTopo, Decomp** decomp, Decomp_Sync** sync )
{
	RangeSet*	rSet;
	unsigned	nLocals, *locals;
	unsigned	nRemotes, *remotes;

	/* Extract local indices. */
	*decomp = Decomp_New();
	Decomp_SetComm( *decomp, self->comm );
	locals = NULL;
	RangeSet_GetIndices( claimed, &nLocals, &locals );
	Decomp_SetLocals( *decomp, nLocals, locals );
	FreeArray( locals );

	/* Build a set of remotes. */
	rSet = RangeSet_New();
	RangeSet_SetIndices( rSet, nDomains, domains );
	RangeSet_Subtraction( rSet, claimed );

	/* Set remote indices. */
	*sync = Decomp_Sync_New();
	Decomp_Sync_SetDecomp( *sync, *decomp );
	Decomp_Sync_SetCommTopology( *sync, commTopo );
	remotes = NULL;
	RangeSet_GetIndices( rSet, &nRemotes, &remotes );
	Decomp_Sync_SetRemotes( *sync, nRemotes, remotes );
	FreeArray( remotes );

	/* Destroy rental set. */
	FreeObject( rSet );
}
