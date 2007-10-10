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
** $Id: CommTopology.c 3584 2006-05-16 11:11:07Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <mpi.h>

#include <StGermain/StGermain.h>
#include <StgDomain/Geometry/Geometry.h>
#include "Mesh.h"


/* Textual name of this class */
const Type CommTopology_Type = "CommTopology";


/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/

CommTopology* CommTopology_New() {
	return _CommTopology_New( sizeof(CommTopology), 
				  CommTopology_Type, 
				  _CommTopology_Delete, 
				  _CommTopology_Print, 
				  NULL );
}

CommTopology* _CommTopology_New( COMMTOPOLOGY_DEFARGS ) {
	CommTopology* self;
	
	/* Allocate memory */
	assert( sizeOfSelf >= sizeof(CommTopology) );
	self = (CommTopology*)_Stg_Class_New( STG_CLASS_PASSARGS );

	/* Virtual info */

	/* CommTopology info */
	_CommTopology_Init( self );

	return self;
}

void _CommTopology_Init( CommTopology* self ) {
	self->comm = MPI_COMM_WORLD;
	self->nIncRanks = 0;
	self->incRanks = NULL;
	self->glMap = UIntMap_New();
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _CommTopology_Delete( void* commTopology ) {
	CommTopology*	self = (CommTopology*)commTopology;

	CommTopology_Destruct( self );
	FreeObject( self->glMap );

	/* Delete the parent. */
	_Stg_Class_Delete( self );
}

void _CommTopology_Print( void* commTopology, Stream* stream ) {
	CommTopology*	self = (CommTopology*)commTopology;
	
	/* Set the Journal for printing informations */
	Stream* commTopologyStream;
	commTopologyStream = Journal_Register( InfoStream_Type, "CommTopologyStream" );

	/* Print parent */
	Journal_Printf( stream, "CommTopology (ptr): (%p)\n", self );
	_Stg_Class_Print( self, stream );
}


/*--------------------------------------------------------------------------------------------------------------------------
** Public Functions
*/

void CommTopology_SetComm( void* commTopology, MPI_Comm comm ) {
	CommTopology*	self = (CommTopology*)commTopology;

	assert( self );

	CommTopology_Destruct( self );
	self->comm = comm;
}

void CommTopology_SetIncidence( void* commTopology, unsigned nIncRanks, unsigned* incRanks ) {
	CommTopology*	self = (CommTopology*)commTopology;

	assert( self );
	assert( !nIncRanks || incRanks );

	CommTopology_Destruct( self );
	CommTopology_AddIncidence( self, nIncRanks, incRanks );
}

void CommTopology_AddIncidence( void* commTopology, unsigned nIncRanks, unsigned* incRanks ) {
	CommTopology*	self = (CommTopology*)commTopology;
	unsigned	r_i;

	assert( self );
	assert( !nIncRanks || incRanks );
	assert( CommTopology_ValidateIncidence( self, nIncRanks, incRanks ) );

	if( nIncRanks ) {
		unsigned*	order;

		for( r_i = 0; r_i < nIncRanks; r_i++ )
			UIntMap_Insert( self->glMap, incRanks[r_i], self->nIncRanks + r_i );

		self->incRanks = ReallocNamedArray( self->incRanks, unsigned, self->nIncRanks + nIncRanks, 
						    "CommTopology::incRanks" );
		memcpy( self->incRanks + self->nIncRanks, incRanks, nIncRanks * sizeof(unsigned) );
		self->nIncRanks += nIncRanks;

		order = AllocArray( unsigned, self->nIncRanks * 2 );
		for( r_i = 0; r_i < self->nIncRanks; r_i++ ) {
			order[2 * r_i] = self->incRanks[r_i];
			order[2 * r_i + 1] = r_i;
		}
		qsort( order, self->nIncRanks, 2 * sizeof(unsigned), CommTopology_CmpRanks );
		self->order = ReallocNamedArray( self->order, unsigned, self->nIncRanks, "CommTopology::order" );
		for( r_i = 0; r_i < self->nIncRanks; r_i++ )
			self->order[r_i] = order[2 * r_i + 1];
		FreeArray( order );
	}
}

MPI_Comm CommTopology_GetComm( void* commTopology ) {
	CommTopology*	self = (CommTopology*)commTopology;

	assert( self );

	return self->comm;
}

unsigned CommTopology_GetCommSize( void* commTopology ) {
	CommTopology*	self = (CommTopology*)commTopology;
	unsigned	nProcs;

	assert( self );

	MPI_Comm_size( self->comm, (int*)&nProcs );

	return nProcs;
}

unsigned CommTopology_GetCommRank( void* commTopology ) {
	CommTopology*	self = (CommTopology*)commTopology;
	unsigned	rank;

	assert( self );

	MPI_Comm_rank( self->comm, (int*)&rank );

	return rank;
}

unsigned CommTopology_GetIncidenceSize( void* commTopology ) {
	CommTopology*	self = (CommTopology*)commTopology;

	assert( self );

	return self->nIncRanks;
}

void CommTopology_GetIncidence( void* commTopology, unsigned* nIncRanks, unsigned** incRanks ) {
	CommTopology*	self = (CommTopology*)commTopology;

	assert( self );
	assert( nIncRanks && incRanks );

	*nIncRanks = self->nIncRanks;
	*incRanks = self->incRanks;
}

unsigned CommTopology_LocalToGlobal( void* commTopology, unsigned local ) {
	CommTopology*	self = (CommTopology*)commTopology;

	assert( self );
	assert( local < self->nIncRanks );
	assert( self->incRanks );
	assert( self->order );

	return self->incRanks[local];
}

Bool CommTopology_GlobalToLocal( void* commTopology, unsigned global, unsigned* local ) {
	CommTopology*	self = (CommTopology*)commTopology;

	assert( self );

	return UIntMap_Map( self->glMap, global, local );
}

void _CommTopology_Allgather( void* commTopology, 
			      unsigned srcSize, void* srcArray, 
			      unsigned** dstSizes, void*** dstArrays, 
			      unsigned itemSize )
{
	CommTopology*	self = (CommTopology*)commTopology;
	unsigned	tag = 669;
	MPI_Status	status;
	unsigned*	nbrSizes;
	Stg_Byte**	nbrArrays;
	unsigned	p_i;

	assert( self );
	assert( !srcSize || srcArray );
	assert( dstSizes );
	assert( dstArrays );

	/* Skip this if we have no neighbours. */
	if( !self->nIncRanks ) {
		*dstSizes = NULL;
		*dstArrays = NULL;
		return;
	}

	/* Allocate base arrays. */
	nbrSizes = AllocArray( unsigned, self->nIncRanks );
	nbrArrays = AllocArray( Stg_Byte*, self->nIncRanks );

	/* Send/recv with each neighbour.  There won't be deadlocks because the neighbouring
	   ranks are ordered from lowest to highest. Start with sizes. */
	for( p_i = 0; p_i < self->nIncRanks; p_i++ ) {
		unsigned	pos = self->order[p_i];
		unsigned	nbr = self->incRanks[pos];
		unsigned	tmpSize;

		MPI_Sendrecv( &srcSize, 1, MPI_UNSIGNED, nbr, tag, 
			      nbrSizes + pos, 1, MPI_UNSIGNED, nbr, tag, 
			      self->comm, &status );

		tmpSize = nbrSizes[pos] * itemSize;

		/* Alloc storage. */
		if( tmpSize )
			nbrArrays[pos] = AllocArray( Stg_Byte, tmpSize );
		else
			nbrArrays[pos] = NULL;

		MPI_Sendrecv( srcArray, srcSize * itemSize, MPI_BYTE, nbr, tag, 
			      nbrArrays[pos], tmpSize, MPI_BYTE, nbr, tag, 
			      self->comm, &status );
	}

	/* Store results. */
	*dstSizes = nbrSizes;
	*dstArrays = (void**)nbrArrays;
}

void _CommTopology_Alltoall( void* commTopology, 
			     unsigned* srcSizes, void** srcArrays, 
			     unsigned** dstSizes, void*** dstArrays, 
			     unsigned itemSize )
{
	CommTopology*	self = (CommTopology*)commTopology;
	unsigned	tag = 669;
	MPI_Status	status;
	unsigned*	nbrSizes;
	Stg_Byte**	nbrArrays;
	unsigned	p_i;

	assert( self );
	assert( !srcSizes || srcArrays );
	assert( dstSizes );
	assert( dstArrays );

	/* Skip this if we have no neighbours. */
	if( !self->nIncRanks ) {
		*dstSizes = NULL;
		*dstArrays = NULL;
		return;
	}

	/* Allocate base array. */
	nbrSizes = AllocArray( unsigned, self->nIncRanks );
	nbrArrays = AllocArray( Stg_Byte*, self->nIncRanks );

	/* Send/recv with each neighbour.  There won't be deadlocks because the neighbouring
	   ranks are ordered from lowest to highest. Start with sizes. */
	for( p_i = 0; p_i < self->nIncRanks; p_i++ ) {
		unsigned	pos = self->order[p_i];
		unsigned	nbr = self->incRanks[pos];
		unsigned	tmpSize;

		MPI_Sendrecv( srcSizes + pos, 1, MPI_UNSIGNED, nbr, tag, 
			      nbrSizes + pos, 1, MPI_UNSIGNED, nbr, tag, 
			      self->comm, &status );

		tmpSize = nbrSizes[pos] * itemSize;

		/* Alloc storage. */
		if( tmpSize )
			nbrArrays[pos] = AllocArray( Stg_Byte, tmpSize );
		else
			nbrArrays[pos] = NULL;

		/* Transfer between current neighbour. */
		MPI_Sendrecv( srcArrays[pos], srcSizes[pos] * itemSize, MPI_BYTE, nbr, tag, 
			      nbrArrays[pos], tmpSize, MPI_BYTE, nbr, tag, 
			      self->comm, &status );
	}

	/* Store results. */
	*dstSizes = nbrSizes;
	*dstArrays = (void**)nbrArrays;
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Private Functions
*/

int CommTopology_CmpRanks( const void* rank0, const void* rank1 ) {
	if( *((unsigned*)rank0) < *((unsigned*)rank1) )
		return -1;
	else if( *((unsigned*)rank0) > *((unsigned*)rank1) )
		return 1;
	else
		return 0;
}

void CommTopology_Destruct( CommTopology* self ) {
	KillArray( self->incRanks );
	self->nIncRanks = 0;
	UIntMap_Clear( self->glMap );
}

#ifndef NDEBUG
Bool CommTopology_ValidateIncidence( CommTopology* self, unsigned nIncRanks, unsigned* incRanks ) {
	RangeSet	*exSet, *newSet;
	unsigned	nProcs, rank;
	unsigned	size;
	unsigned	inc_i;

	assert( self );
	assert( !nIncRanks || incRanks );

	/* Validate basics. */
	MPI_Comm_size( self->comm, (int*)&nProcs );
	MPI_Comm_rank( self->comm, (int*)&rank );
	for( inc_i = 0; inc_i < nIncRanks; inc_i++ ) {
		if( incRanks[inc_i] >= nProcs || incRanks[inc_i] == rank )
			return False;
	}

	/* Ensure no existing overlap. */
	exSet = RangeSet_New();
	newSet = RangeSet_New();
	RangeSet_SetIndices( exSet, self->nIncRanks, self->incRanks );
	RangeSet_SetIndices( newSet, nIncRanks, incRanks );
	RangeSet_Intersection( newSet, exSet );
	size = RangeSet_GetSize( newSet );
	FreeObject( newSet );
	FreeObject( exSet );
	if( size )
		return False;

	return True;
}
#endif
