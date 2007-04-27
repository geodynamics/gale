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
** $Id: Decomp.c 3584 2006-05-16 11:11:07Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <mpi.h>

#include "Base/Base.h"
#include "Mesh.h"


/* Textual name of this class */
const Type Decomp_Type = "Decomp";


/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/

Decomp* Decomp_New( Name name ) {
	return _Decomp_New( sizeof(Decomp), 
			    Decomp_Type, 
			    _Decomp_Delete, 
			    _Decomp_Print, 
			    NULL );
}

Decomp* _Decomp_New( DECOMP_DEFARGS ) {
	Decomp* self;
	
	/* Allocate memory */
	assert( sizeOfSelf >= sizeof(Decomp) );
	self = (Decomp*)_Stg_Class_New( STG_CLASS_PASSARGS );

	/* Virtual info */

	/* Decomp info */
	_Decomp_Init( self );

	return self;
}

void _Decomp_Init( Decomp* self ) {
	self->comm = MPI_COMM_WORLD;

	self->nGlobals = 0;
	self->nLocals = 0;
	self->locals = NULL;

	self->glMap = UIntMap_New();

	self->syncs = List_New();
	List_SetItemSize( self->syncs, sizeof(Decomp_Sync*) );
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _Decomp_Delete( void* decomp ) {
	Decomp*	self = (Decomp*)decomp;

	Decomp_Destruct( self );
	FreeObject( self->glMap );
	FreeObject( self->syncs );

	/* Delete the parent. */
	_Stg_Class_Delete( self );
}

void _Decomp_Print( void* decomp, Stream* stream ) {
	Decomp*	self = (Decomp*)decomp;
	
	/* Set the Journal for printing informations */
	Stream* decompStream;
	decompStream = Journal_Register( InfoStream_Type, "DecompStream" );

	/* Print parent */
	Journal_Printf( stream, "Decomp (ptr): (%p)\n", self );
	_Stg_Class_Print( self, stream );
}


/*--------------------------------------------------------------------------------------------------------------------------
** Public Functions
*/

void Decomp_SetComm( void* decomp, MPI_Comm comm ) {
	Decomp*	self = (Decomp*)decomp;

	assert( self );

	Decomp_Destruct( self );
	self->comm = comm;
}

void Decomp_SetLocals( void* decomp, unsigned nLocals, unsigned* locals ) {
	Decomp*	self = (Decomp*)decomp;

	assert( self );
	assert( !nLocals || locals );
	assert( Decomp_ValidateDomain( self, nLocals, locals ) );

	Decomp_Destruct( self );
	self->nLocals = nLocals;
	if( self->nLocals ) {
		self->locals = AllocNamedArray( unsigned, nLocals, "Decomp::locals" );
		memcpy( self->locals, locals, nLocals * sizeof(unsigned) );
	}

	Decomp_CalcGlobalSize( self );
	Decomp_BuildGLMap( self );
	Decomp_UpdateSyncs( self );
}

MPI_Comm Decomp_GetComm( void* decomp ) {
	Decomp*	self = (Decomp*)decomp;

	assert( self );

	return self->comm;
}

unsigned Decomp_GetGlobalSize( void* decomp ) {
	Decomp*	self = (Decomp*)decomp;

	assert( self );

	return self->nGlobals;
}

unsigned Decomp_GetLocalSize( void* decomp ) {
	Decomp*	self = (Decomp*)decomp;

	assert( self );

	return self->nLocals;
}

void Decomp_GetLocals( void* decomp, unsigned* nLocals, unsigned** locals ) {
	Decomp*	self = (Decomp*)decomp;

	assert( self );
	assert( nLocals );
	assert( locals );

	*nLocals = self->nLocals;
	*locals = self->locals;
}

List* Decomp_GetSyncList( void* decomp ) {
	Decomp*	self = (Decomp*)decomp;

	assert( self );

	return self->syncs;
}

Bool Decomp_GlobalToLocal( void* decomp, unsigned global, unsigned* local ) {
	Decomp*		self = (Decomp*)decomp;

	assert( self );
	assert( global < self->nGlobals );
	assert( local );

	return UIntMap_Map( self->glMap, global, local );
}

unsigned Decomp_LocalToGlobal( void* decomp, unsigned local ) {
	Decomp*	self = (Decomp*)decomp;

	assert( self );
	assert( local < self->nLocals );
	assert( self->locals );

	return self->locals[local];
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Private Functions
*/

void Decomp_BuildGLMap( Decomp* self ) {
	unsigned	l_i;

	assert( self );

	UIntMap_Clear( self->glMap );
	for( l_i = 0; l_i < self->nLocals; l_i++ )
		UIntMap_Insert( self->glMap, self->locals[l_i], l_i );
}

void Decomp_CalcGlobalSize( Decomp* self ) {
	assert( self );

	MPI_Allreduce( &self->nLocals, &self->nGlobals, 1, MPI_UNSIGNED, MPI_SUM, self->comm );
}

void Decomp_UpdateSyncs( Decomp* self ) {
	unsigned	s_i;

	assert( self );

	for( s_i = 0; s_i < List_GetSize( self->syncs ); s_i++ ) {
		Decomp_Sync*	sync;

		sync = *List_Get( self->syncs, s_i, Decomp_Sync* );
		Decomp_Sync_Update( sync );
	}
}

void Decomp_Destruct( Decomp* self ) {
	assert( self );

	Decomp_DestructSyncs( self );

	self->nLocals = 0;
	KillArray( self->locals );
	UIntMap_Clear( self->glMap );
}

void Decomp_DestructSyncs( Decomp* self ) {
	unsigned	s_i;

	assert( self );

	for( s_i = 0; s_i < List_GetSize( self->syncs ); s_i++ ) {
		Decomp_Sync*	sync;

		sync = *List_Get( self->syncs, s_i, Decomp_Sync* );
		Decomp_Sync_SetDecomp( sync, NULL );
	}

	List_Clear( self->syncs );
}

#ifndef NDEBUG
Bool Decomp_ValidateDomain( Decomp* self, unsigned nLocals, unsigned* locals ) {
	unsigned	rank, nProcs;
	unsigned	nBytes;
	Stg_Byte*	bytes;
	unsigned	tag = 5559;
	RangeSet*	lSet;
	RangeSet*	gSet;
	unsigned	netInds;
	unsigned	nGlobals;

	assert( self );
	assert( !nLocals || locals );

	/* Get basic MPI info. */
	MPI_Comm_rank( self->comm, (int*)&rank );
	MPI_Comm_size( self->comm, (int*)&nProcs );

	/* Create a local range set. */
	lSet = RangeSet_New();
	RangeSet_SetIndices( lSet, nLocals, locals );

	/* Create a global range set. */
	gSet = RangeSet_New();

	/* Receive from previous rank. */
	if( rank > 0 ) {
		MPI_Status	mpiStatus;

		MPI_Recv( &nBytes, 1, MPI_UNSIGNED, rank - 1, tag, self->comm, &mpiStatus );
		bytes = AllocArray( Stg_Byte, nBytes );
		MPI_Recv( bytes, nBytes, MPI_BYTE, rank - 1, tag, self->comm, &mpiStatus );
		RangeSet_Unpickle( gSet, nBytes, bytes );
		FreeArray( bytes );
	}

	/* Combine sets. */
	RangeSet_Union( gSet, lSet );

	if( rank < nProcs - 1 ) {
		RangeSet_Pickle( gSet, &nBytes, &bytes );
		MPI_Send( &nBytes, 1, MPI_UNSIGNED, rank + 1, tag, self->comm );
		MPI_Send( bytes, nBytes, MPI_BYTE, rank + 1, tag, self->comm );
		FreeArray( bytes );
	}
	else {
		nGlobals = RangeSet_GetSize( gSet );
		if( RangeSet_GetNumRanges( gSet ) > 1 )
			return False;
	}

	/* Transfer global count to all. */
	MPI_Bcast( &nGlobals, 1, MPI_UNSIGNED, nProcs - 1, self->comm );

	/* Check for overlap. */
	MPI_Allreduce( &lSet->nInds, &netInds, 1, MPI_UNSIGNED, MPI_SUM, self->comm );
	if( netInds != nGlobals )
		return False;

	/* Free the sets. */
	FreeObject( lSet );
	FreeObject( gSet );

	return True;
}
#endif
