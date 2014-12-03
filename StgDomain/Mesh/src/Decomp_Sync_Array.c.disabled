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

#include "types.h"
#include "shortcuts.h"
#include "CommTopology.h"
#include "Decomp.h"
#include "Decomp_Sync.h"
#include "Decomp_Sync_Array.h"


/* Textual name of this class */
const Type Decomp_Sync_Array_Type = "Decomp_Sync_Array";

/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/

Decomp_Sync_Array* Decomp_Sync_Array_New() {
	return _Decomp_Sync_Array_New( sizeof(Decomp_Sync_Array), 
				       Decomp_Sync_Array_Type, 
				       _Decomp_Sync_Array_Delete, 
				       _Decomp_Sync_Array_Print, 
				       NULL );
}

Decomp_Sync_Array* _Decomp_Sync_Array_New( DECOMP_SYNC_ARRAY_DEFARGS ) {
	Decomp_Sync_Array* self;
	
	/* Allocate memory */
	assert( sizeOfSelf >= sizeof(Decomp_Sync_Array) );
	self = (Decomp_Sync_Array*)_Stg_Class_New( STG_CLASS_PASSARGS );

	/* Virtual info */

	/* Decomp_Sync_Array info */
	_Decomp_Sync_Array_Init( self );

	return self;
}

void _Decomp_Sync_Array_Init( Decomp_Sync_Array* self ) {
	self->sync = NULL;

	self->snkArray = NULL;
	self->snkStride = 0;
	self->snkDisps = NULL;
	self->snkSizes = NULL;
	self->snkOffs = NULL;

	self->srcArray = NULL;
	self->srcStride = 0;
	self->srcDisps = NULL;
	self->srcSizes = NULL;
	self->srcOffs = NULL;

	self->itemSize = 0;
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _Decomp_Sync_Array_Delete( void* sync ) {
	Decomp_Sync_Array*	self = (Decomp_Sync_Array*)sync;

	Decomp_Sync_Array_Destruct( self );

	/* Delete the parent. */
	_Stg_Class_Delete( self );
}

void _Decomp_Sync_Array_Print( void* sync, Stream* stream ) {
	Decomp_Sync_Array*	self = (Decomp_Sync_Array*)sync;

	/* Set the Journal for printing informations */
	Stream* syncStream;
	syncStream = Journal_Register( InfoStream_Type, "Decomp_Sync_ArrayStream" );

	/* Print parent */
	Journal_Printf( stream, "Decomp_Sync_Array (ptr): (%p)\n", self );
	_Stg_Class_Print( self, stream );
}


/*--------------------------------------------------------------------------------------------------------------------------
** Public Functions
*/

void Decomp_Sync_Array_SetSync( void* array, Decomp_Sync* sync ) {
	Decomp_Sync_Array*	self = (Decomp_Sync_Array*)array;

	assert( self );

	Decomp_Sync_Array_Destruct( self );

	self->sync = sync;
	if( sync )
		List_Append( sync->arrays, &self );
}

void Decomp_Sync_Array_SetMemory( void* array, 
				  void* localArray, void* remoteArray, 
				  size_t localStride, size_t remoteStride, 
				  size_t itemSize )
{
	Decomp_Sync_Array*	self = (Decomp_Sync_Array*)array;

	/* Sanity checks. */
	assert( self );
	assert( self->sync );

	/* Store information. */
	self->snkArray = localArray;
	self->snkStride = localStride;
	self->srcArray = remoteArray;
	self->srcStride = remoteStride;
	self->itemSize = itemSize;

	/* Build this array. */
	Decomp_Sync_Array_BuildArray( self );
}

void Decomp_Sync_Array_Sync( void* array ) {
	Decomp_Sync_Array*	self = (Decomp_Sync_Array*)array;

	assert( self );

	Decomp_Sync_SyncArray( self->sync, self );
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Private Functions
*/

void Decomp_Sync_Array_BuildArray( Decomp_Sync_Array* self ) {
	Decomp*		decomp;
	CommTopology*	commTopo;
	unsigned	nInc, *inc;

	assert( self );
	assert( self->sync );

	/* Shortcuts. */
	decomp = Decomp_Sync_GetDecomp( self->sync );
	commTopo = Decomp_Sync_GetCommTopology( self->sync );

	/* Extract incidence. */
	CommTopology_GetIncidence( commTopo, &nInc, &inc );

	if( nInc ) {
		/* Determine sink (local) information. */
		if( self->sync->netSnks > 0 ) {
			unsigned*	snkOffs;
			unsigned*	snkSizes;
			unsigned*	snkDisps;
			unsigned	snkInd = 0;
			unsigned	p_i;

			/* Allocate/reallocate memory. */
			snkDisps = Memory_Alloc_Array( unsigned, nInc, "Decomp_Sync_Array::snkDisps" );
			snkSizes = Memory_Alloc_Array( unsigned, nInc, "Decomp_Sync_Array::snkSizes" );
			snkOffs = Memory_Alloc_Array( unsigned, self->sync->netSnks, "Decomp_Sync_Array::snkOffs" );

			/* Calculate offsets and sizes. */
			for( p_i = 0; p_i < nInc; p_i++ ) {
				unsigned	snk_i;

				snkSizes[p_i] = 0;
				for( snk_i = 0; snk_i < self->sync->nSnks[p_i]; snk_i++ ) {
					unsigned	dInd;

					dInd = self->sync->snks[p_i][snk_i];
					snkOffs[snkInd] = dInd * self->snkStride;
					snkSizes[p_i] += self->itemSize;
					snkInd++;
				}
			}

			/* Calculate the displacements. */
			snkDisps[0] = 0;
			for( p_i = 1; p_i < nInc; p_i++ )
				snkDisps[p_i] = snkDisps[p_i - 1] + snkSizes[p_i - 1];

			/* Store arrays. */
			self->snkOffs = snkOffs;
			self->snkDisps = snkDisps;
			self->snkSizes = snkSizes;
		}
		else {
			/* Store null information. */
			self->snkOffs = NULL;
			self->snkDisps = Memory_Alloc_Array( unsigned, nInc, "Decomp_Sync_Array::snkDisps" );
			self->snkSizes = Memory_Alloc_Array( unsigned, nInc, "Decomp_Sync_Array::snkSizes" );
			memset( self->snkDisps, 0, nInc * sizeof(unsigned) );
			memset( self->snkSizes, 0, nInc * sizeof(unsigned) );
		}

		/* Determine source (shadow) information. */
		if( self->sync->netSrcs > 0 ) {
			unsigned*	srcOffs;
			unsigned*	srcSizes;
			unsigned*	srcDisps;
			unsigned	srcInd = 0;
			unsigned	p_i;

			/* Allocate/reallocate memory. */
			srcDisps = Memory_Alloc_Array( unsigned, nInc, "Decomp_Sync_Array::srcDisps" );
			srcSizes = Memory_Alloc_Array( unsigned, nInc, "Decomp_Sync_Array::srcSizes" );
			srcOffs = Memory_Alloc_Array( unsigned, self->sync->netSrcs, "Decomp_Sync_Array::srcOffs" );

			/* Calculate offsets and sizes. */
			for( p_i = 0; p_i < nInc; p_i++ ) {
				unsigned	src_i;

				srcSizes[p_i] = 0;
				for( src_i = 0; src_i < self->sync->nSrcs[p_i]; src_i++ ) {
					unsigned	sInd;

					sInd = self->sync->srcs[p_i][src_i];
					assert( sInd >= Decomp_GetLocalSize( decomp ) );
					sInd -= Decomp_GetLocalSize( decomp );
					srcOffs[srcInd] = sInd * self->srcStride;
					srcSizes[p_i] += self->itemSize;
					srcInd++;
				}
			}

			/* Calculate the displacements. */
			srcDisps[0] = 0;
			for( p_i = 1; p_i < nInc; p_i++ )
				srcDisps[p_i] = srcDisps[p_i - 1] + srcSizes[p_i - 1];

			/* Store arrays. */
			self->srcOffs = srcOffs;
			self->srcDisps = srcDisps;
			self->srcSizes = srcSizes;
		}
		else {
			/* Store null information. */
			self->srcOffs = NULL;
			self->srcDisps = Memory_Alloc_Array( unsigned, nInc, "Decomp_Sync_Array::srcDisps" );
			self->srcSizes = Memory_Alloc_Array( unsigned, nInc, "Decomp_Sync_Array::srcSizes" );
			memset( self->srcDisps, 0, nInc * sizeof(unsigned) );
			memset( self->srcSizes, 0, nInc * sizeof(unsigned) );
		}
	}
	else {
		self->snkOffs = NULL;
		self->snkDisps = NULL;
		self->snkSizes = NULL;
		self->srcOffs = NULL;
		self->srcDisps = NULL;
		self->srcSizes = NULL;
	}
}

void Decomp_Sync_Array_Destruct( Decomp_Sync_Array* self ) {
	if( self->sync )
		List_Remove( self->sync->arrays, &self );

	self->snkArray = NULL;
	self->snkStride = 0;
	self->srcArray = NULL;
	self->srcStride = 0;
	self->itemSize = 0;

	KillArray( self->snkDisps );
	KillArray( self->snkSizes );
	KillArray( self->snkOffs );
	KillArray( self->srcDisps );
	KillArray( self->srcSizes );
	KillArray( self->srcOffs );
	self->itemSize = 0;
}
