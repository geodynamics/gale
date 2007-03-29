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
** $Id: List.c 2038 2004-11-15 07:11:01Z RaquibulHassan $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <Base/Foundation/Foundation.h>
#include <Base/IO/IO.h>

#include "types.h"
#include "MemoryPool.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

const Type MemoryPool_Type = "MemoryPool";

/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/
MemoryPool* _MemoryPool_New(
			SizeT							_sizeOfSelf,
			Type							type,
			Stg_Class_DeleteFunction*			_delete,
			Stg_Class_PrintFunction*			_print,
			Stg_Class_CopyFunction*				_copy,
			int									elementSize,
			int									numElements,
			int									delta
			)
{
	MemoryPool *self = NULL;

	self = (MemoryPool*)_Stg_Class_New(
							_sizeOfSelf,
							type,
							_delete,
							_print,
							_copy);

	self->numChunks = 1;
	self->delta = delta;
	self->elementSize = elementSize;
	self->numElements = numElements;
	self->numInitialElements = numElements;
	self->numElementsFree = numElements;

	return self;
}

MemoryPool* MemoryPool_NewFunc( SizeT elementSize, int numElements, int delta )
{
	
	MemoryPool* self;
	
	self = _MemoryPool_New( 
		   					sizeof(MemoryPool),
							MemoryPool_Type,
							_MemoryPool_DeleteFunc,
							_MemoryPool_PrintFunc,
							NULL,
							elementSize,
							numElements,
							delta
							);

	/* Virtual functions */
	_MemoryPool_Init( self );
	return self;
}

void _MemoryPool_Init( MemoryPool* self ){
	int i = 0;
	
	/* General info */
	
	assert( self );
	_Stg_Class_Init ((Stg_Class*) self);
	
	assert( self->numChunks == 1 );

	self->elements = (char**)Memory_Alloc_Bytes_Unnamed( sizeof(void*)*self->numChunks, "char*" );
	self->elements[self->numChunks-1] = Memory_Alloc_Bytes_Unnamed( self->elementSize * self->numElements, "char" );
	memset( self->elements[self->numChunks-1], 0, self->elementSize * self->numElements );

	self->pool = Memory_Alloc_Bytes_Unnamed( sizeof( char* ) * self->numElements, "char*" );
	memset( self->pool, 0, sizeof(char*) * self->numElements );
	
	for( i=0; i<self->numElements; i++ ){
		self->pool[i] = &(self->elements[self->numChunks-1][i*self->elementSize]);
	}
}

void *MemoryPool_Init( MemoryPool *self )
{
	/*TODO*/
	return NULL;
}

/*----------------------------------------------------------------------------------------------------------------------------------
** Virtual Functions
*/
void _MemoryPool_PrintFunc ( void *memPool, Stream *stream )
{
	MemoryPool *self = NULL;
	
	self = (MemoryPool*) memPool;
	assert( self );
	assert( stream );

	/* print parent */
	_Stg_Class_Print( (void*) self, stream );

	/* general info */
	Journal_Printf( stream, "MemoryPool (ptr): (%p)\n", self );

	/* Virtual Info */

	/* MemoryPool Info */
	Journal_Printf( stream, "\tElementSize\t\t - %d\n", self->elementSize );
	Journal_Printf( stream, "\tNumElememts\t\t - %d\n", self->numElements );
	Journal_Printf( stream, "\tNumElememtsFree\t\t - %d\n", self->numElementsFree );
	Journal_Printf( stream, "\tNumChunks\t\t - %d\n", self->numChunks );
	Journal_Printf( stream, "\tDelta\t\t - %d\n", self->delta );
}

void _MemoryPool_DeleteFunc( void *memPool )
{
	int i = 0;
	MemoryPool *self = NULL;
	
	self = (MemoryPool*)memPool;
	assert (self);
	
	for( i=0; i<self->numChunks; i++ ){
		Memory_Free( self->elements[i] );
	}

	Memory_Free( self->elements );
	Memory_Free( self->pool );
	
	_Stg_Class_Delete( self );
}

/*----------------------------------------------------------------------------------------------------------------------------------
** Public Functions
*/

void *MemoryPool_NewObjectFunc( SizeT elementSize, MemoryPool *memPool )
{
	int index = 0;
	
	assert( elementSize == memPool->elementSize );
	
label:	index = memPool->numElementsFree - 1;

	if( index < 0 ){
		MemoryPool_Extend( memPool );
		goto label;
	}

	return (void*)(memPool->pool[--memPool->numElementsFree]);
}
		
Bool MemoryPool_DeleteObject( MemoryPool *memPool, void *object )
{
	if( object != NULL ){
		int i = 0;
		int valid = 0;
		
		for ( i=0; i<memPool->numChunks; i++ ){
			int numObjs = 0;

			if( i==0 ){
				numObjs = memPool->numInitialElements;
			}
			else{
				numObjs = memPool->delta;
			}

			if(((char*)object >= memPool->elements[i]) &&
			( (char*)object < (memPool->elements[i] + memPool->elementSize * numObjs) ) ){
				valid = 1;
			}
		}

		if( valid ){
			memset( (char*)object, 0, memPool->elementSize );
			memPool->pool[memPool->numElementsFree++] = (char*)object;
			return 1;
		}
		else{
			return 0;
		}
	}
	else{
		return 0;
	}
}

void MemoryPool_Extend( MemoryPool *memPool )
{
	int i = 0;
	char **newPool;

	assert( memPool );

	memPool->numChunks++;

	memPool->elements = (char**)Memory_Realloc( memPool->elements, sizeof(void*)*memPool->numChunks );
	assert( memPool->elements );

	memPool->elements[memPool->numChunks-1] = Memory_Alloc_Bytes_Unnamed( memPool->elementSize * memPool->delta, "int" );
	memset( memPool->elements[memPool->numChunks-1], 0, memPool->elementSize * memPool->delta );

	newPool = Memory_Alloc_Bytes_Unnamed( sizeof(char*) * (memPool->numElements+memPool->delta), "char*" );
	assert( newPool );

	memcpy( newPool+memPool->delta, memPool->pool, sizeof(char*)*memPool->numElements );

	for( i=0; i<memPool->delta; i++ ){
		newPool[i] = &(memPool->elements[memPool->numChunks-1][i*memPool->elementSize]);
	}

	Memory_Free( memPool->pool );
	memPool->pool = newPool;
	memPool->numElements+=memPool->delta;
	memPool->numElementsFree=memPool->delta;
}

