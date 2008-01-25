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
** $Id: Mapping.c 3584 2006-05-16 11:11:07Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <mpi.h>

#include "Base/Foundation/Foundation.h"
#include "Base/IO/IO.h"

#include "Container.h"


/* Textual name of this class */
const Type	Mapping_Type = "Mapping";
unsigned long	Mapping_VacantIndex = (unsigned long)-1;
unsigned long	Mapping_CollidedIndex = (unsigned long)-2;


/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/

Mapping* Mapping_New( Name name ) {
	return _Mapping_New( sizeof(Mapping), 
			     Mapping_Type, 
			     _Mapping_Delete, 
			     _Mapping_Print, 
			     NULL );
}

Mapping* _Mapping_New( MAPPING_DEFARGS ) {
	Mapping*	self;

	/* Allocate memory */
	assert( sizeOfSelf >= sizeof(Mapping) );
	self = (Mapping*)_Stg_Class_New( STG_CLASS_PASSARGS );

	/* Virtual info */

	/* Mapping info */
	_Mapping_Init( self );

	return self;
}

void _Mapping_Init( Mapping* self ) {
	/* This class relies on the size of an unsigned long and a pointer being the same. */
	assert( sizeof(unsigned long) == sizeof(List*) );

	self->maxItems = 0;
	self->nItems = 0;
	self->keySize = 0;
	self->valSize = 0;

	self->hasher = NULL;
	self->items[0] = NULL;
	self->items[1] = NULL;
	self->states = NULL;
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _Mapping_Delete( void* mapping ) {
	Mapping*	self = (Mapping*)mapping;

	Mapping_Destruct( self );

	/* Delete the parent. */
	_Stg_Class_Delete( self );
}

void _Mapping_Print( void* mapping, Stream* stream ) {
	Mapping*	self = (Mapping*)mapping;
	
	/* Set the Journal for printing informations */
	Stream* mappingStream;
	mappingStream = Journal_Register( InfoStream_Type, "MappingStream" );

	/* Print parent */
	Journal_Printf( stream, "Mapping (ptr): (%p)\n", self );
	_Stg_Class_Print( self, stream );
}


/*--------------------------------------------------------------------------------------------------------------------------
** Public Functions
*/

void Mapping_SetTupleSizes( void* mapping, unsigned keySize, unsigned valueSize ) {
	Mapping*	self = (Mapping*)mapping;

	assert( self );

	Mapping_Destruct( self );

	self->keySize = keySize;
	self->valSize = valueSize;
	Mapping_InitHasher( self );
	Mapping_InitItems( self );
}

void Mapping_SetMaxSize( void* mapping, unsigned maxSize ) {
	Mapping*	self = (Mapping*)mapping;

	assert( self );

	Mapping_DestructData( self );

	self->maxItems = maxSize;
	Mapping_InitHasher( self );
	Mapping_InitItems( self );
}

void Mapping_Clear( void* mapping ) {
	Mapping*	self = (Mapping*)mapping;

	assert( self );

	self->nItems = 0;
}

void Mapping_Insert( void* mapping, void* key, void* value ) {
	Mapping*	self = (Mapping*)mapping;
	unsigned	hashInd;

	assert( self );
	assert( key );
	assert( value );

	hashInd = Hasher_Hash( self->hasher, key );
	if( self->states[hashInd] == Mapping_ItemState_Vacant ) {
		memcpy( self->items[0] + hashInd, key, self->keySize );
		memcpy( self->items[1] + hashInd, value, self->valSize );
		self->states[hashInd] = Mapping_ItemState_Occupied;
	}
	else
		Mapping_Collision( self, key, value, hashInd );

	self->nItems++;
}

void Mapping_Remove( void* mapping, void* key ) {
	Mapping*	self = (Mapping*)mapping;

	assert( self );
	assert( key );

	/* TODO */
	abort();
}

Bool Mapping_Map( void* mapping, void* key, void** value ) {
	Mapping*	self = (Mapping*)mapping;
	unsigned	hashInd;

	assert( self );
	assert( key );
	assert( value );

	hashInd = Hasher_Hash( self->hasher, key );
	if( self->states[hashInd] == Mapping_ItemState_Occupied ) {
		*value = (void*)(self->items[1] + hashInd);
		return True;
	}
	else if( self->states[hashInd] == Mapping_ItemState_Collision ) {
		List*		lists[2];
		unsigned	itm_i;

		lists[0] = *(List**)(self->items[0] + self->keySize * hashInd);
		lists[1] = *(List**)(self->items[1] + self->valSize * hashInd);

		for( itm_i = 0; itm_i < List_GetSize( lists[0] ); itm_i++ ) {
			if( memcmp( List_GetItem( lists[0], itm_i ), key, self->keySize ) ) {
				*value = List_GetItem( lists[1], itm_i );
				return True;
			}
		}
	}

	return False;
}

unsigned Mapping_GetSize( void* mapping ) {
	Mapping*	self = (Mapping*)mapping;

	assert( self );

	return self->nItems;
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Private Functions
*/

void Mapping_InitHasher( Mapping* self ) {
	if( self->keySize && self->maxItems ) {
		if( self->keySize == sizeof(unsigned) )
			self->hasher = (Hasher*)NumberHasher_New();
		else
			abort();

		Hasher_SetKeySize( self->hasher, self->keySize );
		Hasher_SetMaxItems( self->hasher, self->maxItems );
	}
}

void Mapping_InitItems( Mapping* self ) {
	if( self->maxItems && self->keySize && self->valSize ) {
		unsigned	tblSize;

		tblSize = Hasher_GetTableSize( self->hasher );
		self->items[0] = Memory_Alloc_Array_Bytes_Unnamed( self->keySize, tblSize, "" );
		self->items[1] = Memory_Alloc_Array_Bytes_Unnamed( self->valSize, tblSize, "" );
		self->states = Memory_Alloc_Array_Unnamed( Mapping_ItemState, tblSize );
		memset( self->states, 0, tblSize * sizeof(Mapping_ItemState) );
	}
}

void Mapping_Collision( Mapping* self, void* key, void* val, unsigned hashInd ) {
	List*	lists[2];

	if( self->states[hashInd] == Mapping_ItemState_Occupied ) {
		lists[0] = List_New();
		lists[1] = List_New();
		List_SetDelta( lists[0], 3 );
		List_SetDelta( lists[1], 3 );
		List_SetItemSize( lists[0], self->keySize );
		List_SetItemSize( lists[1], self->valSize );
		List_Append( lists[0], key );
		List_Append( lists[1], val );

		*(List**)(self->items[0] + self->keySize * hashInd) = lists[0];
		*(List**)(self->items[1] + self->valSize * hashInd) = lists[1];

		self->states[hashInd] = Mapping_ItemState_Collision;
	}
	else {
		lists[0] = *(List**)(self->items[0] + self->keySize * hashInd);
		lists[1] = *(List**)(self->items[1] + self->valSize * hashInd);

		List_Append( lists[0], key );
		List_Append( lists[1], val );
	}
}

void Mapping_DestructData( Mapping* self ) {
	if( self->hasher ) {
		unsigned	itm_i;

		for( itm_i = 0; itm_i < Hasher_GetTableSize( self->hasher ); itm_i++ ) {
			if( self->states[itm_i] != Mapping_ItemState_Collision )
				continue;

			FreeObject( self->items[0] + self->keySize * itm_i );
			FreeObject( self->items[1] + self->valSize * itm_i );
		}
	}

	KillArray( self->items[0] );
	KillArray( self->items[1] );
	KillArray( self->states );

	self->nItems = 0;
	self->maxItems = 0;
}

void Mapping_Destruct( Mapping* self ) {
	Mapping_DestructData( self );

	self->keySize = 0;
	self->valSize = 0;

	KillObject( self->hasher );
}
