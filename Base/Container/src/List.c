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
** $Id: List.c 3584 2006-05-16 11:11:07Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <mpi.h>

#include "Base/Foundation/Foundation.h"
#include "Base/IO/IO.h"

#include "types.h"
#include "List.h"


/* Textual name of this class */
const Type List_Type = "List";


/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/

List* List_New( Name name ) {
	/* Variables set in this function */
	SizeT                      _sizeOfSelf = sizeof(List);
	Type                              type = List_Type;
	Stg_Class_DeleteFunction*      _delete = _List_Delete;
	Stg_Class_PrintFunction*        _print = _List_Print;
	Stg_Class_CopyFunction*          _copy = NULL;

	return _List_New(  LIST_PASSARGS  );
}

List* _List_New(  LIST_DEFARGS  ) {
	List*	self;
	
	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(List) );
	self = (List*)_Stg_Class_New(  STG_CLASS_PASSARGS  );

	/* Virtual info */

	/* List info */
	_List_Init( self );

	return self;
}

void _List_Init( List* self ) {
	self->nItems = 0;
	self->items = NULL;
	self->itemSize = 0;
	self->maxItems = 0;
	self->delta = 10;
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _List_Delete( void* list ) {
	List*	self = (List*)list;

	List_Destruct( self );

	/* Delete the parent. */
	_Stg_Class_Delete( self );
}

void _List_Print( void* list, Stream* stream ) {
	List*	self = (List*)list;
	
	/* Set the Journal for printing informations */
	Stream* listStream;
	listStream = Journal_Register( InfoStream_Type, "ListStream" );

	/* Print parent */
	Journal_Printf( stream, "List (ptr): (%p)\n", self );
	_Stg_Class_Print( self, stream );
}


/*--------------------------------------------------------------------------------------------------------------------------
** Public Functions
*/

void List_SetDelta( void* list, unsigned delta ) {
	List*	self = (List*)list;

	assert( self );
	assert( delta );

	self->delta = delta;
}

void List_SetItemSize( void* list, unsigned itemSize ) {
	List*	self = (List*)list;

	assert( self );
	assert( itemSize );

	List_Destruct( self );

	self->itemSize = itemSize;
}

void List_Clear( void* list ) {
	List*	self = (List*)list;

	assert( self );

	self->nItems = 0;
	KillArray( self->items );
	self->maxItems = 0;
}

void List_Insert( void* list, unsigned index, void* data ) {
	List*		self = (List*)list;
	unsigned	item_i;

	assert( self );
	assert( index <= self->nItems );
	assert( data );
	assert( self->itemSize );

	if( self->nItems == self->maxItems )
		List_Expand( self );

	for( item_i = self->nItems; item_i > index; item_i-- )
		memcpy( self->items + self->itemSize * item_i, 
			self->items + self->itemSize * (item_i - 1), 
			self->itemSize );

	memcpy( self->items + self->itemSize * index, data, self->itemSize );

	self->nItems++;
}

void List_Append( void* list, void* data ) {
	List*	self = (List*)list;

	assert( self );

	List_Insert( self, self->nItems, data );
}

void List_Prepend( void* list, void* data ) {
	List_Insert( list, 0, data );
}

void List_Remove( void* list, void* data ) {
	List*		self = (List*)list;
	unsigned	item_i;

	assert( self );
	assert( data );

	for( item_i = 0; item_i < self->nItems; item_i++ ) {
		if( !memcmp( self->items + self->itemSize * item_i, data, self->itemSize ) )
			break;
	}
	assert( item_i < self->nItems );

	for( item_i++; item_i < self->nItems; item_i++ ) {
		memcpy( self->items + self->itemSize * (item_i - 1), 
			self->items + self->itemSize * item_i, 
			self->itemSize );
	}

	if( --self->nItems % self->delta == 0 )
		List_Contract( self );
}

void* List_GetItem( void* list, unsigned index ) {
	List*	self = (List*)list;

	assert( self );
	assert( index < self->nItems );

	return self->items + self->itemSize * index;
}

unsigned List_GetSize( void* list ) {
	List*	self = (List*)list;

	assert( self );

	return self->nItems;
}

Bool List_Exists( void* list, void* data ) {
	List*		self = (List*)list;
	unsigned	item_i;

	assert( self );

	for( item_i = 0; item_i < self->nItems; item_i++ ) {
		if( !memcmp( self->items + self->itemSize * item_i, data, self->itemSize ) )
			break;
	}

	return (item_i < self->nItems) ? True : False;
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Private Functions
*/

void List_Expand( List* self ) {
	self->maxItems += self->delta;
	if( !self->items )
		self->items = Memory_Alloc_Array_Bytes( self->itemSize, self->maxItems, "", "List::items" );
	else
		self->items = Memory_Realloc_Array_Bytes( self->items, self->itemSize, self->maxItems );
}

void List_Contract( List* self ) {
	if( self->delta > self->maxItems )
		self->maxItems = 0;
	else
		self->maxItems -= self->delta;

	if( !self->maxItems )
		KillArray( self->items );
	else
		self->items = Memory_Realloc_Array_Bytes( self->items, self->itemSize, self->maxItems );
}

void List_Destruct( List* self ) {
	self->nItems = 0;
	KillArray( self->items );
	self->itemSize = 0;
	self->maxItems = 0;
}


