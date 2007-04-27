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
** $Id: Hasher.c 3584 2006-05-16 11:11:07Z PatrickSunter $
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
const Type Hasher_Type = "Hasher";


/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/

Hasher* _Hasher_New( HASHER_DEFARGS ) {
	Hasher*	self;

	/* Allocate memory */
	assert( sizeOfSelf >= sizeof(Hasher) );
	self = (Hasher*)_Stg_Class_New( STG_CLASS_PASSARGS );

	/* Virtual info */
	self->hashFunc = hashFunc;
	self->calcTableSizeFunc = calcTableSizeFunc;

	/* Hasher info */
	_Hasher_Init( self );

	return self;
}

void _Hasher_Init( Hasher* self ) {
	self->keySize = 0;
	self->maxItems = 0;
	self->tableSize = 0;
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _Hasher_Delete( void* hasher ) {
	Hasher*	self = (Hasher*)hasher;

	/* Delete the parent. */
	_Stg_Class_Delete( self );
}

void _Hasher_Print( void* hasher, Stream* stream ) {
	Hasher*	self = (Hasher*)hasher;
	
	/* Set the Journal for printing informations */
	Stream* hasherStream;
	hasherStream = Journal_Register( InfoStream_Type, "HasherStream" );

	/* Print parent */
	Journal_Printf( stream, "Hasher (ptr): (%p)\n", self );
	_Stg_Class_Print( self, stream );
}

void _Hasher_CalcTableSize( void* hasher ) {
	Hasher*	self = (Hasher*)hasher;

	assert( self );

	self->tableSize = self->maxItems;
}


/*--------------------------------------------------------------------------------------------------------------------------
** Public Functions
*/

void Hasher_SetKeySize( void* hasher, unsigned keySize ) {
	Hasher*	self = (Hasher*)hasher;

	assert( self );

	self->keySize = keySize;
}

void Hasher_SetMaxItems( void* hasher, unsigned maxItems ) {
	Hasher*	self = (Hasher*)hasher;

	assert( self );

	self->maxItems = maxItems;
	Hasher_CalcTableSize( self );
}

unsigned Hasher_GetTableSize( void* hasher ) {
	Hasher*	self = (Hasher*)hasher;

	assert( self );

	return self->tableSize;
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Private Functions
*/
