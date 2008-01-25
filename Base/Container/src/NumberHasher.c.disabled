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
** $Id: NumberHasher.c 3584 2006-05-16 11:11:07Z PatrickSunter $
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
const Type NumberHasher_Type = "NumberHasher";


/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/

NumberHasher* NumberHasher_New( Name name ) {
	return _NumberHasher_New( sizeof(NumberHasher), 
				  NumberHasher_Type, 
				  _NumberHasher_Delete, 
				  _NumberHasher_Print, 
				  NULL, 
				  _NumberHasher_Hash, 
				  _NumberHasher_CalcTableSize );
}

NumberHasher* _NumberHasher_New( NUMBERHASHER_DEFARGS ) {
	NumberHasher*	self;

	/* Allocate memory */
	assert( sizeOfSelf >= sizeof(NumberHasher) );
	self = (NumberHasher*)_Hasher_New( HASHER_PASSARGS );

	/* Virtual info */

	/* NumberHasher info */
	_NumberHasher_Init( self );

	return self;
}

void _NumberHasher_Init( NumberHasher* self ) {
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _NumberHasher_Delete( void* hasher ) {
	NumberHasher*	self = (NumberHasher*)hasher;

	/* Delete the parent. */
	_Hasher_Delete( self );
}

void _NumberHasher_Print( void* hasher, Stream* stream ) {
	NumberHasher*	self = (NumberHasher*)hasher;
	
	/* Set the Journal for printing informations */
	Stream* hasherStream;
	hasherStream = Journal_Register( InfoStream_Type, "NumberHasherStream" );

	/* Print parent */
	Journal_Printf( stream, "NumberHasher (ptr): (%p)\n", self );
	_Hasher_Print( self, stream );
}

unsigned _NumberHasher_Hash( void* hasher, void* key ) {
	NumberHasher*	self = (NumberHasher*)hasher;

	assert( self );
	assert( self->keySize == sizeof(unsigned) );

	return *(unsigned*)key % self->tableSize;
}

void _NumberHasher_CalcTableSize( void* hasher ) {
	NumberHasher*	self = (NumberHasher*)hasher;

	assert( self );

	/*
	** Need to transform the maximum number of items into a number suitable for modulus on a number less than
	** the maximum items, probably some kind of prime.
	*/

	/* Expand the table size by a prescribed factor. */
	self->tableSize = (unsigned)((double)self->maxItems * 1.18);

	/* Make table size odd. */
	if( self->tableSize % 2 == 0 )
		self->tableSize++;
}


/*--------------------------------------------------------------------------------------------------------------------------
** Public Functions
*/

/*----------------------------------------------------------------------------------------------------------------------------------
** Private Functions
*/
