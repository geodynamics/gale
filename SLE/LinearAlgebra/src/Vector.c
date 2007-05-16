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
** $Id: Vector.c 3584 2006-05-16 11:11:07Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <mpi.h>
#include "StGermain/StGermain.h"
#include "Discretisation/Discretisation.h"
#include "LinearAlgebra.h"


/* Textual name of this class */
const Type Vector_Type = "Vector";


/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/

Vector* _Vector_New( VECTOR_DEFARGS ) {
	Vector*	self;

	/* Allocate memory */
	assert( sizeOfSelf >= sizeof(Vector) );
	self = (Vector*)_Stg_Component_New( STG_COMPONENT_PASSARGS );

	/* Virtual info */
	self->setCommFunc = setCommFunc;
	self->setGlobalSizeFunc = setGlobalSizeFunc;
	self->setLocalSizeFunc = setLocalSizeFunc;
	self->addEntriesFunc = addEntriesFunc;
	self->insertEntriesFunc = insertEntriesFunc;
	self->setScalarFunc = setScalarFunc;
	self->zeroFunc = zeroFunc;
	self->assemblyBeginFunc = assemblyBeginFunc;
	self->assemblyEndFunc = assemblyEndFunc;

	self->addFunc = addFunc;
	self->addScaledFunc = addScaledFunc;
	self->scaleAddFunc = scaleAddFunc;
	self->subtractFunc = subtractFunc;
	self->scaleFunc = scaleFunc;
	self->dotProductFunc = dotProductFunc;
	self->l2NormFunc = l2NormFunc;
	self->lInfNormFunc = lInfNormFunc; 
	self->pointwiseMultiplyFunc = pointwiseMultiplyFunc;
	self->pointwiseDivideFunc = pointwiseDivideFunc;
	self->reciprocalFunc = reciprocalFunc;

	self->getSizeFunc = getSizeFunc;
	self->getLocalSizeFunc = getLocalSizeFunc;
	self->getArrayFunc = getArrayFunc;
	self->restoreArrayFunc = restoreArrayFunc;
	self->duplicateFunc = duplicateFunc;
	self->copyEntriesFunc = copyEntriesFunc;
	self->viewFunc = viewFunc;

	/* Vector info */
	_Vector_Init( self );

	return self;
}

void _Vector_Init( Vector* self ) {
	assert( self && Stg_CheckType( self, Vector ) );

	self->comm = MPI_COMM_WORLD;
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _Vector_Delete( void* vector ) {
	Vector*	self = (Vector*)vector;

	assert( self && Stg_CheckType( self, Vector ) );

	/* Delete the parent. */
	_Stg_Component_Delete( self );
}

void _Vector_Print( void* vector, Stream* stream ) {
	Vector*	self = (Vector*)vector;
	
	/* Set the Journal for printing informations */
	Stream* vectorStream;
	vectorStream = Journal_Register( InfoStream_Type, "VectorStream" );

	assert( self && Stg_CheckType( self, Vector ) );

	/* Print parent */
	Journal_Printf( stream, "Vector (ptr): (%p)\n", self );
	_Stg_Component_Print( self, stream );
}

void _Vector_Construct( void* vector, Stg_ComponentFactory* cf, void* data ) {
	Vector*	self = (Vector*)vector;

	assert( self && Stg_CheckType( self, Vector ) );
	assert( cf );
}

void _Vector_Build( void* vector, void* data ) {
}

void _Vector_Initialise( void* vector, void* data ) {
}

void _Vector_Execute( void* vector, void* data ) {
}

void _Vector_Destroy( void* vector, void* data ) {
}

void _Vector_SetComm( void* vector, MPI_Comm comm ) {
	Vector*	self = (Vector*)vector;

	assert( self && Stg_CheckType( self, Vector ) );

	self->comm = comm;
}

void _Vector_View( void* vector, void* _stream ) {
	Vector*		self = (Vector*)vector;
	Stream*		stream = (Stream*)_stream;
	unsigned	size;
	double*		array;
	unsigned	entry_i;

	assert( self && Stg_CheckType( self, Vector ) );

	if( !stream )
		stream = Journal_Register( Info_Type, "tmp" );

	size = Vector_GetLocalSize( self );
	Vector_GetArray( self, &array );
	Journal_Printf( stream, "%s = [\n", self->name );
	for( entry_i = 0; entry_i < size; entry_i++ ) 
		Journal_Printf( stream , "\t%u: \t %.12g\n", entry_i, array[entry_i] );
	Journal_Printf( stream, "];\n");
	Vector_RestoreArray( self, &array );
}


/*--------------------------------------------------------------------------------------------------------------------------
** Public Functions
*/

void Vector_Dump( void* vector, const char* filename ) {
	Vector*		self = (Vector*)vector;
	FILE*		fp;
	unsigned	size;
	double*		array;
	unsigned	entry_i, b_i;

	assert( self && Stg_CheckType( self, Vector ) );

	insist( fp = fopen( filename, "w" ), != NULL );

	size = Vector_GetLocalSize( self );
	Vector_GetArray( self, &array );
	for( entry_i = 0; entry_i < size; entry_i++ ) {
		for( b_i = 0; b_i < sizeof(double) / sizeof(unsigned); b_i++ )
			fprintf( fp, "%x", ((unsigned*)(array + entry_i))[b_i] );
		fprintf( fp, "\n" );
	}
	Vector_RestoreArray( self, &array );

	fclose( fp );
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Private Functions
*/
