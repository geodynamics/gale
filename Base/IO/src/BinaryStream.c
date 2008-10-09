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
** $Id:  $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>

#include "Base/Foundation/Foundation.h"

#include "types.h"
#include "shortcuts.h"
#include "JournalFile.h"
#include "CFile.h"
#include "Stream.h"
#include "BinaryStream.h"

#include "Base/IO/mpirecord/mpimessaging.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include <stdarg.h>   /* Subsequent files need this for variable argument lists. */

#include "Journal.h"


const Type BinaryStream_Type = "BinaryStream";


Stream* BinaryStream_New( Name name )
{
	return (Stream*)_BinaryStream_New( sizeof(BinaryStream), BinaryStream_Type, _BinaryStream_Delete, _BinaryStream_Print, _Stream_Copy, 
		name, _BinaryStream_Printf, _BinaryStream_Write, _BinaryStream_Dump, _BinaryStream_SetFile );
}

BinaryStream* _BinaryStream_New(
	SizeT			_sizeOfSelf, 
	Type			type, 
	Stg_Class_DeleteFunction*	_delete, 
	Stg_Class_PrintFunction* 	_print,
	Stg_Class_CopyFunction*	_copy, 
	Name			name,
	Stream_PrintfFunction*	_printf, 
	Stream_WriteFunction*	_write, 
	Stream_DumpFunction*	_dump,
	Stream_SetFileFunction*	_setFile )
{
	BinaryStream* self;
	
	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(BinaryStream) );
	self = (BinaryStream*)_Stream_New( _sizeOfSelf, type, _delete, _print, _copy, name, 
		_printf, _write, _dump, _setFile );
	
	_BinaryStream_Init( self );
	
	return self;
}

void _BinaryStream_Init( BinaryStream* self )
{
	self->defaultFileType = CFileBinary_New;
}
	
void _BinaryStream_Delete( void* cStream )
{
	BinaryStream* self = (BinaryStream*)cStream;
	
	/* Stg_Class_Delete parent */
	_Stream_Delete( self );
}

void _BinaryStream_Print( void* BinaryStream, Stream* stream ) {
	/* TODO? */
		
}
	
SizeT _BinaryStream_Printf( Stream* stream, char *fmt, va_list args )
{
	/* TODO? */
	return False;
}
	
SizeT _BinaryStream_Write( Stream* stream, void *data, SizeT elem_size, SizeT num_elems )
{
	BinaryStream* self = (BinaryStream*)stream;
	
	return fwrite( data, elem_size, num_elems, self->_file->fileHandle );
}
	
Bool _BinaryStream_Dump( Stream* stream, void *data )
{
	/* No specific dumping mechanism, can create in derived classes */
	return False;
}

Bool _BinaryStream_SetFile( Stream* stream, JournalFile* file )
{
	if ( file->type == CFile_Type )
	{
		stream->_file = file;
		return True;
	}
	return False;
}

SizeT BinaryStream_WriteAllProcessors( Name filename, void *data, SizeT elem_size, SizeT num_elems, MPI_Comm comm ) {
	Stream*    stream = Journal_Register( BinaryStream_Type, BinaryStream_Type );
	MPI_Status status;
	int        rank;
	int        nproc;
	int        confirmation = 0;
	const int         FINISHED_WRITING_TAG = 100;	
	MPI_Comm_rank( comm, &rank );
	MPI_Comm_size( comm, &nproc );

	/* wait for go-ahead from process ranked lower than me, to avoid competition writing to file */
	if ( rank != 0 ) {
		MPI_Recv( &confirmation, 1, MPI_INT, rank - 1, FINISHED_WRITING_TAG, comm, &status );
	}	
        /* open the file */
	if ( rank == 0 ) {
		Stream_RedirectFile( stream, filename );
	}
	else {
		Stream_AppendFile( stream, filename );
	}	

        /* write the data */
	Stream_Write( stream, data, elem_size, num_elems );

	/* close the file */
	Stream_CloseFile( stream);	
	/* send go-ahead from process ranked lower than me, to avoid competition writing to file */
	if ( rank != nproc - 1 ) {
		MPI_Ssend( &confirmation, 1, MPI_INT, rank + 1, FINISHED_WRITING_TAG, comm );
	}	

	return;
}
