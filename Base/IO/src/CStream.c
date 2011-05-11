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
** $Id: CStream.c 3570 2006-05-15 03:28:02Z AlanLo $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include "Base/Foundation/Foundation.h"

#include "types.h"
#include "JournalFile.h"
#include "CFile.h"
#include "Stream.h"
#include "CStream.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include <stdarg.h>   /* Subsequent files need this for variable argument lists. */

#include "Journal.h"


const Type CStream_Type = "CStream";


Stream* CStream_New( Name name )
{
	/* Variables set in this function */
	SizeT                      _sizeOfSelf = sizeof(CStream);
	Type                              type = CStream_Type;
	Stg_Class_DeleteFunction*      _delete = _CStream_Delete;
	Stg_Class_PrintFunction*        _print = _CStream_Print;
	Stg_Class_CopyFunction*          _copy = _Stream_Copy;
	Stream_PrintfFunction*         _printf = _CStream_Printf;
	Stream_WriteFunction*           _write = _CStream_Write;
	Stream_DumpFunction*             _dump = _CStream_Dump;
	Stream_SetFileFunction*       _setFile = _CStream_SetFile;

	return (Stream*)_CStream_New(  CSTREAM_PASSARGS  );
}

void CStream_Init( CStream* self, Name name )
{
	
}


CStream* _CStream_New(  CSTREAM_DEFARGS  )
{
	CStream* self;
	
	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(CStream) );
	self = (CStream*)_Stream_New(  STREAM_PASSARGS  );
	
	_CStream_Init( self );
	
	return self;
}

void _CStream_Init( CStream* self )
{
	self->defaultFileType = CFile_New;

}
	
void _CStream_Delete( void* cStream )
{
	CStream* self = (CStream*)cStream;
	
	/* Stg_Class_Delete parent */
	_Stream_Delete( self );
}

void _CStream_Print( void* cStream, Stream* stream ) {

	CStream* self = (CStream*)cStream;
	
	/* General info */
	Journal_Printf( stream, "CStream (ptr): %p\n", cStream );
	
	/* Print parent */
	_Stream_Print( self, stream );
		
}
	
SizeT _CStream_Printf( Stream* stream, const char *fmt, va_list args )
{
	CStream* self = (CStream*)stream;
	SizeT    printResult;

	if ( self->_file == NULL )
	{
		return 0;
	}

	printResult = vfprintf( (FILE*) self->_file->fileHandle, fmt, args );
	return printResult;
}
	
SizeT _CStream_Write( Stream* stream, const void *data, SizeT elem_size, SizeT num_elems )
{
	CStream* self = (CStream*)stream;
	return fwrite( data, elem_size, num_elems, (FILE*) (self->_file->fileHandle) );
}
	
Bool _CStream_Dump( Stream* stream, const void *data )
{
	/* Traditional C does not have a dumping function. Hence, CStream performs no operation here. */
	return False;
}

Bool _CStream_SetFile( Stream* stream, JournalFile* file )
{
	if ( file->type == CFile_Type )
	{
		stream->_file = file;
		return True;
	}
	return False;
}





