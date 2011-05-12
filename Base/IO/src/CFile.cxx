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
** $Id: CFile.c 3462 2006-02-19 06:53:24Z WalterLandry $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include "Base/Foundation/Foundation.h"

#include "types.h"
#include "shortcuts.h"
#include "JournalFile.h"
#include "CFile.h"
#include "Stream.h"
#include "Journal.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>


const Type CFile_Type = "CFile";


JournalFile* CFile_New()
{
	/* Variables set in this function */
	SizeT                      _sizeOfSelf = sizeof(CFile);
	Type                              type = CFile_Type;
	Stg_Class_DeleteFunction*      _delete = _CFile_Delete;
	Stg_Class_PrintFunction*        _print = _CFile_Print;
	Stg_Class_CopyFunction*          _copy = NULL;
	Bool                            binary = False;

	return (JournalFile*)_CFile_New(  CFILE_PASSARGS  );
}

JournalFile* CFileBinary_New()
{
	/* Variables set in this function */
	SizeT                      _sizeOfSelf = sizeof(CFile);
	Type                              type = CFile_Type;
	Stg_Class_DeleteFunction*      _delete = _CFile_Delete;
	Stg_Class_PrintFunction*        _print = _CFile_Print;
	Stg_Class_CopyFunction*          _copy = NULL;
	Bool                            binary = True;

	return (JournalFile*)_CFile_New(  CFILE_PASSARGS  );
}


JournalFile* CFile_New2( Name fileName )
{
	JournalFile* result = CFile_New();

	if ( !JournalFile_Open( result, fileName ) )
	{
		/* File could not be opened successfully. Return cleanly. */
		Stg_Class_Delete( result );
		result = NULL;
	}
	
	return result;
}

CFile* _CFile_New(  CFILE_DEFARGS  )
{
	CFile* self;
	
	self = (CFile*)_Stg_Class_New(  STG_CLASS_PASSARGS  );
	self->binary = binary;
	
	_CFile_Init( self );
	
	return self;
}
	
void CFile_Init( CFile* self )
{
	/* Set virtual info. */
	self->_sizeOfSelf = sizeof(CFile);
	self->type = CFile_Type;
	self->_delete = _CFile_Delete;
	self->_print = _CFile_Print;
	self->_copy = NULL;
	
	_CFile_Init( self );
}
	
void _CFile_Init( CFile* self )
{
	_JournalFile_Init( (JournalFile*)self, (JournalFile_OpenFunction*)_CFile_Open,
		(JournalFile_AppendFunction*)_CFile_Append, _CFile_Close, _CFile_Flush );
}
	
void _CFile_Delete( void* cfile )
{
	CFile* self = (CFile*)cfile;
	
	_JournalFile_Delete( self );
}
void _CFile_Print( void* cfile, Stream* stream )
{
	CFile* self = (CFile*)cfile;
	
	_JournalFile_Print( self, stream );
}

	
Bool _CFile_Open( void* file, Name const fileName )
{
	CFile* self = (CFile*) file;
	FILE* filePtr;
	
	if(!self->binary)
		filePtr = fopen( fileName, "w" );
	else
		filePtr = fopen( fileName, "wb" );
	
	if ( filePtr == NULL )
	{
		return False;
	}

	self->fileHandle = (void*) filePtr;
	
	return True;	
}
	
Bool _CFile_Append( void* file, Name const fileName )
{
	CFile* self = (CFile*) file;
	FILE* filePtr;
	
	if(!self->binary)
		filePtr = fopen( fileName, "a" );
	else
		filePtr = fopen( fileName, "ab" );
	
	if ( filePtr == NULL )
	{
		return False;
	}

	self->fileHandle = (void*) filePtr;
	
	return True;	
}
Bool _CFile_Close( void* file )
{
	CFile* self = (CFile*) file;
	if ( self->fileHandle != NULL )
	{
          return fclose(  (FILE*) self->fileHandle ) == 0 ? True : False;
	}
	return False;
}

Bool _CFile_Flush( void* file )
{
	CFile* self = (CFile*) file;
	if ( self->fileHandle != NULL )
	{
          return fflush( (FILE*) self->fileHandle ) == 0 ? True : False;
	}
	return False;
}



