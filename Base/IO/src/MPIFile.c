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
** $Id: MPIFile.c 2745 2005-03-05 08:12:18Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>

#include "Base/Foundation/Foundation.h"

#include "types.h"
#include "shortcuts.h"
#include "JournalFile.h"
#include "MPIFile.h"
#include "Stream.h"
#include "Journal.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>


const Type MPIFile_Type = "MPIFile";


JournalFile* MPIFile_New()
{
	/* Variables set in this function */
	SizeT                      _sizeOfSelf = sizeof(MPIFile);
	Type                              type = MPIFile_Type;
	Stg_Class_DeleteFunction*      _delete = _MPIFile_Delete;
	Stg_Class_PrintFunction*        _print = _MPIFile_Print;
	Stg_Class_CopyFunction*          _copy = NULL;

	return (JournalFile*)_MPIFile_New(  MPIFILE_PASSARGS  );
}

JournalFile* MPIFile_New2( const char* const fileName )
{
	JournalFile* result = MPIFile_New();

	if ( !JournalFile_Open( result, fileName ) )
	{
		/* File could not be opened successfully. Return cleanly. */
		Stg_Class_Delete( result );
		result = NULL;
	}
	
	return result;
}

MPIFile* _MPIFile_New(  MPIFILE_DEFARGS  )
{
	MPIFile* self;
	
	self = (MPIFile*)_Stg_Class_New(  STG_CLASS_PASSARGS  );
	
	_MPIFile_Init( self );
	
	return self;
}
	
void MPIFile_Init( MPIFile* self )
{
	/* Set virtual info. */
	self->_sizeOfSelf = sizeof(MPIFile);
	self->type = MPIFile_Type;
	self->_delete = _MPIFile_Delete;
	self->_print = _MPIFile_Print;
	self->_copy = NULL;
	
	_MPIFile_Init( self );
}
	
void _MPIFile_Init( MPIFile* self )
{
	_JournalFile_Init( (JournalFile*)self, _MPIFile_Open, _MPIFile_Append, _MPIFile_Close, _MPIFile_Flush );
}
	
void _MPIFile_Delete( void* cfile )
{
	MPIFile* self = (MPIFile*)cfile;
	
	_JournalFile_Delete( self );
}
void _MPIFile_Print( void* cfile, Stream* stream )
{
	MPIFile* self = (MPIFile*)cfile;
	
	_JournalFile_Print( self, stream );
}

	
Bool _MPIFile_Open( void* file, const char* const fileName )
{
	MPIFile*   self = (MPIFile*) file;
	int        fileOpenResult = 0;

	/* Remove the file */
        remove( fileName );

	fileOpenResult = MPI_File_open( MPI_COMM_WORLD, (char*)fileName,
				MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_EXCL, 
				MPI_INFO_NULL, &(self->mpiFile) );

	if (fileOpenResult != MPI_SUCCESS) {
		char         errorString[2000];
		int          errorStringLength = 0;
		Stream*      errorStream = Journal_Register( Error_Type, MPIFile_Type );
		int          myRank = 0;

		MPI_Comm_rank( MPI_COMM_WORLD, &myRank );
		MPI_Error_string( fileOpenResult, errorString, &errorStringLength);
		Journal_Printf( errorStream, "%3d: %s\n", myRank, errorString );
		MPI_Abort(MPI_COMM_WORLD, fileOpenResult );
	}

	self->fileHandle = &(self->mpiFile);
	return True;	
}
	
Bool _MPIFile_Append( void* file, const char* const fileName )
{
	MPIFile*   self = (MPIFile*) file;
	int        fileOpenResult = 0;

	fileOpenResult = MPI_File_open( MPI_COMM_WORLD, (char*)fileName, MPI_MODE_WRONLY, MPI_INFO_NULL,
		&(self->mpiFile) );

	if (fileOpenResult != MPI_SUCCESS) {
		char         errorString[2000];
		int          errorStringLength = 0;
		Stream*      errorStream = Journal_Register( Error_Type, MPIFile_Type );
		int          myRank = 0;

		MPI_Comm_rank( MPI_COMM_WORLD, &myRank );
		MPI_Error_string( fileOpenResult, errorString, &errorStringLength);
		Journal_Printf( errorStream, "%3d: %s\n", myRank, errorString );
		MPI_Abort(MPI_COMM_WORLD, fileOpenResult );
	}

	self->fileHandle = &(self->mpiFile);
	return True;	
}


Bool _MPIFile_Close( void* file )
{
	MPIFile* self = (MPIFile*) file;
	
	MPI_File_close( &(self->mpiFile) );
	
	return True;
}

Bool _MPIFile_Flush( void* file )
{

	/* No flush */
	
	return True;
}


