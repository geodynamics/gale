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
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/** \file
** <b>Role:</b>
**	Encapsulation of C-style printing.
**
** <b>Assumptions</b>
**	None
**
** <b>Comments</b>
**	None
**
** <b>Description</b>
**	A Wrapper class for C output functions such as fprintf() and fwrite().
**
** Comments:
**
** $Id:  $
**
**/

#ifndef __Base_IO_BinaryStream_h__
#define __Base_IO_BinaryStream_h__
	
	/** Textual name for BinaryStream class. */
	extern const Type BinaryStream_Type;
	
	
	/** \def __BinaryStream See BinaryStream. */
	#define __BinaryStream \
		/* General info */ \
		__Stream
	struct BinaryStream { __BinaryStream };


	/** Create a new BinaryStream */
	Stream* BinaryStream_New( Name name );

	/** Inits a BinaryStream. */
	void _BinaryStream_Init( BinaryStream* self );

	/** Constructor interface. */
	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define BINARYSTREAM_DEFARGS \
                STREAM_DEFARGS

	#define BINARYSTREAM_PASSARGS \
                STREAM_PASSARGS

	BinaryStream* _BinaryStream_New(  BINARYSTREAM_DEFARGS  );

	/** Stg_Class_Delete interface. */
	void _BinaryStream_Delete( void* cStream );
	
	/** Print interface. */
	void _BinaryStream_Print( void* cStream, Stream* stream );


	/** Printf() implementation. */
	SizeT _BinaryStream_Printf( Stream* stream, char *fmt, va_list args );
	
	/** Write() implementation. */
	SizeT _BinaryStream_Write( Stream* stream, void *data, SizeT elem_size, SizeT num_elems );
	
	/** Dump() implementation. Performs no operation for BinaryStreams. */
	Bool _BinaryStream_Dump( Stream* stream, void *data );
	
	/** SetFile() implementation. */
	Bool _BinaryStream_SetFile( Stream* stream, JournalFile* file );
	
	SizeT BinaryStream_WriteAllProcessors( Name filename, void *data, SizeT elem_size, SizeT num_elems, MPI_Comm communicator ) ;
	
#endif /* __IO_BinaryStreamFile_h__ */




