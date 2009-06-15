/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
** Copyright (c) 2005, Monash Cluster Computing 
** All rights reserved.
** Redistribution and use in source and binary forms, with or without modification,
** are permitted provided that the following conditions are met:
**
** 		* Redistributions of source code must retain the above copyright notice, 
** 			this list of conditions and the following disclaimer.
** 		* Redistributions in binary form must reproduce the above copyright 
**			notice, this list of conditions and the following disclaimer in the 
**			documentation and/or other materials provided with the distribution.
** 		* Neither the name of the Monash University nor the names of its contributors 
**			may be used to endorse or promote products derived from this software 
**			without specific prior written permission.
**
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
** THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR 
** PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS 
** BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
** CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
** SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) 
** HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT 
** LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT 
** OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
**
**
** Contact:
*%		Cecile Duboz - Cecile.Duboz@sci.monash.edu.au
*%
** Contributors:
*+		Cecile Duboz
*+		Robert Turnbull
*+		Alan Lo
*+		Louis Moresi
*+		David Stegman
*+		David May
*+		Stevan Quenette
*+		Patrick Sunter
*+		Greg Watson
*+
** $Id: OutputPNG.c 740 2007-10-11 08:05:31Z SteveQuenette $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
#ifdef HAVE_LIBPNG

#include <mpi.h>
#include <png.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>

#include <glucifer/Base/Base.h>

#include "types.h"
#include "OutputPNG.h"

#include <assert.h>
#include <string.h>


/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type lucOutputPNG_Type = "lucOutputPNG";

/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
lucOutputPNG* _lucOutputPNG_New( 
		SizeT                                              sizeOfSelf,
		Type                                               type,
		Stg_Class_DeleteFunction*                          _delete,
		Stg_Class_PrintFunction*                           _print,
		Stg_Class_CopyFunction*                            _copy, 
		Stg_Component_DefaultConstructorFunction*          _defaultConstructor,
		Stg_Component_ConstructFunction*                   _construct,
		Stg_Component_BuildFunction*                       _build,
		Stg_Component_InitialiseFunction*                  _initialise,
		Stg_Component_ExecuteFunction*                     _execute,
		Stg_Component_DestroyFunction*                     _destroy,
		lucOutputFormat_OutputFunction*                    _output,
		Name                                               name ) 
{
	lucOutputPNG*					self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( sizeOfSelf >= sizeof(lucOutputPNG) );
	self = (lucOutputPNG*) _lucOutputFormat_New( 
			sizeOfSelf,
			type, 
			_delete,
			_print,
			_copy,
			_defaultConstructor,
			_construct,
			_build,
			_initialise,
			_execute,
			_destroy,
			_output,
			name );
	
	return self;
}

void _lucOutputPNG_Init( 
		lucOutputPNG*                                                self )
{
}

void _lucOutputPNG_Delete( void* outputFormat ) {
	lucOutputPNG*  self = (lucOutputPNG*)outputFormat;

	_lucOutputFormat_Delete( self );
}

void _lucOutputPNG_Print( void* outputFormat, Stream* stream ) {
	lucOutputPNG*  self = (lucOutputPNG*)outputFormat;

	_lucOutputFormat_Print( self, stream );
}

void* _lucOutputPNG_Copy( void* outputFormat, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap) {
	lucOutputPNG*  self = (lucOutputPNG*)outputFormat;
	lucOutputPNG* newOutputFormat;

	newOutputFormat = _lucOutputFormat_Copy( self, dest, deep, nameExt, ptrMap );

	/* TODO */
	abort();

	return (void*) newOutputFormat;
}


void* _lucOutputPNG_DefaultNew( Name name ) {
	return (void*) _lucOutputPNG_New(
		sizeof(lucOutputPNG),
		lucOutputPNG_Type,
		_lucOutputPNG_Delete,
		_lucOutputPNG_Print,
		NULL,
		_lucOutputPNG_DefaultNew,
		_lucOutputPNG_Construct,
		_lucOutputPNG_Build,
		_lucOutputPNG_Initialise,
		_lucOutputPNG_Execute,
		_lucOutputPNG_Destroy,
		_lucOutputPNG_Output,
		name );
}

void _lucOutputPNG_Construct( void* outputFormat, Stg_ComponentFactory* cf, void* data ){
	lucOutputPNG*  self = (lucOutputPNG*)outputFormat;

	/* Construct Parent */
	lucOutputFormat_InitAll( self, "png" );

	_lucOutputPNG_Init( self );
}

void _lucOutputPNG_Build( void* outputFormat, void* data ) {}
void _lucOutputPNG_Initialise( void* outputFormat, void* data ) {}
void _lucOutputPNG_Execute( void* outputFormat, void* data ) {}
void _lucOutputPNG_Destroy( void* outputFormat, void* data ) {}

/* Define png_jmpbuf() in case we are using a pre-1.0.6 version of libpng */
#ifndef png_jmpbuf
	#define png_jmpbuf(png_ptr) png_ptr->jmpbuf
#endif


void lucImagePNG_Write(png_structp png_ptr, png_bytep data, png_size_t length) {
	Stream* stream = (Stream*) png_get_io_ptr(png_ptr);
	
	Journal_Write( stream, (void*) data, 1, length );
}

void _lucOutputPNG_Output( void* outputFormat, lucWindow* window, AbstractContext* context, lucPixel* pixelData ) {
	lucOutputPNG* self       = (lucOutputPNG*) outputFormat;
	Pixel_Index   width        = window->width;
	Pixel_Index   height       = window->height;
	png_bytep     pixels       = (png_bytep) pixelData;
	int           rowStride    = (width * 3 + 3) & ~0x3;
	Stream*       stream       = lucOutputFormat_OpenStream( self, window, context );
	png_bytep*    row_pointers = Memory_Alloc_Array( png_bytep, height, "Row Pointers" );
	png_structp   pngWrite;
	png_infop     pngInfo;
	Pixel_Index   pixel_I;
	int           result;
	
	for ( pixel_I = 0 ; pixel_I < height ; pixel_I++ )
		row_pointers[pixel_I] = (png_bytep) &pixels[rowStride * (height - pixel_I - 1)];

	pngWrite = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
	Journal_Firewall( pngWrite != NULL, lucError, "Cannot create PNG write struct.\n" );

	pngInfo = png_create_info_struct(pngWrite);
	Journal_Firewall( pngInfo != NULL, lucError, "Cannot create PNG info struct.\n" );

	result = setjmp(png_jmpbuf(pngWrite));
	Journal_Firewall( result == 0, lucError, "In func %s: setjmp failed.\n", __func__ );

	png_set_write_fn(pngWrite, (void*) stream, lucImagePNG_Write, NULL);
	png_set_compression_level(pngWrite, Z_BEST_COMPRESSION);
	
	png_set_IHDR(pngWrite, pngInfo,
		width, height,
		8,
		PNG_COLOR_TYPE_RGB,
		PNG_INTERLACE_NONE,
		PNG_COMPRESSION_TYPE_DEFAULT,
		PNG_FILTER_TYPE_DEFAULT);

	png_write_info(pngWrite, pngInfo);
	
	png_write_image(pngWrite, row_pointers);
	png_write_end(pngWrite, pngInfo);

	/* Clean Up */
	png_destroy_info_struct(pngWrite, &pngInfo);
	png_destroy_write_struct(&pngWrite, NULL);
	Memory_Free( row_pointers );
	Stream_CloseFile( stream );	/* Release this file. Otherwise too many files will be opened at a time. */
}

#endif /* HAVE_LIBPNG */

