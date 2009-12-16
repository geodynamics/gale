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
lucOutputPNG* _lucOutputPNG_New(  LUCOUTPUTPNG_DEFARGS  ) 
{
	lucOutputPNG*					self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( _sizeOfSelf >= sizeof(lucOutputPNG) );
	self = (lucOutputPNG*) _lucOutputFormat_New(  LUCOUTPUTFORMAT_PASSARGS  );
	
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
	/* Variables set in this function */
	SizeT                                              _sizeOfSelf = sizeof(lucOutputPNG);
	Type                                                      type = lucOutputPNG_Type;
	Stg_Class_DeleteFunction*                              _delete = _lucOutputPNG_Delete;
	Stg_Class_PrintFunction*                                _print = _lucOutputPNG_Print;
	Stg_Class_CopyFunction*                                  _copy = NULL;
	Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _lucOutputPNG_DefaultNew;
	Stg_Component_ConstructFunction*                    _construct = _lucOutputPNG_AssignFromXML;
	Stg_Component_BuildFunction*                            _build = _lucOutputPNG_Build;
	Stg_Component_InitialiseFunction*                  _initialise = _lucOutputPNG_Initialise;
	Stg_Component_ExecuteFunction*                        _execute = _lucOutputPNG_Execute;
	Stg_Component_DestroyFunction*                        _destroy = _lucOutputPNG_Destroy;
	lucOutputFormat_OutputFunction*                        _output = _lucOutputPNG_Output;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

	return (void*) _lucOutputPNG_New(  LUCOUTPUTPNG_PASSARGS  );
}

void _lucOutputPNG_AssignFromXML( void* outputFormat, Stg_ComponentFactory* cf, void* data ){
	lucOutputPNG*  self = (lucOutputPNG*)outputFormat;

	/* Construct Parent */
   self->extension = "png";
	_lucOutputFormat_AssignFromXML( outputFormat, cf, data);

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

void _lucOutputPNG_Output( void* outputFormat, lucWindow* window, AbstractContext* context, void* pixelData ) {
	lucOutputPNG* self       = (lucOutputPNG*) outputFormat;
	Pixel_Index   width        = window->width;
	Pixel_Index   height       = window->height;
	png_bytep     pixels       = (png_bytep) pixelData;
	int           rowStride    = width * 3; /* Don't pad lines! pack alignment is set to 1 */
	Stream*       stream       = lucOutputFormat_OpenStream( self, window, context );
	png_bytep*    row_pointers = Memory_Alloc_Array( png_bytep, height, "Row Pointers" );
	png_structp   pngWrite;
	png_infop     pngInfo;
	Pixel_Index   pixel_I;
	int           result;
   int colour_type = PNG_COLOR_TYPE_RGB; 
	
   if (self->transparent) {
	   rowStride    = width * 4; /* Don't pad lines! pack alignment is set to 1 */
      colour_type = PNG_COLOR_TYPE_RGBA;
   }

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
		colour_type,
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

#endif /* HAVE_PNG */



