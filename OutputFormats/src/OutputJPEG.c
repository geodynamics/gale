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
** $Id: OutputJPEG.c 740 2007-10-11 08:05:31Z SteveQuenette $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifdef HAVE_JPEG

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>

#include <glucifer/Base/Base.h>

#include "types.h"
#include "OutputJPEG.h"

#include <assert.h>
#include <string.h>

#include <jpeglib.h>

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type lucOutputJPEG_Type = "lucOutputJPEG";

/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
lucOutputJPEG* _lucOutputJPEG_New( 
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
	lucOutputJPEG*					self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( sizeOfSelf >= sizeof(lucOutputJPEG) );
	self = (lucOutputJPEG*) _lucOutputFormat_New( 
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

void _lucOutputJPEG_Init( 
		lucOutputJPEG*                                                self,
		int                                                           quality )
{
	self->quality = quality;

	assert ( quality >= 0 && quality <= 100 );
}

void _lucOutputJPEG_Delete( void* outputFormat ) {
	lucOutputJPEG*  self = (lucOutputJPEG*)outputFormat;

	_lucOutputFormat_Delete( self );
}

void _lucOutputJPEG_Print( void* outputFormat, Stream* stream ) {
	lucOutputJPEG*  self = (lucOutputJPEG*)outputFormat;

	_lucOutputFormat_Print( self, stream );
}

void* _lucOutputJPEG_Copy( void* outputFormat, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap) {
	lucOutputJPEG*  self = (lucOutputJPEG*)outputFormat;
	lucOutputJPEG* newOutputFormat;

	newOutputFormat = _lucOutputFormat_Copy( self, dest, deep, nameExt, ptrMap );

	/* TODO */
	abort();

	return (void*) newOutputFormat;
}


void* _lucOutputJPEG_DefaultNew( Name name ) {
	return (void*) _lucOutputJPEG_New(
		sizeof(lucOutputJPEG),
		lucOutputJPEG_Type,
		_lucOutputJPEG_Delete,
		_lucOutputJPEG_Print,
		NULL,
		_lucOutputJPEG_DefaultNew,
		_lucOutputJPEG_Construct,
		_lucOutputJPEG_Build,
		_lucOutputJPEG_Initialise,
		_lucOutputJPEG_Execute,
		_lucOutputJPEG_Destroy,
		_lucOutputJPEG_Output,
		name );
}

void _lucOutputJPEG_Construct( void* outputFormat, Stg_ComponentFactory* cf, void* data ){
	lucOutputJPEG*  self = (lucOutputJPEG*)outputFormat;

	/* Construct Parent */
	lucOutputFormat_InitAll( self, "jpeg" );

	_lucOutputJPEG_Init( 
			self,
			Stg_ComponentFactory_GetInt( cf, self->name, "quality", 93 ) );
}

void _lucOutputJPEG_Build( void* outputFormat, void* data ) {}
void _lucOutputJPEG_Initialise( void* outputFormat, void* data ) {}
void _lucOutputJPEG_Execute( void* outputFormat, void* data ) {}
void _lucOutputJPEG_Destroy( void* outputFormat, void* data ) {}

void _lucOutputJPEG_Output( void* outputFormat, lucWindow* window, AbstractContext* context, lucPixel* pixelData ) {
	lucOutputJPEG*              self         = (lucOutputJPEG*) outputFormat;
	Pixel_Index                 width        = window->width;
	Pixel_Index                 height       = window->height;
	int                         rowStride    = (width * 3 + 3) & ~0x3;
	unsigned char*              pixels       = (unsigned char*) pixelData;
	FILE*                       file         = lucOutputFormat_OpenFile( self, window, context, "wb" );
	struct jpeg_compress_struct cinfo;
	struct jpeg_error_mgr       jerr;
    JSAMPROW                    row;

	cinfo.err = jpeg_std_error(&jerr);
	jpeg_create_compress(&cinfo);

	jpeg_stdio_dest(&cinfo, file);

	cinfo.image_width = width;
	cinfo.image_height = height;
	cinfo.input_components = 3;
	cinfo.in_color_space = JCS_RGB;

	jpeg_set_defaults(&cinfo);
	jpeg_set_quality(&cinfo, self->quality, TRUE);

	jpeg_start_compress(&cinfo, TRUE);

	while (cinfo.next_scanline < cinfo.image_height) {
		row = &pixels[rowStride * (cinfo.image_height - cinfo.next_scanline - 1)];
		jpeg_write_scanlines(&cinfo, &row, 1);
	}
	
	jpeg_finish_compress(&cinfo);
	fclose(file);
	jpeg_destroy_compress(&cinfo);	
}

#endif

