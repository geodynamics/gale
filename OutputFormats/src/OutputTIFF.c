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
** $Id: OutputTIFF.c 628 2006-10-12 08:23:07Z SteveQuenette $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifdef HAVE_TIFF

#include <mpi.h>
#include <StGermain/StGermain.h>

#include <glucifer/Base/Base.h>

#include "types.h"
#include "OutputTIFF.h"

#include <assert.h>
#include <string.h>

#include <tiff.h>
#include <tiffio.h>

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type lucOutputTIFF_Type = "lucOutputTIFF";

/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
lucOutputTIFF* _lucOutputTIFF_New( 
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
	lucOutputTIFF*					self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( sizeOfSelf >= sizeof(lucOutputTIFF) );
	self = (lucOutputTIFF*) _lucOutputFormat_New( 
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

void _lucOutputTIFF_Init( 
		lucOutputTIFF*                                                self )
{
}

void _lucOutputTIFF_Delete( void* outputFormat ) {
	lucOutputTIFF*  self = (lucOutputTIFF*)outputFormat;

	_lucOutputFormat_Delete( self );
}

void _lucOutputTIFF_Print( void* outputFormat, Stream* stream ) {
	lucOutputTIFF*  self = (lucOutputTIFF*)outputFormat;

	_lucOutputFormat_Print( self, stream );
}

void* _lucOutputTIFF_Copy( void* outputFormat, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap) {
	lucOutputTIFF*  self = (lucOutputTIFF*)outputFormat;
	lucOutputTIFF* newOutputFormat;

	newOutputFormat = _lucOutputFormat_Copy( self, dest, deep, nameExt, ptrMap );

	/* TODO */
	abort();

	return (void*) newOutputFormat;
}


void* _lucOutputTIFF_DefaultNew( Name name ) {
	return (void*) _lucOutputTIFF_New(
		sizeof(lucOutputTIFF),
		lucOutputTIFF_Type,
		_lucOutputTIFF_Delete,
		_lucOutputTIFF_Print,
		NULL,
		_lucOutputTIFF_DefaultNew,
		_lucOutputTIFF_Construct,
		_lucOutputTIFF_Build,
		_lucOutputTIFF_Initialise,
		_lucOutputTIFF_Execute,
		_lucOutputTIFF_Destroy,
		_lucOutputTIFF_Output,
		name );
}

void _lucOutputTIFF_Construct( void* outputFormat, Stg_ComponentFactory* cf, void* data ){
	lucOutputTIFF*  self = (lucOutputTIFF*)outputFormat;

	/* Construct Parent */
	lucOutputFormat_InitAll( self, "tiff" );

	_lucOutputTIFF_Init( self );
}

void _lucOutputTIFF_Build( void* outputFormat, void* data ) {}
void _lucOutputTIFF_Initialise( void* outputFormat, void* data ) {}
void _lucOutputTIFF_Execute( void* outputFormat, void* data ) {}
void _lucOutputTIFF_Destroy( void* outputFormat, void* data ) {}

void _lucOutputTIFF_Output( void* outputFormat, lucWindow* window, AbstractContext* context, lucPixel* pixelData ) {
	lucOutputTIFF*              self         = (lucOutputTIFF*) outputFormat;
	Pixel_Index                 width        = window->width;
	Pixel_Index                 height       = window->height;
	Pixel_Index                 line_I;
	TIFF*                       file;
	Name                        filename;
	lucPixel*                   linePtr;

	/* Open File */
	filename = lucOutputFormat_GetImageFilename( self, window, context );
	file = TIFFOpen(filename, "w");
	Journal_Firewall( file != NULL, lucError, "Cannot Open File %s\n", filename );
	Memory_Free( filename );

	TIFFSetField(file, TIFFTAG_IMAGEWIDTH, (uint32) width);
	TIFFSetField(file, TIFFTAG_IMAGELENGTH, (uint32) height);
	TIFFSetField(file, TIFFTAG_BITSPERSAMPLE, 8);
	TIFFSetField(file, TIFFTAG_COMPRESSION, COMPRESSION_PACKBITS );
	TIFFSetField(file, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
	TIFFSetField(file, TIFFTAG_SAMPLESPERPIXEL, 3);
	TIFFSetField(file, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
	TIFFSetField(file, TIFFTAG_ROWSPERSTRIP, 1);
	TIFFSetField(file, TIFFTAG_IMAGEDESCRIPTION, window->name );
	
	linePtr = pixelData;
	for ( line_I = height - 1;  line_I != (Pixel_Index) -1 ;  line_I--) {
		if (TIFFWriteScanline(file, linePtr, line_I, 0) < 0) {
			TIFFClose(file);
			return;
		}
		linePtr = (lucPixel*)((ArithPointer)linePtr + (ArithPointer)(width * sizeof(lucPixel)));
	}
	TIFFClose(file);
}

#endif

