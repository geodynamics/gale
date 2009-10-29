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
** $Id: Arrhenius.c 78 2005-11-29 11:58:21Z RobertTurnbull $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifdef HAVE_TIFF

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>

#include <glucifer/Base/Base.h>

#include "types.h"
#include "InputTIFF.h"

#include <assert.h>
#include <string.h>

#ifdef HAVE_TIFF
	#include <tiffio.h>
#endif

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type lucInputTIFF_Type = "lucInputTIFF";

/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
lucInputTIFF* _lucInputTIFF_New( 
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
		lucInputFormat_InputFunction*                      _Input,
		Name                                               name ) 
{
	lucInputTIFF*					self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( sizeOfSelf >= sizeof(lucInputTIFF) );
	self = (lucInputTIFF*) _lucInputFormat_New( 
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
			_Input,
			name );
	
	return self;
}

void _lucInputTIFF_Init( 
		lucInputTIFF*                                                self )
{
}

void _lucInputTIFF_Delete( void* InputFormat ) {
	lucInputTIFF*  self = (lucInputTIFF*)InputFormat;

	_lucInputFormat_Delete( self );
}

void _lucInputTIFF_Print( void* InputFormat, Stream* stream ) {
	lucInputTIFF*  self = (lucInputTIFF*)InputFormat;

	_lucInputFormat_Print( self, stream );
}

void* _lucInputTIFF_Copy( void* InputFormat, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap) {
	lucInputTIFF*  self = (lucInputTIFF*)InputFormat;
	lucInputTIFF* newInputFormat;

	newInputFormat = _lucInputFormat_Copy( self, dest, deep, nameExt, ptrMap );

	/* TODO */
	abort();

	return (void*) newInputFormat;
}


void* _lucInputTIFF_DefaultNew( Name name ) {
	return (void*) _lucInputTIFF_New(
		sizeof(lucInputTIFF),
		lucInputTIFF_Type,
		_lucInputTIFF_Delete,
		_lucInputTIFF_Print,
		NULL,
		_lucInputTIFF_DefaultNew,
		_lucInputTIFF_AssignFromXML,
		_lucInputTIFF_Build,
		_lucInputTIFF_Initialise,
		_lucInputTIFF_Execute,
		_lucInputTIFF_Destroy,
		_lucInputTIFF_Input,
		name );
}

void _lucInputTIFF_AssignFromXML( void* InputFormat, Stg_ComponentFactory* cf, void* data ){
	lucInputTIFF*  self = (lucInputTIFF*)InputFormat;

	/* Construct Parent */
	lucInputFormat_InitAll( self, "tiff" );

	_lucInputTIFF_Init( self );
}

void _lucInputTIFF_Build( void* InputFormat, void* data ) {}
void _lucInputTIFF_Initialise( void* InputFormat, void* data ) {}
void _lucInputTIFF_Execute( void* InputFormat, void* data ) {}
void _lucInputTIFF_Destroy( void* InputFormat, void* data ) {}

lucPixel* _lucInputTIFF_Input( void* inputFormat, Name imageName, Pixel_Index *width, Pixel_Index* height ){

	/* Using Sam Leffler's libtiff library 
	 * http://www.remotesensing.org/libtiff/ */
	TIFFRGBAImage  img;
	uint32*        raster;
	size_t         npixels;
	int            hasABGR = 0;
 	TIFF*          tif;
	char           emsg[1024];
	int            i;
	unsigned char* cp;
	lucPixel*      pixelData;
	
	lucInputTIFF*   self         = (lucInputTIFF*)inputFormat;
	
	tif = TIFFOpen(imageName, "r");
	Journal_Firewall( tif != NULL, Journal_MyStream( Error_Type, self ),
			"Error in func '%s' for %s '%s' - Cannot open '%s'\n", __func__, self->type, self->name, imageName );
  	
	if (TIFFRGBAImageBegin(&img, tif, 0,emsg)){
		npixels = img.width*img.height; 
		raster = (uint32 *)_TIFFmalloc(npixels*sizeof(uint32)); 
		if (raster != NULL){ 
			if (TIFFRGBAImageGet(&img, raster, img.width, img.height) == 0){ 
				TIFFError(imageName, emsg); 
				abort(); 
			} 
		} 
		TIFFRGBAImageEnd(&img); 
	} 
	else { 
		TIFFError(imageName, emsg); 
		abort();
	}
    
	self->imageWidth = img.width; 
	self->imageHeight = img.height; 
	*width = img.width;
	*height = img.height;
	pixelData = Memory_Alloc_Array( lucPixel, self->imageWidth * self->imageHeight, "pixel data" );
	
	/* code based upon http://www.opengl.org/developers/code/mjktips/libtiff/showtiff.c */
	/* If cannot directly display ABGR format, we need to reverse the component ordering in each pixel. :-( */
	if (!hasABGR) { 
		for (i = 0; i < npixels; i++) { 
			register unsigned char *cp = (unsigned char *) &raster[i]; 
			int t; 
			
			t = cp[3]; 
			cp[3] = cp[0]; 
			cp[0] = t; 
			t = cp[2]; 
			cp[2] = cp[1]; 
			cp[1] = t; 
		} 
	}
  
	for (i = 0; i < npixels; i++) {	  
		cp = (unsigned char *) &raster[i]; 
		pixelData[i][0] = cp[3] ;		
		pixelData[i][1]=  cp[2]; 
		pixelData[i][2] = cp[1];
	}
	
	return pixelData;
}

#endif

