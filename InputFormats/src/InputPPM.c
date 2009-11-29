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


#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>

#include <glucifer/Base/Base.h>

#include "types.h"
#include "InputPPM.h"

#include <assert.h>
#include <string.h>

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type lucInputPPM_Type = "lucInputPPM";

/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
lucInputPPM* _lucInputPPM_New(  LUCINPUTPPM_DEFARGS  ) 
{
	lucInputPPM*					self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( _sizeOfSelf >= sizeof(lucInputPPM) );
	self = (lucInputPPM*) _lucInputFormat_New(  LUCINPUTFORMAT_PASSARGS  );	
	
	return self;
}

void _lucInputPPM_Init( 
		lucInputPPM*                                                self )
{
}

void _lucInputPPM_Delete( void* InputFormat ) {
	lucInputPPM*  self = (lucInputPPM*)InputFormat;

	_lucInputFormat_Delete( self );
}

void _lucInputPPM_Print( void* InputFormat, Stream* stream ) {
	lucInputPPM*  self = (lucInputPPM*)InputFormat;

	_lucInputFormat_Print( self, stream );
}

void* _lucInputPPM_Copy( void* InputFormat, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap) {
	lucInputPPM*  self = (lucInputPPM*)InputFormat;
	lucInputPPM* newInputFormat;

	newInputFormat = _lucInputFormat_Copy( self, dest, deep, nameExt, ptrMap );

	/* TODO */
	abort();

	return (void*) newInputFormat;
}


void* _lucInputPPM_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                              _sizeOfSelf = sizeof(lucInputPPM);
	Type                                                      type = lucInputPPM_Type;
	Stg_Class_DeleteFunction*                              _delete = _lucInputPPM_Delete;
	Stg_Class_PrintFunction*                                _print = _lucInputPPM_Print;
	Stg_Class_CopyFunction*                                  _copy = NULL;
	Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _lucInputPPM_DefaultNew;
	Stg_Component_ConstructFunction*                    _construct = _lucInputPPM_AssignFromXML;
	Stg_Component_BuildFunction*                            _build = _lucInputPPM_Build;
	Stg_Component_InitialiseFunction*                  _initialise = _lucInputPPM_Initialise;
	Stg_Component_ExecuteFunction*                        _execute = _lucInputPPM_Execute;
	Stg_Component_DestroyFunction*                        _destroy = _lucInputPPM_Destroy;
	lucInputFormat_InputFunction*                           _input = _lucInputPPM_Input;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = ZERO;

	return (void*) _lucInputPPM_New(  LUCINPUTPPM_PASSARGS  );
}

void _lucInputPPM_AssignFromXML( void* InputFormat, Stg_ComponentFactory* cf, void* data ){
	lucInputPPM*  self = (lucInputPPM*)InputFormat;

	/* Construct Parent */
	lucInputFormat_InitAll( self, "ppm" );

	_lucInputPPM_Init( self );
}

void _lucInputPPM_Build( void* InputFormat, void* data ) {}
void _lucInputPPM_Initialise( void* InputFormat, void* data ) {}
void _lucInputPPM_Execute( void* InputFormat, void* data ) {}
void _lucInputPPM_Destroy( void* InputFormat, void* data ) {}

lucPixel* _lucInputPPM_Input( void* inputFormat, Name imageName, Pixel_Index *width, Pixel_Index* height ) {
	lucInputPPM*  self         = (lucInputPPM*)inputFormat;
	FILE*         imageFile;
	int           i,j;
	lucPixel*     pixelData;
	Bool          readTag = False;
	Bool          readWidth = False;
	Bool          readHeight = False;
	Bool          readColourCount = False;
	char*         charPtr;
	char          stringBuffer[241];
	int           ppmType;
	int           colourCount;
	
	imageFile = fopen( imageName, "r" );
	Journal_Firewall( imageFile != NULL, Journal_MyStream( Error_Type, self ),
			"Error in func '%s' for %s '%s' - Cannot open '%s'\n", __func__, self->type, self->name, imageName );

	while ( !readTag || ! readWidth || ! readHeight || !readColourCount ) {
		/* Read in a new line from file */
		charPtr = fgets( stringBuffer, 240, imageFile );
		assert ( charPtr );

		for ( charPtr = stringBuffer ; charPtr < stringBuffer + 240 ; charPtr++ ) {
			/* Check if we should go to a new line - this will happen for comments, line breaks and terminator characters */
			if ( *charPtr == '#' || *charPtr == '\n' || *charPtr == '\0' )
				break;

			/* Check if this is a space - if this is the case, then go to next line */
			if ( *charPtr == ' ' || *charPtr == '\t' )
				continue;

			if ( !readTag ) {
				sscanf( charPtr, "P%d", &ppmType );
				readTag = True;
			}
			else if ( !readWidth ) {
				sscanf( charPtr, "%u", width );
				readWidth = True;
			}
			else if ( !readHeight ) {
				sscanf( charPtr, "%u", height );
				readHeight = True;
			}
			else if ( !readColourCount ) {
				sscanf( charPtr, "%d", &colourCount );
				readColourCount = True;
			}

			/* Go to next white space */
			charPtr = strpbrk( charPtr, " \t" );

			/* If there are no more characters in line then go to next line */
			if ( charPtr == NULL )
				break;
		}
	}

	/* Only allow PPM images of type P6 and with 256 colours */
	assert( ppmType == 6 );
	assert( colourCount == 255 );

	/* Set width and height on object - TODO - Fix this hack, this shouldn't need to happen */
	self->imageWidth  = *width;
	self->imageHeight = *height;
	
	pixelData = Memory_Alloc_Array( lucPixel, self->imageWidth * self->imageHeight, "pixel data" );

	for ( j = self->imageHeight - 1 ; j >= 0 ; j--) {
	    for ( i = 0 ; i < self->imageWidth ; i++) {
		     fread( &pixelData[ self->imageWidth * j + i ], sizeof(lucPixel), 1, imageFile );
		}
	}				
	fclose( imageFile );
		
	return pixelData;
}


