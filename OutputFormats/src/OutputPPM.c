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
** $Id: OutputPPM.c 740 2007-10-11 08:05:31Z SteveQuenette $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>

#include <glucifer/Base/Base.h>

#include "types.h"
#include "OutputPPM.h"

#include <assert.h>
#include <string.h>

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type lucOutputPPM_Type = "lucOutputPPM";

/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
lucOutputPPM* _lucOutputPPM_New(  LUCOUTPUTPPM_DEFARGS  ) 
{
	lucOutputPPM*					self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( _sizeOfSelf >= sizeof(lucOutputPPM) );
	self = (lucOutputPPM*) _lucOutputFormat_New(  LUCOUTPUTFORMAT_PASSARGS  );
	
	return self;
}

void _lucOutputPPM_Init( 
		lucOutputPPM*                                                self )
{
}

void _lucOutputPPM_Delete( void* outputFormat ) {
	lucOutputPPM*  self = (lucOutputPPM*)outputFormat;

	_lucOutputFormat_Delete( self );
}

void _lucOutputPPM_Print( void* outputFormat, Stream* stream ) {
	lucOutputPPM*  self = (lucOutputPPM*)outputFormat;

	_lucOutputFormat_Print( self, stream );
}

void* _lucOutputPPM_Copy( void* outputFormat, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap) {
	lucOutputPPM*  self = (lucOutputPPM*)outputFormat;
	lucOutputPPM* newOutputFormat;

	newOutputFormat = _lucOutputFormat_Copy( self, dest, deep, nameExt, ptrMap );

	/* TODO */
	abort();

	return (void*) newOutputFormat;
}


void* _lucOutputPPM_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                              _sizeOfSelf = sizeof(lucOutputPPM);
	Type                                                      type = lucOutputPPM_Type;
	Stg_Class_DeleteFunction*                              _delete = _lucOutputPPM_Delete;
	Stg_Class_PrintFunction*                                _print = _lucOutputPPM_Print;
	Stg_Class_CopyFunction*                                  _copy = NULL;
	Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _lucOutputPPM_DefaultNew;
	Stg_Component_ConstructFunction*                    _construct = _lucOutputPPM_AssignFromXML;
	Stg_Component_BuildFunction*                            _build = _lucOutputPPM_Build;
	Stg_Component_InitialiseFunction*                  _initialise = _lucOutputPPM_Initialise;
	Stg_Component_ExecuteFunction*                        _execute = _lucOutputPPM_Execute;
	Stg_Component_DestroyFunction*                        _destroy = _lucOutputPPM_Destroy;
	lucOutputFormat_OutputFunction*                        _output = _lucOutputPPM_Output;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = ZERO;

	return (void*) _lucOutputPPM_New(  LUCOUTPUTPPM_PASSARGS  );
}

void _lucOutputPPM_AssignFromXML( void* outputFormat, Stg_ComponentFactory* cf, void* data ){
	lucOutputPPM*  self = (lucOutputPPM*)outputFormat;

	/* Construct Parent */
	lucOutputFormat_InitAll( self, "ppm" );

	_lucOutputPPM_Init( self );
}

void _lucOutputPPM_Build( void* outputFormat, void* data ) {}
void _lucOutputPPM_Initialise( void* outputFormat, void* data ) {}
void _lucOutputPPM_Execute( void* outputFormat, void* data ) {}
void _lucOutputPPM_Destroy( void* outputFormat, void* data ) {}

void _lucOutputPPM_Output( void* outputFormat, lucWindow* window, AbstractContext* context, lucPixel* pixelData ) {
	lucOutputPPM*  self         = (lucOutputPPM*)outputFormat;
	Pixel_Index    windowWidth  = window->width;
	Pixel_Index    windowHeight = window->height;
	int            i, j;
	Stream*        stream;
	
	/* Open file */
	stream = lucOutputFormat_OpenStream( self, window, context );
	Journal_Firewall( stream != NULL, lucError, "Can't open file\n");

	/* Write header for PPM */
	Journal_Printf( stream, "P6\n%d %d\n255\n", windowWidth, windowHeight);
	
	/* Write RGB info */
	/* Top to bottom */
	for ( j = windowHeight - 1 ; j >= 0 ; j--) 
		for ( i = 0 ; i < windowWidth ; i++) 
			Journal_Write( stream,  &pixelData[ windowWidth * j + i ], sizeof(lucPixel), 1 );

	/* Bottom to Top */
	/* Journal_Write( stream, pixelData, sizeof(lucPixel), windowWidth * windowHeight ); */

	Stream_CloseFile( stream );
}


