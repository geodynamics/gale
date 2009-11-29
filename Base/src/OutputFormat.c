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
** $Id: OutputFormat.c 740 2007-10-11 08:05:31Z SteveQuenette $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>

#include "types.h"
#include "OutputFormat.h"
#include "ColourMap.h"
#include "Window.h"
#include "DrawingObject_Register.h"
#include "DrawingObject.h"
#include "Camera.h"
#include "Init.h"
#include "Window.h"

#include <assert.h>

#ifndef MASTER
	#define MASTER 0
#endif

const Type lucOutputFormat_Type = "lucOutputFormat";

lucOutputFormat* _lucOutputFormat_New(  LUCOUTPUTFORMAT_DEFARGS  )
{
	lucOutputFormat*    self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( _sizeOfSelf >= sizeof(lucOutputFormat) );
	/* The following terms are parameters that have been passed into this function but are being set before being passed onto the parent */
	/* This means that any values of these parameters that are passed into this function are not passed onto the parent function
	   and so should be set to ZERO in any children of this class. */
	nameAllocationType = NON_GLOBAL;

	self = (lucOutputFormat*) _Stg_Component_New(  STG_COMPONENT_PASSARGS  );

	self->_output = _output;

	return self;
}

void _lucOutputFormat_Init( 
		lucOutputFormat*                                   self, 
		Name                                               extension )
{
	self->extension     = StG_Strdup( extension );
}

void lucOutputFormat_InitAll( 
		void*                                              outputFormat,
		Name                                               extension )
{
	lucOutputFormat* self        = outputFormat;

	_lucOutputFormat_Init( self, extension );
}

	
void _lucOutputFormat_Delete( void* outputFormat ) {
	lucOutputFormat* self        = outputFormat;

	Memory_Free( self->extension );
	
	_Stg_Component_Delete( self );
}

void _lucOutputFormat_Print( void* outputFormat, Stream* stream ) {
	lucOutputFormat*          self        = outputFormat;
	
	Journal_Printf( stream, "lucOutputFormat: %s\n", self->name );

	/* Print Parent */
	_Stg_Component_Print( self, stream );

	Journal_PrintString( stream, self->extension );
}

void* _lucOutputFormat_Copy( void* outputFormat, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	lucOutputFormat* self        = outputFormat;
	lucOutputFormat* newOutputFormat;

	newOutputFormat = _Stg_Component_Copy( self, dest, deep, nameExt, ptrMap );

	newOutputFormat->extension = StG_Strdup( self->extension );

	return (void*) newOutputFormat;
}

void _lucOutputFormat_AssignFromXML( void* outputFormat, Stg_ComponentFactory* cf, void* data ) {
	lucOutputFormat* self        = outputFormat;

	self->context = Stg_ComponentFactory_ConstructByKey( cf, self->name, "Context", AbstractContext, False, data );
	if( !self->context ) 
		self->context = Stg_ComponentFactory_ConstructByName( cf, "context", AbstractContext, True, data );
}
void _lucOutputFormat_Build( void* outputFormat, void* data ) { }
void _lucOutputFormat_Initialise( void* outputFormat, void* data ) { }
void _lucOutputFormat_Execute( void* outputFormat, void* data ) { }
void _lucOutputFormat_Destroy( void* outputFormat, void* data ) { }


void lucOutputFormat_Output( void* outputFormat, lucWindow* window, AbstractContext* context, lucPixel* pixelData ) {
	lucOutputFormat*        self               = (lucOutputFormat*) outputFormat;

	if ( context->rank != MASTER )
		return;

	self->_output( self, window, context, pixelData );
}


Name lucOutputFormat_GetImageFilename( void* outputFormat, lucWindow* window, void* _context ) {
	lucOutputFormat* self       = (lucOutputFormat*) outputFormat;
	AbstractContext* context    = (AbstractContext*) _context;
	Stream*          infoStream = Journal_MyStream( Info_Type, self );
	Name             filename;

	if ( lucWindow_HasStereoCamera( window ) ) 
		Stg_asprintf( &filename, "%s/%s.%05d.%s.%s", 
				context->outputPath, 
				window->name, 
				context->timeStep, 
				window->currStereoBuffer == lucLeft ? "left" : "right", 
				self->extension );
	else 
		Stg_asprintf( &filename, "%s/%s.%05d.%s", context->outputPath, window->name, context->timeStep, self->extension );
	
	Journal_Printf( infoStream, "Creating %s file: %s\n", self->extension, filename );

	return filename;
}

FILE* lucOutputFormat_OpenFile( void* outputFormat, lucWindow* window, void* _context, const char *mode ) {
	lucOutputFormat* self       = (lucOutputFormat*) outputFormat;
	AbstractContext* context    = (AbstractContext*) _context;
	Stream*          error      = Journal_MyStream( Error_Type, self );
	Name             filename;
	FILE*            file;

	filename = lucOutputFormat_GetImageFilename( self, window, context );
	file = fopen( filename, mode );

	Journal_Firewall( file != NULL, error, "In func %s: Cannot open file %s.\n", __func__, filename );

	Memory_Free( filename );
	return file;
}

Stream* lucOutputFormat_OpenStream( void* outputFormat, lucWindow* window, void* context ) {
	lucOutputFormat* self       = (lucOutputFormat*) outputFormat;
	Stream*          stream     = Journal_MyStream( Dump_Type, self );
	Stream*          error      = Journal_MyStream( Error_Type, self );
	Name             filename;
	Bool             result;

	filename = lucOutputFormat_GetImageFilename( self, window, context );
	result = Stream_RedirectFile( stream, filename );

	Journal_Firewall( result, error, "In func %s: Cannot open file %s.\n", __func__, filename );

	Memory_Free( filename );
	
	/* Setup stream */
	Stream_Enable( stream, True );
	Stream_SetAutoFlush( stream, True );

	return stream;
}


