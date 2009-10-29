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
** $Id: RenderingEngineTest.c 740 2007-10-11 08:05:31Z SteveQuenette $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>

#include <glucifer/Base/Base.h>

#include <assert.h>
#include <string.h>
#include "RenderingEngineTest.h"

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type lucRenderingEngineTest_Type = "lucRenderingEngineTest";

/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
lucRenderingEngineTest* _lucRenderingEngineTest_New( 
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
		lucRenderingEngine_RenderFunction*                 _render,
		lucRenderingEngine_GetPixelDataFunction*           _getPixelData,
		lucRenderingEngine_CompositeViewportFunction*      _compositeViewport,
		Name                                               name ) 
{
	lucRenderingEngineTest*					self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( sizeOfSelf >= sizeof(lucRenderingEngineTest) );
	self = (lucRenderingEngineTest*) _lucRenderingEngine_New( 
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
			_render,
			_getPixelData,
			_compositeViewport,
			name );
	
	return self;
}

void _lucRenderingEngineTest_Init( 
		lucRenderingEngineTest*                                      self  ) 
{
	/* Initial malloc of memory */
	self->buffer = Memory_Alloc_Array( lucPixel, 1, "buffer" );
}

void _lucRenderingEngineTest_Delete( void* renderingEngine ) {
	lucRenderingEngineTest*  self = (lucRenderingEngineTest*)renderingEngine;

	_lucRenderingEngine_Delete( self );
}

void _lucRenderingEngineTest_Print( void* renderingEngine, Stream* stream ) {
	lucRenderingEngineTest*  self = (lucRenderingEngineTest*)renderingEngine;

	_lucRenderingEngine_Print( self, stream );
}

void* _lucRenderingEngineTest_Copy( void* renderingEngine, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap) {
	lucRenderingEngineTest*  self = (lucRenderingEngineTest*)renderingEngine;
	lucRenderingEngineTest* newRenderingEngine;

	newRenderingEngine = _lucRenderingEngine_Copy( self, dest, deep, nameExt, ptrMap );

	/* TODO */
	abort();

	return (void*) newRenderingEngine;
}


void* _lucRenderingEngineTest_DefaultNew( Name name ) {
	return (void*) _lucRenderingEngineTest_New(
		sizeof(lucRenderingEngineTest),
		lucRenderingEngineTest_Type,
		_lucRenderingEngineTest_Delete,
		_lucRenderingEngineTest_Print,
		NULL,
		_lucRenderingEngineTest_DefaultNew,
		_lucRenderingEngineTest_AssignFromXML,
		_lucRenderingEngineTest_Build,
		_lucRenderingEngineTest_Initialise,
		_lucRenderingEngineTest_Execute,
		_lucRenderingEngineTest_Destroy,
		_lucRenderingEngineTest_Render,
		_lucRenderingEngineTest_GetPixelData,
		_lucRenderingEngineTest_CompositeViewport,
		name );
}

void _lucRenderingEngineTest_AssignFromXML( void* renderingEngine, Stg_ComponentFactory* cf, void* data ){
	lucRenderingEngineTest*  self = (lucRenderingEngineTest*)renderingEngine;

	/* Construct Parent */
	_lucRenderingEngine_AssignFromXML( self, cf, data );
	
	_lucRenderingEngineTest_Init( self );
}

void _lucRenderingEngineTest_Build( void* renderingEngine, void* data ) {}
void _lucRenderingEngineTest_Initialise( void* renderingEngine, void* data ) {}
void _lucRenderingEngineTest_Execute( void* renderingEngine, void* data ) {}
void _lucRenderingEngineTest_Destroy( void* renderingEngine, void* data ) {}


void _lucRenderingEngineTest_Render( void* renderingEngine, lucWindow* window, AbstractContext* context ) {
	lucRenderingEngineTest* self              = (lucRenderingEngineTest*) renderingEngine;
	Pixel_Index           width  = window->width;
	Pixel_Index           height = window->height;
	Pixel_Index           horizontal_I;
	Pixel_Index           vertical_I;
	unsigned char*        pixel;

	Journal_DPrintfL( lucDebug, 2, "In func: %s for %s '%s'\n", __func__, self->type, self->name );
	Stream_Indent( lucDebug );

	self->buffer = Memory_Realloc_Array( self->buffer, lucPixel, width * height );
	memset( self->buffer, 0, width * height * sizeof(lucPixel) );
	
	for ( vertical_I = 0 ; vertical_I < height ; vertical_I++ ) {
		for ( horizontal_I = 0 ; horizontal_I < width ; horizontal_I++ ) {
			pixel = (unsigned char*) &self->buffer[ horizontal_I + vertical_I * width ];

			pixel[0] = (unsigned char) (255.0/(double)height * (double)vertical_I);
			pixel[1] = (unsigned char) (255.0 - 255.0/(double)width * (double)horizontal_I);
			pixel[2] = (unsigned char) (255.0/(double)width * (double)horizontal_I);
		}
	}
	
	Stream_UnIndent( lucDebug );
	Journal_DPrintfL( lucDebug, 2, "Leaving func %s\n", __func__ );
}

void _lucRenderingEngineTest_GetPixelData( void* renderingEngine, lucWindow* window, lucPixel* buffer ) {
	lucRenderingEngineTest* self              = (lucRenderingEngineTest*) renderingEngine;
	Pixel_Index             width             = window->width;
	Pixel_Index             height            = window->height;

	memcpy( buffer, self->buffer, sizeof(lucPixel) * width * height );
}

void _lucRenderingEngineTest_CompositeViewport( 
		void*                                              renderingEngine, 
		lucViewportInfo*                                   viewportInfo, 
		AbstractContext*                                   context, 
		Bool                                               broadcast ) 
{
}

void RenderingEngineTest_Register( AbstractContext* context ) {
	Stg_ComponentRegister* componentRegister = Stg_ComponentRegister_Get_ComponentRegister();
	Stg_ComponentRegister_Add( componentRegister, lucRenderingEngineTest_Type,     "0", _lucRenderingEngineTest_DefaultNew );
	RegisterParent( lucRenderingEngineTest_Type, lucRenderingEngine_Type );
}
