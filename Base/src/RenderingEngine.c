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
** $Id: RenderingEngine.c 740 2007-10-11 08:05:31Z SteveQuenette $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>

#include "types.h"
#include "ColourMap.h"
#include "X11Colours.h"
#include "Window.h"
#include "RenderingEngine.h"
#include "Init.h"

#include <assert.h>

const Type lucRenderingEngine_Type = "lucRenderingEngine";

lucRenderingEngine* _lucRenderingEngine_New(
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
		lucRenderingEngine_ClearFunction*				   _clear,
		lucRenderingEngine_GetPixelDataFunction*           _getPixelData,
		lucRenderingEngine_CompositeViewportFunction*      _compositeViewport,
		Name                                               name )
{
	lucRenderingEngine*    self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( sizeOfSelf >= sizeof(lucRenderingEngine) );
	self = (lucRenderingEngine*) _Stg_Component_New( 
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
			name, 
			NON_GLOBAL );

	self->_render            = _render;
	self->_clear			 = _clear;
	self->_getPixelData      = _getPixelData;
	self->_compositeViewport = _compositeViewport;
	
	return self;
}

void _lucRenderingEngine_Init( lucRenderingEngine* self ) {
	self->isConstructed = True;
}

void lucRenderingEngine_InitAll( void* renderingEngine ) {
	lucRenderingEngine* self        = renderingEngine;

	_lucRenderingEngine_Init( self );
}

void _lucRenderingEngine_Delete( void* renderingEngine ) {
	lucRenderingEngine* self        = renderingEngine;
	
	_Stg_Component_Delete( self );
}

void _lucRenderingEngine_Print( void* renderingEngine, Stream* stream ) {
	lucRenderingEngine* self        = renderingEngine;
	
	Journal_Printf( stream, "lucRenderingEngine: %s\n", self->name );

	Stream_Indent( stream );

	/* Print Parent */
	_Stg_Component_Print( self, stream );

	Journal_PrintPointer( stream, self->_getPixelData );
	
	Stream_UnIndent( stream );
}

void* _lucRenderingEngine_Copy( void* renderingEngine, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	lucRenderingEngine* self        = renderingEngine;
	lucRenderingEngine* newRenderingEngine;

	newRenderingEngine = _Stg_Component_Copy( self, dest, deep, nameExt, ptrMap );

	return (void*) newRenderingEngine;
}

void _lucRenderingEngine_AssignFromXML( void* renderingEngine, Stg_ComponentFactory* cf, void* data ) {
	lucRenderingEngine*        self            = (lucRenderingEngine*) renderingEngine ;

	self->context = Stg_ComponentFactory_ConstructByKey( cf, self->name, "Context", AbstractContext, False, data );
	if( !self->context ) 
		self->context = Stg_ComponentFactory_ConstructByName( cf, "context", AbstractContext, True, data );

	_lucRenderingEngine_Init( self );
}

void _lucRenderingEngine_Build( void* camera, void* data ) { }
void _lucRenderingEngine_Initialise( void* camera, void* data ) { }
void _lucRenderingEngine_Execute( void* camera, void* data ) { }
void _lucRenderingEngine_Destroy( void* camera, void* data ) { }

void lucRenderingEngine_Render( void* renderingEngine, lucWindow* window, AbstractContext* context ) {
	lucRenderingEngine*   self       = (lucRenderingEngine*) renderingEngine ;
	
	lucDebug_PrintFunctionBegin( self, 2 );

	self->_render( self, window, context );
	
	lucDebug_PrintFunctionEnd( self, 2 );
}

void lucRenderingEngine_Clear( void* renderingEngine,  lucWindow* window, Bool clearAll ) {
	lucRenderingEngine*   self       = (lucRenderingEngine*) renderingEngine ;
	
	lucDebug_PrintFunctionBegin( self, 2 );

	self->_clear( self, window, clearAll );
	
	lucDebug_PrintFunctionEnd( self, 2 );
}

void lucRenderingEngine_GetPixelData( void* renderingEngine, lucWindow* window, lucPixel* pixelData ) {
	lucRenderingEngine*   self       = (lucRenderingEngine*) renderingEngine ;
	
	lucDebug_PrintFunctionBegin( self, 2 );

	self->_getPixelData( self, window, pixelData );

	lucDebug_PrintFunctionEnd( self, 2 );
}

void lucRenderingEngine_CompositeViewport( 
		void*                                              renderingEngine, 
		lucViewportInfo*                                   viewportInfo, 
		AbstractContext*                                   context, 
		Bool                                               broadcast )
{
	lucRenderingEngine*   self       = (lucRenderingEngine*) renderingEngine ;
	
	lucDebug_PrintFunctionBegin( self, 2 );

	if ( context->nproc ==1 ) {
		Journal_DPrintfL( lucDebug, 2, "Running in serial - No need to composite.\n" );
	}
	else {
		self->_compositeViewport( self, viewportInfo, context, broadcast );
	}

	lucDebug_PrintFunctionEnd( self, 2 );
}
