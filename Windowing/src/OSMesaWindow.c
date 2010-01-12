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
** $Id: OSMesaWindow.c 740 2007-10-11 08:05:31Z SteveQuenette $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifdef HAVE_OSMESA

#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>

#include <glucifer/Base/Base.h>

#include <glu.h>
#include <osmesa.h>
#include "types.h"
#include "OSMesaWindow.h"

#include <assert.h>

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type lucOSMesaWindow_Type = "lucOSMesaWindow";

/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
lucOSMesaWindow* _lucOSMesaWindow_New(  LUCOSMESAWINDOW_DEFARGS  ) 
{
	lucOSMesaWindow*					self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( _sizeOfSelf >= sizeof(lucOSMesaWindow) );
	self = (lucOSMesaWindow*) _lucWindow_New(  LUCWINDOW_PASSARGS  );
	
	return self;
}

void _lucOSMesaWindow_Init( lucOSMesaWindow* self ) {
}

void _lucOSMesaWindow_Delete( void* window ) {
	_lucWindow_Delete( window );
}

void _lucOSMesaWindow_Print( void* window, Stream* stream ) {
	lucOSMesaWindow*  self = (lucOSMesaWindow*)window;

	_lucWindow_Print( self, stream );
	Journal_PrintPointer( stream, self->pixelBuffer );
	Journal_PrintPointer( stream, self->osMesaContext );
}

void* _lucOSMesaWindow_Copy( void* window, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap) {
	lucOSMesaWindow*  self = (lucOSMesaWindow*)window;
	lucOSMesaWindow* newWindow;

	newWindow = _lucWindow_Copy( self, dest, deep, nameExt, ptrMap );

	/* TODO */
	abort();

	return (void*) newWindow;
}


void* _lucOSMesaWindow_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                              _sizeOfSelf = sizeof(lucOSMesaWindow);
	Type                                                      type = lucOSMesaWindow_Type;
	Stg_Class_DeleteFunction*                              _delete = _lucOSMesaWindow_Delete;
	Stg_Class_PrintFunction*                                _print = _lucOSMesaWindow_Print;
	Stg_Class_CopyFunction*                                  _copy = NULL;
	Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _lucOSMesaWindow_DefaultNew;
	Stg_Component_ConstructFunction*                    _construct = _lucOSMesaWindow_AssignFromXML;
	Stg_Component_BuildFunction*                            _build = _lucOSMesaWindow_Build;
	Stg_Component_InitialiseFunction*                  _initialise = _lucOSMesaWindow_Initialise;
	Stg_Component_ExecuteFunction*                        _execute = _lucOSMesaWindow_Execute;
	Stg_Component_DestroyFunction*                        _destroy = _lucOSMesaWindow_Destroy;
	lucWindow_DisplayFunction*                      _displayWindow = lucWindow_Display;
	lucWindow_EventsWaitingFunction*                _eventsWaiting = lucWindow_EventsWaiting;
	lucWindow_EventProcessorFunction*              _eventProcessor = lucWindow_EventProcessor;
	lucWindow_ResizeFunction*                        _resizeWindow = lucWindow_Resize;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

	return (void*) _lucOSMesaWindow_New(  LUCOSMESAWINDOW_PASSARGS  );
}

void _lucOSMesaWindow_AssignFromXML( void* window, Stg_ComponentFactory* cf, void* data ){
	lucOSMesaWindow*  self = (lucOSMesaWindow*)window;

	/* Construct Parent */
	_lucWindow_AssignFromXML( self, cf, data );
	
	_lucOSMesaWindow_Init( self );
}

void _lucOSMesaWindow_Build( void* window, void* data ) {
	/* Run the parent function to build window... */
	_lucWindow_Build(window, data);	
}

void _lucOSMesaWindow_Initialise( void* window, void* data ) {
	lucOSMesaWindow*     self      = (lucOSMesaWindow*)window;

	/* Init OSMesa display buffer */
	self->pixelBuffer = Memory_Alloc_Array( lucAlphaPixel, self->width * self->height, "OSMesa pixelBuffer" );
    /* 24 bit depth, 1 bit stencil, no accum */
	self->osMesaContext = OSMesaCreateContextExt( OSMESA_RGBA, 24, 1, 0, NULL); 

	OSMesaMakeCurrent( self->osMesaContext, self->pixelBuffer, GL_UNSIGNED_BYTE, self->width, self->height );

	/* Run the parent function to init window... */
	_lucWindow_Initialise(window, data);	
}

void _lucOSMesaWindow_Execute( void* window, void* data ) {
    /* Make sure we are using the correct context when more than one are created */
	lucOSMesaWindow*     self      = (lucOSMesaWindow*)window;
	OSMesaMakeCurrent( self->osMesaContext, self->pixelBuffer, GL_UNSIGNED_BYTE, self->width, self->height );

	/* Run the parent function to execute window... */
	_lucWindow_Execute(window, data);	
}

void _lucOSMesaWindow_Destroy( void* window, void* data ) {
	lucOSMesaWindow*     self      = (lucOSMesaWindow*)window;

	/* Moved from _Delete */
	OSMesaDestroyContext( self->osMesaContext );
	Memory_Free( self->pixelBuffer );

	/* Run the parent function to destroy window... */
	_lucWindow_Destroy(window, data);	
}

#endif /* HAVE_OSMESA */


