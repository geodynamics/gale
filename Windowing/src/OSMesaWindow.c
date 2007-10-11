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

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>

#include <glucifer/Base/Base.h>

#include <osmesa.h>
#include "types.h"
#include "OSMesaWindow.h"

#include <assert.h>

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type lucOSMesaWindow_Type = "lucOSMesaWindow";

/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
lucOSMesaWindow* _lucOSMesaWindow_New( 
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
		Name                                               name ) 
{
	lucOSMesaWindow*					self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( sizeOfSelf >= sizeof(lucOSMesaWindow) );
	self = (lucOSMesaWindow*) _lucWindow_New( 
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
			name );
	
	return self;
}

void _lucOSMesaWindow_Init( lucOSMesaWindow* self ) {
	self->pixelBuffer = Memory_Alloc_Array( lucAlphaPixel, self->width * self->height, "pixelBuffer" );
	self->osMesaContext = OSMesaCreateContext( OSMESA_RGBA, NULL );
}

void _lucOSMesaWindow_Delete( void* window ) {
	lucOSMesaWindow*  self = (lucOSMesaWindow*)window;

	OSMesaDestroyContext( self->osMesaContext );
	Memory_Free( self->pixelBuffer );

	_lucWindow_Delete( self );
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
	return (void*) _lucOSMesaWindow_New(
		sizeof(lucOSMesaWindow),
		lucOSMesaWindow_Type,
		_lucOSMesaWindow_Delete,
		_lucOSMesaWindow_Print,
		NULL,
		_lucOSMesaWindow_DefaultNew,
		_lucOSMesaWindow_Construct,
		_lucOSMesaWindow_Build,
		_lucOSMesaWindow_Initialise,
		_lucOSMesaWindow_Execute,
		_lucOSMesaWindow_Destroy,
		name );
}

void _lucOSMesaWindow_Construct( void* window, Stg_ComponentFactory* cf, void* data ){
	lucOSMesaWindow*  self = (lucOSMesaWindow*)window;

	/* Construct Parent */
	_lucWindow_Construct( self, cf, data );
	
	_lucOSMesaWindow_Init( self );
}

void _lucOSMesaWindow_Build( void* window, void* data ) {}
void _lucOSMesaWindow_Initialise( void* window, void* data ) {}

void _lucOSMesaWindow_Execute( void* window, void* data ) {
	lucOSMesaWindow*     self      = (lucOSMesaWindow*)window;

	lucDebug_PrintFunctionBegin( self, 1 );
	
	lucWindow_SetViewportNeedsToSetupFlag( self, True );
	lucWindow_SetViewportNeedsToDrawFlag( self, True );

	OSMesaMakeCurrent( self->osMesaContext, self->pixelBuffer, GL_UNSIGNED_BYTE, self->width, self->height );
	_lucWindow_SetupGLRasterFont( self );
	
	lucWindow_Draw( self, data );
	lucWindow_Dump( self, data );
	lucWindow_CleanUp( self, data );
	
	lucDebug_PrintFunctionEnd( self, 1 );
}

void _lucOSMesaWindow_Destroy( void* window, void* data ) {}

#endif /* HAVE_OSMESA */
