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
** $Id: SDLWindow.c 628 2006-10-12 08:23:07Z SteveQuenette $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifdef HAVE_SDL

#include <mpi.h>
#include <StGermain/StGermain.h>

#include <glucifer/Base/Base.h>

#include "types.h"
#include "SDLWindow.h"

#include <assert.h>
#include <SDL/SDL.h>

#ifndef MASTER
	#define MASTER 0
#endif

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type lucSDLWindow_Type = "lucSDLWindow";

/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
lucSDLWindow* _lucSDLWindow_New( 
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
	lucSDLWindow*					self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( sizeOfSelf >= sizeof(lucSDLWindow) );
	self = (lucSDLWindow*) _lucWindow_New( 
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

void _lucSDLWindow_Init( lucSDLWindow*                                      self ) {
}

void _lucSDLWindow_Delete( void* window ) {
	lucSDLWindow*  self = (lucSDLWindow*)window;

	_lucWindow_Delete( self );
}

void _lucSDLWindow_Print( void* window, Stream* stream ) {
	lucSDLWindow*  self = (lucSDLWindow*)window;

	_lucWindow_Print( self, stream );
}

void* _lucSDLWindow_Copy( void* window, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap) {
	lucSDLWindow*  self = (lucSDLWindow*)window;
	lucSDLWindow* newWindow;

	newWindow = _lucWindow_Copy( self, dest, deep, nameExt, ptrMap );

	/* TODO */
	abort();

	return (void*) newWindow;
}


void* _lucSDLWindow_DefaultNew( Name name ) {
	return (void*) _lucSDLWindow_New(
		sizeof(lucSDLWindow),
		lucSDLWindow_Type,
		_lucSDLWindow_Delete,
		_lucSDLWindow_Print,
		NULL,
		_lucSDLWindow_DefaultNew,
		_lucSDLWindow_Construct,
		_lucSDLWindow_Build,
		_lucSDLWindow_Initialise,
		_lucSDLWindow_Execute,
		_lucSDLWindow_Destroy,
		name );
}

void _lucSDLWindow_Construct( void* window, Stg_ComponentFactory* cf, void* data ){
	lucSDLWindow*  self = (lucSDLWindow*)window;

	/* Construct Parent */
	_lucWindow_Construct( self, cf, data );
	
	_lucSDLWindow_Init( self );
} 

void _lucSDLWindow_Build( void* window, void* data ) {}
void _lucSDLWindow_Initialise( void* window, void* data ) {}

void _lucSDLWindow_Execute( void* window, void* data ) {
	lucSDLWindow*     self      = (lucSDLWindow*)window;
	AbstractContext*  context   = (AbstractContext*) data;

	lucDebug_PrintFunctionBegin( self, 1 );
	
	lucWindow_SetViewportNeedsToSetupFlag( self, True );
	lucWindow_SetViewportNeedsToDrawFlag( self, True );

	/* Initialise SDL */
	if( SDL_Init( SDL_INIT_VIDEO ) < 0 ) {
		Journal_Printf( lucError, "In func %s: Unable to initialize SDL: %s\n", __func__, SDL_GetError() );
		abort();
	}

	/* Create a OpenGL screen */
	if( SDL_SetVideoMode( self->width, self->height, 0, SDL_OPENGL ) == NULL ) {
		Journal_Printf( lucError, "In func %s: Unable to create OpenGL screen: %s\n", __func__, SDL_GetError() );
		abort();
	}

	_lucWindow_SetupGLRasterFont( self );
	
	/* Draw Window */
	lucWindow_Draw( self, context );
	SDL_GL_SwapBuffers();
	if ( self->interactive ) {
		lucWindow_InteractionHelpMessage( self, Journal_MyStream( Info_Type, self ) );
		lucSDLWindow_EventLoop( self, context );
	}
	
	/* Dump image */
	lucWindow_Dump( self, context );
	lucWindow_CleanUp( self, context );
	
	SDL_Quit();
}

void _lucSDLWindow_Destroy( void* window, void* data ) {}

void lucSDLWindow_EventLoop( void* window, AbstractContext* context) {
	lucSDLWindow*     self      = (lucSDLWindow*)window;
	SDL_Event         event;
	Pixel_Index       startx     = 0;
	Pixel_Index       starty     = 0;
	int               button     = 0;
	char              keyPressed;
	Bool              buttonDown = False;
	
	lucWindow_BeginEventLoop( self );
	while ( !self->quitEventLoop ) {
		if( SDL_WaitEvent( &event ) ) {
			switch( event.type ) {
				case SDL_QUIT:
					return ;
				case SDL_VIDEORESIZE:
					lucWindow_SetViewportNeedsToDrawFlag( self, True );
					break;
				case SDL_MOUSEMOTION:
					if (!buttonDown) break;
					lucWindow_MouseMotion(self, button , event.motion.x, self->height - event.motion.y, startx, starty);
					startx = event.motion.x;
					starty = self->height - event.motion.y;
					break;
				case SDL_KEYDOWN:
					keyPressed = event.key.keysym.sym;
					lucWindow_KeyboardEvent( self, keyPressed, startx, starty);
					break;
				case SDL_MOUSEBUTTONDOWN: 
					button = event.button.button;
					startx = event.button.x;
					starty = self->height - event.button.y;
					lucWindow_MouseClick( self, button, event.type, startx, starty);
					buttonDown = True;
					break;
				case SDL_MOUSEBUTTONUP:
					buttonDown = False;
					break;
				default:
					lucWindow_SetViewportNeedsToDrawFlag( self, True );
			}

			/* check to see whether we should continue interactivity */
			MPI_Bcast( &self->quitEventLoop, 1, MPI_INT, MASTER, context->communicator );
			if ( self->quitEventLoop )
				break;
			
			lucWindow_Draw( window, context );
			SDL_GL_SwapBuffers();
		}
	}

}


#endif

