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
** $Id: SDLWindow.c 740 2007-10-11 08:05:31Z SteveQuenette $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifdef HAVE_SDL

#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>

#include <glucifer/Base/Base.h>

#include "types.h"
#include "SDLWindow.h"
#include <stdlib.h>
#include <signal.h>
#include <assert.h>
#include <gl.h>
#include <glu.h>

#ifdef HAVE_OSMESA
	#include <osmesa.h>
#endif

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type lucSDLWindow_Type = "lucSDLWindow";

/* Globals required to store maximum window size as SDL only allows a single actual window to be open */
int SDL_widthMax = 0;
int SDL_heightMax = 0;
SDL_Surface *screen = NULL;
int SDL_useCount = 0;

/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
lucSDLWindow* _lucSDLWindow_New( 
		SizeT                                           sizeOfSelf,
		Type                                            type,
		Stg_Class_DeleteFunction*                       _delete,
		Stg_Class_PrintFunction*                        _print,
		Stg_Class_CopyFunction*                         _copy, 
		Stg_Component_DefaultConstructorFunction*       _defaultConstructor,
		Stg_Component_ConstructFunction*                _construct,
		Stg_Component_BuildFunction*                    _build,
		Stg_Component_InitialiseFunction*               _initialise,
		Stg_Component_ExecuteFunction*                  _execute,
		Stg_Component_DestroyFunction*                  _destroy,
		lucWindow_DisplayFunction						_displayWindow,	
		lucWindow_EventsWaitingFunction*				_eventsWaiting,	
		lucWindow_EventProcessorFunction				_eventProcessor,	
		lucWindow_ResizeFunction						_resizeWindow,	
		Name                                            name ) 
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
			_displayWindow,	
			_eventsWaiting,
			_eventProcessor,
			_resizeWindow,	
			name );
	
	return self;
}


void _lucSDLWindow_Delete( void* window ) {
	_lucWindow_Delete( window );
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
		_lucSDLWindow_Display,	
		_lucSDLWindow_EventsWaiting,
		_lucSDLWindow_EventProcessor, 
		_lucSDLWindow_Resize,	
		name );
}

void _lucSDLWindow_Construct( void* window, Stg_ComponentFactory* cf, void* data ) {
	lucSDLWindow*  self = (lucSDLWindow*)window;

	/* Construct Parent */
	_lucWindow_Construct( self, cf, data );
} 

void _lucSDLWindow_Build( void* window, void* data ) {
	lucSDLWindow*     self      = (lucSDLWindow*)window;

	/* Run the parent function to build window... */
	_lucWindow_Build(window, data);	

    /* Save largest window dimensions required */
    if (self->width > SDL_widthMax) SDL_widthMax = self->width; 
    if (self->height > SDL_heightMax) SDL_heightMax = self->height; 
}

void _lucSDLWindow_Initialise( void* window, void* data ) {
	lucSDLWindow*     self      = (lucSDLWindow*)window;

	/* Initialise SDL Video subsystem */
    if( SDL_Init( SDL_INIT_VIDEO | SDL_INIT_TIMER ) < 0 ) { 
        Journal_Printf( lucError, "In func %s: Unable to initialize SDL: %s\n", __func__, SDL_GetError() );
        abort();
    }

    putenv("SDL_VIDEO_CENTERED=1");
    
    const SDL_VideoInfo *pSDLVideoInfo = SDL_GetVideoInfo();

    if( !pSDLVideoInfo )
    {
        Journal_Printf( lucError, "In func %s: SDL_GetVideoInfo() failed. SDL Error: %s\n", __func__, SDL_GetError() );
        SDL_Quit();
        exit(1);
    }

    /*** SDL will use OSMesa as the OpenGL implementation if it is present, by copying the OSMesa output
     *** to the SDL display. This allows SDL on-screen and OSMesa off-screen rendering available in the same binary.
     *** For this to work, OSMesa must be linked without/before any other OpenGL implementations. */
  #ifdef HAVE_OSMESA
    Journal_Printf( Journal_MyStream( Info_Type, self ), "*** Using OSMesa library for OpenGL graphics in SDL Window ***.\n" );
    self->osMesaContext = OSMesaCreateContextExt( GL_RGBA, 16, 1, 0, NULL );    /* 16 bit depth, 1 bit stencil */
    if (!self->osMesaContext) {
        Journal_Printf( lucError, "In func %s: OSMesaCreateContext failed!\n", __func__);
        abort();
    }
    if( pSDLVideoInfo->hw_available ) 	/* Hardware surfaces enabled? */
        self->sdlFlags = SDL_RESIZABLE | SDL_HWSURFACE | SDL_DOUBLEBUF;
    else
        self->sdlFlags = SDL_RESIZABLE | SDL_SWSURFACE;
  #else
    self->sdlFlags = SDL_OPENGL | SDL_RESIZABLE;
    /* set opengl attributes */
    SDL_GL_SetAttribute(SDL_GL_RED_SIZE,        8);
    SDL_GL_SetAttribute(SDL_GL_GREEN_SIZE,      8);
    SDL_GL_SetAttribute(SDL_GL_BLUE_SIZE,       8);
    SDL_GL_SetAttribute(SDL_GL_ALPHA_SIZE,      8);
    SDL_GL_SetAttribute(SDL_GL_DEPTH_SIZE,      16);
    SDL_GL_SetAttribute(SDL_GL_STENCIL_SIZE,    1);
    SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER,    1);
  #endif

    /* Create display */	
    lucSDLWindow_CreateWindow(self);
    
    if (self->interactive && self->isMaster)
        /* Install 1sec idle timer */
        self->timer = SDL_AddTimer(1000, lucSDLWindow_IdleTimer, self);

    /* NOTE: we still want Ctrl-C to work, so we undo the SDL redirections */
    signal(SIGINT, SIG_DFL);
    signal(SIGQUIT, SIG_DFL);

	/* Run the parent function to init window... */
	_lucWindow_Initialise(window, data);	

	/* Refresh display */
	//_lucSDLWindow_Display(window);
}

void _lucSDLWindow_Execute( void* window, void* data ) {
	lucSDLWindow*   self = (lucSDLWindow*)window;

    /* Update title */
    SDL_WM_SetCaption( self->title, NULL );
  #ifdef HAVE_OSMESA
    OSMesaMakeCurrent( self->osMesaContext, self->buffer->pixels, GL_UNSIGNED_BYTE, self->width, self->height );
    Journal_DPrintfL( lucDebug, 2, "OSMesa make current %d,%d\n", self->width, self->height);
  #else
    /* Clear background */
		self->renderingEngine->_clear(self, window, False);
  #endif
	/* Run the parent function to execute the window... */
	_lucWindow_Execute(window, data);	
}

void _lucSDLWindow_Destroy( void* window, void* data ) {
	lucSDLWindow*   self = (lucSDLWindow*)window;

	/* Run the parent function to destroy window... */
	_lucWindow_Destroy(window, data);

    lucSDLWindow_DeleteWindow(window);

  #ifdef HAVE_OSMESA
    /* destroy the context */
    if (self->osMesaContext) OSMesaDestroyContext( self->osMesaContext );
    if (!self->interactive || !self->isMaster) return;  /* Already quit sdl */
  #endif

	/* Shut down SDL */
	SDL_Quit();
}

/* Window Virtuals */
void _lucSDLWindow_Display( void* window ) {
	lucSDLWindow*        self = (lucSDLWindow*) window; 

	/* Run the parent function to display window... */
	lucWindow_Display(window);	

  #ifdef HAVE_OSMESA
	/* Render in SDL using OSMesa output buffer */
    if (self->interactive)
    {
        SDL_PixelFormat *fmt = self->screen->format;
        lucColour *c = &self->backgroundColour;
        SDL_FillRect(self->screen, NULL, SDL_MapRGBA(fmt, (Uint8)(c->red * 255), (Uint8)(c->green * 255), 
                                                          (Uint8)(c->blue * 255), (Uint8)(c->opacity * 255)));
        Journal_DPrintfL( lucDebug, 2, "SDL BLIT SURFACE src %d,%d to dst %d,%d\n\n", self->buffer->w, self->buffer->h, self->screen->w, self->screen->h);
    	SDL_BlitSurface(self->buffer,NULL,self->screen,NULL);
	    SDL_Flip(self->screen);
    }
  #else	
	/* Swap buffers */
	SDL_GL_SwapBuffers();
  #endif
}

int _lucSDLWindow_EventsWaiting( void* window )
{
	/* Check for events without removing from queue */
	int numevents = 0;
	SDL_Event events[10];
	SDL_PumpEvents();
	return SDL_PeepEvents(events, 10, SDL_PEEKEVENT, SDL_ALLEVENTS);
}

Bool _lucSDLWindow_EventProcessor( void* window ) {
	lucSDLWindow*   self = (lucSDLWindow*)window;
	char            keyPressed;
	static int      button = 0;
	static Bool     buttonDown = False;
	Bool			redisplay = True;
	SDL_Event       event;
	
	/* Wait for next event */
	SDL_WaitEvent( &event );

	switch( event.type ) {
		case SDL_QUIT:
			lucWindow_ToggleApplicationQuit( window );
			break;	
		case SDL_VIDEORESIZE:
			redisplay = lucWindow_SetSize( self, event.resize.w, event.resize.h);
			break;
		case SDL_KEYDOWN:
			keyPressed = event.key.keysym.sym;
			int xpos, ypos;
			SDL_GetMouseState(&xpos, &ypos);
			lucWindow_KeyboardEvent( self, keyPressed, xpos, self->height - ypos);
			break;
		case SDL_MOUSEMOTION:
			if (buttonDown)
				lucWindow_MouseMotion(self, button, event.motion.x, self->height - event.motion.y);
			else
				redisplay = False;
			break;
		case SDL_MOUSEBUTTONDOWN: 
			buttonDown = True;
			button = event.button.button;
			lucWindow_MouseClick( self, button, event.type, event.button.x, self->height - event.button.y);
			break;
		case SDL_MOUSEBUTTONUP:
			buttonDown = False;
			break;
		case SDL_USEREVENT:
			/* Timer event */
			if (!self->interactive) 
			{
				/* Interactive mode switched off */
                Journal_DPrintfL( lucDebug, 2, "Interactive mode OFF\n");
				lucWindow_SetViewportNeedsToDrawFlag( self, True );
	    		SDL_RemoveTimer(self->timer);
                #ifdef HAVE_OSMESA
                    lucSDLWindow_DeleteWindow(window);
                #endif
				self->quitEventLoop = True;
                self->resized = True;
			}
			else
			{
				/* idle timeout check */
				lucWindow_IdleCheck(self);
		        redisplay = False;
			}
			break;
		case SDL_ACTIVEEVENT:
			if (event.active.state == SDL_APPACTIVE && event.active.gain == 1)	/* Restored from icon */
				lucWindow_SetViewportNeedsToDrawFlag( self, True );
			else
				redisplay = False;
			break;
		case SDL_VIDEOEXPOSE:
		default:	
			redisplay = False;	/* No change to display, don't redraw */
	}
	
	/* Reset idle timer */
	lucWindow_IdleReset(self);

	/* Returns true if display needs refresh */
	return redisplay;
}

void _lucSDLWindow_Resize( void* window ) {
	lucSDLWindow*        self = (lucSDLWindow*) window; 

    /* Free existing window data structures */
    lucSDLWindow_DeleteWindow(self);

    /* Recreate in new dimensions */
    lucSDLWindow_CreateWindow(self);
  #ifdef HAVE_OSMESA
    if (!self->interactive) return; /* Interactive switched off, no need to redo viewports */
  #endif
	/* Run the parent function to resize window viewports... */
	lucWindow_Resize(window);	
}

/* Timer callback */
Uint32 lucSDLWindow_IdleTimer(Uint32 interval, void* param) {
	lucSDLWindow*        self = (lucSDLWindow*) param; 

    /* Create a user event and post */
    SDL_Event event;
    
    event.type = SDL_USEREVENT;
    event.user.code = 1;
    event.user.data1 = 0;
    event.user.data2 = 0;
    
    SDL_PushEvent(&event);

	return interval;
}

void lucSDLWindow_CreateWindow( void* window ) {

	lucSDLWindow* self = (lucSDLWindow*)window;
    Journal_DPrintfL( lucDebug, 2, "*** Create window %d,%d (%s)\n", self->width, self->height, self->name);
    if (self->width > SDL_widthMax) SDL_widthMax = self->width; 
    if (self->height > SDL_heightMax) SDL_heightMax = self->height; 

  #ifdef HAVE_OSMESA
    /* SDL interprets each pixel as a 32-bit number, so our masks depend on the byte order */
	Uint32 rmask, gmask, bmask, amask;
	#if SDL_BYTEORDER == SDL_BIG_ENDIAN 
		rmask = 0xff000000; gmask = 0x00ff0000;	bmask = 0x0000ff00;	amask = 0x000000ff;
	#else
		rmask = 0x000000ff;	gmask = 0x0000ff00;	bmask = 0x00ff0000;	amask = 0xff000000;
	#endif

	/* buffer for display */
	self->buffer = SDL_CreateRGBSurface(SDL_SWSURFACE, self->width, self->height, 32, rmask, gmask,  bmask, amask);
    SDL_SetAlpha(self->buffer, 0, 0);
         
	/* Bind the buffer to the context and make it current */
   	if (!OSMesaMakeCurrent( self->osMesaContext, self->buffer->pixels, GL_UNSIGNED_BYTE, self->width, self->height )) {
		Journal_Printf( lucError, "In func %s: OSMesaMakeCurrent failed!\n", __func__);
		abort();
   	}
	OSMesaPixelStore(OSMESA_Y_UP,0);
    if (!self->interactive || !self->isMaster) return;  /* No SDL window required */
  #endif

    /* Create our rendering surface */
    if (screen == NULL || screen->w < SDL_widthMax || screen->h < SDL_heightMax) 
    {
        if (screen != NULL) SDL_FreeSurface(self->screen);
        
        screen = SDL_SetVideoMode( SDL_widthMax, SDL_heightMax, 32, self->sdlFlags );
        Journal_DPrintfL( lucDebug, 2, "SDL SET VIDEO %d,%d\n\n", SDL_widthMax, SDL_heightMax);
        if (!screen)
        {
            Journal_Printf( lucError, "In func %s: Call to SDL_SetVideoMode() failed! - SDL_Error: %s\n", __func__, SDL_GetError() );
            SDL_Quit();
            abort();
        }
    }
    self->screen = screen;
    SDL_useCount++;;
}

void lucSDLWindow_DeleteWindow(void *window) {
	lucSDLWindow* self = (lucSDLWindow*)window;

  #ifdef HAVE_OSMESA
    /* free the image buffer */
    if (self->buffer) SDL_FreeSurface(self->buffer);
  #endif

    /* Decrement usage count, free if no longer needed */
    if (self->screen) 
    {
        SDL_useCount--;
        if (SDL_useCount < 1)
        { 
            SDL_FreeSurface(self->screen);
            SDL_Quit();
            screen = NULL;
        }
    }

    self->buffer = NULL;
    self->screen = NULL;
	Journal_DPrintfL( lucDebug, 2, "Delete Window %s in %s\n", self->name, __func__);
}

#endif
