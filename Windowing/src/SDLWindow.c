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

#ifdef HAVE_OSMESA
	#include <osmesa.h>
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
		lucWindow_DisplayFunction						   _displayWindow,	
		lucWindow_EventsWaitingFunction*				   _eventsWaiting,	
		lucWindow_EventProcessorFunction				   _eventProcessor,	
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
			_displayWindow,	
			_eventsWaiting,
			_eventProcessor,
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
		name );
}

void _lucSDLWindow_Construct( void* window, Stg_ComponentFactory* cf, void* data ) {
	lucSDLWindow*  self = (lucSDLWindow*)window;

	/* Construct Parent */
	_lucWindow_Construct( self, cf, data );
} 

void _lucSDLWindow_Build( void* window, void* data ) {
	/* Run the parent function to build window... */
	_lucWindow_Build(window, data);	
}

void _lucSDLWindow_Initialise( void* window, void* data ) {
	
	lucSDLWindow*     self      = (lucSDLWindow*)window;

	/* Initialise SDL Video subsystem */
	putenv("SDL_VIDEO_CENTERED=1");
	if( SDL_Init( SDL_INIT_VIDEO | SDL_INIT_TIMER ) < 0 ) { 
		Journal_Printf( lucError, "In func %s: Unable to initialize SDL: %s\n", __func__, SDL_GetError() );
		abort();
	}

    SDL_WM_SetCaption( self->title, NULL );
	/* Hide if not using interactive mode */
	if (!self->interactive) SDL_WM_IconifyWindow();
	
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
	   	self->sdlFlags = SDL_RESIZABLE | SDL_HWPALETTE;
		self->osMesaContext = OSMesaCreateContextExt( GL_RGBA, 16, 0, 0, NULL );
		if (!self->osMesaContext) {
			Journal_Printf( lucError, "In func %s: OSMesaCreateContext failed!\n", __func__);
			abort();
		}
	#else
	   	self->sdlFlags = SDL_OPENGL | SDL_GL_DOUBLEBUFFER | SDL_RESIZABLE | SDL_HWPALETTE;
	#endif

	if( pSDLVideoInfo->hw_available ) // Hardware surfaces enabled?
		self->sdlFlags |= SDL_HWSURFACE;
	else
		self->sdlFlags |= SDL_SWSURFACE;
	if( pSDLVideoInfo->blit_hw ) // Hardware supported blitting?
		self->sdlFlags |= SDL_HWACCEL;

	/* Resize/init the display */	
	self->buffer = NULL;
	lucSDLWindow_Resize(self, self->width, self->height);
	
	/* Install 1sec idle timer */
	self->timer = SDL_AddTimer(1000, lucSDLWindow_IdleTimer, self);
	
	/* NOTE: we still want Ctrl-C to work, so we undo the SDL redirections */
    signal(SIGINT, SIG_DFL);
    signal(SIGQUIT, SIG_DFL);
	
	/* Run the parent function to init window... */
	_lucWindow_Initialise(window, data);	
}

void _lucSDLWindow_Execute( void* window, void* data ) {

	/* Run the parent function to execute the window... */
	_lucWindow_Execute(window, data);	
}

void _lucSDLWindow_Destroy( void* window, void* data ) {
	lucSDLWindow*        self = (lucSDLWindow*) window; 

	/* Run the parent function to destroy window... */
	_lucWindow_Destroy(window, data);	

	#ifdef HAVE_OSMESA
	   /* free the image buffer */
	 	SDL_FreeSurface(self->buffer);

	   /* destroy the context */
	   OSMesaDestroyContext( self->osMesaContext );
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
		/* Render to SDL using OSMesa output buffer */
		SDL_Rect oclip;

		/* make backup of clipping area */
		SDL_GetClipRect(self->screen,&oclip);

		/* clip to full screen */
		SDL_SetClipRect(self->screen,NULL);
		SDL_BlitSurface(self->buffer,NULL,self->screen,NULL);
		SDL_SetClipRect(self->screen,&oclip);
		SDL_Flip(self->screen);
	#else	
		/* Swap buffers */
		SDL_GL_SwapBuffers();
	#endif
}

int _lucSDLWindow_EventsWaiting( void* window )
{
	/* Check for events without removing from queue */
	return SDL_PollEvent(NULL);
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
			lucSDLWindow_Resize(window, event.resize.w, event.resize.h);
			_lucWindow_Initialise(window, self->context);	/* Reset font stuff */
			lucWindow_Resize( self, event.resize.w, event.resize.h);
			break;
		case SDL_KEYDOWN:
			keyPressed = event.key.keysym.sym;
			int xpos, ypos;
			SDL_GetMouseState(&xpos, &ypos);
			lucWindow_KeyboardEvent( self, keyPressed, xpos, ypos);
			break;
		case SDL_MOUSEMOTION:
			if (buttonDown){
				int x = event.motion.x;
				int y = self->height - event.motion.y;
				int dx, dy;
				dx = x - self->startx;
				dy = y - self->starty;
				if (dx * dx + dy * dy > 25)	/* Process if movement magnitude > 5 */
				{
					lucWindow_MouseMotion(self, button , x, y);
					break;
				}
			}
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
		case SDL_ACTIVEEVENT:
			if (event.active.state == SDL_APPACTIVE)	/* Restored from icon */
				lucWindow_SetViewportNeedsToDrawFlag( self, True );
			else
				redisplay = False;
			break;
		case SDL_VIDEOEXPOSE:
		default:	
			redisplay = False;	/* No change to display, don't redraw */
	}

	if (!self->interactive) {
		/* No longer interactive? Drop timer, minimize window and quit event loop */
		SDL_RemoveTimer(self->timer);
		SDL_WM_IconifyWindow();
		self->quitEventLoop = True;
	 }
	 
	/* Reset idle timer */
	lucWindow_IdleReset(self);

	/* Returns true if display needs refresh */
	return redisplay;
}

void lucSDLWindow_Resize( void* window, Pixel_Index width, Pixel_Index height ) {
	lucSDLWindow*     self      = (lucSDLWindow*)window;

    /* Create our rendering surface */
    self->screen = SDL_SetVideoMode( width, height, 32, self->sdlFlags );
    if( !self->screen)
    {
		Journal_Printf( lucError, "In func %s: Call to SDL_SetVideoMode() failed! - SDL_Error: %s\n", __func__, SDL_GetError() );
        SDL_Quit();
		abort();
	}

	#ifdef HAVE_OSMESA
	    /* SDL interprets each pixel as a 32-bit number, so our masks must depend
		   on the endianness (byte order) of the machine */
		Uint32 rmask, gmask, bmask, amask;
		#if SDL_BYTEORDER == SDL_BIG_ENDIAN /* Is this the default anyway? */
			rmask = 0xff000000;
			gmask = 0x00ff0000;
			bmask = 0x0000ff00;
			amask = 0x000000ff;
		#else
			rmask = 0x000000ff;
			gmask = 0x0000ff00;
			bmask = 0x00ff0000;
			amask = 0xff000000;
		#endif

		/* buffer for display */
		if (self->buffer != NULL) SDL_FreeSurface(self->buffer); /* free the image buffer */
		self->buffer = SDL_CreateRGBSurface(SDL_SWSURFACE, width, height, 32, rmask, gmask,  bmask, amask);
         
		/* Bind the buffer to the context and make it current */
	   	if (!OSMesaMakeCurrent( self->osMesaContext, self->buffer->pixels, GL_UNSIGNED_BYTE, width, height )) {
			Journal_Printf( lucError, "In func %s: OSMesaMakeCurrent failed!\n", __func__);
			abort();
	   	}
		OSMesaPixelStore(OSMESA_Y_UP,0);
	#endif
	
}

/* Timer callback */
Uint32 lucSDLWindow_IdleTimer(Uint32 interval, void* param) {
	lucSDLWindow*        self = (lucSDLWindow*) param; 
	/* idle timeout check */
	lucWindow_IdleCheck(self);
	
	/*
    SDL_Event event;
    SDL_UserEvent userevent;
  
    userevent.type = SDL_USEREVENT;
    userevent.code = 0;
    userevent.data1 = NULL;
    userevent.data2 = NULL;

    event.type = SDL_USEREVENT;
    event.user = userevent;
    SDL_PushEvent(&event);
	
    return(interval);*/
}

#endif
