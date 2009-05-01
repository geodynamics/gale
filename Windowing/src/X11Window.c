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
** $Id: X11Window.c 787 2008-08-26 07:57:24Z JulianGiordani $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
 
#ifdef HAVE_X11

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>

#include <glucifer/Base/Base.h>

#include "types.h"
#include "X11Window.h"

#include <assert.h>
#ifdef __APPLE__
#include <CoreServices/CoreServices.h>
#endif

#ifndef MASTER
	#define MASTER 0
#endif

#include <signal.h>
#include <sys/time.h>
lucWindow* parent;	/* Need to save in global so signal handler for idle timer can access */

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type lucX11Window_Type = "lucX11Window";

/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
lucX11Window* _lucX11Window_New( 
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
	lucX11Window*					self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( sizeOfSelf >= sizeof(lucX11Window) );
	self = (lucX11Window*) _lucWindow_New( 
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

void _lucX11Window_Init( 
		lucX11Window*                                      self,
		Name                                               xFontName,
		Name                                               host,
		unsigned int                                       displayNumber,
		unsigned int                                       displayScreen )
{

	/* Setup display name */
	Stg_asprintf( &self->displayName, "%s:%u.%u", host, displayNumber, displayScreen );
	self->host          = StG_Strdup( host );
	self->displayNumber = displayNumber;
	self->displayScreen = displayScreen;

	/* Get font from X11 */
	self->xFontName = xFontName;
}

void _lucX11Window_Delete( void* window ) {
	_lucWindow_Delete( window );
}

void _lucX11Window_Print( void* window, Stream* stream ) {
	lucX11Window*  self = (lucX11Window*)window;

	_lucWindow_Print( self, stream );
}

void* _lucX11Window_Copy( void* window, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap) {
	lucX11Window*  self = (lucX11Window*)window;
	lucX11Window* newWindow;

	newWindow = _lucWindow_Copy( self, dest, deep, nameExt, ptrMap );

	/* TODO */
	abort();

	return (void*) newWindow;
}


void* _lucX11Window_DefaultNew( Name name ) {
	return (void*) _lucX11Window_New(
		sizeof(lucX11Window),
		lucX11Window_Type,
		_lucX11Window_Delete,
		_lucX11Window_Print,
		NULL,
		_lucX11Window_DefaultNew,
		_lucX11Window_Construct,
		_lucX11Window_Build,
		_lucX11Window_Initialise,
		_lucX11Window_Execute,
		_lucX11Window_Destroy,
		_lucX11Window_Display,	
		_lucX11Window_EventsWaiting,
		_lucX11Window_EventProcessor, 
		name );
}

void _lucX11Window_Construct( void* window, Stg_ComponentFactory* cf, void* data ){
	lucX11Window*  self = (lucX11Window*)window;

	/* Construct Parent */
	_lucWindow_Construct( self, cf, data ); 
				
	_lucX11Window_Init( 
			self,
			Stg_ComponentFactory_GetString( cf, self->name, "xFontName", "-adobe-courier-bold-r-normal--14-140-75-75-m-90-iso8859-1" ),
			Stg_ComponentFactory_GetString( cf, self->name, "host", "localhost" ),
			Stg_ComponentFactory_GetUnsignedInt( cf, self->name, "displayNumber", 0 ),
			Stg_ComponentFactory_GetUnsignedInt( cf, self->name, "displayScreen", 0 ) );
			
			
}

void _lucX11Window_Build( void* window, void* data ) {
	/* Run the parent function to build window... */
	_lucWindow_Build(window, data);	
}

void _lucX11Window_Initialise( void* window, void* data ) {
	lucX11Window*  self = (lucX11Window*)window;

	lucDebug_PrintFunctionBegin( self, 1 );

	XSetErrorHandler(lucX11Window_Error);
	
	if ( self->interactive && self->isMaster) {
		Journal_DPrintf( lucDebug, "Opening Interactive window.\n");
		lucX11Window_CreateInteractiveWindow( self ); 
	}
	else {
		Journal_DPrintf( lucDebug, "Opening background window.\n");
		lucX11Window_CreateBackgroundWindow( self );
	}

	/* Setup timer */
	parent = (lucWindow*)window;	/* Save parent ref for signal handler... */
	signal(SIGALRM, lucX11Window_Timer);
	struct itimerval timerval;
	timerval.it_value.tv_sec	= 1;
	timerval.it_value.tv_usec	= 0;
	timerval.it_interval = timerval.it_value;
	setitimer(ITIMER_REAL, &timerval, NULL);

	/* Run the parent function to init window... */
	_lucWindow_Initialise(window, data);

	/* call here to override default font display lists */
	lucX11Window_SetupFonts( self );

	lucDebug_PrintFunctionEnd( self, 1 );
}

void _lucX11Window_Execute( void* window, void* data ) {
	/* Run the parent function to execute the window... */
	_lucWindow_Execute(window, data);	
}

void _lucX11Window_Destroy( void* window, void* data ) {
	lucX11Window*  self = (lucX11Window*)window;

	/* Run the parent function to destroy window... */
	_lucWindow_Destroy(window, data);	

	if (self->win) 
		lucX11Window_CloseInteractiveWindow( self );
	else
		lucX11Window_CloseBackgroundWindow( self );

	/* Close glx window */
	XFree( self->vi );
	self->vi = 0;

	glXDestroyContext( self->display,  self->glxcontext);
	self->glxcontext = 0;

	XSetCloseDownMode( self->display,  DestroyAll);
	XCloseDisplay( self->display );
	self->display = 0;

	/* Moved from _Delete */
	Memory_Free( self->displayName );
	Memory_Free( self->host );

	lucDebug_PrintFunctionEnd( self, 1 );
}

/* Window Virtuals */
void _lucX11Window_Display( void* window ) {
	lucX11Window*        self = (lucX11Window*) window; 

	/* Run the parent function to display window... */
	lucWindow_Display(window);	

	/* Swap buffers if interactive */
	if ( self->isMaster && self->interactive && self->doubleBuffer) 
		glXSwapBuffers(self->display, self->win);
}

int _lucX11Window_EventsWaiting( void* window )
{
	lucX11Window*  self = (lucX11Window*)window;
	return XPending(self->display);

}

Bool _lucX11Window_EventProcessor( void* window ) {
	lucX11Window*  self = (lucX11Window*)window;
	KeySym         ks;
	XEvent         event;
	Atom           wmDeleteWindow;
	static unsigned int button = 0;
	static int visible = True;
	static int hiding = False;
	int width, height;

	lucDebug_PrintFunctionBegin( self, 1 );
	/* Avoid busy wait loop by using blocking event check 
	   (signal callback wakes every few seconds to check idle timeout) */
	
	XNextEvent(self->display, &event);

	/* Reset idle timer */
	lucWindow_IdleReset(window);
	
	switch (event.type) {
		case ButtonRelease:
		case ButtonPress: 
			button = event.xbutton.button;
    		lucWindow_MouseClick( self, button, event.type, event.xmotion.x, self->height - event.xmotion.y);
	    	break;
		case MotionNotify:
			lucWindow_MouseMotion( self, button , event.xmotion.x, self->height - event.xmotion.y);
    		break;
		case KeyPress:
			ks = XLookupKeysym((XKeyEvent *) & event, 0);
			lucWindow_KeyboardEvent( self, ks, event.xkey.x, self->height - event.xkey.y);
			break;
    	case ClientMessage:
			if (event.xclient.data.l[0] == self->wmDeleteWindow)
				lucWindow_ToggleApplicationQuit( window );
			break;
   		case MapNotify:
			/* Window shown */
			if (!visible) lucWindow_SetViewportNeedsToDrawFlag( window, True );
			visible = True;
			break;
   		case UnmapNotify:
			/* Window hidden, iconized */
			visible = False;
			if (hiding) 
			{
				lucWindow_SetViewportNeedsToDrawFlag( window, True );
				self->_displayWindow( self );
				self->interactive = False;
			}
			break;
   		case ConfigureNotify:
		{
			/* Notification of window actions */
			width = event.xconfigure.width;
			height = event.xconfigure.height;
			/* Only resize viewports if window dimensions changed */
			if (width != self->width || height != self->height)
				lucWindow_Resize( self, width, height);
    		break;
		}
	}

	if (!self->interactive && !hiding)
	{
		/* Hide window, must now wait for UnmapNotify */
		hiding = True; 
		self->interactive = True;   /* Reset until events processed */
		XUnmapWindow(self->display, self->win);	 
	}

	/* Returns true if event processed */
	return True;
}

void lucX11Window_SetupFonts( void* window ) {
	lucX11Window*  self      = (lucX11Window*)window;
	
	lucDebug_PrintFunctionBegin( self, 2 );

	/* Tell GLX to use this font */
	glXUseXFont( self->font->fid, 32, 96, 2000+32 );
	
	lucDebug_PrintFunctionEnd( self, 2 );
}

Bool lucX11Window_CreateDisplay( void* window )  {
	lucX11Window*  self      = (lucX11Window*)window;
	static int configuration[] = { GLX_DOUBLEBUFFER, GLX_RGBA, GLX_DEPTH_SIZE, 16,
			GLX_RED_SIZE, 1, GLX_GREEN_SIZE, 1, GLX_BLUE_SIZE, 1, None};
	static int alphaConfiguration[] = { GLX_DOUBLEBUFFER, GLX_RGBA, GLX_DEPTH_SIZE, 16,
			GLX_RED_SIZE, 1, GLX_GREEN_SIZE, 1, GLX_BLUE_SIZE, 1, GLX_ALPHA_SIZE, 1, None};
	
	/*********************** Create Display ******************************/
	self->display = XOpenDisplay(NULL);
	if (self->display == NULL) {
		Journal_Printf( lucError, "In func %s: Function XOpenDisplay(NULL) returned NULL\n", __func__);

		/* Second Try */
		self->display = XOpenDisplay(self->displayName);
		if (self->display == NULL) {
			Journal_Printf( lucError, "In func %s: Function XOpenDisplay(%s) didn't work.\n", __func__ , self->displayName);
		
			/* Third Try */
			self->display = XOpenDisplay(":0.0");
			if (self->display == NULL) {
				Journal_Printf( lucError, "In func %s: Function XOpenDisplay(\":0.0\") returned NULL\n", __func__);

				return False;
			}
		}
	}
	
	/* Check to make sure display we've just opened has a glx extension */
	if (!glXQueryExtension(self->display, NULL, NULL)) {
		Journal_Printf( lucError,"In func %s: X server has no OpenGL GLX extension\n", __func__);
		return False;
	}

	/* find an OpenGL-capable display - trying different configurations if nessesary */
	self->vi = glXChooseVisual(self->display, DefaultScreen(self->display), &alphaConfiguration[0]);
	self->doubleBuffer = True;
	if (self->vi == NULL) {
		Journal_Printf( lucError, "In func %s: Couldn't open RGBA Double Buffer display\n", __func__);
		self->vi = glXChooseVisual( self->display, DefaultScreen( self->display), &alphaConfiguration[1]);
		self->doubleBuffer = False;
	}
	if (self->vi == NULL) {
		Journal_Printf( lucError, "In func %s: Couldn't open RGBA display\n", __func__);
		self->vi = glXChooseVisual(self->display, DefaultScreen(self->display), &configuration[0]);
		self->doubleBuffer = True;
	}
	if (self->vi == NULL) {
		Journal_Printf( lucError, "In func %s: Couldn't open Double Buffer display\n", __func__);
		self->vi = glXChooseVisual(self->display, DefaultScreen(self->display), &configuration[1]);
		self->doubleBuffer = False;
	}
	if (self->vi == NULL) {
		Journal_Printf( lucError, "In func %s: Couldn't open display\n", __func__);
		return False;
	}

	/* Create an OpenGL rendering context */
	self->glxcontext = glXCreateContext(self->display, self->vi, NULL, True); //Enable Direct Rendering
	if (self->glxcontext == NULL) {
		Journal_Printf(  lucError, "In func %s: Could not create GLX rendering context.\n", __func__);
		return False;
	}
	
	self->font = XLoadQueryFont( self->display, self->xFontName );
		
	if (! lucX11Window_FindFont( window) ){
		Journal_Printf(  lucError, "In func %s: Attempt to request a different font failed.\n", __func__);
		Journal_Printf(  lucError, "In func %s: Can not find any suitable font \n", __func__ );
		return True;
	}

	return True;
}

Bool lucX11Window_FindFont( void* window )  {
	lucX11Window*  self      = (lucX11Window*)window;
  	Index font_I = 0;
	Name fontNames[] = {
			"-misc-fixed-bold-r-normal--14-130-75-75-c-70-iso10646-1",
			"-misc-fixed-bold-r-normal--14-130-75-75-c-70-iso8859-1",
			"-misc-fixed-bold-r-normal--14-130-75-75-c-70-iso8859-14",
			"-misc-fixed-bold-r-normal--14-130-75-75-c-70-iso8859-15",
			"-misc-fixed-bold-r-normal--14-130-75-75-c-70-iso8859-2",
			"-misc-fixed-bold-r-normal--14-130-75-75-c-70-iso8859-5",
			"-misc-fixed-bold-r-normal--14-130-75-75-c-70-iso8859-7",
			"-misc-fixed-bold-r-normal--14-130-75-75-c-70-iso8859-8",
			"-misc-fixed-bold-r-normal--14-130-75-75-c-70-iso8859-9",
			"-misc-fixed-medium-o-normal--13-120-75-75-c-70-iso10646-1",
			"-misc-fixed-medium-o-normal--13-120-75-75-c-70-iso8859-1",
			"-misc-fixed-medium-o-normal--13-120-75-75-c-70-iso8859-14",
			"-misc-fixed-medium-o-normal--13-120-75-75-c-70-iso8859-15",
			"-misc-fixed-medium-o-normal--13-120-75-75-c-70-iso8859-2",
			"-misc-fixed-medium-o-normal--13-120-75-75-c-70-iso8859-5",
			"-misc-fixed-medium-o-normal--13-120-75-75-c-70-iso8859-7",
			"-misc-fixed-medium-o-normal--13-120-75-75-c-70-iso8859-9",
			"-misc-fixed-medium-o-normal--13-120-75-75-c-80-iso10646-1",
			"-misc-fixed-medium-o-normal--13-120-75-75-c-80-iso8859-1",
			"-misc-fixed-medium-o-normal--13-120-75-75-c-80-iso8859-14",
			"-misc-fixed-medium-o-normal--13-120-75-75-c-80-iso8859-15",
			"-misc-fixed-medium-o-normal--13-120-75-75-c-80-iso8859-2",
			"-misc-fixed-medium-o-normal--13-120-75-75-c-80-iso8859-5",
			"-misc-fixed-medium-o-normal--13-120-75-75-c-80-iso8859-7",
			"-misc-fixed-medium-o-normal--13-120-75-75-c-80-iso8859-9"
		};
		
	for ( font_I = 0 ; font_I < 25 ; font_I++ ) {
		/* Get font from X11 */
		self->font = XLoadQueryFont( self->display, fontNames[font_I] );
		
		if( self->font != NULL ){
			Stg_asprintf(&self->xFontName, fontNames[font_I]);
			return True;
		}
	}
	if( self->font == NULL )
		return False;
	return True;
}

void lucX11Window_CreateBackgroundWindow( void* window )  {
	lucX11Window*  self         = (lucX11Window*)window;
	
	lucDebug_PrintFunctionBegin( self, 1 );

	Journal_Firewall( 
		lucX11Window_CreateDisplay( self ),
		lucError,
		"Error in func '%s' for %s '%s': Cannot create display.\n", __func__, self->type, self->name );
	
	/* Create Pixmap Window */
	self->pmap = XCreatePixmap(
			self->display,
			RootWindow(self->display, self->vi->screen ), 
			self->width, 
			self->height,
			self->vi->depth );

	self->glxpmap = glXCreateGLXPixmap( self->display,  self->vi, self->pmap );

	glXMakeCurrent( self->display, self->glxpmap, self->glxcontext);
	
	lucDebug_PrintFunctionEnd( self, 1 );
}


void lucX11Window_CreateInteractiveWindow( void* window ) {
	lucX11Window*        self         = (lucX11Window*)window;
	Colormap             cmap;
	XSetWindowAttributes swa;
	XWMHints *           wmHints;
	XSizeHints *		 sHints;
	
	lucDebug_PrintFunctionBegin( self, 1 );
		
	Journal_Firewall( 
		lucX11Window_CreateDisplay( self ),
		lucError,
		"Error in func '%s' for %s '%s': Cannot create display.\n", __func__, self->type, self->name );

	/* Create Colourmap */
	cmap = lucX11Window_GetShareableColormap( self );
	swa.colormap = cmap;
	swa.border_pixel = 0;
	swa.event_mask = StructureNotifyMask | ButtonPressMask | ButtonMotionMask | KeyPressMask;

	/* Setup window manager hints */
	sHints = XAllocSizeHints();
	wmHints = XAllocWMHints();

    if ( sHints && wmHints) {
        sHints->min_width = 32;
        sHints->min_height = 32;
        sHints->max_height = 4096;
        sHints->max_width = 4096;

		sHints->flags = PMaxSize | PMinSize | USPosition;
        /* Center */
	    sHints->x = (DisplayWidth(self->display, self->vi->screen) - self->width) / 2;
		sHints->y = (DisplayHeight(self->display, self->vi->screen) - self->height) / 2;
    

		/* Create X window */
		self->win = XCreateWindow(
			self->display, 
			RootWindow(self->display, self->vi->screen),
			sHints->x,
			sHints->y, 
			self->width,
			self->height, 
			0,
			self->vi->depth, 
			InputOutput, 
			self->vi->visual, 
			CWBorderPixel | CWColormap | CWEventMask, 
			&swa);

		wmHints->initial_state = NormalState;
		wmHints->flags = StateHint;

		XTextProperty title;
		XStringListToTextProperty(&self->title, 1, &title); 	/* argv, argc, normal_hints, wm_hints, class_hints */
		XSetWMProperties(self->display, self->win, &title, &title, NULL, 0, sHints, wmHints, NULL);

		self->wmDeleteWindow = XInternAtom(self->display, "WM_DELETE_WINDOW", True);
		XSetWMProtocols(self->display, self->win, &self->wmDeleteWindow, 1);
	
		glXMakeCurrent( self->display, self->win, self->glxcontext);
	
		XMapWindow( self->display, self->win ); /* Show the window */
	}
	else
		abort();
	
	if (sHints) XFree(sHints);
	if (wmHints) XFree(wmHints);
	
	lucDebug_PrintFunctionEnd( self, 1 );
}	

Colormap lucX11Window_GetShareableColormap( lucX11Window* self ) {
	Status             status;
	XStandardColormap* standardCmaps;
	Colormap           cmap;
	int                i;
	int                numCmaps;
	
	lucDebug_PrintFunctionBegin( self, 2 );
	
	/* be lazy; using DirectColor too involved for this example */
	if (self->vi->class != TrueColor)
		Journal_Printf( lucError, "No support for non-TrueColor visual.");
	
	/* if no standard colormap but TrueColor, just make an unshared one */
	status = XmuLookupStandardColormap(self->display, self->vi->screen, self->vi->visualid,
		self->vi->depth, XA_RGB_DEFAULT_MAP, /* replace */ False, /* retain */ True);
		
	if (status == 1) {
		status = XGetRGBColormaps(self->display, RootWindow(self->display, self->vi->screen),
			&standardCmaps, &numCmaps, XA_RGB_DEFAULT_MAP);
	
		if (status == 1)
			for (i = 0; i < numCmaps; i++)
				if (standardCmaps[i].visualid == self->vi->visualid) {
					cmap = standardCmaps[i].colormap;
					XFree(standardCmaps);
					lucDebug_PrintFunctionEnd( self, 2 );
					return cmap;
				}
	}
	cmap = XCreateColormap(self->display, RootWindow(self->display, self->vi->screen), self->vi->visual, AllocNone);
	
	lucDebug_PrintFunctionEnd( self, 2 );
	return cmap;
}

void lucX11Window_CloseInteractiveWindow( lucX11Window* self ) {
	
	lucDebug_PrintFunctionBegin( self, 1 );

	Journal_DPrintf( lucDebug, "Calling XDestroyWindow.\n");

	XDestroyWindow( self->display , self->win );
	self->win = 0;

	lucDebug_PrintFunctionEnd( self, 1 );
}
	
void lucX11Window_CloseBackgroundWindow( lucX11Window* self ) {

	lucDebug_PrintFunctionBegin( self, 1 );

	glXDestroyGLXPixmap(self->display, self->glxpmap);
	self->glxpmap = 0;
	XFreePixmap(self->display, self->pmap);
	self->pmap = 0;

	lucDebug_PrintFunctionEnd( self, 1 );
}

void lucX11Window_Timer( int x)
{
	if (parent->interactive)
		lucWindow_IdleCheck(parent);
	else
	{
		/* Remove timer */
		struct itimerval timerval;
		timerval.it_value.tv_sec = 0;
		timerval.it_value.tv_usec = 0;
		timerval.it_interval = timerval.it_value;
		setitimer(ITIMER_REAL, &timerval, NULL);
	}
}

int lucX11Window_Error(Display* display, XErrorEvent* error)
{
	char error_str[256];
	XGetErrorText(display, error->error_code, error_str, 256);
	Journal_DPrintf( lucError, "X11 Error: %d -> %s\n", error->error_code, error_str);
}

#endif
