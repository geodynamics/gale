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
** $Id: X11Window.c 740 2007-10-11 08:05:31Z SteveQuenette $
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
#include <Power.h>
#endif

#ifndef MASTER
	#define MASTER 0
#endif

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
			name );
	
	return self;
}

void _lucX11Window_Init( 
		lucX11Window*                                      self,
		Name                                               xFontName,
		Name                                               host,
		unsigned int                                       displayNumber,
		unsigned int                                       displayScreen,
		Bool                                               isTimedOut,
		double                                             maxIdleTime ) 
{

	/* Setup display name */
	Stg_asprintf( &self->displayName, "%s:%u.%u", host, displayNumber, displayScreen );
	self->host          = StG_Strdup( host );
	self->displayNumber = displayNumber;
	self->displayScreen = displayScreen;
	self->isTimedOut    = isTimedOut;
	self->maxIdleTime   = maxIdleTime;

	/* Get font from X11 */
	self->xFontName = xFontName;
}

void _lucX11Window_Delete( void* window ) {
	lucX11Window*  self = (lucX11Window*)window;

	lucX11Window_CloseBackgroundWindow( self );
	lucX11Window_CloseInteractiveWindow( self );

	Memory_Free( self->displayName );
	Memory_Free( self->host );

	_lucWindow_Delete( self );
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
			Stg_ComponentFactory_GetUnsignedInt( cf, self->name, "displayScreen", 0 ),
			Stg_ComponentFactory_GetBool( cf, self->name, "isTimedOut", True ),
			Stg_ComponentFactory_GetDouble( cf, self->name, "maxIdleTime", 600.0 ) );
			
			
}

void _lucX11Window_Build( void* window, void* data ) {}
void _lucX11Window_Initialise( void* window, void* data ) {}

void _lucX11Window_Execute( void* window, void* data ) {
	lucX11Window*     self      = (lucX11Window*)window;
	AbstractContext*  context   = (AbstractContext*) data;
	Bool              iAmMaster = (context->rank == MASTER);

	lucDebug_PrintFunctionBegin( self, 1 );
	
	lucWindow_SetViewportNeedsToSetupFlag( self, True );
	lucWindow_SetViewportNeedsToDrawFlag( self, True );

	#ifdef __APPLE__
	UpdateSystemActivity(IdleActivity);
	#endif

	if ( self->interactive ) {
		if (iAmMaster) {
			Journal_DPrintf( lucDebug, "Opening Interactive window.\n");
			lucX11Window_CreateInteractiveWindow( self ); 
		}
		else {
			Journal_DPrintf( lucDebug, "Using background window on slave processor.\n");
			lucX11Window_CreateBackgroundWindow( self );
		}

		lucX11Window_MakeCurrent( self );
		lucX11Window_EventLoop( self, context );
	}
	else {
		/* Open window */
		lucX11Window_CreateBackgroundWindow( self );
		
		lucX11Window_MakeCurrent( self );
		lucWindow_Draw( self, context );
	}
	
	/* Dump image */
	lucWindow_Dump( self, context );

	lucWindow_CleanUp( self, context );
	lucX11Window_CloseInteractiveWindow( self );

	lucDebug_PrintFunctionEnd( self, 1 );
}

void _lucX11Window_Destroy( void* window, void* data ) {}

Bool lucX11Window_EventLoop( void* window, AbstractContext* context) {
	lucX11Window*  self             = (lucX11Window*)window;
	Display*       dpy              = self->display;
	KeySym         ks;
	XEvent         event;
	Atom           wmDeleteWindow;
	unsigned int   button           = 0;
	Pixel_Index    startx           = 0;
	Pixel_Index    starty           = 0;
	Bool           iAmMaster        = (context->rank == MASTER);
	Bool           handledEvent     = False;
	double         maxIdleTime      = MPI_Wtime() + self->maxIdleTime;

	lucDebug_PrintFunctionBegin( self, 1 );
	
	lucWindow_BeginEventLoop( self );

	while( 1 ) {
		if ( iAmMaster ) {
			while ( XEventsQueued(dpy, QueuedAfterFlush ) ) {
				handledEvent = True;
				XNextEvent(dpy, &event);
				switch (event.type) {
					case ConfigureNotify:
						//Reshape( windowingContextExt, event.xconfigure.width, event.xconfigure.height);
						/* fall through... */
					case Expose:
						lucWindow_SetViewportNeedsToDrawFlag( self, True );
						break;
					case ButtonPress: case ButtonRelease:
						button = event.xbutton.button;
						startx = event.xbutton.x;
						starty = self->height - event.xbutton.y;
						lucWindow_MouseClick( self, button, event.type, startx, starty);
						break;
					case MotionNotify:
						lucWindow_MouseMotion( self, button , event.xmotion.x, self->height - event.xmotion.y, startx, starty);
						startx = event.xmotion.x;
						starty = self->height - event.xmotion.y;
						break;
					case KeyPress:
						ks = XLookupKeysym((XKeyEvent *) & event, 0);
						lucWindow_KeyboardEvent( self, ks, event.xkey.x, self->height - event.xkey.y);
						break;
					case ClientMessage:
						if (event.xclient.data.l[0] == wmDeleteWindow)
							lucWindow_QuitEventLoop( self );
						break;
					default:
						lucWindow_SetViewportNeedsToDrawFlag( self, True );
				}
			}

			if(self->isTimedOut){
				if ( MPI_Wtime() > maxIdleTime ) {
					Journal_Printf( lucError, 
							"Error in func '%s'" 
							" - Interactive window '%s' open for too long (over %g seconds) without interaction.\n",
							__func__, self->name, self->maxIdleTime );
					abort();
				}
			}
		}

		MPI_Bcast( &handledEvent, 1, MPI_INT, MASTER, context->communicator );
		if ( !handledEvent )
			continue;

		MPI_Bcast( &self->quitEventLoop, 1, MPI_INT, MASTER, context->communicator );
		if ( self->quitEventLoop )
			break;

		Journal_DPrintfL( lucDebug, 2, "In func %s: Event handled - Now redrawing.\n", __func__ );
		lucWindow_Draw( self, context );
		lucX11Window_SwapBuffers( self, context );
		handledEvent = False;
		maxIdleTime = MPI_Wtime() + self->maxIdleTime;
	}

	if( self->toggleApplicationQuit ) {
		context->gracefulQuit = True;
	}

	lucDebug_PrintFunctionEnd( self, 1 );
	return True;
}

void lucX11Window_SetupFonts( void* window ) {
	lucX11Window*  self      = (lucX11Window*)window;
	
	lucDebug_PrintFunctionBegin( self, 2 );

	/* Tell GLX to use this font */
	glXUseXFont( self->font->fid, 32, 96, 2000+32 );
	glListBase(2000);
	
	lucDebug_PrintFunctionEnd( self, 2 );
}

void lucX11Window_SwapBuffers( void* window, AbstractContext* context ) {
	lucX11Window*  self      = (lucX11Window*)window;
	
	lucDebug_PrintFunctionBegin( self, 2 );

	if (self->doubleBuffer && context->rank == MASTER) 
		glXSwapBuffers(self->display, self->win);
	
	lucDebug_PrintFunctionEnd( self, 2 );
}



void lucX11Window_MakeCurrent( void* window )  {
	lucX11Window*  self         = (lucX11Window*)window;
	
	lucDebug_PrintFunctionBegin( self, 1 );

	if ( self->backgroundWindowOpen ) {
		glXMakeCurrent( self->display, self->glxpmap, self->glxcontext);
	}
	else if ( self->interactiveWindowOpen ){
		glXMakeCurrent( self->display, self->win, self->glxcontext);
		XMapWindow( self->display, self->win );
		glXSwapBuffers( self->display, self->win );
	}

	lucX11Window_SetupFonts( self );
	
	lucDebug_PrintFunctionEnd( self, 1 );
}

Bool lucX11Window_CreateDisplay( void* window )  {
	lucX11Window*  self      = (lucX11Window*)window;
	static int configuration[] = { GLX_DOUBLEBUFFER, GLX_RGBA, GLX_DEPTH_SIZE, 16,
			GLX_RED_SIZE, 1, GLX_GREEN_SIZE, 1, GLX_BLUE_SIZE, 1, None};
	static int alphaConfiguration[] = { GLX_DOUBLEBUFFER, GLX_RGBA, GLX_DEPTH_SIZE, 16,
			GLX_RED_SIZE, 1, GLX_GREEN_SIZE, 1, GLX_BLUE_SIZE, 1, GLX_ALPHA_SIZE, 1, None};
	Display *display;
	XVisualInfo* vi;
	GLXContext   glxcontext;
	
	/*********************** Create Display ******************************/
	display = XOpenDisplay(NULL);
	if (display == NULL) {
		Journal_Printf( lucError, "In func %s: Function XOpenDisplay(NULL) returned NULL\n", __func__);

		/* Second Try */
		display = XOpenDisplay(self->displayName);
		if (display == NULL) {
			Journal_Printf( lucError, "In func %s: Function XOpenDisplay(%s) didn't work.\n", __func__ , self->displayName);
		
			/* Third Try */
			display = XOpenDisplay(":0.0");
			if (display == NULL) {
				Journal_Printf( lucError, "In func %s: Function XOpenDisplay(\":0.0\") returned NULL\n", __func__);

				return False;
			}
		}
	}
	/* Store display on context */
	self->display = display;
	
	/* Check to make sure display we've just opened has a glx extension */
	if (!glXQueryExtension(display, NULL, NULL)) {
		Journal_Printf( lucError,"In func %s: X server has no OpenGL GLX extension\n", __func__);
		return False;
	}

	/* find an OpenGL-capable display - trying different configurations if nessesary */
	vi = glXChooseVisual(display, DefaultScreen(display), &alphaConfiguration[0]);
	self->doubleBuffer = True;
	if (vi == NULL) {
		Journal_Printf( lucError, "In func %s: Couldn't open RGBA Double Buffer display\n", __func__);
		vi = glXChooseVisual( self->display, DefaultScreen( self->display), &alphaConfiguration[1]);
		self->doubleBuffer = False;
	}
	if (vi == NULL) {
		Journal_Printf( lucError, "In func %s: Couldn't open RGBA display\n", __func__);
		vi = glXChooseVisual(display, DefaultScreen(display), &configuration[0]);
		self->doubleBuffer = True;
	}
	if (vi == NULL) {
		Journal_Printf( lucError, "In func %s: Couldn't open Double Buffer display\n", __func__);
		vi = glXChooseVisual(display, DefaultScreen(display), &configuration[1]);
		self->doubleBuffer = False;
	}
	if (vi == NULL) {
		Journal_Printf( lucError, "In func %s: Couldn't open display\n", __func__);
		return False;
	}

	/* Store visual on context */
	self->vi = vi;

	/* Create an OpenGL rendering context */
	glxcontext = glXCreateContext(display, vi, NULL, False);
	if (glxcontext == NULL) {
		Journal_Printf(  lucError, "In func %s: Could not create GLX rendering context.\n", __func__);
		return False;
	}
	/* Store rendering context on our context extension */
	self->glxcontext = glxcontext;
	
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

	lucX11Window_CloseInteractiveWindow( self );

	if ( self->backgroundWindowOpen ) {
		Journal_DPrintf( lucDebug, "Background window already open - returning.\n");
	}
	else {
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
		
		self->backgroundWindowOpen = True;
	}
	
	lucDebug_PrintFunctionEnd( self, 1 );
}


void lucX11Window_CreateInteractiveWindow( void* window ) {
	lucX11Window*        self         = (lucX11Window*)window;
	Pixel_Index          windowWidth  = self->width;
	Pixel_Index          windowHeight = self->height;
	Colormap             cmap;
	XSetWindowAttributes swa;
	XWMHints *           wmHints;
	Atom                 wmDeleteWindow;
	Display*             display;
	XVisualInfo*         vi;
	
	lucDebug_PrintFunctionBegin( self, 1 );
		
	Journal_Firewall( 
		lucX11Window_CreateDisplay( self ),
		lucError,
		"Error in func '%s' for %s '%s': Cannot create display.\n", __func__, self->type, self->name );

	display = self->display;
	vi = self->vi;

	/* Create Colourmap */
	cmap = lucX11Window_GetShareableColormap( self );
	swa.colormap = cmap;
	swa.border_pixel = 0;
	swa.event_mask = ExposureMask | StructureNotifyMask | ButtonPressMask | ButtonMotionMask | KeyPressMask;

	/* Create X window */
	self->win = XCreateWindow(
			display, 
			RootWindow(display, vi->screen),
			0,
			0, 
			windowWidth,
			windowHeight, 
			0,
			vi->depth, 
			InputOutput, 
			vi->visual, 
			CWBorderPixel | CWColormap | CWEventMask, 
			&swa);
		
	XSetStandardProperties(display, self->win, self->name, self->name, None, NULL, 0, NULL);
	
	wmHints = XAllocWMHints();
	wmHints->initial_state = NormalState;
	wmHints->flags = StateHint;
	
	XSetWMHints(display, self->win, wmHints);
	
	wmDeleteWindow = XInternAtom(display, "WM_DELETE_WINDOW", False);
	XSetWMProtocols(display, self->win, &wmDeleteWindow, 1);

	self->interactiveWindowOpen = True;

	lucWindow_InteractionHelpMessage( self, Journal_MyStream( Info_Type, self ) );
	
	lucDebug_PrintFunctionEnd( self, 1 );
}	

Colormap lucX11Window_GetShareableColormap( lucX11Window* self ) {
	Display*           display = self->display;
	XVisualInfo*       vi      =  self->vi;
	Status             status;
	XStandardColormap* standardCmaps;
	Colormap           cmap;
	int                i;
	int                numCmaps;
	
	lucDebug_PrintFunctionBegin( self, 2 );
	
	/* be lazy; using DirectColor too involved for this example */
	if (vi->class != TrueColor)
		Journal_Printf( lucError, "No support for non-TrueColor visual.");
	
	/* if no standard colormap but TrueColor, just make an unshared one */
	status = XmuLookupStandardColormap(display, vi->screen, vi->visualid,
		vi->depth, XA_RGB_DEFAULT_MAP, /* replace */ False, /* retain */ True);
		
	if (status == 1) {
		status = XGetRGBColormaps(display, RootWindow(display, vi->screen),
			&standardCmaps, &numCmaps, XA_RGB_DEFAULT_MAP);
	
		if (status == 1)
			for (i = 0; i < numCmaps; i++)
				if (standardCmaps[i].visualid == vi->visualid) {
					cmap = standardCmaps[i].colormap;
					XFree(standardCmaps);
					lucDebug_PrintFunctionEnd( self, 2 );
					return cmap;
				}
	}
	cmap = XCreateColormap(display, RootWindow(display, vi->screen), vi->visual, AllocNone);
	
	lucDebug_PrintFunctionEnd( self, 2 );
	return cmap;
}

void lucX11Window_CloseInteractiveWindow( lucX11Window* self ) {
	
	lucDebug_PrintFunctionBegin( self, 1 );

	if ( ! self->interactiveWindowOpen ) {
		Journal_DPrintf( lucDebug, "No interactive window open - returning.\n");
	}
	else {
		Journal_DPrintf( lucDebug, "Calling XDestroyWindow.\n");

		XDestroyWindow( self->display , self->win );

		self->win = 0;

		lucX11Window_CloseDisplay( self );
		self->interactiveWindowOpen = False;
	}
	
	lucDebug_PrintFunctionEnd( self, 1 );
}
	
void lucX11Window_CloseBackgroundWindow( lucX11Window* self ) {

	lucDebug_PrintFunctionBegin( self, 1 );

	if ( ! self->backgroundWindowOpen ) {
		Journal_DPrintf( lucDebug, "No background window open - returning.\n");
	}
	else {
		glXDestroyGLXPixmap(self->display, self->glxpmap);
		self->glxpmap = 0;
		XFreePixmap(self->display, self->pmap);
		self->pmap = 0;
		lucX11Window_CloseDisplay( self );
		self->backgroundWindowOpen = False;
	}
	
	lucDebug_PrintFunctionEnd( self, 1 );
}

void lucX11Window_CloseDisplay( lucX11Window* self ) {
	
	lucDebug_PrintFunctionBegin( self, 1 );
			
	XFree( self->vi );
	self->vi = 0;

	glXDestroyContext( self->display,  self->glxcontext);
	self->glxcontext = 0;

	XSetCloseDownMode( self->display,  DestroyAll);
	XCloseDisplay( self->display );
	self->display = 0;
	
	lucDebug_PrintFunctionEnd( self, 1 );
}

#endif
