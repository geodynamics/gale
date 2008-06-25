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
** $Id: Arrhenius.c 78 2005-11-29 11:58:21Z RobertTurnbull $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
#ifdef HAVE_CARBON

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>

#include <Carbon/Carbon.h>
#include <AGL/agl.h>
#include <AGL/agl.h>
#include <OpenGL/OpenGL.h>

#include <glucifer/Base/Base.h>

#include "types.h"
#include "CarbonWindow.h"

#include <assert.h>

#ifndef MASTER
	#define MASTER 0
#endif
	
static pascal OSStatus lucCarbonWindow_EventHandler(EventHandlerCallRef nextHandler, EventRef event, void *data) ;
Bool lucCarbonWindow_IsPointOutsideWindow( lucCarbonWindow* self, Point* point ) ;
void lucCarbonWindow_GetPixelIndicies( lucCarbonWindow* self, Point* point, Pixel_Index *xPos, Pixel_Index *yPos ) ;

/* this is a secret undocumented Mac function that allows code that isn't part of a bundle to be a foreground operation 
 * we need to give the prototype here because the function isn't in the header files. */
void CPSEnableForegroundOperation( void* psn ); 

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type lucCarbonWindow_Type = "lucCarbonWindow";

/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
lucCarbonWindow* _lucCarbonWindow_New( 
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
	lucCarbonWindow*					self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( sizeOfSelf >= sizeof(lucCarbonWindow) );
	self = (lucCarbonWindow*) _lucWindow_New( 
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

void _lucCarbonWindow_Init( 
		lucCarbonWindow*                                   self,
		Pixel_Index                                        offsetX,
		Pixel_Index                                        offsetY,
		double                                             maxIdleTime ) 
{
	int titleLength;

	/* Title */
	titleLength = strlen( self->name );
	self->title = Memory_Alloc_Array( unsigned char, titleLength + 2, "title + extras" );
	self->title[ 0 ] = (unsigned char) titleLength;
	strcpy( (char*) (self->title + 1), self->name );
	
	self->offsetX = offsetX;
	self->offsetY = offsetY;

	self->fontID = 0;
	self->fontSize = 12;
}

void _lucCarbonWindow_Delete( void* window ) {
	lucCarbonWindow*  self = (lucCarbonWindow*)window;

	Memory_Free( self->title );

	_lucWindow_Delete( self );
}

void _lucCarbonWindow_Print( void* window, Stream* stream ) {
	lucCarbonWindow*  self = (lucCarbonWindow*)window;

	_lucWindow_Print( self, stream );
}

void* _lucCarbonWindow_Copy( void* window, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap) {
	lucCarbonWindow*  self = (lucCarbonWindow*)window;
	lucCarbonWindow* newWindow;

	newWindow = _lucWindow_Copy( self, dest, deep, nameExt, ptrMap );

	/* TODO */
	abort();

	return (void*) newWindow;
}


void* _lucCarbonWindow_DefaultNew( Name name ) {
	return (void*) _lucCarbonWindow_New(
		sizeof(lucCarbonWindow),
		lucCarbonWindow_Type,
		_lucCarbonWindow_Delete,
		_lucCarbonWindow_Print,
		NULL,
		_lucCarbonWindow_DefaultNew,
		_lucCarbonWindow_Construct,
		_lucCarbonWindow_Build,
		_lucCarbonWindow_Initialise,
		_lucCarbonWindow_Execute,
		_lucCarbonWindow_Destroy,
		name );
}

void _lucCarbonWindow_Construct( void* window, Stg_ComponentFactory* cf, void* data ){
	lucCarbonWindow*  self = (lucCarbonWindow*)window;

	/* Construct Parent */
	_lucWindow_Construct( self, cf , data);
	
	_lucCarbonWindow_Init( 
			self,
			Stg_ComponentFactory_GetUnsignedInt( cf, self->name, "offsetX", 50 ),
			Stg_ComponentFactory_GetUnsignedInt( cf, self->name, "offsetY", 50 ),
			Stg_ComponentFactory_GetDouble( cf, self->name, "maxIdleTime", 600.0 ) );
}

void _lucCarbonWindow_Build( void* window, void* data ) {}
void _lucCarbonWindow_Initialise( void* window, void* data ) {}

void _lucCarbonWindow_Execute( void* carbonWindow, void* data ) {
	lucCarbonWindow*        self       = (lucCarbonWindow*) carbonWindow; 
	
	lucWindow_SetViewportNeedsToSetupFlag( self, True );
	lucWindow_SetViewportNeedsToDrawFlag( self, True );
	
	self->currentContext = (AbstractContext*) data;
	
	if (self->interactive) {
		_lucCarbonWindow_ExecuteInteractive( self, data );
	}
	else {
		_lucCarbonWindow_ExecuteOffscreen( self, data );
	}
}


void _lucCarbonWindow_ExecuteOffscreen( void* carbonWindow, void* data ) {
	lucCarbonWindow*        self       = (lucCarbonWindow*) carbonWindow; 

	lucCarbonWindow_CreateBackgroundWindow( self, data );

	/***** Perform offscreen drawing *****/
	/* Setup fonts */
	_lucWindow_SetupGLRasterFont( self );

	lucWindow_Draw( self, self->currentContext );	
	
	/* Dump image */
	lucWindow_Dump( self, data );

	/* Clean Up */
	lucWindow_CleanUp( self, data );
	
	lucCarbonWindow_DestroyWindow( self, data );
}


	
void _lucCarbonWindow_ExecuteInteractive( void* carbonWindow, void* data ) {
	lucCarbonWindow*        self = (lucCarbonWindow*) carbonWindow; 
	AbstractContext*        context;
	Bool                    iAmMaster;

	lucDebug_PrintFunctionBegin( self, 1 );

	context   = self->currentContext;
	iAmMaster = ( context->rank == MASTER );

	/* Create the window...  */
	lucCarbonWindow_CreateWindow( self, data );

	/* Loop for interactivity */
	lucWindow_BeginEventLoop( self );
	while ( !self->quitEventLoop ) {
		if ( iAmMaster ) {
			if (self->windowIsVisible)
				SetEventLoopTimerNextFireTime( self->timer, 0.05);

			RunApplicationEventLoop();
		}

		/* Broadcast information about loop */
		MPI_Bcast( &self->quitEventLoop,        1, MPI_INT, MASTER, context->communicator );
		MPI_Bcast( &self->interactive,          1, MPI_INT, MASTER, context->communicator );
		MPI_Bcast( &self->windowIsVisible,      1, MPI_INT, MASTER, context->communicator );
		MPI_Bcast( &self->hackNonInteractive,   1, MPI_INT, MASTER, context->communicator );

		if (self->hackNonInteractive) {
			lucCarbonWindow_Draw(self);
			break;
		}

		if (self->windowIsVisible) {
			lucCarbonWindow_Draw(self);
		}

		if ( !self->interactive )
			break;
	}

	/* Dump image */
	if ( iAmMaster )
		lucWindow_Dump( self, data );

	/* Clean Up */
	lucWindow_CleanUp( self, data );

	lucCarbonWindow_DestroyWindow( self, data );

	lucDebug_PrintFunctionEnd( self, 1 );

}

void _lucCarbonWindow_Destroy( void* window, void* data ) {}

void lucCarbonWindow_Draw( void* window ) {
	lucCarbonWindow*        self        = (lucCarbonWindow*) window;
	AbstractContext*        context     = self->currentContext;
	Bool                    iAmMaster   = ( context->rank == MASTER );
	
	lucDebug_PrintFunctionBegin( self, 2 );

	if ( self->windowIsInteractive )
		aglSetCurrentContext(self->graphicsContext);

	/* Setup fonts */
	_lucWindow_SetupGLRasterFont( self );
#if 0
	aglUseFont( self->graphicsContext, self->fontID, bold, self->fontSize, 32, 96, 2000 );
	glListBase(2000);
#endif

	lucWindow_Draw( self, self->currentContext );
	if ( iAmMaster )
		aglSwapBuffers(self->graphicsContext);
	
	lucDebug_PrintFunctionEnd( self, 2 );
}



void lucCarbonWindow_CreateWindow( void* carbonWindow, void* data ) {
	lucCarbonWindow*        self      = (lucCarbonWindow*) carbonWindow; 
	AbstractContext*        context   = self->currentContext;
	Bool                    iAmMaster = ( context->rank == MASTER );

	if ( iAmMaster && self->interactive ) {
		lucCarbonWindow_CreateInteractiveWindow( self, data );
	}
	else {
		lucCarbonWindow_CreateBackgroundWindow( self, data );
	}
}

void lucCarbonWindow_CreateInteractiveWindow( void* carbonWindow, void* data ) {
	lucCarbonWindow*        self = (lucCarbonWindow*) carbonWindow; 
	AGLPixelFormat          format;		/* OpenGL pixel format */
	WindowPtr               window;		/* Window */
	int                     winattrs;	/* Window attributes */
	Rect                    rect;		/* Rectangle definition */
	EventHandlerUPP	        handler;	/* Event handler */
	EventLoopTimerUPP       thandler;	/* Timer handler */
	EventLoopTimerRef       timer;		/* Timer for animating the window */
	ProcessSerialNumber     psn;		/* Process serial number */
	static EventTypeSpec    events[] =	/* Events we are interested in... */
			{
			  { kEventClassMouse, kEventMouseDown },
			  { kEventClassMouse, kEventMouseUp },
			  { kEventClassMouse, kEventMouseMoved },
			  { kEventClassMouse, kEventMouseDragged },
			  { kEventClassMouse, kEventMouseWheelMoved },
			  { kEventClassKeyboard, kEventRawKeyDown },
			  { kEventClassWindow, kEventWindowDrawContent },
			  { kEventClassWindow, kEventWindowShown },
			  { kEventClassWindow, kEventWindowHidden },
			  { kEventClassWindow, kEventWindowActivated },
			  { kEventClassWindow, kEventWindowDeactivated },
			  { kEventClassWindow, kEventWindowClose },
			  { kEventClassWindow, kEventWindowBoundsChanged }
			};
	static GLint 		attributes[] =	/* OpenGL attributes */
			{
			  AGL_RGBA,
			  AGL_GREEN_SIZE, 1,
			  AGL_DOUBLEBUFFER,
			  AGL_DEPTH_SIZE, 16,
			  AGL_NONE
			};

	lucDebug_PrintFunctionBegin( self, 1 );

	Journal_Firewall( 
			self->currentContext->rank == MASTER,
			Journal_MyStream( Error_Type, self ),
			"Error in func %s for %s '%s':\n"
			"Only the MASTER processor should be calling this function.\n", 
			__func__, self->type, self->name );

	/* Create the window...  */
	self->graphicsContext = NULL;
	self->windowIsVisible = False;

	SetRect(&rect, (int)self->offsetX, (int) self->offsetY, (int) (self->width + self->offsetX), (int) (self->height + self->offsetY) );

	winattrs = kWindowStandardHandlerAttribute | kWindowCloseBoxAttribute |
		     kWindowCollapseBoxAttribute | kWindowFullZoomAttribute |
		     kWindowResizableAttribute | kWindowLiveResizeAttribute;
	winattrs &= GetAvailableWindowAttributes(kDocumentWindowClass);

	CreateNewWindow(kDocumentWindowClass, winattrs, &rect, &window);
	SetWTitle(window, self->title);

	handler = NewEventHandlerUPP(lucCarbonWindow_EventHandler);
	InstallWindowEventHandler(window, handler,
				    sizeof(events) / sizeof(events[0]),
				    events, self, 0L);
	thandler = NewEventLoopTimerUPP((void (*)(EventLoopTimerRef, void *)) lucCarbonWindow_IdleFunc );
	InstallEventLoopTimer(GetMainEventLoop(), 0, 0, thandler, 0, &timer);

	/* This hack with the 'if' condition is to make the CarbonWindowing code not crash if it starts with an interactive window then starts to use a non-interactive window - if it was originally interactive it must continue to use the same windowing code, but we will just arrange so that it doesn't show up for the user */
	if ( ! self->hackNonInteractive ) {
		GetCurrentProcess(&psn);
		/* this is a secret undocumented Mac function that allows code that isn't part of a bundle to be a foreground operation */
		CPSEnableForegroundOperation( &psn ); 
		SetFrontProcess( &psn );

		/* Window operations - Not all of these are nessesary */
		DrawGrowIcon(window);
		ShowWindow(window);
		SetUserFocusWindow( window );
		BringToFront( window );
		ActivateWindow( window, true );
		SelectWindow(window);

		lucWindow_InteractionHelpMessage( self, Journal_MyStream( Info_Type, self ) );
	}

	/* Create the OpenGL context and bind it to the window.  */
	format = aglChoosePixelFormat(NULL, 0, attributes);
	self->graphicsContext = aglCreateContext(format, NULL);
	assert( self->graphicsContext );
	aglDestroyPixelFormat(format);
	
	aglSetDrawable(self->graphicsContext, GetWindowPort(window));

	self->windowIsInteractive = True;
	self->timer = timer;
	self->timerHandler = thandler;
	self->handler = handler;
	self->window = window;
}

/* Steps taken from http://gemma.apple.com/documentation/GraphicsImaging/Conceptual/OpenGL-MacProgGuide/OpenGLProg_MacOSX.pdf */

void lucCarbonWindow_CreateBackgroundWindow( void* carbonWindow, void* data ) {
	lucCarbonWindow*        self       = (lucCarbonWindow*) carbonWindow; 

	GLint                   numPixelFormats;
	CGLContextObj           contextObj;
	lucAlphaPixel*          memBuffer;
	CGLPixelFormatObj       pixelFormatObj;
	CGLPBufferObj           pBuffer;
	CGLError 				cglReturnValue;
	
	
	/* Step 1:
	 * Sets up an array of pixel format attributes 
	 * - an offscreen drawable object and a color buffer with a size of 32 bytes. 
	 * Note that the list must be terminated by NULL. */
	
/*	CGLPixelFormatAttribute attribs[] = 	{
		kCGLPFARemotePBuffer,
		kCGLPFAOffScreen,
		kCGLPFAColorSize, 32,
		NULL
	} ; */
	
	/* Note: the apple document referred to above was recently updated */
	
	CGLPixelFormatAttribute attribs[] = 	{
		kCGLPFAOffScreen,
		kCGLPFAColorSize, 32,
		0
	} ;
	

	/* Step 2 - Creates a pixel format object that has the specified renderer and buffer attributes. */
	cglReturnValue = CGLChoosePixelFormat (attribs, &pixelFormatObj, &numPixelFormats);
		
	if(cglReturnValue != kCGLNoError) 
		Journal_DPrintf( lucDebug,"Error setting up background window (Step 2): %s\n", CGLErrorString(cglReturnValue));
		
		

	/* Step 3 - Creates a CGL context using the newly created pixel format object. */
	cglReturnValue = CGLCreateContext (pixelFormatObj, NULL, &contextObj); /* Step 3 */
	if(cglReturnValue != kCGLNoError) 
		Journal_DPrintf( lucDebug,"Error setting up background window (Step 3): %s\n", CGLErrorString(cglReturnValue));	

	cglReturnValue = CGLDestroyPixelFormat (pixelFormatObj);
	if(cglReturnValue != kCGLNoError) 
		Journal_DPrintf( lucDebug,"Error setting up background window (Step 3): %s\n", CGLErrorString(cglReturnValue));	
	

	/* Step 4a - Sets the current context to the newly created offscreen CGL context. */
	cglReturnValue = CGLSetCurrentContext (contextObj);
	if(cglReturnValue != kCGLNoError) 
		Journal_DPrintf( lucDebug,"Error setting up background window (Step 4a): %s\n", CGLErrorString(cglReturnValue));
		
	/* Step 4b - Create Pixel Buffer */
	cglReturnValue = CGLCreatePBuffer( self->width, self->height, GL_TEXTURE_RECTANGLE_EXT, GL_RGBA, 0, &pBuffer );
	if(cglReturnValue != kCGLNoError) 
		Journal_DPrintf( lucDebug,"Error setting up background window (Step 4b): %s\n", CGLErrorString(cglReturnValue));
	
	cglReturnValue = CGLSetPBuffer( contextObj, pBuffer, 0, 0, 0 );
	if(cglReturnValue != kCGLNoError) 
		Journal_DPrintf( lucDebug,"Error setting up background window (Step 4b): %s\n", CGLErrorString(cglReturnValue));

	/* Step 5 - Allocates memory for the offscreen drawable object. */
	memBuffer = Memory_Alloc_Array( lucAlphaPixel, self->width * self->height, "Image Buffer" );

	/* Step 6 - 
	 * Binds the CGL context to the newly allocated offscreen memory buffer. 
	 * You need to specify the width and height of the offscreen buffer (in pixels), 
	 * the number of bytes per row, and a pointer to the block of memory you want to render the context into.
	 * The number of bytes per row must be at least the width times the bytes per pixels. */
	
	cglReturnValue = CGLSetOffScreen (contextObj, (GLsizei) self->width, (GLsizei) self->height, (GLsizei) self->width * sizeof(lucAlphaPixel), memBuffer);
	if(cglReturnValue != kCGLNoError) 
		Journal_DPrintf( lucDebug,"Error setting up background window (Step 6): %s", CGLErrorString(cglReturnValue));

	/* Set pointer to graphics context on my object */
	self->windowIsInteractive = False;
	self->graphicsContext = contextObj;
}


void lucCarbonWindow_DestroyWindow( void* carbonWindow, void* data ) {
	lucCarbonWindow*        self      = (lucCarbonWindow*) carbonWindow; 

	if ( self->windowIsInteractive ) {
		RemoveEventLoopTimer( self->timer );
		DisposeEventLoopTimerUPP( self->timerHandler ); 
		DisposeEventHandlerUPP( self->handler ); 
		DisposeWindow( self->window );

		aglSetDrawable( self->graphicsContext, 0 );
		aglSetCurrentContext( 0 );
		aglDestroyContext( self->graphicsContext );
		self->graphicsContext = NULL;
	}
	
	else {
				
		CGLContextObj contextObj = (CGLContextObj) self->graphicsContext;
		CGLPBufferObj pBuffer;
		CGLError cglReturnValue;
		GLenum face;
		GLint level;
		GLint screen;
	
		cglReturnValue = CGLGetPBuffer( contextObj, &pBuffer, &face, &level, &screen ) ;		
		if(cglReturnValue != kCGLNoError)
			Journal_DPrintf( lucDebug,"Error cleaning up background window: %s", CGLErrorString(cglReturnValue));
		
		
			CGLDestroyPBuffer (pBuffer);
			CGLSetCurrentContext (NULL);
			CGLClearDrawable (contextObj);
			CGLDestroyContext (contextObj);		
	}
}

static pascal OSStatus lucCarbonWindow_EventHandler(EventHandlerCallRef nextHandler, EventRef event, void* userData) {
	lucCarbonWindow*        self = (lucCarbonWindow*) userData;
	UInt32                  kind;			/* Kind of event */
	Rect                    rect;			/* New window size */
	EventMouseButton        button;			/* Mouse button */
	static Point            point;			/* Mouse position */
	static Pixel_Index      startx      = 0;
	static Pixel_Index      starty      = 0;
	static lucMouseButton   whichButton = 0;
	static lucMouseState    mouseState  = 0;
	EventClass              eventClass;
	Pixel_Index             xPos, yPos;
	SInt32					delta;

	assert( self );
	
	kind = GetEventKind(event);
	eventClass = GetEventClass( event );

	if (eventClass == kEventClassWindow) {
		switch (kind) {
			case kEventWindowDrawContent:
				if (self->windowIsVisible && self->graphicsContext) 
					lucCarbonWindow_Draw(self); 
				break; 
			case kEventWindowBoundsChanged:
			       GetEventParameter(event, kEventParamCurrentBounds, typeQDRectangle, NULL, sizeof(Rect), NULL, &rect); 
			       
			       if (self->graphicsContext)
				      aglUpdateContext(self->graphicsContext);

			       self->offsetX = (Pixel_Index) rect.left;
			       self->offsetY = (Pixel_Index) rect.top; 
			       lucWindow_Resize( self, (Pixel_Index) (rect.right - rect.left), (Pixel_Index) (rect.bottom - rect.top) );
			       
			       if (self->windowIsVisible && self->graphicsContext)
				      lucCarbonWindow_Draw(self); 
			       break; 
			case kEventWindowShown: 
			       self->windowIsVisible = True; 
			       
			       if (self->graphicsContext) 
				       lucCarbonWindow_Draw(self); 
			       break; 
			case kEventWindowHidden: 
			       self->windowIsVisible = False; 
			       break; 
			
			case kEventWindowClose: 
				   lucWindow_QuitEventLoop( self );
			       break; 
		} 
	} 
	else if ( eventClass == kEventClassKeyboard ) {
		if ( kind == kEventRawKeyDown ) {
			char key;
			GetEventParameter(event, kEventParamKeyMacCharCodes, typeChar, NULL, sizeof(char), NULL, &key); 
			lucCarbonWindow_GetPixelIndicies( self, &point, &xPos, &yPos );

			lucWindow_KeyboardEvent( self, key, xPos, yPos );
			/* HACK */
			if ( key == 'i' ) {
				self->interactive = true;
				self->hackNonInteractive = true;
				self->quitEventLoop = true;
			}
		}
	}
	else if ( eventClass == kEventClassMouse ) { 
		switch (kind) { 
			case kEventMouseMoved:
				GetEventParameter(event, kEventParamMouseLocation, typeQDPoint, NULL, sizeof(Point), NULL, &point); 
				return (CallNextEventHandler(nextHandler, event)); 			
			case kEventMouseDown:  case kEventMouseUp:
				GetEventParameter(event, kEventParamMouseButton, typeMouseButton, NULL, sizeof(EventMouseButton), NULL, &button); 
				GetEventParameter(event, kEventParamMouseLocation, typeQDPoint, NULL, sizeof(Point), NULL, &point); 
				
				/* Make sure that mouse click was within window */
				if ( lucCarbonWindow_IsPointOutsideWindow( self, &point )) 
					return (CallNextEventHandler(nextHandler, event)); 
					
				lucCarbonWindow_GetPixelIndicies( self, &point, &startx, &starty );					

				whichButton = ( button == kEventMouseButtonTertiary ? lucMiddleButton :
						button == kEventMouseButtonSecondary ? lucRightButton : lucLeftButton );

				mouseState = (kind == kEventMouseDown ? lucButtonPress : lucButtonRelease );

				lucWindow_MouseClick( self, whichButton, mouseState, startx, starty);
				break; 
			case kEventMouseDragged : 
				GetEventParameter(event, kEventParamMouseLocation, typeQDPoint, NULL, sizeof(Point), NULL, &point); 
				
				/* Make sure that mouse click was within window */
				if ( lucCarbonWindow_IsPointOutsideWindow( self, &point )) 
					return (CallNextEventHandler(nextHandler, event)); 
				
				lucCarbonWindow_GetPixelIndicies( self, &point, &xPos, &yPos );
				lucWindow_MouseMotion(self, whichButton, xPos, yPos, startx, starty);
				startx = xPos;
				starty = yPos;
				break; 
			case kEventMouseWheelMoved:
				GetEventParameter(event, kEventParamMouseLocation, typeQDPoint, NULL, sizeof(Point), NULL, &point); 
				GetEventParameter(event, kEventParamMouseWheelDelta, typeLongInteger, NULL, sizeof(delta), NULL, &delta );
				
				/* Make sure that mouse click was within window */
				if ( lucCarbonWindow_IsPointOutsideWindow( self, &point )) 
					return (CallNextEventHandler(nextHandler, event)); 
					
				lucCarbonWindow_GetPixelIndicies( self, &point, &startx, &starty );					

				whichButton = ( delta > 0 ? lucWheelUp : lucWheelDown );

				lucWindow_MouseClick( self, whichButton, lucButtonPress, startx, starty);
				break;
			default : 
				return (CallNextEventHandler(nextHandler, event)); 
		} 
	} 
	else {
		Journal_Firewall( False, Journal_MyStream( Error_Type, self ), 
				"In func '%s' - Cannot understand event class type '%d'\n", __func__, eventClass );
	}
	
	/* Return whether we handled the event...  */ 
	return  noErr;
}

void lucCarbonWindow_IdleFunc() {
  QuitApplicationEventLoop();
}

Bool lucCarbonWindow_IsPointOutsideWindow( lucCarbonWindow* self, Point* point ) {
	return point->v < self->offsetY || (point->v > (self->offsetY + self->height - 8) && point->h > (self->offsetX + self->width - 8));
}

void lucCarbonWindow_GetPixelIndicies( lucCarbonWindow* self, Point* point, Pixel_Index *xPos, Pixel_Index *yPos ) {
	*xPos = (Pixel_Index) point->h - self->offsetX;
	*yPos = self->height - ((Pixel_Index) point->v - self->offsetY);
}

#endif
