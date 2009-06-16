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

#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>

#include <AGL/agl.h>
#include <OpenGL/OpenGL.h>

#include <glucifer/Base/Base.h>

#include "types.h"
#include "CarbonWindow.h"

#include <assert.h>

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
		lucWindow_DisplayFunction*						   _displayWindow,	
		lucWindow_EventsWaitingFunction*				   _eventsWaiting,	
		lucWindow_EventProcessorFunction*				   _eventProcessor,	
		lucWindow_ResizeFunction*						   _resizeWindow,	
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
			_displayWindow,		
			_eventsWaiting,		
			_eventProcessor,	
			_resizeWindow,		
			name );
	
	return self;
}

void _lucCarbonWindow_Init( 
		lucCarbonWindow*                                   self,
		Pixel_Index                                        offsetX,
		Pixel_Index                                        offsetY,
		double                                             maxIdleTime ) 
{
	self->offsetX = offsetX;
	self->offsetY = offsetY;
}

void _lucCarbonWindow_Delete( void* window ) {
	_lucWindow_Delete( window );
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
		_lucCarbonWindow_Display,
		_lucCarbonWindow_EventsWaiting,
		_lucCarbonWindow_EventProcessor,
		_lucCarbonWindow_Resize,
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

void _lucCarbonWindow_Build( void* window, void* data ) {
	/* Run the parent function to build window... */
	_lucWindow_Build(window, data);	
}

void _lucCarbonWindow_Initialise( void* window, void* data ) {
	/* OK: Moved from ExecuteInteractive/Offscreen */
	lucCarbonWindow*     self      = (lucCarbonWindow*)window;
		
	/* Create the window...  */
	if ( self->isMaster && self->interactive ) {
		lucCarbonWindow_CreateInteractiveWindow( self );
	}
	else {
		lucCarbonWindow_CreateBackgroundWindow( self );
	}	

	/* Run the parent function to init window... */
	_lucWindow_Initialise(window, data);	
	
	/* Run the parent function to display window... */
	lucWindow_Display(window);	
}


void _lucCarbonWindow_Execute( void* window, void* data ) {

	/* Post a dummy event, hack to ensure waiting events are processed in continuous mode */
	if (_lucCarbonWindow_EventsWaiting(window) == 0)
	{
		EventRef dummyEvent;
		OSStatus        err;
		err = CreateEvent(NULL, kEventClassWindow, kEventWindowDrawContent, GetCurrentEventTime(), kEventAttributeNone, &dummyEvent);
		if (err == noErr) err = PostEventToQueue(GetMainEventQueue(), dummyEvent, kEventPriorityHigh);
	}
	
	/* Run the parent function to execute window... */
	_lucWindow_Execute(window, data);	
}

void _lucCarbonWindow_Destroy( void* window, void* data ) {
	lucCarbonWindow*        self = (lucCarbonWindow*) window; 

	/* Run the parent function to destroy window... */
	_lucWindow_Destroy(window, data);	
	
	/* Destroy the window...  */
	if ( self->window )
		lucCarbonWindow_DestroyInteractiveWindow( self );
	else 
		lucCarbonWindow_DestroyBackgroundWindow( self );
}

/* Window Virtuals */
void _lucCarbonWindow_Display( void* window ) {
	lucCarbonWindow*        self = (lucCarbonWindow*) window; 
	
	/* Run the parent function to display window... */
	lucWindow_Display(window);	

	/* Swap buffers if interactive */
	if ( self->isMaster && self->interactive)
		aglSwapBuffers(self->graphicsContext);
}

int _lucCarbonWindow_EventsWaiting( void* window ) {
	EventQueueRef evq = GetMainEventQueue();
	return GetNumEventsInQueue(evq);
}

Bool _lucCarbonWindow_EventProcessor( void* window ) {

	RunApplicationEventLoop();
	/* Returns true when redisplay required */
	return True;
}

void _lucCarbonWindow_Resize( void* window ) {
	lucCarbonWindow*        self = (lucCarbonWindow*) window; 

    if (self->interactive) 
	{
		if (self->isMaster)
		{
			/* Reset context after resize */
			aglSetCurrentContext (self->graphicsContext);
			aglUpdateContext (self->graphicsContext);
		}
		else
		{
			/* Master window resized? Create new background window of required size */
			lucCarbonWindow_DestroyBackgroundWindow( self );
			lucCarbonWindow_CreateBackgroundWindow( self );
		}
	}
	
	/* Run the parent function to resize window... */
	lucWindow_Resize(window);	
}

void lucCarbonWindow_CreateInteractiveWindow( void* window ) {
	lucCarbonWindow*        self = (lucCarbonWindow*) window; 
	AGLPixelFormat          format;		/* OpenGL pixel format */
	int                     winattrs;	/* Window attributes */
	Rect                    rect;		/* Rectangle definition */
	ProcessSerialNumber     psn;		/* Process serial number */
	static EventTypeSpec    events[] =	/* Events we are interested in... */
			{
			  { kEventClassMouse, kEventMouseDown },
			  { kEventClassMouse, kEventMouseUp },
			  { kEventClassMouse, kEventMouseDragged },
			  { kEventClassMouse, kEventMouseWheelMoved },
			  { kEventClassKeyboard, kEventRawKeyDown },
			  { kEventClassWindow, kEventWindowDrawContent },
			  { kEventClassWindow, kEventWindowClose },
			  { kEventClassWindow, kEventWindowBoundsChanged },
			  { kEventClassWindow, kEventWindowResizeCompleted},
			};
	static GLint 		attributes[] =	/* OpenGL attributes */
			{
			  AGL_RGBA,
			  AGL_GREEN_SIZE, 1,
			  AGL_DOUBLEBUFFER,
			  AGL_DEPTH_SIZE, 16,
              AGL_STENCIL_SIZE, 1,
			  AGL_NONE
			};

	lucDebug_PrintFunctionBegin( self, 1 );

	/* Create the window...  */
	self->graphicsContext = NULL;

	SetRect(&rect, (int)self->offsetX, (int) self->offsetY, (int) (self->width + self->offsetX), (int) (self->height + self->offsetY) );

	winattrs = kWindowStandardHandlerAttribute | kWindowCloseBoxAttribute |
		     kWindowCollapseBoxAttribute | kWindowFullZoomAttribute |
		     kWindowResizableAttribute | kWindowLiveResizeAttribute;
	winattrs &= GetAvailableWindowAttributes(kDocumentWindowClass);

	CreateNewWindow(kDocumentWindowClass, winattrs, &rect, &self->window);
	SetWTitle(self->window, (unsigned char*)self->title);
//	SetWindowTitleWithCFString(self->window, CFSTR(self->title));

	self->handler = NewEventHandlerUPP(lucCarbonWindow_EventHandler);
	InstallWindowEventHandler(self->window, self->handler, sizeof(events) / sizeof(events[0]), events, self, 0L);

    self->timerHandler = NewEventLoopIdleTimerUPP(lucCarbonWindow_IdleTimer);
    InstallEventLoopIdleTimer(GetMainEventLoop(), kEventDurationSecond * 2, kEventDurationSecond * 1, self->timerHandler, self, &self->timer);

	GetCurrentProcess(&psn);
	/* this is a secret undocumented Mac function that allows code that isn't part of a bundle to be a foreground operation */
	CPSEnableForegroundOperation( &psn ); 
	SetFrontProcess( &psn );

	/* Show Window */
	ShowWindow(self->window);
	SelectWindow(self->window);

	/* Create the OpenGL context and bind it to the window.  */
	format = aglChoosePixelFormat(NULL, 0, attributes);
	self->graphicsContext = aglCreateContext(format, NULL);
	assert( self->graphicsContext );
	aglDestroyPixelFormat(format);
	
	aglSetDrawable(self->graphicsContext, GetWindowPort(self->window));
	aglSetCurrentContext(self->graphicsContext);

}

/* Steps taken from http://gemma.apple.com/documentation/GraphicsImaging/Conceptual/OpenGL-MacProgGuide/OpenGLProg_MacOSX.pdf */

void lucCarbonWindow_CreateBackgroundWindow( void* window ) {
	lucCarbonWindow*        self       = (lucCarbonWindow*) window; 

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
	
	/* Note: the apple document referred to above was recently updated */
	
	CGLPixelFormatAttribute attribs[] = 	{
		kCGLPFAOffScreen,
		kCGLPFAColorSize, 32,
		kCGLPFADepthSize, 16,
		kCGLPFAStencilSize, 1,
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
	self->graphicsContext = contextObj;
	
	self->window = 0;
}

void lucCarbonWindow_DestroyInteractiveWindow( void* window ) {
	lucCarbonWindow*        self      = (lucCarbonWindow*) window; 

	DisposeEventHandlerUPP( self->handler ); 
    if (self->timer != NULL) RemoveEventLoopTimer( self->timer );
    DisposeEventLoopIdleTimerUPP( self->timerHandler ); 
	DisposeWindow( self->window );

	aglSetDrawable( self->graphicsContext, 0 );
	aglSetCurrentContext( 0 );
	aglDestroyContext( self->graphicsContext );
	self->graphicsContext = NULL;
}

void lucCarbonWindow_DestroyBackgroundWindow( void* window ) {
	lucCarbonWindow*        self      = (lucCarbonWindow*) window; 

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

static pascal OSStatus lucCarbonWindow_EventHandler(EventHandlerCallRef nextHandler, EventRef event, void* userData) {
	lucCarbonWindow*        self = (lucCarbonWindow*) userData;
	UInt32                  kind;			/* Kind of event */
	Rect                    rect;			/* New window size */
	EventMouseButton        button;			/* Mouse button */
	static Point            point;			/* Mouse position */
	static lucMouseButton   whichButton = 0;
	static lucMouseState    mouseState  = 0;
	static Pixel_Index      width, height;
	EventClass              eventClass;
	Pixel_Index             xPos, yPos;
	SInt32					delta;

	assert( self );
	
	kind = GetEventKind(event);
	eventClass = GetEventClass( event );

	if (eventClass == kEventClassWindow) {
		switch (kind) {
			case kEventWindowDrawContent:
				break; 
			case kEventWindowResizeCompleted:
				lucWindow_SetSize( self, width, height );
				break;
			case kEventWindowBoundsChanged:
			{
				GetEventParameter(event, kEventParamCurrentBounds, typeQDRectangle, NULL, sizeof(Rect), NULL, &rect); 
				width =  (Pixel_Index) (rect.right - rect.left);
				height = (Pixel_Index) (rect.bottom - rect.top);
				self->offsetX = (Pixel_Index) rect.left;
				self->offsetY = (Pixel_Index) rect.top; 
				break; 
			}
			case kEventWindowClose: 
				lucWindow_ToggleApplicationQuit(self);
				self->quitEventLoop = true;
			    break; 
		} 
	} 
	else if ( eventClass == kEventClassKeyboard ) {
		if ( kind == kEventRawKeyDown ) {
			char key;
			GetEventParameter(event, kEventParamKeyMacCharCodes, typeChar, NULL, sizeof(char), NULL, &key); 
			lucCarbonWindow_GetPixelIndicies( self, &point, &xPos, &yPos );

			lucWindow_KeyboardEvent( self, key, xPos, yPos );

			if ( !self->interactive )
			{
				/* Hide Window */
				HideWindow(self->window);
				self->quitEventLoop = true;
			}
		}
	}
	else if ( eventClass == kEventClassMouse ) { 
		switch (kind) { 
			case kEventMouseDown:  case kEventMouseUp:
				GetEventParameter(event, kEventParamMouseButton, typeMouseButton, NULL, sizeof(EventMouseButton), NULL, &button); 
				GetEventParameter(event, kEventParamMouseLocation, typeQDPoint, NULL, sizeof(Point), NULL, &point); 
				
				/* Make sure that mouse click was within window */
				if ( lucCarbonWindow_IsPointOutsideWindow( self, &point )) 
					return (CallNextEventHandler(nextHandler, event)); 
					
				lucCarbonWindow_GetPixelIndicies( self, &point, &xPos, &yPos );

				whichButton = ( button == kEventMouseButtonTertiary ? lucMiddleButton :
						button == kEventMouseButtonSecondary ? lucRightButton : lucLeftButton );

				mouseState = (kind == kEventMouseDown ? lucButtonPress : lucButtonRelease );

				lucWindow_MouseClick( self, whichButton, mouseState, xPos, yPos);
				break; 
			case kEventMouseDragged : 
				GetEventParameter(event, kEventParamMouseLocation, typeQDPoint, NULL, sizeof(Point), NULL, &point); 
				
				/* Make sure that mouse click was within window */
				if ( lucCarbonWindow_IsPointOutsideWindow( self, &point )) 
					return (CallNextEventHandler(nextHandler, event)); 
				
				lucCarbonWindow_GetPixelIndicies( self, &point, &xPos, &yPos );
				lucWindow_MouseMotion(self, whichButton, xPos, yPos); 
				break; 
			case kEventMouseWheelMoved:
				GetEventParameter(event, kEventParamMouseLocation, typeQDPoint, NULL, sizeof(Point), NULL, &point); 
				GetEventParameter(event, kEventParamMouseWheelDelta, typeLongInteger, NULL, sizeof(delta), NULL, &delta );
				
				/* Make sure that mouse click was within window */
				if ( lucCarbonWindow_IsPointOutsideWindow( self, &point )) 
					return (CallNextEventHandler(nextHandler, event)); 

				lucCarbonWindow_GetPixelIndicies( self, &point, &xPos, &yPos );

				whichButton = ( delta > 0 ? lucWheelUp : lucWheelDown );

				lucWindow_MouseClick( self, whichButton, lucButtonPress, xPos, yPos);
				break;
			default : 
				return (CallNextEventHandler(nextHandler, event)); 
		} 
	} 
	else {
		Journal_Firewall( False, Journal_MyStream( Error_Type, self ), 
				"In func '%s' - Cannot understand event class type '%d'\n", __func__, eventClass );
	}

	/* Break out of event processing if no more events */
	if (self->quitEventLoop || _lucCarbonWindow_EventsWaiting(self) == 0)
		QuitApplicationEventLoop();

	return  noErr;
}

pascal void lucCarbonWindow_IdleTimer(EventLoopTimerRef inTimer, EventLoopIdleTimerMessage inState, void * inUserData) {
		
	lucCarbonWindow*        self = (lucCarbonWindow*) inUserData;

	/* idle timeout check */
	if (inState == kEventLoopIdleTimerIdling)
	{
		if (self->interactive)
			lucWindow_IdleCheck(inUserData);
		else
		{
			RemoveEventLoopTimer( self->timer );
			self->timer = NULL;
		}

		/* Drop out of event loop if in continuous mode */
		if (self->continuous) QuitApplicationEventLoop();
	}
	else
	{
		/* Reset idle timer */
		lucWindow_IdleReset(inUserData);
	}
}

Bool lucCarbonWindow_IsPointOutsideWindow( lucCarbonWindow* self, Point* point ) {
	return point->v < self->offsetY || (point->v > (self->offsetY + self->height - 8) && point->h > (self->offsetX + self->width - 8));
}

void lucCarbonWindow_GetPixelIndicies( lucCarbonWindow* self, Point* point, Pixel_Index *xPos, Pixel_Index *yPos ) {
	*xPos = (Pixel_Index) point->h - self->offsetX;
	*yPos = self->height - ((Pixel_Index) point->v - self->offsetY);
}

#endif
