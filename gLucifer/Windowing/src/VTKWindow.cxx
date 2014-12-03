#ifdef HAVE_VTK

#include <vtkRenderer.h>
#include <vtkXRenderWindowInteractor.h>
#include <vtkRenderWindow.h>

#ifdef __cplusplus
extern "C" {
#endif

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <glucifer/Base/Base.h>
#include "types.h"


#include "VTKWindow.h"

#include <assert.h>
#ifdef __APPLE__
#include <Power.h>
#endif


#ifndef MASTER
	#define MASTER 0
#endif

/* Textual name of this class */
const Type lucVTKWindow_Type = "lucVTKWindow";

/* Creation implementation / Virtual constructor */
lucVTKWindow* _lucVTKWindow_New(  LUCVTKWINDOW_DEFARGS  ) 
{
	lucVTKWindow*					self;

	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(lucVTKWindow) );
	self = (lucVTKWindow*) _lucWindow_New(  LUCWINDOW_PASSARGS  );
	
	return self;
}

void _lucVTKWindow_Init( 
		lucVTKWindow*                                      self,
		Name                                               xFontName,
		Name                                               host,
		unsigned int                                       displayNumber,
		unsigned int                                       displayScreen,
		double                                             maxIdleTime ) 
{

	/* Setup display name */
	Stg_asprintf( &self->displayName, "%s:%u.%u", host, displayNumber, displayScreen );
	self->host          = StG_Strdup( host );
	self->displayNumber = displayNumber;
	self->displayScreen = displayScreen;
	self->maxIdleTime   = maxIdleTime;

	/* Get font from VTK */
	//self->xFontName = xFontName;
}

void _lucVTKWindow_Delete( void* window ) {
	lucVTKWindow*  self = (lucVTKWindow*)window;

	lucVTKWindow_CloseBackgroundWindow( self );
	lucVTKWindow_CloseInteractiveWindow( self );

	Memory_Free( self->displayName );
	Memory_Free( self->host );

	_lucWindow_Delete( self );
}

void _lucVTKWindow_Print( void* window, Stream* stream ) {
	lucVTKWindow*  self = (lucVTKWindow*)window;

	_lucWindow_Print( self, stream );
}

void* _lucVTKWindow_Copy( void* window, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap) {
	lucVTKWindow*  self = (lucVTKWindow*)window;
	lucVTKWindow* newWindow;

	newWindow = (lucVTKWindow*)_lucWindow_Copy( self, dest, deep, nameExt, ptrMap );

	/* TODO */
	abort();

	return (void*) newWindow;
}


void* _lucVTKWindow_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                              _sizeOfSelf = sizeof(lucVTKWindow);
	Type                                                      type = lucVTKWindow_Type;
	Stg_Class_DeleteFunction*                              _delete = _lucVTKWindow_Delete;
	Stg_Class_PrintFunction*                                _print = _lucVTKWindow_Print;
	Stg_Class_CopyFunction*                                  _copy = NULL;
	Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _lucVTKWindow_DefaultNew;
	Stg_Component_ConstructFunction*                    _construct = _lucVTKWindow_Construct;
	Stg_Component_BuildFunction*                            _build = _lucVTKWindow_Build;
	Stg_Component_InitialiseFunction*                  _initialise = _lucVTKWindow_Initialise;
	Stg_Component_ExecuteFunction*                        _execute = _lucVTKWindow_Execute;
	Stg_Component_DestroyFunction*                        _destroy = _lucVTKWindow_Destroy;
	lucWindow_DisplayFunction*                      _displayWindow = lucWindow_Display;
	lucWindow_EventsWaitingFunction*                _eventsWaiting = lucWindow_EventsWaiting;
	lucWindow_EventProcessorFunction*              _eventProcessor = lucWindow_EventProcessor;
	lucWindow_ResizeFunction*                        _resizeWindow = lucWindow_Resize;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

	return (void*) _lucVTKWindow_New(  LUCVTKWINDOW_PASSARGS  );
}

void _lucVTKWindow_Construct( void* window, Stg_ComponentFactory* cf ){
	lucVTKWindow*  self = (lucVTKWindow*)window;

	/* Construct Parent */
	_lucWindow_Construct( self, cf , data); 
				
	_lucVTKWindow_Init( 
			self,
			Stg_ComponentFactory_GetString( cf, self->name, "xFontName", "-adobe-courier-bold-r-normal--14-140-75-75-m-90-iso8859-1" ),
			Stg_ComponentFactory_GetString( cf, self->name, "host", "localhost" ),
			Stg_ComponentFactory_GetUnsignedInt( cf, self->name, "displayNumber", 0 ),
			Stg_ComponentFactory_GetUnsignedInt( cf, self->name, "displayScreen", 0 ),
			Stg_ComponentFactory_GetDouble( cf, self->name, "maxIdleTime", 600.0 ) );
			
			
}

void _lucVTKWindow_Build( void* window, void* data ) {}
void _lucVTKWindow_Initialise( void* window, void* data ) {}

void _lucVTKWindow_Execute( void* window, void* data ) {
	lucVTKWindow*     self      = (lucVTKWindow*)window;
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
			lucVTKWindow_CreateInteractiveWindow( self ); 
		}
		else {
			Journal_DPrintf( lucDebug, "Using background window on slave processor.\n");
			lucVTKWindow_CreateBackgroundWindow( self );
		}

		lucVTKWindow_MakeCurrent( self );
		lucVTKWindow_EventLoop( self, context );
	}
	else {
		/* Open window */
		lucVTKWindow_CreateBackgroundWindow( self );
		
		lucVTKWindow_MakeCurrent( self );
		lucWindow_Draw( self, context );
	}
	
	/* Dump image */
	lucWindow_Dump( self, context );

	lucWindow_CleanUp( self, context );
	lucVTKWindow_CloseInteractiveWindow( self );

	lucDebug_PrintFunctionEnd( self, 1 );
}

void _lucVTKWindow_Destroy( void* window, void* data ) {}

Bool lucVTKWindow_EventLoop( void* window, AbstractContext* context) {
	lucVTKWindow*  self             = (lucVTKWindow*)window;
	/*Display*       dpy              = self->display;
	KeySym         ks;
	XEvent         event;
	Atom           wmDeleteWindow;
	Bool           loop             = True;
	unsigned int   button           = 0;
	Pixel_Index    startx           = 0;
	Pixel_Index    starty           = 0;
	Bool           iAmMaster        = (context->rank == MASTER);
	Bool           handledEvent     = False;
	double         maxIdleTime      = MPI_Wtime() + self->maxIdleTime;

	lucDebug_PrintFunctionBegin( self, 1 );

	while (True) {
		if ( iAmMaster ) {
			while ( XEventsQueued(dpy, QueuedAfterFlush ) ) {
				handledEvent = True;
				XNextEvent(dpy, &event);
				switch (event.type) {
					case ConfigureNotify:
						//Reshape( windowingContextExt, event.xconfigure.width, event.xconfigure.height);
					
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
						loop = lucWindow_KeyboardEvent( self, ks, event.xkey.x, self->height - event.xkey.y);
						break;
					case ClientMessage:
						if (event.xclient.data.l[0] == wmDeleteWindow)
						loop = False;
						break;
					default:
						lucWindow_SetViewportNeedsToDrawFlag( self, True );
				}
			}

			if ( MPI_Wtime() > maxIdleTime ) {
				Journal_Printf( lucError, 
						"Error in func '%s'" 
						" - Interactive window '%s' open for too long (over %g seconds) without interaction.\n",
						__func__, self->name, self->maxIdleTime );
				abort();
			}
		}

		MPI_Bcast( &handledEvent, 1, MPI_INT, MASTER, context->communicator );
		if ( !handledEvent )
			continue;

		MPI_Bcast( &loop, 1, MPI_INT, MASTER, context->communicator );
		if ( loop == False )
			break;

		Journal_DPrintfL( lucDebug, 2, "In func %s: Event handled - Now redrawing.\n", __func__ );
		lucWindow_Draw( self, context );
		lucVTKWindow_SwapBuffers( self, context );
		handledEvent = False;
		maxIdleTime = MPI_Wtime() + self->maxIdleTime;
		
	}
*/
	lucDebug_PrintFunctionEnd( self, 1 );
	return True;
}

void lucVTKWindow_SetupFonts( void* window ) {
	lucVTKWindow*  self      = (lucVTKWindow*)window;
	
	lucDebug_PrintFunctionBegin( self, 2 );

	/* Tell GLX to use this font */
//	glXUseXFont( self->font->fid, 32, 96, 2000+32 );
	//glListBase(2000);
	
	lucDebug_PrintFunctionEnd( self, 2 );
}

void lucVTKWindow_SwapBuffers( void* window, AbstractContext* context ) {
	lucVTKWindow*  self      = (lucVTKWindow*)window;
	
	lucDebug_PrintFunctionBegin( self, 2 );

	//if (self->doubleBuffer && context->rank == MASTER) 
		//glXSwapBuffers(self->display, self->win);
	
	lucDebug_PrintFunctionEnd( self, 2 );
}



void lucVTKWindow_MakeCurrent( void* window )  {
	lucVTKWindow*  self         = (lucVTKWindow*)window;
	
	lucDebug_PrintFunctionBegin( self, 1 );

	//if ( self->backgroundWindowOpen ) {
		//glXMakeCurrent( self->display, self->glxpmap, self->glxcontext);
	//}
	//else if ( self->interactiveWindowOpen ){
		//glXMakeCurrent( self->display, self->win, self->glxcontext);
		//XMapWindow( self->display, self->win );
		//glXSwapBuffers( self->display, self->win );
	//}

	lucVTKWindow_SetupFonts( self );
	
	lucDebug_PrintFunctionEnd( self, 1 );
}

Bool lucVTKWindow_CreateDisplay( void* window )  {
	lucVTKWindow*  self         = (lucVTKWindow*)window;
	
	vtkRenderWindow* rendererWindow = (vtkRenderWindow*) self->rendererWindow;
	
	rendererWindow = vtkRenderWindow::New();
	//window->DoubleBufferOn();
	rendererWindow->SetSize( self->width, self->height );   
	rendererWindow->LineSmoothingOn();
	
	rendererWindow->SetSubFrames(self->subframes);
	
   rendererWindow->SetWindowName( self->title );

	if (self->interactive) {
		rendererWindow->OffScreenRenderingOff();
		if (self->fullscreen)
	    		rendererWindow->FullScreenOn();
		else 
	    		rendererWindow->FullScreenOff();
	}
	else 
		rendererWindow->OffScreenRenderingOn();
		
   vtkXRenderWindowInteractor* interactor = (vtkXRenderWindowInteractor*) self->interactor ;
	 interactor = vtkXRenderWindowInteractor::New();
	 rendererWindow->SetInteractor(interactor);    

	return True;
}

Bool lucVTKWindow_FindFont( void* window )  {
	
	return True;
}

void lucVTKWindow_CreateBackgroundWindow( void* window )  {
	//lucVTKWindow*  self         = (lucVTKWindow*)window;
	
}


void lucVTKWindow_CreateInteractiveWindow( void* window ) {
	
}	

Colormap lucVTKWindow_GetShareableColormap( lucVTKWindow* self ) {
	Colormap cmap;
	return cmap;
}

void lucVTKWindow_CloseInteractiveWindow( lucVTKWindow* self ) {
	
	
}
	
void lucVTKWindow_CloseBackgroundWindow( lucVTKWindow* self ) {


}

void lucVTKWindow_CloseDisplay( lucVTKWindow* self ) {
	

}

#ifdef __cplusplus
}
#endif

#endif


