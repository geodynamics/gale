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
lucVTKWindow* _lucVTKWindow_New( 
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
		lucWindow_DisplayFunction*						_displayWindow,	
		lucWindow_EventsWaitingFunction*				_eventsWaiting,	
		lucWindow_EventProcessorFunction*				_eventProcessor,	
		lucWindow_ResizeFunction*						_resizeWindow,	
		Name                                            name ) 
{
	lucVTKWindow*					self;

	/* Allocate memory */
	assert( sizeOfSelf >= sizeof(lucVTKWindow) );
	self = (lucVTKWindow*) _lucWindow_New( 
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
	return (void*) _lucVTKWindow_New(
		sizeof(lucVTKWindow),
		lucVTKWindow_Type,
		_lucVTKWindow_Delete,
		_lucVTKWindow_Print,
		NULL,
		_lucVTKWindow_DefaultNew,
		_lucVTKWindow_Construct,
		_lucVTKWindow_Build,
		_lucVTKWindow_Initialise,
		_lucVTKWindow_Execute,
		_lucVTKWindow_Destroy,
		lucWindow_Display,	/* Use parent class default implementations */
		lucWindow_EventsWaiting,
		lucWindow_EventProcessor,
		lucWindow_Resize,
		name );
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
