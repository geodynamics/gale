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
** $Id: Window.c 791 2008-09-01 02:09:06Z JulianGiordani $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>

#include "types.h"
#include "ColourMap.h"
#include "X11Colours.h"
#include "Window.h"
#include "Viewport.h"
#include "Camera.h"
#include "OutputFormat.h"
#include "OutputFormat_Register.h"
#include "Init.h"
#include "RenderingEngine.h"
#include "Light.h"
#include "Light_Register.h"

#include "WindowInteraction.h"
#include "WindowInteraction_Register.h"

#include "ViewportInfo.h"

#include <ctype.h>
#include <assert.h>
#include <string.h>

/* ASCII Characters */
#define ESCAPE 27

#ifndef MASTER
	#define MASTER 0
#endif

const Type lucWindow_Type = "lucWindow";

MPI_Datatype lucWindow_MPI_Datatype;
MPI_Datatype lucViewportInfo_MPI_Datatype;

lucWindow* _lucWindow_New(
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
	lucWindow*    self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( sizeOfSelf >= sizeof(lucWindow) );
	self = (lucWindow*) _Stg_Component_New( 
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
			name, 
			NON_GLOBAL );

	/* Virtual functions */
	self->_displayWindow = _displayWindow;
	self->_eventsWaiting = _eventsWaiting;
	self->_eventProcessor = _eventProcessor;
	self->_resizeWindow = _resizeWindow;
	
	return self;
}

void _lucWindow_Init( 
		lucWindow*                                         self,
		lucRenderingEngine*                                renderingEngine,
		lucViewportInfo*                                   viewportInfoList,
		Viewport_Index                                     viewportCount,
		lucOutputFormat**                                  outputFormatList,
		OutputFormat_Index                                 outputFormatCount,
		lucWindowInteraction**                             windowInteractionList, 
		WindowInteraction_Index                            windowInteractionCount,
		AbstractContext*                                   context,
		Pixel_Index                                        width,
		Pixel_Index                                        height,
		Name                                               backgroundColourName,
		Bool                                               interactive,
		Bool                                               continuous,
		Bool                                               isTimedOut,
		double                                             maxIdleTime ) 
{
	OutputFormat_Index   outputFormat_I;
	WindowInteraction_Index windowInteraction_I;

	self->renderingEngine = renderingEngine;
	self->width = width;
	self->height = height;
    self->resized = False;
	self->interactive = interactive;
	self->continuous = continuous; 

	self->viewportInfoList = Memory_Alloc_Array( lucViewportInfo, viewportCount, "viewport info Array" );
	memcpy( self->viewportInfoList, viewportInfoList, viewportCount * sizeof( lucViewportInfo ) );
	self->viewportCount = viewportCount;
	
	/* Setup output format stuff */
	self->outputFormat_Register = lucOutputFormat_Register_New();
	for ( outputFormat_I = 0 ; outputFormat_I < outputFormatCount ; outputFormat_I++ )
		lucOutputFormat_Register_Add( self->outputFormat_Register, outputFormatList[ outputFormat_I ] );
		
	/* Setup window interaction stuff */
	self->windowInteraction_Register = lucWindowInteraction_Register_New();
	self->defaultWindowInteraction = lucWindowInteraction_New( "defaultWindowInteraction" );
	lucWindowInteraction_Register_Add( self->windowInteraction_Register, self->defaultWindowInteraction );
	for ( windowInteraction_I = 0 ; windowInteraction_I < windowInteractionCount ; windowInteraction_I++ )
		lucWindowInteraction_Register_Add( self->windowInteraction_Register, windowInteractionList[ windowInteraction_I ] );
	

	lucColour_FromString( &self->backgroundColour, backgroundColourName );

	self->currStereoBuffer = lucLeft;

	/* Get window to 'execute' at each 'Dump' entry point, and force the context not to re-execute itself. */
	//context->hasExecuted = True;
	EP_AppendClassHook( Context_GetEntryPoint( context, AbstractContext_EP_DumpClass ), self->_execute, self );
	
	self->isTimedOut    = isTimedOut;
	self->maxIdleTime   = maxIdleTime;
	self->idleTime = 0;
	self->startx = 0;
	self->starty = 0;
	
	self->title = Memory_Alloc_Array( char, 100, "title string" );
	strcpy(self->title, " gLucifer Interactive Output ");
}
		
void _lucWindow_Delete( void* window ) {
	lucWindow*     self      = (lucWindow*)window;

	Stg_Class_Delete( self->outputFormat_Register );
	Stg_Class_Delete( self->defaultWindowInteraction );
	Stg_Class_Delete( self->windowInteraction_Register );

	_Stg_Component_Delete( window );
}

void _lucWindow_Print( void* window, Stream* stream ) {
	lucWindow*          self        = window;
	
	Journal_Printf( stream, "lucWindow: %s\n", self->name );

	Stream_Indent( stream );

	/* Print Parent */
	_Stg_Component_Print( self, stream );


	Journal_PrintValue( stream, self->width );
	Journal_PrintValue( stream, self->height );
	Journal_PrintValue( stream, self->interactive );
	Journal_PrintValue( stream, self->continuous );
	Journal_Printf( stream, "self->currStereoBuffer = ");
	switch ( self->currStereoBuffer ) {
		case lucLeft:
			Journal_Printf( stream, "lucLeft\n" ); break;
		case lucRight:
			Journal_Printf( stream, "lucRight\n" ); break;
	}
	Stream_UnIndent( stream );
}

void* _lucWindow_Copy( void* window, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	lucWindow* self        = window;
	lucWindow* newWindow;

	newWindow = _Stg_Component_Copy( self, dest, deep, nameExt, ptrMap );

	/* TODO */
	abort();

	return (void*) newWindow;
}

void* _lucWindow_DefaultNew( Name name ) {
	return _lucWindow_New( 
			sizeof( lucWindow ),
			lucWindow_Type,
			_lucWindow_Delete,
			_lucWindow_Print,
			_lucWindow_Copy,
			_lucWindow_DefaultNew,
			_lucWindow_Construct,
			_lucWindow_Build,
			_lucWindow_Initialise,
			_lucWindow_Execute,
			_lucWindow_Destroy,
			lucWindow_Display,
			lucWindow_EventsWaiting,
			lucWindow_EventProcessor,
			lucWindow_Resize,
			name );
}

void _lucWindow_Construct( void* window, Stg_ComponentFactory* cf, void* data ) {
	lucWindow*               self        = window;
	lucViewportInfo*         viewportInfoList;
	Viewport_Index           viewportCount;
	lucOutputFormat**        outputFormatList;
	OutputFormat_Index       outputFormatCount;
	AbstractContext*         context;
	lucRenderingEngine*      renderingEngine;
	Pixel_Index              width;
	Pixel_Index              height;
	lucWindowInteraction**   windowInteractionList;
	WindowInteraction_Index  windowInteractionCount;
	Bool                     interactiveDefault;
	Bool                     interactive;
	Bool                     continuousDefault;
	Bool                     continuous;
	
	width = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, "width", 400 );
	height = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, "height", 400 );

	/* Get information about whether this window is interactive or not
	 * All lucWindow objects will check in the root dictionary first to see if interactivity in general is turned on
	 * Specific lucWindow objects can override this parameter in their own component structs */
	interactiveDefault = Stg_ComponentFactory_GetRootDictBool( cf, "interactive", False );
	interactive = Stg_ComponentFactory_GetBool( cf, self->name, "interactive", interactiveDefault );

	/* Get information about whether this window is in continuous mode or not */
	continuousDefault = Stg_ComponentFactory_GetRootDictBool( cf, "continuous", False );
	continuous = Stg_ComponentFactory_GetBool( cf, self->name, "continuous", continuousDefault );

	/* Grab information about what viewports are going to be plotting in this window */
	viewportInfoList = lucWindow_ConstructViewportInfoList( self, cf, width, height, &viewportCount, data );

	/* Grab a list of different output formats for this window to be dumped in */
	outputFormatList = Stg_ComponentFactory_ConstructByList( 
		cf, 
		self->name, 
		"OutputFormat", 
		Stg_ComponentFactory_Unlimited, 
		lucOutputFormat, 
		False, 
		&outputFormatCount, 
		data );
			
	/* Grab a list of interactions for the user to interact with this window */
	windowInteractionList = Stg_ComponentFactory_ConstructByList( 
		cf, 
		self->name, 
		"WindowInteraction", 
		Stg_ComponentFactory_Unlimited, 
		lucWindowInteraction, 
		False, 
		&windowInteractionCount,
		data );
		
	/* The window needs information about the context so that it can attach itself 
	 * onto the AbstractContext_EP_DumpClass entry point. */
	context = Stg_ComponentFactory_ConstructByName( cf, "context", AbstractContext, True, data ); 

	renderingEngine = Stg_ComponentFactory_ConstructByKey( cf, self->name, "RenderingEngine", lucRenderingEngine, True, data );

	_lucWindow_Init( 
			self,
			renderingEngine,
			viewportInfoList,
			viewportCount,
			outputFormatList,
			outputFormatCount,
			windowInteractionList,
			windowInteractionCount,
			context,
			width,
			height,
			Stg_ComponentFactory_GetString( cf, self->name, "backgroundColour", "white" ),
			interactive,
			continuous,
			Stg_ComponentFactory_GetBool( cf, self->name, "isTimedOut", True ),
			Stg_ComponentFactory_GetDouble( cf, self->name, "maxIdleTime", 600.0 ) 
			);
		
	/* Free Memory */
	Memory_Free(viewportInfoList); 
	if ( windowInteractionList )
		Memory_Free(windowInteractionList); 
	if ( outputFormatList )
		Memory_Free(outputFormatList); 

}

void _lucWindow_Build( void* window, void* data ) {
	/* Save context and master flag */
	lucWindow*  self = (lucWindow*)window;
	self->context   = (AbstractContext*) data;
	self->isMaster = (self->context->rank == MASTER);
}

void _lucWindow_Initialise( void* window, void* data ) {
	lucWindow*           self            = (lucWindow*)window;

	lucSetupRasterFont();

	/* Flag display lists must be created and objects drawn */
	lucWindow_SetViewportNeedsToSetupFlag( self, True );
	lucWindow_SetViewportNeedsToDrawFlag( self, True );
}

void _lucWindow_Execute( void* window, void* data ) { 
	/* Display graphics and allow GUI interaction if enabled */
	lucWindow*     self         = (lucWindow*) window ;
	lucDebug_PrintFunctionBegin( self, 1 );

	/* Reset idle timer */
	lucWindow_IdleReset(self);	

	/* Flag viewports need to re-render new information */
	lucWindow_SetViewportNeedsToSetupFlag( self, True );
	lucWindow_SetViewportNeedsToDrawFlag( self, True );

	/* Draw Window (Call virtual to display) 
	 initial output for background & interactive modes */
	self->_displayWindow( self );

	/* Interactive mode? Enter event loop processing */
	if ( self->interactive ) {
		Bool redisplay = False;
		lucWindow_InteractionHelpMessage( self, Journal_MyStream( Info_Type, self ) );
		
		/* Clear quit flags */
		self->quitEventLoop = False;
		self->toggleApplicationQuit = False;
		
		while ( !self->quitEventLoop ) {
            /* Only master processes events */
            int events = 0;
			if (self->isMaster)
            {  
                /* Check for events */
                events = self->_eventsWaiting(self);
                if (self->continuous && events == 0)
                {
                    /* Continuous mode and no events waiting, quit loop */
                    self->quitEventLoop = True; 
                    redisplay = True;
                }
                else	
                    /* Call virtual to wait for and process events */
                    if (self->_eventProcessor(self) || events > 1) redisplay = True;
            }

        	/* Broadcast information about event loop*/ 
        	MPI_Bcast( &events, 1, MPI_INT, MASTER, self->context->communicator );
        	MPI_Bcast( &redisplay, 1, MPI_INT, MASTER, self->context->communicator );
      
			/* Still events to process? delay redisplay until queue empty */
			if (events <= 1)
			{
				/* Redraw Window (Call virtual to display) */
				if (redisplay) self->_displayWindow( self );
				redisplay = False;
			}
		}

		/* Close application if requested */
		if( self->toggleApplicationQuit ) {
			self->context->gracefulQuit = True;
		}
	}
	
	/* Dump outputs to disk */
	if ( self->isMaster ) lucWindow_Dump( self, self->context );
	
	/* Stop idle timeout */
	self->idleTime = 0;
	
	lucDebug_PrintFunctionEnd( self, 1 );
}

void _lucWindow_Destroy( void* window, void* data ) {
	lucWindow*     self      = (lucWindow*)window;

	lucWindow_CleanUp( window, data );
	Memory_Free(self->title);
}

void lucWindow_Display( void* window )
{
	/* Default display function... to be called from derived class display function */
	lucWindow*     self      = (lucWindow*)window;

	lucWindow_Broadcast( window, 0, MPI_COMM_WORLD );

    /* Resize viewports if necessary */
    if (self->resized) self->_resizeWindow(window);
    
	lucRenderingEngine_Render( self->renderingEngine, self, self->context );
}

int lucWindow_EventsWaiting( void* window )
{
	/* Dummy events waiting function... */
	return 0;
}

Bool lucWindow_EventProcessor( void *window)
{
	/* Dummy event processor function... */
	lucWindow*     self      = (lucWindow*)window;
	self->quitEventLoop = True;
	/* Returns true redisplay required */
	return False;
}

void lucWindow_Resize( void* window ) {
    /* Default window resize function */
	lucWindow*       self      = (lucWindow*) window;
	Viewport_Index   viewport_I;
	Viewport_Index   viewportCount = self->viewportCount;
	Viewport_Index   horizontalCount = viewportCount / self->verticalCount;
	lucViewportInfo* viewportInfo;

	for ( viewport_I = 0 ; viewport_I < viewportCount ; viewport_I++ ) {
		viewportInfo = &self->viewportInfoList[ viewport_I ];
        int vertical_I = viewport_I / horizontalCount;
        int horizontal_I = viewport_I % horizontalCount;
		viewportInfo->width  = div( self->width, horizontalCount ).quot;
		viewportInfo->height = div( self->height, self->verticalCount ).quot;
		viewportInfo->startx = horizontal_I * viewportInfo->width;
		viewportInfo->starty = (self->verticalCount - 1 - vertical_I) * viewportInfo->height;
	}

    self->resized = False;
     _lucWindow_Initialise(self, self->context);	/* Reset font stuff */
}

void lucWindow_Dump( void* window, AbstractContext* context ) {
	lucWindow*     self         = (lucWindow*) window ;
	Pixel_Index    width        = self->width;
	Pixel_Index    height       = self->height;
	lucPixel*      imageBuffer  = NULL;
	Stream*        errorStream  = Journal_MyStream( Error_Type, self );

	lucDebug_PrintFunctionBegin( self, 1 );

	/* Block until all rendering complete */
	glFinish();  

	/* Allocate Memory */
	imageBuffer = Memory_Alloc_Array( lucPixel, width * height, "Pixels" );
	Journal_Firewall( imageBuffer != NULL, errorStream, "In func %s: Cannot allocate array.", __func__ );

	/* Grab Pixels from window */
	lucRenderingEngine_GetPixelData( self->renderingEngine, self, imageBuffer );

	/* Output in different formats that the user gives */
	lucOutputFormat_Register_OutputAll( self->outputFormat_Register, self, context, imageBuffer );
	
	/* Free memory */
	Memory_Free( imageBuffer );
	lucDebug_PrintFunctionEnd( self, 1 );
}

void lucWindow_CleanUp( void* window, void* context ) {
	lucWindow*     self      = (lucWindow*) window;
	Viewport_Index viewport_I;
	Viewport_Index viewportCount = self->viewportCount;
	lucViewport*   viewport;

	for ( viewport_I = 0 ; viewport_I < viewportCount ; viewport_I++ ) {
		viewport = self->viewportInfoList[ viewport_I ].viewport;
		
		lucViewport_CleanUp( viewport, context );
	}
}

/* Window event processing - mouse & keyboard */
void lucWindow_MouseMotion( void* window, lucMouseButton button, Pixel_Index xpos, Pixel_Index ypos) {
	lucWindow*            self      = (lucWindow*) window;
	Index                 windowInteraction_I;
	Index                 windowInteractionCount = lucWindowInteraction_Register_GetCount( self->windowInteraction_Register ); 
	lucWindowInteraction* windowInteraction;

	lucDebug_PrintFunctionBegin( self, 2 );

	for ( windowInteraction_I = 0 ; windowInteraction_I < windowInteractionCount ; windowInteraction_I++ ) {
		windowInteraction = lucWindowInteraction_Register_GetByIndex( self->windowInteraction_Register, windowInteraction_I );
		lucWindowInteraction_MouseMotion( windowInteraction, window, button, xpos, ypos, self->startx, self->starty );
	}
	self->startx = xpos;
	self->starty = ypos;

	lucDebug_PrintFunctionEnd( self, 2 );
}

void lucWindow_MouseClick( void* window, lucMouseButton button, lucMouseState state, Pixel_Index xpos, Pixel_Index ypos) {
	lucWindow*            self      = (lucWindow*) window;
	Index                 windowInteraction_I;
	Index                 windowInteractionCount = lucWindowInteraction_Register_GetCount( self->windowInteraction_Register ); 
	lucWindowInteraction* windowInteraction;
	
	lucDebug_PrintFunctionBegin( self, 2 );
	for ( windowInteraction_I = 0 ; windowInteraction_I < windowInteractionCount ; windowInteraction_I++ ) {
		windowInteraction = lucWindowInteraction_Register_GetByIndex( self->windowInteraction_Register, windowInteraction_I );
		lucWindowInteraction_MouseClick( windowInteraction, window, button, state, xpos, ypos );
	}
	self->startx = xpos;
	self->starty = ypos;

	lucDebug_PrintFunctionEnd( self, 2 );
}

void lucWindow_KeyboardEvent( void* window, char key, Pixel_Index xpos, Pixel_Index ypos ) {
	lucWindow*              self      = (lucWindow*) window;
	WindowInteraction_Index windowInteraction_I;
	WindowInteraction_Index windowInteractionCount = lucWindowInteraction_Register_GetCount( self->windowInteraction_Register ); 
	lucWindowInteraction*   windowInteraction;
	
	lucDebug_PrintFunctionBegin( self, 2 );

	for ( windowInteraction_I = 0 ; windowInteraction_I < windowInteractionCount; windowInteraction_I++ ) {
		windowInteraction = lucWindowInteraction_Register_GetByIndex( self->windowInteraction_Register, windowInteraction_I );
		lucWindowInteraction_KeyboardEvent( windowInteraction, window, key, xpos, ypos );
	}
	self->startx = xpos;
	self->starty = ypos;

	lucDebug_PrintFunctionEnd( self, 2 );
}

Bool lucWindow_HasStereoCamera( lucWindow* self ) {
	Viewport_Index viewport_I;
	Viewport_Index viewportCount = self->viewportCount;
	lucViewport*   viewport;

	for ( viewport_I = 0 ; viewport_I < viewportCount ; viewport_I++ ) {
		viewport = self->viewportInfoList[ viewport_I ].viewport;
		
		if (viewport->camera->stereoType != lucMono )
			return True;
	}
	return False;
}

void lucWindow_CheckCameraFlag( void* window ) {
	lucWindow*            self          = (lucWindow*) window;
	Viewport_Index        viewport_I;
	Viewport_Index        viewportCount = self->viewportCount;
	lucViewport*          viewport;
	lucViewportInfo*      viewportInfo;
	
	lucDebug_PrintFunctionBegin( self, 2 );

	for ( viewport_I = 0 ; viewport_I < viewportCount ; viewport_I++ ) {
		viewportInfo = &self->viewportInfoList[ viewport_I ];
		viewport = viewportInfo->viewport;

		if ( viewport->camera->needsToDraw )
			viewportInfo->needsToDraw = True;

	}

	/* Now that camera has passed message to viewport - the camera's flag is no longer nessesary */
	for ( viewport_I = 0 ; viewport_I < viewportCount ; viewport_I++ ) {
		viewport = self->viewportInfoList[ viewport_I ].viewport;

		viewport->camera->needsToDraw = False;
	}
	
	lucDebug_PrintFunctionEnd( self, 2 );
}


void lucWindow_CheckLightFlag( void* window ) {
	lucWindow*            self          = (lucWindow*) window;
	Viewport_Index        viewport_I;
	Viewport_Index        viewportCount = self->viewportCount;
	lucViewport*          viewport;
	lucViewportInfo*      viewportInfo;
	
	lucLight*                   light;
	Light_Index                 light_I = 0;
	Light_Index                 lightCount = 0;

	lucDebug_PrintFunctionBegin( self, 2 );
	
	/* Loop through lights that are registered on the window */
	
	for ( viewport_I = 0 ; viewport_I < viewportCount ; viewport_I++ ) {
	
		viewportInfo = &self->viewportInfoList[ viewport_I ];
		viewport = viewportInfo->viewport;

		/* Finds the number of lights, and if one has been modified, then redraw */
		lightCount   = lucLight_Register_GetCount( viewport->light_Register );
		for ( light_I = 0 ; light_I < lightCount ; light_I++ ) {
			light = lucLight_Register_GetByIndex( viewport->light_Register, light_I );

			if ( light->needsToDraw )
				viewportInfo->needsToDraw = True;

		}
	}

	/* Now that lights has passed message to viewport - the lights' flag is no longer nessesary */
	for ( viewport_I = 0 ; viewport_I < viewportCount ; viewport_I++ ) {
		viewport = self->viewportInfoList[ viewport_I ].viewport;

		lightCount   = lucLight_Register_GetCount( viewport->light_Register );
		for ( light_I = 0 ; light_I < lightCount ; light_I++ ) {
			light = lucLight_Register_GetByIndex( viewport->light_Register, light_I );						 light->needsToDraw = False;
		}
	}
	
	lucDebug_PrintFunctionEnd( self, 2 );
}

void lucWindow_SetViewportNeedsToSetupFlag( void* window, Bool flag ) {
	lucWindow*              self         = (lucWindow*) window ;
	Viewport_Index          viewport_I;
	
	lucDebug_PrintFunctionBegin( self, 2 );
	
	for ( viewport_I = 0 ; viewport_I < self->viewportCount ; viewport_I++ ) {
		lucViewport_SetNeedsToSetupFlag( self->viewportInfoList[ viewport_I ].viewport, flag );
	}
	lucDebug_PrintFunctionEnd( self, 2 );
}

void lucWindow_SetViewportNeedsToDrawFlag( void* window, Bool flag ) {
	lucWindow*              self         = (lucWindow*) window ;
	Viewport_Index          viewport_I;

	lucDebug_PrintFunctionBegin( self, 2 );
	for ( viewport_I = 0 ; viewport_I < self->viewportCount ; viewport_I++ ) {
		self->viewportInfoList[ viewport_I ].needsToDraw = flag;
	}
	lucDebug_PrintFunctionEnd( self, 2 );
}

lucViewportInfo* lucWindow_GetViewportInfoByPixel( void* window, Pixel_Index xPixel, Pixel_Index yPixel ) {
	lucWindow*              self         = (lucWindow*) window ;
	Viewport_Index          viewport_I;
	
	for ( viewport_I = 0 ; viewport_I < self->viewportCount ; viewport_I++ ) {
		if ( lucViewportInfo_IsPixelInside( &self->viewportInfoList[ viewport_I ], xPixel, yPixel ) )
			return &self->viewportInfoList[ viewport_I ];
	}

	return NULL;	
}

lucViewport* lucWindow_GetViewportByPixel( void* window, Pixel_Index xPixel, Pixel_Index yPixel ) {
	lucWindow*              self         = (lucWindow*) window ;
	lucViewportInfo*        viewportInfo;

	viewportInfo = lucWindow_GetViewportInfoByPixel( self, xPixel, yPixel );
	
	if ( viewportInfo )
		return viewportInfo->viewport;
	else 
		return NULL;
}


void lucWindow_Broadcast( void* window, int rootRank, MPI_Comm comm ) {
	lucWindow*         self      = (lucWindow*) window;
	Viewport_Index     viewport_I;
	
	lucDebug_PrintFunctionBegin( self, 2 );

	MPI_Bcast( self, 1, lucWindow_MPI_Datatype, rootRank, comm );

	/* Broadcast Viewport Infos */
	for ( viewport_I = 0 ; viewport_I < self->viewportCount ; viewport_I++ ) {
		lucViewportInfo_Broadcast( &self->viewportInfoList[ viewport_I ], rootRank, comm );
	}
	lucDebug_PrintFunctionEnd( self, 2 );
}

void lucViewportInfo_Broadcast( lucViewportInfo* self, int rootRank, MPI_Comm comm ) {
	Journal_DPrintfL( lucDebug, 2, "In %s.\n", __func__ );

	lucViewport_Broadcast( self->viewport, rootRank, comm );

	MPI_Bcast( self, 1, lucViewportInfo_MPI_Datatype, rootRank, comm );
}

Bool lucViewportInfo_IsPixelInside( lucViewportInfo* viewportInfo, Pixel_Index xPixel, Pixel_Index yPixel ) {
	if ( xPixel < viewportInfo->startx ) 
		return False;
	if ( yPixel < viewportInfo->starty ) 
		return False;
	
	if ( xPixel >= viewportInfo->startx + viewportInfo->width ) 
		return False;
	if ( yPixel >= viewportInfo->starty + viewportInfo->height ) 
		return False;

	return True;
}	

lucViewportInfo* lucWindow_ConstructViewportInfoList( 
		lucWindow* self, 
		Stg_ComponentFactory* cf, 
		Pixel_Index width, 
		Pixel_Index height, 
		Viewport_Index* viewportCount,
		void* data ) 
{
	Viewport_Index          vertical_I;
	Viewport_Index          horizontal_I;
	Viewport_Index          horizontalCount;
	Viewport_Index          total_I                 = 0;
	Pixel_Index             viewportHeight;
	Pixel_Index             viewportWidth;
	Dictionary_Entry_Value* list;
	char*                   charPtr;
	char*                   horizontalVP_String;
	const char*             breakChars              = "\n\t ,;";
	lucViewportInfo*        viewportInfoList        = Memory_Alloc_Array( lucViewportInfo, 1, "viewportInfoArray" );
	lucViewportInfo*        currViewportInfo;
	Dictionary*             dictionary              = Dictionary_GetDictionary( cf->componentDict, self->name );
	Name                    viewportName;
	char*                   bufferPtr;

	*viewportCount = 0;
	
	list = Dictionary_Get( dictionary, "Viewport" );
	Journal_Firewall( list != NULL, lucError, "Cannot Find 'Viewport' for %s '%s'.\n", self->type, self->name );
		
	self->verticalCount = Dictionary_Entry_Value_GetCount( list );

	/* Calculate viewport height */
	viewportHeight = div( height, self->verticalCount ).quot;

	for ( vertical_I = 0 ; vertical_I < self->verticalCount ; vertical_I++ ) {
		horizontalVP_String = StG_Strdup( Dictionary_Entry_Value_AsString( Dictionary_Entry_Value_GetElement( list, vertical_I ) ) );
	
		/* Find number of horizontal layers */
		horizontalCount = 1;
		charPtr = strpbrk( horizontalVP_String, breakChars );
		while ( charPtr != NULL ) {
			charPtr = strpbrk( charPtr + 1 , breakChars );
			horizontalCount++;
		}
		
		/* Sum up total number of viewports */
		*viewportCount += horizontalCount;

		/* Reallocate */
		viewportInfoList = Memory_Realloc_Array( viewportInfoList, lucViewportInfo, *viewportCount );
		
		/* Calculate viewport width */
		viewportWidth = div( width, horizontalCount ).quot;

		/* Read String to get colour map */
		charPtr = strtok_r( horizontalVP_String, breakChars, &bufferPtr );
		for ( horizontal_I = 0 ; horizontal_I < horizontalCount ; horizontal_I++ ) {
			currViewportInfo = &viewportInfoList[ total_I ];

			/* Find viewport */
			viewportName = StG_Strdup( charPtr );
			currViewportInfo->viewport = Stg_ComponentFactory_ConstructByName( cf, viewportName, lucViewport, True, data ) ;
			Memory_Free( viewportName );

			/* Setup viewport dimensions */
			currViewportInfo->startx = horizontal_I * viewportWidth;
			currViewportInfo->starty = (self->verticalCount - 1 - vertical_I) * viewportHeight;
			currViewportInfo->width  = viewportWidth;
			currViewportInfo->height = viewportHeight;

			charPtr = strtok_r( NULL, breakChars, &bufferPtr );

			total_I++;
		}

		Memory_Free( horizontalVP_String );
	}	

	Journal_Firewall( total_I == *viewportCount, lucError, 
			"Something went wrong in %s for %s '%s' - Incorrectly counted number of viewports.\n", 
			__func__, self->type, self->name );

	return viewportInfoList;
}

void lucWindow_InteractionHelpMessage( void* window, Stream* stream ) {
	lucWindow*              self      = (lucWindow*) window;
	WindowInteraction_Index windowInteraction_I;
	WindowInteraction_Index windowInteractionCount = lucWindowInteraction_Register_GetCount( self->windowInteraction_Register ); 
	lucWindowInteraction*   windowInteraction;

	Journal_Printf( stream, "--------------------------------------------------------------------\n" );
	Journal_Printf( stream, "Help for %s '%s':\n", self->type, self->name );

	Journal_Printf( stream, "Window Interactions registered on this window are:\n" );
	Stream_Indent( stream );
	for ( windowInteraction_I = 0 ; windowInteraction_I < windowInteractionCount ; windowInteraction_I++ ) {
		windowInteraction = lucWindowInteraction_Register_GetByIndex( self->windowInteraction_Register, windowInteraction_I );
		Journal_Printf( stream, "'%s' (type %s)\n", windowInteraction->name, windowInteraction->type );
	}
	Stream_UnIndent( stream );
	
	Journal_Printf( stream, "Mouse Interaction:\n" );
	Stream_Indent( stream );
	for ( windowInteraction_I = 0 ; windowInteraction_I < windowInteractionCount ; windowInteraction_I++ ) {
		windowInteraction = lucWindowInteraction_Register_GetByIndex( self->windowInteraction_Register, windowInteraction_I );
		lucWindowInteraction_MouseMessage( windowInteraction, stream );
	}
	Stream_UnIndent( stream );	

	Journal_Printf( stream, "Keyboard Interaction:\n" );
	Stream_Indent( stream );
	for ( windowInteraction_I = 0 ; windowInteraction_I < windowInteractionCount ; windowInteraction_I++ ) {
		windowInteraction = lucWindowInteraction_Register_GetByIndex( self->windowInteraction_Register, windowInteraction_I );
		lucWindowInteraction_KeyboardMessage( windowInteraction, stream );
	}
	Stream_UnIndent( stream );	
	Journal_Printf( stream, "--------------------------------------------------------------------\n" );
}
	

#define lucWindow_TypesCount 7
void lucWindow_Create_MPI_Datatype() {
	MPI_Datatype        typeList[lucWindow_TypesCount]     = { MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_FLOAT, MPI_INT };
	int                 blocklen[lucWindow_TypesCount]     = {1, 1, 1, 1, 1, 4, 1};
	MPI_Aint            displacement[lucWindow_TypesCount];
	lucWindow           window;

	displacement[0] = GetOffsetOfMember( window, width );
	displacement[1] = GetOffsetOfMember( window, height );
	displacement[2] = GetOffsetOfMember( window, interactive );
	displacement[3] = GetOffsetOfMember( window, quitEventLoop );
	displacement[4] = GetOffsetOfMember( window, resized );
	displacement[5] = GetOffsetOfMember( window, backgroundColour );
	displacement[6] = GetOffsetOfMember( window, currStereoBuffer );
	
	MPI_Type_struct( lucWindow_TypesCount, blocklen, displacement, typeList, &lucWindow_MPI_Datatype );
	MPI_Type_commit( & lucWindow_MPI_Datatype );
}


#define lucViewportInfo_TypesCount 5
void lucViewportInfo_Create_MPI_Datatype() {
	MPI_Datatype        typeList[lucViewportInfo_TypesCount]     = { MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT };
	int                 blocklen[lucViewportInfo_TypesCount]     = {1, 1, 1, 1, 1};
	MPI_Aint            displacement[lucViewportInfo_TypesCount];
	lucViewportInfo     viewportInfo;

	displacement[0] = GetOffsetOfMember( viewportInfo, startx );
	displacement[1] = GetOffsetOfMember( viewportInfo, starty );
	displacement[2] = GetOffsetOfMember( viewportInfo, width );
	displacement[3] = GetOffsetOfMember( viewportInfo, height );
	displacement[4] = GetOffsetOfMember( viewportInfo, needsToDraw );
	
	MPI_Type_struct( lucViewportInfo_TypesCount, blocklen, displacement, typeList, &lucViewportInfo_MPI_Datatype );
	MPI_Type_commit( & lucViewportInfo_MPI_Datatype );
}

void lucWindow_ChangeInteractiveMode( void* window ) {
	lucWindow*              self         = (lucWindow*) window ;
	self->interactive = !self->interactive;
}


void lucWindow_ChangeContinuousMode( void* window ) {
	lucWindow*              self         = (lucWindow*) window ;
	self->continuous = !self->continuous;
}

void lucWindow_QuitEventLoop( void* window ) {
	lucWindow*              self         = (lucWindow*) window ;

	self->quitEventLoop = True;
}

void lucWindow_ToggleApplicationQuit( void* window ) {
	lucWindow*              self         = (lucWindow*) window ;

	lucWindow_QuitEventLoop( self );
	self->toggleApplicationQuit = True;
}

Bool lucWindow_SetSize( void* window, Pixel_Index newWidth, Pixel_Index newHeight ) {
	lucWindow*       self      = (lucWindow*) window;

	if (newWidth < 32 || newHeight < 32) return False;
    if (newWidth == self->width && newHeight == self->height) return False; /* No resize neccessary */

	self->width = newWidth;
	self->height = newHeight;

    self->resized = True;
  	Journal_DPrintfL( lucDebug, 2, "Window resized %d x %d\n", self->width, self->height);

    return True;
}

void lucWindow_IdleReset(void *window) {
	/* Update Idle timer */
	lucWindow* self = (lucWindow*)window;
	self->idleTime = MPI_Wtime();
}

void lucWindow_IdleCheck(void *window) {
	lucWindow* self = (lucWindow*)window;
	if (self->isTimedOut && self->idleTime > 0){
		if ( MPI_Wtime() > self->idleTime + self->maxIdleTime ) {
			Journal_Printf( lucError, "Error in func '%s' - Interactive window '%s' open for too long (over %g seconds) without interaction.\n",
							__func__, self->name, self->maxIdleTime );
			abort();
		}
	}
}


