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
** $Id: Window.c 658 2007-02-05 00:58:07Z CecileDuboz $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>

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

const Type lucWindow_Type = "lucWindow";

MPI_Datatype lucWindow_MPI_Datatype;
MPI_Datatype lucViewportInfo_MPI_Datatype;

lucWindow* _lucWindow_New(
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
		Bool                                               interactive )
{
	OutputFormat_Index   outputFormat_I;
	WindowInteraction_Index windowInteraction_I;

	self->renderingEngine = renderingEngine;
	self->width = width;
	self->height = height;
	self->interactive = interactive;

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
}
		
void _lucWindow_Delete( void* window ) {
	lucWindow* self        = window;

	Stg_Class_Delete( self->outputFormat_Register );
	Stg_Class_Delete( self->defaultWindowInteraction );
	Stg_Class_Delete( self->windowInteraction_Register );
		
	_Stg_Component_Delete( self );
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
	
	width = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, "width", 400 );
	height = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, "height", 400 );

	/* Get information about whether this window is interactive or not
	 * All lucWindow objects will check in the root dictionary first to see if interactivity in general is turned on
	 * Specific lucWindow objects can override this parameter in their own component structs */
	interactiveDefault = Stg_ComponentFactory_GetRootDictBool( cf, "interactive", False );
	interactive = Stg_ComponentFactory_GetBool( cf, self->name, "interactive", interactiveDefault );

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
		
	/* The window needs information about the context so that it can attach itself onto the AbstractContext_EP_DumpClass entry 
	 * point. */
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
			interactive );
		
	/* Free Memory */
	Memory_Free(viewportInfoList); 
	if ( windowInteractionList )
		Memory_Free(windowInteractionList); 
	if ( outputFormatList )
		Memory_Free(outputFormatList); 

}

void _lucWindow_Build( void* window, void* data ) { }
void _lucWindow_Initialise( void* window, void* data ) { }
void _lucWindow_Execute( void* window, void* data ) { 
	lucWindow*     self         = (lucWindow*) window ;

	lucWindow_SetViewportNeedsToSetupFlag( self, True );
	lucWindow_SetViewportNeedsToDrawFlag( self, True );

	lucWindow_Draw( self, data );
	lucWindow_Dump( self, data );
	lucWindow_CleanUp( self, data );
}
void _lucWindow_Destroy( void* window, void* data ) { }

void lucWindow_Draw( void* window, AbstractContext* context ) {
	lucWindow*     self         = (lucWindow*) window ;

	lucRenderingEngine_Render( self->renderingEngine, self, context );
}

void lucWindow_Dump( void* window, AbstractContext* context ) {
	lucWindow*     self         = (lucWindow*) window ;
	Pixel_Index    width        = self->width;
	Pixel_Index    height       = self->height;
	lucPixel*      imageBuffer  = NULL;
	Stream*        errorStream  = Journal_MyStream( Error_Type, self );
	
	lucDebug_PrintFunctionBegin( self, 1 );

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

void lucWindow_MouseMotion( void* window, lucMouseButton button, Pixel_Index xpos, Pixel_Index ypos, Pixel_Index startx, Pixel_Index starty) {
	lucWindow*            self      = (lucWindow*) window;
	Index                 windowInteraction_I;
	Index                 windowInteractionCount = lucWindowInteraction_Register_GetCount( self->windowInteraction_Register ); 
	lucWindowInteraction* windowInteraction;
	
	lucDebug_PrintFunctionBegin( self, 2 );

	for ( windowInteraction_I = 0 ; windowInteraction_I < windowInteractionCount ; windowInteraction_I++ ) {
		windowInteraction = lucWindowInteraction_Register_GetByIndex( self->windowInteraction_Register, windowInteraction_I );
		lucWindowInteraction_MouseMotion( windowInteraction, window, button, xpos, ypos, startx, starty );
	}

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
	Viewport_Index          verticalCount;
	Viewport_Index          horizontalCount;
	Viewport_Index          vertical_I;
	Viewport_Index          horizontal_I;
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
		
	verticalCount = Dictionary_Entry_Value_GetCount( list );

	/* Calculate viewport height */
	viewportHeight = div( height, verticalCount ).quot;

	for ( vertical_I = 0 ; vertical_I < verticalCount ; vertical_I++ ) {
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
			currViewportInfo->starty = (verticalCount - 1 - vertical_I) * viewportHeight;
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
	

#define lucWindow_TypesCount 5
void lucWindow_Create_MPI_Datatype() {
	MPI_Datatype        typeList[lucWindow_TypesCount]     = { MPI_INT, MPI_INT, MPI_INT, MPI_FLOAT, MPI_INT };
	int                 blocklen[lucWindow_TypesCount]     = {1, 1, 1, 4, 1};
	MPI_Aint            displacement[lucWindow_TypesCount];
	lucWindow           window;

	displacement[0] = GetOffsetOfMember( window, width );
	displacement[1] = GetOffsetOfMember( window, height );
	displacement[2] = GetOffsetOfMember( window, interactive );
	displacement[3] = GetOffsetOfMember( window, backgroundColour );
	displacement[4] = GetOffsetOfMember( window, currStereoBuffer );
	
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

#ifdef HAVE_GL

#include <gl.h>

void _lucWindow_SetupGLRasterFont( void* window ) {
	GLuint i, j;

	GLubyte rasters[][13] = {
	{0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00}, 
	{0x00, 0x00, 0x18, 0x18, 0x00, 0x00, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18}, 
	{0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x36, 0x36, 0x36, 0x36}, 
	{0x00, 0x00, 0x00, 0x66, 0x66, 0xff, 0x66, 0x66, 0xff, 0x66, 0x66, 0x00, 0x00}, 
	{0x00, 0x00, 0x18, 0x7e, 0xff, 0x1b, 0x1f, 0x7e, 0xf8, 0xd8, 0xff, 0x7e, 0x18}, 
	{0x00, 0x00, 0x0e, 0x1b, 0xdb, 0x6e, 0x30, 0x18, 0x0c, 0x76, 0xdb, 0xd8, 0x70}, 
	{0x00, 0x00, 0x7f, 0xc6, 0xcf, 0xd8, 0x70, 0x70, 0xd8, 0xcc, 0xcc, 0x6c, 0x38}, 
	{0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x18, 0x1c, 0x0c, 0x0e}, 
	{0x00, 0x00, 0x0c, 0x18, 0x30, 0x30, 0x30, 0x30, 0x30, 0x30, 0x30, 0x18, 0x0c}, 
	{0x00, 0x00, 0x30, 0x18, 0x0c, 0x0c, 0x0c, 0x0c, 0x0c, 0x0c, 0x0c, 0x18, 0x30}, 
	{0x00, 0x00, 0x00, 0x00, 0x99, 0x5a, 0x3c, 0xff, 0x3c, 0x5a, 0x99, 0x00, 0x00}, 
	{0x00, 0x00, 0x00, 0x18, 0x18, 0x18, 0xff, 0xff, 0x18, 0x18, 0x18, 0x00, 0x00}, 
	{0x00, 0x00, 0x30, 0x18, 0x1c, 0x1c, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00}, 
	{0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0xff, 0x00, 0x00, 0x00, 0x00, 0x00}, 
	{0x00, 0x00, 0x00, 0x38, 0x38, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00}, 
	{0x00, 0x60, 0x60, 0x30, 0x30, 0x18, 0x18, 0x0c, 0x0c, 0x06, 0x06, 0x03, 0x03}, 
	{0x00, 0x00, 0x3c, 0x66, 0xc3, 0xe3, 0xf3, 0xdb, 0xcf, 0xc7, 0xc3, 0x66, 0x3c}, 
	{0x00, 0x00, 0x7e, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x78, 0x38, 0x18}, 
	{0x00, 0x00, 0xff, 0xc0, 0xc0, 0x60, 0x30, 0x18, 0x0c, 0x06, 0x03, 0xe7, 0x7e}, 
	{0x00, 0x00, 0x7e, 0xe7, 0x03, 0x03, 0x07, 0x7e, 0x07, 0x03, 0x03, 0xe7, 0x7e}, 
	{0x00, 0x00, 0x0c, 0x0c, 0x0c, 0x0c, 0x0c, 0xff, 0xcc, 0x6c, 0x3c, 0x1c, 0x0c}, 
	{0x00, 0x00, 0x7e, 0xe7, 0x03, 0x03, 0x07, 0xfe, 0xc0, 0xc0, 0xc0, 0xc0, 0xff}, 
	{0x00, 0x00, 0x7e, 0xe7, 0xc3, 0xc3, 0xc7, 0xfe, 0xc0, 0xc0, 0xc0, 0xe7, 0x7e}, 
	{0x00, 0x00, 0x30, 0x30, 0x30, 0x30, 0x18, 0x0c, 0x06, 0x03, 0x03, 0x03, 0xff}, 
	{0x00, 0x00, 0x7e, 0xe7, 0xc3, 0xc3, 0xe7, 0x7e, 0xe7, 0xc3, 0xc3, 0xe7, 0x7e}, 
	{0x00, 0x00, 0x7e, 0xe7, 0x03, 0x03, 0x03, 0x7f, 0xe7, 0xc3, 0xc3, 0xe7, 0x7e}, 
	{0x00, 0x00, 0x00, 0x38, 0x38, 0x00, 0x00, 0x38, 0x38, 0x00, 0x00, 0x00, 0x00}, 
	{0x00, 0x00, 0x30, 0x18, 0x1c, 0x1c, 0x00, 0x00, 0x1c, 0x1c, 0x00, 0x00, 0x00}, 
	{0x00, 0x00, 0x06, 0x0c, 0x18, 0x30, 0x60, 0xc0, 0x60, 0x30, 0x18, 0x0c, 0x06}, 
	{0x00, 0x00, 0x00, 0x00, 0xff, 0xff, 0x00, 0xff, 0xff, 0x00, 0x00, 0x00, 0x00}, 
	{0x00, 0x00, 0x60, 0x30, 0x18, 0x0c, 0x06, 0x03, 0x06, 0x0c, 0x18, 0x30, 0x60}, 
	{0x00, 0x00, 0x18, 0x00, 0x00, 0x18, 0x18, 0x0c, 0x06, 0x03, 0xc3, 0xc3, 0x7e}, 
	{0x00, 0x00, 0x3f, 0x60, 0xcf, 0xdb, 0xd3, 0xdd, 0xc3, 0x7e, 0x00, 0x00, 0x00}, 
	{0x00, 0x00, 0xc3, 0xc3, 0xc3, 0xc3, 0xff, 0xc3, 0xc3, 0xc3, 0x66, 0x3c, 0x18}, 
	{0x00, 0x00, 0xfe, 0xc7, 0xc3, 0xc3, 0xc7, 0xfe, 0xc7, 0xc3, 0xc3, 0xc7, 0xfe}, 
	{0x00, 0x00, 0x7e, 0xe7, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xe7, 0x7e}, 
	{0x00, 0x00, 0xfc, 0xce, 0xc7, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc7, 0xce, 0xfc}, 
	{0x00, 0x00, 0xff, 0xc0, 0xc0, 0xc0, 0xc0, 0xfc, 0xc0, 0xc0, 0xc0, 0xc0, 0xff}, 
	{0x00, 0x00, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xfc, 0xc0, 0xc0, 0xc0, 0xff}, 
	{0x00, 0x00, 0x7e, 0xe7, 0xc3, 0xc3, 0xcf, 0xc0, 0xc0, 0xc0, 0xc0, 0xe7, 0x7e}, 
	{0x00, 0x00, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xff, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3}, 
	{0x00, 0x00, 0x7e, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x7e}, 
	{0x00, 0x00, 0x7c, 0xee, 0xc6, 0x06, 0x06, 0x06, 0x06, 0x06, 0x06, 0x06, 0x06}, 
	{0x00, 0x00, 0xc3, 0xc6, 0xcc, 0xd8, 0xf0, 0xe0, 0xf0, 0xd8, 0xcc, 0xc6, 0xc3}, 
	{0x00, 0x00, 0xff, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0}, 
	{0x00, 0x00, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xdb, 0xff, 0xff, 0xe7, 0xc3}, 
	{0x00, 0x00, 0xc7, 0xc7, 0xcf, 0xcf, 0xdf, 0xdb, 0xfb, 0xf3, 0xf3, 0xe3, 0xe3}, 
	{0x00, 0x00, 0x7e, 0xe7, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xe7, 0x7e}, 
	{0x00, 0x00, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xfe, 0xc7, 0xc3, 0xc3, 0xc7, 0xfe}, 
	{0x00, 0x00, 0x3f, 0x6e, 0xdf, 0xdb, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0x66, 0x3c}, 
	{0x00, 0x00, 0xc3, 0xc6, 0xcc, 0xd8, 0xf0, 0xfe, 0xc7, 0xc3, 0xc3, 0xc7, 0xfe}, 
	{0x00, 0x00, 0x7e, 0xe7, 0x03, 0x03, 0x07, 0x7e, 0xe0, 0xc0, 0xc0, 0xe7, 0x7e}, 
	{0x00, 0x00, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0xff}, 
	{0x00, 0x00, 0x7e, 0xe7, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3}, 
	{0x00, 0x00, 0x18, 0x3c, 0x3c, 0x66, 0x66, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3}, 
	{0x00, 0x00, 0xc3, 0xe7, 0xff, 0xff, 0xdb, 0xdb, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3}, 
	{0x00, 0x00, 0xc3, 0x66, 0x66, 0x3c, 0x3c, 0x18, 0x3c, 0x3c, 0x66, 0x66, 0xc3}, 
	{0x00, 0x00, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x3c, 0x3c, 0x66, 0x66, 0xc3}, 
	{0x00, 0x00, 0xff, 0xc0, 0xc0, 0x60, 0x30, 0x7e, 0x0c, 0x06, 0x03, 0x03, 0xff}, 
	{0x00, 0x00, 0x3c, 0x30, 0x30, 0x30, 0x30, 0x30, 0x30, 0x30, 0x30, 0x30, 0x3c}, 
	{0x00, 0x03, 0x03, 0x06, 0x06, 0x0c, 0x0c, 0x18, 0x18, 0x30, 0x30, 0x60, 0x60}, 
	{0x00, 0x00, 0x3c, 0x0c, 0x0c, 0x0c, 0x0c, 0x0c, 0x0c, 0x0c, 0x0c, 0x0c, 0x3c}, 
	{0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0xc3, 0x66, 0x3c, 0x18}, 
	{0xff, 0xff, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00}, 
	{0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x18, 0x38, 0x30, 0x70}, 
	{0x00, 0x00, 0x7f, 0xc3, 0xc3, 0x7f, 0x03, 0xc3, 0x7e, 0x00, 0x00, 0x00, 0x00}, 
	{0x00, 0x00, 0xfe, 0xc3, 0xc3, 0xc3, 0xc3, 0xfe, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0}, 
	{0x00, 0x00, 0x7e, 0xc3, 0xc0, 0xc0, 0xc0, 0xc3, 0x7e, 0x00, 0x00, 0x00, 0x00}, 
	{0x00, 0x00, 0x7f, 0xc3, 0xc3, 0xc3, 0xc3, 0x7f, 0x03, 0x03, 0x03, 0x03, 0x03}, 
	{0x00, 0x00, 0x7f, 0xc0, 0xc0, 0xfe, 0xc3, 0xc3, 0x7e, 0x00, 0x00, 0x00, 0x00}, 
	{0x00, 0x00, 0x30, 0x30, 0x30, 0x30, 0x30, 0xfc, 0x30, 0x30, 0x30, 0x33, 0x1e}, 
	{0x7e, 0xc3, 0x03, 0x03, 0x7f, 0xc3, 0xc3, 0xc3, 0x7e, 0x00, 0x00, 0x00, 0x00}, 
	{0x00, 0x00, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xc3, 0xfe, 0xc0, 0xc0, 0xc0, 0xc0}, 
	{0x00, 0x00, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x00, 0x00, 0x18, 0x00}, 
	{0x38, 0x6c, 0x0c, 0x0c, 0x0c, 0x0c, 0x0c, 0x0c, 0x0c, 0x00, 0x00, 0x0c, 0x00}, 
	{0x00, 0x00, 0xc6, 0xcc, 0xf8, 0xf0, 0xd8, 0xcc, 0xc6, 0xc0, 0xc0, 0xc0, 0xc0}, 
	{0x00, 0x00, 0x7e, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x78}, 
	{0x00, 0x00, 0xdb, 0xdb, 0xdb, 0xdb, 0xdb, 0xdb, 0xfe, 0x00, 0x00, 0x00, 0x00}, 
	{0x00, 0x00, 0xc6, 0xc6, 0xc6, 0xc6, 0xc6, 0xc6, 0xfc, 0x00, 0x00, 0x00, 0x00}, 
	{0x00, 0x00, 0x7c, 0xc6, 0xc6, 0xc6, 0xc6, 0xc6, 0x7c, 0x00, 0x00, 0x00, 0x00}, 
	{0xc0, 0xc0, 0xc0, 0xfe, 0xc3, 0xc3, 0xc3, 0xc3, 0xfe, 0x00, 0x00, 0x00, 0x00}, 
	{0x03, 0x03, 0x03, 0x7f, 0xc3, 0xc3, 0xc3, 0xc3, 0x7f, 0x00, 0x00, 0x00, 0x00}, 
	{0x00, 0x00, 0xc0, 0xc0, 0xc0, 0xc0, 0xc0, 0xe0, 0xfe, 0x00, 0x00, 0x00, 0x00}, 
	{0x00, 0x00, 0xfe, 0x03, 0x03, 0x7e, 0xc0, 0xc0, 0x7f, 0x00, 0x00, 0x00, 0x00}, 
	{0x00, 0x00, 0x1c, 0x36, 0x30, 0x30, 0x30, 0x30, 0xfc, 0x30, 0x30, 0x30, 0x00}, 
	{0x00, 0x00, 0x7e, 0xc6, 0xc6, 0xc6, 0xc6, 0xc6, 0xc6, 0x00, 0x00, 0x00, 0x00}, 
	{0x00, 0x00, 0x18, 0x3c, 0x3c, 0x66, 0x66, 0xc3, 0xc3, 0x00, 0x00, 0x00, 0x00}, 
	{0x00, 0x00, 0xc3, 0xe7, 0xff, 0xdb, 0xc3, 0xc3, 0xc3, 0x00, 0x00, 0x00, 0x00}, 
	{0x00, 0x00, 0xc3, 0x66, 0x3c, 0x18, 0x3c, 0x66, 0xc3, 0x00, 0x00, 0x00, 0x00}, 
	{0xc0, 0x60, 0x60, 0x30, 0x18, 0x3c, 0x66, 0x66, 0xc3, 0x00, 0x00, 0x00, 0x00}, 
	{0x00, 0x00, 0xff, 0x60, 0x30, 0x18, 0x0c, 0x06, 0xff, 0x00, 0x00, 0x00, 0x00}, 
	{0x00, 0x00, 0x0f, 0x18, 0x18, 0x18, 0x38, 0xf0, 0x38, 0x18, 0x18, 0x18, 0x0f}, 
	{0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18, 0x18}, 
	{0x00, 0x00, 0xf0, 0x18, 0x18, 0x18, 0x1c, 0x0f, 0x1c, 0x18, 0x18, 0x18, 0xf0}, 
	{0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x06, 0x8f, 0xf1, 0x60, 0x00, 0x00, 0x00} 
	};

	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

	for (i = 0,j = ' ' ; j < '~' ; i++,j++) {
		glNewList( 2000 + j, GL_COMPILE);
		glBitmap(8, 13, 0.0, 2.0, 10.0, 0.0, rasters[i]);
		glEndList();
	}
	glListBase(2000);

}	

void lucWindow_ChangeInteractiveMode( void* window ) {
	lucWindow*              self         = (lucWindow*) window ;
	self->interactive = !self->interactive;
}


void lucWindow_BeginEventLoop( void* window ) {
	lucWindow*              self         = (lucWindow*) window ;

	self->quitEventLoop = False;
	self->toggleApplicationQuit = False;
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

#endif /* HAVE_GL */
