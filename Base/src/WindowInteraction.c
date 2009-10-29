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

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>

#include "types.h"
#include "ColourMap.h"
#include "X11Colours.h"
#include "ViewportInfo.h"
#include "WindowInteraction.h"
#include "WindowInteraction_Register.h"
#include "Viewport.h"
#include "Camera.h"
#include "OutputFormat.h"
#include "OutputFormat_Register.h"
#include "Init.h"
#include "RenderingEngine.h"
#include "Window.h"


#include <ctype.h>
#include <assert.h>
#include <string.h>

/* ASCII Characters */
#define ESCAPE 27

const Type lucWindowInteraction_Type = "lucWindowInteraction";

lucWindowInteraction* lucWindowInteraction_New( Name name ) {
	lucWindowInteraction* self = (lucWindowInteraction*) _lucWindowInteraction_DefaultNew( name );

	lucWindowInteraction_InitAll( self );

	return self;
}

lucWindowInteraction* _lucWindowInteraction_New(
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
		lucWindowInteraction_MouseMotionFunction*          _mouseMotion,
		lucWindowInteraction_MouseClickFunction*           _mouseClick,
		lucWindowInteraction_MouseMessageFunction*         _mouseMessage,
		lucWindowInteraction_KeyboardEventFunction*        _keyboardEvent,
		lucWindowInteraction_KeyboardMessageFunction*      _keyboardMessage,
		Name                                               name )
{
	lucWindowInteraction*    self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( sizeOfSelf >= sizeof(lucWindowInteraction) );
	self = (lucWindowInteraction*) _Stg_Component_New( 
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

	/* Set virtual functions specific to this sub-class here: */
	self->_mouseMotion     = _mouseMotion;
	self->_mouseClick      = _mouseClick;
	self->_mouseMessage    = _mouseMessage;
	self->_keyboardEvent   = _keyboardEvent;
	self->_keyboardMessage = _keyboardMessage;

	return self;
}

void _lucWindowInteraction_Init( lucWindowInteraction*                              self ) {
}

void lucWindowInteraction_InitAll( void* windowInteractor ) {
	lucWindowInteraction* self        = windowInteractor;

	_lucWindowInteraction_Init( self );
}
		
void _lucWindowInteraction_Delete( void* windowInteractor ) {
	lucWindowInteraction* self        = windowInteractor;
	
	_Stg_Component_Delete( self );
}

void _lucWindowInteraction_Print( void* windowInteractor, Stream* stream ) {
	lucWindowInteraction*          self        = windowInteractor;
	
	Journal_Printf( stream, "lucWindowInteraction: %s\n", self->name );

	Stream_Indent( stream );

	/* Print Parent */
	_Stg_Component_Print( self, stream );

	Stream_UnIndent( stream );
}

void* _lucWindowInteraction_Copy( void* windowInteractor, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	lucWindowInteraction* self        = windowInteractor;
	lucWindowInteraction* newWindowInteraction;

	newWindowInteraction = _Stg_Component_Copy( self, dest, deep, nameExt, ptrMap );

	/* TODO */
	abort();

	return (void*) newWindowInteraction;
}

void* _lucWindowInteraction_DefaultNew( Name name ) {
	return _lucWindowInteraction_New( 
			sizeof( lucWindowInteraction ),
			lucWindowInteraction_Type,
			_lucWindowInteraction_Delete,
			_lucWindowInteraction_Print,
			_lucWindowInteraction_Copy,
			_lucWindowInteraction_DefaultNew,
			_lucWindowInteraction_AssignFromXML,
			_lucWindowInteraction_Build,
			_lucWindowInteraction_Initialise,
			_lucWindowInteraction_Execute,
			_lucWindowInteraction_Destroy,
			_lucWindowInteraction_MouseMotion,
			_lucWindowInteraction_MouseClick,
			_lucWindowInteraction_MouseMessage,
			_lucWindowInteraction_KeyboardEvent,
			_lucWindowInteraction_KeyboardMessage,
			name );
}

void _lucWindowInteraction_AssignFromXML( void* windowInteractor, Stg_ComponentFactory* cf, void* data ) {
	lucWindowInteraction*          self        = windowInteractor;
	
	self->context = Stg_ComponentFactory_ConstructByKey( cf, self->name, "Context", AbstractContext, False, data );
	if( !self->context ) 
		self->context = Stg_ComponentFactory_ConstructByName( cf, "context", AbstractContext, True, data );
	
	_lucWindowInteraction_Init( self );
}

void _lucWindowInteraction_Build( void* windowInteractor, void* data ) { }
void _lucWindowInteraction_Initialise( void* windowInteractor, void* data ) { }
void _lucWindowInteraction_Execute( void* windowInteractor, void* data ) { }
void _lucWindowInteraction_Destroy( void* windowInteractor, void* data ) { }

/* Wrappers to virtual functions */
void lucWindowInteraction_MouseMotion( void* windowInteraction, lucWindow* window, lucMouseButton button, Pixel_Index xpos, Pixel_Index ypos, Pixel_Index startx, Pixel_Index starty) {
	lucWindowInteraction*     self      = (lucWindowInteraction*) windowInteraction;

	self->_mouseMotion( self, window, button, xpos, ypos, startx, starty );
}

void lucWindowInteraction_MouseClick( void* windowInteraction, lucWindow* window, lucMouseButton button, lucMouseState state, Pixel_Index xpos, Pixel_Index ypos) {
	lucWindowInteraction*     self      = (lucWindowInteraction*) windowInteraction;

	self->_mouseClick( self, window, button, state, xpos, ypos );
}

void lucWindowInteraction_MouseMessage( void* windowInteraction, Stream* stream ) {
	lucWindowInteraction*     self      = (lucWindowInteraction*) windowInteraction;

	self->_mouseMessage( self, stream );
}

void lucWindowInteraction_KeyboardEvent( void* windowInteraction, lucWindow* window, char key, Pixel_Index xpos, Pixel_Index ypos) {
	lucWindowInteraction*     self      = (lucWindowInteraction*) windowInteraction;

	self->_keyboardEvent( self, window, key, xpos, ypos );
}

void lucWindowInteraction_KeyboardMessage( void* windowInteraction, Stream* stream ) {
	lucWindowInteraction*     self      = (lucWindowInteraction*) windowInteraction;

	self->_keyboardMessage( self, stream );
}

/* Default Virtual Function Implementations */
void _lucWindowInteraction_MouseMotion( void* windowInteractor, lucWindow* window, lucMouseButton button, Pixel_Index xpos, Pixel_Index ypos, Pixel_Index startx, Pixel_Index starty) {
	lucWindowInteraction*     self      = (lucWindowInteraction*) windowInteractor;
	lucViewport*   viewport;
	lucCamera*     camera;
	
	lucDebug_PrintFunctionBegin( self, 2 );

	viewport = lucWindow_GetViewportByPixel( window, xpos, ypos );
	if (viewport == NULL) {
		lucDebug_PrintFunctionEnd( self, 2 );
		return;
	}
	camera = viewport->camera;
	
	switch (button) {
		case lucLeftButton:
			lucCamera_RotateAroundUpDirection(  camera, ((double)startx - (double)xpos) * M_PI/180.0 );
			lucCamera_RotateTowardsUpDirection( camera, ((double)starty - (double)ypos) * M_PI/180.0 );
			break;
		case lucRightButton: 
		    lucCamera_ChangeFocalPoint( camera, startx, starty, xpos, ypos);
			break;
		case lucMiddleButton:
			lucCamera_Zoom( camera, 1.0 + ((double)starty - (double)ypos)/10.0 );
			break;
		default:
			break;
	}
	lucDebug_PrintFunctionEnd( self, 2 );
}

void _lucWindowInteraction_MouseClick( void* windowInteractor, lucWindow* window, lucMouseButton button, lucMouseState state, Pixel_Index xpos, Pixel_Index ypos) {
	lucWindowInteraction*     self      = (lucWindowInteraction*) windowInteractor;
	lucViewport*   viewport;
	
	lucDebug_PrintFunctionBegin( self, 2 );

	viewport = lucWindow_GetViewportByPixel( window, xpos, ypos );
	if (viewport == NULL) 
		return;

	switch (state) {
		case lucButtonPress:
			break;
		case lucButtonRelease:
			break;
		default:
			break;
	}

	switch (button) {
		case lucLeftButton:
			break;
		case lucRightButton:
			break;
		case lucMiddleButton:
			break;
		case lucWheelUp:
			lucCamera_Zoom( viewport->camera, 0.9 );
			break;
		case lucWheelDown:
			lucCamera_Zoom( viewport->camera, 1.1 );
			break;
	}

	lucDebug_PrintFunctionEnd( self, 2 );
}	

void _lucWindowInteraction_MouseMessage( void* windowInteractor, Stream* stream ) {
	Journal_Printf( stream, "Left Button:                  Rotate camera around model.\n" );
	Journal_Printf( stream, "Middle Button:                Zoom camera in or out of the model.\n" );
	Journal_Printf( stream, "Right Button:                 Move camera's focal point.\n" );
}

void _lucWindowInteraction_KeyboardEvent( void* windowInteractor, lucWindow* window, char key, Pixel_Index xpos, Pixel_Index ypos) {
	lucWindowInteraction* self      = (lucWindowInteraction*) windowInteractor;
	lucViewport*          viewport = lucWindow_GetViewportByPixel( window, xpos, ypos );
	
	lucDebug_PrintFunctionBegin( self, 2 );

	/* Make Case Insensitive */
	/*key = tolower( key )*/;

	Journal_DPrintfL( lucDebug, 0, "Key pressed '%c' (%d)\n", key, key );

	switch (key) {
		case ESCAPE:
			Journal_Printf( Journal_MyStream( Info_Type, self ), "Escape key pressed - Exiting program.\n" );
			lucDebug_PrintFunctionEnd( self, 2 );
			lucWindow_ToggleApplicationQuit( window );
			break;
		case 'q': case ' ':
			lucDebug_PrintFunctionEnd( self, 2 );
			lucWindow_QuitEventLoop( window );
			break;
		case 'h': 
			lucWindow_InteractionHelpMessage( window, Journal_MyStream( Info_Type, window ) );     
			break;
		case 'i': 
			lucWindow_ChangeInteractiveMode( window );     
			Journal_Printf( Journal_MyStream( Info_Type, self ), "Interactivity for %s '%s' is %s.\n",
					window->type, window->name, window->interactive ? "True" : "False" );
			break;
		case 'r': 
			if (viewport)
				lucViewport_Reset( viewport );
			break;
		case 's': 
			if (viewport)
				lucCamera_Pickle( viewport->camera, Journal_MyStream( Info_Type, viewport->camera ) );
			break;
		case '[':
			if (viewport)
				lucCamera_Zoom( viewport->camera, 1.1 );
			break;
		case ']':
			if (viewport)
				lucCamera_Zoom( viewport->camera, 0.9 );
			break;
		case 'o':
			lucWindow_ChangeContinuousMode( window );
			Journal_Printf( Journal_MyStream( Info_Type, self ), "Continuous for %s '%s' is %s.\n",
					window->type, window->name, window->continuous ? "True" : "False" );
			lucWindow_QuitEventLoop( window );
			break;
		
		default:
			break;
	}

	lucDebug_PrintFunctionEnd( self, 2 );
}

void _lucWindowInteraction_KeyboardMessage( void* windowInteractor, Stream* stream ) {
	Journal_Printf( stream, "escape:                       Exit out of program.\n" );
	Journal_Printf( stream, "q, space:                     Close window and continue running program.\n" );
	Journal_Printf( stream, "h:                            Print this help message.\n" );
	Journal_Printf( stream, "i:                            Toggle Window Interactivity.\n" );
	Journal_Printf( stream, "r:                            Reset camera for viewport under cursor.\n" );
	Journal_Printf( stream, "s:                            Output information for camera associated with viewport under cursor.\n" );
	Journal_Printf( stream, "[:                            Zoom out with camera associated with viewport under cursor.\n" );
	Journal_Printf( stream, "]:                            Zoom in with camera associated with viewport under cursor.\n" );
	Journal_Printf( stream, "o:                            Continuous mode on/off (keep rendering new frames while interacting).\n" );
}
