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
** $Id: ColourBarInteraction.c 510 2006-02-17 04:33:32Z RobertTurnbull $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>

#include <glucifer/Base/Base.h>
#include <glucifer/Windowing/Windowing.h>
#include <glucifer/RenderingEngines/RenderingEngines.h>
#include <glucifer/DrawingObjects/DrawingObjects.h>

#include "types.h"
#include "ColourBarInteraction.h"

#include <assert.h>

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type lucColourBarInteraction_Type = "lucColourBarInteraction";

/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
lucColourBarInteraction* _lucColourBarInteraction_New( 
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
	lucColourBarInteraction*					self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( sizeOfSelf >= sizeof(lucColourBarInteraction) );
	self = (lucColourBarInteraction*) _lucWindowInteraction_New( 
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
			_mouseMotion,
			_mouseClick, 
			_mouseMessage, 
			_keyboardEvent,
			_keyboardMessage,
			name );
	
	return self;
}

void _lucColourBarInteraction_Init( 
		lucColourBarInteraction*                                      self  ) 
{
}

void _lucColourBarInteraction_Delete( void* ColourBarInteraction ) {
	lucColourBarInteraction*  self = (lucColourBarInteraction*)ColourBarInteraction;

	_lucWindowInteraction_Delete( self );
}

void _lucColourBarInteraction_Print( void* ColourBarInteraction, Stream* stream ) {
	lucColourBarInteraction*  self = (lucColourBarInteraction*)ColourBarInteraction;

	_lucWindowInteraction_Print( self, stream );
}

void* _lucColourBarInteraction_Copy( void* ColourBarInteraction, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap) {
	lucColourBarInteraction*  self = (lucColourBarInteraction*)ColourBarInteraction;
	lucColourBarInteraction* newColourBarInteraction;

	newColourBarInteraction = _lucWindowInteraction_Copy( self, dest, deep, nameExt, ptrMap );

	/* TODO */
	abort();

	return (void*) newColourBarInteraction;
}


void* _lucColourBarInteraction_DefaultNew( Name name ) {
	return (void*) _lucColourBarInteraction_New(
		sizeof(lucColourBarInteraction),
		lucColourBarInteraction_Type,
		_lucColourBarInteraction_Delete,
		_lucColourBarInteraction_Print,
		NULL,
		_lucColourBarInteraction_DefaultNew,
		_lucColourBarInteraction_AssignFromXML,
		_lucColourBarInteraction_Build,
		_lucColourBarInteraction_Initialise,
		_lucColourBarInteraction_Execute,
		_lucColourBarInteraction_Destroy,		
		_lucColourBarInteraction_MouseMotion,
		_lucColourBarInteraction_MouseClick,
		_lucColourBarInteraction_MouseMessage,
		_lucColourBarInteraction_KeyboardEvent,
		_lucColourBarInteraction_KeyboardMessage,
		name );
}

void _lucColourBarInteraction_AssignFromXML( void* ColourBarInteraction, Stg_ComponentFactory* cf, void* data ){
	lucColourBarInteraction*  self = (lucColourBarInteraction*)ColourBarInteraction;

	/* Construct Parent */
	_lucWindowInteraction_AssignFromXML( self, cf, data );
	
	_lucColourBarInteraction_Init( self );
}

void _lucColourBarInteraction_Build( void* renderingEngine, void* data ) {}
void _lucColourBarInteraction_Initialise( void* renderingEngine, void* data ) {}
void _lucColourBarInteraction_Execute( void* renderingEngine, void* data ) {}
void _lucColourBarInteraction_Destroy( void* renderingEngine, void* data ) {}

/* This component doesn't use any mouse interaction */
void _lucColourBarInteraction_MouseMotion( void* windowInteraction, lucWindow* window, lucMouseButton button, Pixel_Index xpos, Pixel_Index ypos, Pixel_Index startx, Pixel_Index starty) {}
void _lucColourBarInteraction_MouseClick( void* windowInteraction, lucWindow* window, lucMouseButton button, lucMouseState state, Pixel_Index xpos, Pixel_Index ypos) { } 
void _lucColourBarInteraction_MouseMessage( void* windowInteraction, Stream* stream ) { }


void _lucColourBarInteraction_KeyboardEvent( void* windowInteraction, lucWindow* window, char key, Pixel_Index xpos, Pixel_Index ypos) {
	lucColourBarInteraction*   self = (lucColourBarInteraction*) windowInteraction;
	lucViewportInfo*            viewportInfo;
	lucViewport*                viewport;
	Coord                       coord;
	lucDrawingObject*           object;
	Stream*                     stream = Journal_MyStream( Info_Type, self );
	DrawingObject_Index         object_I;
	DrawingObject_Index         objectCount;
	double                      value;
	lucColourBar*               colourBar;
	Pixel_Index                 length;
	Pixel_Index                 height;
	Pixel_Index                 pixel_I;
	int                         startPos[2];


	/* This function works when the key pressed is 'c' */
	if ( key != 'c' )
		return;
	
	/* Find which viewport the user clicked on */
	viewportInfo = lucWindow_GetViewportInfoByPixel( window, xpos, ypos );
	if ( viewportInfo == NULL )
		return; /* If the user hasn't clicked on a viewport at all, then return */

	viewport = viewportInfo->viewport;

	/* Get spatial coordinate that the user clicked on */
	lucViewportInfo_GetCoordFromPixel( viewportInfo, xpos, ypos, coord );

	/* Loop through lucScalarFields that are registered on the viewport */
	objectCount   = lucDrawingObject_Register_GetCount( viewport->drawingObject_Register );
	for ( object_I = 0 ; object_I < objectCount ; object_I++ ) {
		object = lucDrawingObject_Register_GetByIndex( viewport->drawingObject_Register, object_I );

		/* Check if this drawing object is a colourBar type */
		if ( Stg_Class_IsInstance( object, lucColourBar_Type ) ) {
			colourBar = (lucColourBar*) object;
			length        = (Pixel_Index) ((double) viewportInfo->width * colourBar->lengthFactor);
			height         = colourBar->height;

			/* Finds if the position is in the colourBar */		
			startPos[0] = (viewportInfo->width - length)/2 + viewportInfo->startx;
			startPos[1] = colourBar->margin + viewportInfo->starty;

			if ( (xpos >= startPos[0]) && ( xpos <= startPos[0] + length) && (ypos >= startPos[1]) && (ypos <= startPos[1]+height)){
				/* finds the value of the pixel */
				pixel_I =  xpos - startPos[0];
				value = colourBar->colourMap->minimum + (double) pixel_I * (colourBar->colourMap->maximum - colourBar->colourMap->minimum) / (double) length;
			
		    	Journal_Printf( stream, "Value for this colour is  '%.4E' \n ",value );
			}
			else Journal_Printf( stream, "The mouse should be positionned on the colourBar \n " );
		}
	}
}

void _lucColourBarInteraction_KeyboardMessage( void* windowInteraction, Stream* stream ) {
	Journal_Printf( stream,
			"c:                            The value of ColourBar's colour under the cursor will be printed to screen.\n" );
}


