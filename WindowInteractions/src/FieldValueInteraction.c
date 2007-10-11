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
** $Id: FieldValueInteraction.c 740 2007-10-11 08:05:31Z SteveQuenette $
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
#include "FieldValueInteraction.h"

#include <assert.h>

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type lucFieldValueInteraction_Type = "lucFieldValueInteraction";

/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
lucFieldValueInteraction* _lucFieldValueInteraction_New( 
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
	lucFieldValueInteraction*					self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( sizeOfSelf >= sizeof(lucFieldValueInteraction) );
	self = (lucFieldValueInteraction*) _lucWindowInteraction_New( 
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

void _lucFieldValueInteraction_Init( 
		lucFieldValueInteraction*                                      self  ) 
{
}

void _lucFieldValueInteraction_Delete( void* FieldValueInteraction ) {
	lucFieldValueInteraction*  self = (lucFieldValueInteraction*)FieldValueInteraction;

	_lucWindowInteraction_Delete( self );
}

void _lucFieldValueInteraction_Print( void* FieldValueInteraction, Stream* stream ) {
	lucFieldValueInteraction*  self = (lucFieldValueInteraction*)FieldValueInteraction;

	_lucWindowInteraction_Print( self, stream );
}

void* _lucFieldValueInteraction_Copy( void* FieldValueInteraction, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap) {
	lucFieldValueInteraction*  self = (lucFieldValueInteraction*)FieldValueInteraction;
	lucFieldValueInteraction* newFieldValueInteraction;

	newFieldValueInteraction = _lucWindowInteraction_Copy( self, dest, deep, nameExt, ptrMap );

	/* TODO */
	abort();

	return (void*) newFieldValueInteraction;
}


void* _lucFieldValueInteraction_DefaultNew( Name name ) {
	return (void*) _lucFieldValueInteraction_New(
		sizeof(lucFieldValueInteraction),
		lucFieldValueInteraction_Type,
		_lucFieldValueInteraction_Delete,
		_lucFieldValueInteraction_Print,
		NULL,
		_lucFieldValueInteraction_DefaultNew,
		_lucFieldValueInteraction_Construct,
		_lucFieldValueInteraction_Build,
		_lucFieldValueInteraction_Initialise,
		_lucFieldValueInteraction_Execute,
		_lucFieldValueInteraction_Destroy,		
		_lucFieldValueInteraction_MouseMotion,
		_lucFieldValueInteraction_MouseClick,
		_lucFieldValueInteraction_MouseMessage,
		_lucFieldValueInteraction_KeyboardEvent,
		_lucFieldValueInteraction_KeyboardMessage,
		name );
}

void _lucFieldValueInteraction_Construct( void* FieldValueInteraction, Stg_ComponentFactory* cf, void* data ){
	lucFieldValueInteraction*  self = (lucFieldValueInteraction*)FieldValueInteraction;

	/* Construct Parent */
	_lucWindowInteraction_Construct( self, cf, data );
	
	_lucFieldValueInteraction_Init( self );
}

void _lucFieldValueInteraction_Build( void* renderingEngine, void* data ) {}
void _lucFieldValueInteraction_Initialise( void* renderingEngine, void* data ) {}
void _lucFieldValueInteraction_Execute( void* renderingEngine, void* data ) {}
void _lucFieldValueInteraction_Destroy( void* renderingEngine, void* data ) {}

/* This component doesn't use any mouse interaction */
void _lucFieldValueInteraction_MouseMotion( void* windowInteraction, lucWindow* window, lucMouseButton button, Pixel_Index xpos, Pixel_Index ypos, Pixel_Index startx, Pixel_Index starty) {}
void _lucFieldValueInteraction_MouseClick( void* windowInteraction, lucWindow* window, lucMouseButton button, lucMouseState state, Pixel_Index xpos, Pixel_Index ypos) { } 
void _lucFieldValueInteraction_MouseMessage( void* windowInteraction, Stream* stream ) { }


void _lucFieldValueInteraction_KeyboardEvent( void* windowInteraction, lucWindow* window, char key, Pixel_Index xpos, Pixel_Index ypos) {
	lucFieldValueInteraction*   self = (lucFieldValueInteraction*) windowInteraction;
	lucViewportInfo*            viewportInfo;
	lucViewport*                viewport;
	Coord                       coord;
	lucDrawingObject*           object;
	lucScalarFieldCrossSection* scalarField;
	Stream*                     stream = Journal_MyStream( Info_Type, self );
	DrawingObject_Index         object_I;
	DrawingObject_Index         objectCount;
	double                      value;
	FieldVariable*              fieldVariable;

	/* This function works when the key pressed is 'f' */
	if ( key != 'f' )
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

		/* Check if this drawing object is a scalar field */
		if ( Stg_Class_IsInstance( object, lucScalarFieldCrossSection_Type ) ) {
			scalarField = (lucScalarFieldCrossSection*) object;
			fieldVariable = scalarField->fieldVariable;

			/* Get Value of the field, and store which processor it is on */
			if ( FieldVariable_InterpolateValueAt( fieldVariable, coord, &value ) == LOCAL ) {
				Journal_Printf( stream, "Field '%s' (type %s) at coord (%g, %g, %g) is %g\n", 
						fieldVariable->name, fieldVariable->type, 
						coord[ I_AXIS ], coord[ J_AXIS ], coord[ K_AXIS ], value );
			}
		}
	}

	/* Loop through other FieldVariables that are registered on this object */
}

void _lucFieldValueInteraction_KeyboardMessage( void* windowInteraction, Stream* stream ) {
	Journal_Printf( stream,
			"f:                            The value of fields under the cursor with be printed to screen.\n" );
}


