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
** $Id: DrawingObject_Register.c 740 2007-10-11 08:05:31Z SteveQuenette $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>

#include "types.h"
#include "ViewportInfo.h"
#include "Viewport.h"
#include "ColourMap.h"
#include "Window.h"
#include "RenderingEngine.h"
#include "DrawingObject_Register.h"
#include "DrawingObject.h"
#include "Light_Register.h"


const Type lucDrawingObject_Register_Type = "lucDrawingObject_Register";


lucDrawingObject_Register*	lucDrawingObject_Register_New( void ) {
	lucDrawingObject_Register* self;
	
	self = (lucDrawingObject_Register*) _NamedObject_Register_New(
		sizeof(lucDrawingObject_Register),
		lucDrawingObject_Register_Type,
		_NamedObject_Register_Delete,
		_NamedObject_Register_Print,
		_NamedObject_Register_Copy );

	return self;
}

void lucDrawingObject_Register_SetNeedsToSetupFlag( void* drawingObject_Register, Bool flag ) {
	lucDrawingObject_Register* self           = (lucDrawingObject_Register*) drawingObject_Register;
	DrawingObject_Index        object_I;
	DrawingObject_Index        objectCount    = lucDrawingObject_Register_GetCount( self );
	lucDrawingObject*          object;

	for ( object_I = 0 ; object_I < objectCount ; object_I++ ) {
		object = lucDrawingObject_Register_GetByIndex( self, object_I );

		object->needsToSetup = flag;
	}
}

void lucDrawingObject_Register_DrawAll( void* drawingObject_Register, lucWindow* window, lucViewportInfo* viewportInfo, void* context, Bool compositeEachDraw ) {
	lucDrawingObject_Register* self          = (lucDrawingObject_Register*) drawingObject_Register;
	DrawingObject_Index        object_I;
	DrawingObject_Index        objectCount   = lucDrawingObject_Register_GetCount( self );
	lucDrawingObject*          object;
	lucViewport*               viewport      =  (lucViewport*) viewportInfo->viewport;
	lucLight_Register*         lightRegister =  (lucLight_Register*) (viewport->light_Register);

	for ( object_I = 0 ; object_I < objectCount ; object_I++ ) {
		object = lucDrawingObject_Register_GetByIndex( self, object_I );
		lucLight_Register_EnableAll( lightRegister );
		lucDrawingObject_Draw( object, window, viewportInfo, context );

		if ( compositeEachDraw )
			lucRenderingEngine_CompositeViewport( window->renderingEngine, viewportInfo, context, True );
	}

	if ( !compositeEachDraw )
		lucRenderingEngine_CompositeViewport( window->renderingEngine, viewportInfo, context, False );
}

void lucDrawingObject_Register_CleanUpAll( void* drawingObject_Register, void* context ) {
	lucDrawingObject_Register* self          = (lucDrawingObject_Register*) drawingObject_Register;
	DrawingObject_Index        object_I;
	DrawingObject_Index        objectCount   = lucDrawingObject_Register_GetCount( self );
	lucDrawingObject*          object;

	for ( object_I = 0 ; object_I < objectCount ; object_I++ ) {
		object = lucDrawingObject_Register_GetByIndex( self, object_I );
		lucDrawingObject_CleanUp( object, context );
	}
}
