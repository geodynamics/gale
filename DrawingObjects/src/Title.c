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
** $Id: Title.c 510 2006-02-17 04:33:32Z RobertTurnbull $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>

#include <glucifer/Base/Base.h>
#include <glucifer/RenderingEngines/RenderingEngines.h>

#ifdef HAVE_GL2PS
	#include <gl2ps.h>
#endif


#include "types.h"
#include "OpenGLDrawingObject.h"
#include "Title.h"

#include <assert.h>
#ifdef HAVE_OPENGL_FRAMEWORK
	#include <OpenGL/gl.h>
	#include <OpenGL/glu.h>
#else
	#include <gl.h>
	#include <glu.h>
#endif
#include <string.h>

#ifndef MASTER
	#define MASTER 0
#endif



const Type lucTitle_Type = "lucTitle";

lucTitle* lucTitle_New( 
		Name                                               name,
    		char *                                             titleString,
		lucColour                                          colour,
		int                                                yPos )
{
	lucTitle* self = (lucTitle*) _lucTitle_DefaultNew( name );

	lucTitle_InitAll( self, titleString, colour, yPos );

	return self;
}

lucTitle* _lucTitle_New(
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
		lucDrawingObject_SetupFunction*                    _setup,
		lucDrawingObject_DrawFunction*                     _draw,
		lucDrawingObject_CleanUpFunction*                  _cleanUp,
		lucOpenGLDrawingObject_BuildDisplayListFunction*   _buildDisplayList,
	
		Name                                               name )
{
	lucTitle*    self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( sizeOfSelf >= sizeof(lucTitle) );
	self = (lucTitle*)  _lucOpenGLDrawingObject_New( 
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
			_setup,
			_draw,
			_cleanUp,
			_buildDisplayList,
			name );
	
	
	return self;
}

void lucTitle_Init(		
		lucTitle*                                          self,
		char*                                              titleString,
		lucColour                                          colour,
		int                                                yPos )
{
	
	self->titleString = Memory_Alloc_Array( char, (strlen(titleString) + 1), "self->titleString" );
	strcpy( self->titleString, titleString );
	memcpy( &(self->colour), &colour, sizeof(lucColour) );
	self->yPos = yPos;
}

void lucTitle_InitAll( 
		void*                                              title,
		char*                                              titleString,
		lucColour                                          colour,
		int                                                yPos )
{
	lucTitle* self        = title;

	/* TODO Init parent */
	lucTitle_Init( self, titleString, colour, yPos );
}

void _lucTitle_Delete( void* drawingObject ) {
	lucTitle*  self = (lucTitle*)drawingObject;

	Memory_Free( self->titleString );
	_lucOpenGLDrawingObject_Delete( self );
}

void _lucTitle_Print( void* drawingObject, Stream* stream ) {
	lucTitle*  self = (lucTitle*)drawingObject;

	_lucOpenGLDrawingObject_Print( self, stream );
}

void* _lucTitle_Copy( void* title, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	lucTitle* self        = title;
	lucTitle* newTitle;

	newTitle = _Stg_Component_Copy( self, dest, deep, nameExt, ptrMap );

	newTitle->titleString = StG_Strdup( self->titleString );
	memcpy( &(newTitle->colour),       &(self->colour),       sizeof(lucColour) );

	
	return (void*) newTitle;
}

void* _lucTitle_DefaultNew( Name name ) {
	return _lucTitle_New( 
			sizeof( lucTitle ),
			lucTitle_Type,
			_lucTitle_Delete,
			_lucTitle_Print,
			_lucTitle_Copy,
			_lucTitle_DefaultNew,
			_lucTitle_Construct,
			_lucTitle_Build,
			_lucTitle_Initialise,
			_lucTitle_Execute,
			_lucTitle_Destroy,		
		        _lucTitle_Setup,
			_lucTitle_Draw,
	                _lucTitle_CleanUp,
	 		_lucTitle_BuildDisplayList,
			name );
}

void _lucTitle_Construct( void* title, Stg_ComponentFactory* cf, void* data ) {
	lucTitle*             self               = (lucTitle*) title;
        Name colourName;
	
	/* Get Stereo Type */
         /* Construct Parent */
	_lucDrawingObject_Construct( self, cf, data );

	colourName  = Stg_ComponentFactory_GetString( cf, self->name, "colour", "Black") ;
	
	lucColour_FromString( &self->colour, colourName );

	lucTitle_InitAll( self, 
			Stg_ComponentFactory_GetString( cf, self->name, "string", "" ),
		        self->colour,
			Stg_ComponentFactory_GetInt( cf, self->name, "yPos", 20 ) );
			
}

void _lucTitle_Build( void* Title, void* data ) { }
void _lucTitle_Initialise( void* Title, void* data ) { }
void _lucTitle_Execute( void* Title, void* data ) { }
void _lucTitle_Destroy( void* Title, void* data ) { }

void _lucTitle_Setup( void* drawingObject, void* _context ) {
	lucTitle*       self            = (lucTitle*)drawingObject;
	_lucOpenGLDrawingObject_Setup( self, _context );
}

void _lucTitle_Draw( void* drawingObject, lucWindow* window, lucViewportInfo* viewportInfo, void* _context ) {
	lucTitle*    self            = (lucTitle*)drawingObject;
    lucViewport* viewport        = viewportInfo->viewport;
	char*        title;
	
	/* Set up 2D Viewer the size of the viewport */
	lucViewport2d(True, viewportInfo);

	/* Set title */
	title = strlen( self->titleString ) ? self->titleString : viewport->name;

	/* Print white first for contrast */
	lucColour white = {1.0f, 1.0f, 1.0f, 1.0f};
	lucColour_SetOpenGLColour( &white );
	lucPrint(viewportInfo->width/2 - lucStringWidth(title)/2 - 1, self->yPos - 1, title );

	/* Set the colour */
	lucColour_SetOpenGLColour( &self->colour );

	/* Print Title */
	lucPrint(viewportInfo->width/2 - lucStringWidth(title)/2, self->yPos, title );

	/* Restore the viewport */
	lucViewport2d(False, viewportInfo);
}

void _lucTitle_CleanUp( void* drawingObject, void* _context ) {
	lucTitle*      self            = (lucTitle*)drawingObject;
	
	_lucOpenGLDrawingObject_CleanUp( self, _context );

}

void _lucTitle_BuildDisplayList( void* drawingObject, void* _context ) {
	}





