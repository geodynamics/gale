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
** $Id: OpenGLDrawingObject.c 791 2008-09-01 02:09:06Z JulianGiordani $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>

#include <glucifer/Base/Base.h>
#include <glucifer/RenderingEngines/RenderingEngines.h>

#include "types.h"
#include "OpenGLDrawingObject.h"

#include <assert.h>
#include <gl.h>
#include <glu.h>
#include <string.h>

#ifndef MASTER
	#define MASTER 0
#endif

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type lucOpenGLDrawingObject_Type = "lucOpenGLDrawingObject";

/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
lucOpenGLDrawingObject* _lucOpenGLDrawingObject_New( 
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
	lucOpenGLDrawingObject*					self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( sizeOfSelf >= sizeof(lucOpenGLDrawingObject) );
	self = (lucOpenGLDrawingObject*) _lucDrawingObject_New( 
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
			name );

	self->_buildDisplayList = _buildDisplayList;
	
	return self;
}

void _lucOpenGLDrawingObject_Init( 
		lucOpenGLDrawingObject*                            self )
{
    self->displayList = 0;
}

void _lucOpenGLDrawingObject_Delete( void* drawingObject ) {
	lucOpenGLDrawingObject*  self = (lucOpenGLDrawingObject*)drawingObject;

	_lucDrawingObject_Delete( self );
}

void _lucOpenGLDrawingObject_Print( void* drawingObject, Stream* stream ) {
	lucOpenGLDrawingObject*  self = (lucOpenGLDrawingObject*)drawingObject;

	_lucDrawingObject_Print( self, stream );

	Journal_PrintValue( stream, self->displayList );
}

void* _lucOpenGLDrawingObject_Copy( void* drawingObject, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap) {
	lucOpenGLDrawingObject*  self = (lucOpenGLDrawingObject*)drawingObject;
	lucOpenGLDrawingObject* newDrawingObject;

	newDrawingObject = _lucDrawingObject_Copy( self, dest, deep, nameExt, ptrMap );

	newDrawingObject->displayList = self->displayList;

	return (void*) newDrawingObject;
}



void _lucOpenGLDrawingObject_AssignFromXML( void* drawingObject, Stg_ComponentFactory* cf, void* data ){
	lucOpenGLDrawingObject*         self               = (lucOpenGLDrawingObject*)drawingObject;

	/* Construct Parent */
	_lucDrawingObject_AssignFromXML( self, cf, data );

	_lucOpenGLDrawingObject_Init( self );
}

void _lucOpenGLDrawingObject_Build( void* drawingObject, void* data ) {}
void _lucOpenGLDrawingObject_Initialise( void* drawingObject, void* data ) {}
void _lucOpenGLDrawingObject_Execute( void* drawingObject, void* data ) {}
void _lucOpenGLDrawingObject_Destroy( void* drawingObject, void* data ) {}

/* Drawing Object implementations */
void _lucOpenGLDrawingObject_Setup( void* drawingObject, void* context ) {
	lucOpenGLDrawingObject*           self            = (lucOpenGLDrawingObject*)drawingObject;

    /* Generate a display list id */
	if (!self->displayList)	self->displayList = glGenLists( 1 );

	/* Create/replace OpenGL display list */
	glNewList( self->displayList, GL_COMPILE);

	/* Run the virtual function for building the display list - 
	 * this should contain as much of the opengl drawing primitives as possible */
	lucOpenGLDrawingObject_BuildDisplayList( self, context );
	
	/* Tell OpenGL that we've finished creating the list now */
	glEndList();	

}
	
void _lucOpenGLDrawingObject_Draw( void* drawingObject, lucWindow* window, lucViewportInfo* viewportInfo, void* _context ) {
	lucOpenGLDrawingObject*           self          = (lucOpenGLDrawingObject*)drawingObject;

	/* We should make sure that the rendering engine for this window is the lucRenderingEngineGL */
	Journal_Firewall( 
			Stg_Class_IsInstance( window->renderingEngine, lucRenderingEngineGL_Type ),
			self->errorStream,
			"Error for %s '%s' - This class only works with rendering engines of type %s.\n"
			"%s '%s' is using a rendering engine of type %s. Please correct this.\n", 
			self->type, self->name, lucRenderingEngineGL_Type, window->type, window->name, window->renderingEngine->type );

	/* All that we need to do to visualise this object now is to call the display list
	 * this should have been created in the setup phase */
	if (self->displayList) glCallList( self->displayList );
}

void _lucOpenGLDrawingObject_CleanUp( void* drawingObject, void* _context ) {
	lucOpenGLDrawingObject*           self          = (lucOpenGLDrawingObject*)drawingObject;
	
	if (self->displayList) glDeleteLists( self->displayList, 1 );
    self->displayList = 0;
}

/* Wrappers for virtual functions */
void lucOpenGLDrawingObject_BuildDisplayList( void* drawingObject, void* context ) {
	lucOpenGLDrawingObject*           self            = (lucOpenGLDrawingObject*)drawingObject;

	self->_buildDisplayList( self, context );
}

/* HACK - a function to check whether a field is an FeVariable or not before it does an FeVariable_SyncShadowValues */
void lucOpenGLDrawingObject_SyncShadowValues( void* drawingObject, void* field ) {
	if ( field && Stg_Class_IsInstance( field, FeVariable_Type ) )
		FeVariable_SyncShadowValues( field );
}
