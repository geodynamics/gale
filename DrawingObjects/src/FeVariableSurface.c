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
** $Id: FeVariableSurface.c 694 2007-05-22 03:17:41Z LukeHodkinson $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgFEM/StgFEM.h>

#include <glucifer/Base/Base.h>
#include <glucifer/RenderingEngines/RenderingEngines.h>

#include "types.h"
#include "OpenGLDrawingObject.h"
#include "FeVariableSurface.h"

#include <assert.h>
#include <gl.h>
#include <glu.h>
#include <string.h>

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type lucFeVariableSurface_Type = "lucFeVariableSurface";

/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
lucFeVariableSurface* _lucFeVariableSurface_New( 
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
	lucFeVariableSurface*					self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( sizeOfSelf >= sizeof(lucFeVariableSurface) );
	self = (lucFeVariableSurface*) _lucOpenGLDrawingObject_New( 
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

void _lucFeVariableSurface_Init( 
		lucFeVariableSurface*                                        self,
		FieldVariable*                                               feVariable,
		lucColourMap*                                                colourMap,
		Name                                                         colourName,
		Bool                                                         wireframe,
		float                                                        lineWidth,
		float                                                        scaleHeight )
{
	self->feVariable  = feVariable;
	self->colourMap   = colourMap;
	lucColour_FromString( &self->colour, colourName );
	self->wireframe   = wireframe;
	self->lineWidth   = lineWidth;
	self->scaleHeight = scaleHeight;

	assert( Stg_Class_IsInstance( feVariable, FeVariable_Type ) );
}

void _lucFeVariableSurface_Delete( void* drawingObject ) {
	lucFeVariableSurface*  self = (lucFeVariableSurface*)drawingObject;

	_lucOpenGLDrawingObject_Delete( self );
}

void _lucFeVariableSurface_Print( void* drawingObject, Stream* stream ) {
	lucFeVariableSurface*  self = (lucFeVariableSurface*)drawingObject;

	_lucOpenGLDrawingObject_Print( self, stream );
}

void* _lucFeVariableSurface_Copy( void* drawingObject, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap) {
	lucFeVariableSurface*  self = (lucFeVariableSurface*)drawingObject;
	lucFeVariableSurface* newDrawingObject;

	newDrawingObject = _lucOpenGLDrawingObject_Copy( self, dest, deep, nameExt, ptrMap );

	/* TODO */
	abort();

	return (void*) newDrawingObject;
}


void* _lucFeVariableSurface_DefaultNew( Name name ) {
	return (void*) _lucFeVariableSurface_New(
		sizeof(lucFeVariableSurface),
		lucFeVariableSurface_Type,
		_lucFeVariableSurface_Delete,
		_lucFeVariableSurface_Print,
		NULL,
		_lucFeVariableSurface_DefaultNew,
		_lucFeVariableSurface_Construct,
		_lucFeVariableSurface_Build,
		_lucFeVariableSurface_Initialise,
		_lucFeVariableSurface_Execute,
		_lucFeVariableSurface_Destroy,
		_lucFeVariableSurface_Setup,
		_lucFeVariableSurface_Draw,
		_lucFeVariableSurface_CleanUp,
		_lucFeVariableSurface_BuildDisplayList,
		name );
}

void _lucFeVariableSurface_Construct( void* drawingObject, Stg_ComponentFactory* cf, void* data ){
	lucFeVariableSurface*  self = (lucFeVariableSurface*)drawingObject;
	FieldVariable*         feVariable;
	lucColourMap*          colourMap;

	/* Construct Parent */
	_lucDrawingObject_Construct( self, cf, data );

	feVariable    =  Stg_ComponentFactory_ConstructByKey( cf, self->name, "FeVariable", FieldVariable, True,  data );
	colourMap     =  Stg_ComponentFactory_ConstructByKey( cf, self->name, "ColourMap",  lucColourMap,  False, data );
	
	_lucFeVariableSurface_Init( 
			self, 
			feVariable,
			colourMap,
			Stg_ComponentFactory_GetString( cf, self->name, "colour", "black" ),
			Stg_ComponentFactory_GetBool( cf, self->name, "wireframe", False ),
			(float) Stg_ComponentFactory_GetDouble( cf, self->name, "lineWidth", 1.0 ),
			(float) Stg_ComponentFactory_GetDouble( cf, self->name, "scaleHeight", 0.0 ) );
}

void _lucFeVariableSurface_Build( void* drawingObject, void* data ) {}
void _lucFeVariableSurface_Initialise( void* drawingObject, void* data ) {}
void _lucFeVariableSurface_Execute( void* drawingObject, void* data ) {}
void _lucFeVariableSurface_Destroy( void* drawingObject, void* data ) {}

void _lucFeVariableSurface_Setup( void* drawingObject, void* _context ) {
	lucFeVariableSurface*       self            = (lucFeVariableSurface*)drawingObject;
	_lucOpenGLDrawingObject_Setup( self, _context );
}
	
void _lucFeVariableSurface_Draw( void* drawingObject, lucWindow* window, lucViewportInfo* viewportInfo, void* _context ) {
	lucFeVariableSurface*       self            = (lucFeVariableSurface*)drawingObject;
	FeVariable*                    feVariable         = (FeVariable*) self->feVariable;
	lucColourMap*                  colourMap          = self->colourMap;
	if ( colourMap )
		lucColourMap_CalibrateFromFieldVariable( colourMap, feVariable );
	else
		lucColour_SetOpenGLColour( &self->colour );
	_lucOpenGLDrawingObject_Draw( self, window, viewportInfo, _context );
}


void _lucFeVariableSurface_CleanUp( void* drawingObject, void* _context ) {
	lucFeVariableSurface*       self            = (lucFeVariableSurface*)drawingObject;
	_lucOpenGLDrawingObject_CleanUp( self, _context );
}
	
void _lucFeVariableSurface_BuildDisplayList( void* drawingObject, void* _context ) {
	lucFeVariableSurface*          self               = (lucFeVariableSurface*)drawingObject;
	FeVariable*                    feVariable         = (FeVariable*) self->feVariable;
	FeMesh*    		       mesh               = feVariable->feMesh;
	lucColourMap*                  colourMap          = self->colourMap;
	Element_LocalIndex             lElement_I;
	Element_LocalIndex             elementLocalCount  = FeMesh_GetElementLocalSize( mesh );
	Element_NodeIndex              eNode_I;
	Element_NodeIndex              elementNodeCount, *elementNodes;
	Node_LocalIndex                lNode_I;
	double                         nodeValue;
	double                         height;

	FeVariable_SyncShadowValues( feVariable );

	/* Give option to draw surface as wireframe */
	if (self->wireframe) 
		glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );
	else 
		glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );

	glNormal3f( 0.0, 0.0, 1.0 ); /* TODO - FIX */
	glLineWidth( self->lineWidth );
	glDisable( GL_LIGHTING );

	if ( colourMap )
		lucColourMap_CalibrateFromFieldVariable( colourMap, feVariable );
	else
		lucColour_SetOpenGLColour( &self->colour );

	for ( lElement_I = 0 ; lElement_I < elementLocalCount ; lElement_I++ ) {
		FeMesh_GetElementNodes( mesh, lElement_I, &elementNodeCount, &elementNodes );

		glBegin( GL_POLYGON );
		for ( eNode_I = 0 ; eNode_I < elementNodeCount ; eNode_I++ ) {
			lNode_I = elementNodes[ eNode_I ];
			/* Get Value at node */
			nodeValue = FeVariable_GetScalarAtNode( feVariable, lNode_I );

			/* Change Colour */
			if ( colourMap )
				lucColourMap_SetOpenGLColourFromValue( colourMap, nodeValue );

			/* Set Height */
			height = nodeValue * self->scaleHeight;

			/* Plot Vertex */
			glVertex3d( mesh->verts[ lNode_I ][ I_AXIS ], mesh->verts[ lNode_I ][ J_AXIS ], height );
		}
		glEnd();
	}
	glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
	glEnable( GL_LIGHTING );	
}
