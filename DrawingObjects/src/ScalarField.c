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
** $Id: ScalarField.c 786 2008-08-18 13:55:47Z LukeHodkinson $
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
#include "ScalarFieldCrossSection.h"
#include "ScalarField.h"

#include <assert.h>
#include <gl.h>
#include <glu.h>
#include <string.h>

#ifndef MASTER
	#define MASTER 0
#endif

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type lucScalarField_Type = "lucScalarField";

/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
lucScalarField* _lucScalarField_New( 
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
	lucScalarField*					self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( sizeOfSelf >= sizeof(lucScalarField) );
	self = (lucScalarField*) _lucScalarFieldCrossSection_New( 
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

void _lucScalarField_Init( 
		lucScalarField*                                              self,
		Bool                                                         cullFace )
{
	self->cullFace = cullFace;
}

void _lucScalarField_Delete( void* drawingObject ) {
	lucScalarField*  self = (lucScalarField*)drawingObject;

	_lucScalarFieldCrossSection_Delete( self );
}

void _lucScalarField_Print( void* drawingObject, Stream* stream ) {
	lucScalarField*  self = (lucScalarField*)drawingObject;

	_lucScalarFieldCrossSection_Print( self, stream );
}

void* _lucScalarField_Copy( void* drawingObject, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap) {
	lucScalarField*  self = (lucScalarField*)drawingObject;
	lucScalarField* newDrawingObject;

	newDrawingObject = _lucScalarFieldCrossSection_Copy( self, dest, deep, nameExt, ptrMap );

	/* TODO */
	abort();

	return (void*) newDrawingObject;
}


void* _lucScalarField_DefaultNew( Name name ) {
	return (void*) _lucScalarField_New(
		sizeof(lucScalarField),
		lucScalarField_Type,
		_lucScalarField_Delete,
		_lucScalarField_Print,
		NULL,
		_lucScalarField_DefaultNew,
		_lucScalarField_Construct,
		_lucScalarField_Build,
		_lucScalarField_Initialise,
		_lucScalarField_Execute,
		_lucScalarField_Destroy,
		_lucScalarField_Setup,
		_lucScalarField_Draw,
		_lucScalarField_CleanUp,
		_lucScalarField_BuildDisplayList,
		name );
}

void _lucScalarField_Construct( void* drawingObject, Stg_ComponentFactory* cf, void* data ){
	lucScalarField*  self = (lucScalarField*)drawingObject;

	/* Construct Parent */
	_lucScalarFieldCrossSection_Construct( self, cf, data );

	_lucScalarField_Init( 
			self, 
			Stg_ComponentFactory_GetBool( cf, self->name, "cullFace", True ) );

        self->useMesh = Stg_ComponentFactory_GetBool( cf, self->name, "useMesh", False );
}

void _lucScalarField_Build( void* drawingObject, void* data ) {
	lucScalarField*  self = (lucScalarField*)drawingObject;

	/* Call parent function */
	_lucScalarFieldCrossSection_Build( self, data );
}
void _lucScalarField_Initialise( void* drawingObject, void* data ) {
	lucScalarField*  self = (lucScalarField*)drawingObject;

	/* Call parent function */
	_lucScalarFieldCrossSection_Initialise( self, data );
}
void _lucScalarField_Execute( void* drawingObject, void* data ) {}
void _lucScalarField_Destroy( void* drawingObject, void* data ) {}

void _lucScalarField_Setup( void* drawingObject, void* _context ) {
	lucScalarField*          self          = (lucScalarField*)drawingObject;
	
	_lucScalarFieldCrossSection_Setup( self, _context );
}
	
void _lucScalarField_Draw( void* drawingObject, lucWindow* window, lucViewportInfo* viewportInfo, void* _context ) {
	lucScalarField*          self          = (lucScalarField*)drawingObject;
	
	_lucScalarFieldCrossSection_Draw( self, window, viewportInfo, _context );
}

void _lucScalarField_CleanUp( void* drawingObject, void* _context ) {
	lucScalarField*          self          = (lucScalarField*)drawingObject;
	
	_lucScalarFieldCrossSection_CleanUp( self, _context );
}

void _lucScalarField_BuildDisplayList( void* drawingObject, void* _context ) {
	lucScalarField*          self          = (lucScalarField*)drawingObject;
	DomainContext*   context       = (DomainContext*) _context;
	FieldVariable*           fieldVariable = self->fieldVariable;
	Coord                    min;
	Coord                    max;
	Dimension_Index          dim_I;

	/* Scale Colour Map */
	FieldVariable_GetMinAndMaxGlobalCoords( fieldVariable, min, max );

	/* Crop the size of the cros-section that you wish to draw */
	for ( dim_I = 0 ; dim_I < fieldVariable->dim ; dim_I++ ) {
		min[ dim_I ] = MAX( self->minCropValues[ dim_I ], min[ dim_I ]);
		max[ dim_I ] = MIN( self->maxCropValues[ dim_I ], max[ dim_I ]);
	}

	
	if (context->dim == 2) {
           if( self->useMesh )
              lucScalarField_DrawWithMesh( self );
           else
              lucScalarFieldCrossSection_DrawCrossSection( self, 0.0, K_AXIS );
	}
	else {
		if ( self->cullFace ) 
			glEnable(GL_CULL_FACE);
	
		glFrontFace(GL_CCW);
		lucScalarFieldCrossSection_DrawCrossSection( self, min[ I_AXIS ], I_AXIS );
		lucScalarFieldCrossSection_DrawCrossSection( self, max[ J_AXIS ], J_AXIS );
		lucScalarFieldCrossSection_DrawCrossSection( self, min[ K_AXIS ], K_AXIS );
	
		glFrontFace(GL_CW);
		lucScalarFieldCrossSection_DrawCrossSection( self, max[ I_AXIS ], I_AXIS );
		lucScalarFieldCrossSection_DrawCrossSection( self, min[ J_AXIS ], J_AXIS );
		lucScalarFieldCrossSection_DrawCrossSection( self, max[ K_AXIS ], K_AXIS );

		glDisable(GL_CULL_FACE);
	}
}

void lucScalarField_DrawWithMesh( lucScalarField* self ) {
   FeVariable* var = (FeVariable*)self->fieldVariable;
   FeMesh* mesh = var->feMesh;
   lucColourMap* cmap = self->colourMap;
   IArray* inc;
   double value;
   int* nodes, nElements, curNode;
   int nodeMap[4] = {0, 1, 3, 2};
   double xi[3], vertex[3];
   int ii, jj, kk;

   lucOpenGLDrawingObject_SyncShadowValues( self, self->fieldVariable );
   glDisable( GL_LIGHTING );
   glNormal3f( 0.0, 0.0, 1.0 );
   glBegin( GL_QUADS );
   nElements = FeMesh_GetElementLocalSize( mesh );
   inc = IArray_New();
   for( ii = 0; ii < nElements; ii++ ) {
      for( jj = 0; jj < 10; jj++ ) {
         for( kk = 0; kk < 10; kk++ ) {

            xi[0] = -1.0 + ((double)kk / 10.0) * 2.0;
            xi[1] = -1.0 + ((double)jj / 10.0) * 2.0;
            FeVariable_InterpolateWithinElement( var, ii, xi, &value );
            lucColourMap_SetOpenGLColourFromValue( cmap, value );
            FeMesh_CoordLocalToGlobal( mesh, ii, xi, vertex );
            glVertex2dv( vertex );

            xi[0] = -1.0 + ((double)(kk + 1) / 10.0) * 2.0;
            xi[1] = -1.0 + ((double)jj / 10.0) * 2.0;
            FeVariable_InterpolateWithinElement( var, ii, xi, &value );
            lucColourMap_SetOpenGLColourFromValue( cmap, value );
            FeMesh_CoordLocalToGlobal( mesh, ii, xi, vertex );
            glVertex2dv( vertex );

            xi[0] = -1.0 + ((double)(kk + 1) / 10.0) * 2.0;
            xi[1] = -1.0 + ((double)(jj + 1) / 10.0) * 2.0;
            FeVariable_InterpolateWithinElement( var, ii, xi, &value );
            lucColourMap_SetOpenGLColourFromValue( cmap, value );
            FeMesh_CoordLocalToGlobal( mesh, ii, xi, vertex );
            glVertex2dv( vertex );

            xi[0] = -1.0 + ((double)kk / 10.0) * 2.0;
            xi[1] = -1.0 + ((double)(jj + 1) / 10.0) * 2.0;
            FeVariable_InterpolateWithinElement( var, ii, xi, &value );
            lucColourMap_SetOpenGLColourFromValue( cmap, value );
            FeMesh_CoordLocalToGlobal( mesh, ii, xi, vertex );
            glVertex2dv( vertex );
         }
      }

#if 0
      FeMesh_GetElementNodes( mesh, ii, inc );
      assert( IArray_GetSize( inc ) == 4 );
      nodes = IArray_GetPtr( inc );
      for( jj = 0; jj < 4; jj++ ) {
         curNode = nodes[nodeMap[jj]];
         FeVariable_GetValueAtNode( var, curNode, &value );
         lucColourMap_SetOpenGLColourFromValue( cmap, value );
         glVertex2dv( Mesh_GetVertex( mesh, curNode ) );
      }
#endif
   }
   NewClass_Delete( inc );
   glEnd();
   glEnable(GL_LIGHTING);
}
