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
** $Id: ScalarField.c 791 2008-09-01 02:09:06Z JulianGiordani $
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
lucScalarField* _lucScalarField_New(  LUCSCALARFIELD_DEFARGS  ) 
{
	lucScalarField*					self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( _sizeOfSelf >= sizeof(lucScalarField) );
	self = (lucScalarField*) _lucScalarFieldCrossSection_New(  LUCSCALARFIELDCROSSSECTION_PASSARGS  );
	
	return self;
}

void _lucScalarField_Init(lucScalarField* self, Bool useMesh)
{
   self->useMesh = useMesh;
}

void _lucScalarField_Delete( void* drawingObject ) {
	lucScalarField*  self = (lucScalarField*)drawingObject;

	_lucScalarFieldCrossSection_Delete( self );
}

void _lucScalarField_Print( void* drawingObject, Stream* stream ) {
	lucScalarField*  self = (lucScalarField*)drawingObject;

	_lucScalarFieldCrossSection_Print( self, stream );
}

void* _lucScalarField_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                                     _sizeOfSelf = sizeof(lucScalarField);
	Type                                                             type = lucScalarField_Type;
	Stg_Class_DeleteFunction*                                     _delete = _lucScalarField_Delete;
	Stg_Class_PrintFunction*                                       _print = _lucScalarField_Print;
	Stg_Class_CopyFunction*                                         _copy = NULL;
	Stg_Component_DefaultConstructorFunction*         _defaultConstructor = _lucScalarField_DefaultNew;
	Stg_Component_ConstructFunction*                           _construct = _lucScalarField_AssignFromXML;
	Stg_Component_BuildFunction*                                   _build = _lucScalarFieldCrossSection_Build;
	Stg_Component_InitialiseFunction*                         _initialise = _lucScalarFieldCrossSection_Initialise;
	Stg_Component_ExecuteFunction*                               _execute = _lucScalarFieldCrossSection_Execute;
	Stg_Component_DestroyFunction*                               _destroy = _lucScalarFieldCrossSection_Destroy;
	lucDrawingObject_SetupFunction*                                _setup = _lucScalarFieldCrossSection_Setup;
	lucDrawingObject_DrawFunction*                                  _draw = _lucOpenGLDrawingObject_Draw;
	lucDrawingObject_CleanUpFunction*                            _cleanUp = _lucOpenGLDrawingObject_CleanUp;
	lucOpenGLDrawingObject_BuildDisplayListFunction*    _buildDisplayList = _lucScalarField_BuildDisplayList;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

	return (void*) _lucScalarField_New(  LUCSCALARFIELD_PASSARGS  );
}

void _lucScalarField_AssignFromXML( void* drawingObject, Stg_ComponentFactory* cf, void* data ){
	lucScalarField*  self = (lucScalarField*)drawingObject;

	/* Construct Parent */
	_lucScalarFieldCrossSection_AssignFromXML( self, cf, data );

	_lucScalarField_Init(self, Stg_ComponentFactory_GetBool( cf, self->name, "useMesh", False ));
}

void _lucScalarField_BuildDisplayList( void* drawingObject, void* _context ) {
	lucScalarField*          self          = (lucScalarField*)drawingObject;
	DomainContext*   context       = (DomainContext*) _context;

	if (context->dim == 2) 
   {
       if( self->useMesh )
          lucScalarField_DrawWithMesh( self );
       else {
          lucScalarFieldCrossSection_DrawCrossSection( lucCrossSection_Set(self, 0.0, K_AXIS, False), GL_CCW);
       }
	}
	else 
   {
		/* Cross sections at minimums, default winding for faces */
		lucScalarFieldCrossSection_DrawCrossSection( lucCrossSection_Set(self, 0.0, I_AXIS, True), GL_CCW);
		lucScalarFieldCrossSection_DrawCrossSection( lucCrossSection_Set(self, 0.0, J_AXIS, True), GL_CCW);
		lucScalarFieldCrossSection_DrawCrossSection( lucCrossSection_Set(self, 0.0, K_AXIS, True), GL_CCW);

		/* Cross sections at maximums, reverse winding for faces */
		lucScalarFieldCrossSection_DrawCrossSection( lucCrossSection_Set(self, 1.0, I_AXIS, True), GL_CW);
		lucScalarFieldCrossSection_DrawCrossSection( lucCrossSection_Set(self, 1.0, J_AXIS, True), GL_CW);
		lucScalarFieldCrossSection_DrawCrossSection( lucCrossSection_Set(self, 1.0, K_AXIS, True), GL_CW);
	}
}

void lucScalarField_DrawWithMesh( lucScalarField* self ) {
   FeVariable* var = (FeVariable*)self->fieldVariable;
   FeMesh* mesh = var->feMesh;
   lucColourMap* cmap = self->colourMap;
   IArray* inc;
   double value;
   unsigned nElements;
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


