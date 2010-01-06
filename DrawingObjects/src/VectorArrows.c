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
** $Id: VectorArrows.c 791 2008-09-01 02:09:06Z JulianGiordani $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>

#include <glucifer/Base/Base.h>
#include <glucifer/RenderingEngines/RenderingEngines.h>

#include "types.h"
#include "OpenGLDrawingObject.h"
#include "VectorArrowCrossSection.h"
#include "VectorArrows.h"

#include <assert.h>
#include <gl.h>
#include <glu.h>
#include <string.h>

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type lucVectorArrows_Type = "lucVectorArrows";

/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
lucVectorArrows* _lucVectorArrows_New(  LUCVECTORARROWS_DEFARGS  ) 
{
	lucVectorArrows*					self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( _sizeOfSelf >= sizeof(lucVectorArrows) );
	self = (lucVectorArrows*) _lucVectorArrowCrossSection_New(  LUCVECTORARROWCROSSSECTION_PASSARGS  );
	
	return self;
}

void _lucVectorArrows_Init( lucVectorArrows*                                             self ) {
}

void _lucVectorArrows_Delete( void* drawingObject ) {
	lucVectorArrows*  self = (lucVectorArrows*)drawingObject;

	_lucVectorArrowCrossSection_Delete( self );
}

void _lucVectorArrows_Print( void* drawingObject, Stream* stream ) {
	lucVectorArrows*  self = (lucVectorArrows*)drawingObject;

	_lucVectorArrowCrossSection_Print( self, stream );
}

void* _lucVectorArrows_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                                     _sizeOfSelf = sizeof(lucVectorArrows);
	Type                                                             type = lucVectorArrows_Type;
	Stg_Class_DeleteFunction*                                     _delete = _lucVectorArrows_Delete;
	Stg_Class_PrintFunction*                                       _print = _lucVectorArrows_Print;
	Stg_Class_CopyFunction*                                         _copy = NULL;
	Stg_Component_DefaultConstructorFunction*         _defaultConstructor = _lucVectorArrows_DefaultNew;
	Stg_Component_ConstructFunction*                           _construct = _lucVectorArrows_AssignFromXML;
	Stg_Component_BuildFunction*                                   _build = _lucVectorArrows_Build;
	Stg_Component_InitialiseFunction*                         _initialise = _lucVectorArrows_Initialise;
	Stg_Component_ExecuteFunction*                               _execute = _lucVectorArrows_Execute;
	Stg_Component_DestroyFunction*                               _destroy = _lucVectorArrows_Destroy;
	lucDrawingObject_SetupFunction*                                _setup = _lucOpenGLDrawingObject_Setup;
	lucDrawingObject_DrawFunction*                                  _draw = _lucOpenGLDrawingObject_Draw;
	lucDrawingObject_CleanUpFunction*                            _cleanUp = _lucOpenGLDrawingObject_CleanUp;
	lucOpenGLDrawingObject_BuildDisplayListFunction*    _buildDisplayList = _lucVectorArrows_BuildDisplayList;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

	return (void*) _lucVectorArrows_New(  LUCVECTORARROWS_PASSARGS  );
}

void _lucVectorArrows_AssignFromXML( void* drawingObject, Stg_ComponentFactory* cf, void* data ){
	lucVectorArrows* self = (lucVectorArrows*)drawingObject;

	/* Construct Parent */
	_lucVectorArrowCrossSection_AssignFromXML( self, cf, data );
	
	_lucVectorArrows_Init( self );
}

void _lucVectorArrows_Build( void* drawingObject, void* data ) {}
void _lucVectorArrows_Initialise( void* drawingObject, void* data ) {}
void _lucVectorArrows_Execute( void* drawingObject, void* data ) {}
void _lucVectorArrows_Destroy( void* drawingObject, void* data ) {}

void _lucVectorArrows_BuildDisplayList( void* drawingObject, void* _context ) {
	lucVectorArrows*       self            = (lucVectorArrows*)drawingObject;
	DomainContext* context         = (DomainContext*) _context;
	Dimension_Index        dim             = context->dim;

	if ( dim == 2 ) 
   {
      _lucVectorArrowCrossSection_DrawCrossSection( lucCrossSection_Set(self, 0.0, K_AXIS, False), dim);
	}
	else 
   {
		double dz = 1/(double)self->resolution[ K_AXIS ];
      self->axis = K_AXIS;
      self->interpolate = True;
		for ( self->value = 0.0 ; self->value < 1.0+dz ; self->value += dz) {
		   _lucVectorArrowCrossSection_DrawCrossSection( self, dim);
		}
	}
}


