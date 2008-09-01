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
#ifdef HAVE_OPENGL_FRAMEWORK
	#include <OpenGL/gl.h>
	#include <OpenGL/glu.h>
#else
	#include <gl.h>
	#include <glu.h>
#endif
#include <string.h>

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type lucVectorArrows_Type = "lucVectorArrows";

/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
lucVectorArrows* _lucVectorArrows_New( 
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
	lucVectorArrows*					self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( sizeOfSelf >= sizeof(lucVectorArrows) );
	self = (lucVectorArrows*) _lucVectorArrowCrossSection_New( 
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

void* _lucVectorArrows_Copy( void* drawingObject, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap) {
	lucVectorArrows*  self = (lucVectorArrows*)drawingObject;
	lucVectorArrows* newDrawingObject;

	newDrawingObject = _lucVectorArrowCrossSection_Copy( self, dest, deep, nameExt, ptrMap );

	/* TODO */
	abort();

	return (void*) newDrawingObject;
}


void* _lucVectorArrows_DefaultNew( Name name ) {
	return (void*) _lucVectorArrows_New(
		sizeof(lucVectorArrows),
		lucVectorArrows_Type,
		_lucVectorArrows_Delete,
		_lucVectorArrows_Print,
		NULL,
		_lucVectorArrows_DefaultNew,
		_lucVectorArrows_Construct,
		_lucVectorArrows_Build,
		_lucVectorArrows_Initialise,
		_lucVectorArrows_Execute,
		_lucVectorArrows_Destroy,
		_lucVectorArrows_Setup,
		_lucVectorArrows_Draw,
		_lucVectorArrows_CleanUp,
		_lucVectorArrows_BuildDisplayList,
		name );
}

void _lucVectorArrows_Construct( void* drawingObject, Stg_ComponentFactory* cf, void* data ){
	lucVectorArrows* self = (lucVectorArrows*)drawingObject;

	/* Construct Parent */
	_lucVectorArrowCrossSection_Construct( self, cf, data );
	
	_lucVectorArrows_Init( self );
}

void _lucVectorArrows_Build( void* drawingObject, void* data ) {}
void _lucVectorArrows_Initialise( void* drawingObject, void* data ) {}
void _lucVectorArrows_Execute( void* drawingObject, void* data ) {}
void _lucVectorArrows_Destroy( void* drawingObject, void* data ) {}

void _lucVectorArrows_Setup( void* drawingObject, void* _context ) {
	lucVectorArrows*       self            = (lucVectorArrows*)drawingObject;
	
	_lucVectorArrowCrossSection_Setup( self, _context );
}

void _lucVectorArrows_Draw( void* drawingObject, lucWindow* window, lucViewportInfo* viewportInfo, void* _context ) {
	lucVectorArrows*       self            = (lucVectorArrows*)drawingObject;

	_lucVectorArrowCrossSection_Draw( self, window, viewportInfo, _context );
}

void _lucVectorArrows_CleanUp( void* drawingObject, void* _context ) {
	lucVectorArrows*       self            = (lucVectorArrows*)drawingObject;
	
	_lucVectorArrowCrossSection_CleanUp( self, _context );
}
	
void _lucVectorArrows_BuildDisplayList( void* drawingObject, void* _context ) {
	lucVectorArrows*       self            = (lucVectorArrows*)drawingObject;
	DomainContext* context         = (DomainContext*) _context;
	Dimension_Index        dim             = context->dim;

	if ( dim == 2 ) {
		_lucVectorArrowCrossSection_DrawCrossSection( self, dim, 0.0, K_AXIS );
	}
	else {
		Coord             globalMin;
		Coord             globalMax;
		double            dz;
		double            depth;

		FieldVariable_GetMinAndMaxGlobalCoords( self->vectorVariable, globalMin, globalMax );
	
		dz = (globalMax[K_AXIS] - globalMin[K_AXIS])/(double)self->resolution[ K_AXIS ];

		for ( depth = globalMin[ K_AXIS ] + dz * 0.5 ; depth < globalMax[ K_AXIS ] ; depth += dz) {
			_lucVectorArrowCrossSection_DrawCrossSection( self, dim, depth, K_AXIS );
		}
	}
}
