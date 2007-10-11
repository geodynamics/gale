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
** $Id: Eigenvectors.c 740 2007-10-11 08:05:31Z SteveQuenette $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>

#include <glucifer/Base/Base.h>
#include <glucifer/RenderingEngines/RenderingEngines.h>

#include "types.h"
#include "OpenGLDrawingObject.h"
#include "EigenvectorsCrossSection.h"
#include "Eigenvectors.h"

#include <assert.h>
#include <gl.h>
#include <glu.h>
#include <string.h>

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type lucEigenvectors_Type = "lucEigenvectors";

/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
lucEigenvectors* _lucEigenvectors_New( 
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
	lucEigenvectors*					self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( sizeOfSelf >= sizeof(lucEigenvectors) );
	self = (lucEigenvectors*) _lucEigenvectorsCrossSection_New( 
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

void _lucEigenvectors_Init( lucEigenvectors*                                             self ) {
}

void _lucEigenvectors_Delete( void* drawingObject ) {
	lucEigenvectors*  self = (lucEigenvectors*)drawingObject;

	_lucEigenvectorsCrossSection_Delete( self );
}

void _lucEigenvectors_Print( void* drawingObject, Stream* stream ) {
	lucEigenvectors*  self = (lucEigenvectors*)drawingObject;

	_lucEigenvectorsCrossSection_Print( self, stream );
}

void* _lucEigenvectors_Copy( void* drawingObject, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap) {
	lucEigenvectors*  self = (lucEigenvectors*)drawingObject;
	lucEigenvectors* newDrawingObject;

	newDrawingObject = _lucEigenvectorsCrossSection_Copy( self, dest, deep, nameExt, ptrMap );

	/* TODO */
	abort();

	return (void*) newDrawingObject;
}


void* _lucEigenvectors_DefaultNew( Name name ) {
	return (void*) _lucEigenvectors_New(
		sizeof(lucEigenvectors),
		lucEigenvectors_Type,
		_lucEigenvectors_Delete,
		_lucEigenvectors_Print,
		NULL,
		_lucEigenvectors_DefaultNew,
		_lucEigenvectors_Construct,
		_lucEigenvectors_Build,
		_lucEigenvectors_Initialise,
		_lucEigenvectors_Execute,
		_lucEigenvectors_Destroy,
		_lucEigenvectors_Setup,
		_lucEigenvectors_Draw,
		_lucEigenvectors_CleanUp,
		_lucEigenvectors_BuildDisplayList,
		name );
}

void _lucEigenvectors_Construct( void* drawingObject, Stg_ComponentFactory* cf, void* data ){
	lucEigenvectors* self = (lucEigenvectors*)drawingObject;

	/* Construct Parent */
	_lucEigenvectorsCrossSection_Construct( self, cf, data );
	
	_lucEigenvectors_Init( self );
}

void _lucEigenvectors_Build( void* drawingObject, void* data ) {}
void _lucEigenvectors_Initialise( void* drawingObject, void* data ) {}
void _lucEigenvectors_Execute( void* drawingObject, void* data ) {}
void _lucEigenvectors_Destroy( void* drawingObject, void* data ) {}

void _lucEigenvectors_Setup( void* drawingObject, void* _context ) {
	lucEigenvectors*       self            = (lucEigenvectors*)drawingObject;
	
	_lucEigenvectorsCrossSection_Setup( self, _context );
}

void _lucEigenvectors_Draw( void* drawingObject, lucWindow* window, lucViewportInfo* viewportInfo, void* _context ) {
	lucEigenvectors*       self            = (lucEigenvectors*)drawingObject;

	_lucEigenvectorsCrossSection_Draw( self, window, viewportInfo, _context );
}

void _lucEigenvectors_CleanUp( void* drawingObject, void* _context ) {
	lucEigenvectors*       self            = (lucEigenvectors*)drawingObject;
	
	_lucEigenvectorsCrossSection_CleanUp( self, _context );
}
	
void _lucEigenvectors_BuildDisplayList( void* drawingObject, void* _context ) {
	lucEigenvectors*       self            = (lucEigenvectors*)drawingObject;
	DomainContext* context         = (DomainContext*) _context;
	Dimension_Index        dim             = context->dim;

	if ( dim == 2 ) {
		_lucEigenvectorsCrossSection_DrawCrossSection( self, dim, 0.0, K_AXIS );
	}
	else {
		Coord             globalMin;
		Coord             globalMax;
		double            dz;
		double            depth;

		FieldVariable_GetMinAndMaxGlobalCoords( self->tensorField, globalMin, globalMax );
	
		dz = (globalMax[K_AXIS] - globalMin[K_AXIS])/(double)self->resolution[ K_AXIS ];

		for ( depth = globalMin[ K_AXIS ] + dz * 0.5 ; depth < globalMax[ K_AXIS ] ; depth += dz) {
			_lucEigenvectorsCrossSection_DrawCrossSection( self, dim, depth, K_AXIS );
		}
	}
}
