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
** $Id: ScalarField.c 564 2006-05-12 07:36:25Z RobertTurnbull $
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
#include "ScalarFieldOnMeshCrossSection.h"
#include "ScalarFieldOnMesh.h"

#include <assert.h>
#include <gl.h>
#include <glu.h>
#include <string.h>


/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type lucScalarFieldOnMesh_Type = "lucScalarFieldOnMesh";

/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
lucScalarFieldOnMesh* _lucScalarFieldOnMesh_New( 
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
	lucScalarFieldOnMesh*					self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( sizeOfSelf >= sizeof(lucScalarFieldOnMesh) );
	self = (lucScalarFieldOnMesh*) _lucScalarFieldOnMeshCrossSection_New( 
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

void _lucScalarFieldOnMesh_Init( 
		lucScalarFieldOnMesh*                                              self,
		Bool                                                         cullFace )
{
	self->cullFace = cullFace;
}

void _lucScalarFieldOnMesh_Delete( void* drawingObject ) {
	lucScalarFieldOnMesh*  self = (lucScalarFieldOnMesh*)drawingObject;

	_lucScalarFieldOnMeshCrossSection_Delete( self );
}

void _lucScalarFieldOnMesh_Print( void* drawingObject, Stream* stream ) {
	lucScalarFieldOnMesh*  self = (lucScalarFieldOnMesh*)drawingObject;

	_lucScalarFieldOnMeshCrossSection_Print( self, stream );
}

void* _lucScalarFieldOnMesh_Copy( void* drawingObject, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap) {
	lucScalarFieldOnMesh*  self = (lucScalarFieldOnMesh*)drawingObject;
	lucScalarFieldOnMesh* newDrawingObject;

	newDrawingObject = _lucScalarFieldOnMeshCrossSection_Copy( self, dest, deep, nameExt, ptrMap );

	/* TODO */
	abort();

	return (void*) newDrawingObject;
}


void* _lucScalarFieldOnMesh_DefaultNew( Name name ) {
	return (void*) _lucScalarFieldOnMesh_New(
		sizeof(lucScalarFieldOnMesh),
		lucScalarFieldOnMesh_Type,
		_lucScalarFieldOnMesh_Delete,
		_lucScalarFieldOnMesh_Print,
		NULL,
		_lucScalarFieldOnMesh_DefaultNew,
		_lucScalarFieldOnMesh_AssignFromXML,
		_lucScalarFieldOnMesh_Build,
		_lucScalarFieldOnMesh_Initialise,
		_lucScalarFieldOnMesh_Execute,
		_lucScalarFieldOnMesh_Destroy,
		_lucScalarFieldOnMesh_Setup,
		_lucScalarFieldOnMesh_Draw,
		_lucScalarFieldOnMesh_CleanUp,
		_lucScalarFieldOnMesh_BuildDisplayList,
		name );
}

void _lucScalarFieldOnMesh_AssignFromXML( void* drawingObject, Stg_ComponentFactory* cf, void* data ){
	lucScalarFieldOnMesh*  self = (lucScalarFieldOnMesh*)drawingObject;

	/* Construct Parent */
	_lucScalarFieldOnMeshCrossSection_AssignFromXML( self, cf, data );

	_lucScalarFieldOnMesh_Init( 
			self, 
			Stg_ComponentFactory_GetBool( cf, self->name, "cullFace", True ) );
}

void _lucScalarFieldOnMesh_Build( void* drawingObject, void* data ) {
	lucScalarFieldOnMesh*  self = (lucScalarFieldOnMesh*)drawingObject;

	/* Call parent function */
	_lucScalarFieldOnMeshCrossSection_Build( self, data );
}
void _lucScalarFieldOnMesh_Initialise( void* drawingObject, void* data ) {
	lucScalarFieldOnMesh*  self = (lucScalarFieldOnMesh*)drawingObject;

	/* Call parent function */
	_lucScalarFieldOnMeshCrossSection_Initialise( self, data );
}
void _lucScalarFieldOnMesh_Execute( void* drawingObject, void* data ) {}
void _lucScalarFieldOnMesh_Destroy( void* drawingObject, void* data ) {}

void _lucScalarFieldOnMesh_Setup( void* drawingObject, void* _context ) {
	lucScalarFieldOnMesh*          self          = (lucScalarFieldOnMesh*)drawingObject;
	
	_lucScalarFieldOnMeshCrossSection_Setup( self, _context );
}
	
void _lucScalarFieldOnMesh_Draw( void* drawingObject, lucWindow* window, lucViewportInfo* viewportInfo, void* _context ) {
	lucScalarFieldOnMesh*          self          = (lucScalarFieldOnMesh*)drawingObject;
	
	_lucScalarFieldOnMeshCrossSection_Draw( self, window, viewportInfo, _context );
}

void _lucScalarFieldOnMesh_CleanUp( void* drawingObject, void* _context ) {
	lucScalarFieldOnMesh*          self          = (lucScalarFieldOnMesh*)drawingObject;
	
	_lucScalarFieldOnMeshCrossSection_CleanUp( self, _context );
}

void _lucScalarFieldOnMesh_BuildDisplayList( void* drawingObject, void* _context ) {
	lucScalarFieldOnMesh*  self          = (lucScalarFieldOnMesh*)drawingObject;
	FeVariable*            fieldVariable = (FeVariable*) self->fieldVariable;
	Mesh*                  mesh          = (Mesh*) fieldVariable->feMesh;
	Grid*                  vertGrid;

	vertGrid = *(Grid**)ExtensionManager_Get( mesh->info, mesh, self->vertexGridHandle );
	
	if (fieldVariable->dim == 2) {
		lucScalarFieldOnMeshCrossSection_DrawCrossSection( self, 0, K_AXIS );
	}
	else {
		if ( self->cullFace ) 
			glEnable(GL_CULL_FACE);
	
		glFrontFace(GL_CCW);
		lucScalarFieldOnMeshCrossSection_DrawCrossSection( self, 0, I_AXIS );
		lucScalarFieldOnMeshCrossSection_DrawCrossSection( self, vertGrid->sizes[ J_AXIS ] - 1, J_AXIS );
		lucScalarFieldOnMeshCrossSection_DrawCrossSection( self, 0, K_AXIS );
	
		glFrontFace(GL_CW);
		lucScalarFieldOnMeshCrossSection_DrawCrossSection( self, vertGrid->sizes[ I_AXIS ] - 1, I_AXIS );
		lucScalarFieldOnMeshCrossSection_DrawCrossSection( self, 0, J_AXIS );
		lucScalarFieldOnMeshCrossSection_DrawCrossSection( self, vertGrid->sizes[ K_AXIS ] - 1, K_AXIS );

		glDisable(GL_CULL_FACE);
	}
}

