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
lucScalarFieldOnMesh* _lucScalarFieldOnMesh_New(  LUCSCALARFIELDONMESH_DEFARGS  ) 
{
	lucScalarFieldOnMesh*					self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( _sizeOfSelf >= sizeof(lucScalarFieldOnMesh) );
	self = (lucScalarFieldOnMesh*) _lucScalarFieldOnMeshCrossSection_New(  LUCSCALARFIELDONMESHCROSSSECTION_PASSARGS  );
	
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

void* _lucScalarFieldOnMesh_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                                     _sizeOfSelf = sizeof(lucScalarFieldOnMesh);
	Type                                                             type = lucScalarFieldOnMesh_Type;
	Stg_Class_DeleteFunction*                                     _delete = _lucScalarFieldOnMesh_Delete;
	Stg_Class_PrintFunction*                                       _print = _lucScalarFieldOnMesh_Print;
	Stg_Class_CopyFunction*                                         _copy = NULL;
	Stg_Component_DefaultConstructorFunction*         _defaultConstructor = _lucScalarFieldOnMesh_DefaultNew;
	Stg_Component_ConstructFunction*                           _construct = _lucScalarFieldOnMesh_AssignFromXML;
	Stg_Component_BuildFunction*                                   _build = _lucScalarFieldOnMesh_Build;
	Stg_Component_InitialiseFunction*                         _initialise = _lucScalarFieldOnMesh_Initialise;
	Stg_Component_ExecuteFunction*                               _execute = _lucScalarFieldOnMesh_Execute;
	Stg_Component_DestroyFunction*                               _destroy = _lucScalarFieldOnMesh_Destroy;
	lucDrawingObject_SetupFunction*                                _setup = _lucScalarFieldOnMeshCrossSection_Setup;
	lucDrawingObject_DrawFunction*                                  _draw = _lucOpenGLDrawingObject_Draw;
	lucDrawingObject_CleanUpFunction*                            _cleanUp = _lucOpenGLDrawingObject_CleanUp;
	lucOpenGLDrawingObject_BuildDisplayListFunction*    _buildDisplayList = _lucScalarFieldOnMesh_BuildDisplayList;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

	return (void*) _lucScalarFieldOnMesh_New(  LUCSCALARFIELDONMESH_PASSARGS  );
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

void _lucScalarFieldOnMesh_BuildDisplayList( void* drawingObject, void* _context ) {
	lucScalarFieldOnMesh*  self          = (lucScalarFieldOnMesh*)drawingObject;
	FeVariable*            fieldVariable = (FeVariable*) self->fieldVariable;
	Mesh*                  mesh          = (Mesh*) fieldVariable->feMesh;
	Grid*                  vertGrid;

	vertGrid = *(Grid**)ExtensionManager_Get( mesh->info, mesh, self->vertexGridHandle );
	
	if (fieldVariable->dim == 2) {
		lucScalarFieldOnMeshCrossSection_DrawCrossSection( lucCrossSection_Set(self, 0.0, K_AXIS, False));
	}
	else {
	   glEnable(GL_LIGHTING);

		if ( self->cullFace ) 
			glEnable(GL_CULL_FACE);
	
		glFrontFace(GL_CCW);
		lucScalarFieldOnMeshCrossSection_DrawCrossSection( lucCrossSection_Set(self, 0.0, I_AXIS, False));
		lucScalarFieldOnMeshCrossSection_DrawCrossSection( lucCrossSection_Set(self, vertGrid->sizes[ J_AXIS ] - 1, J_AXIS, False));
		lucScalarFieldOnMeshCrossSection_DrawCrossSection( lucCrossSection_Set(self, 0.0, K_AXIS, False));

		glFrontFace(GL_CW);
		lucScalarFieldOnMeshCrossSection_DrawCrossSection( lucCrossSection_Set(self, vertGrid->sizes[ I_AXIS ] - 1, I_AXIS, False));
		lucScalarFieldOnMeshCrossSection_DrawCrossSection( lucCrossSection_Set(self, 0.0, J_AXIS, False));
		lucScalarFieldOnMeshCrossSection_DrawCrossSection( lucCrossSection_Set(self, vertGrid->sizes[ K_AXIS ] - 1, K_AXIS, False));


		glFrontFace(GL_CCW);
		//glDisable(GL_CULL_FACE);
	}
}



