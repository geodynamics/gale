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
** $Id: FieldVariableBorder.c 791 2008-09-01 02:09:06Z JulianGiordani $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>

#include <glucifer/Base/Base.h>
#include <glucifer/RenderingEngines/RenderingEngines.h>

#include "types.h"
#include "OpenGLDrawingObject.h"
#include "FieldVariableBorder.h"

#include <assert.h>
#include <gl.h>
#include <glu.h>
#include <string.h>

#ifndef MASTER
	#define MASTER 0
#endif

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type lucFieldVariableBorder_Type = "lucFieldVariableBorder";

/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
lucFieldVariableBorder* _lucFieldVariableBorder_New(  LUCFIELDVARIABLEBORDER_DEFARGS  ) 
{
	lucFieldVariableBorder*					self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( _sizeOfSelf >= sizeof(lucFieldVariableBorder) );
	self = (lucFieldVariableBorder*) _lucOpenGLDrawingObject_New(  LUCOPENGLDRAWINGOBJECT_PASSARGS  );
	
	return self;
}

void _lucFieldVariableBorder_Init( 
		lucFieldVariableBorder*                                      self,
		FieldVariable*                                               fieldVariable,
		Name                                                         colourName,
		float                                                        lineWidth ) 
{
	self->fieldVariable = fieldVariable;
	
	lucColour_FromString( &self->colour, colourName );
	self->lineWidth = lineWidth;
}

void _lucFieldVariableBorder_Delete( void* drawingObject ) {
	lucFieldVariableBorder*  self = (lucFieldVariableBorder*)drawingObject;

	_lucOpenGLDrawingObject_Delete( self );
}

void _lucFieldVariableBorder_Print( void* drawingObject, Stream* stream ) {
	lucFieldVariableBorder*  self = (lucFieldVariableBorder*)drawingObject;

	_lucOpenGLDrawingObject_Print( self, stream );
}

void* _lucFieldVariableBorder_Copy( void* drawingObject, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap) {
	lucFieldVariableBorder*  self = (lucFieldVariableBorder*)drawingObject;
	lucFieldVariableBorder* newDrawingObject;

	newDrawingObject = _lucOpenGLDrawingObject_Copy( self, dest, deep, nameExt, ptrMap );

	/* TODO */
	abort();

	return (void*) newDrawingObject;
}


void* _lucFieldVariableBorder_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                                     _sizeOfSelf = sizeof(lucFieldVariableBorder);
	Type                                                             type = lucFieldVariableBorder_Type;
	Stg_Class_DeleteFunction*                                     _delete = _lucFieldVariableBorder_Delete;
	Stg_Class_PrintFunction*                                       _print = _lucFieldVariableBorder_Print;
	Stg_Class_CopyFunction*                                         _copy = NULL;
	Stg_Component_DefaultConstructorFunction*         _defaultConstructor = _lucFieldVariableBorder_DefaultNew;
	Stg_Component_ConstructFunction*                           _construct = _lucFieldVariableBorder_AssignFromXML;
	Stg_Component_BuildFunction*                                   _build = _lucFieldVariableBorder_Build;
	Stg_Component_InitialiseFunction*                         _initialise = _lucFieldVariableBorder_Initialise;
	Stg_Component_ExecuteFunction*                               _execute = _lucFieldVariableBorder_Execute;
	Stg_Component_DestroyFunction*                               _destroy = _lucFieldVariableBorder_Destroy;
	lucDrawingObject_SetupFunction*                                _setup = _lucFieldVariableBorder_Setup;
	lucDrawingObject_DrawFunction*                                  _draw = _lucFieldVariableBorder_Draw;
	lucDrawingObject_CleanUpFunction*                            _cleanUp = _lucFieldVariableBorder_CleanUp;
	lucOpenGLDrawingObject_BuildDisplayListFunction*    _buildDisplayList = _lucFieldVariableBorder_BuildDisplayList;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

	return (void*) _lucFieldVariableBorder_New(  LUCFIELDVARIABLEBORDER_PASSARGS  );
}

void _lucFieldVariableBorder_AssignFromXML( void* drawingObject, Stg_ComponentFactory* cf, void* data ){
	lucFieldVariableBorder*  self = (lucFieldVariableBorder*)drawingObject;
	FieldVariable*           fieldVariable;

	/* Construct Parent */
	_lucOpenGLDrawingObject_AssignFromXML( self, cf, data );

	fieldVariable =  Stg_ComponentFactory_ConstructByKey( cf, self->name, "FieldVariable", FieldVariable, True, data );
	
	_lucFieldVariableBorder_Init( 
			self, 
			fieldVariable, 
			Stg_ComponentFactory_GetString( cf, self->name, "colour", "Black" ),
			(float) Stg_ComponentFactory_GetDouble( cf, self->name, "lineWidth", 0.5 ) );

}

void _lucFieldVariableBorder_Build( void* drawingObject, void* data ) {}
void _lucFieldVariableBorder_Initialise( void* drawingObject, void* data ) {}
void _lucFieldVariableBorder_Execute( void* drawingObject, void* data ) {}
void _lucFieldVariableBorder_Destroy( void* drawingObject, void* data ) {}

void _lucFieldVariableBorder_Setup( void* drawingObject, void* _context ) {
	lucFieldVariableBorder*       self            = (lucFieldVariableBorder*)drawingObject;
	_lucOpenGLDrawingObject_Setup( self, _context );
}
	
void _lucFieldVariableBorder_Draw( void* drawingObject, lucWindow* window, lucViewportInfo* viewportInfo, void* _context ) {
	lucFieldVariableBorder*       self            = (lucFieldVariableBorder*)drawingObject;
	_lucOpenGLDrawingObject_Draw( self, window, viewportInfo, _context );
}


void _lucFieldVariableBorder_CleanUp( void* drawingObject, void* _context ) {
	lucFieldVariableBorder*       self            = (lucFieldVariableBorder*)drawingObject;
	_lucOpenGLDrawingObject_CleanUp( self, _context );
}	

void _lucFieldVariableBorder_BuildDisplayList( void* drawingObject, void* _context ) {
	lucFieldVariableBorder*  self    = (lucFieldVariableBorder*)drawingObject;
	DomainContext*   context = (DomainContext*) _context;
	Dimension_Index          dim     = context->dim;
	FieldVariable*           field   = self->fieldVariable;
	Coord                    min;
	Coord                    max;

	/* THis needs to be first so all processors call it - collective comms */
	/* Grab values from field variable */
	FieldVariable_GetMinAndMaxGlobalCoords( field, min, max );

	/* Initialise OpenGL stuff */
	glDisable(GL_LIGHTING);

	/* Line smoothing will not work correctly with depth testing enabled*/
	glDisable(GL_DEPTH_TEST);
	glEnable(GL_LINE_SMOOTH);
	glHint (GL_LINE_SMOOTH_HINT, GL_NICEST);
	glEnable(GL_BLEND);
	glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glPolygonMode ( GL_FRONT_AND_BACK, GL_LINE ) ;

	lucColour_SetOpenGLColour( &self->colour );
	glLineWidth( self->lineWidth );

	if (dim == 2) 
		glRectd( min[ I_AXIS ], min[ J_AXIS ], max[ I_AXIS ], max[ J_AXIS ]);
	else if ( dim == 3 ) {
		glBegin( GL_QUAD_STRIP );
			glVertex3d( min[ I_AXIS ], min[ J_AXIS ], min[ K_AXIS ]);
			glVertex3d( min[ I_AXIS ], max[ J_AXIS ], min[ K_AXIS ]);

			glVertex3d( max[ I_AXIS ], min[ J_AXIS ], min[ K_AXIS ]);
			glVertex3d( max[ I_AXIS ], max[ J_AXIS ], min[ K_AXIS ]);

			glVertex3d( max[ I_AXIS ], min[ J_AXIS ], max[ K_AXIS ]);
			glVertex3d( max[ I_AXIS ], max[ J_AXIS ], max[ K_AXIS ]);

			glVertex3d( min[ I_AXIS ], min[ J_AXIS ], max[ K_AXIS ]);
			glVertex3d( min[ I_AXIS ], max[ J_AXIS ], max[ K_AXIS ]);
			
			glVertex3d( min[ I_AXIS ], min[ J_AXIS ], min[ K_AXIS ]);
			glVertex3d( min[ I_AXIS ], max[ J_AXIS ], min[ K_AXIS ]);
		glEnd();
	}
	
	/* Clean up OpenGL stuff */
 	glEnable(GL_DEPTH_TEST);
	glPolygonMode ( GL_FRONT_AND_BACK, GL_FILL ) ;
	glEnable( GL_LIGHTING );
}


