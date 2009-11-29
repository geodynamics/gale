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
** $Id: Axis.c 510 2006-02-17 04:33:32Z RobertTurnbull $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>

#include <glucifer/Base/Base.h>
#include <glucifer/RenderingEngines/RenderingEngines.h>

#ifdef HAVE_GL2PS
	#include <gl2ps.h>
#endif

#include "types.h"
#include "OpenGLDrawingObject.h"
#include "Axis.h"

#include <assert.h>
#include <gl.h>
#include <glu.h>
#include <string.h>

#ifndef MASTER
	#define MASTER 0
#endif

const Type lucAxis_Type = "lucAxis";

lucAxis* lucAxis_New( 
		Name                                                  name,
		Coord                                                 origin,
    		float 				                      length,
		lucColour                                             colourX,
		lucColour                                             colourY,
		lucColour                                             colourZ)
{
	lucAxis* self = (lucAxis*) _lucAxis_DefaultNew( name );

	lucAxis_InitAll( self, origin, length, colourX, colourY, colourZ);

	return self;
}

lucAxis* _lucAxis_New(  LUCAXIS_DEFARGS  )
{
	lucAxis*    self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( _sizeOfSelf >= sizeof(lucAxis) );
	self = (lucAxis*)  _lucOpenGLDrawingObject_New(  LUCOPENGLDRAWINGOBJECT_PASSARGS  );
	
	
	return self;
}

void lucAxis_Init(		
		lucAxis*                                         self,
		Coord                                            origin,
		float 				                 length,
		lucColour                                        colourX,
		lucColour                                        colourY,
		lucColour                                        colourZ) 
{
	
	self->length = length;
	memcpy( self->origin, origin, sizeof(Coord) );	
	memcpy( &(self->colourX), &colourX, sizeof(lucColour) );	
	memcpy( &(self->colourY), &colourY, sizeof(lucColour) );	
	memcpy( &(self->colourZ), &colourZ, sizeof(lucColour) );	

	
}

void lucAxis_InitAll( 
		void*                                              axis,
		Coord                                              origin,
	        float 				                   length,
		lucColour                                          colourX,
		lucColour                                          colourY,
		lucColour                                          colourZ)
{
	lucAxis* self        = axis;

	/* TODO Init parent */
	lucAxis_Init( self, origin, length, colourX, colourY, colourZ );
}

void _lucAxis_Delete( void* drawingObject ) {
	lucAxis*  self = (lucAxis*)drawingObject;

	_lucOpenGLDrawingObject_Delete( self );
}

void _lucAxis_Print( void* drawingObject, Stream* stream ) {
	lucAxis*  self = (lucAxis*)drawingObject;

	_lucOpenGLDrawingObject_Print( self, stream );
}

void* _lucAxis_Copy( void* axis, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	lucAxis* self        = axis;
	lucAxis* newAxis;

	newAxis = _Stg_Component_Copy( self, dest, deep, nameExt, ptrMap );

	newAxis->length = self->length;
	memcpy( &(newAxis->colourX),       &(self->colourX),       sizeof(lucColour) );
	memcpy( &(newAxis->colourY),       &(self->colourY),       sizeof(lucColour) );
	memcpy( &(newAxis->colourZ),       &(self->colourZ),       sizeof(lucColour) );
	memcpy( newAxis->origin,       self->origin,       sizeof(Coord) );
	
	return (void*) newAxis;
}

void* _lucAxis_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                                     _sizeOfSelf = sizeof( lucAxis );
	Type                                                             type = lucAxis_Type;
	Stg_Class_DeleteFunction*                                     _delete = _lucAxis_Delete;
	Stg_Class_PrintFunction*                                       _print = _lucAxis_Print;
	Stg_Class_CopyFunction*                                         _copy = _lucAxis_Copy;
	Stg_Component_DefaultConstructorFunction*         _defaultConstructor = _lucAxis_DefaultNew;
	Stg_Component_ConstructFunction*                           _construct = _lucAxis_AssignFromXML;
	Stg_Component_BuildFunction*                                   _build = _lucAxis_Build;
	Stg_Component_InitialiseFunction*                         _initialise = _lucAxis_Initialise;
	Stg_Component_ExecuteFunction*                               _execute = _lucAxis_Execute;
	Stg_Component_DestroyFunction*                               _destroy = _lucAxis_Destroy;
	lucDrawingObject_SetupFunction*                                _setup = _lucAxis_Setup;
	lucDrawingObject_DrawFunction*                                  _draw = _lucAxis_Draw;
	lucDrawingObject_CleanUpFunction*                            _cleanUp = _lucAxis_CleanUp;
	lucOpenGLDrawingObject_BuildDisplayListFunction*    _buildDisplayList = _lucAxis_BuildDisplayList;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = ZERO;

	return _lucAxis_New(  LUCAXIS_PASSARGS  );
}

void _lucAxis_AssignFromXML( void* axis, Stg_ComponentFactory* cf, void* data ) {
	lucAxis*             self               = (lucAxis*) axis;
        Name colourNameX;
	Name colourNameY;	
	Name colourNameZ;		
	
	Coord origin;
	
	/* Get Stereo Type */
         /* Construct Parent */
	_lucDrawingObject_AssignFromXML( self, cf, data );

	colourNameX  = Stg_ComponentFactory_GetString( cf, self->name, "colourX", "Red") ;
	colourNameY  = Stg_ComponentFactory_GetString( cf, self->name, "colourY", "Green") ;
	colourNameZ  = Stg_ComponentFactory_GetString( cf, self->name, "colourZ", "Blue") ;
	
	lucColour_FromString( &self->colourX, colourNameX );	
	lucColour_FromString( &self->colourY, colourNameY );
	lucColour_FromString( &self->colourZ, colourNameZ );
	
	origin[I_AXIS]  = Stg_ComponentFactory_GetDouble( cf, self->name, "originX", -0.05 );
	origin[J_AXIS]  = Stg_ComponentFactory_GetDouble( cf, self->name, "originY", -0.05 );
	origin[K_AXIS]  = Stg_ComponentFactory_GetDouble( cf, self->name, "originZ", -0.05 );
	
       	lucAxis_InitAll( self, 
	                origin,
			Stg_ComponentFactory_GetDouble( cf, self->name, "length", 0.2 ),
		        self->colourX,
			self->colourY,
			self->colourZ);
			
}

void _lucAxis_Build( void* Axis, void* data ) { }
void _lucAxis_Initialise( void* Axis, void* data ) { }
void _lucAxis_Execute( void* Axis, void* data ) { }
void _lucAxis_Destroy( void* Axis, void* data ) { }

void _lucAxis_Setup( void* drawingObject, void* _context ) {
	lucAxis*       self            = (lucAxis*)drawingObject;
	_lucOpenGLDrawingObject_Setup( self, _context );
}
void _lucAxis_Draw( void* drawingObject, lucWindow* window, lucViewportInfo* viewportInfo, void* _context ) {
	lucAxis*         self     = (lucAxis*)drawingObject;
        lucViewport*     viewport = viewportInfo->viewport;
	DomainContext*   context  = (DomainContext*) _context;
	Dimension_Index  dim      = context->dim;
        double rodLength          = 0.0;
	double arrowHeadLength    = 0.0;
	double textSpacing        = 0.0;
		
	/* Initialise OpenGL stuff */
	glShadeModel(GL_SMOOTH);
	glDisable(GL_LIGHTING);

	glEnable(GL_LINE_SMOOTH);
	glEnable(GL_BLEND);	
	
	/* Disable lighting because we don't want a 3D effect */
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

	/* The rodlength is the total length of the arrow line.
	   By default it is 0.25 */
	rodLength = self->length;
	 
	/* The tip of the arrow hea starts at rodLength. 
	   The size of the arrow is a fifth of the total length */
	arrowHeadLength = rodLength/5.0;

	textSpacing = 0; //arrowHeadLength;
	 
	if (dim == 2) {
	        /* Drawing the X axis, default is the RED color */
		lucColour_SetOpenGLColour( &self->colourX );

		glBegin( GL_LINES );
			glVertex2f( self->origin[I_AXIS], self->origin[J_AXIS] ); 
			glVertex2f( self->origin[I_AXIS] + rodLength , self->origin[J_AXIS]  );
		glEnd(); 
		glBegin(GL_TRIANGLES);
			glVertex2f( self->origin[I_AXIS] + rodLength, self->origin[J_AXIS] );
			glVertex2f( self->origin[I_AXIS] + rodLength - arrowHeadLength, self->origin[J_AXIS] - arrowHeadLength/2.0);
			glVertex2f( self->origin[I_AXIS] + rodLength - arrowHeadLength, self->origin[J_AXIS] + arrowHeadLength/2.0);
		glEnd();
		lucPrint(self->origin[I_AXIS] + rodLength + textSpacing, self->origin[J_AXIS], "X");
		
		/* Drawing the Y axis, default is the GREEN color */
		lucColour_SetOpenGLColour( &self->colourY );

		glBegin( GL_LINES );
			glVertex2f( self->origin[I_AXIS], self->origin[J_AXIS] ); 
			glVertex2f( self->origin[I_AXIS], self->origin[J_AXIS] + rodLength );
		glEnd();	
		glBegin(GL_TRIANGLES);
			glVertex2f( self->origin[I_AXIS], self->origin[J_AXIS] + rodLength );
			glVertex2f( self->origin[I_AXIS] + arrowHeadLength/2.0, self->origin[J_AXIS] + rodLength - arrowHeadLength);
			glVertex2f( self->origin[I_AXIS] - arrowHeadLength/2.0, self->origin[J_AXIS] + rodLength - arrowHeadLength);
		glEnd();
		lucPrint(self->origin[I_AXIS], self->origin[J_AXIS] + rodLength + arrowHeadLength, "Y");
	}
	else if ( dim == 3 ) {
		/* Drawing the X axis, by default using the RED color */
		lucColour_SetOpenGLColour( &self->colourX );

		glBegin(GL_TRIANGLES);
			glVertex3f( self->origin[I_AXIS] + rodLength, self->origin[J_AXIS], self->origin[K_AXIS] );
			glVertex3f( self->origin[I_AXIS] + rodLength - arrowHeadLength, 
				    self->origin[J_AXIS] - arrowHeadLength/2.0, self->origin[K_AXIS] );
			glVertex3f( self->origin[I_AXIS] + rodLength - arrowHeadLength,
				    self->origin[J_AXIS] + arrowHeadLength/2.0,
				    self->origin[K_AXIS] );
		glEnd();

		glBegin( GL_LINES );
			glVertex3f( self->origin[I_AXIS], self->origin[J_AXIS] , self->origin[K_AXIS] ); 
			glVertex3f( self->origin[I_AXIS] + rodLength, self->origin[J_AXIS] , self->origin[K_AXIS] );
		glEnd(); 
		
		lucPrint3d( self->origin[I_AXIS] + rodLength + textSpacing, self->origin[J_AXIS], self->origin[K_AXIS], "X");
		
		/* Drawing the Y axis, by default using the GREEN color */
		lucColour_SetOpenGLColour( &self->colourY );

		glBegin(GL_TRIANGLES);
			glVertex3f( self->origin[I_AXIS], self->origin[J_AXIS] + rodLength, self->origin[K_AXIS]  );
			glVertex3f( self->origin[I_AXIS] + arrowHeadLength/2.0, self->origin[J_AXIS] + rodLength -arrowHeadLength, 
				    self->origin[K_AXIS] );
			glVertex3f( self->origin[I_AXIS] - arrowHeadLength/2.0, self->origin[J_AXIS] + rodLength -arrowHeadLength, 
				    self->origin[K_AXIS] );
		glEnd();

		glBegin( GL_LINES );
			glVertex3f(  self->origin[I_AXIS], self->origin[J_AXIS] , self->origin[K_AXIS]  ); 
			glVertex3f(  self->origin[I_AXIS], self->origin[J_AXIS] + rodLength , self->origin[K_AXIS]  );
		glEnd();

		lucPrint3d( self->origin[I_AXIS], self->origin[J_AXIS]+ rodLength + arrowHeadLength, self->origin[K_AXIS], "Y");
		
		
		/* Drawing the Z axis, by default using the BLUE color */
		lucColour_SetOpenGLColour( &self->colourZ );
		glBegin(GL_TRIANGLES);
			glVertex3f( self->origin[I_AXIS], self->origin[J_AXIS] , self->origin[K_AXIS] + rodLength );
			glVertex3f( self->origin[I_AXIS] + arrowHeadLength/2.0, self->origin[J_AXIS] , 
				    self->origin[K_AXIS] + rodLength - arrowHeadLength );
			glVertex3f( self->origin[I_AXIS] - arrowHeadLength/2.0, self->origin[J_AXIS], 
				    self->origin[K_AXIS] + rodLength -arrowHeadLength );
		glEnd();

		glBegin( GL_LINES );
			glVertex3f(  self->origin[I_AXIS], self->origin[J_AXIS] , self->origin[K_AXIS] ); 
			glVertex3f( self->origin[I_AXIS], self->origin[J_AXIS] , self->origin[K_AXIS] + rodLength );
		glEnd(); 

		lucPrint3d( self->origin[I_AXIS], self->origin[J_AXIS] , self->origin[K_AXIS] + rodLength + textSpacing, "Z");
	}
	/* Put back settings */
	glEnable(GL_LIGHTING);
}

void _lucAxis_CleanUp( void* drawingObject, void* _context ) {
	lucAxis*      self            = (lucAxis*)drawingObject;
	
	_lucOpenGLDrawingObject_CleanUp( self, _context );

}

void _lucAxis_BuildDisplayList( void* drawingObject, void* _context ) {
}







