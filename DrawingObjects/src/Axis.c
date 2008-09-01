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
#ifdef HAVE_OPENGL_FRAMEWORK
	#include <OpenGL/gl.h>
	#include <OpenGL/glu.h>
#else
	#include <gl.h>
	#include <glu.h>
#endif
#include <string.h>

#ifndef MASTER
	#define MASTER 0
#endif

const Type lucAxis_Type = "lucAxis";

lucAxis* lucAxis_New( 
		Name                                                  name,
		Coord                                                 origin,
    		float 				                      scale,
		lucColour                                             colourX,
		lucColour                                             colourY,
		lucColour                                             colourZ)
{
	lucAxis* self = (lucAxis*) _lucAxis_DefaultNew( name );

	lucAxis_InitAll( self, origin, scale, colourX, colourY, colourZ);

	return self;
}

lucAxis* _lucAxis_New(
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
	lucAxis*    self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( sizeOfSelf >= sizeof(lucAxis) );
	self = (lucAxis*)  _lucOpenGLDrawingObject_New( 
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

void lucAxis_Init(		
		lucAxis*                                         self,
		Coord                                            origin,
		float 				                 scale,
		lucColour                                        colourX,
		lucColour                                        colourY,
		lucColour                                        colourZ) 
{
	
	self->scale = scale;
	memcpy( self->origin, origin, sizeof(Coord) );	
	memcpy( &(self->colourX), &colourX, sizeof(lucColour) );	
	memcpy( &(self->colourY), &colourY, sizeof(lucColour) );	
	memcpy( &(self->colourZ), &colourZ, sizeof(lucColour) );	

	
}

void lucAxis_InitAll( 
		void*                                              axis,
		Coord                                              origin,
	        float 				                   scale,
		lucColour                                          colourX,
		lucColour                                          colourY,
		lucColour                                          colourZ)
{
	lucAxis* self        = axis;

	/* TODO Init parent */
	lucAxis_Init( self, origin, scale, colourX, colourY, colourZ );
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

	newAxis->scale = self->scale;
	memcpy( &(newAxis->colourX),       &(self->colourX),       sizeof(lucColour) );
	memcpy( &(newAxis->colourY),       &(self->colourY),       sizeof(lucColour) );
	memcpy( &(newAxis->colourZ),       &(self->colourZ),       sizeof(lucColour) );
	memcpy( newAxis->origin,       self->origin,       sizeof(Coord) );
	
	return (void*) newAxis;
}

void* _lucAxis_DefaultNew( Name name ) {
	return _lucAxis_New( 
			sizeof( lucAxis ),
			lucAxis_Type,
			_lucAxis_Delete,
			_lucAxis_Print,
			_lucAxis_Copy,
			_lucAxis_DefaultNew,
			_lucAxis_Construct,
			_lucAxis_Build,
			_lucAxis_Initialise,
			_lucAxis_Execute,
			_lucAxis_Destroy,		
		        _lucAxis_Setup,
			_lucAxis_Draw,
	                _lucAxis_CleanUp,
	 		_lucAxis_BuildDisplayList,
			name );
}

void _lucAxis_Construct( void* axis, Stg_ComponentFactory* cf, void* data ) {
	lucAxis*             self               = (lucAxis*) axis;
        Name colourNameX;
	Name colourNameY;	
	Name colourNameZ;		
	
	Coord origin;
	
	/* Get Stereo Type */
         /* Construct Parent */
	_lucDrawingObject_Construct( self, cf, data );

	colourNameX  = Stg_ComponentFactory_GetString( cf, self->name, "colourX", "Red") ;
	colourNameY  = Stg_ComponentFactory_GetString( cf, self->name, "colourY", "Green") ;
	colourNameZ  = Stg_ComponentFactory_GetString( cf, self->name, "colourZ", "Blue") ;
	
	lucColour_FromString( &self->colourX, colourNameX );	
	lucColour_FromString( &self->colourY, colourNameY );
	lucColour_FromString( &self->colourZ, colourNameZ );
	
	origin[I_AXIS]  = Stg_ComponentFactory_GetDouble( cf, self->name, "originX", -0.25 );
	origin[J_AXIS]  = Stg_ComponentFactory_GetDouble( cf, self->name, "originY", 0.0 );
	origin[K_AXIS]  = Stg_ComponentFactory_GetDouble( cf, self->name, "originZ", 0.0 );
	
       	lucAxis_InitAll( self, 
	                origin,
			Stg_ComponentFactory_GetDouble( cf, self->name, "scale", 1.0),
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
	lucAxis*      self            = (lucAxis*)drawingObject;
        lucViewport* viewport = viewportInfo->viewport;
	DomainContext*   context = (DomainContext*) _context;
	Dimension_Index          dim     = context->dim;

        double rodLength = 0.0;
	double arrowLength = 0.0;

		
	/* Initialise OpenGL stuff */
	glShadeModel(GL_SMOOTH);
	glDisable(GL_LIGHTING);

	glEnable(GL_LINE_SMOOTH);
	glEnable(GL_BLEND);	
	
	/* Disable lighting because we don't want a 3D effect */
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

	/* The rodlength is the total length of the arrow line.
	   By default it is 0.25, but can be scaled */
         rodLength =  self->scale*0.25;
	 
	 /* The tip of the arrow hea starts at rofLength. The size of the arrow
	 is a fifth of the total length
	 ---->
	 */
	 arrowLength = rodLength/5.0;
	 
	 if (dim == 2) {
	        /* Drawing the X axis, default is the RED color */
		glColor3f(self->colourX.red, self->colourX.green, self->colourX.blue);

		glBegin( GL_LINES );
			glVertex2f( self->origin[I_AXIS], self->origin[J_AXIS] ); 
			glVertex2f( self->origin[I_AXIS] + self->scale*0.25 , self->origin[J_AXIS]  );
		glEnd(); 
		glBegin(GL_TRIANGLES);
			glVertex2f( self->origin[I_AXIS] + self->scale*0.25, self->origin[J_AXIS] );
			glVertex2f( self->origin[I_AXIS] + self->scale*0.25 - arrowLength, self->origin[J_AXIS] -arrowLength/2.0);
			glVertex2f( self->origin[I_AXIS] + self->scale*0.25 - arrowLength, self->origin[J_AXIS] +arrowLength/2.0);
		glEnd();
		glRasterPos2f( self->origin[I_AXIS] + self->scale*0.25 + self->scale*0.05, self->origin[J_AXIS] );
		lucPrintString( "X");
		
		/* Drawing the Y axis, default is the GREEN color */
		glColor3f(self->colourY.red, self->colourY.green, self->colourY.blue);

		glBegin( GL_LINES );
			glVertex2f( self->origin[I_AXIS], self->origin[J_AXIS] ); 
			glVertex2f( self->origin[I_AXIS], self->origin[J_AXIS] + self->scale*0.25 );
		glEnd();	
		glBegin(GL_TRIANGLES);
			glVertex2f( self->origin[I_AXIS], self->origin[J_AXIS] + self->scale*0.25 );
			glVertex2f( self->origin[I_AXIS] + arrowLength/2.0, self->origin[J_AXIS] + self->scale*0.25 -arrowLength);
			glVertex2f( self->origin[I_AXIS] - arrowLength/2.0, self->origin[J_AXIS] + self->scale*0.25 -arrowLength);
		glEnd();
		glRasterPos2f( self->origin[I_AXIS], self->origin[J_AXIS]+ self->scale*0.25 + self->scale*0.05 );	
		lucPrintString( "Y");
	}
	else if ( dim == 3 ) {
		/* Drawing the X axis, by default using the RED color */
		glColor3f(self->colourX.red, self->colourX.green, self->colourX.blue);

		glBegin(GL_TRIANGLES);
			glVertex3f( self->origin[I_AXIS] + self->scale*0.25, self->origin[J_AXIS], self->origin[K_AXIS] );
			glVertex3f( self->origin[I_AXIS] + self->scale*0.25 - arrowLength, 
				    self->origin[J_AXIS] -arrowLength/2.0, self->origin[K_AXIS] );
			glVertex3f( self->origin[I_AXIS] + self->scale*0.25 - arrowLength,
				    self->origin[J_AXIS] +arrowLength/2.0,
				    self->origin[K_AXIS] );
		glEnd();

		glBegin( GL_LINES );
			glVertex3f( self->origin[I_AXIS], self->origin[J_AXIS] , self->origin[K_AXIS] ); 
			glVertex3f( self->origin[I_AXIS] + self->scale*0.25, self->origin[J_AXIS] , self->origin[K_AXIS] );
		glEnd(); 
		
		glRasterPos3f( self->origin[I_AXIS] + self->scale*0.25 + self->scale*0.05, self->origin[J_AXIS], self->origin[K_AXIS] );
		lucPrintString( "X");
		
		/* Drawing the X axis, by default using the GREEN color */
		glColor3f(self->colourY.red, self->colourY.green, self->colourY.blue);

		glBegin(GL_TRIANGLES);
			glVertex3f( self->origin[I_AXIS], self->origin[J_AXIS] + self->scale*0.25, self->origin[K_AXIS]  );
			glVertex3f( self->origin[I_AXIS] + arrowLength/2.0, self->origin[J_AXIS] + self->scale*0.25 -arrowLength, 
				    self->origin[K_AXIS] );
			glVertex3f( self->origin[I_AXIS] - arrowLength/2.0, self->origin[J_AXIS] + self->scale*0.25 -arrowLength, 
				    self->origin[K_AXIS] );
		glEnd();

		glBegin( GL_LINES );
			glVertex3f(  self->origin[I_AXIS], self->origin[J_AXIS] , self->origin[K_AXIS]  ); 
			glVertex3f(  self->origin[I_AXIS], self->origin[J_AXIS] + self->scale*0.25 , self->origin[K_AXIS]  );
		glEnd();
		
		glRasterPos3f( self->origin[I_AXIS], self->origin[J_AXIS]+ self->scale*0.25 + self->scale*0.05, self->origin[K_AXIS] );
		lucPrintString( "Y");
		
		
		/* Drawing the X axis, by default using the BLUE color */
		glColor3f(self->colourZ.red, self->colourZ.green, self->colourZ.blue);
		glBegin(GL_TRIANGLES);
			glVertex3f( self->origin[I_AXIS], self->origin[J_AXIS] , self->origin[K_AXIS] + self->scale*0.25  );
			glVertex3f( self->origin[I_AXIS] + arrowLength/2.0, self->origin[J_AXIS] , 
				    self->origin[K_AXIS] + self->scale*0.25 - arrowLength );
			glVertex3f( self->origin[I_AXIS] - arrowLength/2.0, self->origin[J_AXIS], 
				    self->origin[K_AXIS] + self->scale*0.25 -arrowLength );
		glEnd();

		glBegin( GL_LINES );
			glVertex3f(  self->origin[I_AXIS], self->origin[J_AXIS] , self->origin[K_AXIS] ); 
			glVertex3f( self->origin[I_AXIS], self->origin[J_AXIS] , self->origin[K_AXIS] + self->scale*0.25 );
		glEnd(); 
		glRasterPos3f( self->origin[I_AXIS], self->origin[J_AXIS] , self->origin[K_AXIS] + self->scale*0.25+ self->scale*0.05 );
		lucPrintString("Z");
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





