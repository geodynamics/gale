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
** $Id: ColourBar.c 785 2008-08-18 13:55:00Z LukeHodkinson $
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
#include "ColourBar.h"

#include <assert.h>
#include <gl.h>
#include <glu.h>
#include <string.h>

#ifndef MASTER
	#define MASTER 0
#endif

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type lucColourBar_Type = "lucColourBar";

/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
lucColourBar* _lucColourBar_New( 
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
		Name                                               name ) 
{
	lucColourBar*					self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( sizeOfSelf >= sizeof(lucColourBar) );
	self = (lucColourBar*) _lucDrawingObject_New( 
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
			name );
	
	return self;
}

void _lucColourBar_Init( 
		lucColourBar*                                                self,
		lucColourMap*                                                colourMap,
		double                                                       lengthFactor,    
		Pixel_Index                                                  height,    
		Pixel_Index                                                  margin,
		float                                                        borderWidth,
		int                                                          precision,
		Bool                                                         scientific,
		int                                                          ticks,
		Bool                                                         printTickValue,
		float                                                        scaleValue )
{
	self->colourMap = colourMap;
	self->lengthFactor = lengthFactor;
	self->height = height;
	self->margin = margin;
	self->borderWidth = borderWidth;
	self->precision = precision;
	self->scientific = scientific;
	self->ticks = ticks;
	self->printTickValue = printTickValue;
	self->scaleValue = scaleValue;

}

void _lucColourBar_Delete( void* drawingObject ) {
	lucColourBar*  self = (lucColourBar*)drawingObject;

	_lucDrawingObject_Delete( self );
}

void _lucColourBar_Print( void* drawingObject, Stream* stream ) {
	lucColourBar*  self = (lucColourBar*)drawingObject;

	_lucDrawingObject_Print( self, stream );
}

void* _lucColourBar_Copy( void* drawingObject, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap) {
	lucColourBar*  self = (lucColourBar*)drawingObject;
	lucColourBar* newDrawingObject;

	newDrawingObject = _lucDrawingObject_Copy( self, dest, deep, nameExt, ptrMap );

	/* TODO */
	abort();

	return (void*) newDrawingObject;
}


void* _lucColourBar_DefaultNew( Name name ) {
	return (void*) _lucColourBar_New(
		sizeof(lucColourBar),
		lucColourBar_Type,
		_lucColourBar_Delete,
		_lucColourBar_Print,
		NULL,
		_lucColourBar_DefaultNew,
		_lucColourBar_Construct,
		_lucColourBar_Build,
		_lucColourBar_Initialise,
		_lucColourBar_Execute,
		_lucColourBar_Destroy,
		_lucColourBar_Setup,
		_lucColourBar_Draw,
		_lucColourBar_CleanUp,
		name );
}

void _lucColourBar_Construct( void* drawingObject, Stg_ComponentFactory* cf, void* data ){
	lucColourBar*  self = (lucColourBar*)drawingObject;
	lucColourMap*    colourMap;

	/* Construct Parent */
	_lucDrawingObject_Construct( self, cf, data );

	colourMap     =  Stg_ComponentFactory_ConstructByKey(  cf,  self->name,  "ColourMap", lucColourMap, True, data ) ;
	
	_lucColourBar_Init( 
			self, 
			colourMap,
			Stg_ComponentFactory_GetDouble( cf, self->name, "lengthFactor", 0.8 ),
			Stg_ComponentFactory_GetUnsignedInt( cf, self->name, "height",  10 ),
			Stg_ComponentFactory_GetUnsignedInt( cf, self->name, "margin",  16 ),
			(float) Stg_ComponentFactory_GetDouble( cf, self->name, "borderWidth", 1 ) ,
			Stg_ComponentFactory_GetUnsignedInt( cf, self->name, "precision", 2 ) ,
			Stg_ComponentFactory_GetBool(cf, self->name, "scientific", False ),

			Stg_ComponentFactory_GetUnsignedInt(cf, self->name, "ticks", 0 ) ,
			Stg_ComponentFactory_GetBool(cf, self->name, "printTickValue", False ),
			(float) Stg_ComponentFactory_GetDouble(cf, self->name, "scaleValue", 1.0 ) 
			);
}

void _lucColourBar_Build( void* drawingObject, void* data ) {}
void _lucColourBar_Initialise( void* drawingObject, void* data ) {}
void _lucColourBar_Execute( void* drawingObject, void* data ) {}
void _lucColourBar_Destroy( void* drawingObject, void* data ) {}

void _lucColourBar_Setup( void* drawingObject, void* _context ) {
}

void _lucColourBar_WithPrecision(char *string, Bool scientific, int precision, float scaleValue, double value){

if(precision > 5 ) precision = 2;

/* For display purpose, scales the printed values if needed */
value = scaleValue * value;

	if(scientific){
		if( precision == 1 ){
			sprintf(string, "%.1e", value);
			return;
			}
		if( precision == 2 ){
			sprintf(string, "%.2e", value);
			return;
			}
		if( precision == 3 ){
			sprintf(string, "%.3e", value);
			return;
			}
		if( precision == 4 ){
			sprintf(string, "%.4e", value);
			return;
			}
		if( precision == 5 ){
			sprintf(string, "%.5e", value);
			return;
			}
	}
	else{
		if( precision == 1 ){
			sprintf(string, "%.1g", value);
			return;
			}
		if( precision == 2 ){
			sprintf(string, "%.2g", value);
			return;
			}
		if( precision == 3 ){
			sprintf(string, "%.3g", value);
			return;
			}
		if( precision == 4 ){
			sprintf(string, "%.4g", value);
			return;
			}
		if( precision == 5 ){
			sprintf(string, "%.5g", value);
			return;
			}
	}

}
	
void _lucColourBar_Draw( void* drawingObject, lucWindow* window, lucViewportInfo* viewportInfo, void* _context ) {
	lucColourBar*            self          = (lucColourBar*)drawingObject;
	lucColourMap*            colourMap     = self->colourMap;
	Pixel_Index              length        = (Pixel_Index) ((double) viewportInfo->width * self->lengthFactor);
	Pixel_Index              height         = self->height;
	Pixel_Index              pixel_I;
	double                   value;
	double                   tickValue;
	int                      startPos[2];
	GLint                    rasterPos[2];
	char                     string[20];
	int                      stringWidth;
	int i = 0;
	AbstractContext*         context       = (AbstractContext*) _context;

	
	/* Only get master to draw colour bar */
	if ( context->rank != MASTER )
		return;

	/* Set up 2D Viewer the size of the viewport */
	
	
	//glPushMatrix();
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	gluOrtho2D((GLfloat) 0.0, (GLfloat) viewportInfo->width, (GLfloat) 0.0, (GLfloat) viewportInfo->height );
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	/* Disable lighting because we don't want a 3D effect */
	glDisable(GL_LIGHTING);

	startPos[0] = (viewportInfo->width - length)/2;
	startPos[1] = self->margin;
	
	/* Draw Colour Bar */
	for ( pixel_I = 0 ; pixel_I < length ; pixel_I++ ) {
           if( colourMap->logScale ) {
              value = log10(colourMap->minimum) +
                 (double)pixel_I * (log10(colourMap->maximum) - log10(colourMap->minimum)) / (double)length;
              value = pow( 10.0, value );
           }
           else
              value = colourMap->minimum + (double) pixel_I * (colourMap->maximum - colourMap->minimum) / (double) length;

		lucColourMap_SetOpenGLColourFromValue( colourMap, value );
	
		glRecti( startPos[0] + pixel_I, startPos[1], startPos[0] + pixel_I + 1 , startPos[1] + height );
	}

	/* Draw Box around colour bar */
	glLineWidth( self->borderWidth );
	lucColour_SetComplimentaryOpenGLColour( &window->backgroundColour );
	glBegin( GL_LINE_LOOP );
		glVertex2i( startPos[0], startPos[1] );
		glVertex2i( startPos[0] + length, startPos[1] );
		glVertex2i( startPos[0] + length, startPos[1] + height );
		glVertex2i( startPos[0] , startPos[1] + height );
	glEnd();

	/* Write scale */
	if ( !self->scientific && fabs(colourMap->minimum) < 1.0e-5 )
		sprintf( string, "0" );
	else{
		 _lucColourBar_WithPrecision(string, self->scientific, self->precision,  self->scaleValue, colourMap->minimum);
	}
	stringWidth = lucStringWidth( string );

	rasterPos[0] = startPos[0] - (int) (0.5 * (float)stringWidth);
	rasterPos[1] = startPos[1] - 13;
	if (rasterPos[0] < 0)
		rasterPos[0] = 1;

	glRasterPos2iv( rasterPos );
	lucPrintString( string );
	
	if ( !self->scientific && fabs(colourMap->maximum) < 1.0e-5 )
		sprintf( string, "0" );
	else{
		_lucColourBar_WithPrecision( string, self->scientific, self->precision, self->scaleValue, colourMap->maximum);
	}
	stringWidth = lucStringWidth( string );

	rasterPos[0] = startPos[0] + length - (int) (0.5 * (float)stringWidth);
	glRasterPos2iv( rasterPos );
	lucPrintString( string );

	/* Write ticks */
	for(i = 1; i< self->ticks; i++){
		/* Computes the tick position */
	       	rasterPos[0] = startPos[0]+ i*(length/ self->ticks);
        	rasterPos[1] = startPos[1] - 13;

                /* Draws the tick */
		glLineWidth( self->borderWidth );
		lucColour_SetComplimentaryOpenGLColour( &window->backgroundColour );
		glBegin(GL_LINES);
			glVertex2i( startPos[0]+ i*(length/ self->ticks), startPos[1]-5 );
			glVertex2i( startPos[0] + i*(length/ self->ticks), startPos[1] );
		glEnd();

		
		/* Computse the tick value */
                if( colourMap->logScale ) {
                   tickValue = log10(colourMap->minimum) + 
                      ((double)i * (log10(colourMap->maximum) - log10(colourMap->minimum)) / (double)self->ticks);
                   tickValue = pow( 10.0, tickValue );
                }
                else
                   tickValue = colourMap->minimum + ( i*(colourMap->maximum - colourMap->minimum)/self->ticks );
		
		if(self->printTickValue)  {
			_lucColourBar_WithPrecision( string, self->scientific, self->precision, self->scaleValue, tickValue);
			
			stringWidth = lucStringWidth( string );
			rasterPos[0] -= (int) (0.5 * (float)stringWidth);
			glRasterPos2iv( rasterPos );
			lucPrintString( string );
		}

	}

	/* Put back settings */
	glEnable(GL_LIGHTING);
	glPopMatrix();
	
	
	/*Set back the viewport to what it should be to render any other object */
	/* If this is not done, than any object displayed after the colour bar will not appear,*/
	/* because the projection matrix and lookAt point have been altered */
	lucViewportInfo_SetOpenGLCamera( viewportInfo );
}

void _lucColourBar_CleanUp( void* drawingObject, void* _context ) {
}
