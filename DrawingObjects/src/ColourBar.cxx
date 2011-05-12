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
** $Id: ColourBar.c 791 2008-09-01 02:09:06Z JulianGiordani $
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
lucColourBar* _lucColourBar_New(  LUCCOLOURBAR_DEFARGS  ) 
{
	lucColourBar*					self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( _sizeOfSelf >= sizeof(lucColourBar) );
	self = (lucColourBar*) _lucDrawingObject_New(  LUCDRAWINGOBJECT_PASSARGS  );
	
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
		float                                                        scaleValue,
        double                                                       tickValues[] )
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
    int i;
    for (i = 0; i < 10; i++)    
        self->tickValues[i] = tickValues[i];
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
	/* Variables set in this function */
	SizeT                                              _sizeOfSelf = sizeof(lucColourBar);
	Type                                                      type = lucColourBar_Type;
	Stg_Class_DeleteFunction*                              _delete = _lucColourBar_Delete;
	Stg_Class_PrintFunction*                                _print = _lucColourBar_Print;
	Stg_Class_CopyFunction*                                  _copy = NULL;
	Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _lucColourBar_DefaultNew;
	Stg_Component_ConstructFunction*                    _construct = _lucColourBar_AssignFromXML;
	Stg_Component_BuildFunction*                            _build = _lucColourBar_Build;
	Stg_Component_InitialiseFunction*                  _initialise = _lucColourBar_Initialise;
	Stg_Component_ExecuteFunction*                        _execute = _lucColourBar_Execute;
	Stg_Component_DestroyFunction*                        _destroy = _lucColourBar_Destroy;
	lucDrawingObject_SetupFunction*                         _setup = _lucColourBar_Setup;
	lucDrawingObject_DrawFunction*                           _draw = _lucColourBar_Draw;
	lucDrawingObject_CleanUpFunction*                     _cleanUp = _lucColourBar_CleanUp;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

	return (void*) _lucColourBar_New(  LUCCOLOURBAR_PASSARGS  );
}

void _lucColourBar_AssignFromXML( void* drawingObject, Stg_ComponentFactory* cf, void* data ){
	lucColourBar*   self = (lucColourBar*)drawingObject;
	lucColourMap*   colourMap;
    unsigned int    i, defaultTicks, ticks;
    double          tickValues[11];
    char            tickLabel[10];

	/* Construct Parent */
	_lucDrawingObject_AssignFromXML( self, cf, data );

	colourMap     =  Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"ColourMap", lucColourMap, True, data  ) ;

    /* Default to 0 tick marks for linear, 1 for fixed centre, 2 for logarithmic scale */
    defaultTicks = 0;
    if (colourMap->centreOnFixedValue) defaultTicks = 1;
    if (colourMap->logScale) defaultTicks = 2;
	ticks = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, (Dictionary_Entry_Key)"ticks", defaultTicks );
    if (ticks > 9) ticks = 9;

    /* Load any provided intermediate tick values (tick1-9) */
    for (i = 1; i < ticks+1; i++ ) {
        sprintf(tickLabel, "tick%d", i);
		tickValues[i] = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)tickLabel, colourMap->maximum + 1 );
	}
    tickValues[0] = colourMap->minimum;
    tickValues[ticks+1] = colourMap->maximum;

	_lucColourBar_Init( 
			self, 
			colourMap,
			Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"lengthFactor", 0.8  ),
			Stg_ComponentFactory_GetUnsignedInt( cf, self->name, (Dictionary_Entry_Key)"height", 10  ),
			Stg_ComponentFactory_GetUnsignedInt( cf, self->name, (Dictionary_Entry_Key)"margin", 16  ),
			(float) Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"borderWidth", 1  ) ,
			Stg_ComponentFactory_GetUnsignedInt( cf, self->name, (Dictionary_Entry_Key)"precision", 2  ) ,
			Stg_ComponentFactory_GetBool( cf, self->name, (Dictionary_Entry_Key)"scientific", False  ),
   			ticks,
			Stg_ComponentFactory_GetBool( cf, self->name, (Dictionary_Entry_Key)"printTickValue", True  ),
			(float) Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"scaleValue", 1.0  ),
            tickValues
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
	int                      	startx, starty;
	char                     string[20];
	int i = 0;
	/* AbstractContext*         context = (AbstractContext*) _context; */

	/* Only get master to draw colour bar - /
	if ( context->rank != MASTER ) return; */

	/* Set up 2D Viewer the size of the viewport */
	lucViewport2d(True, viewportInfo);

	lucSetFontCharset(FONT_SMALL);
	
	startx = (viewportInfo->width - length)/2;
	starty = viewportInfo->height - self->margin - self->height;

    /* Write scale */
    lucColour_SetComplimentaryOpenGLColour( &window->backgroundColour );
    glLineWidth( self->borderWidth );

    /* Update min/max end ticks */
    self->tickValues[0] = colourMap->minimum;
    self->tickValues[self->ticks+1] = colourMap->maximum;

    for (i = 0; i < self->ticks+2; i++){
        /* Calculate tick position */
        float scaledPos;
        if (self->tickValues[i] > colourMap->maximum)
        {
            /* First get scaled position 0-1 */
            if (colourMap->logScale)
            {
                /* Space ticks based on a logarithmic scale of log(1) to log(11)
                   shows non-linearity while keeping ticks spaced apart enough to read labels */
                float tickpos = 1.0 + (float)i * (10.0 / (self->ticks+1));
                scaledPos = (log10(tickpos) / log10(11));
            }
            else
                /* Default linear scale evenly spaced ticks */
                scaledPos = (float)i / (self->ticks+1);

            /* Compute the tick value */
            if (colourMap->logScale ) {
                /* Reverse calc to find Value tick value at calculated position 0-1: */
                tickValue = log10(colourMap->minimum) + scaledPos 
                          * (log10(colourMap->maximum) - log10(colourMap->minimum));
                tickValue = pow( 10.0, tickValue );
            }
            else
            {
                /* Reverse scale calc and find value of tick at position 0-1 */
                if (colourMap->centreOnFixedValue && colourMap->centringValue > 0)
                {
                    /* Using fixed centre value, even linear scales either side */
                    if (scaledPos > 0.5)
                        tickValue = (colourMap->maximum - colourMap->centringValue)
                                  * (scaledPos - 0.5) / 0.5 + colourMap->centringValue;
                    else
                        tickValue = (colourMap->centringValue - colourMap->minimum)
                                  * scaledPos / 0.5 + colourMap->minimum;
                }
                else
                    /* Using even linear scale */
                    tickValue = colourMap->minimum + scaledPos * (colourMap->maximum - colourMap->minimum);
            }
        }
        else
        {
            /* User specified value */
            tickValue = self->tickValues[i];
            /* Calculate scaled position from value */
            scaledPos = lucColourMap_ScaleValue(colourMap, tickValue);
        }

        /* Calculate pixel position */
        int xpos = startx + length * scaledPos;

        /* Draws the tick */
        glBegin(GL_LINES);
            glVertex2i( xpos, starty+5+self->height );
            glVertex2i( xpos, starty+self->height);
        glEnd();
        
 
        /* Always print end values, print others if flag set */	
        if (self->printTickValue || i == 0 || i == self->ticks+1) 
        {
            if ( !self->scientific && fabs(tickValue) < 1.0e-5 )
                sprintf( string, "0" );
            else
                _lucColourBar_WithPrecision( string, self->scientific, self->precision, self->scaleValue, tickValue);

            lucPrint(xpos - (int) (0.5 * (float)lucStringWidth(string)),  starty + 13, string );
        }
    }

    /* Draw Colour Bar */
    for ( pixel_I = 0 ; pixel_I < length ; pixel_I++ ) {
        value = ((float)pixel_I / length);
        lucColourMap_SetOpenGLColourFromScaledValue( colourMap, value);
        glRecti( startx + pixel_I, starty, startx + pixel_I + 1 , starty + height );
    }

	/* Draw Box around colour bar */
	lucColour_SetComplimentaryOpenGLColour( &window->backgroundColour );
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	glRecti(startx, starty, startx + length, starty + height);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

	lucSetFontCharset(FONT_DEFAULT);

	/* Restore the viewport */
	lucViewport2d(False, viewportInfo);
}

void _lucColourBar_CleanUp( void* drawingObject, void* _context ) {
}


