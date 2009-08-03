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
** $Id: ColourMap.c 784 2008-08-18 13:54:31Z LukeHodkinson $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include "types.h"
#include "ColourMap.h"
#include "X11Colours.h"

#include <string.h>
#include <assert.h>

const Type lucColourMap_Type = "lucColourMap";

lucColourMap* lucColourMap_New( 
		Name                                               name,
		char*                                              colourMapString, 
		double                                             minimum,
		double                                             maximum,
		Bool                                               logScale,
		Bool                                               dynamicRange,
		Bool											   centreOnFixedValue,
		double											   centringValue	 )
{
	lucColourMap* self = (lucColourMap*) _lucColourMap_DefaultNew( name );

	lucColourMap_InitAll( self, colourMapString, minimum, maximum, logScale, dynamicRange, centreOnFixedValue, centringValue );

	return self;
}

lucColourMap* _lucColourMap_New(
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
		Name                                               name )
{
	lucColourMap*    self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( sizeOfSelf >= sizeof(lucColourMap) );
	self = (lucColourMap*) _Stg_Component_New( 
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
			name, 
			NON_GLOBAL );
	
	return self;
}

void _lucColourMap_Init( 
		lucColourMap*                                      self, 
		char*                                              _colourMapString, 
		double                                             minimum,
		double                                             maximum,
		Bool                                               logScale,
		Bool                                               dynamicRange,
		Bool											   centreOnFixedValue,
		double											   centringValue	 )
{	
	char*        colourMapString = StG_Strdup( _colourMapString );
	Colour_Index colourCount;
	Colour_Index colour_I;
	char*        charPointer;
	lucColour*    colour;
	char*        breakChars = " \t\n;,"; 

	self->minimum      = minimum;
	self->maximum      = maximum;
	self->logScale     = logScale;
	self->dynamicRange = dynamicRange;
	self->centreOnFixedValue = centreOnFixedValue;

  if( centringValue == maximum || centringValue == minimum )
    self->centringValue = 0.5*(maximum-minimum) + minimum;
  else
    self->centringValue = centringValue;
	
	/* Find number of colours */
	colourCount = 1;
	charPointer = strpbrk( colourMapString, breakChars );
	while ( charPointer != NULL ) {
		charPointer = strpbrk( charPointer + 1 , breakChars );
		colourCount++;
	}
	self->colourCount = colourCount;

	/* Allocate space for colour map */
	self->colourList = Memory_Alloc_Array( lucColour , colourCount, "Colour map");

	/* Initialise the space in the colour map */
	for( colour_I = 0 ; colour_I < colourCount ; colour_I++ ) {
		self->colourList[colour_I].red = 0;
		self->colourList[colour_I].green = 0;
		self->colourList[colour_I].blue = 0;
		self->colourList[colour_I].opacity = 0;
	}

	/* Read String to get colour map */
	charPointer = strtok( colourMapString, breakChars );
	for ( colour_I = 0 ; colour_I < colourCount ; colour_I++ ) {
		colour = lucColourMap_GetColourFromList( self, colour_I );
		lucColour_FromString( colour, charPointer );
		charPointer = strtok( NULL, breakChars );
	}

	Memory_Free( colourMapString );
}

void lucColourMap_InitAll( 
		void*                                              colourMap, 
		char*                                              colourMapString, 
		double                                             minimum,
		double                                             maximum,
		Bool                                               logScale,
		Bool                                               dynamicRange, 
		Bool											   centreOnFixedValue,
		double											   centringValue	 )
{
	lucColourMap* self        = colourMap;

	_lucColourMap_Init( self, colourMapString, minimum, maximum, logScale, dynamicRange, centreOnFixedValue, centringValue );
}
		
void _lucColourMap_Delete( void* colourMap ) {
	lucColourMap* self        = colourMap;

	Memory_Free( self->colourList );

	_Stg_Component_Delete( self );
}

void _lucColourMap_Print( void* colourMap, Stream* stream ) {
	lucColourMap* self        = colourMap;
	Colour_Index  colour_I;
	lucColour*    colour;
	lucColour     black;
	lucColour     white;
	lucColour     colourFromValue;
	Index         charCount = 100;
	Index         char_I;
	
	black.red = black.green = black.blue = 0.0;
	white.red = white.green = white.blue = 1.0;

	Journal_Printf( stream, "lucColourMap: %s\n", self->name );

	/* Print Parent */
	_Stg_Component_Print( self, stream );

	Journal_PrintValue( stream, self->colourCount );
	for ( colour_I = 0 ; colour_I < self->colourCount ; colour_I++ ) {
		colour = lucColourMap_GetColourFromList( self, colour_I );
		Journal_Printf( stream, "\tColour %u: Red - %6.3g, Green - %6.3g, Blue - %6.3g, Opacity - %6.3g\n", 
				colour_I, colour->red, colour->green, colour->blue, colour->opacity );
	}
	
	Journal_PrintValue( stream, self->minimum );
	Journal_PrintValue( stream, self->maximum );
	Journal_PrintBool( stream, self->logScale );
	Journal_PrintBool( stream, self->dynamicRange );

	/* Print Colour Map using terminal colours */
	for ( char_I = 0 ; char_I < charCount ; char_I++ ) {
		double value = (double) char_I / (double) (charCount - 1);
		lucColourMap_GetColourFromValue( self, value, &colourFromValue );
		lucColour_SetTerminalColours( &black, &colourFromValue, stream );
		Journal_Printf( stream, " " );
	}
	lucColour_SetTerminalColours( &white, &black, stream );
	Journal_Printf( stream, "\n" );

	/* Print Scale For Colour Map */
	Journal_Printf( stream, "%.4g", self->minimum );
	for ( char_I = 0 ; char_I < charCount ; char_I++ ) 
		Journal_Printf( stream, " " );
	Journal_Printf( stream, "%.4g\n", self->maximum );

}

void* _lucColourMap_Copy( void* colourMap, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	lucColourMap* self        = colourMap;
	lucColourMap* newColourMap;

	newColourMap = _Stg_Component_Copy( self, dest, deep, nameExt, ptrMap );

	newColourMap->colourCount         = self->colourCount;
	newColourMap->minimum             = self->minimum;
	newColourMap->maximum             = self->maximum;
	newColourMap->logScale            = self->logScale;
	newColourMap->dynamicRange        = self->dynamicRange;
	newColourMap->centreOnFixedValue  = self->centreOnFixedValue;
	newColourMap->centringValue       = self->centringValue;

	if (deep)
		memcpy( newColourMap->colourList, self->colourList, self->colourCount * sizeof(lucColour) );
	else 
		newColourMap->colourList = self->colourList;

	return (void*) newColourMap;
}


void* _lucColourMap_DefaultNew( Name name ) {
	return _lucColourMap_New( 
			sizeof( lucColourMap ),
			lucColourMap_Type,
			_lucColourMap_Delete,
			_lucColourMap_Print,
			_lucColourMap_Copy,
			_lucColourMap_DefaultNew,
			_lucColourMap_Construct,
			_lucColourMap_Build,
			_lucColourMap_Initialise,
			_lucColourMap_Execute,
			_lucColourMap_Destroy,
			name );
}

void _lucColourMap_Construct( void* colourMap, Stg_ComponentFactory* cf, void* data ) {
	lucColourMap* self             = (lucColourMap*) colourMap;

	_lucColourMap_Init( 
			self, 
			Stg_ComponentFactory_GetString( cf, self->name, "colours", "Blue;White;Red"), 
			Stg_ComponentFactory_GetDouble( cf, self->name, "minimum", 0.0 ),
			Stg_ComponentFactory_GetDouble( cf, self->name, "maximum", 1.0 ),
			Stg_ComponentFactory_GetBool( cf, self->name, "logScale", False ),
			Stg_ComponentFactory_GetBool( cf, self->name, "dynamicRange", False ),
			Stg_ComponentFactory_GetBool( cf, self->name, "centreOnFixedValue", False ),
			Stg_ComponentFactory_GetDouble( cf, self->name, "centringValue", 0.0 )
	);

    self->discrete = Stg_ComponentFactory_GetBool( cf, self->name, "discrete", False );
}

void _lucColourMap_Build( void* colourMap, void* data ) { }
void _lucColourMap_Initialise( void* colourMap, void* data ) { }
void _lucColourMap_Execute( void* colourMap, void* data ) { }
void _lucColourMap_Destroy( void* colourMap, void* data ) { }

void lucColourMap_GetColourFromValue( void* colourMap, double value, lucColour* colour ) {
	lucColourMap* self        = colourMap;

    /* Scale value to range [0,1] */
    float scaledValue = lucColourMap_ScaleValue(colourMap, value);

    /* Convert scaled value to colour */
    lucColourMap_GetColourFromScaledValue(colourMap, scaledValue, colour );

    /* Check for invalid range - set colour to invisible */
    if (self->maximum == self->minimum) colour->opacity = 0;
}

float lucColourMap_ScaleValue( void* colourMap, double value ) {
	lucColourMap* self        = colourMap;
	float         scaledValue;
	float 	      max, min, centre, sampleValue;

    /* To get a log scale, transform each value to log10(value) */
	if (self->logScale == True) {
		max 	    = log10(self->maximum);
		min  	    = log10(self->minimum);
		sampleValue = log10(value);
	}
	else {
		max 	    = self->maximum;
		min  	    = self->minimum;
		centre 	    = self->centringValue;
		sampleValue = value;
	}
	
    /* Scale value so that it is between 0 and 1 */
	if (self->logScale || !self->centreOnFixedValue) {
        scaledValue = (sampleValue - min) / (max - min);
    }
    else
    {
    	/* Scale value so that it is between 0 and 1, taking into account
    		the fact that centringValue should be at a scaled value 0.5 */
    	if (sampleValue > centre)
    		scaledValue = 0.5 + 0.5 * (sampleValue - centre)/(max - centre);
    	else
    		scaledValue = 0.5 * (sampleValue - min) / (centre - min); 
    }

    return scaledValue;
}

void lucColourMap_GetColourFromScaledValue( void* colourMap, float scaledValue, lucColour* colour ) 
{
	lucColourMap* self        = colourMap;
	Colour_Index  colourBelow_I;
	lucColour*    colourBelow;
	lucColour*    colourAbove;
	float         remainder;
	Colour_Index  colourCount = self->colourCount;

    /* Check within range */
    if (scaledValue <= 0.0 || colourCount == 1) {
		memcpy( colour, lucColourMap_GetColourFromList( self, 0 ), sizeof(lucColour) );
		return;
	}
	if (scaledValue >= 1.0) {
		memcpy( colour, lucColourMap_GetColourFromList( self, colourCount - 1 ), sizeof(lucColour) );
		return;
	}

	/* Discrete colourmap option does not interpolate between colours */
    if( !self->discrete ) {
       colourBelow_I = (Colour_Index) ( ( colourCount - 1 ) * scaledValue );
       colourBelow   = lucColourMap_GetColourFromList( self, colourBelow_I );
       colourAbove   = lucColourMap_GetColourFromList( self, colourBelow_I + 1 );

       remainder = (float)( colourCount - 1 ) * scaledValue - (float) colourBelow_I;

       /* Do linear interpolation between colours */
       colour->red     = ( colourAbove->red     - colourBelow->red     ) * remainder + colourBelow->red;
       colour->green   = ( colourAbove->green   - colourBelow->green   ) * remainder + colourBelow->green;
       colour->blue    = ( colourAbove->blue    - colourBelow->blue    ) * remainder + colourBelow->blue;
       colour->opacity = ( colourAbove->opacity - colourBelow->opacity ) * remainder + colourBelow->opacity;
    }
    else {
       colourBelow_I = (Colour_Index) ( (float)colourCount * scaledValue );
       colourBelow = lucColourMap_GetColourFromList( self, colourBelow_I );
       colour->red = colourBelow->red;
       colour->green = colourBelow->green;
       colour->blue = colourBelow->blue;
       colour->opacity = colourBelow->opacity;
    }
}

void lucColourMap_SetMinMax( void* colourMap, double min, double max ) {
	lucColourMap* self        = colourMap;
	double        tolerance   = 1e-10;

	/* Shift max and min if they are too close /
	if (fabs(min - max) < tolerance) {	
		max += 0.5 * tolerance;
		min -= 0.5 * tolerance;
	}
    Removed, caused attempt to draw colourBar with incorrect values
    now checked in colour calculations */

	/* Copy to colour map */
	self->minimum = min;
	self->maximum = max;
	
    /* If a centringValue has been imposed on a dynamic problem it should
        force the max / min to contain it correctly */
            
    if (Num_Approx(self->centringValue, 0.0)) self->centringValue = 0.5 * (max - min) + min;

    if (self->centringValue < min)
        min = self->centringValue - tolerance;
    
    if (self->centringValue > max)	
        max = self->centringValue + tolerance;
		
}

void lucColourMap_CalibrateFromVariable( void* colourMap, void* _variable ) {
	lucColourMap* self        = colourMap;
	Variable*     variable    = (Variable*)_variable;
	Index         array_I;
	Index         arrayCount  = *variable->arraySizePtr;
	double        value;
	double        localMin    = 1.0e99;
	double        localMax    = -1.0e99;
	double        globalMax;
	double        globalMin;

	if ( !self->dynamicRange ) return;

	for ( array_I = 0 ; array_I < arrayCount ; array_I++ ){
		/* Get scalar value from particle */
		value = Variable_GetValueDouble( variable, array_I ) ;

		if ( value < localMin ) 
			localMin = value;
		else if ( value > localMax ) 
			localMax = value;
	}
	
	MPI_Allreduce( &localMin, &globalMin, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD );
	MPI_Allreduce( &localMax, &globalMax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD );

	lucColourMap_SetMinMax( self, globalMin, globalMax );
}

void lucColourMap_CalibrateFromFieldVariable( void* colourMap, void* _fieldVariable ) {
	lucColourMap*   self          = colourMap;
	FieldVariable*  fieldVariable = (FieldVariable*)_fieldVariable;

	if ( !self->dynamicRange ) return;

	lucColourMap_SetMinMax( 
			self, 
			FieldVariable_GetMinGlobalFieldMagnitude( fieldVariable ), 
			FieldVariable_GetMaxGlobalFieldMagnitude( fieldVariable ) );
}

void lucColourMap_CalibrateFromSwarmVariable( void* colourMap, void* swarmVariable ) {
	lucColourMap*   self          = colourMap;

	if ( !self->dynamicRange ) return;

	lucColourMap_SetMinMax( 
			self, 
			SwarmVariable_GetMinGlobalMagnitude( swarmVariable ), 
			SwarmVariable_GetMaxGlobalMagnitude( swarmVariable ) );
}

void lucColour_FromHSV( lucColour* self, float hue, float saturation, float value, float opacity ) {
	int   Hi     = (int) ( hue/60.0 );
	float f      = hue/60.0 - (float) Hi;
	float p      = value * ( 1.0 - saturation );
	float q      = value * ( 1.0 - saturation * f );
	float t      = value * ( 1.0 - saturation * ( 1.0 - f ));

	switch ( Hi ){
		case 0:
			self->red = value;   self->green = t;       self->blue = p;     break;
		case 1:
			self->red = q;       self->green = value;   self->blue = p;     break;
		case 2:
			self->red = p;       self->green = value;   self->blue = t;     break;
		case 3:
			self->red = p;       self->green = q;       self->blue = value; break;
		case 4:
			self->red = t;       self->green = p;       self->blue = value; break;
		case 5:
			self->red = value;   self->green = p;       self->blue = q;     break;
	}

	self->opacity = opacity;
}
	
	

void lucColour_FromString( lucColour* self, char* string ) {
	char* charPointer;
	float opacity;

	lucColour_FromX11ColourName( self, string );

	/* Get Opacity From String */
	/* Opacity must be read in after the ":" of the name of the colour */
	charPointer = strchr( string, ':' );

	if (charPointer != NULL) {
		/* Return full opactity (non-transparent) if no opacity is set */
		if (sscanf( charPointer + 1, "%f", &opacity )	!= 1) 
			opacity = 1.0;
	}
	else 
		opacity = 1.0;
	self->opacity = opacity;
}

typedef enum {
	Terminal_Black,
	Terminal_Red,
	Terminal_Green,
	Terminal_Yellow,
	Terminal_Blue,
	Terminal_Magenta,
	Terminal_Cyan,
	Terminal_Grey,
	Terminal_White 
} Terminal_Colour;

Terminal_Colour lucColour_GetClosestTerminalColour( lucColour* self ) {
	Bool hasRed   = ( self->red   > 0.5 );
	Bool hasGreen = ( self->green > 0.5 );
	Bool hasBlue  = ( self->blue  > 0.5 );
	
	/* Shades */
	if ( ! hasRed && ! hasGreen && ! hasBlue ) 
		return Terminal_Black;

	if ( hasRed && hasGreen && hasBlue ) 
		return Terminal_White;
	
	/* Primary Colours */
	if ( hasRed && ! hasGreen && ! hasBlue ) 
		return Terminal_Red;

	if ( ! hasRed && hasGreen && ! hasBlue ) 
		return Terminal_Green;
	
	if ( ! hasRed && ! hasGreen && hasBlue ) 
		return Terminal_Blue;

	/* Secondary Colours */
	if ( hasRed && ! hasGreen && hasBlue ) 
		return Terminal_Magenta;

	if ( ! hasRed && hasGreen && hasBlue ) 
		return Terminal_Cyan;
	
	if ( hasRed && hasGreen && ! hasBlue ) 
		return Terminal_Yellow;

	abort();
	return Terminal_Black;
}

void lucColour_SetTerminalColours( lucColour* textColour, lucColour* backgroundColour, Stream* stream ) {
	int reset = 0;
	
	/* Command is the control command to the terminal */
	Journal_Printf( stream, "%c[%d;%d;%dm", 
			0x1B, 
			reset, 
			lucColour_GetClosestTerminalColour( textColour ) + 30, 
			lucColour_GetClosestTerminalColour( backgroundColour ) + 40 );
}
