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
** $Id: ScalarFieldCrossSection.c 791 2008-09-01 02:09:06Z JulianGiordani $
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
#include "ScalarFieldCrossSection.h"

#include <assert.h>
#ifdef HAVE_OPENGL_FRAMEWORK
	#include <OpenGL/gl.h>
	#include <OpenGL/glu.h>
#else
	#include <gl.h>
	#include <glu.h>
#endif
#include <string.h>
#include <ctype.h>

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type lucScalarFieldCrossSection_Type = "lucScalarFieldCrossSection";

/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
lucScalarFieldCrossSection* _lucScalarFieldCrossSection_New( 
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
	lucScalarFieldCrossSection*					self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( sizeOfSelf >= sizeof(lucScalarFieldCrossSection) );
	self = (lucScalarFieldCrossSection*) _lucOpenGLDrawingObject_New( 
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

void _lucScalarFieldCrossSection_Init( 
		lucScalarFieldCrossSection*                                  self,
		Name                                                         fieldVariableName,
		lucColourMap*                                                colourMap,
		IJK                                                          resolution,
		double                                                       crossSectionValue,
		Axis                                                         crossSectionAxis,
		XYZ                                                          minCropValues,
		XYZ                                                          maxCropValues ) 
{
//	self->fieldVariable = fieldVariable;
	self->fieldVariableName = fieldVariableName;
	self->colourMap = colourMap;
	memcpy( self->resolution, resolution, sizeof(IJK) );
	self->crossSectionValue = crossSectionValue;
	self->crossSectionAxis = crossSectionAxis;
	memcpy( self->minCropValues, minCropValues, sizeof(XYZ) );
	memcpy( self->maxCropValues, maxCropValues, sizeof(XYZ) );
}

void _lucScalarFieldCrossSection_Delete( void* drawingObject ) {
	lucScalarFieldCrossSection*  self = (lucScalarFieldCrossSection*)drawingObject;

	_lucOpenGLDrawingObject_Delete( self );
}

void _lucScalarFieldCrossSection_Print( void* drawingObject, Stream* stream ) {
	lucScalarFieldCrossSection*  self = (lucScalarFieldCrossSection*)drawingObject;

	_lucOpenGLDrawingObject_Print( self, stream );
}

void* _lucScalarFieldCrossSection_Copy( void* drawingObject, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap) {
	lucScalarFieldCrossSection*  self = (lucScalarFieldCrossSection*)drawingObject;
	lucScalarFieldCrossSection* newDrawingObject;

	newDrawingObject = _lucOpenGLDrawingObject_Copy( self, dest, deep, nameExt, ptrMap );

	/* TODO */
	abort();

	return (void*) newDrawingObject;
}


void* _lucScalarFieldCrossSection_DefaultNew( Name name ) {
	return (void*) _lucScalarFieldCrossSection_New(
		sizeof(lucScalarFieldCrossSection),
		lucScalarFieldCrossSection_Type,
		_lucScalarFieldCrossSection_Delete,
		_lucScalarFieldCrossSection_Print,
		NULL,
		_lucScalarFieldCrossSection_DefaultNew,
		_lucScalarFieldCrossSection_Construct,
		_lucScalarFieldCrossSection_Build,
		_lucScalarFieldCrossSection_Initialise,
		_lucScalarFieldCrossSection_Execute,
		_lucScalarFieldCrossSection_Destroy,
		_lucScalarFieldCrossSection_Setup,
		_lucScalarFieldCrossSection_Draw,
		_lucScalarFieldCrossSection_CleanUp,
		_lucScalarFieldCrossSection_BuildDisplayList,
		name );
}

void _lucScalarFieldCrossSection_Construct( void* drawingObject, Stg_ComponentFactory* cf, void* data ){
	lucScalarFieldCrossSection*     self = (lucScalarFieldCrossSection*)drawingObject;
	lucColourMap*    colourMap;
	Index            defaultResolution;
	IJK              resolution;
	char             axisChar;
	double           value               = 0.0;
	Axis             axis                = 0;
	Name             crossSectionName;
	Name             fieldVariableName;
	XYZ              minCropValues;
	XYZ              maxCropValues;

	/* Construct Parent */
	_lucOpenGLDrawingObject_Construct( self, cf, data );

	fieldVariableName = Stg_ComponentFactory_GetString( cf, self->name, "FieldVariable", "defaultName" );

	/* This variable is now construct in build phase.
	   fieldVariable =  Stg_ComponentFactory_ConstructByKey( cf, self->name, "FieldVariable", FieldVariable, True ) ;
	*/


	colourMap = Stg_ComponentFactory_ConstructByKey( cf, self->name, "ColourMap", lucColourMap, True, data ) ;

	defaultResolution = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, "resolution", 128 );
	resolution[ I_AXIS ] = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, "resolutionX", defaultResolution );
	resolution[ J_AXIS ] = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, "resolutionY", defaultResolution );
	resolution[ K_AXIS ] = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, "resolutionZ", defaultResolution );
			
	crossSectionName = Stg_ComponentFactory_GetString( cf, self->name, "crossSection", "" );
	if ( sscanf( crossSectionName, "%c=%lf", &axisChar, &value ) == 2 ) {
		if ( toupper( axisChar ) >= 'X' )
			axis = toupper( axisChar ) - 'X';
	}

	/* Get Values with which to crop the cross section */
	minCropValues[ I_AXIS ] = Stg_ComponentFactory_GetDouble( cf, self->name, "minCropX", -HUGE_VAL );
	minCropValues[ J_AXIS ] = Stg_ComponentFactory_GetDouble( cf, self->name, "minCropY", -HUGE_VAL );
	minCropValues[ K_AXIS ] = Stg_ComponentFactory_GetDouble( cf, self->name, "minCropZ", -HUGE_VAL );
	maxCropValues[ I_AXIS ] = Stg_ComponentFactory_GetDouble( cf, self->name, "maxCropX", +HUGE_VAL );
	maxCropValues[ J_AXIS ] = Stg_ComponentFactory_GetDouble( cf, self->name, "maxCropY", +HUGE_VAL );
	maxCropValues[ K_AXIS ] = Stg_ComponentFactory_GetDouble( cf, self->name, "maxCropZ", +HUGE_VAL );
	
	_lucScalarFieldCrossSection_Init( 
			self, 
			fieldVariableName,
			colourMap,
			resolution,
			value,
			axis,
			minCropValues,
			maxCropValues );
}

void _lucScalarFieldCrossSection_Build( void* drawingObject, void* data ) {
	lucScalarFieldCrossSection*     self        = (lucScalarFieldCrossSection*)drawingObject;
	AbstractContext*                context     = Stg_CheckType( data, AbstractContext );
	Stg_ComponentFactory*           cf          = context->CF;
	Stream*                         errorStream = Journal_Register( Error_Type, self->type );

	/* HACK - Get pointer to FieldVariable in build phase just to let FieldVariables be created in plugins */
	self->fieldVariable = Stg_ComponentFactory_ConstructByName( cf, self->fieldVariableName, FieldVariable, True, data ); 
	Journal_Firewall( self->fieldVariable->fieldComponentCount == 1, errorStream,
		"Error - in %s(): provided FieldVariable \"%s\" has %u components - but %s Component "
		"can only visualise FieldVariables with 1 component. Did you mean to visualise the "
		"magnitude of the given field?\n", __func__, self->fieldVariable->name,
		self->fieldVariable->fieldComponentCount, self->type );
}

void _lucScalarFieldCrossSection_Initialise( void* drawingObject, void* data ) {}
void _lucScalarFieldCrossSection_Execute( void* drawingObject, void* data ) {}
void _lucScalarFieldCrossSection_Destroy( void* drawingObject, void* data ) {}

void _lucScalarFieldCrossSection_Setup( void* drawingObject, void* _context ) {
	lucScalarFieldCrossSection*       self            = (lucScalarFieldCrossSection*)drawingObject;

	lucColourMap_CalibrateFromFieldVariable( self->colourMap, self->fieldVariable );
	_lucOpenGLDrawingObject_Setup( self, _context );
}
	
void _lucScalarFieldCrossSection_Draw( void* drawingObject, lucWindow* window, lucViewportInfo* viewportInfo, void* _context ) {
	lucScalarFieldCrossSection*       self            = (lucScalarFieldCrossSection*)drawingObject;
	_lucOpenGLDrawingObject_Draw( self, window, viewportInfo, _context );
}


void _lucScalarFieldCrossSection_CleanUp( void* drawingObject, void* _context ) {
	lucScalarFieldCrossSection*       self            = (lucScalarFieldCrossSection*)drawingObject;
	_lucOpenGLDrawingObject_CleanUp( self, _context );
}

void _lucScalarFieldCrossSection_BuildDisplayList( void* drawingObject, void* _context ) {
	lucScalarFieldCrossSection*       self            = (lucScalarFieldCrossSection*)drawingObject;
	lucScalarFieldCrossSection_DrawCrossSection( self, self->crossSectionValue, self->crossSectionAxis );
}

#define FUDGE_FACTOR 0.0001

void lucScalarFieldCrossSection_DrawCrossSection( void* drawingObject, double crossSectionValue, Axis axis ) {
	lucScalarFieldCrossSection*       self            = (lucScalarFieldCrossSection*)drawingObject;
	FieldVariable* fieldVariable = self->fieldVariable;
	Axis           aAxis;
	Axis           bAxis;
	Coord          min;
	Coord          max;
	Coord          pos;
	Coord          interpolationCoord;
	float          normal[3];
	Index          aResolution;
	Index          bResolution;
	Index          aIndex;
	Index          bIndex;
	double         aLength;
	double         bLength;
	Dimension_Index dim_I;

	/* Ensure the field is synchronised. */
	lucOpenGLDrawingObject_SyncShadowValues( self, self->fieldVariable );
	
	glDisable(GL_LIGHTING);
	
	aAxis = ( axis == I_AXIS ? J_AXIS : I_AXIS );
	bAxis = ( axis == K_AXIS ? J_AXIS : K_AXIS );
	
	aResolution = self->resolution[ aAxis ];
	bResolution = self->resolution[ bAxis ];
	
	Journal_DPrintfL( self->debugStream, 2, 
			"%s called on field %s, with res[A] as %u, res[B] as %u, axis of cross section as %d, crossSectionValue as %g\n",
			__func__, fieldVariable->name, aResolution, bResolution, axis, crossSectionValue );
	
	FieldVariable_GetMinAndMaxLocalCoords( fieldVariable, min, max );
	/* Crop the size of the cros-section that you wish to draw */
	for ( dim_I = 0 ; dim_I < fieldVariable->dim ; dim_I++ ) {
		min[ dim_I ] = MAX( self->minCropValues[ dim_I ], min[ dim_I ]);
		max[ dim_I ] = MIN( self->maxCropValues[ dim_I ], max[ dim_I ]);
	}

	/* Create normal */
	normal[axis]  = 1.0;
	normal[aAxis] = 0.0;
	normal[bAxis] = 0.0;
	glNormal3fv( normal );

	/* Find position of cross - section */
	pos[axis] = crossSectionValue;

	aLength = (max[aAxis] - min[aAxis])/(double)aResolution;
	bLength = (max[bAxis] - min[bAxis])/(double)bResolution;

	Journal_DPrintfL( self->debugStream, 2, "Calculated aLength as %g, bLength as %g\n", aLength, bLength );

	/* Plot a number of tiles with a colour map to represent scalar field */
	/* OpenGL will interpolate colours within tile */
	for ( aIndex = 0 ; aIndex < aResolution + 1 ; aIndex++ ) {
		glBegin(GL_QUAD_STRIP);
		for ( bIndex = 0 ; bIndex < bResolution + 2 ; bIndex++ ) {
			/* Get position */
			pos[ aAxis ] = min[ aAxis ] + (double)aIndex * aLength;
			pos[ bAxis ] = min[ bAxis ] + (double)bIndex * bLength;

			memcpy( interpolationCoord, pos, sizeof(Coord) );

			if ( pos[ bAxis ] >= max[ bAxis ] ) { 
				pos[ bAxis ] = max[ bAxis ];
				interpolationCoord[ bAxis ] = max[ bAxis ] - FUDGE_FACTOR/bLength;
			}
			if (pos[ aAxis ] >= max[ aAxis ]) {
				pos[ aAxis ] = max[ aAxis ];
				interpolationCoord[ aAxis ] = max[ aAxis ] - FUDGE_FACTOR/aLength;
			}
			
			/* Plot bottom left corner of coloured tile */
			lucScalarFieldCrossSection_PlotColouredVertex( self, interpolationCoord, pos );

			/* Plot top left corner of coloured tile */
			pos[ aAxis ] += aLength;
			if (pos[ aAxis ] >= max[ aAxis ]) {
				pos[ aAxis ] = max[ aAxis ];
				interpolationCoord[ aAxis ] = max[ aAxis ] - FUDGE_FACTOR/aLength;
			}
			else
				interpolationCoord[ aAxis ] = pos[ aAxis ];
			
			lucScalarFieldCrossSection_PlotColouredVertex( self, interpolationCoord, pos );
		}
		glEnd();
	}

	glEnable(GL_LIGHTING);
}

Bool lucScalarFieldCrossSection_PlotColouredVertex( void* drawingObject, Coord interpolationCoord, Coord plotCoord ) {
	lucScalarFieldCrossSection*       self            = (lucScalarFieldCrossSection*)drawingObject;
	Bool result;
	FieldVariable* fieldVariable = self->fieldVariable;
	lucColourMap*  cmap          = self->colourMap;
	double         quantity;

	Journal_DPrintfL( self->debugStream, 3, "%s called at interpolationCoord (%g,%g,%g)  - ",
			__func__, interpolationCoord[0], interpolationCoord[1], interpolationCoord[2] );

	/* Interpolate value to this position */
	result = FieldVariable_InterpolateValueAt( fieldVariable, interpolationCoord, &quantity );
	if ( LOCAL == result || SHADOW == result ) {
		/* Get colour for this value from colour map */
		lucColourMap_SetOpenGLColourFromValue( cmap, quantity );

		Journal_DPrintfL( self->debugStream, 3, "%s is %g there.\n", fieldVariable->name, quantity );
		
		Journal_DPrintfL( self->debugStream, 3, "Plotting At Vertex (%g,%g,%g).\n", 
				plotCoord[0], plotCoord[1], plotCoord[2]  );

		/* Plot Vertex */
		glVertex3dv(plotCoord);

		return True;
	}

	Journal_DPrintfL( self->debugStream, 3, "%s NOT FOUND THERE.\n", fieldVariable->name );
	/* If value could not be interpolated return warning */
	return False;
}
