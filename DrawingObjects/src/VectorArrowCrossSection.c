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
** $Id: VectorArrowCrossSection.c 791 2008-09-01 02:09:06Z JulianGiordani $
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
#include "VectorArrowCrossSection.h"

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

#ifndef MASTER
	#define MASTER 0
#endif

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type lucVectorArrowCrossSection_Type = "lucVectorArrowCrossSection";

/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
lucVectorArrowCrossSection* _lucVectorArrowCrossSection_New( 
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
	lucVectorArrowCrossSection*					self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( sizeOfSelf >= sizeof(lucVectorArrowCrossSection) );
	self = (lucVectorArrowCrossSection*) _lucOpenGLDrawingObject_New( 
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

void _lucVectorArrowCrossSection_Init( 
		lucVectorArrowCrossSection*                                  self,
		FieldVariable*                                               vectorVariable,
		Name                                                         colourName,
		IJK                                                          resolution,
		double                                                       arrowHeadSize,
		double                                                       maximum,
		Bool                                                         dynamicRange,
		double                                                       lengthScale,
		float                                                        lineWidth,
		double                                                       crossSectionValue,
		Axis                                                         crossSectionAxis ) 
{
	Stream* errorStream   = Journal_MyStream( Error_Type, self );

	self->vectorVariable = vectorVariable;
	lucColour_FromString( &self->colour, colourName );
	memcpy( self->resolution, resolution, sizeof(IJK) );
	self->arrowHeadSize = arrowHeadSize;
	Journal_Firewall( ( arrowHeadSize <= 1 && arrowHeadSize >= 0 ), errorStream,
			"Error in %s:\narrowHeadSize given for %s was not in the range [0, 1]. " 
			"Please use an arrowHeadSize within this range\n", __func__, self->name );
	self->maximum = maximum;
	self->dynamicRange = dynamicRange;
	self->lengthScale = lengthScale;
	self->lineWidth = lineWidth;
	
	self->crossSectionValue = crossSectionValue;
	self->crossSectionAxis = crossSectionAxis;
}

void _lucVectorArrowCrossSection_Delete( void* drawingObject ) {
	lucVectorArrowCrossSection*  self = (lucVectorArrowCrossSection*)drawingObject;

	_lucOpenGLDrawingObject_Delete( self );
}

void _lucVectorArrowCrossSection_Print( void* drawingObject, Stream* stream ) {
	lucVectorArrowCrossSection*  self = (lucVectorArrowCrossSection*)drawingObject;

	_lucOpenGLDrawingObject_Print( self, stream );
}

void* _lucVectorArrowCrossSection_Copy( void* drawingObject, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap) {
	lucVectorArrowCrossSection*  self = (lucVectorArrowCrossSection*)drawingObject;
	lucVectorArrowCrossSection* newDrawingObject;

	newDrawingObject = _lucOpenGLDrawingObject_Copy( self, dest, deep, nameExt, ptrMap );

	/* TODO */
	abort();

	return (void*) newDrawingObject;
}


void* _lucVectorArrowCrossSection_DefaultNew( Name name ) {
	return (void*) _lucVectorArrowCrossSection_New(
		sizeof(lucVectorArrowCrossSection),
		lucVectorArrowCrossSection_Type,
		_lucVectorArrowCrossSection_Delete,
		_lucVectorArrowCrossSection_Print,
		NULL,
		_lucVectorArrowCrossSection_DefaultNew,
		_lucVectorArrowCrossSection_Construct,
		_lucVectorArrowCrossSection_Build,
		_lucVectorArrowCrossSection_Initialise,
		_lucVectorArrowCrossSection_Execute,
		_lucVectorArrowCrossSection_Destroy,
		_lucVectorArrowCrossSection_Setup,
		_lucVectorArrowCrossSection_Draw,
		_lucVectorArrowCrossSection_CleanUp,
		_lucVectorArrowCrossSection_BuildDisplayList,
		name );
}

void _lucVectorArrowCrossSection_Construct( void* drawingObject, Stg_ComponentFactory* cf, void* data ){
	lucVectorArrowCrossSection* self = (lucVectorArrowCrossSection*)drawingObject;
	FieldVariable*   vectorVariable;
	Index            defaultResolution;
	IJK              resolution;
	double           value               = 0.0;
	Axis             axis                = 0;
	char             axisChar;
	Name             crossSectionName;

	/* Construct Parent */
	_lucOpenGLDrawingObject_Construct( self, cf, data );

	vectorVariable =  Stg_ComponentFactory_ConstructByKey( cf, self->name, "VectorVariable", FieldVariable, True, data );

	defaultResolution = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, "resolution", 8 );
	resolution[ I_AXIS ] = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, "resolutionX", defaultResolution );
	resolution[ J_AXIS ] = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, "resolutionY", defaultResolution );
	resolution[ K_AXIS ] = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, "resolutionZ", defaultResolution );
			
	crossSectionName = Stg_ComponentFactory_GetString( cf, self->name, "crossSection", "" );
	if ( sscanf( crossSectionName, "%c=%lf", &axisChar, &value ) == 2 ) {
		if ( toupper( axisChar ) >= 'X' )
			axis = toupper( axisChar ) - 'X';
	}
	
	_lucVectorArrowCrossSection_Init( 
			self, 
			vectorVariable,
			Stg_ComponentFactory_GetString( cf, self->name, "colour", "black" ),
			resolution,
			Stg_ComponentFactory_GetDouble( cf, self->name, "arrowHeadSize", 0.3 ),
			Stg_ComponentFactory_GetDouble( cf, self->name, "maximum", 1.0 ),
			Stg_ComponentFactory_GetBool( cf, self->name, "dynamicRange", True ),
			Stg_ComponentFactory_GetDouble( cf, self->name, "lengthScale", 0.3 ),
			(float) Stg_ComponentFactory_GetDouble( cf, self->name, "lineWidth", 1.0 ),
			value,
			axis );
}

void _lucVectorArrowCrossSection_Build( void* drawingObject, void* data ) {}
void _lucVectorArrowCrossSection_Initialise( void* drawingObject, void* data ) {}
void _lucVectorArrowCrossSection_Execute( void* drawingObject, void* data ) {}
void _lucVectorArrowCrossSection_Destroy( void* drawingObject, void* data ) {}

void _lucVectorArrowCrossSection_Setup( void* drawingObject, void* _context ) {
	lucVectorArrowCrossSection*       self            = (lucVectorArrowCrossSection*)drawingObject;

	_lucOpenGLDrawingObject_Setup( self, _context );
}
	
void _lucVectorArrowCrossSection_Draw( void* drawingObject, lucWindow* window, lucViewportInfo* viewportInfo, void* _context ) {
	lucVectorArrowCrossSection*       self            = (lucVectorArrowCrossSection*)drawingObject;

	_lucOpenGLDrawingObject_Draw( self, window, viewportInfo, _context );
}


void _lucVectorArrowCrossSection_CleanUp( void* drawingObject, void* _context ) {
	lucVectorArrowCrossSection*       self            = (lucVectorArrowCrossSection*)drawingObject;

	_lucOpenGLDrawingObject_CleanUp( self, _context );
}

void _lucVectorArrowCrossSection_BuildDisplayList( void* drawingObject, void* _context ) {
	lucVectorArrowCrossSection*       self            = (lucVectorArrowCrossSection*)drawingObject;
	DomainContext*            context         = (DomainContext*) _context;

	_lucVectorArrowCrossSection_DrawCrossSection( self, context->dim, self->crossSectionValue, self->crossSectionAxis );
}

void _lucVectorArrowCrossSection_DrawCrossSection( void* drawingObject, Dimension_Index dim, double crossSectionValue, Axis axis ) {
	lucVectorArrowCrossSection*  self           = (lucVectorArrowCrossSection*)drawingObject;
	FieldVariable*    vectorVariable = self->vectorVariable;
	Axis              aAxis          = (axis == I_AXIS ? J_AXIS : I_AXIS);
	Axis              bAxis          = (axis == K_AXIS ? J_AXIS : K_AXIS);
	Coord             pos;
	XYZ               vector;
	Coord             globalMin;
	Coord             globalMax;
	Coord             localMin;
	Coord             localMax;
	double            dA, dB;
	double            scaleValue;
	
	Journal_DPrintf( self->debugStream, "In %s():\n", __func__ );
	Stream_Indent( self->debugStream );

	lucOpenGLDrawingObject_SyncShadowValues( self, vectorVariable );

	if ( True == self->dynamicRange ) {
		scaleValue = 1 / FieldVariable_GetMaxGlobalFieldMagnitude( vectorVariable );
		Journal_DPrintf( self->debugStream, "Dynamic range enabled -> scale value set to "
			"1 / field maximum %g\n", scaleValue );
	}
	else {
		scaleValue = 1 / self->maximum; 
		Journal_DPrintf( self->debugStream, "Dynamic range disabled -> scale value set to "
			"1 / maximum from XML of %g\n", scaleValue );
	}

	scaleValue = scaleValue * self->lengthScale;

	Journal_DPrintf( self->debugStream, "Specified lengthScale for arrows is %.2f -> "
		"final scale value is %g\n", self->lengthScale, scaleValue );
	
	lucColour_SetOpenGLColour( &self->colour );

	FieldVariable_GetMinAndMaxGlobalCoords( vectorVariable, globalMin, globalMax );
	FieldVariable_GetMinAndMaxLocalCoords( vectorVariable, localMin, localMax );

	glLineWidth(self->lineWidth);
	
	dA = (globalMax[ aAxis ] - globalMin[ aAxis ])/(double)self->resolution[ aAxis ];
	dB = (globalMax[ bAxis ] - globalMin[ bAxis ])/(double)self->resolution[ bAxis ];
	
	pos[axis] = crossSectionValue;

	for ( pos[ aAxis ] = globalMin[ aAxis ] + dA * 0.5 ; pos[ aAxis ] < globalMax[ aAxis ] ; pos[ aAxis ] += dA ) {
		for ( pos[ bAxis ] = globalMin[ bAxis ] + dB * 0.5 ; pos[ bAxis ] < globalMax[ bAxis ] ; pos[ bAxis ] += dB ) {

			if ( pos[ aAxis ] < localMin[ aAxis ] || pos[ aAxis ] >= localMax[ aAxis ] )
				continue;
			if ( pos[ bAxis ] < localMin[ bAxis ] || pos[ bAxis ] >= localMax[ bAxis ] )
				continue;

			/* Get Value of Vector */
			if ( FieldVariable_InterpolateValueAt( vectorVariable, pos, vector ) == LOCAL ) {
				/* Draw vector */
				luc_DrawVector( dim, pos, vector, scaleValue, self->arrowHeadSize );
			}
		}
	}
	Stream_UnIndent( self->debugStream );
}
