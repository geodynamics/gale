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
** $Id: Contour.c 628 2006-10-12 08:23:07Z SteveQuenette $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#include <mpi.h>
#include <StGermain/StGermain.h>

#include <glucifer/Base/Base.h>
#include <glucifer/RenderingEngines/RenderingEngines.h>

#include "types.h"
#include "OpenGLDrawingObject.h"
#include "Contour.h"

#include <assert.h>
#include <gl.h>
#include <glu.h>
#include <string.h>

#ifndef MASTER
	#define MASTER 0
#endif

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type lucContour_Type = "lucContour";

/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
lucContour* _lucContour_New( 
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
	lucContour*					self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( sizeOfSelf >= sizeof(lucContour) );
	self = (lucContour*) _lucOpenGLDrawingObject_New( 
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

void _lucContour_Init( 
		lucContour*                                         self,
		FieldVariable*                                      fieldVariable,
		lucColourMap*                                       colourMap,
		Name                                                colourName,
		IJK                                                 resolution,
		double                                              lineWidth,
		Index                                               isovalueCount,
		double*                                             isovalueList,
		double                                              interval )
{
	self->fieldVariable = fieldVariable;
	self->colourMap     = colourMap;
	lucColour_FromString( &self->colour, colourName );
	memcpy( self->resolution, resolution, sizeof(IJK) );
	self->lineWidth = lineWidth;
	self->isovalueCount = isovalueCount;

	self->isovalueList = Memory_Alloc_Array( double, isovalueCount, "isovalue list" );
	memcpy( self->isovalueList, isovalueList, isovalueCount * sizeof(double) );
	self->interval = interval;
}

void _lucContour_Delete( void* drawingObject ) {
	lucContour*  self = (lucContour*)drawingObject;

	Memory_Free( self->isovalueList );

	_lucOpenGLDrawingObject_Delete( self );
}

void _lucContour_Print( void* drawingObject, Stream* stream ) {
	lucContour*  self = (lucContour*)drawingObject;

	_lucOpenGLDrawingObject_Print( self, stream );
}

void* _lucContour_Copy( void* drawingObject, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap) {
	lucContour*  self = (lucContour*)drawingObject;
	lucContour* newDrawingObject;

	newDrawingObject = _lucOpenGLDrawingObject_Copy( self, dest, deep, nameExt, ptrMap );

	/* TODO */
	abort();

	return (void*) newDrawingObject;
}


void* _lucContour_DefaultNew( Name name ) {
	return (void*) _lucContour_New(
		sizeof(lucContour),
		lucContour_Type,
		_lucContour_Delete,
		_lucContour_Print,
		NULL,
		_lucContour_DefaultNew,
		_lucContour_Construct,
		_lucContour_Build,
		_lucContour_Initialise,
		_lucContour_Execute,
		_lucContour_Destroy,
		_lucContour_Setup,
		_lucContour_Draw,
		_lucContour_CleanUp,
		_lucContour_BuildDisplayList,
		name );
}

void _lucContour_Construct( void* drawingObject, Stg_ComponentFactory* cf, void* data ){
	lucContour*      self = (lucContour*)drawingObject;
	Index            defaultResolution;
	FieldVariable*   fieldVariable;
	lucColourMap*    colourMap;
	IJK              resolution;

	/* Construct Parent */
	_lucOpenGLDrawingObject_Construct( self, cf, data );

	fieldVariable =  Stg_ComponentFactory_ConstructByKey( cf, self->name, "FieldVariable", FieldVariable, True,  data );
	colourMap     =  Stg_ComponentFactory_ConstructByKey( cf, self->name, "ColourMap",     lucColourMap,  False, data );

	defaultResolution = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, "resolution", 8 );
	resolution[ I_AXIS ] = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, "resolutionX", defaultResolution );
	resolution[ J_AXIS ] = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, "resolutionY", defaultResolution );
	resolution[ K_AXIS ] = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, "resolutionZ", defaultResolution );
			
	_lucContour_Init( 
			self, 
			fieldVariable,
			colourMap,
			Stg_ComponentFactory_GetString( cf, self->name, "colour", "black" ),
			resolution,
			(float) Stg_ComponentFactory_GetDouble( cf, self->name, "lineWidth", 1.0 ),
			0,
			NULL,
			Stg_ComponentFactory_GetDouble( cf, self->name, "interval", -1.0 ) ) ;
}

void _lucContour_Build( void* drawingObject, void* data ) {}
void _lucContour_Initialise( void* drawingObject, void* data ) {}
void _lucContour_Execute( void* drawingObject, void* data ) {}
void _lucContour_Destroy( void* drawingObject, void* data ) {}

void _lucContour_Setup( void* drawingObject, void* _context ) {
	lucContour*            self            = (lucContour*)drawingObject;
	_lucOpenGLDrawingObject_Setup( self, _context );
}

void _lucContour_Draw( void* drawingObject, lucWindow* window, lucViewportInfo* viewportInfo, void* _context ) {
	lucContour*            self            = (lucContour*)drawingObject;
	_lucOpenGLDrawingObject_Draw( self, window, viewportInfo, _context );
}

void _lucContour_CleanUp( void* drawingObject, void* _context ) {
	lucContour*            self            = (lucContour*)drawingObject;
	_lucOpenGLDrawingObject_CleanUp( self, _context );
}
	
void _lucContour_BuildDisplayList( void* drawingObject, void* _context ) {
	lucContour*            self            = (lucContour*)drawingObject;
	FieldVariable*         fieldVariable   = self->fieldVariable;
	double                 isovalue;
	double                 interval        = self->interval;
	double                 minIsovalue     = FieldVariable_GetMinGlobalFieldMagnitude( fieldVariable );
	double                 maxIsovalue     = FieldVariable_GetMaxGlobalFieldMagnitude( fieldVariable );
	lucColourMap*          colourMap       = self->colourMap;
	Index                  isovalue_I;
	Coord min, max;
	
	glLineWidth(self->lineWidth);

	FieldVariable_GetMinAndMaxLocalCoords( fieldVariable, min, max );
	
	lucColour_SetOpenGLColour( &self->colour );

	/* Draw static isovalues */
	for ( isovalue_I = 0 ; isovalue_I < self->isovalueCount ; isovalue_I++ ) {
		isovalue = self->isovalueList[isovalue_I];
	
		if ( colourMap )
			lucColourMap_SetOpenGLColourFromValue( colourMap, isovalue );

		lucContour_DrawContour( self, isovalue, 0.0, K_AXIS, min, max ); 
	}
	
	/* Draw isovalues at interval */
	if ( interval <= 0.0 ) 
		return;

	for ( isovalue = minIsovalue ; isovalue < maxIsovalue ; isovalue += interval ) {
		if ( colourMap )
			lucColourMap_SetOpenGLColourFromValue( colourMap, isovalue );
		lucContour_DrawContour( self, isovalue, 0.0, K_AXIS, min, max ); 
	}
}


#define LEFT   0
#define RIGHT  1
#define BOTTOM 2
#define TOP    3

void lucContour_DrawContour( 
		void*                                             drawingObject,
		double                                            isovalue,
		double                                            planeHeight,
		Axis                                              planeAxis,
		Coord                                             min,
		Coord                                             max )
{
	lucContour*            self            = (lucContour*)drawingObject;
	FieldVariable*         fieldVariable   = self->fieldVariable;
	Axis                   aAxis           = ( planeAxis == I_AXIS ? J_AXIS : I_AXIS );
	Axis                   bAxis           = ( planeAxis == K_AXIS ? J_AXIS : K_AXIS );
	unsigned int           elementType;
	unsigned int           i, j;
	Coord                  pos;
	double **              array;
	Index                  resolutionA     = self->resolution[ aAxis ];
	Index                  resolutionB     = self->resolution[ bAxis ];
	double                 dA, dB;
	
	/* Find position of cross - section */
	pos[planeAxis] = planeHeight;
	
	/* Calculate number of points in direction A and B */
	dA = (max[ aAxis ] - min[ aAxis ])/(double) (resolutionA - 1);
	dB = (max[ bAxis ] - min[ bAxis ])/(double) (resolutionB - 1);

	array = Memory_Alloc_2DArray( double , resolutionA, resolutionB, "Field Values");
	for ( i = 0, pos[ aAxis ] = min[ aAxis ] ; i < resolutionA ; i++, pos[aAxis] += dA ) {
		for ( j = 0, pos[bAxis] = min[ bAxis ] ; j < resolutionB ; j++, pos[bAxis] += dB ) {
			if (pos[aAxis] > max[ aAxis ]) 
				pos[aAxis] = max[ aAxis ];
			if (pos[bAxis] > max[ bAxis ]) 
				pos[bAxis] = max[ bAxis ];

			/* Interpolate value to point */
			FieldVariable_InterpolateValueAt( fieldVariable, pos, &array[i][j] );
		}
	}
		
	/* Initialise OpenGL stuff */
	glDisable(GL_LIGHTING);
	glEnable(GL_LINE_SMOOTH);
	glEnable(GL_BLEND);
	glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
	glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);

	glBegin(GL_LINES);

	/* Start marching rectangles */
	for ( i = 0, pos[ aAxis ] = min[ aAxis ] ; i < resolutionA - 1 ; i++, pos[ aAxis ] += dA ) {
		for ( j = 0, pos[ bAxis ] = min[ bAxis ] ; j < resolutionB - 1 ; j++, pos[ bAxis ] += dB ) {
			/* Assign a unique number to the square type from 0 to 15 */
			elementType = 0;
			if (array[i][j]     > isovalue) 	elementType += 1;
			if (array[i+1][j]   > isovalue) 	elementType += 2;
			if (array[i][j+1]   > isovalue) 	elementType += 4;
			if (array[i+1][j+1] > isovalue) 	elementType += 8;

			switch ( elementType ) {
				case 0:
					/*  @@  */
					/*  @@  */
					break;
				case 1:		
					/*  @@  */
					/*  #@  */
					lucContour_PlotPoint( LEFT,   isovalue, array[i][j], array[i+1][j], array[i][j+1], array[i+1][j+1], pos, dA, dB, planeAxis );
					lucContour_PlotPoint( BOTTOM, isovalue, array[i][j], array[i+1][j], array[i][j+1], array[i+1][j+1], pos, dA, dB, planeAxis );
					break;
				case 2:
					/*  @@  */
					/*  @#  */	
					lucContour_PlotPoint( RIGHT, isovalue, array[i][j], array[i+1][j], array[i][j+1], array[i+1][j+1], pos, dA, dB, planeAxis );
					lucContour_PlotPoint( BOTTOM, isovalue, array[i][j], array[i+1][j], array[i][j+1], array[i+1][j+1], pos, dA, dB, planeAxis );
					break;
				case 3:
					/*  @@  */
					/*  ##  */	
					lucContour_PlotPoint( LEFT, isovalue, array[i][j], array[i+1][j], array[i][j+1], array[i+1][j+1], pos, dA, dB, planeAxis );
					lucContour_PlotPoint( RIGHT, isovalue, array[i][j], array[i+1][j], array[i][j+1], array[i+1][j+1], pos, dA, dB, planeAxis );
					break;
				case 4:
					/*  #@  */
					/*  @@  */
					lucContour_PlotPoint( LEFT  , isovalue, array[i][j], array[i+1][j], array[i][j+1], array[i+1][j+1], pos, dA, dB, planeAxis );
					lucContour_PlotPoint( TOP   , isovalue, array[i][j], array[i+1][j], array[i][j+1], array[i+1][j+1], pos, dA, dB, planeAxis );
					break;
				case 5:
					/*  #@  */
					/*  #@  */
					lucContour_PlotPoint( TOP   , isovalue, array[i][j], array[i+1][j], array[i][j+1], array[i+1][j+1], pos, dA, dB, planeAxis );
					lucContour_PlotPoint( BOTTOM, isovalue, array[i][j], array[i+1][j], array[i][j+1], array[i+1][j+1], pos, dA, dB, planeAxis );
					break;
				case 6:
					/*  #@  */
					/*  @#  */
					lucContour_PlotPoint( LEFT, isovalue, array[i][j], array[i+1][j], array[i][j+1], array[i+1][j+1], pos, dA, dB, planeAxis );
					lucContour_PlotPoint( TOP , isovalue, array[i][j], array[i+1][j], array[i][j+1], array[i+1][j+1], pos, dA, dB, planeAxis );

					lucContour_PlotPoint( RIGHT , isovalue, array[i][j], array[i+1][j], array[i][j+1], array[i+1][j+1], pos, dA, dB, planeAxis );
					lucContour_PlotPoint( BOTTOM, isovalue, array[i][j], array[i+1][j], array[i][j+1], array[i+1][j+1], pos, dA, dB, planeAxis );
					break;
				case 7:
					/*  #@  */
					/*  ##  */
					lucContour_PlotPoint( TOP, isovalue, array[i][j], array[i+1][j], array[i][j+1], array[i+1][j+1], pos, dA, dB, planeAxis );
					lucContour_PlotPoint( RIGHT, isovalue, array[i][j], array[i+1][j], array[i][j+1], array[i+1][j+1], pos, dA, dB, planeAxis );
					break;
				case 8:
					/*  @#  */
					/*  @@  */
					lucContour_PlotPoint( TOP, isovalue, array[i][j], array[i+1][j], array[i][j+1], array[i+1][j+1], pos, dA, dB, planeAxis );
					lucContour_PlotPoint( RIGHT, isovalue, array[i][j], array[i+1][j], array[i][j+1], array[i+1][j+1], pos, dA, dB, planeAxis );
					break;
				case 9:
					/*  @#  */
					/*  #@  */
					lucContour_PlotPoint( TOP, isovalue, array[i][j], array[i+1][j], array[i][j+1], array[i+1][j+1], pos, dA, dB, planeAxis );
					lucContour_PlotPoint( RIGHT, isovalue, array[i][j], array[i+1][j], array[i][j+1], array[i+1][j+1], pos, dA, dB, planeAxis );

					lucContour_PlotPoint( BOTTOM, isovalue, array[i][j], array[i+1][j], array[i][j+1], array[i+1][j+1], pos, dA, dB, planeAxis );
					lucContour_PlotPoint( LEFT, isovalue, array[i][j], array[i+1][j], array[i][j+1], array[i+1][j+1], pos, dA, dB, planeAxis );
					break;
				case 10:
					/*  @#  */
					/*  @#  */
					lucContour_PlotPoint( TOP, isovalue, array[i][j], array[i+1][j], array[i][j+1], array[i+1][j+1], pos, dA, dB, planeAxis );
					lucContour_PlotPoint( BOTTOM, isovalue, array[i][j], array[i+1][j], array[i][j+1], array[i+1][j+1], pos, dA, dB, planeAxis );
					break;
				case 11:
					/*  @#  */
					/*  ##  */
					lucContour_PlotPoint( TOP, isovalue, array[i][j], array[i+1][j], array[i][j+1], array[i+1][j+1], pos, dA, dB, planeAxis );
					lucContour_PlotPoint( LEFT, isovalue, array[i][j], array[i+1][j], array[i][j+1], array[i+1][j+1], pos, dA, dB, planeAxis );
					break;
				case 12:
					/*  ##  */
					/*  @@  */
					lucContour_PlotPoint( LEFT, isovalue, array[i][j], array[i+1][j], array[i][j+1], array[i+1][j+1], pos, dA, dB, planeAxis );
					lucContour_PlotPoint( RIGHT, isovalue, array[i][j], array[i+1][j], array[i][j+1], array[i+1][j+1], pos, dA, dB, planeAxis );
					break;
				case 13:
					/*  ##  */
					/*  #@  */
					lucContour_PlotPoint( RIGHT, isovalue, array[i][j], array[i+1][j], array[i][j+1], array[i+1][j+1], pos, dA, dB, planeAxis );
					lucContour_PlotPoint( BOTTOM, isovalue, array[i][j], array[i+1][j], array[i][j+1], array[i+1][j+1], pos, dA, dB, planeAxis );
					break;
				case 14:
					/*  ##  */
					/*  @#  */
					lucContour_PlotPoint( LEFT, isovalue, array[i][j], array[i+1][j], array[i][j+1], array[i+1][j+1], pos, dA, dB, planeAxis );
					lucContour_PlotPoint( BOTTOM, isovalue, array[i][j], array[i+1][j], array[i][j+1], array[i+1][j+1], pos, dA, dB, planeAxis );
					break;
				case 15:
					/*  ##  */
					/*  ##  */
					break;
				default:
					fprintf(stderr, "In func %s: Cannot find element %d.\n", __func__, elementType );
					abort();
			}
			
		}
	}
	glEnd();
	glEnable(GL_LIGHTING);

	/* Clean up */
	Memory_Free(array);
}

void lucContour_PlotPoint( char edge, double isovalue, double leftBtm, double rightBtm, double leftTop, double rightTop , Coord pos, double dA, double dB, Axis planeAxis ) {
	Axis            aAxis           = ( planeAxis == I_AXIS ? J_AXIS : I_AXIS );
	Axis            bAxis           = ( planeAxis == K_AXIS ? J_AXIS : K_AXIS );
	Coord           vertex;
	
	vertex[planeAxis] = pos[planeAxis];

	switch (edge) {
		case BOTTOM:
			vertex[aAxis] = pos[aAxis] + dA * (isovalue - leftBtm)/(rightBtm - leftBtm) ; 
			vertex[bAxis] = pos[bAxis];
			break;
		case TOP:	
			vertex[aAxis] = pos[aAxis] + dA * (isovalue - leftTop)/(rightTop - leftTop); 
			vertex[bAxis] = pos[bAxis] + dB;
			break;
		case LEFT:	
			vertex[aAxis] = pos[aAxis];
			vertex[bAxis] = pos[bAxis] + dB * (isovalue - leftBtm)/(leftTop - leftBtm); 
			break;
		case RIGHT:	
			vertex[aAxis] = pos[aAxis] + dA;
			vertex[bAxis] = pos[bAxis] + dB * (isovalue - rightBtm)/(rightTop - rightBtm ); 
			break;
	}
	glVertex3dv(vertex);
}
