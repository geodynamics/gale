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
** $Id: EigenvectorsCrossSection.c 791 2008-09-01 02:09:06Z JulianGiordani $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>

#include <glucifer/Base/Base.h>
#include <glucifer/RenderingEngines/RenderingEngines.h>

#include "types.h"
#include "OpenGLDrawingObject.h"
#include "EigenvectorsCrossSection.h"

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
const Type lucEigenvectorsCrossSection_Type = "lucEigenvectorsCrossSection";

/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
lucEigenvectorsCrossSection* _lucEigenvectorsCrossSection_New( 
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
	lucEigenvectorsCrossSection*					self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( sizeOfSelf >= sizeof(lucEigenvectorsCrossSection) );
	self = (lucEigenvectorsCrossSection*) _lucOpenGLDrawingObject_New( 
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

void _lucEigenvectorsCrossSection_Init( 
		lucEigenvectorsCrossSection*                                 self,
		FieldVariable*                                               tensorField,
		Dimension_Index                                              dim,
		Name                                                         leastColourName,
		Name                                                         middleColourName,
		Name                                                         greatestColourName,
		IJK                                                          resolution,
		double                                                       arrowHeadSize,
		double                                                       lengthScale,
		float                                                        lineWidth,
		Bool 				                             useEigenValue,
		double                                                       notEigenValue,
		Bool 							     plotEigenVector,
		Bool 							     plotEigenValue,
		double 							     scaleEigenValue,
		Name                                                         leastColourForNegativeName,
		Name                                                         middleColourForNegativeName,
		Name                                                         greatestColourForNegativeName,
		double                                                       crossSectionValue,
		Axis                                                         crossSectionAxis ) 
{
	Stream* errorStream         = Journal_MyStream( Error_Type, self );
	self->tensorField = tensorField;
	if ( dim == 2 ) {
		lucColour_FromString( &self->colour[0], leastColourName );
		lucColour_FromString( &self->colour[1], greatestColourName );
		lucColour_FromString( &self->colourForNegative[0], leastColourForNegativeName );
		lucColour_FromString( &self->colourForNegative[1], greatestColourForNegativeName );

	}
	else {
		lucColour_FromString( &self->colour[0], leastColourName );
		lucColour_FromString( &self->colour[1], middleColourName );
		lucColour_FromString( &self->colour[2], greatestColourName );

	        lucColour_FromString( &self->colourForNegative[0], leastColourForNegativeName );
		lucColour_FromString( &self->colourForNegative[1], middleColourForNegativeName );
		lucColour_FromString( &self->colourForNegative[2], greatestColourForNegativeName );
	}

	memcpy( self->resolution, resolution, sizeof(IJK) );
	self->arrowHeadSize = arrowHeadSize;
	Journal_Firewall( ( arrowHeadSize <= 1 && arrowHeadSize >= 0 ), errorStream,
			"Error in %s:\narrowHeadSize given for %s was not in the range [0, 1]. " 
			"Please use an arrowHeadSize within this range\n", __func__, self->name );
	self->lengthScale = lengthScale;
	self->lineWidth = lineWidth;
	
	
	self->useEigenValue = useEigenValue;
	self->notEigenValue = notEigenValue;

	self->plotEigenVector = plotEigenVector;
	self->plotEigenValue = plotEigenValue;
	self->scaleEigenValue = scaleEigenValue;

	self->crossSectionValue = crossSectionValue;
	self->crossSectionAxis = crossSectionAxis;
}

void _lucEigenvectorsCrossSection_Delete( void* drawingObject ) {
	lucEigenvectorsCrossSection*  self = (lucEigenvectorsCrossSection*)drawingObject;

	_lucOpenGLDrawingObject_Delete( self );
}

void _lucEigenvectorsCrossSection_Print( void* drawingObject, Stream* stream ) {
	lucEigenvectorsCrossSection*  self = (lucEigenvectorsCrossSection*)drawingObject;

	_lucOpenGLDrawingObject_Print( self, stream );
}

void* _lucEigenvectorsCrossSection_Copy( void* drawingObject, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap) {
	lucEigenvectorsCrossSection*  self = (lucEigenvectorsCrossSection*)drawingObject;
	lucEigenvectorsCrossSection* newDrawingObject;

	newDrawingObject = _lucOpenGLDrawingObject_Copy( self, dest, deep, nameExt, ptrMap );

	/* TODO */
	abort();

	return (void*) newDrawingObject;
}


void* _lucEigenvectorsCrossSection_DefaultNew( Name name ) {
	return (void*) _lucEigenvectorsCrossSection_New(
		sizeof(lucEigenvectorsCrossSection),
		lucEigenvectorsCrossSection_Type,
		_lucEigenvectorsCrossSection_Delete,
		_lucEigenvectorsCrossSection_Print,
		NULL,
		_lucEigenvectorsCrossSection_DefaultNew,
		_lucEigenvectorsCrossSection_Construct,
		_lucEigenvectorsCrossSection_Build,
		_lucEigenvectorsCrossSection_Initialise,
		_lucEigenvectorsCrossSection_Execute,
		_lucEigenvectorsCrossSection_Destroy,
		_lucEigenvectorsCrossSection_Setup,
		_lucEigenvectorsCrossSection_Draw,
		_lucEigenvectorsCrossSection_CleanUp,
		_lucEigenvectorsCrossSection_BuildDisplayList,
		name );
}

void _lucEigenvectorsCrossSection_Construct( void* drawingObject, Stg_ComponentFactory* cf, void* data ){
	lucEigenvectorsCrossSection* self = (lucEigenvectorsCrossSection*)drawingObject;
	FieldVariable*   tensorField;
	Index            defaultResolution;
	IJK              resolution;
	double           value               = 0.0;
	Axis             axis                = 0;
	char             axisChar;
	Name             crossSectionName;

	/* Construct Parent */
	_lucOpenGLDrawingObject_Construct( self, cf, data );

	tensorField =  Stg_ComponentFactory_ConstructByKey( cf, self->name, "TensorField", FieldVariable, True, data );

	defaultResolution = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, "resolution", 8 );
	resolution[ I_AXIS ] = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, "resolutionX", defaultResolution );
	resolution[ J_AXIS ] = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, "resolutionY", defaultResolution );
	resolution[ K_AXIS ] = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, "resolutionZ", defaultResolution );
			
	crossSectionName = Stg_ComponentFactory_GetString( cf, self->name, "crossSection", "" );
	if ( sscanf( crossSectionName, "%c=%lf", &axisChar, &value ) == 2 ) {
		if ( toupper( axisChar ) >= 'X' )
			axis = toupper( axisChar ) - 'X';
	}
	
	_lucEigenvectorsCrossSection_Init( 
			self, 
			tensorField,
			Stg_ComponentFactory_GetRootDictUnsignedInt( cf, "dim", 2 ),
			Stg_ComponentFactory_GetString( cf, self->name, "leastColour", "black" ),
			Stg_ComponentFactory_GetString( cf, self->name, "middleColour", "black" ),
			Stg_ComponentFactory_GetString( cf, self->name, "greatestColour", "black" ),
			resolution,
			Stg_ComponentFactory_GetDouble( cf, self->name, "arrowHeadSize", 0.3 ),
			Stg_ComponentFactory_GetDouble( cf, self->name, "lengthScale", 0.3 ),
			(float) Stg_ComponentFactory_GetDouble( cf, self->name, "lineWidth", 1.0 ),
			Stg_ComponentFactory_GetBool( cf, self->name, "useEigenValue", True ),
			Stg_ComponentFactory_GetDouble( cf, self->name, "notEigenValue", 0.3 ),
			Stg_ComponentFactory_GetBool( cf, self->name, "plotEigenVector", True ),
			Stg_ComponentFactory_GetBool( cf, self->name, "plotEigenValue", False ),
			Stg_ComponentFactory_GetDouble( cf, self->name, "scaleEigenValue", 1.0 ),
			Stg_ComponentFactory_GetString( cf, self->name, "leastColourForNegative", "black" ),
			Stg_ComponentFactory_GetString( cf, self->name, "middleColourForNegative", "black" ),
			Stg_ComponentFactory_GetString( cf, self->name, "greatestColourForNegative", "black" ),
			value, 
			axis );
}

void _lucEigenvectorsCrossSection_Build( void* drawingObject, void* data ) {}
void _lucEigenvectorsCrossSection_Initialise( void* drawingObject, void* data ) {}
void _lucEigenvectorsCrossSection_Execute( void* drawingObject, void* data ) {}
void _lucEigenvectorsCrossSection_Destroy( void* drawingObject, void* data ) {}

void _lucEigenvectorsCrossSection_Setup( void* drawingObject, void* _context ) {
	lucEigenvectorsCrossSection*       self            = (lucEigenvectorsCrossSection*)drawingObject;

	_lucOpenGLDrawingObject_Setup( self, _context );
}
	
void _lucEigenvectorsCrossSection_Draw( void* drawingObject, lucWindow* window, lucViewportInfo* viewportInfo, void* _context ) {
	lucEigenvectorsCrossSection*       self            = (lucEigenvectorsCrossSection*)drawingObject;

	_lucOpenGLDrawingObject_Draw( self, window, viewportInfo, _context );
}


void _lucEigenvectorsCrossSection_CleanUp( void* drawingObject, void* _context ) {
	lucEigenvectorsCrossSection*       self            = (lucEigenvectorsCrossSection*)drawingObject;

	_lucOpenGLDrawingObject_CleanUp( self, _context );
}

void _lucEigenvectorsCrossSection_BuildDisplayList( void* drawingObject, void* _context ) {
	lucEigenvectorsCrossSection*       self            = (lucEigenvectorsCrossSection*)drawingObject;
	DomainContext*            context         = (DomainContext*) _context;

	_lucEigenvectorsCrossSection_DrawCrossSection( self, context->dim, self->crossSectionValue, self->crossSectionAxis );
}

void _lucEigenvectorsCrossSection_DrawCrossSection( void* drawingObject, Dimension_Index dim, double crossSectionValue, Axis axis ) {
	lucEigenvectorsCrossSection*  self           = (lucEigenvectorsCrossSection*)drawingObject;
	FieldVariable*    tensorField    = self->tensorField;
	Axis              aAxis          = (axis == I_AXIS ? J_AXIS : I_AXIS);
	Axis              bAxis          = (axis == K_AXIS ? J_AXIS : K_AXIS);
	Coord             pos;
	SymmetricTensor   tensor;
	Coord             globalMin;
	Coord             globalMax;
	Coord             localMin;
	Coord             localMax;
	double            dA, dB;
	Eigenvector       eigenvectorList[3];
	Dimension_Index   dim_I;
	
	FieldVariable_GetMinAndMaxGlobalCoords( tensorField, globalMin, globalMax );
	FieldVariable_GetMinAndMaxLocalCoords( tensorField, localMin, localMax );

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

			/* Get Value of Tensor at this point in space */
			if ( FieldVariable_InterpolateValueAt( tensorField, pos, tensor ) == LOCAL ) {
				SymmetricTensor_CalcAllEigenvectors( tensor, dim, eigenvectorList );

                                if(self->plotEigenVector){
					for ( dim_I = 0 ; dim_I < dim ; dim_I++ ) {
						
						lucColour_SetOpenGLColour( &self->colour[ dim_I ] );
						if(self->useEigenValue){
						     luc_DrawVector( dim, pos, eigenvectorList[ dim_I ].vector, 
								(eigenvectorList[ dim_I ].eigenvalue * self->scaleEigenValue), self->arrowHeadSize );
						}
						else{
							luc_DrawVector( dim, pos, eigenvectorList[ dim_I ].vector,
									 self->notEigenValue, self->arrowHeadSize );
						}
					}
				}
				if(self->plotEigenValue){
					GLfloat pointSize = 0;

					for ( dim_I = 0 ; dim_I < dim ; dim_I++ ) {
					        /* The EigenValue can be negative.... Got to attribute a potential */
						/* colour for negative values, one for each dim as well */
						if ( eigenvectorList[ dim_I ].eigenvalue >= 0) {
						        pointSize = eigenvectorList[ dim_I ].eigenvalue * self->scaleEigenValue;
							lucColour_SetOpenGLColour( &self->colour[ dim_I ] );
	                                        }
						else {
						        lucColour_SetOpenGLColour( &self->colourForNegative[ dim_I ] );
							pointSize = - eigenvectorList[ dim_I ].eigenvalue * self->scaleEigenValue;
	                                        }
						glPointSize( pointSize );
					
						glBegin(GL_POINTS);
							if (dim == 2)
					      			glVertex3f( (GLfloat)pos[0], (GLfloat)pos[1], 0.001 );
							else 
						    		glVertex3dv( pos );
						glEnd();
					
					}
				}

			}
		}
	}
}
