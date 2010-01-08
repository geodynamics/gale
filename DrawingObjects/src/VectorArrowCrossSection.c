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
#include "CrossSection.h"
#include "VectorArrowCrossSection.h"

#include <assert.h>
#include <gl.h>
#include <glu.h>
#include <string.h>
#include <ctype.h>

#ifndef MASTER
	#define MASTER 0
#endif

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type lucVectorArrowCrossSection_Type = "lucVectorArrowCrossSection";

/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
lucVectorArrowCrossSection* _lucVectorArrowCrossSection_New(  LUCVECTORARROWCROSSSECTION_DEFARGS  ) 
{
	lucVectorArrowCrossSection*					self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( _sizeOfSelf >= sizeof(lucVectorArrowCrossSection) );
	self = (lucVectorArrowCrossSection*) _lucCrossSection_New(  LUCCROSSSECTION_PASSARGS  );
	
	return self;
}

void _lucVectorArrowCrossSection_Init( 
		lucVectorArrowCrossSection*                                  self,
		IJK                                                          resolution,
		double                                                       arrowHeadSize,
		double                                                       maximum,
		Bool                                                         dynamicRange,
		double                                                       lengthScale,
		float                                                        lineWidth)
{
	Stream* errorStream   = Journal_MyStream( Error_Type, self );

	memcpy( self->resolution, resolution, sizeof(IJK) );
	self->arrowHeadSize = arrowHeadSize;
	Journal_Firewall( ( arrowHeadSize <= 1 && arrowHeadSize >= 0 ), errorStream,
			"Error in %s:\narrowHeadSize given for %s was not in the range [0, 1]. " 
			"Please use an arrowHeadSize within this range\n", __func__, self->name );
	self->maximum = maximum;
	self->dynamicRange = dynamicRange;
	self->lengthScale = lengthScale;
	self->lineWidth = lineWidth;
}

void _lucVectorArrowCrossSection_Delete( void* drawingObject ) {
	lucVectorArrowCrossSection*  self = (lucVectorArrowCrossSection*)drawingObject;
	_lucCrossSection_Delete( self );
}

void _lucVectorArrowCrossSection_Print( void* drawingObject, Stream* stream ) {
	lucVectorArrowCrossSection*  self = (lucVectorArrowCrossSection*)drawingObject;

	_lucCrossSection_Print( self, stream );
}

void* _lucVectorArrowCrossSection_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                                     _sizeOfSelf = sizeof(lucVectorArrowCrossSection);
	Type                                                             type = lucVectorArrowCrossSection_Type;
	Stg_Class_DeleteFunction*                                     _delete = _lucVectorArrowCrossSection_Delete;
	Stg_Class_PrintFunction*                                       _print = _lucVectorArrowCrossSection_Print;
	Stg_Class_CopyFunction*                                         _copy = NULL;
	Stg_Component_DefaultConstructorFunction*         _defaultConstructor = _lucVectorArrowCrossSection_DefaultNew;
	Stg_Component_ConstructFunction*                           _construct = _lucVectorArrowCrossSection_AssignFromXML;
	Stg_Component_BuildFunction*                                   _build = _lucVectorArrowCrossSection_Build;
	Stg_Component_InitialiseFunction*                         _initialise = _lucVectorArrowCrossSection_Initialise;
	Stg_Component_ExecuteFunction*                               _execute = _lucVectorArrowCrossSection_Execute;
	Stg_Component_DestroyFunction*                               _destroy = _lucVectorArrowCrossSection_Destroy;
	lucDrawingObject_SetupFunction*                                _setup = _lucOpenGLDrawingObject_Setup;
	lucDrawingObject_DrawFunction*                                  _draw = _lucCrossSection_Draw;
	lucDrawingObject_CleanUpFunction*                            _cleanUp = _lucOpenGLDrawingObject_CleanUp;
	lucOpenGLDrawingObject_BuildDisplayListFunction*    _buildDisplayList = _lucVectorArrowCrossSection_BuildDisplayList;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

	return (void*) _lucVectorArrowCrossSection_New(  LUCVECTORARROWCROSSSECTION_PASSARGS  );
}

void _lucVectorArrowCrossSection_AssignFromXML( void* drawingObject, Stg_ComponentFactory* cf, void* data ){
	lucVectorArrowCrossSection* self = (lucVectorArrowCrossSection*)drawingObject;
	Index            defaultResolution;
	IJK              resolution;

	/* Construct Parent */
	_lucCrossSection_AssignFromXML( self, cf, data );
   strcpy(self->fieldVariableName, "VectorVariable");

	defaultResolution = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, "resolution", 8 );
	resolution[ I_AXIS ] = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, "resolutionX", defaultResolution );
	resolution[ J_AXIS ] = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, "resolutionY", defaultResolution );
	resolution[ K_AXIS ] = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, "resolutionZ", defaultResolution );
			
	_lucVectorArrowCrossSection_Init( 
			self, 
			resolution,
			Stg_ComponentFactory_GetDouble( cf, self->name, "arrowHeadSize", 0.3 ),
			Stg_ComponentFactory_GetDouble( cf, self->name, "maximum", 1.0 ),
			Stg_ComponentFactory_GetBool( cf, self->name, "dynamicRange", True ),
			Stg_ComponentFactory_GetDouble( cf, self->name, "lengthScale", 0.3 ),
			(float) Stg_ComponentFactory_GetDouble( cf, self->name, "lineWidth", 1.0 ));
}

void _lucVectorArrowCrossSection_Build( void* drawingObject, void* data ) {
   /* Build field variable in parent */
   _lucCrossSection_Build(drawingObject, data);
}

void _lucVectorArrowCrossSection_Initialise( void* drawingObject, void* data ) {}
void _lucVectorArrowCrossSection_Execute( void* drawingObject, void* data ) {}
void _lucVectorArrowCrossSection_Destroy( void* drawingObject, void* data ) {}

void _lucVectorArrowCrossSection_BuildDisplayList( void* drawingObject, void* _context ) {
	lucVectorArrowCrossSection*       self            = (lucVectorArrowCrossSection*)drawingObject;
	DomainContext*            context         = (DomainContext*) _context;

	_lucVectorArrowCrossSection_DrawCrossSection( self, context->dim );
}

void _lucVectorArrowCrossSection_DrawCrossSection( void* drawingObject, Dimension_Index dim ) {
	lucVectorArrowCrossSection*  self           = (lucVectorArrowCrossSection*)drawingObject;
	FieldVariable*    vectorVariable = self->fieldVariable;
	Coord             pos;
	XYZ               vector;
	Coord             globalMin;
	Coord             globalMax;
	Coord             localMin;
	Coord             localMax;
	double            dA, dB;
	double            scaleValue;
	Stream*           errorStream = Journal_Register( Error_Type, self->type );
	
	Journal_DPrintf( self->debugStream, "In %s():\n", __func__ );
	Stream_Indent( self->debugStream );

	Journal_Firewall( vectorVariable->fieldComponentCount == vectorVariable->dim, errorStream,
		"Error - in %s(): provided FieldVariable \"%s\" has %u components - but %s Component "
		"can only visualse FieldVariables with %d components.\n", __func__, vectorVariable->name,
		vectorVariable->fieldComponentCount, self->type, vectorVariable->dim );

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
	
	dA = (globalMax[ self->axis1 ] - globalMin[ self->axis1 ])/(double)self->resolution[ self->axis1 ];
	dB = (globalMax[ self->axis2 ] - globalMin[ self->axis2 ])/(double)self->resolution[ self->axis2 ];
	
	pos[self->axis] = lucCrossSection_GetValue(self, globalMin[self->axis], globalMax[self->axis]);
	Journal_DPrintf( self->debugStream, "-- Drawing cross section on axis %d at value %lf\n", self->axis, pos[self->axis]);

	for ( pos[ self->axis1 ] = globalMin[ self->axis1 ] + dA * 0.5 ; pos[ self->axis1 ] < globalMax[ self->axis1 ] ; pos[ self->axis1 ] += dA ) {
		for ( pos[ self->axis2 ] = globalMin[ self->axis2 ] + dB * 0.5 ; pos[ self->axis2 ] < globalMax[ self->axis2 ] ; pos[ self->axis2 ] += dB ) {

			if ( pos[ self->axis1 ] < localMin[ self->axis1 ] || pos[ self->axis1 ] >= localMax[ self->axis1 ] )
				continue;
			if ( pos[ self->axis2 ] < localMin[ self->axis2 ] || pos[ self->axis2 ] >= localMax[ self->axis2 ] )
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


