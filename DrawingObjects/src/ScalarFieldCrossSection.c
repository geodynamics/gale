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
#include "CrossSection.h"
#include "ScalarFieldCrossSection.h"

#include <assert.h>
#include <gl.h>
#include <glu.h>
#include <string.h>
#include <ctype.h>

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type lucScalarFieldCrossSection_Type = "lucScalarFieldCrossSection";

/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
lucScalarFieldCrossSection* _lucScalarFieldCrossSection_New(  LUCSCALARFIELDCROSSSECTION_DEFARGS  ) 
{
	lucScalarFieldCrossSection*					self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( _sizeOfSelf >= sizeof(lucScalarFieldCrossSection) );
	self = (lucScalarFieldCrossSection*) _lucCrossSection_New(  LUCCROSSSECTION_PASSARGS  );
	
	return self;
}

void _lucScalarFieldCrossSection_Init( 
		lucScalarFieldCrossSection*                                  self,
		lucColourMap*                                                colourMap,
		IJK                                                          resolution,
		XYZ                                                          minCropValues,
		XYZ                                                          maxCropValues,
      Bool                                                        cullFace ) 
{
	self->colourMap = colourMap;
	memcpy( self->resolution, resolution, sizeof(IJK) );
	memcpy( self->minCropValues, minCropValues, sizeof(XYZ) );
	memcpy( self->maxCropValues, maxCropValues, sizeof(XYZ) );
	self->cullFace = cullFace;
}

void _lucScalarFieldCrossSection_Delete( void* drawingObject ) {
	lucScalarFieldCrossSection*  self = (lucScalarFieldCrossSection*)drawingObject;

	_lucCrossSection_Delete( self );
}

void _lucScalarFieldCrossSection_Print( void* drawingObject, Stream* stream ) {
	lucScalarFieldCrossSection*  self = (lucScalarFieldCrossSection*)drawingObject;

	_lucCrossSection_Print( self, stream );
}

void* _lucScalarFieldCrossSection_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                                     _sizeOfSelf = sizeof(lucScalarFieldCrossSection);
	Type                                                             type = lucScalarFieldCrossSection_Type;
	Stg_Class_DeleteFunction*                                     _delete = _lucScalarFieldCrossSection_Delete;
	Stg_Class_PrintFunction*                                       _print = _lucScalarFieldCrossSection_Print;
	Stg_Class_CopyFunction*                                         _copy = NULL;
	Stg_Component_DefaultConstructorFunction*         _defaultConstructor = _lucScalarFieldCrossSection_DefaultNew;
	Stg_Component_ConstructFunction*                           _construct = _lucScalarFieldCrossSection_AssignFromXML;
	Stg_Component_BuildFunction*                                   _build = _lucScalarFieldCrossSection_Build;
	Stg_Component_InitialiseFunction*                         _initialise = _lucScalarFieldCrossSection_Initialise;
	Stg_Component_ExecuteFunction*                               _execute = _lucScalarFieldCrossSection_Execute;
	Stg_Component_DestroyFunction*                               _destroy = _lucScalarFieldCrossSection_Destroy;
	lucDrawingObject_SetupFunction*                                _setup = _lucScalarFieldCrossSection_Setup;
	lucDrawingObject_DrawFunction*                                  _draw = _lucCrossSection_Draw;
	lucDrawingObject_CleanUpFunction*                            _cleanUp = _lucOpenGLDrawingObject_CleanUp;
	lucOpenGLDrawingObject_BuildDisplayListFunction*    _buildDisplayList = _lucScalarFieldCrossSection_BuildDisplayList;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

	return (void*) _lucScalarFieldCrossSection_New(  LUCSCALARFIELDCROSSSECTION_PASSARGS  );
}

void _lucScalarFieldCrossSection_AssignFromXML( void* drawingObject, Stg_ComponentFactory* cf, void* data ){
	lucScalarFieldCrossSection*     self = (lucScalarFieldCrossSection*)drawingObject;
	lucColourMap*    colourMap;
	Index            defaultResolution;
	IJK              resolution;
	XYZ              minCropValues;
	XYZ              maxCropValues;
   Bool              cullFace;

	/* Construct Parent */
	_lucCrossSection_AssignFromXML( self, cf, data );

	colourMap = Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"ColourMap", lucColourMap, True, data  ) ;

	defaultResolution = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, (Dictionary_Entry_Key)"resolution", 128  );
	resolution[ I_AXIS ] = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, (Dictionary_Entry_Key)"resolutionX", defaultResolution  );
	resolution[ J_AXIS ] = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, (Dictionary_Entry_Key)"resolutionY", defaultResolution  );
	resolution[ K_AXIS ] = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, (Dictionary_Entry_Key)"resolutionZ", defaultResolution  );
			
	/* Get Values with which to crop the cross section */
	minCropValues[ I_AXIS ] = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"minCropX", -HUGE_VAL  );
	minCropValues[ J_AXIS ] = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"minCropY", -HUGE_VAL  );
	minCropValues[ K_AXIS ] = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"minCropZ", -HUGE_VAL  );
	maxCropValues[ I_AXIS ] = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"maxCropX", +HUGE_VAL  );
	maxCropValues[ J_AXIS ] = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"maxCropY", +HUGE_VAL  );
	maxCropValues[ K_AXIS ] = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"maxCropZ", +HUGE_VAL  );

   cullFace = Stg_ComponentFactory_GetBool( cf, self->name, (Dictionary_Entry_Key)"cullFace", True  );
	
	_lucScalarFieldCrossSection_Init( 
			self, 
			colourMap,
			resolution,
			minCropValues,
			maxCropValues,
         cullFace );
}

void _lucScalarFieldCrossSection_Build( void* drawingObject, void* data ) {
	lucScalarFieldCrossSection*     self        = (lucScalarFieldCrossSection*)drawingObject;
	Stream*                         errorStream = Journal_Register( Error_Type, (Name)self->type  );

   /* Build field variable in parent */
   _lucCrossSection_Build(self, data);

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
	
void _lucScalarFieldCrossSection_BuildDisplayList( void* drawingObject, void* _context ) {
	lucScalarFieldCrossSection*       self            = (lucScalarFieldCrossSection*)drawingObject;
	lucScalarFieldCrossSection_DrawCrossSection( self, GL_CCW );
}

#define FUDGE_FACTOR 0.0001

void lucScalarFieldCrossSection_DrawCrossSection( void* drawingObject, int direction) {
	lucScalarFieldCrossSection*       self            = (lucScalarFieldCrossSection*)drawingObject;
	FieldVariable* fieldVariable = self->fieldVariable;
	Coord          min, globalMin;
	Coord          max, globalMax;
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

   /* set polygon face winding */
	glFrontFace(direction); 
   /* Visible from both sides? */
	if ( self->cullFace )
      glEnable(GL_CULL_FACE);
	else
      glDisable(GL_CULL_FACE);

   glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	if (fieldVariable->dim == 2) 
      glDisable(GL_LIGHTING);
   else
      glEnable(GL_LIGHTING);

	aResolution = self->resolution[ self->axis1 ];
	bResolution = self->resolution[ self->axis2 ];
	
	FieldVariable_GetMinAndMaxGlobalCoords( fieldVariable, globalMin, globalMax );
	FieldVariable_GetMinAndMaxLocalCoords( fieldVariable, min, max );

	/* Crop the size of the cros-section that you wish to draw */
	for ( dim_I = 0 ; dim_I < fieldVariable->dim ; dim_I++ ) {
		min[ dim_I ] = MAX( self->minCropValues[ dim_I ], min[ dim_I ]);
		max[ dim_I ] = MIN( self->maxCropValues[ dim_I ], max[ dim_I ]);
	}

	/* Find position of cross section */
	pos[self->axis] = lucCrossSection_GetValue(self, globalMin[self->axis], globalMax[self->axis]);
	Journal_DPrintfL( self->debugStream, 2, 
			"%s called on field %s, with res[A] as %u, res[B] as %u, axis of cross section as %d, crossSection value as %g\n",
			__func__, fieldVariable->name, aResolution, bResolution, self->axis, pos[self->axis]);

	/* Create normal */
	normal[self->axis]  = 1.0;
	normal[self->axis1] = 0.0;
	normal[self->axis2] = 0.0;
   /* Flip normal to face opposite direction if cross-section past mid-point in 3d */
   if (fieldVariable->dim == 3 && pos[self->axis] >= globalMin[self->axis] + 0.5 * globalMax[self->axis] - globalMin[self->axis]) 
      normal[self->axis] = -1.0;
	glNormal3fv( normal );

	aLength = (max[self->axis1] - min[self->axis1])/(double)aResolution;
	bLength = (max[self->axis2] - min[self->axis2])/(double)bResolution;

	Journal_DPrintfL( self->debugStream, 2, "Calculated aLength as %g, bLength as %g\n", aLength, bLength );

	/* Plot a number of tiles with a colour map to represent scalar field */
	/* OpenGL will interpolate colours within tile */
	for ( aIndex = 0 ; aIndex < aResolution + 1 ; aIndex++ ) {
		glBegin(GL_QUAD_STRIP);
		for ( bIndex = 0 ; bIndex < bResolution + 2 ; bIndex++ ) {
			/* Get position */
			pos[ self->axis1 ] = min[ self->axis1 ] + (double)aIndex * aLength;
			pos[ self->axis2 ] = min[ self->axis2 ] + (double)bIndex * bLength;

			memcpy( interpolationCoord, pos, sizeof(Coord) );

			if ( pos[ self->axis2 ] >= max[ self->axis2 ] ) { 
				pos[ self->axis2 ] = max[ self->axis2 ];
				interpolationCoord[ self->axis2 ] = max[ self->axis2 ] - FUDGE_FACTOR/bLength;
			}
			if (pos[ self->axis1 ] >= max[ self->axis1 ]) {
				pos[ self->axis1 ] = max[ self->axis1 ];
				interpolationCoord[ self->axis1 ] = max[ self->axis1 ] - FUDGE_FACTOR/aLength;
			}
			
			/* Plot bottom left corner of coloured tile */
			lucScalarFieldCrossSection_PlotColouredVertex( self, interpolationCoord, pos );

			/* Plot top left corner of coloured tile */
			pos[ self->axis1 ] += aLength;
			if (pos[ self->axis1 ] >= max[ self->axis1 ]) {
				pos[ self->axis1 ] = max[ self->axis1 ];
				interpolationCoord[ self->axis1 ] = max[ self->axis1 ] - FUDGE_FACTOR/aLength;
			}
			else
				interpolationCoord[ self->axis1 ] = pos[ self->axis1 ];
			
			lucScalarFieldCrossSection_PlotColouredVertex( self, interpolationCoord, pos );
		}
		glEnd();
	}

   glEnable(GL_LIGHTING);
   glEnable(GL_CULL_FACE);
   glFrontFace(GL_CCW);
}

Bool lucScalarFieldCrossSection_PlotColouredVertex( void* drawingObject, Coord interpolationCoord, Coord plotCoord ) {
	lucScalarFieldCrossSection*       self            = (lucScalarFieldCrossSection*)drawingObject;
	Bool result;
	FieldVariable* fieldVariable = self->fieldVariable;
	lucColourMap*  cmap          = self->colourMap;
	double         quantity;

//   static double ttime1, ttime2;
//   double time2, time1 = MPI_Wtime();

	Journal_DPrintfL( self->debugStream, 3, "%s called at interpolationCoord (%g,%g,%g)  - ",
			__func__, interpolationCoord[0], interpolationCoord[1], interpolationCoord[2] );

	/* Interpolate value to this position */
	result = FieldVariable_InterpolateValueAt( fieldVariable, interpolationCoord, &quantity );
//   ttime1 += (MPI_Wtime() - time1);
//   time2 = MPI_Wtime();

	if ( LOCAL == result || SHADOW == result ) {
		/* Get colour for this value from colour map */
		lucColourMap_SetOpenGLColourFromValue( cmap, quantity );

		Journal_DPrintfL( self->debugStream, 3, "%s is %g there.\n", fieldVariable->name, quantity );
		Journal_DPrintfL( self->debugStream, 3, "Plotting At Vertex (%g,%g,%g).\n", plotCoord[0], plotCoord[1], plotCoord[2]  );

		/* Plot Vertex */
		glVertex3dv(plotCoord);

//      ttime2 += (MPI_Wtime() - time2);
//		fprintf(stderr, "Plotting At Vertex (%g,%g,%g). (ttime1 %f ttime2 %f)\n", 
//				plotCoord[0], plotCoord[1], plotCoord[2] , ttime1, ttime2);

		return True;
	}

//   ttime2 += (MPI_Wtime() - time2);
//   fprintf( stderr, "%s NOT FOUND THERE. (ttime1 %f ttime2 %f)\n", fieldVariable->name, ttime1, ttime2);

   if (self->context->nproc == 1)
	   Journal_DPrintfL( self->errorStream, 3, 
      "%s Could not interpolate value at %f,%f,%f, returned OTHER_PROC yet running on 1 processor!\n", 
      fieldVariable->name, interpolationCoord[0], interpolationCoord[1], interpolationCoord[2] );

	Journal_DPrintfL( self->debugStream, 3, "%s InterpolateValueAt: NOT FOUND.\n", fieldVariable->name );
	/* If value could not be interpolated return warning */
	return False;
}


