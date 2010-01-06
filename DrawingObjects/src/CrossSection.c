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
** $Id: CrossSection.c 791 2008-09-01 02:09:06Z JulianGiordani $
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
#include "CrossSection.h"

#include <assert.h>
#include <gl.h>
#include <glu.h>
#include <string.h>
#include <ctype.h>

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type lucCrossSection_Type = "lucCrossSection";

/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
lucCrossSection* _lucCrossSection_New(  LUCCROSSSECTION_DEFARGS  ) 
{
	lucCrossSection*					self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( _sizeOfSelf >= sizeof(lucCrossSection) );
	self = (lucCrossSection*) _lucOpenGLDrawingObject_New(  LUCOPENGLDRAWINGOBJECT_PASSARGS  );
	
	return self;
}

void _lucCrossSection_Init( 
		lucCrossSection*        self,
		Name                    colourName,
		Name                    fieldVariableName,
      double                  value, 
      Axis                    axis, 
      Bool                    interpolate)
{
	lucColour_FromString( &self->colour, colourName );
	self->value = value;
	self->axis = axis;
	self->interpolate = interpolate;
}

void _lucCrossSection_Delete( void* drawingObject ) {
	lucCrossSection*  self = (lucCrossSection*)drawingObject;

	_lucOpenGLDrawingObject_Delete( self );
}

void _lucCrossSection_Print( void* drawingObject, Stream* stream ) {
	lucCrossSection*  self = (lucCrossSection*)drawingObject;
	_lucOpenGLDrawingObject_Print( self, stream );
}

void* _lucCrossSection_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                                     _sizeOfSelf = sizeof(lucCrossSection);
	Type                                                             type = lucCrossSection_Type;
	Stg_Class_DeleteFunction*                                     _delete = _lucCrossSection_Delete;
	Stg_Class_PrintFunction*                                       _print = _lucCrossSection_Print;
	Stg_Class_CopyFunction*                                         _copy = NULL;
	Stg_Component_DefaultConstructorFunction*         _defaultConstructor = _lucCrossSection_DefaultNew;
	Stg_Component_ConstructFunction*                           _construct = _lucCrossSection_AssignFromXML;
	Stg_Component_BuildFunction*                                   _build = _lucCrossSection_Build;
	Stg_Component_InitialiseFunction*                         _initialise = _lucCrossSection_Initialise;
	Stg_Component_ExecuteFunction*                               _execute = _lucCrossSection_Execute;
	Stg_Component_DestroyFunction*                               _destroy = _lucCrossSection_Destroy;
	lucDrawingObject_SetupFunction*                                _setup = _lucOpenGLDrawingObject_Setup;
	lucDrawingObject_DrawFunction*                                  _draw = _lucCrossSection_Draw;
	lucDrawingObject_CleanUpFunction*                            _cleanUp = _lucOpenGLDrawingObject_CleanUp;
	lucOpenGLDrawingObject_BuildDisplayListFunction*    _buildDisplayList = _lucCrossSection_BuildDisplayList;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

	return (void*) _lucCrossSection_New(  LUCCROSSSECTION_PASSARGS  );
}

void _lucCrossSection_AssignFromXML( void* drawingObject, Stg_ComponentFactory* cf, void* data ){
	lucCrossSection*     self = (lucCrossSection*)drawingObject;
   Name crossSectionStr;
   char axisChar;
   char crossSectionVal[20];
   char modifierChar = ' ';
   double value = 0.0;
   Axis axis = 0;
   Bool interpolate = False;

	/* Construct Parent */
	_lucOpenGLDrawingObject_AssignFromXML( self, cf, data );

   /* Read the cross section string specification */
   crossSectionStr = Stg_ComponentFactory_GetString( cf, self->name, "crossSection", "z=min");

   /* axis=value    : draw at this exact value on axis
    * axis=min      : draw at minimum of range on axis
    * axis=max      : draw at maximum of range on axis
    * axis=value%   : draw at interpolated percentage value of range on axis
    * Axis is a single character, one of [xyzXYZ] */

   /* Parse the input string */
   if ( sscanf( crossSectionStr, "%c=%s", &axisChar, crossSectionVal ) == 2 ) 
   {
      /* Axis X/Y/Z */
      if ( toupper( axisChar ) >= 'X' )
		   axis = toupper( axisChar ) - 'X';   /* x=0 y=1 z=2 */

    	if (sscanf( crossSectionVal, "%lf%c", &value, &modifierChar) >= 1)
      {
         /* Found a numeric value  + optional modifier character */
         //fprintf(stderr, "CROSS SECTION VALUE %lf on Axis %c\n",self->value, axisChar);

         /* Interpolate cross section using percentage value */
         if (modifierChar == '%')
         {
            /* Interpolate between max and min value using provided value as percentage */
            interpolate = True;
            //fprintf(stderr, "PERCENTAGE %lf %% CROSS SECTION on Axis %c\n", self->value, axisChar);
            value *= 0.01;
         }
      }
      /* Max or Min specified? */
      else if (strcmp(crossSectionVal, "min") == 0) 
      {
         value = 0.0;
         interpolate = True;
         //fprintf(stderr, "MIN CROSS SECTION AT %lf on Axis %c\n", self->value, axisChar);
      }
      else if (strcmp(crossSectionVal, "max") == 0) 
      {
         value = 1.0;
         interpolate = True;
         //fprintf(stderr, "MAX CROSS SECTION AT %lf on Axis %c\n", self->value, axisChar);
      }
	}

	_lucCrossSection_Init( 
			self, 
			Stg_ComponentFactory_GetString( cf, self->name, "colour", "black" ),
         Stg_ComponentFactory_GetString( cf, self->name, "FieldVariable", "defaultName" ),
         value,
         axis,
         interpolate	);
}

void _lucCrossSection_Build( void* drawingObject, void* data )
{
	lucCrossSection* self    = (lucCrossSection*)drawingObject;
	AbstractContext* context = self->context;
	Stg_ComponentFactory* cf = context->CF;
	
	/* HACK - Get pointer to FieldVariable in build phase just to let FieldVariables be created in plugins */
	self->fieldVariable =  Stg_ComponentFactory_ConstructByKey( cf, self->name, "FieldVariable", FieldVariable, False, data );
 	Stg_Component_Build( self->fieldVariable, data, False );
}

void _lucCrossSection_Initialise( void* drawingObject, void* data ) {}
void _lucCrossSection_Execute( void* drawingObject, void* data ) {}
void _lucCrossSection_Destroy( void* drawingObject, void* data ) {}

void _lucCrossSection_Draw( void* drawingObject, lucWindow* window, lucViewportInfo* viewportInfo, void* _context ) {
	lucCrossSection* self = (lucCrossSection*)drawingObject;
	/* Ensure the field is synchronised. */
	lucOpenGLDrawingObject_SyncShadowValues( self, self->fieldVariable );

   /* Call parent Draw */
	_lucOpenGLDrawingObject_Draw( self, window, viewportInfo, _context );
}
	
/* Default cross-section object allows drawing a cut plane at a specified coord on any axis */
void _lucCrossSection_BuildDisplayList( void* drawingObject, void* _context ) {
	lucCrossSection*  self          = (lucCrossSection*)drawingObject;
   double   plane[4][3], min[3], max[3];
	Axis     axis, aAxis, bAxis;

	/* Get Axis Directions */
   axis = self->axis;
	aAxis = ( axis == I_AXIS ? J_AXIS : I_AXIS );
	bAxis = ( axis == K_AXIS ? J_AXIS : K_AXIS );
	
   FieldVariable_GetMinAndMaxGlobalCoords(self->fieldVariable, min, max );

   /* Fixed value on chosen axis */
   plane[0][axis] = plane[1][axis] = plane[2][axis] = plane[3][axis] = self->value;
   /* Max and min values on other axis */
   plane[0][aAxis] = min[aAxis];
   plane[1][aAxis] = min[aAxis];
   plane[2][aAxis] = max[aAxis];
   plane[3][aAxis] = max[aAxis];

   plane[0][bAxis] = max[bAxis];
   plane[1][bAxis] = min[bAxis];
   plane[2][bAxis] = min[bAxis];
   plane[3][bAxis] = max[bAxis];

	lucColour_SetOpenGLColour( &self->colour );

   glEnable(GL_LIGHTING);
   glDisable(GL_CULL_FACE);
   //glEnable(GL_CULL_FACE);

   /* Create normal aligned to axis */
   double normal[3]  = {0.0, 0.0, 0.0};
   normal[axis] = 1.0;
   glNormal3dv( normal );
   glBegin(GL_QUADS);
      glVertex3dv(plane[0]);
      glVertex3dv(plane[1]);
      glVertex3dv(plane[2]);
      glVertex3dv(plane[3]);
   glEnd();

}

/* Returns the cross section value, interpolating where necessary */
double lucCrossSection_GetValue(void* crossSection, double min, double max) 
{
   lucCrossSection* self = (lucCrossSection*)crossSection;
   if (self->interpolate)
      /* Interpolation factor 0-1 provided to determine cross-section value */
	   return min + self->value * (max - min);
   else
      /* Exact value provided */
	   return self->value;
}

/* Function to set all cross section parameters and return self for use in passing cross-sections to functions */
lucCrossSection* lucCrossSection_Set(void* crossSection, double val, Axis axis, Bool interpolate)
{
   lucCrossSection* self = (lucCrossSection*)crossSection;
   self->value = val;
   self->axis = axis;
   self->interpolate = interpolate;
   return self;
}

