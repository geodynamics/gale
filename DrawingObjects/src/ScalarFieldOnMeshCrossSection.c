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
** $Id: ScalarFieldCrossSection.c 568 2006-06-02 06:21:50Z RobertTurnbull $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgFEM/StgFEM.h>

#ifdef GLUCIFER_USE_PICELLERATOR
	#include <StgFEM/StgFEM.h>
	#include <PICellerator/PICellerator.h>
#endif

#include <glucifer/Base/Base.h>
#include <glucifer/RenderingEngines/RenderingEngines.h>

#include "types.h"
#include "OpenGLDrawingObject.h"
#include "ScalarFieldOnMeshCrossSection.h"

#include <assert.h>
#include <gl.h>
#include <glu.h>
#include <string.h>
#include <ctype.h>

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type lucScalarFieldOnMeshCrossSection_Type = "lucScalarFieldOnMeshCrossSection";

/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
lucScalarFieldOnMeshCrossSection* _lucScalarFieldOnMeshCrossSection_New( 
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
	lucScalarFieldOnMeshCrossSection*					self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( sizeOfSelf >= sizeof(lucScalarFieldOnMeshCrossSection) );
	self = (lucScalarFieldOnMeshCrossSection*) _lucOpenGLDrawingObject_New( 
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

void _lucScalarFieldOnMeshCrossSection_Init( 
		lucScalarFieldOnMeshCrossSection*                                  self,
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

void _lucScalarFieldOnMeshCrossSection_Delete( void* drawingObject ) {
	lucScalarFieldOnMeshCrossSection*  self = (lucScalarFieldOnMeshCrossSection*)drawingObject;

	_lucOpenGLDrawingObject_Delete( self );
}

void _lucScalarFieldOnMeshCrossSection_Print( void* drawingObject, Stream* stream ) {
	lucScalarFieldOnMeshCrossSection*  self = (lucScalarFieldOnMeshCrossSection*)drawingObject;

	_lucOpenGLDrawingObject_Print( self, stream );
}

void* _lucScalarFieldOnMeshCrossSection_Copy( void* drawingObject, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap) {
	lucScalarFieldOnMeshCrossSection*  self = (lucScalarFieldOnMeshCrossSection*)drawingObject;
	lucScalarFieldOnMeshCrossSection* newDrawingObject;

	newDrawingObject = _lucOpenGLDrawingObject_Copy( self, dest, deep, nameExt, ptrMap );

	/* TODO */
	abort();

	return (void*) newDrawingObject;
}


void* _lucScalarFieldOnMeshCrossSection_DefaultNew( Name name ) {
	return (void*) _lucScalarFieldOnMeshCrossSection_New(
		sizeof(lucScalarFieldOnMeshCrossSection),
		lucScalarFieldOnMeshCrossSection_Type,
		_lucScalarFieldOnMeshCrossSection_Delete,
		_lucScalarFieldOnMeshCrossSection_Print,
		NULL,
		_lucScalarFieldOnMeshCrossSection_DefaultNew,
		_lucScalarFieldOnMeshCrossSection_Construct,
		_lucScalarFieldOnMeshCrossSection_Build,
		_lucScalarFieldOnMeshCrossSection_Initialise,
		_lucScalarFieldOnMeshCrossSection_Execute,
		_lucScalarFieldOnMeshCrossSection_Destroy,
		_lucScalarFieldOnMeshCrossSection_Setup,
		_lucScalarFieldOnMeshCrossSection_Draw,
		_lucScalarFieldOnMeshCrossSection_CleanUp,
		_lucScalarFieldOnMeshCrossSection_BuildDisplayList,
		name );
}

void _lucScalarFieldOnMeshCrossSection_Construct( void* drawingObject, Stg_ComponentFactory* cf, void* data ){
	lucScalarFieldOnMeshCrossSection*     self = (lucScalarFieldOnMeshCrossSection*)drawingObject;
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
	
	/* This variable is now constructed in the build phase 	
	fieldVariable =  Stg_ComponentFactory_ConstructByKey( cf, self->name, "FieldVariable", FieldVariable, True ) ;
	*/
	
	colourMap = Stg_ComponentFactory_ConstructByKey( cf, self->name, "ColourMap", lucColourMap, True, data );

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
	
	_lucScalarFieldOnMeshCrossSection_Init( 
			self, 
			fieldVariableName,
			colourMap,
			resolution,
			value,
			axis,
			minCropValues,
			maxCropValues );
}

void _lucScalarFieldOnMeshCrossSection_Build( void* drawingObject, void* data ) {
	lucScalarFieldOnMeshCrossSection*     self    = (lucScalarFieldOnMeshCrossSection*)drawingObject;
	AbstractContext*                context = Stg_CheckType( data, AbstractContext );
	Stg_ComponentFactory*           cf      = context->CF;

	/* HACK - Get pointer to FieldVariable in build phase just to let FieldVariables be created in plugins */
	self->fieldVariable = Stg_ComponentFactory_ConstructByName( cf, self->fieldVariableName, FieldVariable, True, 0 /* dummy */ );
}

void _lucScalarFieldOnMeshCrossSection_Initialise( void* drawingObject, void* data ) {}
void _lucScalarFieldOnMeshCrossSection_Execute( void* drawingObject, void* data ) {}
void _lucScalarFieldOnMeshCrossSection_Destroy( void* drawingObject, void* data ) {}

void _lucScalarFieldOnMeshCrossSection_Setup( void* drawingObject, void* _context ) {
	lucScalarFieldOnMeshCrossSection*       self            = (lucScalarFieldOnMeshCrossSection*)drawingObject;

	lucColourMap_CalibrateFromFieldVariable( self->colourMap, self->fieldVariable );
	_lucOpenGLDrawingObject_Setup( self, _context );
}
	
void _lucScalarFieldOnMeshCrossSection_Draw( void* drawingObject, lucWindow* window, lucViewportInfo* viewportInfo, void* _context ) {
	lucScalarFieldOnMeshCrossSection*       self            = (lucScalarFieldOnMeshCrossSection*)drawingObject;
	_lucOpenGLDrawingObject_Draw( self, window, viewportInfo, _context );
}


void _lucScalarFieldOnMeshCrossSection_CleanUp( void* drawingObject, void* _context ) {
	lucScalarFieldOnMeshCrossSection*       self            = (lucScalarFieldOnMeshCrossSection*)drawingObject;
	_lucOpenGLDrawingObject_CleanUp( self, _context );
}

void _lucScalarFieldOnMeshCrossSection_BuildDisplayList( void* drawingObject, void* _context ) {
	lucScalarFieldOnMeshCrossSection*       self            = (lucScalarFieldOnMeshCrossSection*)drawingObject;
	lucScalarFieldOnMeshCrossSection_DrawCrossSection( self, self->crossSectionValue, self->crossSectionAxis );
}

#define FUDGE_FACTOR 0.0001

Bool lucScalarFieldOnMeshCrossSection_IsInList(int nb, int *list, int nbList){
	int i;
	for (i=0; i<nbList; i++){
		if(nb == list[i])
			return True;
	}
	return False;
}

void lucScalarFieldOnMeshCrossSection_DrawCrossSection( void* drawingObject, double crossSectionValue, Axis axis ) {
	lucScalarFieldOnMeshCrossSection*       self            = (lucScalarFieldOnMeshCrossSection*)drawingObject;
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
	
	Node_LocalIndex      node_lI;
	Node_Index           elNode_I;
	Element_LocalIndex   element_lI;
	unsigned	 nIncVerts, *incVerts;
	IArray		*inc;

	
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

        /* Plots element per element from node to node in order to respect the 
	geometry of a deforming mesh */

	Mesh *mesh = (Mesh*)( ( (FeVariable*)(self->fieldVariable) )->feMesh);
  #ifdef TEST
	/* Get min max of fieldVariable... This is for test purpose */
 
	float minVal = FieldVariable_GetMinGlobalFieldMagnitude( fieldVariable );
	float maxVal =  FieldVariable_GetMaxGlobalFieldMagnitude( fieldVariable );
	printf(" minVal is %f maxval is %f\n", minVal, maxVal);
	printf("mesh->nodeLocalCount is %d\n", mesh->nodeLocalCount);
	float ninterval = (maxVal - minVal) / (mesh->nodeLocalCount);
	float myVal = minVal;

	double* displayVal = Memory_Alloc(double, mesh->nodeLocalCount, "ScalarField_mesh->nodeLocalCount");
	int i;
	int k;

	int minList[15]= {0,5,10,15,20, 2, 7, 12,17,22, 4,9,14,19,24};
	int maxList[10] = {1,6,11,16,21,3,8,13,18,23};

	
	for(i = 0; i< mesh->nodeLocalCount; i++){
	
		if( lucScalarFieldOnMeshCrossSection_IsInList(i, minList, 15) ){
			//displayVal[i] = minVal + i*ninterval;
			displayVal[i] = minVal ;
		}
		else{	
			displayVal[i] = maxVal;
		}

	}
	#endif

	
	 /* TODO: The elements could have nodes insides...... !!!*/
	#define NORMAL_LUCSCALARFIELD_BEHAVIOUR
	#ifdef NORMAL_LUCSCALARFIELD_BEHAVIOUR
	inc = IArray_New();
	for ( element_lI = 0; element_lI < Mesh_GetLocalSize( mesh, Mesh_GetDimSize( mesh ) ); element_lI++ ) {
		Mesh_GetIncidence( mesh, Mesh_GetDimSize( mesh ), element_lI, MT_VERTEX, inc );
		nIncVerts = IArray_GetSize( inc );
		incVerts = IArray_GetPtr( inc );

		/* Normal general case, the element can have whatever number of nodes */
		glBegin(GL_POLYGON);
			for ( elNode_I=0; elNode_I < nIncVerts; elNode_I++ ) {
				node_lI = incVerts[elNode_I];
				//printf("elNode_I is %d node_lI is %d  displayVal[node_lI] is %f \n", elNode_I, node_lI, displayVal[node_lI]);
			
				// Should be there when no test
				memcpy( interpolationCoord, mesh->verts[node_lI], sizeof(Coord) );
				
				// Plot the node 
				lucScalarFieldOnMeshCrossSection_PlotColouredVertex( self, interpolationCoord, mesh->verts[node_lI]);

				#ifdef TEST
				// Plot the node 
			//	lucScalarFieldOnMeshCrossSection_PlotTestColouredVertex( self, displayVal[node_lI], mesh->nodeCoord[node_lI]);
				#endif
			}
		glEnd();
	}
	#endif

     	#ifdef LUC_SCALAR_FIELD_PLOTS_SQUAD_STRIP
	Node_LocalIndex      node_0;
	Node_LocalIndex      node_1;
	Node_LocalIndex      node_2;
	Node_LocalIndex      node_3;

	/** testing the quad_strip way --- ASSUMES only 4 nodes per element */
	/* checking that there is 4 nodes per element - If yes, the quad strips can be used */
	Mesh_GetIncidence( mesh, Mesh_GetDimSize( mesh ), element_lI, MT_VERTEX, inc );
	nIncVerts = IArray_GetSize( inc );
	incVerts = IArray_GetPtr( inc );
	if( nIncVerts == 4 ){
		glBegin(GL_QUAD_STRIP);
			/* The nodes are for 10, 10 msh for instance instance 0,1,10,9 We want to display 0,9 
			and 1,10 for the quads, so display node 0,3,1,2 */
			
			node_0 = incVerts[0];
			node_1 = incVerts[3];
			node_2 = incVerts[1];
			node_3 = incVerts[2];
							
			memcpy( interpolationCoord, mesh->verts[node_0], sizeof(Coord) );
			lucScalarFieldOnMeshCrossSection_PlotColouredVertex( self, interpolationCoord, mesh->verts[node_0]);
			
			memcpy( interpolationCoord, mesh->verts[node_1], sizeof(Coord) ); 
			lucScalarFieldOnMeshCrossSection_PlotColouredVertex( self, interpolationCoord, mesh->verts[node_1]);
			
			memcpy( interpolationCoord, mesh->verts[node_2], sizeof(Coord) );
			lucScalarFieldOnMeshCrossSection_PlotColouredVertex( self, interpolationCoord, mesh->verts[node_2]);
			
			memcpy( interpolationCoord, mesh->verts[node_3], sizeof(Coord) );
			lucScalarFieldOnMeshCrossSection_PlotColouredVertex( self, interpolationCoord, mesh->verts[node_3]);
		
		glEnd();
      	}
	#endif
	glEnable(GL_LIGHTING);

	NewClass_Delete( inc );
}




Bool lucScalarFieldOnMeshCrossSection_PlotTestColouredVertex( void* drawingObject, double quantity, Coord plotCoord ) {
	lucScalarFieldOnMeshCrossSection*       self            = (lucScalarFieldOnMeshCrossSection*)drawingObject;
	FieldVariable* fieldVariable = self->fieldVariable;
	lucColourMap*  cmap          = self->colourMap;
	
	
//	if ( LOCAL == FieldVariable_InterpolateValueAt( fieldVariable, interpolationCoord, &quantity )) {
		/* Get colour for this value from colour map */
		//printf("quantity is : %f\n ", quantity);
		
		lucColourMap_SetOpenGLColourFromValue( cmap, quantity );

	//	Journal_DPrintfL( self->debugStream, 3, "%s is %g there.\n", fieldVariable->name, quantity );
		
	//	Journal_DPrintfL( self->debugStream, 3, "Plotting At Vertex (%g,%g,%g).\n", 
	//			plotCoord[0], plotCoord[1], plotCoord[2]  );

		/* Plot Vertex */
		glVertex3dv(plotCoord);

		return True;
//	}

	Journal_DPrintfL( self->debugStream, 3, "%s NOT FOUND THERE.\n", fieldVariable->name );
	/* If value could not be interpolated return warning */
	return False;
}


Bool lucScalarFieldOnMeshCrossSection_PlotColouredVertex( void* drawingObject, Coord interpolationCoord, Coord plotCoord ) {
	lucScalarFieldOnMeshCrossSection*       self            = (lucScalarFieldOnMeshCrossSection*)drawingObject;
	FieldVariable* fieldVariable = self->fieldVariable;
	lucColourMap*  cmap          = self->colourMap;
	double         quantity;

	Journal_DPrintfL( self->debugStream, 3, "%s called at interpolationCoord (%g,%g,%g)  - ",
			__func__, interpolationCoord[0], interpolationCoord[1], interpolationCoord[2] );


	/* Interpolate value to this position */
	if ( LOCAL == FieldVariable_InterpolateValueAt( fieldVariable, interpolationCoord, &quantity )) {
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
