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
#include <StgDomain/StgDomain.h>
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
		lucScalarFieldOnMeshCrossSection*                            self,
		Name                                                         fieldVariableName,
		lucColourMap*                                                colourMap,
		Node_Index                                                   crossSection_I,
		Axis                                                         crossSectionAxis,
		XYZ                                                          minCropValues,
		XYZ                                                          maxCropValues ) 
{
//	self->fieldVariable = fieldVariable;
	self->fieldVariableName = fieldVariableName;
	self->colourMap = colourMap;
	self->crossSection_I = crossSection_I;
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
		_lucScalarFieldOnMeshCrossSection_AssignFromXML,
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

void _lucScalarFieldOnMeshCrossSection_AssignFromXML( void* drawingObject, Stg_ComponentFactory* cf, void* data ){
	lucScalarFieldOnMeshCrossSection*     self = (lucScalarFieldOnMeshCrossSection*)drawingObject;
	lucColourMap*    colourMap;
	Index            defaultResolution;
	char             axisChar;
	Node_Index       value               = 0;
	Axis             axis                = 0;
	Name             crossSectionName;
	Name             fieldVariableName;
	XYZ              minCropValues;
	XYZ              maxCropValues;

	/* Construct Parent */
	_lucOpenGLDrawingObject_AssignFromXML( self, cf, data );

	fieldVariableName = Stg_ComponentFactory_GetString( cf, self->name, "FieldVariable", "defaultName" );
	
	/* This variable is now constructed in the build phase 	
	fieldVariable =  Stg_ComponentFactory_ConstructByKey( cf, self->name, "FieldVariable", FieldVariable, True ) ;
	*/
	
	colourMap = Stg_ComponentFactory_ConstructByKey( cf, self->name, "ColourMap", lucColourMap, True, data );

	crossSectionName = Stg_ComponentFactory_GetString( cf, self->name, "crossSection", "" );
	if ( sscanf( crossSectionName, "%c=%d", &axisChar, &value ) == 2 ) {
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
			value,
			axis,
			minCropValues,
			maxCropValues );
}

void _lucScalarFieldOnMeshCrossSection_Build( void* drawingObject, void* data ) {
	lucScalarFieldOnMeshCrossSection*     self    = (lucScalarFieldOnMeshCrossSection*)drawingObject;
	AbstractContext*                context = Stg_CheckType( data, AbstractContext );
	Stg_ComponentFactory*           cf      = context->CF;
	FeVariable*                     feVariable;
	Mesh*                           mesh;
	Stream*                         errorStream = Journal_Register( Error_Type, self->type );
	

	/* HACK - Get pointer to FieldVariable in build phase just to let FieldVariables be created in plugins */
	feVariable =  Stg_ComponentFactory_ConstructByName( cf, self->fieldVariableName, FeVariable, True, 0 /* dummy */ );
	self->fieldVariable = (FieldVariable*) feVariable;

	Journal_Firewall( self->fieldVariable->fieldComponentCount == 1, errorStream,
		"Error - in %s(): provided FieldVariable \"%s\" has %u components - but %s Component "
		"can only visualise FieldVariables with 1 component. Did you mean to visualise the "
		"magnitude of the given field?\n", __func__, self->fieldVariable->name,
		self->fieldVariable->fieldComponentCount, self->type );

	Stg_Component_Build( feVariable, data, False );
	mesh    = (Mesh*) feVariable->feMesh;

	/* Store the Vertex Grid */
	self->vertexGridHandle = ExtensionManager_GetHandle( mesh->info, "vertexGrid" );
	if ( self->vertexGridHandle == (ExtensionInfo_Index)-1 )

	Journal_Firewall( self->vertexGridHandle != (ExtensionInfo_Index)-1, errorStream,
		"Error - in %s(): provided FieldVariable \"%s\" doesn't have a Vertex Grid.\n"
		"Try visualising with lucScalarField instead.\n", __func__, self->fieldVariable->name );
		
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
	lucScalarFieldOnMeshCrossSection_DrawCrossSection( self, self->crossSection_I, self->crossSectionAxis );
}

void lucScalarFieldOnMeshCrossSection_DrawCrossSection( void* drawingObject, Node_LocalIndex crossSection_I, Axis axis ) {
	lucScalarFieldOnMeshCrossSection*       self            = (lucScalarFieldOnMeshCrossSection*)drawingObject;
	FeVariable*          fieldVariable = (FeVariable*) self->fieldVariable;
	Mesh*                mesh          = (Mesh*) fieldVariable->feMesh;
	Axis                 aAxis;
	Axis                 bAxis;
	Grid*                vertGrid;
	IJK                  node_ijk;
	float                normal[3];
	Dimension_Index      dim_I;
	Node_GlobalIndex     node_gI;
	Node_DomainIndex     node_dI_1, node_dI_2;
	Node_DomainIndex     nDomainNodes;
	unsigned	     nIncVerts, *incVerts;
	IArray*              inc;

	glDisable(GL_LIGHTING);
	
	/* Get Axis Directions */
	aAxis = ( axis == I_AXIS ? J_AXIS : I_AXIS );
	bAxis = ( axis == K_AXIS ? J_AXIS : K_AXIS );
	
	vertGrid = *(Grid**)ExtensionManager_Get( mesh->info, mesh, self->vertexGridHandle );
	
	Journal_DPrintfL( self->debugStream, 2, 
			"%s called on field %s, with axis of cross section as %d, crossSection_I as %d\n",
			__func__, fieldVariable->name, axis, crossSection_I );
	
	/* Create normal */
	normal[axis]  = 1.0;
	normal[aAxis] = 0.0;
	normal[bAxis] = 0.0;
	glNormal3fv( normal );

	nDomainNodes = Mesh_GetDomainSize( mesh, MT_VERTEX );

	/* Plots cross section */
	node_ijk[ axis ] = crossSection_I;
	for ( node_ijk[ aAxis ] = 0 ; node_ijk[ aAxis ] < vertGrid->sizes[ aAxis ] - 1 ; node_ijk[ aAxis ]++ ) {
		glBegin(GL_QUAD_STRIP);
		for ( node_ijk[ bAxis ] = 0 ; node_ijk[ bAxis ] < vertGrid->sizes[ bAxis ] ; node_ijk[ bAxis ]++ ) {
			node_gI = Grid_Project( vertGrid, node_ijk );
			/* Get Local Node Index */
			if( !Mesh_GlobalToDomain( mesh, MT_VERTEX, node_gI, &node_dI_1 ) || node_dI_1 >= nDomainNodes ){
				continue;
			}
			
			node_ijk[ aAxis ]++;
			node_gI = Grid_Project( vertGrid, node_ijk );
			/* Get Local Node Index */
			if( !Mesh_GlobalToDomain( mesh, MT_VERTEX, node_gI, &node_dI_2 ) || node_dI_2 >= nDomainNodes ){
				continue;
			}
			lucScalarFieldOnMeshCrossSection_PlotColouredNode( self, node_dI_1 );
			lucScalarFieldOnMeshCrossSection_PlotColouredNode( self, node_dI_2 );
			node_ijk[ aAxis ]--;

			/* TODO Cropping */
		}
		glEnd();
	}
	glEnable(GL_LIGHTING);
}

void lucScalarFieldOnMeshCrossSection_PlotColouredNode( void* drawingObject, Node_LocalIndex lNode_I ) {
	lucScalarFieldOnMeshCrossSection*       self            = (lucScalarFieldOnMeshCrossSection*)drawingObject;
	FeVariable*    fieldVariable = (FeVariable*) self->fieldVariable;
	lucColourMap*  cmap          = self->colourMap;
	double         quantity;

	/* Get colour for vertex */
	FeVariable_GetValueAtNode( fieldVariable, lNode_I, &quantity );
	lucColourMap_SetOpenGLColourFromValue( cmap, quantity );
	
	/* Plot vertex */
	if ( fieldVariable->dim == 2 )
		glVertex2dv( fieldVariable->feMesh->verts[lNode_I] );
	else 
		glVertex3dv( fieldVariable->feMesh->verts[lNode_I] );
}
