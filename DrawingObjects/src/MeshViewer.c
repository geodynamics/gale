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
** $Id: Arrhenius.c 78 2005-11-29 11:58:21Z RobertTurnbull $
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
#include "MeshViewer.h"

#include <assert.h>
#ifdef HAVE_OPENGL_FRAMEWORK
	#include <OpenGL/gl.h>
	#include <OpenGL/glu.h>
#else
	#include <gl.h>
	#include <glu.h>
#endif
#include <string.h>


/* Textual name of this class - This is a global pointer which is used for 
   times when you need to refer to class and not a particular instance of a class */
const Type lucMeshViewer_Type = "lucMeshViewer";


/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
lucMeshViewer* _lucMeshViewer_New( 
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
	lucMeshViewer*					self;

	/* Call private constructor of parent - this will set virtual functions of 
	   parent and continue up the hierarchy tree. At the beginning of the tree 
	   it will allocate memory of the size of object and initialise all the 
	   memory to zero. */
	assert( sizeOfSelf >= sizeof(lucMeshViewer) );
	self = (lucMeshViewer*) _lucOpenGLDrawingObject_New( 
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

void _lucMeshViewer_Init( 
                lucMeshViewer*                                            self,
		Mesh*                                                     mesh,
		Name                                                      colourName,
		lucColour                                                 localColour,
		lucColour		                                  shadowColour, 
		lucColour		                                  vacantColour)
{
	self->mesh  = mesh;
	lucColour_FromString( &self->colour, colourName );
	memcpy( &(self->localColour), &localColour, sizeof(lucColour) );
	memcpy( &(self->shadowColour), &shadowColour, sizeof(lucColour) );
	memcpy( &(self->vacantColour), &vacantColour, sizeof(lucColour) );
	assert( Stg_Class_IsInstance( mesh, Mesh_Type ) );

	self->renderEdges = NULL;
}

void _lucMeshViewer_Delete( void* drawingObject ) {
	lucMeshViewer*  self = (lucMeshViewer*)drawingObject;

	_lucOpenGLDrawingObject_Delete( self );
}

void _lucMeshViewer_Print( void* drawingObject, Stream* stream ) {
	lucMeshViewer*  self = (lucMeshViewer*)drawingObject;

	_lucOpenGLDrawingObject_Print( self, stream );
}

void* _lucMeshViewer_Copy( void* drawingObject, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap) {
	lucMeshViewer*  self = (lucMeshViewer*)drawingObject;
	lucMeshViewer* newDrawingObject;
	newDrawingObject = _lucOpenGLDrawingObject_Copy( self, dest, deep, nameExt, ptrMap );
	memcpy( &(newDrawingObject->colour),       &(self->colour),       sizeof(lucColour) );
	memcpy( &(newDrawingObject->localColour),       &(self->localColour),       sizeof(lucColour) );
	memcpy( &(newDrawingObject->shadowColour),       &(self->shadowColour),       sizeof(lucColour) );
	memcpy( &(newDrawingObject->vacantColour),       &(self->vacantColour),       sizeof(lucColour) );

	/* TODO */
	abort();

	return (void*) newDrawingObject;
}


void* _lucMeshViewer_DefaultNew( Name name ) {
	return (void*) _lucMeshViewer_New(
		sizeof(lucMeshViewer),
		lucMeshViewer_Type,
		_lucMeshViewer_Delete,
		_lucMeshViewer_Print,
		NULL,
		_lucMeshViewer_DefaultNew,
		_lucMeshViewer_Construct,
		_lucMeshViewer_Build,
		_lucMeshViewer_Initialise,
		_lucMeshViewer_Execute,
		_lucMeshViewer_Destroy,
		_lucMeshViewer_Setup,
		_lucMeshViewer_Draw,
		_lucMeshViewer_CleanUp,
		_lucMeshViewer_BuildDisplayList,
		name );
}

void _lucMeshViewer_Construct( void* drawingObject, Stg_ComponentFactory* cf, void* data ){
	lucMeshViewer*         self = (lucMeshViewer*)drawingObject;
	Mesh*                  mesh;
	Name localColourName;
	Name shadowColourName;
	Name vacantColourName;
	
	/* Construct Parent */
	_lucOpenGLDrawingObject_Construct( self, cf, data );
	
	mesh = Stg_ComponentFactory_ConstructByKey( cf, self->name, "Mesh", Mesh, True, data );
	localColourName = Stg_ComponentFactory_GetString( cf, self->name, "localColour", "Black") ;
	shadowColourName = Stg_ComponentFactory_GetString( cf, self->name, "shadowColour", "Blue");
	vacantColourName = Stg_ComponentFactory_GetString( cf, self->name, "vacantColour", "Grey");
	self->nodeNumbers = Stg_ComponentFactory_GetBool( cf, self->name, "nodeNumbers", False);
	self->elementNumbers = Stg_ComponentFactory_GetBool( cf, self->name, "elementNumbers", False);
	self->displayNodes = Stg_ComponentFactory_GetBool( cf, self->name, "displayNodes", False);

	lucColour_FromString( &self->localColour, localColourName );
	lucColour_FromString( &self->shadowColour, shadowColourName );
	lucColour_FromString( &self->vacantColour, vacantColourName );
   
	_lucMeshViewer_Init( 
			self, 
		        mesh,
			Stg_ComponentFactory_GetString( cf, self->name, "colour", "black" ),
			self->localColour,
			self->shadowColour,
	                self->vacantColour);
}

void _lucMeshViewer_Build( void* drawingObject, void* data ) {
}

void _lucMeshViewer_Initialise( void* drawingObject, void* data ) {
	lucMeshViewer*	self = (lucMeshViewer*)drawingObject;

	assert( self );

	if( Mesh_HasIncidence( self->mesh, MT_EDGE, MT_VERTEX ) )
		self->renderEdges = lucMeshViewer_RenderEdges_WithInc;
	else {
		lucMeshViewer_BuildEdges( self );
		self->renderEdges = lucMeshViewer_RenderEdges;
	}
}

void _lucMeshViewer_Execute( void* drawingObject, void* data ) {
}

void _lucMeshViewer_Destroy( void* drawingObject, void* data ) {
}

void _lucMeshViewer_Setup( void* drawingObject, void* _context ) {
	lucMeshViewer*          self                = (lucMeshViewer*)drawingObject;
	
	_lucOpenGLDrawingObject_Setup( self, _context );
	 lucMeshViewer_UpdateVariables( self );
	
}

void lucMeshViewer_UpdateVariables( void* drawingObject ) {
}

void _lucMeshViewer_Draw( void* drawingObject, lucWindow* window, lucViewportInfo* viewportInfo, void* _context ) {
	lucMeshViewer*          self          = (lucMeshViewer*)drawingObject;
	lucCamera*               camera        = viewportInfo->viewport->camera;
	XYZ                      normal;
	
	StGermain_VectorSubtraction( normal, camera->coord, camera->focalPoint, 3 );
	glNormal3dv(normal);
	_lucOpenGLDrawingObject_Draw( self, window, viewportInfo, _context );
}

void _lucMeshViewer_CleanUp( void* drawingObject, void* context ) {
	lucMeshViewer*          self          = (lucMeshViewer*)drawingObject;
	
	_lucOpenGLDrawingObject_CleanUp( self, context );
}

void _lucMeshViewer_BuildDisplayList( void* drawingObject, void* _context ) {
	lucMeshViewer*	self = (lucMeshViewer*)drawingObject;
	lucColour	colour;

	/* Initialise colour value */
	memcpy( &colour, &self->colour, sizeof(lucColour) );
	lucColour_SetOpenGLColour( &colour );

	glPointSize( 1.0 );

	/* Plot the mesh */
	lucMeshViewer_Render( drawingObject );
}

#if 0
void lucMeshViewer_RenderGlobal( void* drawingObject ) {
	lucMeshViewer*		self = (lucMeshViewer*)drawingObject;
	unsigned	edge_I;
	unsigned	vert_I;
	
	glColor3f( self->localColour.red, self->localColour.green, self->localColour.blue );
	
	/* Render vertices */
	glBegin( GL_POINTS );
	for( vert_I = 0; vert_I < self->vertCnt * 3; vert_I += 3 ) {
		glVertex3dv( &self->verts[vert_I] );
	}
	glEnd();
	
	/* Render edges */
	glBegin( GL_LINES );
	for( edge_I = 0; edge_I < self->edgeCnt * 2; edge_I += 2 ) {
		unsigned	vert_I = self->edges[edge_I] * 3;
		unsigned	vert_J = self->edges[edge_I + 1] * 3;
		
		glVertex3dv( &self->verts[vert_I] );
		glVertex3dv( &self->verts[vert_J] );
	}
	glEnd();
}

void lucMeshViewer_PrintAllElementsNumber( void* drawingObject, Partition_Index rank ) {
	lucMeshViewer*	     self = (lucMeshViewer*)drawingObject;
	Coord                avgCoord;
	Coord                offset;
	char                 elementNumString[100];
	Dimension_Index      dim_I;
	Node_LocalIndex      node_lI;
	Node_Index           elNode_I;
	Element_LocalIndex   element_lI;


	glDisable(GL_LIGHTING); /*if the lighting is not disabled, the colour won't appear for the numbers*/
	glColor3f( self->localColour.red, self->localColour.green, self->localColour.blue );

	/* Prints the element numbers */
	offset[0] = -0.01;
	offset[1] = -0.01;
	offset[2] = 0;
	for ( element_lI = 0; element_lI < self->mesh->elementLocalCount; element_lI++ ) 
	{	
		sprintf( elementNumString, "el%u", element_lI );

		for ( dim_I=0; dim_I < 3; dim_I++) {	
			avgCoord[dim_I] = 0;
		}
		for ( elNode_I=0; elNode_I < self->mesh->elementNodeCountTbl[element_lI]; elNode_I++ ) {
			node_lI = self->mesh->elementNodeTbl[element_lI][elNode_I];
			for ( dim_I=0; dim_I < ((HexaEL*)(self->mesh->layout->elementLayout))->dim; dim_I++) {	
				avgCoord[dim_I] += self->mesh->nodeCoord[node_lI][dim_I];
			}	
		}
		for ( dim_I=0; dim_I < ((HexaEL*)(self->mesh->layout->elementLayout))->dim; dim_I++) {	
			avgCoord[dim_I] /= (double)self->mesh->elementNodeCountTbl[element_lI];
		}

		if ( ((HexaEL*)(self->mesh->layout->elementLayout))->dim == 2) {
			glRasterPos2f( (float)avgCoord[0] + offset[0], (float)avgCoord[1] + offset[1] );		
		}	
		else {  
			glRasterPos3f( (float)avgCoord[0] + offset[0], (float)avgCoord[1] + offset[1],
				(float)avgCoord[2] + offset[2] );
		}	

		lucPrintString( elementNumString );
	}
	glEnable(GL_LIGHTING);


}

void lucMeshViewer_PrintAllNodesNumber( void* drawingObject, Partition_Index rank ) {
	lucMeshViewer*		self = (lucMeshViewer*)drawingObject;
	double*              coord;
	char                 nodeNumString[100];
	Node_LocalIndex      node_lI;

	glDisable(GL_LIGHTING); /*if the lighting is not disabled, the colour won't appear for the numbers*/
	glColor3f( self->localColour.red, self->localColour.green, self->localColour.blue );

	/* Prints the node numbers */
	for ( node_lI = 0; node_lI < self->mesh->nodeLocalCount; node_lI++ ) 
	{	
		sprintf( nodeNumString, "nl%u", node_lI );
		coord = self->mesh->nodeCoord[node_lI];
		if ( ((HexaEL*)(self->mesh->layout->elementLayout))->dim == 2)
			glRasterPos2f( (float)coord[0] + 0.015, (float)coord[1] + 0.015 );		
		else   
			glRasterPos3f( (float)coord[0] + 0.015, (float)coord[1] + 0.015, (float)coord[2] + 0.015 );

		lucPrintString( nodeNumString );
	}
	glEnable(GL_LIGHTING);

}

void lucMeshViewer_PrintNodeNumber( void* drawingObject, Coord coord, int* nodeNumber ) {
	lucMeshViewer*	     self = (lucMeshViewer*)drawingObject;
	char                 nodeNumString[100];

	unsigned dim  = ((HexaEL*)(self->mesh->layout->elementLayout))->dim ;

	lucMeshViewer_ClosestNode(self, coord, nodeNumber);
	
	glDisable(GL_LIGHTING); /*if the lighting is not disabled, the colour won't appear for the numbers*/
	glColor3f( self->localColour.red, self->localColour.green, self->localColour.blue );

	/* Prints the node numbers */
	sprintf( nodeNumString, "nl%u", *nodeNumber );
	if (dim == 2)
		glRasterPos2f( (float)coord[0] + 0.015, (float)coord[1] + 0.015 );		
	else   
		glRasterPos3f( (float)coord[0] + 0.015, (float)coord[1] + 0.015, (float)coord[2] + 0.015 );

	lucPrintString( nodeNumString );
}

void lucMeshViewer_PrintElementNumber( void* drawingObject, Coord coord, int* elementNumber ) {
	lucMeshViewer*	     self = (lucMeshViewer*)drawingObject;
	Coord                avgCoord;
	Coord                offset;
	char                 elementNumString[100];
	Dimension_Index      dim_I;
	Node_LocalIndex      node_lI;
	Node_Index           elNode_I;
	Element_LocalIndex   element_lI;

        lucMeshViewer_FindElementNumber(drawingObject, coord, elementNumber);

	glDisable(GL_LIGHTING); /*if the lighting is not disabled, the colour won't appear for the numbers*/
	glColor3f( self->localColour.red, self->localColour.green, self->localColour.blue );

	/* Prints the element numbers */
	offset[0] = -0.01;
	offset[1] = -0.01;
	offset[2] = 0;
	element_lI = *elementNumber;
	
	sprintf( elementNumString, "el%u", element_lI );

	for ( dim_I=0; dim_I < 3; dim_I++) {	
		avgCoord[dim_I] = 0;
	}
	for ( elNode_I=0; elNode_I < self->mesh->elementNodeCountTbl[element_lI]; elNode_I++ ) {
		node_lI = self->mesh->elementNodeTbl[element_lI][elNode_I];
		for ( dim_I=0; dim_I < ((HexaEL*)(self->mesh->layout->elementLayout))->dim; dim_I++) {	
			avgCoord[dim_I] += self->mesh->nodeCoord[node_lI][dim_I];
		}	
	}
	for ( dim_I=0; dim_I < ((HexaEL*)(self->mesh->layout->elementLayout))->dim; dim_I++) {	
		avgCoord[dim_I] /= (double)self->mesh->elementNodeCountTbl[element_lI];
	}

	if ( ((HexaEL*)(self->mesh->layout->elementLayout))->dim == 2) {
		glRasterPos2f( (float)avgCoord[0] + offset[0], (float)avgCoord[1] + offset[1] );		
	}	
	else {  
		glRasterPos3f( (float)avgCoord[0] + offset[0], (float)avgCoord[1] + offset[1],
			(float)avgCoord[2] + offset[2] );
	}	

	lucPrintString( elementNumString );
	glEnable(GL_LIGHTING);
}
#endif

void lucMeshViewer_RenderLocal( void* drawingObject ) {
	lucMeshViewer*	self = (lucMeshViewer*)drawingObject;
	Mesh*		mesh;
	vertexFuncType*	vertexFunc;

	assert( self );
	assert( Mesh_GetDomainSize( self->mesh, MT_VERTEX ) );
	assert( self->renderEdges );

	/* Shortcuts. */
	mesh = self->mesh;

	/* Pick the correct dimension. */
	if( Mesh_GetDimSize( mesh ) == 3 )
		vertexFunc = glVertex3dv;
	else
		vertexFunc = glVertex2dv;

	/* Set color. */
	glColor3f( self->localColour.red, self->localColour.green, self->localColour.blue );

	/* Render vertices. */
	if(self->displayNodes){
		unsigned	nVerts;
		unsigned	v_i;

		nVerts = Mesh_GetLocalSize( mesh, MT_VERTEX );
		glPointSize( 5 );
		glBegin( GL_POINTS );
		for( v_i = 0; v_i < nVerts; v_i ++ )
			vertexFunc( Mesh_GetVertex( mesh, v_i ) );
		glEnd();
	}

	/* Render edges */
	self->renderEdges( self, vertexFunc );

#if 0
	/* Prints the element numbers */
	if( self->elementNumbers )
		lucMeshViewer_PrintAllElementsNumber( self );

	/* Prints the node numbers */
	if( self->nodeNumbers )
		lucMeshViewer_PrintAllNodesNumber( self );
#endif
}

#if 0
void lucMeshViewer_RenderShadow( void* drawingObject, Partition_Index rank ) {
	lucMeshViewer*		self = (lucMeshViewer*)drawingObject;
	unsigned	sEdge_I;
	
	assert( rank < self->rankCnt );
	
	if( !self->shadowEdgeCnts[rank] || !self->shadowEdges[rank] ) {
		return;
	}
	
	glColor3f( self->shadowColour.red, self->shadowColour.green, self->shadowColour.blue );
	

	/* Render edges */
	glBegin( GL_LINES );
	for( sEdge_I = 0; sEdge_I < self->shadowEdgeCnts[rank]; sEdge_I++ ) {
		unsigned	  edge_I = self->shadowEdges[rank][sEdge_I] * 2;
		Node_DomainIndex  node1_dI = Mesh_NodeMapGlobalToDomain( self->mesh, self->edges[edge_I] );
		Node_DomainIndex  node2_dI = Mesh_NodeMapGlobalToDomain( self->mesh, self->edges[edge_I + 1] );
		double*         coord1;
		double*         coord2;

		coord1 = self->mesh->nodeCoord[node1_dI];
		coord2 = self->mesh->nodeCoord[node2_dI];
		
		glVertex3dv( coord1 );
		glVertex3dv( coord2 );
	}
	glEnd();
}

void lucMeshViewer_RenderVacant( void* drawingObject, Partition_Index rank ) {
	lucMeshViewer*		self = (lucMeshViewer*)drawingObject;
	unsigned	     lEdge_I;
	unsigned	     edge_I;
	Node_LocalIndex      node1_lI;
	Node_LocalIndex      node2_lI;
	double*              coord1;
	double*              coord2;
	
	assert( rank < self->rankCnt );
	
	assert( rank < self->rankCnt );
	
	if( !self->vacantEdgeCnts || !self->vacantEdgeCnts[rank] || !self->vacantEdges || !self->vacantEdges[rank] ) {
		return;
	}
	
	glColor3f( self->vacantColour.red, self->vacantColour.green, self->vacantColour.blue );
	
	glBegin( GL_LINES );
	for( lEdge_I = 0; lEdge_I < self->localEdgeCnts[rank]; lEdge_I++ ) {
		edge_I = self->localEdges[rank][lEdge_I] * 2;

		node1_lI = Mesh_NodeMapGlobalToLocal( self->mesh, self->edges[edge_I] );
		node2_lI = Mesh_NodeMapGlobalToLocal( self->mesh, self->edges[edge_I + 1] );

		coord1 = self->mesh->nodeCoord[node1_lI];	
		coord2 = self->mesh->nodeCoord[node2_lI];	
		
		glVertex3dv( coord1 );
		glVertex3dv( coord2 );
	}
	glEnd();
}
#endif

void lucMeshViewer_Render( void* drawingObject ) {
	lucMeshViewer*		self = (lucMeshViewer*)drawingObject;

	lucMeshViewer_RenderLocal( self );
}


/*--------------------------------------------------------------------------------------------------------------------------
** Private Member functions
*/

void lucMeshViewer_BuildEdges( lucMeshViewer* self ) {
#if 0
	assert( self );

	nVerts = Mesh_GetLocalSize( );
	done = AllocArray( Bool, nVerts );
#endif
}

void lucMeshViewer_RenderEdges_WithInc( lucMeshViewer* self, vertexFuncType* vertexFunc ) {
	unsigned	nEdges;
	unsigned	nIncVerts, *incVerts;
	IArray*		inc;
	unsigned	e_i;

	assert( self );
	assert( Mesh_GetDomainSize( self->mesh, MT_EDGE ) && 
		Mesh_HasIncidence( self->mesh, MT_EDGE, MT_VERTEX ) );

	nEdges = Mesh_GetLocalSize( self->mesh, MT_EDGE );
	glDisable(GL_LIGHTING);
	glBegin( GL_LINES );
	inc = IArray_New();
	for( e_i = 0; e_i < nEdges; e_i++ ) {
		Mesh_GetIncidence( self->mesh, MT_EDGE, e_i, MT_VERTEX, inc );
		nIncVerts = IArray_GetSize( inc );
		incVerts = IArray_GetPtr( inc );
		assert( nIncVerts == 2 );

		vertexFunc( Mesh_GetVertex( self->mesh, incVerts[0] ) );
		vertexFunc( Mesh_GetVertex( self->mesh, incVerts[1] ) );
	}
	glEnd();
	glEnable(GL_LIGHTING);

	NewClass_Delete( inc );
}

void lucMeshViewer_RenderEdges( lucMeshViewer* self, vertexFuncType* vertexFunc ) {
	unsigned	nEdges, **edges;
	unsigned	e_i;

	assert( self );
	assert( self->nEdges && self->edges );

	nEdges = self->nEdges;
	edges = self->edges;
	glDisable(GL_LIGHTING);
	glBegin( GL_LINES );
	for( e_i = 0; e_i < nEdges; e_i++ ) {
		vertexFunc( Mesh_GetVertex( self->mesh, edges[e_i][0] ) );
		vertexFunc( Mesh_GetVertex( self->mesh, edges[e_i][1] ) );
	}
	glEnd();
	glEnable(GL_LIGHTING);
}
