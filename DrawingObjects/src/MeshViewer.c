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
#include <gl.h>
#include <glu.h>
#include <string.h>

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
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

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
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
}


void _lucMeshViewer_Execute( void* drawingObject, void* data ) {}
void _lucMeshViewer_Destroy( void* drawingObject, void* data ) {}

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
	lucMeshViewer*          self                = (lucMeshViewer*)drawingObject;
	lucColour              colour;

	/*Geometry*		geometry;*/
	Node_GlobalIndex	point_I;
	Edge_Index          edge_I;
	ElementLayout*      elementLayout;
	Partition_Index		rank_I;

	/* Stuff to construct the layout */
	MeshDecomp*		decomp;
	MeshLayout*		meshLayout;
	Dimension_Index		numPartitionedDims;
	Partition_Index		maxRank = 0;
  	Node_Index nodeCount;

	/*TODO*/     
	numPartitionedDims = 2;
	
	meshLayout =  self->mesh->layout;
	decomp = meshLayout->decomp;

	/* Make sure "StoreAll" is set so proc 0 can get all info */
	maxRank = decomp->procsInUse;
	
/*TODO	if( rank == 0 ) {
		glMesh = GLMesh_New();
		GLMesh_BuildFromMesh( glMesh, meshLayout );
	}
*/	

   /* Ensure everything is already freed*/ 
	_lucMeshViewer_CleanMem( self, NULL );
	
	/* Copy the vertices */
	self->vertCnt = meshLayout->nodeLayout->nodeCount;
	self->verts = Memory_Alloc_Array( GLdouble, self->vertCnt * 3, "lucMeshViewer->verts" );

	assert( self->verts );

	nodeCount = self->mesh->nodeLocalCount;
	
	for( point_I = 0; point_I < nodeCount; point_I++ ) {
		double* nodeCoord = Mesh_CoordAt( self->mesh, point_I );
		unsigned	vert_I = point_I * 3;
			
		self->verts[vert_I] = (GLdouble)nodeCoord[0];
		self->verts[vert_I + 1] = (GLdouble)nodeCoord[1];
		self->verts[vert_I + 2] = (GLdouble)nodeCoord[2];

	}
	
	/* Build the edges */
	elementLayout = meshLayout->elementLayout;
	self->edgeCnt = elementLayout->edgeCount;
	self->edges = Memory_Alloc_Array( unsigned, self->edgeCnt * 2, "lucMeshViewer->edges" );
	assert( self->edges );

	for( edge_I = 0; edge_I < elementLayout->edgeCount; edge_I++ ) {
		unsigned	glEdge_I = edge_I * 2;
		Edge		edge;
		
		elementLayout->edgeAt( elementLayout, edge_I, edge );
		self->edges[glEdge_I] = (unsigned)edge[0];
		self->edges[glEdge_I + 1] = (unsigned)edge[1];
	}
	
	/* Build local edge indices */
	self->rankCnt = meshLayout->decomp->procsInUse;
	self->localEdgeCnts = Memory_Alloc_Array( unsigned, self->rankCnt, "lucMeshViewer->localEdgeCnts" );
	memset( self->localEdgeCnts, 0, sizeof(unsigned) * self->rankCnt );
	self->localEdges = Memory_Alloc_Array( unsigned*, self->rankCnt, "lucMeshViewer->localEdges" );
	memset( self->localEdges, 0, sizeof(unsigned*) * self->rankCnt );
	self->shadowEdgeCnts = Memory_Alloc_Array( unsigned, self->rankCnt, "lucMeshViewer->localEdges" );
	memset( self->shadowEdgeCnts, 0, sizeof(unsigned) * self->rankCnt );
	self->shadowEdges = Memory_Alloc_Array( unsigned*, self->rankCnt, "lucMeshViewer->shadowEdges" );
	memset( self->shadowEdges, 0, sizeof(unsigned*) * self->rankCnt );
	self->vacantEdgeCnts = Memory_Alloc_Array( unsigned, self->rankCnt, "lucMeshViewer->vacantEdgeCnts" );
	memset( self->vacantEdgeCnts, 0, sizeof(unsigned) * self->rankCnt );
	self->vacantEdges = Memory_Alloc_Array( unsigned*, self->rankCnt, "lucMeshViewer->vacantEdges" );
	memset( self->vacantEdges, 0, sizeof(unsigned*) * self->rankCnt );
	
	for( rank_I = 0; rank_I < self->rankCnt; rank_I++ ) {
		_lucMeshViewer_BuildLocalEdges( self, meshLayout, rank_I );
		_lucMeshViewer_BuildShadowEdges( self, meshLayout, rank_I );
		_lucMeshViewer_BuildVacantEdges( self, meshLayout, rank_I );
	}

	
	/* Initialise colour value */
	memcpy( &colour, &self->colour, sizeof(lucColour) );
	lucColour_SetOpenGLColour( &colour );

	glPointSize( 1.0 );
	
	/* Plot the mesh */
	lucMeshViewer_RenderRank( drawingObject, 0 );
}

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


void lucMeshViewer_ClosestNode( void* self, Coord crd, int* nodeNumber ) {
	Bool		done;
        Mesh*		mesh = ((lucMeshViewer*)self)->mesh; 
	Coord*		nodeCrds = mesh->nodeCoord;
	unsigned	curNode;
	unsigned        nDims;
	
	nDims  = ((HexaEL*)(mesh->layout->elementLayout))->dim ;

	/* Begin somewhere in the middle. */
	curNode = mesh->nodeLocalCount / 2;

	if(!mesh->nodeNeighbourCountTbl){
		Mesh_ActivateNodeNeighbourTbl( mesh ); 
	}

	/* Loop until we've found closest local node. */
	do {
		unsigned	nNbrs = mesh->nodeNeighbourCountTbl[curNode];
		unsigned*	nbrs = mesh->nodeNeighbourTbl[curNode];
		double		dist;
		double		tmp;
		unsigned	nbr_i, d_i;

		/* Assume we'll be done after this loop. */
		done = True;

		/* Calc distance squared to current node. */
		tmp = nodeCrds[curNode][0] - crd[0];
		dist = tmp * tmp;
		for( d_i = 1; d_i < nDims; d_i++ ) {
			tmp = nodeCrds[curNode][d_i] - crd[d_i];
			dist += tmp * tmp;
		}

		/* Compare to neighbours. */
		for( nbr_i = 0; nbr_i < nNbrs; nbr_i++ ) {
			double	nbrDist;

			/* Just in case... */
			if( nbrs[nbr_i] >= mesh->nodeLocalCount )
				continue;

			tmp = nodeCrds[nbrs[nbr_i]][0] - crd[0];
			nbrDist = tmp * tmp;
			for( d_i = 1; d_i < nDims; d_i++ ) {
				tmp = nodeCrds[nbrs[nbr_i]][d_i] - crd[d_i];
				nbrDist += tmp * tmp;
			}

			if( nbrDist < dist ) {
				curNode = nbrs[nbr_i];
				dist = nbrDist;
				done = False;
			}
		}
	}
	while( !done );

	*nodeNumber = curNode;
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

void lucMeshViewer_FindElementNumber(void* drawingObject, Coord coord, int* elementNumber){
	Mesh*		mesh = ((lucMeshViewer*)drawingObject)->mesh; 
	MeshLayout*		mLayout = mesh->layout;
	ElementLayout*		eLayout = mLayout->elementLayout;
 	Element_DomainIndex	elementCoordIn = (unsigned)-1;
	
	/* locate which mesh element given coord is in : use inclusive upper boundaries to save
		the need to use shadow space if possible */
	if( eLayout->type == ParallelPipedHexaEL_Type ) {
		elementCoordIn = eLayout->elementWithPoint( eLayout, mLayout->decomp, coord, mesh, 
							    INCLUSIVE_UPPER_BOUNDARY, 0, NULL );
	}
	else {
		unsigned	cNode;
 
		/* Find closest node to point. */
		lucMeshViewer_ClosestNode( drawingObject, coord, (int*)&cNode );

		/* Find with hint of incident elements. */
		elementCoordIn = eLayout->elementWithPoint( eLayout, mLayout->decomp, coord, mesh, 
							    INCLUSIVE_UPPER_BOUNDARY, 
							    mesh->nodeElementCountTbl[cNode], mesh->nodeElementTbl[cNode] );

		/* If still no cigar, brute force. */
		if ( elementCoordIn >= mesh->elementDomainCount ) {
			elementCoordIn = eLayout->elementWithPoint( eLayout, mLayout->decomp, coord, mesh, 
								    INCLUSIVE_UPPER_BOUNDARY, 0, NULL );
		}
	}
      	*elementNumber = 	elementCoordIn;
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



void lucMeshViewer_RenderLocal( void* drawingObject, Partition_Index rank ) {
	lucMeshViewer*		self = (lucMeshViewer*)drawingObject;
	unsigned	lEdge_I;
	Node_LocalIndex      node_lI;
	unsigned	     edge_I;
	Node_LocalIndex      node1_lI;
	Node_LocalIndex      node2_lI;
	double*              coord1;
	double*              coord2;
	
	assert( rank < self->rankCnt );
	
	if( !self->localEdgeCnts[rank] || !self->localEdges[rank] ) {
		return;
	}
	
	glColor3f( self->localColour.red, self->localColour.green, self->localColour.blue );
	
	/* Render nodes */
	if(self->displayNodes){
		glPointSize( 5 );
		glBegin( GL_POINTS );
		for( node_lI = 0; node_lI < self->mesh->nodeLocalCount; node_lI ++ ) {
			glVertex3dv( self->mesh->nodeCoord[node_lI] );
		}
		glEnd();
	}
	
	/* Render edges */
	glDisable(GL_LIGHTING);
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
	glEnable(GL_LIGHTING);

	/* Prints the element numbers */
	if(self->elementNumbers)
		lucMeshViewer_PrintAllElementsNumber(self, rank);

	/* Prints the node numbers */
	if(self->nodeNumbers)
		lucMeshViewer_PrintAllNodesNumber(self, rank);
	
}


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


void lucMeshViewer_RenderRank( void* drawingObject, Partition_Index rank ) {
	lucMeshViewer*		self = (lucMeshViewer*)drawingObject;
	Processor_Index         myRank = self->mesh->layout->decomp->rank;

	if ( rank == myRank ) {
		lucMeshViewer_RenderLocal( self, myRank );
		lucMeshViewer_RenderShadow( self, myRank );
        }
	else{
		lucMeshViewer_RenderVacant( self, myRank );
	}
}



/*--------------------------------------------------------------------------------------------------------------------------
** Private Member functions
*/

void _lucMeshViewer_BuildLocalEdges( void* meshViewer, MeshLayout* mesh, Partition_Index rank ) {
	lucMeshViewer*		self = (lucMeshViewer*)meshViewer;
	Index		localElementCnt;
	Index*		localElementSet;
	int		lEl_i;
	
	/* Old code:
	** IndexSet_GetMembers( mesh->decomp->localElementSets[rank], &localElementCnt, &localElementSet );
	**
	** This should really use a StGermain mesh, not the mesh layout. */
	localElementCnt = mesh->decomp->elementLocalCount;
	localElementSet = Memory_Alloc_Array( unsigned, localElementCnt, "localElementSet" );
	for( lEl_i = 0; lEl_i < localElementCnt; lEl_i++ ) {
		localElementSet[lEl_i] = mesh->decomp->elementMapLocalToGlobal( mesh->decomp, lEl_i );
	}
	
	
	if( localElementCnt ) {
		if( self->localEdges[rank] ) {
			Memory_Free( self->localEdges );
		}
		
		self->localEdgeCnts[rank] = ElementLayout_BuildEdgeSubset( mesh->elementLayout, 
									   localElementCnt, 
									   localElementSet, 
									   &self->localEdges[rank] );
	}
	
	Memory_Free( localElementSet );
}


void _lucMeshViewer_BuildShadowEdges( void* meshViewer, MeshLayout* mesh, Partition_Index rank ) {
	lucMeshViewer*		self = (lucMeshViewer*)meshViewer;
	Index		shadowElementCnt;
	Index*		shadowElementSet;
	
	if( !mesh->decomp->shadowElementSets || !mesh->decomp->shadowElementSets[rank] ) {
		return;
	}
	
	IndexSet_GetMembers( mesh->decomp->shadowElementSets[rank], &shadowElementCnt, &shadowElementSet );
	
	if( shadowElementCnt ) {
		IndexSet*	localEdgeSet;
		IndexSet*	shadowEdgeSet;
		Index		lEdge_I;
		Index		sEdge_I;
		
		localEdgeSet = IndexSet_New( mesh->decomp->elementLayout->edgeCount );
		for( lEdge_I = 0; lEdge_I < self->localEdgeCnts[rank]; lEdge_I++ ) {
			IndexSet_Add( localEdgeSet, self->localEdges[rank][lEdge_I] );
		}
		
		if( self->shadowEdges[rank] ) {
			Memory_Free( self->shadowEdges );
		}
		self->shadowEdgeCnts[rank] = ElementLayout_BuildEdgeSubset( mesh->elementLayout, 
									   shadowElementCnt, 
									   shadowElementSet, 
									   &self->shadowEdges[rank] );
		shadowEdgeSet = IndexSet_New( mesh->decomp->elementLayout->edgeCount );
		for( sEdge_I = 0; sEdge_I < self->shadowEdgeCnts[rank]; sEdge_I++ ) {
			if( !IndexSet_IsMember( localEdgeSet, self->shadowEdges[rank][sEdge_I] ) ) {
				IndexSet_Add( shadowEdgeSet, self->shadowEdges[rank][sEdge_I] );
			}
		}
		
		Memory_Free( self->shadowEdges[rank] );
		IndexSet_GetMembers( shadowEdgeSet, &self->shadowEdgeCnts[rank], &self->shadowEdges[rank] );
		Stg_Class_Delete( shadowEdgeSet );
	}
	
	Memory_Free( shadowElementSet );
}


void _lucMeshViewer_BuildVacantEdges( void* meshViewer, MeshLayout* mesh, Partition_Index rank ) {
	lucMeshViewer*		self = (lucMeshViewer*)meshViewer;
	ElementLayout*	elementLayout = mesh->decomp->elementLayout;
	IndexSet*	domainEdgeSet;
	IndexSet*	vacantEdgeSet;
	Index		gEdge_I;
	
	domainEdgeSet = IndexSet_New( elementLayout->edgeCount );
	
	if( self->localEdgeCnts && self->localEdgeCnts[rank] ) {
		unsigned	lEdge_I;
		
		for( lEdge_I = 0; lEdge_I < self->localEdgeCnts[rank]; lEdge_I++ ) {
			IndexSet_Add( domainEdgeSet, self->localEdges[rank][lEdge_I] );
		}
	}
	
	if( self->shadowEdgeCnts && self->shadowEdgeCnts[rank] ) {
		unsigned	sEdge_I;
		
		for( sEdge_I = 0; sEdge_I < self->shadowEdgeCnts[rank]; sEdge_I++ ) {
			IndexSet_Add( domainEdgeSet, self->shadowEdges[rank][sEdge_I] );
		}
	}
	
	vacantEdgeSet = IndexSet_New( elementLayout->edgeCount );
	
	for( gEdge_I = 0; gEdge_I < elementLayout->edgeCount; gEdge_I++ ) {
		if( !IndexSet_IsMember( domainEdgeSet, gEdge_I ) ) {
			IndexSet_Add( vacantEdgeSet, gEdge_I );
		}
	}
	
	IndexSet_GetMembers( vacantEdgeSet, &self->vacantEdgeCnts[rank], &self->vacantEdges[rank] );
	
	Stg_Class_Delete( domainEdgeSet );
	Stg_Class_Delete( vacantEdgeSet );
}

void _lucMeshViewer_CleanMem( void* drawingObject, void* data ) {
	lucMeshViewer*		self = (lucMeshViewer*)drawingObject;
	
	if( self->verts ) {
		Memory_Free( self->verts );
		self->verts = NULL;
	}
	
	if( self->edges ) {
		Memory_Free( self->edges );
		self->edges = NULL;
	}
	
	if( self->localEdgeCnts ) {
		Memory_Free( self->localEdgeCnts );
		self->localEdgeCnts = NULL;
	}
	
	if( self->localEdges ) {
		Partition_Index		rank_I;
		
		for( rank_I = 0; rank_I < self->rankCnt; rank_I++ ) {
			if( self->localEdges[rank_I] ) {
				Memory_Free( self->localEdges[rank_I] );
			}
		}
		Memory_Free( self->localEdges );
		self->localEdges = NULL;
	}
	
	if( self->shadowEdgeCnts ) {
		Memory_Free( self->shadowEdgeCnts );
		self->shadowEdgeCnts = NULL;
	}
	
	if( self->shadowEdges ) {
		Partition_Index		rank_I;
		
		for( rank_I = 0; rank_I < self->rankCnt; rank_I++ ) {
			if( self->shadowEdges[rank_I] ) {
				Memory_Free( self->shadowEdges[rank_I] );
			}
		}
		Memory_Free( self->shadowEdges );
		self->shadowEdges = NULL;
	}
	
	if( self->vacantEdgeCnts ) {
		Memory_Free( self->vacantEdgeCnts );
		self->vacantEdgeCnts = NULL;
	}
	
	if( self->vacantEdges ) {
		Partition_Index		rank_I;
		
		for( rank_I = 0; rank_I < self->rankCnt; rank_I++ ) {
			if( self->vacantEdges[rank_I] ) {
				Memory_Free( self->vacantEdges[rank_I] );
			}
		}
		Memory_Free( self->vacantEdges );
		self->vacantEdges = NULL;
	}
}
