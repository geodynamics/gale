/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003-2006, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street,
**	Melbourne, 3053, Australia.
** Copyright (c) 2005-2006, Monash Cluster Computing, Building 28, Monash University Clayton Campus,
**	Victoria, 3800, Australia
**
** Primary Contributing Organisations:
**	Victorian Partnership for Advanced Computing Ltd, Computational Software Development - http://csd.vpac.org
**	Australian Computational Earth Systems Simulator - http://www.access.edu.au
**	Monash Cluster Computing - http://www.mcc.monash.edu.au
**
** Contributors:
**	Robert Turnbull, Research Assistant, Monash University. (robert.turnbull@sci.monash.edu.au)
**	Patrick D. Sunter, Software Engineer, VPAC. (patrick@vpac.org)
**	Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
**	Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**	David May, PhD Student, Monash University (david.may@sci.monash.edu.au)
**	Vincent Lemiale, Postdoctoral Fellow, Monash University. (vincent.lemiale@sci.monash.edu.au)
**	Julian Giordani, Research Assistant, Monash University. (julian.giordani@sci.monash.edu.au)
**	Louis Moresi, Associate Professor, Monash University. (louis.moresi@sci.monash.edu.au)
**	Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**	Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
**	David Stegman, Postdoctoral Fellow, Monash University. (david.stegman@sci.monash.edu.au)
**	Wendy Sharples, PhD Student, Monash University (wendy.sharples@sci.monash.edu.au)
**
**  This library is free software; you can redistribute it and/or
**  modify it under the terms of the GNU Lesser General Public
**  License as published by the Free Software Foundation; either
**  version 2.1 of the License, or (at your option) any later version.
**
**  This library is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
**  Lesser General Public License for more details.
**
**  You should have received a copy of the GNU Lesser General Public
**  License along with this library; if not, write to the Free Software
**  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
**
** $Id: AdaptiveVoronoi.c 376 2006-10-18 06:58:41Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgFEM/StgFEM.h>

#include "types.h"
#include "DiscreteVoronoi.h"
#include "AdaptiveVoronoi.h"

#include <assert.h>
#include <string.h>
#include <math.h>

#define UNCLAIMED                    ((unsigned int) -1)
#define FINISHED                     ((unsigned int) -1)

/* Cell Vertex Index Macros */
#define LEFT_BOTTOM_FRONT            0
#define RIGHT_BOTTOM_FRONT           1
#define RIGHT_TOP_FRONT              2
#define LEFT_TOP_FRONT               3

/* Cell Vertex Index Macros - For 3D Only */
#define LEFT_BOTTOM_BACK             4
#define RIGHT_BOTTOM_BACK            5
#define RIGHT_TOP_BACK               6
#define LEFT_TOP_BACK                7

/* New Vertex Index Macros */
#define MIDDLE_MIDDLE_FRONT          0
#define LEFT_MIDDLE_FRONT            1
#define RIGHT_MIDDLE_FRONT           2
#define MIDDLE_BOTTOM_FRONT          3
#define MIDDLE_TOP_FRONT             4

/* New Vertex Index Macros - For 3D Only */
#define LEFT_BOTTOM_MIDDLE           5
#define RIGHT_BOTTOM_MIDDLE          6
#define RIGHT_TOP_MIDDLE             7
#define LEFT_TOP_MIDDLE              8
#define MIDDLE_MIDDLE_MIDDLE         9
#define LEFT_MIDDLE_MIDDLE           10
#define RIGHT_MIDDLE_MIDDLE          11
#define MIDDLE_BOTTOM_MIDDLE         12
#define MIDDLE_TOP_MIDDLE            13

#define MIDDLE_MIDDLE_BACK           14
#define LEFT_MIDDLE_BACK             15
#define RIGHT_MIDDLE_BACK            16
#define MIDDLE_BOTTOM_BACK           17
#define MIDDLE_TOP_BACK              18


/* Textual name of this class */
const Type AdaptiveVoronoi_Type = "AdaptiveVoronoi";

/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/

AdaptiveVoronoi* _AdaptiveVoronoi_New(
		SizeT                                              _sizeOfSelf, 
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
		DiscreteVoronoi_CalculateForCellFunction*          _calculate,		
		DiscreteVoronoi_GetParticleIndexFunction*          _getParticleIndex,
		DiscreteVoronoi_GetVolumeFunction*                 _getVolume,
		DiscreteVoronoi_GetCentroidFunction*               _getCentroid,	
		Name                                               name )
{
	AdaptiveVoronoi* self;
	
	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(AdaptiveVoronoi) );
	self = (AdaptiveVoronoi*)_DiscreteVoronoi_New( 
			_sizeOfSelf,
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
			_calculate,
			_getParticleIndex,
			_getVolume,
			_getCentroid,			
			name );

	
	/* General info */

	/* Virtual Info */
	
	return self;
}

void _AdaptiveVoronoi_Init( void* adaptiveVoronoi, Iteration_Index maxIterations ) {
	AdaptiveVoronoi* self = (AdaptiveVoronoi*)adaptiveVoronoi;
	
	self->isConstructed = True;

	self->maxIterations          = maxIterations;
	self->cellVertexCount        = ( self->dim == 2 ? 4 : 8 );
	self->newVertexCountForSplit = ( self->dim == 2 ? 5 : 19 );
}


/*------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _AdaptiveVoronoi_Delete( void* adaptiveVoronoi ) {
	AdaptiveVoronoi* self = (AdaptiveVoronoi*)adaptiveVoronoi;
	
	/* Delete parent */
	_DiscreteVoronoi_Delete( self );

	Memory_Free( self->cells );
	Memory_Free( self->verticies );
	Memory_Free( self->cellList );
}


void _AdaptiveVoronoi_Print( void* adaptiveVoronoi, Stream* stream ) {
	AdaptiveVoronoi* self = (AdaptiveVoronoi*)adaptiveVoronoi;
	
	/* Print parent */
	_DiscreteVoronoi_Print( self, stream );
}

void* _AdaptiveVoronoi_Copy( void* adaptiveVoronoi, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	AdaptiveVoronoi*	self = (AdaptiveVoronoi*)adaptiveVoronoi;
	AdaptiveVoronoi*	newAdaptiveVoronoi;
	
	newAdaptiveVoronoi = (AdaptiveVoronoi*)_DiscreteVoronoi_Copy( self, dest, deep, nameExt, ptrMap );
	
	return (void*)newAdaptiveVoronoi;
}

void* _AdaptiveVoronoi_DefaultNew( Name name ) {
	return (void*) _AdaptiveVoronoi_New(
			sizeof(AdaptiveVoronoi),
			AdaptiveVoronoi_Type,
			_AdaptiveVoronoi_Delete,
			_AdaptiveVoronoi_Print,
			_AdaptiveVoronoi_Copy,
			_AdaptiveVoronoi_DefaultNew,
			_AdaptiveVoronoi_Construct,
			_AdaptiveVoronoi_Build,
			_AdaptiveVoronoi_Initialise,
			_AdaptiveVoronoi_Execute,
			_AdaptiveVoronoi_Destroy,
			_AdaptiveVoronoi_CalculateForCell,
			_AdaptiveVoronoi_GetParticleIndex,
			_AdaptiveVoronoi_GetVolume,
			_AdaptiveVoronoi_GetCentroid,
			name );
}


void _AdaptiveVoronoi_Construct( void* adaptiveVoronoi, Stg_ComponentFactory* cf, void* data ) {
	AdaptiveVoronoi*	     self          = (AdaptiveVoronoi*) adaptiveVoronoi;
	Iteration_Index          maxIterations;

	_DiscreteVoronoi_Construct( self, cf, data );

	maxIterations = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, "maxIterations", 3 );
	
	_AdaptiveVoronoi_Init( self, maxIterations );
}

void _AdaptiveVoronoi_Build( void* adaptiveVoronoi, void* data ) {
	AdaptiveVoronoi*	self = (AdaptiveVoronoi*)adaptiveVoronoi;

	_DiscreteVoronoi_Build( self, data );

	/* Allocate Memory */
	self->cellsAllocatedCount = 5;
	self->cells = Memory_Alloc_Array( AdaptiveVoronoiUnclaimedCell, self->cellsAllocatedCount, "cell list" );

	self->verticiesAllocatedCount = 15;
	self->verticies = Memory_Alloc_Array( AdaptiveVoronoiVertex, self->verticiesAllocatedCount, "AdaptiveVoronoiVertex" );

	self->claimedCellsAllocated = 5;
	self->cellList = Memory_Alloc_Array( AdaptiveVoronoiClaimedCell, self->claimedCellsAllocated, "claimedCellList" );
}

void _AdaptiveVoronoi_Initialise( void* adaptiveVoronoi, void* data ) {
	AdaptiveVoronoi*	self = (AdaptiveVoronoi*)adaptiveVoronoi;
	
	_DiscreteVoronoi_Initialise( self, data );
}
void _AdaptiveVoronoi_Execute( void* adaptiveVoronoi, void* data ) {
	AdaptiveVoronoi*	self = (AdaptiveVoronoi*)adaptiveVoronoi;
	
	_DiscreteVoronoi_Execute( self, data );
}
void _AdaptiveVoronoi_Destroy( void* adaptiveVoronoi, void* data ) {
	AdaptiveVoronoi*	self = (AdaptiveVoronoi*)adaptiveVoronoi;
	
	_DiscreteVoronoi_Destroy( self, data );
}

/*-------------------------------------------------------------------------------------------------------------------------
** Private Functions
*/
void _AdaptiveVoronoi_CalculateForCell( void* adaptiveVoronoi, void* _swarm, Cell_LocalIndex lCell_I ) {
	AdaptiveVoronoi*             self            = (AdaptiveVoronoi*)  adaptiveVoronoi;
	Index                        currentCellIndex;

	/* Initialise Current Info */
	self->swarm                  = (Swarm*) _swarm;
	self->cell_I                 = lCell_I;
	self->claimedCellCount       = 0;
	AdaptiveVoronoi_InitialiseFirstCell( self );

	/* Do AdaptiveVoronoi */
	currentCellIndex = 0;
	while ( currentCellIndex != FINISHED )
		currentCellIndex = AdaptiveVoronoi_Loop( self, currentCellIndex );
}

Particle_InCellIndex _AdaptiveVoronoi_GetParticleIndex( void* discreteVoronoi, Voronoi_CellIndex vCell_I ){
	AdaptiveVoronoi*             self            = (AdaptiveVoronoi*)  discreteVoronoi;

	return self->cellList[ vCell_I ].particle_I;
}

double _AdaptiveVoronoi_GetVolume(void* discreteVoronoi, Voronoi_CellIndex vCell_I ){
	AdaptiveVoronoi*             self       = (AdaptiveVoronoi*)  discreteVoronoi;
	AdaptiveVoronoiClaimedCell*  cell       = &self->cellList[ vCell_I ];

	if (self->dim == 2)
		return StGermain_ConvexQuadrilateralArea( 
			cell->vertex[0], 
			cell->vertex[1], 
			cell->vertex[2], 
			cell->vertex[3],
			self->dim); 
	else 
		return StGermain_ParallelepipedVolume(
			cell->vertex[0], 
			cell->vertex[1], 
			cell->vertex[3], 
			cell->vertex[4] );
}

void _AdaptiveVoronoi_GetCentroid(void* discreteVoronoi, Voronoi_CellIndex vCell_I, Coord centroid ){
	AdaptiveVoronoi*             self            = (AdaptiveVoronoi*)  discreteVoronoi;

	memcpy( centroid, self->cellList[ vCell_I ].centroid, self->dim * sizeof(double) );
}

/*-------------------------------------------------------------------------------------------------------------------------
** Public Functions
*/



Bool AdaptiveVoronoi_NeedSplitCell( AdaptiveVoronoi* self, AdaptiveVoronoiUnclaimedCell* cell ) {
	Vertex_Index           cellVertexCount = self->cellVertexCount;
	Vertex_Index           vertex_I;
	AdaptiveVoronoiVertex* vertex;
	AdaptiveVoronoiVertex* testVertex;

	/* Get First Vertex of Cell to test all other verticies against */
	testVertex = &self->verticies[ cell->vertex_I[ 0 ] ];

	for ( vertex_I = 0 ; vertex_I < cellVertexCount ; vertex_I++ ) {
		vertex = &self->verticies[ cell->vertex_I[ vertex_I ] ];
	
		/* Check if this guy has been claimed */
		if ( vertex->particle_I == UNCLAIMED ) 
			vertex->particle_I = Swarm_FindClosestParticleInCell( self->swarm, self->cell_I, self->dim, vertex->coord, NULL );
					
		/* Check if this vertex is claimed by the same particle as first vertex */
		if ( vertex->particle_I != testVertex->particle_I )
			return True;
		
	}
	return False;
}
	
void AdaptiveVoronoi_InitialiseFirstCell( AdaptiveVoronoi* self ) {
	AdaptiveVoronoiUnclaimedCell* firstCell = &self->cells[0];
	Vertex_Index                  cVertex_I;
	Cell_Index                    cell_I    = self->cell_I;

	/* Initialise data */
	memset( self->cells, 0,     sizeof( AdaptiveVoronoiUnclaimedCell ) * self->cellsAllocatedCount );
	memset( self->verticies, 0, sizeof( AdaptiveVoronoiVertex )        * self->verticiesAllocatedCount );

	/* Set up cell */
	firstCell->level      = 0;
	firstCell->nextCell_I = FINISHED;
	
	for ( cVertex_I = 0 ; cVertex_I < self->cellVertexCount ; cVertex_I++ ) {
		/* Assign Vertex to Cell */
		firstCell->vertex_I[ cVertex_I ] = cVertex_I;

		/* Set up position of cell */
		memcpy( self->verticies[ cVertex_I ].coord, self->swarm->cellPointTbl[cell_I][ cVertex_I ], sizeof(Coord) );

		/* Initialise particle index so we know that this guy is unclaimed */
		self->verticies[ cVertex_I ].particle_I = UNCLAIMED;
	}

	self->cellsCount     = 1;
	self->verticiesCount = self->cellVertexCount;
}


void AdaptiveVoronoi_SplitCell( AdaptiveVoronoi* self, Index parentCell_I ) {
	Cell_Index                      newCell_I;
	Vertex_Index                    cellVertexCount            = self->cellVertexCount;
	double*                         vertexList[8];
	AdaptiveVoronoiUnclaimedCell*   currentCell                = NULL;
	AdaptiveVoronoiUnclaimedCell*   newCells                   = NULL;
	AdaptiveVoronoiUnclaimedCell*   parentCell                 = NULL;
	AdaptiveVoronoiVertex*          newVerticies               = NULL;
	AdaptiveVoronoiVertex*          currentVertex              = NULL;
	AdaptiveVoronoiVertex*          vertex_LEFT_BOTTOM_FRONT   = NULL;
	AdaptiveVoronoiVertex*          vertex_RIGHT_BOTTOM_FRONT  = NULL;
	AdaptiveVoronoiVertex*          vertex_RIGHT_TOP_FRONT     = NULL;
	AdaptiveVoronoiVertex*          vertex_LEFT_TOP_FRONT      = NULL;
	Dimension_Index                 dim                        = self->dim;

	/* Allocate Memory */
	if ( self->cellsAllocatedCount < self->cellsCount + cellVertexCount ) {
		self->cellsAllocatedCount = self->cellsCount + cellVertexCount * 5;
		self->cells = Memory_Realloc_Array( self->cells, AdaptiveVoronoiUnclaimedCell, self->cellsAllocatedCount );
	}
	if ( self->verticiesAllocatedCount < self->verticiesCount + self->newVertexCountForSplit ) {
		self->verticiesAllocatedCount = self->verticiesCount + 5 * self->newVertexCountForSplit;
		self->verticies = Memory_Realloc_Array( self->verticies, AdaptiveVoronoiVertex, self->verticiesAllocatedCount );
	}

	/* Get pointers to memory */
	newCells     = &self->cells[ self->cellsCount ];
	parentCell   = &self->cells[ parentCell_I ];
	newVerticies = &self->verticies[ self->verticiesCount ];

	/* Initialise Memory */
	memset( newCells,     0, sizeof(AdaptiveVoronoiUnclaimedCell) * cellVertexCount );
	memset( newVerticies, 0, sizeof(AdaptiveVoronoiVertex)        * self->newVertexCountForSplit );

	for ( newCell_I = 0 ; newCell_I < cellVertexCount ; newCell_I++ ) {
		currentCell = &newCells[ newCell_I ];
		
		/* Increment Level */
		currentCell->level = parentCell->level + 1;

		/* Set Linked List */
		currentCell->nextCell_I = self->cellsCount + newCell_I + 1;
		
		/* Assign Vertex Which is the Same as the parent cell */
		currentCell->vertex_I[ newCell_I ] = parentCell->vertex_I[ newCell_I ];
	}
	/* Reset last guy to point to parents next cell */
	currentCell->nextCell_I = parentCell->nextCell_I;
	/* Set Parent Cell to point to first child */
	parentCell->nextCell_I  = self->cellsCount;

	/* Grab pointers to verticies of parent cell */
	vertex_LEFT_BOTTOM_FRONT  = &self->verticies[ parentCell->vertex_I[ LEFT_BOTTOM_FRONT  ] ];
	vertex_RIGHT_BOTTOM_FRONT = &self->verticies[ parentCell->vertex_I[ RIGHT_BOTTOM_FRONT ] ];
	vertex_RIGHT_TOP_FRONT    = &self->verticies[ parentCell->vertex_I[ RIGHT_TOP_FRONT    ] ];
	vertex_LEFT_TOP_FRONT     = &self->verticies[ parentCell->vertex_I[ LEFT_TOP_FRONT     ] ];
	
	/* Setup New Verticies for children cells */
	currentVertex = &newVerticies[ MIDDLE_MIDDLE_FRONT ];
	currentVertex->particle_I = UNCLAIMED;
	vertexList[0] = vertex_LEFT_BOTTOM_FRONT->coord;
	vertexList[1] = vertex_RIGHT_BOTTOM_FRONT->coord;
	vertexList[2] = vertex_RIGHT_TOP_FRONT->coord;
	vertexList[3] = vertex_LEFT_TOP_FRONT->coord;
	StGermain_AverageCoord( currentVertex->coord, vertexList, 4, dim );

	currentVertex = &newVerticies[ LEFT_MIDDLE_FRONT ];
	currentVertex->particle_I = UNCLAIMED;
	vertexList[0] = vertex_LEFT_BOTTOM_FRONT->coord;
	vertexList[1] = vertex_LEFT_TOP_FRONT->coord;
	StGermain_AverageCoord( currentVertex->coord, vertexList, 2, dim );

	currentVertex = &newVerticies[ RIGHT_MIDDLE_FRONT ];
	currentVertex->particle_I = UNCLAIMED;
	vertexList[0] = vertex_RIGHT_BOTTOM_FRONT->coord;
	vertexList[1] = vertex_RIGHT_TOP_FRONT->coord;
	StGermain_AverageCoord( currentVertex->coord, vertexList, 2, dim );

	currentVertex = &newVerticies[ MIDDLE_BOTTOM_FRONT ];
	currentVertex->particle_I = UNCLAIMED;
	vertexList[0] = vertex_LEFT_BOTTOM_FRONT->coord;
	vertexList[1] = vertex_RIGHT_BOTTOM_FRONT->coord;
	StGermain_AverageCoord( currentVertex->coord, vertexList, 2, dim );

	currentVertex = &newVerticies[ MIDDLE_TOP_FRONT ];
	currentVertex->particle_I = UNCLAIMED;
	vertexList[0] = vertex_LEFT_TOP_FRONT->coord;
	vertexList[1] = vertex_RIGHT_TOP_FRONT->coord;
	StGermain_AverageCoord( currentVertex->coord, vertexList, 2, dim );

	/* Assign New Verticies */
	currentCell = &newCells[ LEFT_BOTTOM_FRONT ];
	currentCell->vertex_I[ RIGHT_BOTTOM_FRONT ] = self->verticiesCount + MIDDLE_BOTTOM_FRONT ;
	currentCell->vertex_I[ RIGHT_TOP_FRONT    ] = self->verticiesCount + MIDDLE_MIDDLE_FRONT ;
	currentCell->vertex_I[ LEFT_TOP_FRONT     ] = self->verticiesCount + LEFT_MIDDLE_FRONT   ;

	currentCell = &newCells[ RIGHT_BOTTOM_FRONT ];
	currentCell->vertex_I[ LEFT_BOTTOM_FRONT  ] = self->verticiesCount + MIDDLE_BOTTOM_FRONT ;
	currentCell->vertex_I[ LEFT_TOP_FRONT     ] = self->verticiesCount + MIDDLE_MIDDLE_FRONT ;
	currentCell->vertex_I[ RIGHT_TOP_FRONT    ] = self->verticiesCount + RIGHT_MIDDLE_FRONT  ;

	currentCell = &newCells[ LEFT_TOP_FRONT ];
	currentCell->vertex_I[ LEFT_BOTTOM_FRONT  ] = self->verticiesCount + LEFT_MIDDLE_FRONT   ;
	currentCell->vertex_I[ RIGHT_BOTTOM_FRONT ] = self->verticiesCount + MIDDLE_MIDDLE_FRONT ;
	currentCell->vertex_I[ RIGHT_TOP_FRONT    ] = self->verticiesCount + MIDDLE_TOP_FRONT    ;

	currentCell = &newCells[ RIGHT_TOP_FRONT ];
	currentCell->vertex_I[ LEFT_BOTTOM_FRONT  ] = self->verticiesCount + MIDDLE_MIDDLE_FRONT ;
	currentCell->vertex_I[ RIGHT_BOTTOM_FRONT ] = self->verticiesCount + RIGHT_MIDDLE_FRONT  ;
	currentCell->vertex_I[ LEFT_TOP_FRONT     ] = self->verticiesCount + MIDDLE_TOP_FRONT    ;

	if (self->dim == 3){
		AdaptiveVoronoiVertex* vertex_LEFT_BOTTOM_BACK;
		AdaptiveVoronoiVertex* vertex_RIGHT_BOTTOM_BACK;
		AdaptiveVoronoiVertex* vertex_RIGHT_TOP_BACK;
		AdaptiveVoronoiVertex* vertex_LEFT_TOP_BACK;

		/* Grab pointers to verticies of parent cell */
		vertex_LEFT_BOTTOM_BACK  = &self->verticies[ parentCell->vertex_I[ LEFT_BOTTOM_BACK  ] ];
		vertex_RIGHT_BOTTOM_BACK = &self->verticies[ parentCell->vertex_I[ RIGHT_BOTTOM_BACK ] ];
		vertex_RIGHT_TOP_BACK    = &self->verticies[ parentCell->vertex_I[ RIGHT_TOP_BACK    ] ];
		vertex_LEFT_TOP_BACK     = &self->verticies[ parentCell->vertex_I[ LEFT_TOP_BACK     ] ];

		/* Setup New Verticies for children cells */
		currentVertex = &newVerticies[ LEFT_BOTTOM_MIDDLE ];
		currentVertex->particle_I = UNCLAIMED;
		vertexList[0] = vertex_LEFT_BOTTOM_FRONT->coord;
		vertexList[1] = vertex_LEFT_BOTTOM_BACK->coord;
		StGermain_AverageCoord( currentVertex->coord, vertexList, 2, dim );

		currentVertex = &newVerticies[ RIGHT_BOTTOM_MIDDLE ];
		currentVertex->particle_I = UNCLAIMED;
		vertexList[0] = vertex_RIGHT_BOTTOM_FRONT->coord;
		vertexList[1] = vertex_RIGHT_BOTTOM_BACK->coord;
		StGermain_AverageCoord( currentVertex->coord, vertexList, 2, dim );

		currentVertex = &newVerticies[ RIGHT_TOP_MIDDLE ];
		currentVertex->particle_I = UNCLAIMED;
		vertexList[0] = vertex_RIGHT_TOP_FRONT->coord;
		vertexList[1] = vertex_RIGHT_TOP_BACK->coord;
		StGermain_AverageCoord( currentVertex->coord, vertexList, 2, dim );		

		currentVertex = &newVerticies[ LEFT_TOP_MIDDLE ];
		currentVertex->particle_I = UNCLAIMED;
		vertexList[0] = vertex_LEFT_TOP_FRONT->coord;
		vertexList[1] = vertex_LEFT_TOP_BACK->coord;
		StGermain_AverageCoord( currentVertex->coord, vertexList, 2, dim );		

		currentVertex = &newVerticies[ MIDDLE_MIDDLE_MIDDLE ];
		currentVertex->particle_I = UNCLAIMED;
		vertexList[0] = vertex_LEFT_BOTTOM_FRONT->coord;
		vertexList[1] = vertex_RIGHT_BOTTOM_FRONT->coord;
		vertexList[2] = vertex_RIGHT_TOP_FRONT->coord;
		vertexList[3] = vertex_LEFT_TOP_FRONT->coord;
		vertexList[4] = vertex_LEFT_BOTTOM_BACK->coord;
		vertexList[5] = vertex_RIGHT_BOTTOM_BACK->coord;
		vertexList[6] = vertex_RIGHT_TOP_BACK->coord;
		vertexList[7] = vertex_LEFT_TOP_BACK->coord;
		StGermain_AverageCoord( currentVertex->coord, vertexList, 8, dim );				

		currentVertex = &newVerticies[ LEFT_MIDDLE_MIDDLE ];
		currentVertex->particle_I = UNCLAIMED;
		vertexList[0] = vertex_LEFT_BOTTOM_FRONT->coord;
		vertexList[1] = vertex_LEFT_TOP_FRONT->coord;
		vertexList[2] = vertex_LEFT_BOTTOM_BACK->coord;
		vertexList[3] = vertex_LEFT_TOP_BACK->coord;
		StGermain_AverageCoord( currentVertex->coord, vertexList, 4, dim );	
		
		currentVertex = &newVerticies[ RIGHT_MIDDLE_MIDDLE ];
		currentVertex->particle_I = UNCLAIMED;
		vertexList[0] = vertex_RIGHT_BOTTOM_FRONT->coord;
		vertexList[1] = vertex_RIGHT_TOP_FRONT->coord;
		vertexList[2] = vertex_RIGHT_BOTTOM_BACK->coord;
		vertexList[3] = vertex_RIGHT_TOP_BACK->coord;
		StGermain_AverageCoord( currentVertex->coord, vertexList, 4, dim );
		
		currentVertex = &newVerticies[ MIDDLE_BOTTOM_MIDDLE ];
		currentVertex->particle_I = UNCLAIMED;
		vertexList[0] = vertex_LEFT_BOTTOM_FRONT->coord;
		vertexList[1] = vertex_RIGHT_BOTTOM_FRONT->coord;
		vertexList[2] = vertex_LEFT_BOTTOM_BACK->coord;
		vertexList[3] = vertex_RIGHT_BOTTOM_BACK->coord;
		StGermain_AverageCoord( currentVertex->coord, vertexList, 4, dim );		
		
		currentVertex = &newVerticies[ MIDDLE_TOP_MIDDLE ];
		currentVertex->particle_I = UNCLAIMED;
		vertexList[0] = vertex_LEFT_TOP_FRONT->coord;
		vertexList[1] = vertex_RIGHT_TOP_FRONT->coord;
		vertexList[2] = vertex_LEFT_TOP_BACK->coord;
		vertexList[3] = vertex_RIGHT_TOP_BACK->coord;
		StGermain_AverageCoord( currentVertex->coord, vertexList, 4, dim );			

		currentVertex = &newVerticies[ MIDDLE_MIDDLE_BACK ];
		currentVertex->particle_I = UNCLAIMED;
		vertexList[0] = vertex_LEFT_BOTTOM_BACK->coord;
		vertexList[1] = vertex_RIGHT_BOTTOM_BACK->coord;
		vertexList[2] = vertex_RIGHT_TOP_BACK->coord;
		vertexList[3] = vertex_LEFT_TOP_BACK->coord;
		StGermain_AverageCoord( currentVertex->coord, vertexList, 4, dim );			

		currentVertex = &newVerticies[ LEFT_MIDDLE_BACK ];
		currentVertex->particle_I = UNCLAIMED;
		vertexList[0] = vertex_LEFT_BOTTOM_BACK->coord;
		vertexList[1] = vertex_LEFT_TOP_BACK->coord;
		StGermain_AverageCoord( currentVertex->coord, vertexList, 2, dim );				

		currentVertex = &newVerticies[ RIGHT_MIDDLE_BACK ];
		currentVertex->particle_I = UNCLAIMED;
		vertexList[0] = vertex_RIGHT_BOTTOM_BACK->coord;
		vertexList[1] = vertex_RIGHT_TOP_BACK->coord;
		StGermain_AverageCoord( currentVertex->coord, vertexList, 2, dim );			

		currentVertex = &newVerticies[ MIDDLE_BOTTOM_BACK ];
		currentVertex->particle_I = UNCLAIMED;
		vertexList[0] = vertex_LEFT_BOTTOM_BACK->coord;
		vertexList[1] = vertex_RIGHT_BOTTOM_BACK->coord;
		StGermain_AverageCoord( currentVertex->coord, vertexList, 2, dim );	
		
		currentVertex = &newVerticies[ MIDDLE_TOP_BACK ];
		currentVertex->particle_I = UNCLAIMED;
		vertexList[0] = vertex_LEFT_TOP_BACK->coord;
		vertexList[1] = vertex_RIGHT_TOP_BACK->coord;
		StGermain_AverageCoord( currentVertex->coord, vertexList, 2, dim );			
		
		/* Assign New Verticies to children cells - For Front Face Cells */
		currentCell = &newCells[ LEFT_BOTTOM_FRONT ];
		currentCell->vertex_I[ LEFT_BOTTOM_BACK   ] = self->verticiesCount + LEFT_BOTTOM_MIDDLE  ;
		currentCell->vertex_I[ RIGHT_BOTTOM_BACK  ] = self->verticiesCount + MIDDLE_BOTTOM_MIDDLE;
		currentCell->vertex_I[ RIGHT_TOP_BACK     ] = self->verticiesCount + MIDDLE_MIDDLE_MIDDLE;
		currentCell->vertex_I[ LEFT_TOP_BACK      ] = self->verticiesCount + LEFT_MIDDLE_MIDDLE  ;

		currentCell = &newCells[ RIGHT_BOTTOM_FRONT ];
		currentCell->vertex_I[ LEFT_BOTTOM_BACK   ] = self->verticiesCount + MIDDLE_BOTTOM_MIDDLE;
		currentCell->vertex_I[ RIGHT_BOTTOM_BACK  ] = self->verticiesCount + RIGHT_BOTTOM_MIDDLE ;
		currentCell->vertex_I[ RIGHT_TOP_BACK     ] = self->verticiesCount + RIGHT_MIDDLE_MIDDLE ;
		currentCell->vertex_I[ LEFT_TOP_BACK      ] = self->verticiesCount + MIDDLE_MIDDLE_MIDDLE;
		
		currentCell = &newCells[ RIGHT_TOP_FRONT ];
		currentCell->vertex_I[ LEFT_BOTTOM_BACK   ] = self->verticiesCount + MIDDLE_MIDDLE_MIDDLE;
		currentCell->vertex_I[ RIGHT_BOTTOM_BACK  ] = self->verticiesCount + RIGHT_MIDDLE_MIDDLE ;
		currentCell->vertex_I[ RIGHT_TOP_BACK     ] = self->verticiesCount + RIGHT_TOP_MIDDLE    ;		
		currentCell->vertex_I[ LEFT_TOP_BACK      ] = self->verticiesCount + MIDDLE_TOP_MIDDLE   ;		

		currentCell = &newCells[ LEFT_TOP_FRONT ];
		currentCell->vertex_I[ LEFT_BOTTOM_BACK   ] = self->verticiesCount + LEFT_MIDDLE_MIDDLE  ;
		currentCell->vertex_I[ RIGHT_BOTTOM_BACK  ] = self->verticiesCount + MIDDLE_MIDDLE_MIDDLE;
		currentCell->vertex_I[ RIGHT_TOP_BACK     ] = self->verticiesCount + MIDDLE_TOP_MIDDLE   ;
		currentCell->vertex_I[ LEFT_TOP_BACK      ] = self->verticiesCount + LEFT_TOP_MIDDLE     ;

		/* Assign New Verticies to children cells - For Back Face Cells */
		currentCell = &newCells[ LEFT_BOTTOM_BACK ];
		currentCell->vertex_I[ LEFT_BOTTOM_FRONT  ] = self->verticiesCount + LEFT_BOTTOM_MIDDLE  ;
		currentCell->vertex_I[ RIGHT_BOTTOM_FRONT ] = self->verticiesCount + MIDDLE_BOTTOM_MIDDLE;
		currentCell->vertex_I[ RIGHT_TOP_FRONT    ] = self->verticiesCount + MIDDLE_MIDDLE_MIDDLE;
		currentCell->vertex_I[ LEFT_TOP_FRONT     ] = self->verticiesCount + LEFT_MIDDLE_MIDDLE  ;
		currentCell->vertex_I[ RIGHT_BOTTOM_BACK  ] = self->verticiesCount + MIDDLE_BOTTOM_BACK  ;
		currentCell->vertex_I[ RIGHT_TOP_BACK     ] = self->verticiesCount + MIDDLE_MIDDLE_BACK  ;
		currentCell->vertex_I[ LEFT_TOP_BACK      ] = self->verticiesCount + LEFT_MIDDLE_BACK    ;

		currentCell = &newCells[ RIGHT_BOTTOM_BACK ];
		currentCell->vertex_I[ LEFT_BOTTOM_FRONT  ] = self->verticiesCount + MIDDLE_BOTTOM_MIDDLE;
		currentCell->vertex_I[ RIGHT_BOTTOM_FRONT ] = self->verticiesCount + RIGHT_BOTTOM_MIDDLE ;
		currentCell->vertex_I[ RIGHT_TOP_FRONT    ] = self->verticiesCount + RIGHT_MIDDLE_MIDDLE ;
		currentCell->vertex_I[ LEFT_TOP_FRONT     ] = self->verticiesCount + MIDDLE_MIDDLE_MIDDLE;
		currentCell->vertex_I[ LEFT_BOTTOM_BACK   ] = self->verticiesCount + MIDDLE_BOTTOM_BACK  ;
		currentCell->vertex_I[ RIGHT_TOP_BACK     ] = self->verticiesCount + RIGHT_MIDDLE_BACK   ;
		currentCell->vertex_I[ LEFT_TOP_BACK      ] = self->verticiesCount + MIDDLE_MIDDLE_BACK  ;

		currentCell = &newCells[ RIGHT_TOP_BACK ];
		currentCell->vertex_I[ LEFT_BOTTOM_FRONT  ] = self->verticiesCount + MIDDLE_MIDDLE_MIDDLE;
		currentCell->vertex_I[ RIGHT_BOTTOM_FRONT ] = self->verticiesCount + RIGHT_MIDDLE_MIDDLE ;
		currentCell->vertex_I[ RIGHT_TOP_FRONT    ] = self->verticiesCount + RIGHT_TOP_MIDDLE    ;		
		currentCell->vertex_I[ LEFT_TOP_FRONT     ] = self->verticiesCount + MIDDLE_TOP_MIDDLE   ;			
		currentCell->vertex_I[ LEFT_BOTTOM_BACK   ] = self->verticiesCount + MIDDLE_MIDDLE_BACK  ;
		currentCell->vertex_I[ RIGHT_BOTTOM_BACK  ] = self->verticiesCount + RIGHT_MIDDLE_BACK   ;
		currentCell->vertex_I[ LEFT_TOP_BACK      ] = self->verticiesCount + MIDDLE_TOP_BACK     ;

		currentCell = &newCells[ LEFT_TOP_BACK ];
		currentCell->vertex_I[ LEFT_BOTTOM_FRONT  ] = self->verticiesCount + LEFT_MIDDLE_MIDDLE  ;
		currentCell->vertex_I[ RIGHT_BOTTOM_FRONT ] = self->verticiesCount + MIDDLE_MIDDLE_MIDDLE;
		currentCell->vertex_I[ RIGHT_TOP_FRONT    ] = self->verticiesCount + MIDDLE_TOP_MIDDLE   ;
		currentCell->vertex_I[ LEFT_TOP_FRONT     ] = self->verticiesCount + LEFT_TOP_MIDDLE     ;
		currentCell->vertex_I[ LEFT_BOTTOM_BACK   ] = self->verticiesCount + LEFT_MIDDLE_BACK    ;
		currentCell->vertex_I[ RIGHT_BOTTOM_BACK  ] = self->verticiesCount + MIDDLE_MIDDLE_BACK  ;
		currentCell->vertex_I[ RIGHT_TOP_BACK     ] = self->verticiesCount + MIDDLE_TOP_BACK     ;
	}

	/* Add numbers to cells */
	self->cellsCount     += cellVertexCount;
	self->verticiesCount += self->newVertexCountForSplit;
}

void AdaptiveVoronoi_ClaimCell( 
		AdaptiveVoronoi*              self,
		AdaptiveVoronoiUnclaimedCell* cell, 
		Particle_InCellIndex          particle_I )
{
	AdaptiveVoronoiClaimedCell* claimedCell;
	double*                     centroid;
	Dimension_Index             dim             = self->dim;
	Index                       cVertex_I;
	Index                       cellVertexCount = self->cellVertexCount;
	double                      factor          = 1.0/ (double)cellVertexCount;
	AdaptiveVoronoiVertex*      vertex;

	/* Get Pointers */
	claimedCell = AdaptiveVoronoi_NewClaimedCell( self );
	centroid    = claimedCell->centroid;

	/* Find Centroid and store verticies */ 
	centroid[ I_AXIS ] = centroid[ J_AXIS ] = centroid[ K_AXIS ] = 0.0;
	for ( cVertex_I = 0 ; cVertex_I < cellVertexCount ; cVertex_I++ ) {
		vertex = &self->verticies[ cell->vertex_I[ cVertex_I ] ];

		centroid[ I_AXIS ] += vertex->coord[ I_AXIS ];
		centroid[ J_AXIS ] += vertex->coord[ J_AXIS ];
		if (dim == 3)
			centroid[ K_AXIS ] += vertex->coord[ K_AXIS ];
		
		/* Store Verticies */
		memcpy( claimedCell->vertex[ cVertex_I ] , vertex->coord, sizeof(Coord) );
	}
	centroid[ I_AXIS ] *= factor;
	centroid[ J_AXIS ] *= factor;
	if (dim == 3)
		centroid[ K_AXIS ] *= factor;

	/* Claim the cell if it isn't claimed already */
	if ( particle_I == UNCLAIMED ) 
		particle_I = Swarm_FindClosestParticleInCell( self->swarm, self->cell_I, self->dim, centroid, NULL );

	/* Add data to cell */
	claimedCell->particle_I = particle_I;
}

Voronoi_CellIndex AdaptiveVoronoi_Loop( AdaptiveVoronoi* self, Index cell_I ) {
	AdaptiveVoronoiUnclaimedCell* cell = &self->cells[ cell_I ];

	/* If we've split as many times as we are allowed - then just claim cell */
	if (cell->level >= self->maxIterations) 
		AdaptiveVoronoi_ClaimCell(self, cell, UNCLAIMED );
	else {
		/* Check to see whether we should split this cell */
		if (AdaptiveVoronoi_NeedSplitCell( self, cell )) 
			AdaptiveVoronoi_SplitCell( self, cell_I );
		else
			AdaptiveVoronoi_ClaimCell( self, cell, self->verticies[ cell->vertex_I[0] ].particle_I );
	}
	/* reassign pointer, in case there has been a realloc */
	cell = &self->cells[ cell_I ];

	return cell->nextCell_I;
}

void AdaptiveVoronoi_PrintLinkedList( AdaptiveVoronoi* self, Stream* stream ){
	Cell_Index                    cell_I = 0;
	AdaptiveVoronoiUnclaimedCell* cell;

	do {
		cell = &self->cells[ cell_I ];
		Journal_Printf( stream, "cell %u points to cell %u\n", cell_I, cell->nextCell_I );
		cell_I = cell->nextCell_I;
	} while ( cell_I != FINISHED );
	
}


AdaptiveVoronoiClaimedCell* AdaptiveVoronoi_NewClaimedCell( AdaptiveVoronoi* self ) {
	if ( self->claimedCellsAllocated <= self->claimedCellCount ) {
		self->claimedCellsAllocated = self->claimedCellCount + 5;
		self->cellList = Memory_Realloc_Array( self->cellList, AdaptiveVoronoiClaimedCell, self->claimedCellsAllocated );
	}

	self->claimedCellCount++;

	return &self->cellList[ self->claimedCellCount - 1 ];
}
