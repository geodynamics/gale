/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street, Melbourne, 3053, Australia.
**
** Authors:
**	Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**	Patrick D. Sunter, Software Engineer, VPAC. (pds@vpac.org)
**	Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**	Siew-Ching Tan, Software Engineer, VPAC. (siew@vpac.org)
**	Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
**	Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
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
*/
/** \file
**  Role:
** The ParallelDelaunay class computes the constrained delaunay triangulation of a set of points
** in 2 Dimensions along with the voronoi diagram of the point-set.
** Assumptions:
**
** Comments:
**
** The recursive divide and conquer algorithm has been implemented from the pseudo-code
** given in "Primitives for the manipulation of general subdivisions and the computation
** of voronoi diagrams" by Leonidas Guibas and Jorge Stolfi.
** ACM transtactions on graphics, Vol. 4, No. 2, April 1985, Pages 74-123
**
**
** $Id: ParallelDelaunay.c 3462 2006-02-19 06:53:24Z RaquibulHassan $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include "StGermain/StGermain.h"


#include "units.h"
#include "types.h"
#include "QuadEdge.h"
#include "Delaunay.h"
#include "ParallelDelaunay.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>

const Type ParallelDelaunay_Type="ParallelDelaunay";
#define PI 3.1415926535897932384626
#define MERGE_FACTOR 1

	/*--------------------------------------------------------------------------------------------------------------------------
	** Constructors
	*/
	
/** Create a ParallelDelaunay */
ParallelDelaunay* ParallelDelaunay_DefaultNew( Name name )
{
	/* Variables set in this function */
	SizeT                                              _sizeOfSelf = sizeof( ParallelDelaunay );
	Type                                                      type = ParallelDelaunay_Type;
	Stg_Class_DeleteFunction*                              _delete = _ParallelDelaunay_Delete;
	Stg_Class_PrintFunction*                                _print = _ParallelDelaunay_Print;
	Stg_Class_CopyFunction*                                  _copy = _ParallelDelaunay_Copy;
	Stg_Component_DefaultConstructorFunction*  _defaultConstructor = (Stg_Component_DefaultConstructorFunction*)ParallelDelaunay_DefaultNew;
	Stg_Component_ConstructFunction*                    _construct = _ParallelDelaunay_AssignFromXML;
	Stg_Component_BuildFunction*                            _build = _ParallelDelaunay_Build;
	Stg_Component_InitialiseFunction*                  _initialise = _ParallelDelaunay_Initialise;
	Stg_Component_ExecuteFunction*                        _execute = _ParallelDelaunay_Execute;
	Stg_Component_DestroyFunction*                        _destroy = _ParallelDelaunay_Destroy;
	Bool                                                  initFlag = False;
	Dictionary*                                         dictionary = NULL;
	CoordF*                                                  sites = NULL;
	int                                                   numSites = 0;
	int                                                       rank = 0;
	int                                                   numProcs = 0;
	MPI_Comm*                                                 comm = NULL;
	DelaunayAttributes*                                       attr = NULL;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;
	int                       idOffset = ZERO;

	ParallelDelaunay *d = _ParallelDelaunay_New(  PARALLELDELAUNAY_PASSARGS  );

	return d;
}
	
ParallelDelaunay* ParallelDelaunay_New(
	Name						name,
	Dictionary*					dictionary,
	CoordF						*sites,
	int							numSites,
	int							rank,
	int							numProcs,
	MPI_Comm					*comm,
	DelaunayAttributes			*attr )
{
	/* Variables set in this function */
	SizeT                                              _sizeOfSelf = sizeof( ParallelDelaunay );
	Type                                                      type = ParallelDelaunay_Type;
	Stg_Class_DeleteFunction*                              _delete = _ParallelDelaunay_Delete;
	Stg_Class_PrintFunction*                                _print = _ParallelDelaunay_Print;
	Stg_Class_CopyFunction*                                  _copy = _ParallelDelaunay_Copy;
	Stg_Component_DefaultConstructorFunction*  _defaultConstructor = (Stg_Component_DefaultConstructorFunction*)ParallelDelaunay_DefaultNew;
	Stg_Component_ConstructFunction*                    _construct = _ParallelDelaunay_AssignFromXML;
	Stg_Component_BuildFunction*                            _build = _ParallelDelaunay_Build;
	Stg_Component_InitialiseFunction*                  _initialise = _ParallelDelaunay_Initialise;
	Stg_Component_ExecuteFunction*                        _execute = _ParallelDelaunay_Execute;
	Stg_Component_DestroyFunction*                        _destroy = _ParallelDelaunay_Destroy;
	Bool                                                  initFlag = True;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;
	int                       idOffset = ZERO;

	ParallelDelaunay *d = _ParallelDelaunay_New(  PARALLELDELAUNAY_PASSARGS  );
	
	return d;
}

#define MASTER_PROC 0

	/** Initialise a ParallelDelaunay */
void ParallelDelaunay_Init(
	ParallelDelaunay*					self,
	Name						name,
	Dictionary*					dictionary,
	CoordF						*sites,
	int							numSites,
	int							rank,
	int							numProcs,
	MPI_Comm					*comm,
	DelaunayAttributes			*attr )
{
	self->type = ParallelDelaunay_Type;
	self->_sizeOfSelf = sizeof( ParallelDelaunay );
	self->_deleteSelf = False;
	self->dictionary = dictionary;

	self->_delete = _ParallelDelaunay_Delete;
	self->_print = _ParallelDelaunay_Print;
	self->_copy = _ParallelDelaunay_Copy;
	self->_defaultConstructor = (Stg_Component_DefaultConstructorFunction*)ParallelDelaunay_DefaultNew;
	self->_construct = _ParallelDelaunay_AssignFromXML;
	self->_build = _ParallelDelaunay_Build;
	self->_initialise = _ParallelDelaunay_Initialise;
	self->_execute = _ParallelDelaunay_Execute;
	self->_destroy = _ParallelDelaunay_Destroy;

	self->attributes = Memory_Alloc_Unnamed( DelaunayAttributes );
	memcpy( self->attributes, attr, sizeof( DelaunayAttributes ) );
	self->attributes->BuildBoundingTriangle = 0;

	_Stg_Class_Init( (Stg_Class*)self );
	_Stg_Object_Init( (Stg_Object*)self, name, NON_GLOBAL );
	_Stg_Component_Init( (Stg_Component*)self );
	_ParallelDelaunay_Init( self, sites, 0, 0, rank, numProcs, comm, numSites, numSites, dictionary, True );
}

/** Creation implementation */
ParallelDelaunay* _ParallelDelaunay_New(  PARALLELDELAUNAY_DEFARGS  )
{
	DelaunayAttributes* attr_shadowed = attr;
	Bool initFlag_shadowed = initFlag;
	ParallelDelaunay *self = NULL;
	DelaunayAttributes *myAttr = NULL;
	
	myAttr = Memory_Alloc_Unnamed( DelaunayAttributes );
	memcpy( myAttr, attr_shadowed, sizeof( DelaunayAttributes ) );
	myAttr->BuildBoundingTriangle = 0;
	
	assert( _sizeOfSelf >= sizeof(ParallelDelaunay) );
	/* The following terms are parameters that have been passed into this function but are being set before being passed onto the parent */
	/* This means that any values of these parameters that are passed into this function are not passed onto the parent function
	   and so should be set to ZERO in any children of this class. */
	idOffset = 0;
	initFlag = True;
	attr     = myAttr;

	self = (ParallelDelaunay*)_Delaunay_New(  DELAUNAY_PASSARGS  );
	
	_ParallelDelaunay_Init( self, sites, 0, 0, rank, numProcs, comm, numSites, numSites, dictionary, initFlag_shadowed );

	return self;
}

int ParallelDelaunayBtreeCompareFunction( void *a, void *b )
{
	Site *s1, *s2;

	s1 = (Site*)a;
	s2 = (Site*)b;

	if( (*(s1->coord))[0] > (*(s2->coord))[0] ){
		return 1;
	}
	else if( (*(s1->coord))[0] == (*(s2->coord))[0] ){
		if( (*(s1->coord))[1] > (*(s2->coord))[1] ){
			return 1;
		}
		else if( (*(s1->coord))[1] < (*(s2->coord))[1] ){
			return -1;
		}
		else{
			return 0;
		}
	}
	else{
		return -1;
	}
}

#define epsilon 0.0001
#define LOAD_TAG 1
#define DATA_TAG 1<<1
#define NEW_SITE_FACTOR 100
void _ParallelDelaunay_Init( ParallelDelaunay* self, CoordF* points, int leftProc, int rightProc, int rank, int numProcs, MPI_Comm* comm, int numSites, int numInputSites, Dictionary* dictionary, Bool initFlag )
{
	assert( self );

	self->points = points;
	self->leftProc = leftProc;
	self->rightProc = rightProc;
	self->haloSites[0] = NULL;
	self->haloSites[1] = NULL;
	self->localTriangulation = NULL;
	self->rank = rank;
	self->numProcs = numProcs;
	self->comm = comm;

	self->dictionary = dictionary;
	self->numSites = numSites;
	self->numInputSites = numInputSites;
	
    if (initFlag) {
	    self->numHaloSites[0] = 0;
	    self->numHaloSites[1] = 0;
	    self->sitePool = MemoryPool_New( Site, NEW_SITE_FACTOR, NEW_SITE_FACTOR );
	    self->coordPool = MemoryPool_New( CoordF, NEW_SITE_FACTOR, NEW_SITE_FACTOR );
    }
}

	/*--------------------------------------------------------------------------------------------------------------------------
	** Virtual functions
	*/

/** Stg_Class_Delete() implementation */
void _ParallelDelaunay_Delete( void* pd )
{
	ParallelDelaunay *self = (ParallelDelaunay*)pd;
	
	assert( self );

	_Delaunay_Delete( self );
}

/** Stg_Class_Print() implementation */
void _ParallelDelaunay_Print( void* pd, Stream* stream )
{
	ParallelDelaunay *self = ( ParallelDelaunay* )pd;
	
	assert( self );
	assert( stream );

	_Delaunay_Print( (Delaunay*)self, stream );
	Journal_Printf( stream, "ParallelDelaunay (ptr): (%p)\n", self );
}

void *_ParallelDelaunay_Copy( const void* pd, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap )
{
	return NULL;
}

void _ParallelDelaunay_AssignFromXML( void* pd, Stg_ComponentFactory* cf, void* data )
{
	
}

void _ParallelDelaunay_Build( void* pd, void* data ) {
	float _minX, _maxX;
	float _minY, _maxY;
	int numProcs, numSites, i, j, count;
	float stride, start;
	int *alloced = NULL;
	int offset;
	ParallelDelaunay *self = (ParallelDelaunay*)pd;

	assert( self );
	
	numProcs = self->numProcs;
	numSites = self->numSites;

	if( numProcs == 1 ){
		self->leftProc = numProcs;
		self->rightProc = numProcs;
	}
	else{
		if( self->rank == MASTER_PROC ){
			self->leftProc = numProcs;
			self->rightProc = self->rank + 1;
		}
		else if( self->rank == (numProcs-1) ){
			self->leftProc = self->rank - 1;
			self->rightProc = numProcs;
		}
		else{
			self->leftProc = self->rank - 1;
			self->rightProc = self->rank + 1;
		}
	}
	
	self->mapGlobalToLocal = Memory_Alloc_Array_Unnamed( int, numSites );
	self->processorLoad = Memory_Alloc_Array_Unnamed( int, numProcs );
	memset( self->processorLoad, 0, sizeof( int )*numProcs );

	if( self->rank == MASTER_PROC ){
		self->processor = Memory_Alloc_Array_Unnamed( int, numSites );

		alloced = Memory_Alloc_Array_Unnamed( int, numSites );
		memset( alloced, 0, sizeof( int )*numSites );
	
		self->initialOrder = Memory_Alloc_Array_Unnamed( int, numSites );
		memset( self->initialOrder, 0, sizeof( int )*numSites );
		
		Delaunay_FindMinMax( self->sites, self->numSites, &_minX, &_minY, &_maxX, &_maxY );
		Delaunay_SortSites( self->sites, self->numSites );

		for( i=0; i<numSites; i++ ){
			self->initialOrder[i] = self->sites[i].id;
		}
		
		stride = (_maxX - _minX)/((float)numProcs);

		start = _minX;
		for( i=0; i<numProcs; i++ ){
			for( j=0; j<numSites; j++ ){
				if( ((*(self->sites[j].coord))[0] >= start-epsilon) &&
						((*(self->sites[j].coord))[0] <= (start+stride+epsilon)) &&
						(!alloced[j]) ){
					
					alloced[j] = 1;
					self->processorLoad[i]++;
					self->processor[j] = i;
				}
			}
			start+=stride;
		}
		
		/*for( i=0; i<numProcs; i++ ){
			printf( "processorLoad[%d] = %d\n", i, self->processorLoad[i] );
		}*/
		
		for( i=MASTER_PROC+1; i<numProcs; i++ ){
			MPI_Send( &(self->processorLoad[i]), 1, MPI_INT, i, LOAD_TAG, *self->comm );
		}

		self->numLocalSites = self->processorLoad[MASTER_PROC];
		self->localPoints = Memory_Alloc_Array_Unnamed( CoordF, self->numLocalSites );
		
		{
			CoordF **procCoords = NULL;
			int *procCoordCounter = NULL;

			procCoordCounter = (int*)malloc( sizeof(int)*numProcs );
			memset( procCoordCounter, 0, sizeof(int)*numProcs );

			procCoords = (float (**)[3])malloc( sizeof( CoordF* ) * numProcs );
			for( i=MASTER_PROC+1; i<numProcs; i++ ){
                          procCoords[i] = (float (*)[3])malloc( sizeof(CoordF) * self->processorLoad[i] );
				memset( procCoords[i], 0, sizeof( CoordF ) * self->processorLoad[i] );
			}

			count = 0;
			for( i=0; i<numSites; i++ ){
				if( self->processor[i] == MASTER_PROC ){
					memcpy( &(self->localPoints[count++]), self->sites[i].coord, sizeof(CoordF) );
				}
				else{
					memcpy( &(procCoords[self->processor[i]][procCoordCounter[self->processor[i]]++]), self->sites[i].coord,  sizeof(CoordF) );
				}			
			}
			
			for( i=MASTER_PROC+1; i<numProcs; i++ ){
				MPI_Send( procCoords[i], sizeof( CoordF )*self->processorLoad[i], MPI_BYTE, i, DATA_TAG, *self->comm );
			}

			free( procCoordCounter );
			for( i=MASTER_PROC+1; i<numProcs; i++ ){
				free( procCoords[i] );
			}
			free( procCoords );
		}

		Memory_Free( alloced );
	}
	else{
		MPI_Status status;
		
		MPI_Recv( &self->numLocalSites, 1, MPI_INT, MASTER_PROC, LOAD_TAG, *self->comm, &status );

		self->localPoints = Memory_Alloc_Array_Unnamed( CoordF, self->numLocalSites );

		MPI_Recv( self->localPoints, sizeof(CoordF)*self->numLocalSites, MPI_BYTE, MASTER_PROC, DATA_TAG, *self->comm, &status );
	}

	MPI_Bcast( self->processorLoad, numProcs, MPI_INT, MASTER_PROC, *self->comm );
	self->numTotalLocalSites = self->numLocalSites;

	offset = 0;
	for( i=0; i<self->rank; i++ ){
		offset += self->processorLoad[i];
	}
	
	self->localTriangulation = Delaunay_New( "delaunay", self->dictionary, self->localPoints, self->numLocalSites, offset, self->attributes );
	self->localTriangulation->qp = MemoryPool_New( QuadEdge, self->localTriangulation->numSites * (3 + MERGE_FACTOR), 10 );
	Delaunay_SortSites(self->localTriangulation->sites, self->localTriangulation->numSites);
    Delaunay_Recurse(self->localTriangulation, 0, self->localTriangulation->numSites,
			&self->localTriangulation->leftMost, &self->localTriangulation->rightMost);

	for( i=0; i<numSites; i++ ){
		self->mapGlobalToLocal[i] = numSites;
	}
	
	for( i=0; i<self->numLocalSites; i++ ){
		self->mapGlobalToLocal[self->localTriangulation->sites[i].id] = self->localTriangulation->sites[i].id;
	}
	
	if( self->leftProc != self->numProcs ){
		ParallelDelaunayMerge( self, self->comm, self->leftProc );
	}
		
	if( self->rightProc != self->numProcs ){
		ParallelDelaunayMerge( self, self->comm, self->rightProc );
	}
}

void _ParallelDelaunay_Initialise( void* pd, void* data )
{
	
}

void _ParallelDelaunay_Execute( void* pd, void* data )
{
	DelaunayAttributes attr;
	ParallelDelaunay *self = (ParallelDelaunay*)pd;
	int offset, i;
	
	offset = 0;
	for( i=0; i<self->rank; i++ ){
		offset += self->processorLoad[i];
	}

	if( self->leftProc != self->numProcs ){
		ParallelDelaunay_RetriangulateAfterMerge( self, LEFT );
	}

	if( self->rightProc != self->numProcs ){
		ParallelDelaunay_RetriangulateAfterMerge( self, RIGHT );
	}

	Stg_Class_Delete( self->localTriangulation );
	
	memcpy( &attr, self->attributes, sizeof( DelaunayAttributes ) );
	attr.BuildTriangleIndices = 0;
	self->localTriangulation = Delaunay_New( "delaunay", self->dictionary, self->localPoints, self->numTotalLocalSites, offset, &attr );
	Stg_Component_Build( self->localTriangulation, NULL, True );
	
	if( self->attributes->BuildTriangleIndices ){
		ParallelDelaunay_BuildTriangleIndices( self );
	}	
}

void _ParallelDelaunay_Destroy( void* pd, void* data ) { 
   ParallelDelaunay *self = (ParallelDelaunay*)pd;
	
	assert( self );
   Memory_Free( self->localPoints );
	if( self->mappingTable[0] ) Memory_Free( self->mappingTable[0] );
	if( self->mappingTable[1] ) Memory_Free( self->mappingTable[1] );
	Memory_Free( self->mapGlobalToLocal );
	if( self->processor ) Memory_Free( self->processor );
	if( self->initialOrder ) Memory_Free( self->initialOrder );
	Memory_Free( self->processorLoad );
	Memory_Free( self->attributes );

	Stg_Class_Delete( self->localTriangulation );
	Stg_Class_Delete( self->sitePool );
	Stg_Class_Delete( self->coordPool );

   _Delaunay_Destroy(pd, data);
}

#define onCurrentProc( pd, id ) ( id < (pd->processorLoad[pd->rank]+pd->localTriangulation->idOffset) )
#define onLeftProc( pd, id ) ( ( pd->numHaloSites[0]>0 ) && ( id >= (pd->processorLoad[pd->rank]+pd->localTriangulation->idOffset) ) && \
				( id < (pd->processorLoad[pd->rank]+pd->numHaloSites[0]+pd->localTriangulation->idOffset) ) )
#define onRightProc( pd, id ) ( ( pd->numHaloSites[1]>0 ) && ( id >= (pd->processorLoad[pd->rank]+pd->localTriangulation->idOffset+pd->numHaloSites[0]) ) && \
				( id < (pd->numTotalLocalSites+pd->localTriangulation->idOffset) ) )

void ParallelDelaunay_RetriangulateAfterMerge( ParallelDelaunay *pd, Direction d )
{
		BTreeIterator *iter = NULL;
		Site *result = NULL;
		ParallelDelaunay *self = pd;
		int i = 0;
		int numSites;
		int count;
		int offset;

		numSites = self->numSites;
		offset = 0;
		for( i=0; i<self->rank; i++ ){
			offset += self->processorLoad[i];
		}
		
		self->haloSites[d] = BTree_New( ParallelDelaunayBtreeCompareFunction, NULL, NULL, NULL, BTREE_NO_DUPLICATES );
		
		for( i=0; i<self->localTriangulation->qp->numElements; i++ ){
			if( IS_FREE((QuadEdgeRef)&(self->localTriangulation->qp->chunks[0].memory[i*sizeof(QuadEdge)])) ) continue;
	
			if( self->mapGlobalToLocal[((Site*)((QuadEdge*)&(self->localTriangulation->qp->chunks[0].memory[i*sizeof(QuadEdge)]))->data[0])->id] == numSites )
				BTree_InsertNode( self->haloSites[d], ((Site*)((QuadEdge*)&(self->localTriangulation->qp->chunks[0].memory[i*sizeof(QuadEdge)]))->data[0]), sizeof( Site* ) );
			
			if( self->mapGlobalToLocal[((Site*)((QuadEdge*)&(self->localTriangulation->qp->chunks[0].memory[i*sizeof(QuadEdge)]))->data[2])->id] == numSites )
				BTree_InsertNode( self->haloSites[d], ((Site*)((QuadEdge*)&(self->localTriangulation->qp->chunks[0].memory[i*sizeof(QuadEdge)]))->data[2]), sizeof( Site* ) );
		}

		self->localPoints = Memory_Realloc_Array
			( self->localPoints, CoordF, sizeof(CoordF) * (self->numTotalLocalSites + self->haloSites[d]->nodeCount) );
		
		self->mappingTable[d] = Memory_Alloc_Array_Unnamed( int, self->haloSites[d]->nodeCount );
		memset( self->mappingTable[d], 0, sizeof( int ) * self->haloSites[d]->nodeCount );
		
		count=0;
		i = self->numTotalLocalSites;
		iter = BTreeIterator_New( self->haloSites[d] );
		for( result=(Site*)BTreeIterator_First(iter);
				result;
				result=(Site*)BTreeIterator_Next(iter)){
			
			self->mappingTable[d][count++] = result->id;
			memcpy( &(self->localPoints[i++]), result->coord, sizeof( CoordF ) );
		}
		Stg_Class_Delete( self->localTriangulation );

		self->localTriangulation = Delaunay_New( "delaunay", self->dictionary, self->localPoints, self->numLocalSites, offset, self->attributes );
		self->localTriangulation->qp = MemoryPool_New( QuadEdge, self->localTriangulation->numSites * 4, 10 );
		Delaunay_SortSites(self->localTriangulation->sites, self->localTriangulation->numSites);
    	Delaunay_Recurse(self->localTriangulation, 0, self->localTriangulation->numSites,
				&self->localTriangulation->leftMost, &self->localTriangulation->rightMost);

		self->numHaloSites[d] = self->haloSites[d]->nodeCount;
		self->numTotalLocalSites += self->haloSites[d]->nodeCount;
		
		Stg_Class_Delete( self->haloSites[d] );
		Stg_Class_Delete( iter );
}

void ParallelDelaunay_BuildTriangleIndices( ParallelDelaunay *pd )
{
	int i = 0, triCount;
	QuadEdgeRef e = 0, eStart = 0, eOnext = 0, eLnext = 0;
	QuadEdge *edges = NULL;
	int maxEdges = 0;
	int pass = 0;
	unsigned int **triIndices = NULL;
	Delaunay *delaunay = NULL;

	delaunay = pd->localTriangulation;

	assert( delaunay );

	delaunay->triangleIndices = Memory_Alloc_Array_Unnamed( unsigned int*, delaunay->numFaces );
	delaunay->triangleIndices[0] = Memory_Alloc_Array_Unnamed( unsigned int, delaunay->numFaces * 3 );
	memset( delaunay->triangleIndices[0] , 0, sizeof(unsigned int) * delaunay->numFaces * 3 );

	for( i=0; i<delaunay->numFaces; i++ ){
		delaunay->triangleIndices[i] = delaunay->triangleIndices[0]+i*3;
	}
	
	triIndices = delaunay->triangleIndices;
	
	edges = (QuadEdge*)delaunay->qp->chunks[0].memory;
	maxEdges = delaunay->qp->numElements;
		
	for (i = 0; i < maxEdges; i++) {
		edges[i].count = 0;
	}

	triCount = 0;
	for (i = 0; i < maxEdges; i++) {
		
		e = eStart = (QuadEdgeRef)((void*)&(edges[i]));
		
		if( IS_FREE(e) )continue;
		
		do{
			eOnext = ONEXT(e);
			eLnext = LNEXT(e);
			
			if( (COUNT(e)<2) && (COUNT(LNEXT(e))<2) && (COUNT(eOnext)<2) ){
				if( ((((Site*)ORG(eLnext)) == ((Site*)DEST(e)))) &&
				 ((((Site*)DEST(eLnext)) == ((Site*)DEST(eOnext)))) ){
		
					if( onCurrentProc( pd, ((Site*)ORG(e))->id ) ||
						onCurrentProc( pd, ((Site*)DEST(e))->id ) ||
						onCurrentProc( pd, ((Site*)DEST(eOnext))->id ) ){
						
						pass = 0;
						if( !onCurrentProc( pd, ((Site*)ORG(e))->id ) ){
							if( !onRightProc( pd, ((Site*)ORG(e))->id ) ) pass = 1;
						}
						
						if( !onCurrentProc( pd, ((Site*)DEST(e))->id ) ){
							if( !onRightProc( pd, ((Site*)DEST(e))->id ) ) pass = 1;
						}
						
						if( !onCurrentProc( pd, ((Site*)DEST(eOnext))->id ) ){
							if( !onRightProc( pd, ((Site*)DEST(eOnext))->id ) ) pass = 1;
						}
						
						if( !pass ){
							triIndices[triCount][0] = ParallelDelaunay_TranslateLocalToGlobal(pd, ((Site*)ORG(e))->id);
							triIndices[triCount][1] = ParallelDelaunay_TranslateLocalToGlobal(pd, ((Site*)DEST(e))->id);
							triIndices[triCount][2] = ParallelDelaunay_TranslateLocalToGlobal(pd, ((Site*)DEST(eOnext))->id);

							triCount++;
						}
					}

					COUNT(e)++;
					COUNT(LNEXT(e))++;
					COUNT(eOnext)++;
				}
			}
			e = eOnext;
		}while( e != eStart );
	}
	
	delaunay->numTriangles = triCount;
}

#define NEIGHBOURS_TAG 1<<4
#define VORONOI_SIDES_TAG 1<<5
#define VORONOI_AREA_TAG 1<<6
#define NUM_NEIGHBOUR_TAG 1<<7
#define MAX_NEIGHBOURS 100

void ParallelDelaunay_GatherTriangulation( ParallelDelaunay *pd )
{
	int i, j, k, l, count, count1;
	MPI_Status st;
	int numNeighboursSum = 0;
	int stride = 0;
	
	assert( pd );

	if( pd->rank == MASTER_PROC ){
		
		if( pd->attributes->CalculateVoronoiSurfaceArea ){
			pd->voronoiArea = Memory_Alloc_Array_Unnamed( float, pd->numInputSites );
		}

		if( pd->attributes->FindNeighbours ){
			pd->numNeighbours = Memory_Alloc_Array_Unnamed( unsigned int, pd->numInputSites );
		}
		
		count = 0;
		count1 = 0;
		numNeighboursSum = 0;

		{
			float **procVoronoi = NULL;
			int **procNumNeighbours = NULL;
			int *procVoronoiCounter = NULL;
			int *procNumNeighboursCounter = NULL;

			procVoronoiCounter = (int*)malloc( sizeof( int ) * pd->numProcs );
			memset( procVoronoiCounter, 0, sizeof( int ) * pd->numProcs );
			
			procNumNeighboursCounter = (int*)malloc( sizeof( int ) * pd->numProcs );
			memset( procNumNeighboursCounter, 0, sizeof( int ) * pd->numProcs );

			procVoronoi = (float**)malloc( sizeof( float* ) * pd->numProcs ); memset( procVoronoi, 0, sizeof(float*)*pd->numProcs );
			procNumNeighbours = (int**)malloc( sizeof( int* ) * pd->numProcs ); memset( procNumNeighbours, 0, sizeof(int*)*pd->numProcs );

			for( i=MASTER_PROC+1; i<pd->numProcs; i++ ){

				if( pd->attributes->CalculateVoronoiSurfaceArea ){
                                  procVoronoi[i] = (float*)malloc( sizeof(float)*pd->processorLoad[i] );
					memset( procVoronoi[i], 0, sizeof( float ) * pd->processorLoad[i] );
					MPI_Recv( procVoronoi[i], pd->processorLoad[i], MPI_FLOAT, i, VORONOI_AREA_TAG, (*pd->comm), &st );
				}

				if( pd->attributes->FindNeighbours ){
                                  procNumNeighbours[i] = (int*)malloc( sizeof(int)*pd->processorLoad[i] );
					memset( procNumNeighbours[i], 0, sizeof(int) * pd->processorLoad[i] );
					MPI_Recv( procNumNeighbours[i], pd->processorLoad[i], MPI_INT, i, NUM_NEIGHBOUR_TAG, (*pd->comm), &st );
				}
			}

			for( i=0; i<pd->numInputSites; i++ ){
				int proc, order;
				proc = pd->processor[i];
				order = pd->initialOrder[i];

				if( pd->attributes->CalculateVoronoiSurfaceArea ){
					if( pd->processor[i] == MASTER_PROC ){
						memcpy( &(pd->voronoiArea[pd->initialOrder[i]]), &(pd->localTriangulation->voronoiArea[count++]), sizeof( float ) );
					}
					else{
						memcpy( &(pd->voronoiArea[order]), &(procVoronoi[proc][procVoronoiCounter[proc]++]), sizeof(float) );
					}
				}
			
				if( pd->attributes->FindNeighbours ){
					if( pd->processor[i] == MASTER_PROC ){
						memcpy( &(pd->numNeighbours[pd->initialOrder[i]]), &(pd->localTriangulation->numNeighbours[count1++]), sizeof( int ) );
					}
					else{
						memcpy( &(pd->numNeighbours[pd->initialOrder[i]]), &(procNumNeighbours[proc][procNumNeighboursCounter[proc]++]), sizeof( int ) );
					}
					numNeighboursSum += pd->numNeighbours[pd->initialOrder[i]];
				}
			}
			
			for( i=MASTER_PROC+1; i<pd->numProcs; i++ ){
				if( procVoronoi[i] ) free( procVoronoi[i] );
				if( procNumNeighbours[i] ) free( procNumNeighbours[i] );
			}
			free( procVoronoi );
			free( procNumNeighbours );
			free( procVoronoiCounter );
			free( procNumNeighboursCounter );
		}

		if( pd->attributes->FindNeighbours ){
			pd->neighbours = Memory_Alloc_Array_Unnamed( unsigned int*, pd->numInputSites );
			pd->neighbours[0] = Memory_Alloc_Array_Unnamed( unsigned int, numNeighboursSum );
		}
	
		if( pd->attributes->CalculateVoronoiSides ){
			pd->voronoiSides = Memory_Alloc_Array_Unnamed( float*, pd->numInputSites );
			pd->voronoiSides[0] = Memory_Alloc_Array_Unnamed( float, numNeighboursSum );
		}
		
		stride = 0;
		for( j=0; j<pd->numInputSites; j++ ){
			if( pd->attributes->FindNeighbours ){
				pd->neighbours[j] = pd->neighbours[0]+stride;
			}
			
			if( pd->attributes->CalculateVoronoiSides ){
				pd->voronoiSides[j] = pd->voronoiSides[0]+stride;
			}

			if( pd->attributes->FindNeighbours ){
				stride += pd->numNeighbours[j];
			}
		}
		
		{
			float **procVoronoiSides = NULL;
			float **procNeighbours = NULL;
			int *procVoronoiSidesCounter = NULL;
			int *procNeighboursCounter = NULL;
			int *procNumNeighboursCount = NULL;

			procVoronoiSidesCounter = (int*)malloc( sizeof( int ) * pd->numProcs );
			memset( procVoronoiSidesCounter, 0, sizeof( int ) * pd->numProcs );
			
			procNeighboursCounter = (int*)malloc( sizeof( int ) * pd->numProcs );
			memset( procNeighboursCounter, 0, sizeof( int ) * pd->numProcs );
			
			procNumNeighboursCount = (int*)malloc( sizeof( int ) * pd->numProcs );
			memset( procNumNeighboursCount, 0, sizeof( int ) * pd->numProcs );

			procVoronoiSides = (float**)malloc( sizeof( float* ) * pd->numProcs ); memset( procVoronoiSides, 0, sizeof(float*)*pd->numProcs );
			procNeighbours = (float**)malloc( sizeof( unsigned int* ) * pd->numProcs ); memset( procNeighbours, 0, sizeof(unsigned int*)*pd->numProcs );

			for( i=0; i<pd->numInputSites; i++ ){
				procNumNeighboursCount[pd->processor[i]]+=pd->numNeighbours[pd->initialOrder[i]];
			}

			for( i=MASTER_PROC+1; i<pd->numProcs; i++ ){
				if( pd->attributes->CalculateVoronoiSides ){
                                  procVoronoiSides[i] = (float*)malloc( sizeof(float)*procNumNeighboursCount[i] );
					memset( procVoronoiSides[i], 0, sizeof( float ) * procNumNeighboursCount[i] );
					MPI_Recv( procVoronoiSides[i], procNumNeighboursCount[i], MPI_FLOAT, i, VORONOI_SIDES_TAG, (*pd->comm), &st );
				}

				if( pd->attributes->FindNeighbours ){
                                  procNeighbours[i] = (float*)malloc( sizeof(unsigned int)*procNumNeighboursCount[i] );
					memset( procNeighbours[i], 0, sizeof(unsigned int) * procNumNeighboursCount[i] );
					MPI_Recv( procNeighbours[i], procNumNeighboursCount[i], MPI_INT, i, NEIGHBOURS_TAG, (*pd->comm), &st );
				}
			}

			count = 0;
			for(j=0; j<pd->numInputSites; j++){
				int proc, order;
				proc = pd->processor[j];
				order = pd->initialOrder[j];
				
				if( pd->processor[j] == MASTER_PROC ){
					if( pd->attributes->CalculateVoronoiSides ){
						memcpy( (pd->voronoiSides[pd->initialOrder[j]]), (pd->localTriangulation->voronoiSides[count]),
								sizeof(float)*pd->localTriangulation->numNeighbours[count] );
					}
					if( pd->attributes->FindNeighbours ){
						memcpy( (pd->neighbours[pd->initialOrder[j]]), (pd->localTriangulation->neighbours[count]),
								sizeof(int)*pd->localTriangulation->numNeighbours[count] );
						for( i=0; i<(int)(pd->numNeighbours[pd->initialOrder[j]]); i++ ){
							pd->neighbours[pd->initialOrder[j]][i] = 
								pd->initialOrder[ParallelDelaunay_TranslateLocalToGlobal(pd, pd->neighbours[pd->initialOrder[j]][i])];
						}
					}
					count++;
				}
				else{
					if( pd->attributes->CalculateVoronoiSides ){
						memcpy( (pd->voronoiSides[pd->initialOrder[j]]),
								&(procVoronoiSides[proc][procVoronoiSidesCounter[proc]]),
								sizeof(float)*pd->numNeighbours[order] );
						procVoronoiSidesCounter[proc]+=pd->numNeighbours[order];
					}
					if( pd->attributes->FindNeighbours ){
						memcpy( (pd->neighbours[pd->initialOrder[j]]),
								&(procNeighbours[proc][procNeighboursCounter[proc]]),
								sizeof(int)*pd->numNeighbours[order] );
						procNeighboursCounter[proc]+=pd->numNeighbours[order];

						for( i=0; i<(int)(pd->numNeighbours[pd->initialOrder[j]]); i++ ){
							pd->neighbours[pd->initialOrder[j]][i] = pd->initialOrder[pd->neighbours[pd->initialOrder[j]][i]];
						}
					}
				}
			}

			for( i=MASTER_PROC+1; i<pd->numProcs; i++ ){
				if( procVoronoiSides[i] ) free( procVoronoiSides[i] );
				if( procNeighbours[i] ) free( procNeighbours[i] );
			}
			free( procVoronoiSides );
			free( procNeighbours );
			free( procVoronoiSidesCounter );
			free( procNeighboursCounter );			
			free( procNumNeighboursCount );			
		}
	}
	else{
		if( pd->attributes->CalculateVoronoiSurfaceArea ){
			MPI_Send( (pd->localTriangulation->voronoiArea), pd->numLocalSites, MPI_FLOAT, MASTER_PROC, VORONOI_AREA_TAG, (*pd->comm) );
		}
			
		if( pd->attributes->FindNeighbours ){
			MPI_Send( (pd->localTriangulation->numNeighbours), pd->numLocalSites, MPI_INT, MASTER_PROC, NUM_NEIGHBOUR_TAG, (*pd->comm) );
		}
		
		if( pd->attributes->CalculateVoronoiSides ){
			int sum = 0;
			for( i=0; i<pd->numLocalSites; i++ ){
				sum+=pd->localTriangulation->numNeighbours[i];
			}
			MPI_Send( pd->localTriangulation->voronoiSides[0], sum, 
					MPI_FLOAT, MASTER_PROC, VORONOI_SIDES_TAG, *(pd->comm));
		}

		if( pd->attributes->FindNeighbours ){
			int sum = 0;
			for( i=0; i<pd->numLocalSites; i++ ){
				sum+=pd->localTriangulation->numNeighbours[i];
				
				for( j=0; j<(int)(pd->localTriangulation->numNeighbours[i]); j++ ){
					pd->localTriangulation->neighbours[i][j] = ParallelDelaunay_TranslateLocalToGlobal( pd, pd->localTriangulation->neighbours[i][j] );
				}
			}
			MPI_Send( pd->localTriangulation->neighbours[0], sum, MPI_INT, MASTER_PROC, NEIGHBOURS_TAG, *(pd->comm) );
		}
	}
	
	if( pd->attributes->BuildTriangleIndices ){
		unsigned int **triIndices = NULL;
		int globalNumTriangles = 0;
		Delaunay *delaunay = NULL;
		int *triCountArray = NULL;
		int triCount = 0;
		int rank = 0;

		delaunay = pd->localTriangulation;

		assert( delaunay );
	
		triIndices = delaunay->triangleIndices;
		triCount = delaunay->numTriangles;

		MPI_Allreduce( &triCount, &globalNumTriangles, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
	
		pd->numTriangles = globalNumTriangles;
		
		triCountArray = Memory_Alloc_Array_Unnamed( int, pd->numProcs );
		rank = pd->rank;
	
		if( rank == MASTER_PROC ){

			pd->triangleIndices = Memory_Alloc_Array_Unnamed( unsigned int*,  globalNumTriangles );
			pd->triangleIndices[0] = Memory_Alloc_Array_Unnamed( unsigned int, globalNumTriangles*3 );

			if ( pd->attributes->BuildTriangleNeighbours ){
				pd->triangleNeighbours = Memory_Alloc_Array_Unnamed( unsigned int*,  globalNumTriangles );
				pd->triangleNeighbours[0] = Memory_Alloc_Array_Unnamed( unsigned int, globalNumTriangles*3 );
			}

			triCountArray[0] = triCount;
			for( i=MASTER_PROC+1; i<pd->numProcs; i++ ){
				MPI_Recv( &(triCountArray[i]), 1, MPI_INT, i, DATA_TAG, MPI_COMM_WORLD, &st );					
			}

			for( i=0; i<globalNumTriangles; i++ ){
				pd->triangleIndices[i] = pd->triangleIndices[0] + i*3;
	
				if ( pd->attributes->BuildTriangleNeighbours ){
					pd->triangleNeighbours[i] = pd->triangleNeighbours[0] + i*3;

					pd->triangleNeighbours[i][0] = pd->numTriangles;
					pd->triangleNeighbours[i][1] = pd->numTriangles;
					pd->triangleNeighbours[i][2] = pd->numTriangles;
				}
			}

			for(j=0; j<triCountArray[0]; j++){
				memcpy( pd->triangleIndices[j], triIndices[j], sizeof(int)*3 );
			}
		
			stride = j;
			for( i=MASTER_PROC+1; i<pd->numProcs; i++ ){
                          unsigned int *temp = (unsigned*)malloc( sizeof(unsigned int)*triCountArray[i]*3 );
					MPI_Recv( temp, triCountArray[i]*3, MPI_INT, i, DATA_TAG, MPI_COMM_WORLD, &st );					
				for(j=0; j<triCountArray[i]*3; j+=3){
					memcpy( pd->triangleIndices[stride++], &(temp[j]), sizeof(int)*3 );
				}
				free( temp );
			}

			for( i=0; i<pd->numTriangles; i++ ){
				pd->triangleIndices[i][0] = pd->initialOrder[pd->triangleIndices[i][0]];
				pd->triangleIndices[i][1] = pd->initialOrder[pd->triangleIndices[i][1]];
				pd->triangleIndices[i][2] = pd->initialOrder[pd->triangleIndices[i][2]];
			}

			/*for( j=0; j<globalNumTriangles; j++ ){
				printf( "indices[%d] = [%d, %d, %d]\n", j, pd->triangleIndices[j][0], pd->triangleIndices[j][1], pd->triangleIndices[j][2] );
			}*/
			
			stride = 0;
			if ( pd->attributes->BuildTriangleNeighbours ){
				int **nodeToTriangle = NULL;
				int *nodeToTriangleCount = NULL;
				int *triangleNeighboursCount = NULL;
				int triangles[2];
				int counter = 0;
				int found = 0;
				int m = 0;

				nodeToTriangle = Memory_Alloc_Array_Unnamed( int*, pd->numInputSites );
				nodeToTriangle[0] = Memory_Alloc_Array_Unnamed( int, numNeighboursSum );

				nodeToTriangleCount = Memory_Alloc_Array_Unnamed( int, pd->numInputSites );
				memset( nodeToTriangleCount, 0, sizeof( int ) * pd->numInputSites );
				
				triangleNeighboursCount = Memory_Alloc_Array_Unnamed( int, pd->numTriangles );
				memset( triangleNeighboursCount, 0, sizeof( int ) * pd->numTriangles );

				for( i=0; i<pd->numSites; i++ ){
					nodeToTriangle[i] = nodeToTriangle[0] + stride;

					stride += pd->numNeighbours[i];
				}

				for( i=0; i<pd->numTriangles; i++ ){
					nodeToTriangle[pd->triangleIndices[i][0]][nodeToTriangleCount[pd->triangleIndices[i][0]]++] = i;
					nodeToTriangle[pd->triangleIndices[i][1]][nodeToTriangleCount[pd->triangleIndices[i][1]]++] = i;
					nodeToTriangle[pd->triangleIndices[i][2]][nodeToTriangleCount[pd->triangleIndices[i][2]]++] = i;
				}
				
				for( i=0; i<pd->numInputSites; i++ ){
                                  for( j=0; j<(int)(pd->numNeighbours[i]); j++ ){
						triangles[0] = -1;
						triangles[1] = -1;
						counter = 0;
						for( k=0; k<nodeToTriangleCount[i]; k++ ){
							for( l=0; l<3; l++ ){
								if( pd->neighbours[i][j] == pd->triangleIndices[nodeToTriangle[i][k]][l] ){
									triangles[counter++] = nodeToTriangle[i][k];
								}
							}
						}
						
						if( (triangles[0] > -1) ){
							
							found = 0;
							for( m=0; m<3; m++ ){
                                                          if( pd->triangleNeighbours[triangles[0]][m] == (unsigned)(triangles[1]) ) found = 1;
							}
							if( !found ){
								pd->triangleNeighbours[triangles[0]][triangleNeighboursCount[triangles[0]]++] = triangles[1];
							}
						}
						if( (triangles[1] > -1) ){
						
							found = 0;
							for( m=0; m<3; m++ ){
                                                          if( pd->triangleNeighbours[triangles[1]][m] == (unsigned)(triangles[0]) ) found = 1;
							}
							if( !found ){
								pd->triangleNeighbours[triangles[1]][triangleNeighboursCount[triangles[1]]++] = triangles[0];
							}
						}
					}
				}
				
				Memory_Free( nodeToTriangle[0] );
				Memory_Free( nodeToTriangle );
				Memory_Free( nodeToTriangleCount );
				Memory_Free( triangleNeighboursCount );
			}
		}
		else{
			MPI_Send( &triCount, 1, MPI_INT, MASTER_PROC, DATA_TAG, MPI_COMM_WORLD );
			MPI_Send( triIndices[0], 3*triCount, MPI_INT, MASTER_PROC, DATA_TAG, MPI_COMM_WORLD );
		}

		Memory_Free( triCountArray );
	}	
}

	/*--------------------------------------------------------------------------------------------------------------------------
	** Private Member functions
	*/

int ParallelDelaunay_TranslateLocalToGlobal( ParallelDelaunay *self, int id )
{
	if( id >= (self->localTriangulation->idOffset+self->numLocalSites) ){
	
		if( (self->leftProc != self->numProcs) && (id < (self->localTriangulation->idOffset+self->numLocalSites + self->numHaloSites[0])) ){
			
			id = self->mappingTable[0][(id -
						(self->localTriangulation->idOffset+
						 self->numLocalSites))];
		}
		else{

			id = self->mappingTable[1][(id -
						(self->localTriangulation->idOffset+
						 self->numLocalSites+
						 self->numHaloSites[0]))];
		}
		return id;
	}
	else{
		return id;
	}
}

#define EDGE_TAG 102
#define BREAK_TAG 103
typedef struct SitePacket_t{
	int id;
	float xyz[3];
}SitePacket;
void ParallelDelaunaySendEdge( QuadEdgeRef edge, int rank, MPI_Comm *comm, MPI_Request *req )
{
	SitePacket sp[2];
	assert(edge);

	memcpy( sp[0].xyz, ((Site*)ORG(edge))->coord, sizeof( CoordF ) );
	sp[0].id = ((Site*)ORG(edge))->id;
	
	memcpy( sp[1].xyz, ((Site*)DEST(edge))->coord, sizeof( CoordF ) );
	sp[1].id = ((Site*)DEST(edge))->id;
	MPI_Isend( sp, sizeof(SitePacket)*2, MPI_BYTE, rank, EDGE_TAG, *comm, req );
}

QuadEdgeRef ParallelDelaunayRecvEdge( ParallelDelaunay *pd, int rank, MPI_Comm *comm )
{
	QuadEdgeRef edge = 0;
	CoordF *c[2];
	Site *s[2];
	MPI_Status st;
	SitePacket sp[2];
	Delaunay *d = NULL;
	
	assert( pd );
	d = pd->localTriangulation;
	
	edge = MakeQuadEdge( d->qp );
	
	assert(edge);

	s[0] = MemoryPool_NewObject( Site, pd->sitePool );
	s[1] = MemoryPool_NewObject( Site, pd->sitePool );
	
	c[0] = MemoryPool_NewObject( CoordF, pd->coordPool );
	c[1] = MemoryPool_NewObject( CoordF, pd->coordPool );
	
	MPI_Recv( sp, sizeof(SitePacket)*2, MPI_BYTE, rank, EDGE_TAG, *comm, &st );

	memcpy( (c[0]), sp[0].xyz, sizeof( CoordF ) );
	memcpy( (c[1]), sp[1].xyz, sizeof( CoordF ) );

	s[0]->id = sp[0].id;
	s[0]->coord = (c[0]);
	s[1]->id = sp[1].id;
	s[1]->coord = (c[1]);
	
	ORG(edge)=(s[0]);
	DEST(edge)=(s[1]);

	return edge;
}

int ParallelDelaunayListCompareFunction( void *a, void *b )
{
	QuadEdgeRef e1, e2;

	e1 = (QuadEdgeRef)a;
	e2 = (QuadEdgeRef)b;

	if( e1 > e2 ){
		return 1;
	}
	else if( e1 < e2 ){
		return -1;
	}
	else{
		return 0;
	}
}

void ParallelDelaunayListDeleteFunction( void *a )
{

}

QuadEdgeRef ParallelDelaunayFindLowestQuadEdge( ParallelDelaunay *pd, MPI_Comm *comm, int rank )
{
	QuadEdgeRef ldi=0, rdi=0;
	MPI_Request r;
	MPI_Status s;
	int globalBreak, localBreak;
	LinkedList *list = NULL;
	LinkedListIterator *iter;
	QuadEdgeRef result = 0;
	Delaunay *d = NULL;

	assert( pd );
	d = pd->localTriangulation;
	
	list = LinkedList_New( ParallelDelaunayListCompareFunction, NULL, NULL, ParallelDelaunayListDeleteFunction, LINKEDLIST_UNSORTED );
	iter = LinkedListIterator_New( list );
	
	if( rank == pd->leftProc ){
		rdi = d->leftMost;
	}
	else if( rank == pd->rightProc ){
		ldi = d->rightMost;
	}
	else{
		fprintf( stderr, "Failed to find lowest edge on rank %d..!\n Aborting..!\n", rank );
	}

	if( rank == pd->numProcs ) return 0;
	
	localBreak = 0;
	globalBreak = 0;

	while (1) 
	{
		localBreak = 0;
		globalBreak = 0;

		if( rank == pd->leftProc ){
			ParallelDelaunaySendEdge( rdi, rank, comm, &r );
			ldi = ParallelDelaunayRecvEdge( pd, rank, comm );
			MPI_Wait( &r, &s );
			LinkedList_InsertNode( list, (void*)ldi, sizeof( QuadEdgeRef* ) );
			
			if (RightOf((Site*)ORG(ldi), rdi)){
				rdi = ONEXT(SYM(rdi));
				localBreak = 1;
			}
			else{
			}
		}
		
		if( rank == pd->rightProc ){
			ParallelDelaunaySendEdge( ldi, rank, comm, &r );
			rdi = ParallelDelaunayRecvEdge( pd, rank, comm );
			MPI_Wait( &r, &s );
			LinkedList_InsertNode( list, (void*)rdi, sizeof( QuadEdgeRef* ) );
			
			if (LeftOf((Site*)ORG(rdi), ldi)){
				ldi = LNEXT(ldi);
				localBreak = 1;
			}
			else{
			}
		}

		MPI_Isend( &localBreak, 1, MPI_INT, rank, BREAK_TAG, *comm, &r );
		MPI_Recv( &(globalBreak), 1, MPI_INT, rank, BREAK_TAG, *comm, &s );
		MPI_Wait( &r, &s );
		globalBreak |= localBreak;
		
		if( globalBreak == 0 ){
			break;
		}
	}

	for( result=(QuadEdgeRef)((void*)LinkedListIterator_First(iter));
			result != 0;
			result=(QuadEdgeRef)((void*)LinkedListIterator_Next(iter)) ){

		DeleteQuadEdge( d->qp, result );
	}

	Stg_Class_Delete( list );
	Stg_Class_Delete( iter );
	
	if( rank == pd->leftProc ){
		return rdi;
	}
	else if( rank == pd->rightProc ){
		return ldi;
	}
	
	return 0;
}

void ParallelDelaunayMerge( ParallelDelaunay *pd, MPI_Comm *comm, int rank )
{
	QuadEdgeRef lowest = 0, lcand = 0, rcand = 0, basel = 0, baselPrev;
	MPI_Request r;
	MPI_Status s;
	int localBreak=0, globalBreak=0;
	double result = 0.0f;
	Delaunay *d = NULL;
	
	assert( pd );

	d = pd->localTriangulation;

	lowest = ParallelDelaunayFindLowestQuadEdge( pd, comm, rank );
	
	localBreak = 0;
	globalBreak = 0;

	if( rank == pd->numProcs ){
		fprintf( stderr, "Failed to merge rank %d with rank %d..!\nAborting..!\n", pd->rank, rank );
		assert( 0 );
	}
	
	if( rank == pd->leftProc ){
		rcand = lowest;
		
		ParallelDelaunaySendEdge( lowest, rank, comm, &r );
		lcand = ParallelDelaunayRecvEdge( pd, rank, comm );
		MPI_Wait( &r, &s );
		
		basel = MakeQuadEdge( d->qp );
		ORG(basel) = DEST(SYM(rcand));
		DEST(basel) = ORG(lcand);

		DeleteQuadEdge( d->qp, lcand );
		
		CCW( ((Site*)DEST(rcand))->coord, ((Site*)DEST(basel))->coord, ((Site*)ORG(basel))->coord, &result);
		if( result == 0.0f ){
			rcand = OPREV(rcand);
		}
		
		while(1){
			localBreak = 0;
			if (RightOf((Site*)DEST(rcand), basel)){
				while (InCircle((Site*)DEST(basel), (Site*)ORG(basel), (Site*)DEST(rcand), (Site*)DEST(OPREV(rcand)))){
					QuadEdgeRef t = OPREV(rcand);

					DeleteQuadEdge(d->qp, rcand);
					rcand = t;
				}
			}
			
			if (!RightOf((Site*)DEST(rcand), basel)) localBreak = 1;

			MPI_Isend( &localBreak, 1, MPI_INT, rank, BREAK_TAG, *comm, &r );
			MPI_Recv( &(globalBreak), 1, MPI_INT, rank, BREAK_TAG, *comm, &s );
			MPI_Wait( &r, &s );
			globalBreak &= localBreak;

			if( globalBreak ){
				break;
			}

			ParallelDelaunaySendEdge( rcand, rank, comm, &r );
			lcand = ParallelDelaunayRecvEdge( pd, rank, comm );
			MPI_Wait( &r, &s );

			if ( !RightOf((Site*)DEST(lcand), basel) ||
				( RightOf((Site*)DEST(rcand), basel) && 
				InCircle((Site*)DEST(lcand), (Site*)ORG(lcand), (Site*)ORG(rcand), (Site*)DEST(rcand)))){
			
				baselPrev = basel;
				basel = MakeQuadEdge(d->qp);
				ORG(basel) = DEST(rcand);
				DEST(basel) = ORG(SYM(baselPrev));
				
				rcand=LNEXT(rcand);
				DeleteQuadEdge( d->qp, lcand );
			}
			else{
				baselPrev = basel;
				basel = MakeQuadEdge(d->qp);
				ORG(basel) = DEST(SYM(baselPrev));
				DEST(basel) = ORG(SYM(lcand));
			}
		}
	}
	else if( rank == pd->rightProc ){
		lcand = lowest;
		
		ParallelDelaunaySendEdge( lowest, rank, comm, &r );
		rcand = ParallelDelaunayRecvEdge( pd, rank, comm );
		MPI_Wait( &r, &s );
		
		basel = MakeQuadEdge(d->qp);
		ORG(basel) = DEST(SYM(rcand));
		DEST(basel) = ORG(lcand);
		
		DeleteQuadEdge( d->qp, rcand );
		
		CCW( ((Site*)DEST(lcand))->coord, ((Site*)DEST(basel))->coord, ((Site*)ORG(basel))->coord, &result);
		if( result == 0.0f ){
			lcand = ONEXT(lcand);
		}
		
		while(1){
			localBreak = 0;
			if (RightOf((Site*)DEST(lcand), basel)){
				while (InCircle((Site*)DEST(basel), (Site*)ORG(basel), (Site*)DEST(lcand), (Site*)DEST(ONEXT(lcand)))){
					QuadEdgeRef t = ONEXT(lcand);
					
					DeleteQuadEdge(d->qp, lcand);
					lcand = t;
				}
			}
			
			if (!RightOf((Site*)DEST(lcand), basel)) localBreak = 1;

			MPI_Isend( &localBreak, 1, MPI_INT, rank, BREAK_TAG, *comm, &r );
			MPI_Recv( &(globalBreak), 1, MPI_INT, rank, BREAK_TAG, *comm, &s );
			MPI_Wait( &r, &s );
			globalBreak &= localBreak;

			if( globalBreak ){
				break;
			}
			
			ParallelDelaunaySendEdge( lcand, rank, comm, &r );
			rcand = ParallelDelaunayRecvEdge( pd, rank, comm );
			MPI_Wait( &r, &s );

			if ( !RightOf((Site*)DEST(lcand), basel) ||
				( RightOf((Site*)DEST(rcand), basel) && 
				InCircle((Site*)DEST(lcand), (Site*)ORG(lcand), (Site*)ORG(rcand), (Site*)DEST(rcand)))){
				
				baselPrev = basel;
				basel = MakeQuadEdge(d->qp);
				ORG(basel) = DEST(rcand);
				DEST(basel) = ORG(SYM(baselPrev));
			}
			else{
				baselPrev = basel;
				basel = MakeQuadEdge(d->qp);
				ORG(basel) = DEST(SYM(baselPrev));
				DEST(basel) = ORG(SYM(lcand));
				
				lcand = RPREV(lcand);
				DeleteQuadEdge( d->qp, rcand );
			}
		}
	}
}



