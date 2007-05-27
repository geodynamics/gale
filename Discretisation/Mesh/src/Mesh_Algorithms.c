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
** $Id: Mesh_Algorithms.c 3584 2006-05-16 11:11:07Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <mpi.h>

#include "Base/Base.h"
#include "Discretisation/Geometry/Geometry.h"

#include "Mesh.h"


/* Textual name of this class */
const Type Mesh_Algorithms_Type = "Mesh_Algorithms";

/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/

Mesh_Algorithms* Mesh_Algorithms_New( Name name ) {
	return _Mesh_Algorithms_New( sizeof(Mesh_Algorithms), 
				     Mesh_Algorithms_Type, 
				     _Mesh_Algorithms_Delete, 
				     _Mesh_Algorithms_Print, 
				     NULL, 
				     (void* (*)(Name))_Mesh_Algorithms_New, 
				     _Mesh_Algorithms_Construct, 
				     _Mesh_Algorithms_Build, 
				     _Mesh_Algorithms_Initialise, 
				     _Mesh_Algorithms_Execute, 
				     _Mesh_Algorithms_Destroy, 
				     name, 
				     NON_GLOBAL, 
				     _Mesh_Algorithms_SetMesh, 
				     _Mesh_Algorithms_Update, 
				     _Mesh_Algorithms_NearestVertex, 
				     _Mesh_Algorithms_Search, 
				     _Mesh_Algorithms_SearchElements, 
				     _Mesh_Algorithms_GetMinimumSeparation, 
				     _Mesh_Algorithms_GetLocalCoordRange, 
				     _Mesh_Algorithms_GetDomainCoordRange, 
				     _Mesh_Algorithms_GetGlobalCoordRange );
}

Mesh_Algorithms* _Mesh_Algorithms_New( MESH_ALGORITHMS_DEFARGS ) {
	Mesh_Algorithms* self;
	
	/* Allocate memory */
	assert( sizeOfSelf >= sizeof(Mesh_Algorithms) );
	self = (Mesh_Algorithms*)_Stg_Component_New( STG_COMPONENT_PASSARGS );

	/* Virtual info */
	self->setMeshFunc = setMeshFunc;
	self->updateFunc = updateFunc;
	self->nearestVertexFunc = nearestVertexFunc;
	self->searchFunc = searchFunc;
	self->searchElementsFunc = searchElementsFunc;
	self->getMinimumSeparationFunc = getMinimumSeparationFunc;
	self->getLocalCoordRangeFunc = getLocalCoordRangeFunc;
	self->getDomainCoordRangeFunc = getDomainCoordRangeFunc;
	self->getGlobalCoordRangeFunc = getGlobalCoordRangeFunc;

	/* Mesh_Algorithms info */
	_Mesh_Algorithms_Init( self );

	return self;
}

void _Mesh_Algorithms_Init( Mesh_Algorithms* self ) {
	self->nearestVertex = NULL;
	self->search = NULL;
	self->mesh = NULL;
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _Mesh_Algorithms_Delete( void* algorithms ) {
	Mesh_Algorithms*	self = (Mesh_Algorithms*)algorithms;

	/* Delete the parent. */
	_Stg_Component_Delete( self );
}

void _Mesh_Algorithms_Print( void* algorithms, Stream* stream ) {
	Mesh_Algorithms*	self = (Mesh_Algorithms*)algorithms;
	
	/* Set the Journal for printing informations */
	Stream* algorithmsStream;
	algorithmsStream = Journal_Register( InfoStream_Type, "Mesh_AlgorithmsStream" );

	/* Print parent */
	Journal_Printf( stream, "Mesh_Algorithms (ptr): (%p)\n", self );
	_Stg_Component_Print( self, stream );
}

void _Mesh_Algorithms_Construct( void* algorithms, Stg_ComponentFactory* cf, void* data ) {
}

void _Mesh_Algorithms_Build( void* algorithms, void* data ) {
}

void _Mesh_Algorithms_Initialise( void* algorithms, void* data ) {
}

void _Mesh_Algorithms_Execute( void* algorithms, void* data ) {
}

void _Mesh_Algorithms_Destroy( void* algorithms, void* data ) {
}

void _Mesh_Algorithms_SetMesh( void* algorithms, void* mesh ) {
	Mesh_Algorithms*	self = (Mesh_Algorithms*)algorithms;

	assert( self );
	assert( mesh );

	self->mesh = (Mesh*)mesh;
}

void _Mesh_Algorithms_Update( void* algorithms ) {
	Mesh_Algorithms*	self = (Mesh_Algorithms*)algorithms;
	unsigned		nDims;
	unsigned		d_i;

	assert( self );

	if( !self->mesh )
		return;

	if( Mesh_HasIncidence( self->mesh, MT_VERTEX, MT_VERTEX ) )
		self->nearestVertex = Mesh_Algorithms_NearestVertexWithNeighbours;
	else
		self->nearestVertex = Mesh_Algorithms_NearestVertexGeneral;

	nDims = Mesh_GetDimSize( self->mesh );
	for( d_i = 0; d_i < nDims; d_i++ ) {
		if( !Mesh_GetGlobalSize( self->mesh, d_i ) || !Mesh_HasIncidence( self->mesh, nDims, d_i ) )
			break;
	}
	if( d_i == nDims )
		self->search = Mesh_Algorithms_SearchWithFullIncidence;
	else if( Mesh_HasIncidence( self->mesh, MT_VERTEX, Mesh_GetDimSize( self->mesh ) ) )
		self->search = Mesh_Algorithms_SearchWithMinIncidence;
	else
		self->search = Mesh_Algorithms_SearchGeneral;
}

unsigned _Mesh_Algorithms_NearestVertex( void* algorithms, double* point ) {
	Mesh_Algorithms*	self = (Mesh_Algorithms*)algorithms;

	assert( self );
	assert( self->nearestVertex );

	return self->nearestVertex( self, point );
}

Bool _Mesh_Algorithms_Search( void* algorithms, double* point, 
			      MeshTopology_Dim* dim, unsigned* ind )
{
	Mesh_Algorithms*	self = (Mesh_Algorithms*)algorithms;

	assert( self );
	assert( self->search );

	return self->search( self, point, dim, ind );
}

Bool _Mesh_Algorithms_SearchElements( void* algorithms, double* point, 
				      unsigned* elInd )
{
	Mesh_Algorithms*	self = (Mesh_Algorithms*)algorithms;
	Mesh*			mesh;
	unsigned		dim, ind;

	assert( self );
	assert( self->mesh );
	assert( elInd );

	mesh = self->mesh;
	if( Mesh_Algorithms_Search( self, point, &dim, &ind ) ) {
		unsigned	nDims;

		nDims = Mesh_GetDimSize( mesh );
		if( dim != nDims ) {
			unsigned	nInc, *inc;
			unsigned	nLocalEls;
			unsigned	lowest;
			unsigned	global;
			unsigned	inc_i;

			/* Must have required incidence for this to work. */
			assert( Mesh_HasIncidence( mesh, dim, nDims ) );

			nLocalEls = Mesh_GetLocalSize( mesh, nDims );
			Mesh_GetIncidence( mesh, dim, ind, nDims, &nInc, &inc );
			assert( nInc );
			lowest = Mesh_DomainToGlobal( mesh, nDims, inc[0] );
			for( inc_i = 1; inc_i < nInc; inc_i++ ) {
				global = Mesh_DomainToGlobal( mesh, nDims, inc[inc_i] );
				if( global < lowest )
					lowest = global;
			}

			insist( Mesh_GlobalToDomain( mesh, nDims, lowest, elInd), == True );
		}
		else
			*elInd = ind;

		return True;
	}

	return False;
}

double _Mesh_Algorithms_GetMinimumSeparation( void* algorithms, double* perDim ) {
	Mesh_Algorithms*	self = (Mesh_Algorithms*)algorithms;
	Mesh*			mesh;
	unsigned		nDomainEls;
	double			minSep;
	double*			dimSep;
	unsigned		e_i;

	assert( self );
	assert( self->mesh );

	mesh = self->mesh;
	if( perDim )
		dimSep = Memory_Alloc_Array_Unnamed( double, Mesh_GetDimSize( mesh ) );
	else
		dimSep = NULL;

	minSep = HUGE_VAL;
	nDomainEls = Mesh_GetDomainSize( mesh, Mesh_GetDimSize( mesh ) );
	for( e_i = 0; e_i < nDomainEls; e_i++ ) {
		Mesh_ElementType*	elType;
		double			curSep;

		elType = Mesh_GetElementType( mesh, e_i );
		curSep = Mesh_ElementType_GetMinimumSeparation( elType, e_i, dimSep );
		if( curSep < minSep ) {
			minSep = curSep;
			if( perDim )
				memcpy( perDim, dimSep, Mesh_GetDimSize( mesh ) * sizeof(double) );
		}
	}

	return minSep;
}

void _Mesh_Algorithms_GetLocalCoordRange( void* algorithms, double* min, double* max ) {
	Mesh_Algorithms*	self = (Mesh_Algorithms*)algorithms;
	Mesh*			mesh;
	unsigned		nVerts, nEls;
	unsigned*		verts;
	double*			vert;
	unsigned		nDims;
	unsigned		v_i, e_i, d_i;

	assert( self );
	assert( self->mesh );
	assert( min );
	assert( max );

	mesh = self->mesh;
	nDims = Mesh_GetDimSize( mesh );
	nEls = Mesh_GetLocalSize( mesh, nDims );
	memcpy( min, Mesh_GetVertex( mesh, 0 ), nDims * sizeof(double) );
	memcpy( max, Mesh_GetVertex( mesh, 0 ), nDims * sizeof(double) );
	for( e_i = 0; e_i < nEls; e_i++ ) {
		Mesh_GetIncidence( mesh, nDims, e_i, 0, &nVerts, &verts );
		for( v_i = 0; v_i < nVerts; v_i++ ) {
			vert = Mesh_GetVertex( mesh, verts[v_i] );
			for( d_i = 0; d_i < nDims; d_i++ ) {
				if( vert[d_i] < min[d_i] )
					min[d_i] = vert[d_i];
				if( vert[d_i] > max[d_i] )
					max[d_i] = vert[d_i];
			}
		}
	}
}

void _Mesh_Algorithms_GetDomainCoordRange( void* algorithms, double* min, double* max ) {
	Mesh_Algorithms*	self = (Mesh_Algorithms*)algorithms;
	Mesh*			mesh;
	unsigned		nVerts;
	double*			vert;
	unsigned		nDims;
	unsigned		v_i, d_i;

	assert( self );
	assert( self->mesh );
	assert( min );
	assert( max );

	mesh = self->mesh;
	nDims = Mesh_GetDimSize( mesh );
	nVerts = Mesh_GetDomainSize( mesh, MT_VERTEX );
	memcpy( min, Mesh_GetVertex( mesh, 0 ), nDims * sizeof(double) );
	memcpy( max, Mesh_GetVertex( mesh, 0 ), nDims * sizeof(double) );
	for( v_i = 1; v_i < nVerts; v_i++ ) {
		vert = Mesh_GetVertex( mesh, v_i );
		for( d_i = 0; d_i < nDims; d_i++ ) {
			if( vert[d_i] < min[d_i] )
				min[d_i] = vert[d_i];
			if( vert[d_i] > max[d_i] )
				max[d_i] = vert[d_i];
		}
	}
}

void _Mesh_Algorithms_GetGlobalCoordRange( void* algorithms, double* min, double* max ) {
	Mesh_Algorithms*	self = (Mesh_Algorithms*)algorithms;
	Mesh*			mesh;
	unsigned		nDims;
	double			*localMin, *localMax;
	MPI_Comm		comm;
	unsigned		d_i;

	assert( self );
	assert( self->mesh );
	assert( min );
	assert( max );

	mesh = self->mesh;
	nDims = Mesh_GetDimSize( mesh );
	localMin = Memory_Alloc_Array_Unnamed( double, nDims );
	localMax = Memory_Alloc_Array_Unnamed( double, nDims );

	comm = Comm_GetMPIComm( Mesh_GetCommTopology( mesh, MT_VERTEX ) );
	Mesh_Algorithms_GetLocalCoordRange( self, localMin, localMax );
	for( d_i = 0; d_i < Mesh_GetDimSize( mesh ); d_i++ ) {
		MPI_Allreduce( localMin + d_i, min + d_i, 1, MPI_DOUBLE, MPI_MIN, comm );
		MPI_Allreduce( localMax + d_i, max + d_i, 1, MPI_DOUBLE, MPI_MAX, comm );
	}

	FreeArray( localMin );
	FreeArray( localMax );
}


/*--------------------------------------------------------------------------------------------------------------------------
** Public Functions
*/

#define Vec_Sep( nDims, v0, v1 )									\
	(((v0)[0] - (v1)[0]) * ((v0)[0] - (v1)[0]) +							\
	 (((nDims) >= 2) ? (((v0)[1] - (v1)[1]) * ((v0)[1] - (v1)[1]) + 				\
			    (((nDims) == 3) ? (((v0)[2] - (v1)[2]) * ((v0)[2] - (v1)[2])) : 0)) : 0))

unsigned Mesh_Algorithms_NearestVertexWithNeighbours( void* algorithms, double* point ) {
	Mesh_Algorithms*	self = (Mesh_Algorithms*)algorithms;
	Mesh*			mesh;
	unsigned		nDims;
	unsigned		curVert;
	double*			vert;
	double			sep;
	Bool			done;
	unsigned		nNbrs;
	unsigned*		nbrs;
	double			nbrSep;
	unsigned		nbr_i;

	assert( self );
	assert( self->mesh );
	assert( Mesh_HasIncidence( self->mesh, MT_VERTEX, MT_VERTEX ) );

	/* Get dimensionality. */
	mesh = self->mesh;
	nDims = Mesh_GetDimSize( mesh );

	/* Begin somewhere in the middle. */
	curVert = Mesh_GetDomainSize( mesh, MT_VERTEX ) / 2;

	/* Calc distance squared to current node. */
	vert = Mesh_GetVertex( mesh, curVert );
	sep = Vec_Sep( nDims, vert, point );

	/* Loop until we've found closest local node. */
	do {
		/* Get neighbouring vertices. */
		Mesh_GetIncidence( mesh, MT_VERTEX, curVert, MT_VERTEX, &nNbrs, &nbrs );

		/* Assume we'll be done after this loop. */
		done = True;

		/* Compare to neighbours. */
		for( nbr_i = 0; nbr_i < nNbrs; nbr_i++ ) {
			/* Calculate neighbour separation. */
			nbrSep = Vec_Sep( nDims, Mesh_GetVertex( mesh, nbrs[nbr_i] ), point );

			/* Closer? */
			if( nbrSep < sep ) {
				curVert = nbrs[nbr_i];
				sep = nbrSep;
				done = False;
			}
		}
	}
	while( !done );

	return curVert;
}

unsigned Mesh_Algorithms_NearestVertexGeneral( void* algorithms, double* point ) {
	Mesh_Algorithms*	self = (Mesh_Algorithms*)algorithms;
	Mesh*			mesh;
	unsigned		nDims;
	unsigned		nDomainVerts;
	double*			vert;
	unsigned		minVertInd;
	double			curSep, minSep;
	unsigned		v_i;

	assert( self );
	assert( self->mesh );
	assert( Mesh_GetDomainSize( self->mesh, MT_VERTEX ) );

	/* TODO: This is going to be hella slow, need to use some kind of spatial partitioning scheme. */

	mesh = self->mesh;
	nDims = Mesh_GetDimSize( mesh );
	nDomainVerts = Mesh_GetDomainSize( mesh, MT_VERTEX );

	vert = Mesh_GetVertex( mesh, 0 );
	minSep = Vec_Sep( nDims, vert, point );
	minVertInd = 0;

	for( v_i = 1; v_i < nDomainVerts; v_i++ ) {
		vert = Mesh_GetVertex( mesh, v_i );
		curSep = Vec_Sep( nDims, vert, point );
		if( curSep < minSep ) {
			minSep = curSep;
			minVertInd = v_i;
		}
	}

	return minVertInd;
}

Bool Mesh_Algorithms_SearchWithFullIncidence( void* algorithms, double* point, 
					      MeshTopology_Dim* dim, unsigned* ind )
{
	Mesh_Algorithms*	self = (Mesh_Algorithms*)algorithms;
	Mesh*			mesh;
	double			maxCrd[3], minCrd[3];
	unsigned		nDims;
	unsigned		nEls;
	unsigned		nearVert;
	unsigned		nInc, *inc;
	unsigned		e_i, d_i, inc_i;

	assert( self );
	assert( self->mesh );
	assert( Mesh_HasIncidence( self->mesh, MT_VERTEX, Mesh_GetDimSize( self->mesh ) ) );
	assert( dim );
	assert( ind );

	/* Get dimensionality. */
	mesh = self->mesh;
	nDims = Mesh_GetDimSize( mesh );

	/* If outside local range, immediately return false. */
	Mesh_GetDomainCoordRange( mesh, minCrd, maxCrd );
	for( d_i = 0; d_i < nDims; d_i++ ) {
		if( point[d_i] < minCrd[d_i] || point[d_i] > maxCrd[d_i] )
			return False;
	}

	/* Start by locating the closest vertex. */
	nearVert = Mesh_NearestVertex( mesh, point );

	/* Get vertex/element incidence. */
	Mesh_GetIncidence( mesh, MT_VERTEX, nearVert, nDims, 
			   &nInc, &inc );

	/* Search each of these incident elements in turn. */
	for( inc_i = 0; inc_i < nInc; inc_i++ ) {
		if( Mesh_ElementHasPoint( mesh, inc[inc_i], point, dim, ind ) )
			return True;
	}

	/* Brute force, search every element in turn (last resort). */
	nEls = Mesh_GetDomainSize( mesh, nDims );
	for( e_i = 0; e_i < nEls; e_i++ ) {
		if( Mesh_ElementHasPoint( mesh, e_i, point, dim, ind ) )
			return True;
	}

	return False;
}

Bool Mesh_Algorithms_SearchWithMinIncidence( void* algorithms, double* point, 
					      MeshTopology_Dim* dim, unsigned* ind )
{
	Mesh_Algorithms*	self = (Mesh_Algorithms*)algorithms;
	Mesh*			mesh;
	double			maxCrd[3], minCrd[3];
	unsigned		lowest;
	unsigned		nDims;
	unsigned		nEls;
	unsigned		nearVert;
	unsigned		nInc, *inc;
	unsigned		e_i, d_i, inc_i;

	assert( self );
	assert( self->mesh );
	assert( Mesh_HasIncidence( self->mesh, MT_VERTEX, Mesh_GetDimSize( self->mesh ) ) );
	assert( dim );
	assert( ind );

	/* Get dimensionality. */
	mesh = self->mesh;
	nDims = Mesh_GetDimSize( mesh );

	/* If outside local range, immediately return false. */
	Mesh_GetDomainCoordRange( mesh, minCrd, maxCrd );
	for( d_i = 0; d_i < nDims; d_i++ ) {
		if( point[d_i] < minCrd[d_i] || point[d_i] > maxCrd[d_i] )
			return False;
	}

	/* Start by locating the closest vertex. */
	nearVert = Mesh_NearestVertex( mesh, point );

	/* Get vertex/element incidence. */
	Mesh_GetIncidence( mesh, MT_VERTEX, nearVert, nDims, 
			   &nInc, &inc );

	/* Search all of these elements and return the element with lowest global index. */
	lowest = (unsigned)-1;
	for( inc_i = 0; inc_i < nInc; inc_i++ ) {
		if( Mesh_ElementHasPoint( mesh, inc[inc_i], point, dim, ind ) ) {
			unsigned	global;

			global = Mesh_DomainToGlobal( mesh, nDims, inc[inc_i] );
			if( global < lowest )
				lowest = global;
		}
	}
	if( lowest != (unsigned)-1 ) {
		insist( Mesh_GlobalToDomain( mesh, nDims, lowest, ind ), == True );
		*dim = nDims;
		return True;
	}

	/* Brute force, search every element in turn (last resort). */
	lowest = (unsigned)-1;
	nEls = Mesh_GetDomainSize( mesh, nDims );
	for( e_i = 0; e_i < nEls; e_i++ ) {
		if( Mesh_ElementHasPoint( mesh, e_i, point, dim, ind ) ) {
			unsigned	global;

			global = Mesh_DomainToGlobal( mesh, nDims, e_i );
			if( global < lowest )
				lowest = global;
		}
	}
	if( lowest != (unsigned)-1 ) {
		insist( Mesh_GlobalToDomain( mesh, nDims, lowest, ind ), == True );
		*dim = nDims;
		return True;
	}

	return False;
}

Bool Mesh_Algorithms_SearchGeneral( void* algorithms, double* point, 
				    MeshTopology_Dim* dim, unsigned* ind )
{
	Mesh_Algorithms*	self = (Mesh_Algorithms*)algorithms;
	Mesh*			mesh;
	double			maxCrd[3], minCrd[3];
	unsigned		nDims;
	unsigned		nEls;
	unsigned		e_i, d_i;

	assert( self );
	assert( self->mesh );
	assert( dim );
	assert( ind );

	/* Get dimensionality. */
	mesh = self->mesh;
	nDims = Mesh_GetDimSize( mesh );

	/* If outside local range, immediately return false. */
	Mesh_GetDomainCoordRange( mesh, minCrd, maxCrd );
	for( d_i = 0; d_i < nDims; d_i++ ) {
		if( point[d_i] < minCrd[d_i] || point[d_i] > maxCrd[d_i] )
			return False;
	}

	/* Brute force, search every element in turn. */
	nEls = Mesh_GetDomainSize( mesh, nDims );
	for( e_i = 0; e_i < nEls; e_i++ ) {
		if( Mesh_ElementHasPoint( mesh, e_i, point, dim, ind ) )
			return True;
	}

	return False;
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Private Functions
*/
