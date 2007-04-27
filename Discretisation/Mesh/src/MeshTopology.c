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
** $Id: MeshTopology.c 3584 2006-05-16 11:11:07Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <mpi.h>
#include "Base/Base.h"

#include "Discretisation/Geometry/Geometry.h"
#include "Discretisation/Shape/Shape.h"
#include "Mesh.h"


/* Textual name of this class */
const Type MeshTopology_Type = "MeshTopology";


/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/

MeshTopology* MeshTopology_New( Name name ) {
	return _MeshTopology_New( sizeof(MeshTopology), 
				  MeshTopology_Type, 
				  _MeshTopology_Delete, 
				  _MeshTopology_Print, 
				  NULL, 
				  (void* (*)(Name))_MeshTopology_New, 
				  _MeshTopology_Construct, 
				  _MeshTopology_Build, 
				  _MeshTopology_Initialise, 
				  _MeshTopology_Execute, 
				  _MeshTopology_Destroy, 
				  name, 
				  NON_GLOBAL );
}

MeshTopology* _MeshTopology_New( MESHTOPOLOGY_DEFARGS ) {
	MeshTopology* self;
	
	/* Allocate memory */
	assert( sizeOfSelf >= sizeof(MeshTopology) );
	self = (MeshTopology*)_Stg_Component_New( STG_COMPONENT_PASSARGS );

	/* Virtual info */

	/* MeshTopology info */
	_MeshTopology_Init( self );

	return self;
}

void _MeshTopology_Init( MeshTopology* self ) {
	self->nDims = 0;
	self->nTDims = 0;
	self->domains = NULL;
	self->nBndVerts = 0;
	self->bndVerts = NULL;
	self->shadowDepth = 0;
	self->nIncEls = NULL;
	self->incEls = NULL;
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _MeshTopology_Delete( void* topo ) {
	MeshTopology*	self = (MeshTopology*)topo;

	MeshTopology_Destruct( self );

	/* Delete the parent. */
	_Stg_Component_Delete( self );
}

void _MeshTopology_Print( void* topo, Stream* stream ) {
	MeshTopology*	self = (MeshTopology*)topo;
	
	/* Set the Journal for printing informations */
	Stream* topoStream;
	topoStream = Journal_Register( InfoStream_Type, "MeshTopologyStream" );

	/* Print parent */
	Journal_Printf( stream, "MeshTopology (ptr): (%p)\n", self );
	_Stg_Component_Print( self, stream );
}

void _MeshTopology_Construct( void* topo, Stg_ComponentFactory* cf, void* data ) {
}

void _MeshTopology_Build( void* topo, void* data ) {
}

void _MeshTopology_Initialise( void* topo, void* data ) {
}

void _MeshTopology_Execute( void* topo, void* data ) {
}

void _MeshTopology_Destroy( void* topo, void* data ) {
}

/*--------------------------------------------------------------------------------------------------------------------------
** Public Functions
*/

void MeshTopology_SetDimSize( void* topo, unsigned nDims ) {
	MeshTopology*	self = (MeshTopology*)topo;

	assert( self );
	assert( nDims > 0 );

	/* If we're changing dimensions, kill everything and begin again. */
	MeshTopology_Destruct( self );

	/* Set and allocate. */
	self->nDims = nDims;
	self->nTDims = nDims + 1;
	if( nDims ) {
		unsigned	 d_i;

		self->domains = AllocNamedArray( Decomp_Sync*, self->nTDims, "MeshTopology::domains" );
		self->nIncEls = AllocNamedArray2D( unsigned*, self->nTDims, self->nTDims, "MeshTopology::nIncEls" );
		self->incEls = AllocNamedArray2D( unsigned**, self->nTDims, self->nTDims, "MeshTopology::incEls" );
		memset( self->domains, 0, self->nTDims * sizeof(Decomp_Sync*) );
		for( d_i = 0; d_i < self->nTDims; d_i++ ) {
			memset( self->nIncEls[d_i], 0, self->nTDims * sizeof(unsigned*) );
			memset( self->incEls[d_i], 0, self->nTDims * sizeof(unsigned*) );
		}
	}
}

void MeshTopology_SetElements( void* topo, MeshTopology_Dim dim, unsigned nElements, unsigned* elements ) {
	MeshTopology*	self = (MeshTopology*)topo;
	Decomposer*	decomposer;

	assert( self );
	assert( dim < self->nTDims );
	assert( self->domains );
	assert( MeshTopology_ValidateElements( self, dim, nElements, elements ) );

	MeshTopology_ClearElements( self, dim );
	decomposer = Decomposer_New();
	Decomposer_Decompose( decomposer, nElements, elements, 
			      NULL, NULL, self->domains + dim );
	FreeObject( decomposer );

	/* If we just set the vertices, take the union of sink and source elements as the
	   boundary set. */
	if( dim == MT_VERTEX ) {
		Decomp_Sync*	sync;
		CommTopology*	commTopo;
		RangeSet*	bndSet;
		unsigned	nInds, *inds;
		unsigned	p_i;

		sync = MeshTopology_GetSync( self, MT_VERTEX );
		commTopo = MeshTopology_GetCommTopology( self, MT_VERTEX );
		bndSet = RangeSet_New();
		for( p_i = 0; p_i < CommTopology_GetIncidenceSize( commTopo ); p_i++ ) {
			Decomp_Sync_GetSources( sync, p_i, &nInds, &inds );
			RangeSet_AddIndices( bndSet, nInds, inds );
			Decomp_Sync_GetSinks( sync, p_i, &nInds, &inds );
			RangeSet_AddIndices( bndSet, nInds, inds );
		}
		self->bndVerts = NULL;
		RangeSet_GetIndices( bndSet, &self->nBndVerts, &self->bndVerts );
		FreeObject( bndSet );
	}
}

void MeshTopology_SetIncidence( void* topo, MeshTopology_Dim fromDim, MeshTopology_Dim toDim, 
				unsigned* nIncEls, unsigned** incEls )
{
	MeshTopology*	self = (MeshTopology*)topo;
	unsigned	size;
	unsigned	e_i;

	assert( self );
	assert( fromDim < self->nTDims );
	assert( toDim < self->nTDims );
	assert( !nIncEls || incEls );
	assert( Decomp_Sync_GetDomainSize( self->domains[fromDim] ) );
	assert( Decomp_Sync_GetDomainSize( self->domains[toDim] ) );

	/* Clear existing incidence. */
	KillArray( self->nIncEls[fromDim][toDim] );
	if( self->incEls[fromDim][toDim] ) {
		for( e_i = 0; e_i < Decomp_Sync_GetDomainSize( self->domains[fromDim] ); e_i++ )
			FreeArray( self->incEls[fromDim][toDim][e_i] );
		KillArray( self->incEls[fromDim][toDim] );
	}

	if( nIncEls && incEls ) {
		/* Allocate for incoming data. */
		size = Decomp_Sync_GetDomainSize( self->domains[fromDim] );
		self->nIncEls[fromDim][toDim] = AllocNamedArray( unsigned, size, "MeshTopology::nIncEls[][]" );
		self->incEls[fromDim][toDim] = AllocNamedArray( unsigned*, size, "MeshTopology::incEls[][]" );

		/* Copy the lot. */
		memcpy( self->nIncEls[fromDim][toDim], nIncEls, size * sizeof(unsigned) );
		for( e_i = 0; e_i < size; e_i++ ) {
			self->incEls[fromDim][toDim][e_i] = AllocArray( unsigned, nIncEls[e_i] );
			memcpy( self->incEls[fromDim][toDim][e_i], incEls[e_i], 
				self->nIncEls[fromDim][toDim][e_i] * sizeof(unsigned) );
		}
	}
}

void MeshTopology_SetShadowDepth( void* topo, unsigned depth ) {
	Stream*		errorStream = Journal_Register( ErrorStream_Type, "MeshTopology::SetShadowDepth" );
	MeshTopology*	self = (MeshTopology*)topo;
	CommTopology*	commTopo;

	assert( self );
	assert( depth );

	/* Store shadow depth. */
	self->shadowDepth = depth;

	/* Ensure all communication topologies are equivalent. */
	MeshTopology_CommUnion( self );

	/* If there are no neighbouring processors, forget about it. */
	commTopo = MeshTopology_GetCommTopology( self, MT_VERTEX );
	if( !CommTopology_GetIncidenceSize( commTopo ) )
		return;

	/* If we don't have sufficient incidence, we need to error. */
	Journal_Firewall( MeshTopology_HasIncidence( self, self->nDims, MT_VERTEX ) && 
			  MeshTopology_HasIncidence( self, self->nDims, MT_VERTEX ), 
			  errorStream, 
			  "\n" \
			  "********************************************************\n" \
			  "* Error: Topological shadows cannot be automatically   *\n" \
			  "*        built without a minimum of element-vertex and *\n" \
			  "*        vertex-element incidence.                     *\n" \
			  "********************************************************\n" \
			  "\n" );

	/* Construct and exchange shadow elements and incidence. */
	MeshTopology_BuildShadows( self );
}

void MeshTopology_SetSync( void* topo, MeshTopology_Dim dim, Decomp_Sync* sync ) {
	MeshTopology*	self = (MeshTopology*)topo;

	assert( self );
	assert( dim < self->nTDims );
	assert( self->domains );

	MeshTopology_ClearElements( self, dim );
	self->domains[dim] = sync;
	if( sync )
		Stg_Class_AddRef( sync );
}

void MeshTopology_Complete( void* topo ) {
	MeshTopology*	self = (MeshTopology*)topo;
	unsigned	d_i, d_j;

	assert( self );

	/* Complete downwards. */
	for( d_i = 0; d_i < self->nTDims - 1; d_i++ ) {
		for( d_j = d_i + 2; d_j < self->nTDims; d_j++ ) {
			if( !self->nIncEls[d_j][d_i] )
				MeshTopology_Cascade( self, d_j, d_i );
		}
	}

	/* Invert missing up relations. */
	for( d_i = 0; d_i < self->nTDims - 1; d_i++ ) {
		for( d_j = d_i + 1; d_j < self->nTDims; d_j++ ) {
			if( !self->nIncEls[d_i][d_j] )
				MeshTopology_Invert( self, d_i, d_j );
		}
	}

	/* Build neighbourhoods. */
	for( d_i = 0; d_i < self->nTDims; d_i++ ) {
		if( !self->nIncEls[d_i][d_i] )
			MeshTopology_Neighbourhood( self, d_i );
	}
}

void MeshTopology_Invert( void* topo, MeshTopology_Dim fromDim, MeshTopology_Dim toDim ) {
	Stream*		errorStream = Journal_Register( ErrorStream_Type, "MeshTopology::Invert" );
	MeshTopology*	self = (MeshTopology*)topo;
	unsigned	fromSize, toSize;
	unsigned*	invNIncEls;
	unsigned**	invIncEls;
	unsigned*	nIncEls;
	unsigned**	incEls;
	unsigned	e_i;

	/* Sanity check. */
	assert( self );
	assert( fromDim < self->nTDims );
	assert( toDim < self->nTDims );

	Journal_Firewall( self->domains[fromDim] && self->domains[toDim] && 
			  MeshTopology_HasIncidence( self, toDim, fromDim ), 
			  errorStream, 
			  "\n" \
			  "*******************************************************\n" \
			  "* Error: Cannot invert a topological relation without *\n" \
			  "*        pre-existing relation.                       *\n" \
			  "*******************************************************\n" \
			  "\n" );

	/* Shortcuts. */
	fromSize = Decomp_Sync_GetDomainSize( self->domains[fromDim] );
	toSize = Decomp_Sync_GetDomainSize( self->domains[toDim] );
	invNIncEls = self->nIncEls[toDim][fromDim];
	invIncEls = self->incEls[toDim][fromDim];

	/* Allocate some tables. */
	nIncEls = AllocNamedArray( unsigned, fromSize, "MeshTopology::nIncEls[][]" );
	memset( nIncEls, 0, fromSize * sizeof(unsigned) );

	/* Two phase process: count the numbers then allocate and do it again. */
	for( e_i = 0; e_i < toSize; e_i++ ) {
		unsigned	inc_i;

		for( inc_i = 0; inc_i < invNIncEls[e_i]; inc_i++ ) {
			unsigned	elInd = invIncEls[e_i][inc_i];

			nIncEls[elInd]++;
		}
	}

	/* Build up the tables. */
	incEls = AllocNamedComplex2D( unsigned, fromSize, nIncEls, "MeshTopology::incEls[][]" );
	memset( nIncEls, 0, fromSize * sizeof(unsigned) );
	for( e_i = 0; e_i < toSize; e_i++ ) {
		unsigned	inc_i;

		for( inc_i = 0; inc_i < invNIncEls[e_i]; inc_i++ ) {
			unsigned	elInd = invIncEls[e_i][inc_i];

			incEls[elInd][nIncEls[elInd]++] = e_i;
		}
	}

	/* Transfer to permanent storage. */
	MeshTopology_SetIncidence( self, fromDim, toDim, nIncEls, incEls );
	FreeArray( nIncEls );
	FreeArray( incEls );
}

void MeshTopology_Cascade( void* topo, MeshTopology_Dim fromDim, MeshTopology_Dim toDim ) {
	MeshTopology*	self = (MeshTopology*)topo;
	unsigned	maxInc = 0;
	unsigned*	nIncEls;
	unsigned**	incEls;
	unsigned	e_i;

	assert( self );
	assert( fromDim < self->nTDims );
	assert( toDim < self->nTDims );
	assert( fromDim > toDim + 1 );

	/* Determine maximum incidence for any element. */
	for( e_i = 0; e_i < self->domains[fromDim]->nDomains; e_i++ ) {
		unsigned	curInc = 0;
		unsigned	inc_i;

		for( inc_i = 0; inc_i < self->nIncEls[fromDim][fromDim - 1][e_i]; inc_i++ ) {
			unsigned	incEl = self->incEls[fromDim][fromDim - 1][e_i][inc_i];

			curInc += self->nIncEls[fromDim - 1][toDim][incEl];
		}

		if( curInc > maxInc )
			maxInc = curInc;
	}

	/* Allocate. */
	nIncEls = AllocArray( unsigned, self->domains[fromDim]->nDomains );
	incEls = AllocArray2D( unsigned, self->domains[fromDim]->nDomains, maxInc );
	memset( nIncEls, 0, self->domains[fromDim]->nDomains * sizeof(unsigned) );

	/* Determine actual incidence. */
	for( e_i = 0; e_i < self->domains[fromDim]->nDomains; e_i++ ) {
		unsigned	inc_i;

		for( inc_i = 0; inc_i < self->nIncEls[fromDim][fromDim - 1][e_i]; inc_i++ ) {
			unsigned	incEl = self->incEls[fromDim][fromDim - 1][e_i][inc_i];
			unsigned	inc_j;

			for( inc_j = 0; inc_j < self->nIncEls[fromDim - 1][toDim][incEl]; inc_j++ ) {
				unsigned	target = self->incEls[fromDim - 1][toDim][incEl][inc_j];
				unsigned	e_j;

				for( e_j = 0; e_j < nIncEls[e_i]; e_j++ ) {
					if( incEls[e_i][e_j] == target )
						break;
				}
				if( e_j == nIncEls[e_i] )
					incEls[e_i][nIncEls[e_i]++] = target;
			}
		}
	}

	/* Set incidence and kill temporary arrays. */
	MeshTopology_SetIncidence( self, fromDim, toDim, nIncEls, incEls );
	FreeArray( nIncEls );
	FreeArray( incEls );
}

void MeshTopology_Neighbourhood( void* topo, MeshTopology_Dim dim ) {
	MeshTopology*		self = (MeshTopology*)topo;
	unsigned		size;
	MeshTopology_Dim	toDim;
	unsigned*		nNbrs;
	unsigned**		nbrs;
	unsigned		e_i;

	assert( self );
	assert( dim < self->nTDims );
	assert( Decomp_Sync_GetDomainSize( self->domains[dim] ) );

	/* Vertex neighbours search upwards, all others search downwards. */
	size = Decomp_Sync_GetDomainSize( self->domains[dim] );
	toDim = (dim == MT_VERTEX) ? MT_EDGE : MT_VERTEX;

	/* Allocate some space for neighbour counts. */
	nNbrs = AllocArray( unsigned, size ); 

	/* Calculate maximum neighbours for each element. */
	for( e_i = 0; e_i < size; e_i++ ) {
		unsigned	nNodes = self->nIncEls[dim][toDim][e_i];
		unsigned*	nodes = self->incEls[dim][toDim][e_i];
		unsigned	n_i;

		nNbrs[e_i] = 0;
		for( n_i = 0; n_i < nNodes; n_i++ ) {
			unsigned	nCurNbrs = self->nIncEls[toDim][dim][nodes[n_i]];

			nNbrs[e_i] += nCurNbrs - 1;
		}
	}

	/* Allocate for maximum neighbours and clear neighbour counts. */
	nbrs = AllocComplex2D( unsigned, size, nNbrs );
	memset( nNbrs, 0, size * sizeof(unsigned) );

	/* Build neighbours for each element of dimension 'dim'. */
	for( e_i = 0; e_i < size; e_i++ ) {
		unsigned	nNodes = self->nIncEls[dim][toDim][e_i];
		unsigned*	nodes = self->incEls[dim][toDim][e_i];
		unsigned	n_i;

		/* Build a set of unique element nodes' elements. */
		for( n_i = 0; n_i < nNodes; n_i++ ) {
			unsigned	nCurNbrs = self->nIncEls[toDim][dim][nodes[n_i]];
			unsigned	e_j;

			for( e_j = 0; e_j < nCurNbrs; e_j++ ) {
				unsigned	curNbr = self->incEls[toDim][dim][nodes[n_i]][e_j];
				unsigned	nbr_i;

				/* Don't add current element. */
				if( curNbr == e_i )
					continue;

				/* Ensure is unique. */
				for( nbr_i = 0; nbr_i < nNbrs[e_i]; nbr_i++ ) {
					if( nbrs[e_i][nbr_i] == curNbr )
						break;
				}
				if( nbr_i == nNbrs[e_i] )
					nbrs[e_i][nNbrs[e_i]++] = curNbr;
			}
		}
	}

	/* Transfer. */
	MeshTopology_SetIncidence( self, dim, dim, nNbrs, nbrs );
	FreeArray( nNbrs );
	FreeArray( nbrs );
}

unsigned MeshTopology_GetGlobalSize( void* meshTopology, MeshTopology_Dim dim ) {
	MeshTopology*	self = (MeshTopology*)meshTopology;

	assert( self );
	assert( dim < self->nTDims );
	assert( self->domains );

	return self->domains[dim] ? Decomp_Sync_GetGlobalSize( self->domains[dim] ) : 0;
}

unsigned MeshTopology_GetLocalSize( void* meshTopology, MeshTopology_Dim dim ) {
	MeshTopology*	self = (MeshTopology*)meshTopology;

	assert( self );
	assert( dim < self->nTDims );
	assert( self->domains );

	return self->domains[dim] ? Decomp_Sync_GetLocalSize( self->domains[dim] ) : 0;
}

unsigned MeshTopology_GetRemoteSize( void* meshTopology, MeshTopology_Dim dim ) {
	MeshTopology*	self = (MeshTopology*)meshTopology;

	assert( self );
	assert( dim < self->nTDims );
	assert( self->domains );

	return self->domains[dim] ? Decomp_Sync_GetRemoteSize( self->domains[dim] ) : 0;
}

unsigned MeshTopology_GetDomainSize( void* meshTopology, MeshTopology_Dim dim ) {
	MeshTopology*	self = (MeshTopology*)meshTopology;

	assert( self );
	assert( dim < self->nTDims );
	assert( self->domains );

	return self->domains[dim] ? Decomp_Sync_GetDomainSize( self->domains[dim] ) : 0;
}

unsigned MeshTopology_GetSharedSize( void* meshTopology, MeshTopology_Dim dim ) {
	MeshTopology*	self = (MeshTopology*)meshTopology;

	assert( self );
	assert( dim < self->nTDims );
	assert( self->domains );

	return self->domains[dim] ? Decomp_Sync_GetSharedSize( self->domains[dim] ) : 0;
}

void MeshTopology_GetLocalElements( void* meshTopology, MeshTopology_Dim dim, unsigned* nEls, unsigned** els ) {
	MeshTopology*	self = (MeshTopology*)meshTopology;

	assert( self );
	assert( dim < self->nTDims );
	assert( self->domains );
	assert( self->domains[dim] );

	Decomp_GetLocals( self->domains[dim]->decomp, nEls, els );
}

void MeshTopology_GetRemoteElements( void* meshTopology, MeshTopology_Dim dim, unsigned* nEls, unsigned** els ) {
	MeshTopology*	self = (MeshTopology*)meshTopology;

	assert( self );
	assert( dim < self->nTDims );
	assert( self->domains );

	Decomp_Sync_GetRemotes( self->domains[dim], nEls, els );
}

CommTopology* MeshTopology_GetCommTopology( void* meshTopology, MeshTopology_Dim dim ) {
	MeshTopology*	self = (MeshTopology*)meshTopology;

	assert( self );
	assert( dim < self->nTDims );
	assert( self->domains );

	return self->domains[dim] ? Decomp_Sync_GetCommTopology( self->domains[dim] ) : NULL;
}

Decomp_Sync* MeshTopology_GetSync( void* meshTopology, MeshTopology_Dim dim ) {
	MeshTopology*	self = (MeshTopology*)meshTopology;

	assert( self );
	assert( dim < self->nTDims );
	assert( self->domains );

	return self->domains[dim];
}

Bool MeshTopology_GlobalToDomain( void* meshTopology, MeshTopology_Dim dim, unsigned global, unsigned* domain ) {
	MeshTopology*	self = (MeshTopology*)meshTopology;

	assert( self );
	assert( dim < self->nTDims );
	assert( self->domains );
	assert( self->domains[dim] );

	return Decomp_Sync_GlobalToDomain( self->domains[dim], global, domain );
}

unsigned MeshTopology_DomainToGlobal( void* meshTopology, MeshTopology_Dim dim, unsigned domain ) {
	MeshTopology*	self = (MeshTopology*)meshTopology;

	assert( self );
	assert( dim < self->nTDims );
	assert( self->domains );
	assert( self->domains[dim] );

	return Decomp_Sync_DomainToGlobal( self->domains[dim], domain );
}

Bool MeshTopology_DomainToShared( void* meshTopology, MeshTopology_Dim dim, unsigned domain, unsigned* shared ) {
	MeshTopology*	self = (MeshTopology*)meshTopology;

	assert( self );
	assert( dim < self->nTDims );
	assert( self->domains );
	assert( self->domains[dim] );

	return Decomp_Sync_DomainToShared( self->domains[dim], domain, shared );
}

unsigned MeshTopology_SharedToDomain( void* meshTopology, MeshTopology_Dim dim, unsigned shared ) {
	MeshTopology*	self = (MeshTopology*)meshTopology;

	assert( self );
	assert( dim < self->nTDims );
	assert( self->domains );
	assert( self->domains[dim] );

	return Decomp_Sync_SharedToDomain( self->domains[dim], shared );
}

unsigned MeshTopology_GetOwner( void* meshTopology, MeshTopology_Dim dim, unsigned remote ) {
	MeshTopology*	self = (MeshTopology*)meshTopology;

	assert( self );
	assert( dim < self->nTDims );
	assert( self->domains );

	return Decomp_Sync_GetOwner( self->domains[dim], remote );
}

void MeshTopology_GetSharers( void* meshTopology, MeshTopology_Dim dim, unsigned shared, 
			      unsigned* nSharers, unsigned** sharers )
{
	MeshTopology*	self = (MeshTopology*)meshTopology;

	assert( self );
	assert( dim < self->nTDims );
	assert( self->domains );
	assert( self->domains[dim] );

	Decomp_Sync_GetSharers( self->domains[dim], shared, nSharers, sharers );
}

Bool MeshTopology_HasIncidence( void* meshTopology, MeshTopology_Dim fromDim, MeshTopology_Dim toDim ) {
	MeshTopology*	self = (MeshTopology*)meshTopology;

	assert( self );
	assert( fromDim < self->nTDims );
	assert( toDim < self->nTDims );
	assert( self->nIncEls );
	assert( self->incEls );

	return (self->nIncEls[fromDim][toDim] && self->incEls[fromDim][toDim]);
}

unsigned MeshTopology_GetIncidenceSize( void* meshTopology, MeshTopology_Dim fromDim, unsigned fromInd, 
					MeshTopology_Dim toDim )
{
	MeshTopology*	self = (MeshTopology*)meshTopology;

	assert( self );
	assert( self->domains );
	assert( fromDim < self->nTDims );
	assert( toDim < self->nTDims );
	assert( fromInd < Decomp_Sync_GetDomainSize( self->domains[fromDim] ) );
	assert( self->nIncEls );
	assert( self->nIncEls[fromDim][toDim] );

	return self->nIncEls[fromDim][toDim][fromInd];
}

void MeshTopology_GetIncidence( void* meshTopology, MeshTopology_Dim fromDim, unsigned fromInd, MeshTopology_Dim toDim, 
				unsigned* nInc, unsigned** inc )
{
	MeshTopology*	self = (MeshTopology*)meshTopology;

	assert( self );
	assert( self->domains );
	assert( fromDim < self->nTDims );
	assert( toDim < self->nTDims );
	assert( fromInd < Decomp_Sync_GetDomainSize( self->domains[fromDim] ) );
	assert( self->nIncEls );
	assert( self->nIncEls[fromDim][toDim] );
	assert( self->incEls );
	assert( self->incEls[fromDim][toDim] );
	assert( self->incEls[fromDim][toDim][fromInd] );
	assert( nInc );
	assert( inc );

	*nInc = self->nIncEls[fromDim][toDim][fromInd];
	*inc = self->incEls[fromDim][toDim][fromInd];
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Private Functions
*/

void MeshTopology_ClearElements( MeshTopology* self, MeshTopology_Dim dim ) {
	assert( self );
	assert( dim < self->nTDims );
	assert( self->domains );

	if( self->domains[dim] ) {
		MeshTopology_ClearIncidence( self, dim );
		Stg_Class_RemoveRef( self->domains[dim] );
		self->domains[dim] = NULL;
	}

	if( dim == MT_VERTEX ) {
		KillArray( self->bndVerts );
		self->nBndVerts = 0;
	}
}

void MeshTopology_ClearIncidence( MeshTopology* self, MeshTopology_Dim dim ) {
	unsigned	nDomainEls;
	unsigned	d_i;

	assert( self );
	assert( dim < self->nTDims );

	nDomainEls = MeshTopology_GetDomainSize( self, dim );
	for( d_i = 0; d_i < self->nTDims; d_i++ ) {
		KillArray( self->nIncEls[dim][d_i] );
		KillArray2D( nDomainEls, self->incEls[dim][d_i] );
	}
}

void MeshTopology_CommUnion( MeshTopology* self ) {
	Decomp_Sync*	sync;
	CommTopology*	commTopo;
	RangeSet	*allSet, *curSet, *tmpSet;
	unsigned	nIncRanks, *incRanks;
	unsigned	d_i;

	assert( self );

	allSet = RangeSet_New();
	for( d_i = 0; d_i < self->nTDims; d_i++ ) {
		commTopo = MeshTopology_GetCommTopology( self, d_i );
		if( !commTopo )
			continue;
		CommTopology_GetIncidence( commTopo, &nIncRanks, &incRanks );
		RangeSet_AddIndices( allSet, nIncRanks, incRanks );
	}

	curSet = RangeSet_New();
	tmpSet = RangeSet_New();
	for( d_i = 0; d_i < self->nTDims; d_i++ ) {
		sync = MeshTopology_GetSync( self, d_i );
		if( !sync )
			continue;
		RangeSet_Clear( tmpSet );
		RangeSet_Union( tmpSet, allSet );
		commTopo = MeshTopology_GetCommTopology( self, d_i );
		CommTopology_GetIncidence( commTopo, &nIncRanks, &incRanks );
		RangeSet_SetIndices( curSet, nIncRanks, incRanks );
		RangeSet_Subtraction( tmpSet, curSet );
		incRanks = NULL;
		RangeSet_GetIndices( tmpSet, &nIncRanks, &incRanks );
		Decomp_Sync_AddRemoteRanks( sync, nIncRanks, incRanks );
	}

	FreeObject( allSet );
	FreeObject( curSet );
	FreeObject( tmpSet );
}

void MeshTopology_BuildShadows( MeshTopology* self ) {
	CommTopology*	commTopo;
	unsigned	nIncRanks;
	unsigned	**nShdEls, ***shdEls;
	unsigned	**nExtEls, ***extEls;
	unsigned**	nShdInc;
	Stg_Byte***	shdInc;
	unsigned*	nOldIncRanks;
	unsigned	d_i, p_i;

	assert( self );

	commTopo = MeshTopology_GetCommTopology( self, MT_VERTEX );
	nIncRanks = CommTopology_GetIncidenceSize( commTopo );

	nShdEls = AllocArray( unsigned*, self->nTDims );
	shdEls = AllocArray( unsigned**, self->nTDims );
	nExtEls = AllocArray2D( unsigned, self->nTDims, nIncRanks );
	extEls = AllocArray2D( unsigned*, self->nTDims, nIncRanks );

	for( d_i = self->nTDims - 1; d_i < self->nTDims; d_i-- ) {
		nShdEls[d_i] = AllocArray( unsigned, nIncRanks );
		shdEls[d_i] = AllocArray( unsigned*, nIncRanks );
		memset( nShdEls[d_i], 0, nIncRanks * sizeof(unsigned) );
		memset( shdEls[d_i], 0, nIncRanks * sizeof(unsigned*) );
		memset( nExtEls[d_i], 0, nIncRanks * sizeof(unsigned) );
		memset( extEls[d_i], 0, nIncRanks * sizeof(unsigned*) );

		MeshTopology_BuildShadowedElements( self, d_i, nShdEls, shdEls, 
						    nExtEls, extEls );
		MeshTopology_ExchangeShadowedElements( self, d_i, nShdEls[d_i], shdEls[d_i] );
		MeshTopology_ExchangeExternalElements( self, d_i, nExtEls[d_i], extEls[d_i], 
						       nShdEls + d_i, shdEls + d_i );
	}

	for( d_i = 0; d_i < self->nTDims; d_i++ ) {
		for( p_i = 0; p_i < nIncRanks; p_i++ )
			FreeArray( extEls[d_i][p_i] );
	}
	FreeArray( nExtEls ); FreeArray( extEls );

	nOldIncRanks = AllocArray( unsigned, self->nTDims );
	for( d_i = 0; d_i < self->nTDims; d_i++ ) {
		commTopo = MeshTopology_GetCommTopology( self, d_i );
		if( !commTopo )
			continue;

		nOldIncRanks[d_i] = CommTopology_GetIncidenceSize( commTopo );
	}
	MeshTopology_CommUnion( self );
	for( d_i = 0; d_i < self->nTDims; d_i++ ) {
		commTopo = MeshTopology_GetCommTopology( self, d_i );
		if( !commTopo )
			continue;

		nIncRanks = CommTopology_GetIncidenceSize( commTopo );
		if( nIncRanks != nOldIncRanks[d_i] ) {
			nShdEls[d_i] = ReallocArray( nShdEls[d_i], unsigned, nIncRanks );
			shdEls[d_i] = ReallocArray( shdEls[d_i], unsigned*, nIncRanks );
			memset( nShdEls[d_i] + nOldIncRanks[d_i], 0, 
				(nIncRanks - nOldIncRanks[d_i]) * sizeof(unsigned) );
			memset( shdEls[d_i] + nOldIncRanks[d_i], 0, 
				(nIncRanks - nOldIncRanks[d_i]) * sizeof(unsigned*) );
		}
	}
	FreeArray( nOldIncRanks );

	commTopo = MeshTopology_GetCommTopology( self, MT_VERTEX );
	nIncRanks = CommTopology_GetIncidenceSize( commTopo );

	nShdInc = AllocArray2D( unsigned, self->nTDims, nIncRanks );
	shdInc = AllocArray2D( Stg_Byte*, self->nTDims, nIncRanks );

	for( d_i = 0; d_i < self->nTDims; d_i++ ) {
		memset( nShdInc[d_i], 0, nIncRanks * sizeof(unsigned) );
		memset( shdInc[d_i], 0, nIncRanks * sizeof(unsigned*) );

		MeshTopology_BuildShadowedIncidence( self, d_i, nShdEls[d_i], shdEls[d_i], 
						     nShdInc[d_i], shdInc[d_i] );
		MeshTopology_ExchangeShadowedIncidence( self, d_i, nShdInc[d_i], shdInc[d_i] );
	}

	for( d_i = 0; d_i < self->nTDims; d_i++ ) {
		if( !MeshTopology_GetCommTopology( self, d_i ) )
			continue;

		for( p_i = 0; p_i < nIncRanks; p_i++ ) {
			FreeArray( shdEls[d_i][p_i] );
			FreeArray( shdInc[d_i][p_i] );
		}
	}
	FreeArray( nShdEls ); FreeArray2D( self->nTDims, shdEls );
	FreeArray( nShdInc ); FreeArray( shdInc );
}

void MeshTopology_BuildShadowedElements( MeshTopology* self, MeshTopology_Dim dim, 
					 unsigned** nShadowedEls, unsigned*** shadowedEls, 
					 unsigned** nExternalEls, unsigned*** externalEls )
{
	assert( self );

	if( dim == self->nDims ) {
		MeshTopology_BuildTopShadowedElements( self, nShadowedEls, shadowedEls, 
						       nExternalEls, externalEls );
	}
	else {
		MeshTopology_BuildMidShadowedElements( self, dim, nShadowedEls, shadowedEls, 
						       nExternalEls, externalEls );
	}
}

void MeshTopology_BuildTopShadowedElements( MeshTopology* self, unsigned** nShadowedEls, unsigned*** shadowedEls, 
					    unsigned** nExternalEls, unsigned*** externalEls )
{
	CommTopology*	commTopo;
	unsigned	nDomainEls;
	unsigned	nIncRanks;
	IndexSet**	shdSets;
	Decomp_Sync*	sync;
	unsigned	p_i, s_i;

	assert( self );
	assert( nShadowedEls && shadowedEls );
	assert( nExternalEls && externalEls );

	/* Extract communicator incidence. */
	commTopo = MeshTopology_GetCommTopology( self, self->nDims );
	nIncRanks = CommTopology_GetIncidenceSize( commTopo );

	/* Build some index sets. */
	sync = MeshTopology_GetSync( self, MT_VERTEX );
	nDomainEls = Decomp_Sync_GetDomainSize( sync );
	shdSets = AllocArray( IndexSet*, nIncRanks );
	for( p_i = 0; p_i < nIncRanks; p_i++ )
		shdSets[p_i] = IndexSet_New( nDomainEls );

	/* Store connected vertices. */
	for( s_i = 0; s_i < Decomp_Sync_GetSharedSize( sync ); s_i++ ) {
		unsigned	domainInd;
		unsigned	nSharers, *sharers;
		unsigned	sharer_i;

		domainInd = Decomp_Sync_SharedToDomain( sync, s_i );
		Decomp_Sync_GetSharers( sync, s_i, &nSharers, &sharers );
		for( sharer_i = 0; sharer_i < nSharers; sharer_i++ )
			IndexSet_Add( shdSets[sharers[sharer_i]], domainInd );
	}
#if 0
	for( p_i = 0; p_i < nIncRanks; p_i++ ) {
		for( s_i = 0; s_i < sync->nSrcs[p_i]; s_i++ )
			IndexSet_Add( shdSets[p_i], sync->srcs[p_i][s_i] );
		for( s_i = 0; s_i < sync->nSnks[p_i]; s_i++ )
			IndexSet_Add( shdSets[p_i], sync->snks[p_i][s_i] );
	}
#endif

	/* Build range sets of shadowed elements. */
	for( p_i = 0; p_i < nIncRanks; p_i++ ) {
		unsigned	nVerts;
		unsigned*	verts;
		unsigned	e_i;

		/* Extract vertices. */
		IndexSet_GetMembers( shdSets[p_i], &nVerts, &verts );

		/* Replace vertex index set with element index set. */
		shdSets[p_i] = IndexSet_New( MeshTopology_GetLocalSize( self, self->nDims ) );

		/* Build initial set of elements. */
		for( e_i = 0; e_i < nVerts; e_i++ ) {
			unsigned	nIncEls;
			unsigned*	incEls;
			unsigned	inc_i;

			nIncEls = self->nIncEls[MT_VERTEX][self->nDims][verts[e_i]];
			incEls = self->incEls[MT_VERTEX][self->nDims][verts[e_i]];
			for( inc_i = 0; inc_i < nIncEls; inc_i++ ) {
				/* Only store local elements. */
				if( incEls[inc_i] >= MeshTopology_GetLocalSize( self, self->nDims ) )
					continue;

				IndexSet_Add( shdSets[p_i], incEls[inc_i] );
			}
		}

		/* Free vertex array. */
		FreeArray( verts );

		/* Expand the elements to satisfy shadow depth. */
		if( self->shadowDepth > 1 ) {
			IndexSet*	bndIndSet;
			unsigned	d_i;

			bndIndSet = IndexSet_DeepCopy( shdSets[p_i] );
			for( d_i = 1; d_i < self->shadowDepth; d_i++ )
				MeshTopology_ExpandShadows( self, (IndexSet*)shdSets[p_i], bndIndSet );
			FreeObject( bndIndSet );
		}

		/* Get members. */
		IndexSet_GetMembers( shdSets[p_i], nShadowedEls[self->nDims] + p_i, 
				     shadowedEls[self->nDims] + p_i );
		FreeObject( shdSets[p_i] );

		/* Clear externals. */
		nExternalEls[self->nDims][p_i] = 0;
		externalEls[self->nDims][p_i] = NULL;
	}
}

void MeshTopology_ExpandShadows( MeshTopology* self, IndexSet* shadows, IndexSet* boundary ) {
	unsigned	nLocals = MeshTopology_GetLocalSize( self, self->nDims );
	unsigned	nBndEls;
	unsigned*	bndEls;
	unsigned	b_i;

	assert( self );
	assert( shadows );
	assert( boundary );

	IndexSet_GetMembers( boundary, &nBndEls, &bndEls );
	FreeObject( boundary );
	boundary = IndexSet_New( nLocals );

	for( b_i = 0; b_i < nBndEls; b_i++ ) {
		unsigned	nNbrs = self->nIncEls[self->nDims][self->nDims][bndEls[b_i]];
		unsigned*	nbrs = self->incEls[self->nDims][self->nDims][bndEls[b_i]];
		unsigned	n_i;

		for( n_i = 0; n_i < nNbrs; n_i++ ) {
			if( nbrs[n_i] >= nLocals || IndexSet_IsMember( shadows, nbrs[n_i] ) )
				continue;
			IndexSet_Add( boundary, nbrs[n_i] );
		}
	}
	FreeArray( bndEls );

	IndexSet_GetMembers( boundary, &nBndEls, &bndEls );
	for( b_i = 0; b_i < nBndEls; b_i++ )
		IndexSet_Add( shadows, bndEls[b_i] );
	FreeArray( bndEls );
}

void MeshTopology_BuildMidShadowedElements( MeshTopology* self, MeshTopology_Dim dim, 
					    unsigned** nShadowedEls, unsigned*** shadowedEls, 
					    unsigned** nExternalEls, unsigned*** externalEls )
{
	Decomp_Sync*	sync;
	CommTopology	*commTopo, *prevCommTopo;
	unsigned	nIncRanks;
	unsigned	nLocals;
	unsigned	prevDim, prevRank;
	IndexSet	*iSet, *extSet;
	RangeSet*	remSet;
	unsigned	p_i, e_i;

	assert( self );
	assert( dim < self->nTDims );
	assert( nShadowedEls && shadowedEls );
	assert( nShadowedEls[dim] && shadowedEls[dim] );

	if( !self->domains[dim] )
		return;

	sync = MeshTopology_GetSync( self, dim );
	commTopo = MeshTopology_GetCommTopology( self, dim );
	nLocals = MeshTopology_GetLocalSize( self, dim );
	nIncRanks = CommTopology_GetIncidenceSize( commTopo );
	prevDim = self->nDims;
	prevCommTopo = MeshTopology_GetCommTopology( self, prevDim );

	for( p_i = 0; p_i < nIncRanks; p_i++ ) {
		/* Map to previous dimension. */
		prevRank = CommTopology_LocalToGlobal( commTopo, p_i );
		insist( CommTopology_GlobalToLocal( prevCommTopo, prevRank, &prevRank ) );

		/* If there are no elements in this dimension or we have no incidence relations between
		   these levels, skip it. */
		if( MeshTopology_GetLocalSize( self, dim ) == 0 || 
		    !self->nIncEls[prevDim][dim] )
		{
			nShadowedEls[dim][p_i] = 0;
			shadowedEls[dim][p_i] = NULL;
			continue;
		}

		/* Build a set representing all our domains already marked as remote. */
		remSet = RangeSet_New();
		RangeSet_SetIndices( remSet, sync->nSnks[p_i], sync->snks[p_i] );
		RangeSet_AddIndices( remSet, sync->nSrcs[p_i], sync->srcs[p_i] );

		/* We'll need an index set to store indices as we build them. */
		iSet = IndexSet_New( nLocals );
		extSet = IndexSet_New( MeshTopology_GetRemoteSize( self, dim ) );

		/* Build indices on this level based on previous level shadowed elements. */
		for( e_i = 0; e_i < nShadowedEls[prevDim][prevRank]; e_i++ ) {
			unsigned	nIncEls;
			unsigned*	incEls;
			unsigned	inc_i;

			nIncEls = self->nIncEls[prevDim][dim][shadowedEls[prevDim][prevRank][e_i]];
			incEls = self->incEls[prevDim][dim][shadowedEls[prevDim][prevRank][e_i]];
			for( inc_i = 0; inc_i < nIncEls; inc_i++ ) {
				/* Skip elements already marked as remote by target processor. */
				if( RangeSet_HasIndex( remSet, incEls[inc_i] ) )
					continue;

				if( incEls[inc_i] < nLocals )
					IndexSet_Add( iSet, incEls[inc_i] );
				else
					IndexSet_Add( extSet, incEls[inc_i] - nLocals );
			}
		}

		/* Free the remote set. */
		FreeObject( remSet );

		/* Get members. */
		IndexSet_GetMembers( iSet, nShadowedEls[dim] + p_i, shadowedEls[dim] + p_i );
		FreeObject( iSet );
		IndexSet_GetMembers( extSet, nExternalEls[dim] + p_i, externalEls[dim] + p_i );
		FreeObject( extSet );
	}
}

void MeshTopology_BuildShadowedIncidence( MeshTopology* self, MeshTopology_Dim dim, 
					  unsigned* nShadowedEls, unsigned** shadowedEls, 
					  unsigned* nShadowedInc, Stg_Byte** shadowedInc )
{
	CommTopology*	commTopo;
	unsigned	nIncRanks;
	unsigned	nBytes;
	Stg_Byte*	bytes;
	unsigned	p_i, d_i;

	assert( self );
	assert( dim < self->nTDims );
	assert( nShadowedEls && shadowedEls );
	assert( nShadowedInc && shadowedInc );

	if( !self->domains[dim] )
		return;

	commTopo = MeshTopology_GetCommTopology( self, dim );
	nIncRanks = CommTopology_GetIncidenceSize( commTopo );

	for( p_i = 0; p_i < nIncRanks; p_i++ ) {
		nShadowedInc[p_i] = 0;
		shadowedInc[p_i] = NULL;

		for( d_i = 0; d_i < dim; d_i++ ) {
			MeshTopology_PickleIncidence( self, dim, d_i, 
						      nShadowedEls[p_i], shadowedEls[p_i], 
						      &nBytes, &bytes );
			shadowedInc[p_i] = ReallocArray( shadowedInc[p_i], Stg_Byte, 
							 nShadowedInc[p_i] + nBytes + sizeof(unsigned) );
			memcpy( shadowedInc[p_i] + nShadowedInc[p_i], 
				&nBytes, sizeof(unsigned) );
			memcpy( shadowedInc[p_i] + nShadowedInc[p_i] + sizeof(unsigned), 
				bytes, nBytes * sizeof(Stg_Byte) );
			nShadowedInc[p_i] += nBytes + sizeof(unsigned);
			FreeArray( bytes );
		}
	}
}

void MeshTopology_ExchangeShadowedElements( MeshTopology* self, MeshTopology_Dim dim, 
					    unsigned* nShadowedEls, unsigned** shadowedEls )
{
	Decomp_Sync*	sync;
	CommTopology*	commTopo;
	unsigned	nIncRanks, *incRanks, *ranks;
	unsigned	*nShadowEls, **shadowEls;
	unsigned	oldSize, newSize;
	unsigned	d_i, p_i, e_i;

	assert( self );
	assert( dim < self->nTDims );
	assert( nShadowedEls && shadowedEls );

	if( !self->domains[dim] )
		return;

	sync = MeshTopology_GetSync( self, dim );
	commTopo = MeshTopology_GetCommTopology( self, dim );
	CommTopology_GetIncidence( commTopo, &nIncRanks, &incRanks );

	for( p_i = 0; p_i < nIncRanks; p_i++ ) {
		for( e_i = 0; e_i < nShadowedEls[p_i]; e_i++ )
			shadowedEls[p_i][e_i] = Decomp_Sync_DomainToGlobal( sync, shadowedEls[p_i][e_i] );
	}

	oldSize = Decomp_Sync_GetDomainSize( sync );
	nShadowEls = AllocArray( unsigned, nIncRanks );
	shadowEls = AllocArray( unsigned*, nIncRanks );
	MeshTopology_ExchangeElements( self, dim, 
				       nShadowedEls, shadowedEls, 
				       nShadowEls, shadowEls );
	ranks = AllocArray( unsigned, nIncRanks );
	for( p_i = 0; p_i < nIncRanks; p_i++ )
		ranks[p_i] = p_i;
	Decomp_Sync_AddSinks( sync, nIncRanks, ranks, nShadowedEls, shadowedEls );
	Decomp_Sync_AddSources( sync, nIncRanks, ranks, nShadowEls, shadowEls );
	Decomp_Sync_BuildShared( sync );
	FreeArray( ranks );
	FreeArray( nShadowEls );
	FreeArray2D( nIncRanks, shadowEls );
	newSize = Decomp_Sync_GetDomainSize( sync );

	for( p_i = 0; p_i < nIncRanks; p_i++ ) {
		for( e_i = 0; e_i < nShadowedEls[p_i]; e_i++ )
			insist( Decomp_Sync_GlobalToDomain( sync, shadowedEls[p_i][e_i], shadowedEls[p_i] + e_i ) );
	}

	for( d_i = 0; d_i < self->nTDims; d_i++ ) {
		if( self->nIncEls[dim][d_i] ) {
			self->nIncEls[dim][d_i] = ReallocArray( self->nIncEls[dim][d_i], 
								unsigned, 
								newSize );
			memset( self->nIncEls[dim][d_i] + oldSize, 0, 
				(newSize - oldSize) * sizeof(unsigned) );
		}

		if( self->incEls[dim][d_i] ) {
			self->incEls[dim][d_i] = ReallocArray( self->incEls[dim][d_i], 
							       unsigned*, 
							       newSize );
			memset( self->incEls[dim][d_i] + oldSize, 0, 
				(newSize - oldSize) * sizeof(unsigned*) );
		}
	}
}

void MeshTopology_ExchangeExternalElements( MeshTopology* self, MeshTopology_Dim dim, 
					    unsigned* nExternalEls, unsigned** externalEls, 
					    unsigned** nShadowedEls, unsigned*** shadowedEls )
{
	unsigned	tag = 98403;
	Decomp_Sync*	sync;
	Decomp*		decomp;
	CommTopology*	commTopo;
	unsigned	nIncRanks, nOldIncRanks;
	MPI_Comm	comm;
	unsigned	rank, nProcs;
	RangeSet**	shdSets;
	unsigned	nLocals, *locals;
	unsigned	nRemotes, *remotes;
	unsigned	nBytes;
	Stg_Byte*	bytes;
	RangeSet	*extSet, *lSet, *remExtSet;
	unsigned	*nExtEls, **extEls;
	unsigned	p_i, e_i;

	assert( self );
	assert( dim < self->nTDims );
	assert( nExternalEls && externalEls );

	if( !self->domains[dim] )
		return;

	sync = MeshTopology_GetSync( self, dim );
	decomp = Decomp_Sync_GetDecomp( sync );
	commTopo = MeshTopology_GetCommTopology( self, dim );
	nIncRanks = CommTopology_GetIncidenceSize( commTopo );
	comm = Decomp_GetComm( decomp );
	MPI_Comm_size( comm, (int*)&nProcs );
	MPI_Comm_rank( comm, (int*)&rank );

	for( p_i = 0; p_i < nIncRanks; p_i++ ) {
		for( e_i = 0; e_i < nExternalEls[p_i]; e_i++ ) {
			externalEls[p_i][e_i] += Decomp_Sync_GetLocalSize( sync );
			externalEls[p_i][e_i] = Decomp_Sync_DomainToGlobal( sync, externalEls[p_i][e_i] );
		}
	}

	nExtEls = AllocArray( unsigned, nIncRanks );
	extEls = AllocArray( unsigned*, nIncRanks );
	MeshTopology_ExchangeElements( self, dim, 
				       nExternalEls, externalEls, 
				       nExtEls, extEls );
	extSet = RangeSet_New();
	for( p_i = 0; p_i < nIncRanks; p_i++ )
		RangeSet_AddIndices( extSet, nExtEls[p_i], extEls[p_i] );
	FreeArray( nExtEls );
	FreeArray2D( nIncRanks, extEls );
	lSet = RangeSet_New();
	Decomp_GetLocals( decomp, &nLocals, &locals );
	RangeSet_SetIndices( lSet, nLocals, locals );
	remExtSet = RangeSet_New();
	Decomp_Sync_GetRemotes( sync, &nRemotes, &remotes );
	RangeSet_AddIndices( remExtSet, nRemotes, remotes );
	RangeSet_Subtraction( extSet, lSet );
	RangeSet_Subtraction( extSet, remExtSet );

	shdSets = AllocArray( RangeSet*, nIncRanks );
	for( p_i = 0; p_i < nIncRanks; p_i++ ) {
		for( e_i = 0; e_i < (*nShadowedEls)[p_i]; e_i++ )
			(*shadowedEls)[p_i][e_i] = Decomp_Sync_DomainToGlobal( sync, (*shadowedEls)[p_i][e_i] );
		shdSets[p_i] = RangeSet_New();
		RangeSet_SetIndices( shdSets[p_i], (*nShadowedEls)[p_i], (*shadowedEls)[p_i] );
		for( e_i = 0; e_i < (*nShadowedEls)[p_i]; e_i++ )
			insist( Decomp_Sync_GlobalToDomain( sync, (*shadowedEls)[p_i][e_i], (*shadowedEls)[p_i] + e_i) );
	}

	nOldIncRanks = nIncRanks;
	RangeSet_Clear( remExtSet );
	for( p_i = 0; p_i < nProcs; p_i++ ) {
		Stg_Byte	state;
		unsigned	nInds, *inds;
		unsigned	localRank;

		if( rank == p_i )
			RangeSet_Pickle( extSet, &nBytes, &bytes );
		MPI_Bcast( &nBytes, 1, MPI_UNSIGNED, p_i, comm );
		if( rank != p_i )
			bytes = AllocArray( Stg_Byte, nBytes );
		MPI_Bcast( bytes, nBytes, MPI_BYTE, p_i, comm );
		if( rank != p_i ) {
			RangeSet_Unpickle( remExtSet, nBytes, bytes );
			FreeArray( bytes );

			RangeSet_Intersection( remExtSet, lSet );
			if( CommTopology_GlobalToLocal( commTopo, p_i, &localRank ) )
				RangeSet_Subtraction( remExtSet, shdSets[localRank] );
			if( RangeSet_GetSize( remExtSet ) ) {
				state = 1;
				MPI_Gather( &state, 1, MPI_BYTE, NULL, 1, MPI_BYTE, p_i, comm );

				RangeSet_Pickle( remExtSet, &nBytes, &bytes );
				MPI_Send( &nBytes, 1, MPI_UNSIGNED, p_i, tag, comm );
				MPI_Send( bytes, nBytes, MPI_BYTE, p_i, tag, comm );
				FreeArray( bytes );

				inds = NULL;
				RangeSet_GetIndices( remExtSet, &nInds, &inds );
				if( !CommTopology_GlobalToLocal( commTopo, p_i, &localRank ) ) {
					Decomp_Sync_AddRemoteRanks( sync, 1, &p_i );
					insist( CommTopology_GlobalToLocal( commTopo, p_i, &localRank ) );
					nIncRanks = CommTopology_GetIncidenceSize( Decomp_Sync_GetCommTopology( sync ) );
					*nShadowedEls = ReallocArray( *nShadowedEls, unsigned, nIncRanks );
					*shadowedEls = ReallocArray( *shadowedEls, unsigned*, nIncRanks );
					(*nShadowedEls)[localRank] = 0;
					(*shadowedEls)[localRank] = NULL;
				}
				Decomp_Sync_AddSinks( sync, 1, &localRank, &nInds, &inds );

				(*shadowedEls)[localRank] = ReallocArray( (*shadowedEls)[localRank], unsigned, 
									  (*nShadowedEls)[localRank] + nInds );
				for( e_i = 0; e_i < nInds; e_i++ )
					insist( Decomp_Sync_GlobalToDomain( sync, inds[e_i], inds + e_i ) );
				memcpy( (*shadowedEls)[localRank] + (*nShadowedEls)[localRank], inds, 
					nInds * sizeof(unsigned) );
				(*nShadowedEls)[localRank] += nInds;

				FreeArray( inds );
			}
			else {
				state = 0;
				MPI_Gather( &state, 1, MPI_BYTE, NULL, 1, MPI_BYTE, p_i, comm );
			}
		}
		else {
			Stg_Byte*	states;
			unsigned	nActive, *active;
			MPI_Status	status;
			unsigned	newSize, oldSize;
			unsigned	nNewRanks, *newRanks;
			unsigned	p_j, d_i;

			oldSize = Decomp_Sync_GetDomainSize( self->domains[dim] );
			newSize = oldSize;

			FreeArray( bytes );
			state = 0;
			states = AllocArray( Stg_Byte, nProcs );
			MPI_Gather( &state, 1, MPI_BYTE, states, 1, MPI_BYTE, p_i, comm );
			nActive = 0;
			nNewRanks = 0;
			for( p_j = 0; p_j < nProcs; p_j++ ) {
				if( states[p_j] ) {
					nActive++;
					if( !CommTopology_GlobalToLocal( commTopo, p_j, &localRank ) )
						nNewRanks++;
				}
			}
			active = AllocArray( unsigned, nActive );
			newRanks = AllocArray( unsigned, nNewRanks );
			nActive = 0;
			nNewRanks = 0;
			for( p_j = 0; p_j < nProcs; p_j++ ) {
				if( states[p_j] ) {
					active[nActive++] = p_j;
					if( !CommTopology_GlobalToLocal( commTopo, p_j, &localRank ) )
						newRanks[nNewRanks++] = p_j;
				}
			}
			FreeArray( states );

			Decomp_Sync_AddRemoteRanks( sync, nNewRanks, newRanks );
			if( nNewRanks ) {
				nIncRanks = CommTopology_GetIncidenceSize( Decomp_Sync_GetCommTopology( sync ) );
				*nShadowedEls = ReallocArray( *nShadowedEls, unsigned, nIncRanks );
				*shadowedEls = ReallocArray( *shadowedEls, unsigned*, nIncRanks );
				for( p_j = 0; p_j < nNewRanks; p_j++ ) {
					insist( CommTopology_GlobalToLocal( commTopo, newRanks[p_j], &localRank ) );
					(*nShadowedEls)[localRank] = 0;
					(*shadowedEls)[localRank] = NULL;
				}
			}
			FreeArray( newRanks );

			for( p_j = 0; p_j < nActive; p_j++ ) {
				MPI_Recv( &nBytes, 1, MPI_UNSIGNED, active[p_j], tag, comm, &status );
				bytes = AllocArray( Stg_Byte, nBytes );
				MPI_Recv( bytes, nBytes, MPI_BYTE, active[p_j], tag, comm, &status );
				RangeSet_Unpickle( extSet, nBytes, bytes );
				FreeArray( bytes );
				inds = NULL;
				RangeSet_GetIndices( extSet, &nInds, &inds );
				insist( CommTopology_GlobalToLocal( commTopo, active[p_j], &localRank ) );
				Decomp_Sync_AddSources( sync, 1, &localRank, &nInds, &inds );
				FreeArray( inds );

				newSize += nInds;
			}

			FreeArray( active );

			for( d_i = 0; d_i < self->nTDims; d_i++ ) {
				if( self->nIncEls[dim][d_i] ) {
					self->nIncEls[dim][d_i] = ReallocArray( self->nIncEls[dim][d_i], 
										unsigned, 
										newSize );
					memset( self->nIncEls[dim][d_i] + oldSize, 0, 
						(newSize - oldSize) * sizeof(unsigned) );
				}

				if( self->incEls[dim][d_i] ) {
					self->incEls[dim][d_i] = ReallocArray( self->incEls[dim][d_i], 
									       unsigned*, 
									       newSize );
					memset( self->incEls[dim][d_i] + oldSize, 0, 
						(newSize - oldSize) * sizeof(unsigned*) );
				}
			}
		}
	}

	FreeObject( extSet );
	FreeObject( remExtSet );
	FreeObject( lSet );

	for( p_i = 0; p_i < nOldIncRanks; p_i++ )
		FreeObject( shdSets[p_i] );
	FreeArray( shdSets );

	Decomp_Sync_BuildShared( sync );
}

void MeshTopology_ExchangeShadowedIncidence( MeshTopology* self, MeshTopology_Dim dim, 
					     unsigned* nShadowedInc, Stg_Byte** shadowedInc )
{
	CommTopology*	commTopo;
	unsigned	nIncRanks;
	unsigned*	nRemBytes;
	Stg_Byte**	remBytes;
	unsigned	offs;
	unsigned	p_i, d_i;

	assert( self );
	assert( dim < self->nTDims );
	assert( nShadowedInc && shadowedInc );

	if( !self->domains[dim] )
		return;

	commTopo = MeshTopology_GetCommTopology( self, dim );
	nIncRanks = CommTopology_GetIncidenceSize( commTopo );

	CommTopology_Alltoall( commTopo, 
			       nShadowedInc, shadowedInc, 
			       &nRemBytes, &remBytes, 
			       sizeof(Stg_Byte) );

	for( p_i = 0; p_i < nIncRanks; p_i++ ) {
		if( !nRemBytes[p_i] )
			continue;

		offs = 0;
		for( d_i = 0; d_i < dim; d_i++ ) {
			unsigned	nBytes;

			nBytes = ((unsigned*)(remBytes[p_i] + offs))[0];
			MeshTopology_UnpickleIncidence( self, dim, d_i, 
							nBytes, 
							remBytes[p_i] + offs + sizeof(unsigned) );
			offs += nBytes + sizeof(unsigned);
		}
		FreeArray( remBytes[p_i] );
	}

	FreeArray( nRemBytes );
	FreeArray( remBytes );
}

void MeshTopology_ExchangeElements( MeshTopology* self, MeshTopology_Dim dim, 
				    unsigned* nSendElements, unsigned** sendElements, 
				    unsigned* nRecvElements, unsigned** recvElements )
{
	Decomp_Sync*	sync;
	CommTopology*	commTopo;
	unsigned	nIncRanks;
	unsigned	*nSrcBytes, *nSnkBytes;
	Stg_Byte	**srcBytes, **snkBytes;
	RangeSet	*shdSet, *tmpSet;
	unsigned	p_i;

	assert( self );
	assert( dim < self->nTDims );
	assert( nSendElements && sendElements );
	assert( nRecvElements && recvElements );

	/* Shortcuts. */
	sync = MeshTopology_GetSync( self, dim );
	commTopo = MeshTopology_GetCommTopology( self, dim );
	nIncRanks = CommTopology_GetIncidenceSize( commTopo );

	/* Need a shadow set to pickle/unpickle all our sets and one as a temporary. */
	shdSet = RangeSet_New();
	tmpSet = RangeSet_New();

	/* Allocate for byte arrays. */
	nSnkBytes = AllocArray( unsigned, nIncRanks );
	snkBytes = AllocArray( Stg_Byte*, nIncRanks );

	/* Pickle everything and store in byte arrays. */
	for( p_i = 0; p_i < nIncRanks; p_i++ ) {
		RangeSet_SetIndices( shdSet, nSendElements[p_i], sendElements[p_i] );
		RangeSet_Pickle( shdSet, nSnkBytes + p_i, snkBytes + p_i );
	}

	/* No longer need shadow range set. */
	FreeObject( shdSet );

	/* Transfer shadow elements. */
	CommTopology_Alltoall( commTopo, 
			       nSnkBytes, (void**)snkBytes, 
			       &nSrcBytes, (void***)&srcBytes, 
			       sizeof(Stg_Byte) );

	/* Don't need sink bytes anymore. */
	FreeArray( nSnkBytes );
	FreeArray2D( nIncRanks, snkBytes );

	/* Unpickle sourced bytes. */
	for( p_i = 0; p_i < nIncRanks; p_i++ ) {
		RangeSet_Unpickle( tmpSet, nSrcBytes[p_i], srcBytes[p_i] );
		recvElements[p_i] = NULL;
		RangeSet_GetIndices( tmpSet, nRecvElements + p_i, recvElements + p_i );
	}

	/* Don't need source bytes or temporary set. */
	FreeArray( nSrcBytes );
	FreeArray( srcBytes );
	FreeObject( tmpSet );
}

void MeshTopology_PickleIncidence( MeshTopology* self, MeshTopology_Dim fromDim, MeshTopology_Dim toDim, 
				   unsigned nEls, unsigned* els, 
				   unsigned* nBytes, Stg_Byte** bytes )
{
	Decomp_Sync	*fromSync, *toSync;
	unsigned	nEntries, *entries;
	unsigned	curEntry;
	unsigned	e_i;

	assert( self );
	assert( fromDim < self->nTDims );
	assert( toDim < self->nTDims );
	assert( !self->nIncEls[fromDim][toDim] || self->incEls[fromDim][toDim] );
	assert( !nEls || els );
	assert( nBytes );
	assert( bytes );

	/* If we have no incidence in this dimension-dimension tuple, clear results. */
	if( !self->nIncEls[fromDim][toDim] ) {
		*nBytes = 0;
		*bytes = NULL;
		return;
	}

	/*
	** Linearise this incidence relation as follows: 
	** {number of elements, [(global element index, number of incident elements, [incident elements])]}
	*/

	fromSync = self->domains[fromDim];
	toSync = self->domains[toDim];

	/* Count the number of entries. */
	nEntries = 1 + 2 * nEls;
	for( e_i = 0; e_i < nEls; e_i++ )
		nEntries += self->nIncEls[fromDim][toDim][e_i];

	/* Allocate for pickled incidence. */
	entries = AllocArray( unsigned, nEntries );

	/* Pack entries. */
	entries[0] = nEls;
	curEntry = 1;
	for( e_i = 0; e_i < nEls; e_i++ ) {
		unsigned	nIncEls, *incEls;
		unsigned	inc_i;

		nIncEls = self->nIncEls[fromDim][toDim][els[e_i]];
		incEls = self->incEls[fromDim][toDim][els[e_i]];
		entries[curEntry++] = Decomp_Sync_DomainToGlobal( fromSync, els[e_i] );
		entries[curEntry++] = nIncEls;
		for( inc_i = 0; inc_i < nIncEls; inc_i++ )
			entries[curEntry++] = Decomp_Sync_DomainToGlobal( toSync, incEls[inc_i] );
	}

	/* Store results. */
	*nBytes = nEntries * sizeof(unsigned);
	*bytes = (Stg_Byte*)entries;
}

void MeshTopology_UnpickleIncidence( MeshTopology* self, MeshTopology_Dim fromDim, MeshTopology_Dim toDim, 
				     unsigned nBytes, Stg_Byte* bytes )
{
	Decomp_Sync	*fromSync, *toSync;
	unsigned	nEntries, *entries;
	unsigned	curEntry;
	unsigned	nEls;
	unsigned	e_i;

	assert( self );
	assert( fromDim < self->nTDims );
	assert( toDim < self->nTDims );
	assert( !nBytes || bytes );
	assert( nBytes % sizeof(unsigned) == 0 );

	/* If we have an empt byte set return. */
	if( !nBytes )
		return;

	/* Setup entry array. */
	fromSync = self->domains[fromDim];
	toSync = self->domains[toDim];
	nEntries = nBytes / sizeof(unsigned);
	entries = (unsigned*)bytes;

	/* Transfer incidence. */
	nEls = entries[0];
	curEntry = 1;
	for( e_i = 0; e_i < nEls; e_i++ ) {
		unsigned	nIncEls, **incEls;
		unsigned	el;
		unsigned	inc_i;

		/* Convert element index to domain. */
		insist( Decomp_Sync_GlobalToDomain( fromSync, entries[curEntry++], &el ) );

		/* Get the size. */
		nIncEls = entries[curEntry++];
		self->nIncEls[fromDim][toDim][el] = nIncEls;

		/* Resize array. */
		incEls = self->incEls[fromDim][toDim];
		incEls[el] = ReallocArray( incEls[el], unsigned, nIncEls );

		/* Copy incidence, converting to domain indices. */
		for( inc_i = 0; inc_i < nIncEls; inc_i++ )
			insist( Decomp_Sync_GlobalToDomain( toSync, entries[curEntry++], incEls[el] + inc_i ) );
	}
}

void MeshTopology_Destruct( MeshTopology* self ) {
	unsigned	d_i;

	assert( self );

	self->shadowDepth = 0;

	for( d_i = 0; d_i < self->nTDims; d_i++ ) {
		unsigned	d_j;

		if( self->domains && self->domains[d_i] ) {
			Stg_Class_RemoveRef( self->domains[d_i] );
			self->domains[d_i] = NULL;
		}

		for( d_j = 0; d_j < self->nTDims; d_j++ ) {
			if( self->incEls )
				KillArray( self->incEls[d_i][d_j] );
			if( self->nIncEls )
				KillArray( self->nIncEls[d_i][d_j] );
		}
	}
	KillArray( self->incEls );
	KillArray( self->nIncEls );
	KillArray( self->domains );
}

#ifndef NDEBUG
Bool MeshTopology_ValidateElements( MeshTopology* self, MeshTopology_Dim dim, unsigned nEls, unsigned* els ) {
	Decomp_Sync*	sync;
	Decomp*		decomp;
	RangeSet	*exSet, *newSet;
	unsigned	nInds, *inds;

	assert( self );
	assert( !nEls || els );

	sync = self->domains[dim];
	if( !sync )
		return True;

	exSet = RangeSet_New();
	decomp = Decomp_Sync_GetDecomp( sync );
	Decomp_GetLocals( decomp, &nInds, &inds );
	RangeSet_SetIndices( exSet, nInds, inds );
	FreeArray( inds );
	Decomp_Sync_GetRemotes( sync, &nInds, &inds );
	RangeSet_AddIndices( exSet, nInds, inds );

	newSet = RangeSet_New();
	RangeSet_AddIndices( newSet, nEls, els );

	RangeSet_Intersection( newSet, exSet );
	if( RangeSet_GetSize( newSet ) )
		return False;

	return True;
}
#endif

#if 0
void MeshTopology_BuildExternalElements( MeshTopology* self, MeshTopology_Dim dim, 
					 unsigned* nShadowedEls, unsigned** shadowedEls, 
					 unsigned* nExtEls, unsigned** extEls )
{
	Decomp_Sync*	sync;
	Decomp*		decomp;
	unsigned	nLocals;
	unsigned	p_i;

	assert( self );
	assert( dim < self->nTDims );
	assert( nExtEls && extEls );

	sync = MeshTopology_GetSync( self, dim );
	decomp = Decomp_Sync_GetDecomp( sync );
	nLocals = Decomp_GetLocalSize( decomp );
	for( p_i = 0; p_i < nIncRanks; p_i++ ) {
		IndexSet*	remSet;
		unsigned	e_i;

		remSet = IndexSet_New( Decomp_Sync_GetRemoteSize( sync ) );
		for( e_i = 0; e_i < nShadowedEls[p_i]; e_i++ ) {
			if( shadowedEls[p_i][e_i] >= nLocals ) {
				unsigned	remInd;

				remInd = shadowedEls[p_i][e_i] - nLocals;
				IndexSet_Add( remSet, remInd );
			}
		}

		IndexSet_GetMembers( remSet, nExtEls + p_i, extEls + p_i );
		FreeObject( remSet );
	}
}

void MeshTopology_ExchangeExternalElements( MeshTopology* self, MeshTopology_Dim dim, 
					    unsigned* nShadowedEls, unsigned** shadowedEls, 
					    unsigned* nExtEls, unsigned** extEls )
{
	unsigned	*nSplitExts, **splitExts;

	ownSets = AllocArray( IndexSet*, nIncRanks );
	for( p_i = 0; p_i < nIncRanks; p_i++ )
		ownSets[p_i] = IndexSet_New( Decomp_Sync_GetRemoteSize( sync ) );

	nSplitExts = AllocArray( unsigned, nIncRanks );
	splitExts = AllocArray( unsigned*, nIncRanks );

	for( p_i = 0; p_i < nIncRanks; p_i++ ) {
		for( e_i = 0; e_i < nExtEls[p_i]; e_i++ ) {
			unsigned	owner;

			owner = Decomp_Sync_GetOwner( sync, extEls[p_i][e_i] );
			IndexSet_Add( ownSets[owner], extEls[p_i][e_i] );
		}
		for( p_j = 0; p_j < nIncRanks; p_j++ ) {
			IndexSet_GetMembers( ownSets[p_j], nSplitExts + p_j, splitExts + p_j );
			IndexSet_Clear( ownSets[p_j] );
		}

		for( p_j = 0; p_j < nIncRanks; p_j++ ) {
			if( !nSplitExts[p_j] )
				continue;
			offs = nSnkEntries[p_j];
			nSnkEntries[p_j] += nSplitExts[p_j] + 2;
			snkEntries[p_j] = ReallocArray( snkEntries[p_j], unsigned, nSnkEntries[p_j] );
			snkEntries[p_j][offs] = CommTopology_LocalToGlobal( commTopo, p_i );
			snkEntries[p_j][offs + 1] = nSplitExts[p_j];
			memcpy( snkEntries[p_j] + offs + 2, splitExts[p_j], nSplitExts[p_j] * sizeof(unsigned) );
		}

		offs = nSrcEntries[p_i];
		for( p_j = 0; p_j < nIncRanks; p_j++ )
			nSrcEntries[p_i] += nSplitExts[p_j] + (nSplitExts[p_j] ? 2 : 0);
		srcEntries[p_i] = ReallocArray( srcEntries[p_i], unsigned, nSrcEntries[p_i] );
		for( p_j = 0; p_j < nIncRanks; p_j++ ) {
			if( !nSplitExts[p_j] )
				continue;
			srcEntries[p_i][offs++] = CommTopology_LocalToGlobal( commTopo, p_j );
			srcEntries[p_i][offs++] = nSplitExts[p_j];
			memcpy( srcEntries[p_i] + offs, splitExts[p_j], nSplitExts[p_j] * sizeof(unsigned) );
			offs += nSplitExts[p_j];
		}
	}

	CommTopology_Alltoall( commTopo, 
			       nSnkEntries, snkEntries, 
			       nRecvSnkEntries, recvSnkEntries, 
			       sizeof(unsigned) );
	CommTopology_Alltoall( commTopo, 
			       nSrcEntries, srcEntries, 
			       nRecvSrcEntries, recvSrcEntries, 
			       sizeof(unsigned) );
	for( p_i = 0; p_i < nIncRanks; p_i++ ) {
		FreeArray( srcEntries[p_i] );
		FreeArray( snkEntries[p_i] );
	}

	rankSet = RangeSet_New();
	newRanks = NULL;
	for( p_i = 0; p_i < nIncRanks; p_i++ ) {
		unsigned	nItems;

		nItems = recvSnkEntries[p_i][0];
		newRanks = ReallocArray( unsigned, nItems );
		offs = 1;
		for( itm_i = 0; itm_i < nItems; itm_i++ ) {
			unsigned	nEls;

			newRanks[itm_i] = recvSnkEntries[p_i][offs++];
			nEls = recvSnkEntries[p_i][offs++];
			offs += nEls;
		}
		RangeSet_AddIndices( rankSet, nItems, newRanks );

		nItems = recvSrcEntries[p_i][0];
		newRanks = ReallocArray( unsigned, nItems );
		offs = 1;
		for( itm_i = 0; itm_i < nItems; itm_i++ ) {
			unsigned	nEls;

			newRanks[itm_i] = recvSrcEntries[offs++];
			nEls = recvSrcEntries[offs++];
			offs += nEls;
		}
		RangeSet_AddIndices( rankSet, nItems, newRanks );
	}
	FreeArray( newRanks );
	newRanks = NULL;
	RangeSet_GetIndices( rankSet, &nNewRanks, &newRanks );
	Decomp_Sync_AddRemoteRanks( sync, nNewRanks, newRanks );
	FreeArray( newRanks );

	for( p_i = 0; p_i < nIncRanks; p_i++ ) {
	}

	FreeArray2D( nIncRanks, snkEntries );
	FreeArray2D( nIncRanks, srcEntries );
	FreeArray( nSnkEntries );
	FreeArray( nSrcEntries );
}
#endif
