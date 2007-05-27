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
** $Id: CartesianGenerator.c 3584 2006-05-16 11:11:07Z PatrickSunter $
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

#include "types.h"
#include "shortcuts.h"
#include "Grid.h"
#include "Decomp.h"
#include "Sync.h"
#include "MeshTopology.h"
#include "Mesh_ElementType.h"
#include "Mesh_HexType.h"
#include "MeshClass.h"
#include "MeshGenerator.h"
#include "Mesh_Algorithms.h"
#include "Mesh_RegularAlgorithms.h"
#include "CartesianGenerator.h"


/* Textual name of this class */
const Type CartesianGenerator_Type = "CartesianGenerator";

/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/

CartesianGenerator* CartesianGenerator_New( Name name ) {
	return _CartesianGenerator_New( sizeof(CartesianGenerator), 
					CartesianGenerator_Type, 
					_CartesianGenerator_Delete, 
					_CartesianGenerator_Print, 
					NULL, 
					(void* (*)(Name))_CartesianGenerator_New, 
					_CartesianGenerator_Construct, 
					_CartesianGenerator_Build, 
					_CartesianGenerator_Initialise, 
					_CartesianGenerator_Execute, 
					_CartesianGenerator_Destroy, 
					name, 
					NON_GLOBAL, 
					CartesianGenerator_SetDimSize, 
					CartesianGenerator_Generate, 
					_CartesianGenerator_SetTopologyParams, 
					_CartesianGenerator_GenElements, 
					_CartesianGenerator_GenFaces, 
					_CartesianGenerator_GenEdges, 
					_CartesianGenerator_GenVertices, 
					_CartesianGenerator_GenElementVertexInc, 
					_CartesianGenerator_GenVolumeEdgeInc, 
					_CartesianGenerator_GenVolumeFaceInc, 
					_CartesianGenerator_GenFaceVertexInc, 
					_CartesianGenerator_GenFaceEdgeInc, 
					_CartesianGenerator_GenEdgeVertexInc, 
					_CartesianGenerator_GenElementTypes );
}

CartesianGenerator* _CartesianGenerator_New( CARTESIANGENERATOR_DEFARGS ) {
	CartesianGenerator* self;
	
	/* Allocate memory */
	assert( sizeOfSelf >= sizeof(CartesianGenerator) );
	self = (CartesianGenerator*)_MeshGenerator_New( MESHGENERATOR_PASSARGS );

	/* Virtual info */
	self->setTopologyParamsFunc = setTopologyParamsFunc;
	self->genElementsFunc = genElementsFunc;
	self->genFacesFunc = genFacesFunc;
	self->genEdgesFunc = genEdgesFunc;
	self->genVerticesFunc = genVerticesFunc;
	self->genElementVertexIncFunc = genElementVertexIncFunc;
	self->genVolumeEdgeIncFunc = genVolumeEdgeIncFunc;
	self->genVolumeFaceIncFunc = genVolumeFaceIncFunc;
	self->genFaceVertexIncFunc = genFaceVertexIncFunc;
	self->genFaceEdgeIncFunc = genFaceEdgeIncFunc;
	self->genEdgeVertexIncFunc = genEdgeVertexIncFunc;
	self->genElementTypesFunc = genElementTypesFunc;

	/* CartesianGenerator info */
	_CartesianGenerator_Init( self );

	return self;
}

void _CartesianGenerator_Init( CartesianGenerator* self ) {
	Stream*	stream;

	assert( self && Stg_CheckType( self, CartesianGenerator ) );

	stream = Journal_Register( Info_Type, self->type );
	Stream_SetPrintingRank( stream, 0 );

	self->comm = NULL;
	self->regular = True;
	self->maxDecompDims = 0;
	self->minDecomp = NULL;
	self->maxDecomp = NULL;
	self->shadowDepth = 1;
	self->crdMin = NULL;
	self->crdMax = NULL;

	self->vertGrid = NULL;
	self->elGrid = NULL;
	self->procGrid = NULL;
	self->origin = NULL;
	self->range = NULL;
	self->vertOrigin = NULL;
	self->vertRange = NULL;
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _CartesianGenerator_Delete( void* meshGenerator ) {
	CartesianGenerator*	self = (CartesianGenerator*)meshGenerator;

	CartesianGenerator_Destruct( self );

	/* Delete the parent. */
	_MeshGenerator_Delete( self );
}

void _CartesianGenerator_Print( void* meshGenerator, Stream* stream ) {
	CartesianGenerator*	self = (CartesianGenerator*)meshGenerator;
	
	/* Set the Journal for printing informations */
	Stream* meshGeneratorStream;
	meshGeneratorStream = Journal_Register( InfoStream_Type, "CartesianGeneratorStream" );

	assert( self && Stg_CheckType( self, CartesianGenerator ) );

	/* Print parent */
	Journal_Printf( stream, "CartesianGenerator (ptr): (%p)\n", self );
	_MeshGenerator_Print( self, stream );
}

void _CartesianGenerator_Construct( void* meshGenerator, Stg_ComponentFactory* cf, void* data ) {
	CartesianGenerator*	self = (CartesianGenerator*)meshGenerator;
	Dictionary*		dict;
	Dictionary_Entry_Value*	tmp;
	char*			rootKey;
	Dictionary_Entry_Value*	sizeList;
	Dictionary_Entry_Value	*minList, *maxList;
	unsigned		maxDecompDims;
	unsigned		*minDecomp, *maxDecomp;
	double			*crdMin, *crdMax;
	unsigned*		size;
	unsigned		shadowDepth;
	Stream*			stream;
	unsigned		d_i;

	assert( self && Stg_CheckType( self, CartesianGenerator ) );
	assert( cf );

	/* Call parent construct. */
	_MeshGenerator_Construct( self, cf, data );

	/* Rip out the components structure as a dictionary. */
	dict = Dictionary_Entry_Value_AsDictionary( Dictionary_Get( cf->componentDict, self->name ) );

	/* Read the sizes. */
	sizeList = Dictionary_Get( dict, "size" );
	assert( sizeList );
	assert( Dictionary_Entry_Value_GetCount( sizeList ) >= self->nDims );
	size = Memory_Alloc_Array_Unnamed( unsigned, self->nDims );
	for( d_i = 0; d_i < self->nDims; d_i++ ) {
		tmp = Dictionary_Entry_Value_GetElement( sizeList, d_i );
		rootKey = Dictionary_Entry_Value_AsString( tmp );

		if( !Stg_StringIsNumeric( rootKey ) )
			tmp = Dictionary_Get( cf->rootDict, rootKey );
		size[d_i] = Dictionary_Entry_Value_AsUnsignedInt( tmp );
	}

	/* Read decomposition restrictions. */
	maxDecompDims = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, "maxDecomposedDims", 0 );

	minList = Dictionary_Get( dict, "minDecompositions" );
	if( minList ) {
		minDecomp = AllocArray( unsigned, self->nDims );
		for( d_i = 0; d_i < self->nDims; d_i++ ) {
			tmp = Dictionary_Entry_Value_GetElement( minList, d_i );
			rootKey = Dictionary_Entry_Value_AsString( tmp );

			if( !Stg_StringIsNumeric( rootKey ) )
				tmp = Dictionary_Get( cf->rootDict, rootKey );
			minDecomp[d_i] = Dictionary_Entry_Value_AsUnsignedInt( tmp );
		}
	}
	else
		minDecomp = NULL;

	maxList = Dictionary_Get( dict, "maxDecompositions" );
	if( maxList ) {
		maxDecomp = AllocArray( unsigned, self->nDims );
		for( d_i = 0; d_i < self->nDims; d_i++ ) {
			tmp = Dictionary_Entry_Value_GetElement( maxList, d_i );
			rootKey = Dictionary_Entry_Value_AsString( tmp );

			if( !Stg_StringIsNumeric( rootKey ) )
				tmp = Dictionary_Get( cf->rootDict, rootKey );
			maxDecomp[d_i] = Dictionary_Entry_Value_AsUnsignedInt( tmp );
		}
	}
	else
		maxDecomp = NULL;

	/* Initial setup. */
	CartesianGenerator_SetTopologyParams( self, size, maxDecompDims, minDecomp, maxDecomp );

	/* Read geometry. */
	minList = Dictionary_Get( dict, "minCoord" );
	maxList = Dictionary_Get( dict, "maxCoord" );
	if( minList && maxList ) {
		assert( Dictionary_Entry_Value_GetCount( minList ) >= self->nDims );
		assert( Dictionary_Entry_Value_GetCount( maxList ) >= self->nDims );
		crdMin = Memory_Alloc_Array_Unnamed( double, self->nDims );
		crdMax = Memory_Alloc_Array_Unnamed( double, self->nDims );
		for( d_i = 0; d_i < self->nDims; d_i++ ) {
			tmp = Dictionary_Entry_Value_GetElement( minList, d_i );
			rootKey = Dictionary_Entry_Value_AsString( tmp );

			if( !Stg_StringIsNumeric( rootKey ) )
				tmp = Dictionary_Get( cf->rootDict, rootKey );
			crdMin[d_i] = Dictionary_Entry_Value_AsDouble( tmp );

			tmp = Dictionary_Entry_Value_GetElement( maxList, d_i );
			rootKey = Dictionary_Entry_Value_AsString( tmp );

			if( !Stg_StringIsNumeric( rootKey ) )
				tmp = Dictionary_Get( cf->rootDict, rootKey );
			crdMax[d_i] = Dictionary_Entry_Value_AsDouble( tmp );
		}

		/* Initial setup. */
		CartesianGenerator_SetGeometryParams( self, crdMin, crdMax );

		/* Free coordinate arrays. */
		FreeArray( crdMin );
		FreeArray( crdMax );
	}

	/* Read and set shadow depth. */
	shadowDepth = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, "shadowDepth", 1 );
	CartesianGenerator_SetShadowDepth( self, shadowDepth );

	/* Read regular flag. */
	self->regular = Stg_ComponentFactory_GetBool( cf, self->name, "regular", True );

	/* Read a general dictionary flag for which processor to watch. */
	stream = Journal_Register( Info_Type, self->type );
	Stream_SetPrintingRank( stream, Dictionary_GetUnsignedInt_WithDefault( cf->rootDict, "rankToWatch", 0 ) );

	/* Free stuff. */
	FreeArray( size );
	FreeArray( minDecomp );
	FreeArray( maxDecomp );
}

void _CartesianGenerator_Build( void* meshGenerator, void* data ) {
	_MeshGenerator_Build( meshGenerator, data );
}

void _CartesianGenerator_Initialise( void* meshGenerator, void* data ) {
}

void _CartesianGenerator_Execute( void* meshGenerator, void* data ) {
}

void _CartesianGenerator_Destroy( void* meshGenerator, void* data ) {
}

void CartesianGenerator_SetDimSize( void* meshGenerator, unsigned nDims ) {
	CartesianGenerator*	self = (CartesianGenerator*)meshGenerator;

	assert( self && Stg_CheckType( self, CartesianGenerator ) );

	_MeshGenerator_SetDimSize( self, nDims );

	self->minDecomp = AllocArray( unsigned, self->nDims );
	memset( self->minDecomp, 0, nDims * sizeof(unsigned) );
	self->maxDecomp = AllocArray( unsigned, self->nDims );
	memset( self->maxDecomp, 0, nDims * sizeof(unsigned) );
}

void CartesianGenerator_Generate( void* meshGenerator, void* _mesh ) {
	CartesianGenerator*	self = (CartesianGenerator*)meshGenerator;
	Stream*			stream = Journal_Register( Info_Type, self->type );
	Mesh*			mesh = (Mesh*)_mesh;
	Grid**			grid;
	unsigned		*localRange, *localOrigin;

	/* Sanity check. */
	assert( self );
	assert( !self->elGrid || mesh );

	Journal_Printf( stream, "Cartesian generator: '%s'\n", self->name );
	Stream_Indent( stream );

	/* If we havn't been given anything, don't do anything. */
	if( self->elGrid ) {
		unsigned	d_i;

		Journal_Printf( stream, "Target mesh: '%s'\n", mesh->name );
		Journal_Printf( stream, "Global element size: %d", self->elGrid->sizes[0] );
		for( d_i = 1; d_i < self->elGrid->nDims; d_i++ )
			Journal_Printf( stream, "x%d", self->elGrid->sizes[d_i] );
		Journal_Printf( stream, "\n" );
		Journal_Printf( stream, "Local offset of rank %d: %d", stream->_printingRank, self->origin[0] );
		for( d_i = 1; d_i < self->elGrid->nDims; d_i++ )
			Journal_Printf( stream, "x%d", self->origin[d_i] );
		Journal_Printf( stream, "\n" );
		Journal_Printf( stream, "Local range of rank %d: %d", stream->_printingRank, self->range[0] );
		for( d_i = 1; d_i < self->elGrid->nDims; d_i++ )
			Journal_Printf( stream, "x%d", self->range[d_i] );
		Journal_Printf( stream, "\n" );

		/* Fill topological values. */
		MeshTopology_SetComm( mesh->topo, self->comm );
		MeshTopology_SetNumDims( mesh->topo, self->elGrid->nDims );
		CartesianGenerator_GenTopo( self, mesh->topo );

		/* Fill geometric values. */
		CartesianGenerator_GenGeom( self, mesh );

		/* Fill element types. */
		CartesianGenerator_GenElementTypes( self, mesh );
	}
	else {
		MeshTopology_SetNumDims( mesh->topo, 0 );
	}

	/* Add extensions to the mesh and fill with cartesian information. */
	ExtensionManager_Add( mesh->info, "vertexGrid", sizeof(Grid*) );
	grid = (Grid**)ExtensionManager_Get( mesh->info, mesh, 
					     ExtensionManager_GetHandle( mesh->info, "vertexGrid" ) );
	*grid = Grid_New();
	Grid_SetNumDims( *grid, self->vertGrid->nDims );
	Grid_SetSizes( *grid, self->vertGrid->sizes );

	ExtensionManager_Add( mesh->info, "elementGrid", sizeof(Grid*) );
	grid = (Grid**)ExtensionManager_Get( mesh->info, mesh, 
					     ExtensionManager_GetHandle( mesh->info, "elementGrid" ) );
	*grid = Grid_New();
	Grid_SetNumDims( *grid, self->elGrid->nDims );
	Grid_SetSizes( *grid, self->elGrid->sizes );

	ExtensionManager_AddArray( mesh->info, "localOrigin", sizeof(unsigned), Mesh_GetDimSize( mesh ) );
	localOrigin = (unsigned*)ExtensionManager_Get( mesh->info, mesh, 
						       ExtensionManager_GetHandle( mesh->info, "localOrigin" ) );
	memcpy( localOrigin, self->origin, Mesh_GetDimSize( mesh ) * sizeof(unsigned) );

	ExtensionManager_AddArray( mesh->info, "localRange", sizeof(unsigned), Mesh_GetDimSize( mesh ) );
	localRange = (unsigned*)ExtensionManager_Get( mesh->info, mesh, 
						      ExtensionManager_GetHandle( mesh->info, "localRange" ) );
	memcpy( localRange, self->range, Mesh_GetDimSize( mesh ) * sizeof(unsigned) );

	Stream_UnIndent( stream );
}

void _CartesianGenerator_SetTopologyParams( void* meshGenerator, unsigned* sizes, 
					    unsigned maxDecompDims, unsigned* minDecomp, unsigned* maxDecomp )
{
	CartesianGenerator*	self = (CartesianGenerator*)meshGenerator;
	unsigned		d_i;

	/* Sanity check. */
	assert( self && Stg_CheckType( self, CartesianGenerator ) );
	assert( !self->nDims || sizes );
	assert( self->nDims <= 3 );

	/* Kill everything we have, topologically. */
	CartesianGenerator_DestructTopology( self );

	/* Don't continue if we have nothing. */
	if( !self->nDims )
		return;

	/* Set the parameters. */
	self->elGrid = Grid_New();
	Grid_SetNumDims( self->elGrid, self->nDims );
	Grid_SetSizes( self->elGrid, sizes );

	self->vertGrid = Grid_New();
	Grid_SetNumDims( self->vertGrid, self->nDims );
	for( d_i = 0; d_i < self->nDims; d_i++ )
		sizes[d_i]++;
	Grid_SetSizes( self->vertGrid, sizes );
	for( d_i = 0; d_i < self->nDims; d_i++ )
		sizes[d_i]--;

	if( minDecomp )
		memcpy( self->minDecomp, minDecomp, self->nDims * sizeof(unsigned) );
	else
		memset( self->minDecomp, 0, self->nDims * sizeof(unsigned) );
	if( maxDecomp )
		memcpy( self->maxDecomp, maxDecomp, self->nDims * sizeof(unsigned) );
	else
		memset( self->maxDecomp, 0, self->nDims * sizeof(unsigned) );
	self->maxDecompDims = maxDecompDims;

	/* As soon as we know the topology, we can decompose. */
	CartesianGenerator_BuildDecomp( self );
}

void _CartesianGenerator_GenElements( void* meshGenerator, MeshTopology* topo, Grid*** grids ) {
	CartesianGenerator*	self = (CartesianGenerator*)meshGenerator;
	Stream*			stream = Journal_Register( Info_Type, self->type );
	Grid*			grid;
	unsigned		nEls;
	unsigned*		els;
	unsigned*		dimInds;
	unsigned		d_i, e_i;

	assert( self && Stg_CheckType( self, CartesianGenerator ) );
	assert( topo );
	assert( grids );

	Journal_Printf( stream, "Generating elements...\n" );
	Stream_Indent( stream );

	grid = Grid_New();
	Grid_SetNumDims( grid, self->elGrid->nDims );
	Grid_SetSizes( grid, self->range );

	nEls = grid->sizes[0];
	for( d_i = 1; d_i < grid->nDims; d_i++ )
		nEls *= grid->sizes[d_i];
	els = Memory_Alloc_Array_Unnamed( unsigned, nEls );

	dimInds = Memory_Alloc_Array_Unnamed( unsigned, self->elGrid->nDims );

	for( e_i = 0; e_i < nEls; e_i++ ) {
		Grid_Lift( grid, e_i, dimInds );
		for( d_i = 0; d_i < grid->nDims; d_i++ )
			dimInds[d_i] += self->origin[d_i];
		els[e_i] = Grid_Project( self->elGrid, dimInds );
	}

	MeshTopology_SetLocalElements( topo, grid->nDims, nEls, els );
	FreeArray( els );
	FreeArray( dimInds );
	FreeObject( grid );

	MPI_Barrier( self->mpiComm );
	Journal_Printf( stream, "... done.\n" );
	Stream_UnIndent( stream );
}

void _CartesianGenerator_GenVertices( void* meshGenerator, MeshTopology* topo, Grid*** grids ) {
	CartesianGenerator*	self = (CartesianGenerator*)meshGenerator;
	Stream*			stream = Journal_Register( Info_Type, self->type );
	unsigned		rank;
	Grid*			globalGrid;
	Grid*			grid;
	unsigned		nEls;
	unsigned		nLocals, *locals;
	unsigned		nRemotes, *remotes;
	unsigned		*dstArray, *dstCount;
	unsigned		*dimInds, *rankInds;
	unsigned		d_i, e_i;

	assert( self && Stg_CheckType( self, CartesianGenerator ) );
	assert( topo );
	assert( grids );

	Journal_Printf( stream, "Generating vertices...\n" );
	Stream_Indent( stream );

	MPI_Comm_rank( Comm_GetMPIComm( self->comm ), (int*)&rank );

	globalGrid = self->vertGrid;

	grid = Grid_New();
	Grid_SetNumDims( grid, self->elGrid->nDims );
	Grid_SetSizes( grid, self->vertRange );

	nEls = grid->sizes[0];
	for( d_i = 1; d_i < grid->nDims; d_i++ )
		nEls *= grid->sizes[d_i];
	locals = Memory_Alloc_Array_Unnamed( unsigned, nEls );
/*
	remotes = Memory_Alloc_Array_Unnamed( unsigned, nEls );
*/

	dimInds = Memory_Alloc_Array_Unnamed( unsigned, self->elGrid->nDims );
	rankInds = Memory_Alloc_Array_Unnamed( unsigned, self->elGrid->nDims );
	Grid_Lift( self->procGrid, rank, rankInds );

	nLocals = 0;
	nRemotes = 0;
	for( e_i = 0; e_i < nEls; e_i++ ) {
		Grid_Lift( grid, e_i, dimInds );
/*
		for( d_i = 0; d_i < grid->nDims; d_i++ ) {
			if( dimInds[d_i] == 0 && rankInds[d_i] > 0 )
				break;
		}
		if( d_i < grid->nDims ) {
			dstArray = remotes;
			dstCount = &nRemotes;
		}
		else {
*/
			dstArray = locals;
			dstCount = &nLocals;
/*
		}
*/

		for( d_i = 0; d_i < grid->nDims; d_i++ )
			dimInds[d_i] += self->vertOrigin[d_i];

		dstArray[(*dstCount)++] = Grid_Project( globalGrid, dimInds );
	}

/*
	MeshTopology_SetLocalElements( topo, 0, nLocals, locals );
	MeshTopology_SetRemoteElements( topo, 0, nRemotes, remotes );
*/
	MeshTopology_SetElements( topo, 0, nLocals, locals );
	FreeArray( locals );
/*
	FreeArray( remotes );
*/
	FreeArray( dimInds );
	FreeArray( rankInds );
	FreeObject( grid );

	MPI_Barrier( self->mpiComm );
	Journal_Printf( stream, "... done.\n" );
	Stream_UnIndent( stream );
}

void _CartesianGenerator_GenEdges( void* meshGenerator, MeshTopology* topo, Grid*** grids ) {
	CartesianGenerator*	self = (CartesianGenerator*)meshGenerator;
	Stream*			stream = Journal_Register( Info_Type, self->type );

	assert( self && Stg_CheckType( self, CartesianGenerator ) );
	assert( topo );
	assert( self->elGrid->nDims >= 2 );
	assert( self->elGrid->nDims <= 3 );

	Journal_Printf( stream, "Generating edges...\n" );
	Stream_Indent( stream );

	if( self->elGrid->nDims == 2 )
		CartesianGenerator_GenEdges2D( self, topo, grids );
	else
		CartesianGenerator_GenEdges3D( self, topo, grids );

	MPI_Barrier( self->mpiComm );
	Journal_Printf( stream, "... done.\n" );
	Stream_UnIndent( stream );
}

void _CartesianGenerator_GenFaces( void* meshGenerator, MeshTopology* topo, Grid*** grids ) {
  /*
	CartesianGenerator*	self = (CartesianGenerator*)meshGenerator;
	Stream*			stream = Journal_Register( Info_Type, self->type );
	Grid*			globalGrid;
	unsigned		nGlobalEls[2];
	Grid*			grid;
	unsigned		nEls[3];
	unsigned*		els;
	unsigned*		dimInds;
	unsigned		d_i, e_i;

	assert( self && Stg_CheckType( self, CartesianGenerator ) );
	assert( topo && Stg_CheckType( topo, MeshTopology ) );
	assert( grids );
	assert( self->elGrid->nDims == 3 );

	Journal_Printf( stream, "Generating faces...\n" );
	Stream_Indent( stream );

	globalGrid = Grid_New();
	Grid_SetNumDims( globalGrid, self->elGrid->nDims );
	self->elGrid->sizes[2]++;
	Grid_SetSizes( globalGrid, self->elGrid->sizes );
	self->elGrid->sizes[2]--;

	grid = Grid_New();
	Grid_SetNumDims( grid, self->elGrid->nDims );
	self->range[2]++;
	Grid_SetSizes( grid, self->range );
	self->range[2]--;

	nEls[0] = grid->sizes[0];
	nGlobalEls[0] = globalGrid->sizes[0];
	for( d_i = 1; d_i < grid->nDims; d_i++ ) {
		nEls[0] *= grid->sizes[d_i];
		nGlobalEls[0] *= globalGrid->sizes[d_i];
	}
	els = Memory_Alloc_Array_Unnamed( unsigned, nEls[0] );

	dimInds = Memory_Alloc_Array_Unnamed( unsigned, self->elGrid->nDims );

	for( e_i = 0; e_i < nEls[0]; e_i++ ) {
		Grid_Lift( grid, e_i, dimInds );
		for( d_i = 0; d_i < grid->nDims; d_i++ )
			dimInds[d_i] += self->origin[d_i];
		els[e_i] = Grid_Project( globalGrid, dimInds );
	}

	self->elGrid->sizes[1]++;
	Grid_SetSizes( globalGrid, self->elGrid->sizes );
	self->elGrid->sizes[1]--;

	self->range[1]++;
	Grid_SetSizes( grid, self->range );
	self->range[1]--;

	nEls[1] = grid->sizes[0];
	nGlobalEls[1] = globalGrid->sizes[0];
	for( d_i = 1; d_i < grid->nDims; d_i++ ) {
		nEls[1] *= grid->sizes[d_i];
		nGlobalEls[1] *= globalGrid->sizes[d_i];
	}
	els = Memory_Realloc_Array( els, unsigned, nEls[0] + nEls[1] );

	for( e_i = 0; e_i < nEls[1]; e_i++ ) {
		Grid_Lift( grid, e_i, dimInds );
		for( d_i = 0; d_i < grid->nDims; d_i++ )
			dimInds[d_i] += self->origin[d_i];
		els[nEls[0] + e_i] = nGlobalEls[0] + Grid_Project( globalGrid, dimInds );
	}

	self->elGrid->sizes[0]++;
	Grid_SetSizes( globalGrid, self->elGrid->sizes );
	self->elGrid->sizes[0]--;

	self->range[0]++;
	Grid_SetSizes( grid, self->range );
	self->range[0]--;

	nEls[2] = grid->sizes[0];
	for( d_i = 1; d_i < grid->nDims; d_i++ )
		nEls[2] *= grid->sizes[d_i];
	els = Memory_Realloc_Array( els, unsigned, nEls[0] + nEls[1] + nEls[2] );

	for( e_i = 0; e_i < nEls[2]; e_i++ ) {
		Grid_Lift( grid, e_i, dimInds );
		for( d_i = 0; d_i < grid->nDims; d_i++ )
			dimInds[d_i] += self->origin[d_i];
		els[nEls[0] + nEls[1] + e_i] = nGlobalEls[0] + nGlobalEls[1] + Grid_Project( globalGrid, dimInds );
	}

	MeshTopology_SetElements( topo, MT_FACE, nEls[0] + nEls[1] + nEls[2], els );

	FreeArray( dimInds );
	FreeArray( els );
	FreeObject( grid );
	FreeObject( globalGrid );

	MPI_Barrier( self->mpiComm );
	Journal_Printf( stream, "... done.\n" );
	Stream_UnIndent( stream );
  */
}

void _CartesianGenerator_GenElementVertexInc( void* meshGenerator, MeshTopology* topo, Grid*** grids ) {
	CartesianGenerator*	self = (CartesianGenerator*)meshGenerator;
	Stream*			stream = Journal_Register( Info_Type, self->type );
	int			nDomainEls;
	unsigned*		incEls;
	unsigned*		dimInds;
	unsigned		vertsPerEl;
	const Sync*		sync;
	unsigned		e_i;

	assert( self && Stg_CheckType( self, CartesianGenerator ) );
	assert( topo );
	assert( grids );

	Journal_Printf( stream, "Generating element-vertex incidence...\n" );
	Stream_Indent( stream );

	vertsPerEl = (topo->nDims == 1) ? 2 : (topo->nDims == 2) ? 4 : 8;

	sync = MeshTopology_GetDomain( topo, topo->nDims );
	nDomainEls = Sync_GetNumDomains( sync );
	incEls = Memory_Alloc_Array_Unnamed( unsigned, vertsPerEl );
	dimInds = Memory_Alloc_Array_Unnamed( unsigned, topo->nDims );
	for( e_i = 0; e_i < nDomainEls; e_i++ ) {
		unsigned	gInd = Sync_DomainToGlobal( sync, e_i );
		unsigned	curNode = 0;

		Grid_Lift( grids[topo->nDims][0], gInd, dimInds );

		incEls[curNode++] = Grid_Project( grids[0][0], dimInds );

		dimInds[0]++;
		incEls[curNode++] = Grid_Project( grids[0][0], dimInds );
		dimInds[0]--;

		if( topo->nDims >= 2 ) {
			dimInds[1]++;
			incEls[curNode++] = Grid_Project( grids[0][0], dimInds );

			dimInds[0]++;
			incEls[curNode++] = Grid_Project( grids[0][0], dimInds );
			dimInds[0]--;
			dimInds[1]--;

			if( topo->nDims >= 3 ) {
				dimInds[2]++;
				incEls[curNode++] = Grid_Project( grids[0][0], dimInds );

				dimInds[0]++;
				incEls[curNode++] = Grid_Project( grids[0][0], dimInds );
				dimInds[0]--;

				dimInds[1]++;
				incEls[curNode++] = Grid_Project( grids[0][0], dimInds );

				dimInds[0]++;
				incEls[curNode++] = Grid_Project( grids[0][0], dimInds );
				dimInds[0]--;
				dimInds[1]--;
				dimInds[2]--;
			}
		}

		CartesianGenerator_MapToDomain( self, (Sync*)MeshTopology_GetDomain( topo, 0 ), 
						vertsPerEl, incEls );
		MeshTopology_SetIncidence( topo, topo->nDims, e_i, MT_VERTEX, vertsPerEl, incEls );
	}

	FreeArray( incEls );
	FreeArray( dimInds );

	MPI_Barrier( self->mpiComm );
	Journal_Printf( stream, "... done.\n" );
	Stream_UnIndent( stream );
}

void _CartesianGenerator_GenVolumeEdgeInc( void* meshGenerator, MeshTopology* topo, Grid*** grids ) {
/*
	CartesianGenerator*	self = (CartesianGenerator*)meshGenerator;
	Stream*			stream = Journal_Register( Info_Type, self->type );
	unsigned*		nIncEls;
	unsigned**		incEls;
	unsigned*		dimInds;
	unsigned		e_i;

	assert( self && Stg_CheckType( self, CartesianGenerator ) );
	assert( topo && Stg_CheckType( topo, MeshTopology ) );
	assert( grids );
	assert( topo->nDims >= 3 );

	Journal_Printf( stream, "Generating volume-edge incidence...\n" );
	Stream_Indent( stream );

	nIncEls = Memory_Alloc_Array_Unnamed( unsigned, topo->domains[topo->nDims]->nDomains );
	incEls = Memory_Alloc_2DArray_Unnamed( unsigned, topo->domains[topo->nDims]->nDomains, 12 );
	dimInds = Memory_Alloc_Array_Unnamed( unsigned, topo->nDims );
	for( e_i = 0; e_i < topo->domains[MT_VOLUME]->nDomains; e_i++ ) {
		unsigned	gInd = Sync_DomainToGlobal( topo->domains[MT_VOLUME], e_i );

		nIncEls[e_i] = 12;
		Grid_Lift( grids[topo->nDims][0], gInd, dimInds );

		incEls[e_i][0] = Grid_Project( grids[1][0], dimInds );

		dimInds[1]++;
		incEls[e_i][1] = Grid_Project( grids[1][0], dimInds );
		dimInds[1]--;

		incEls[e_i][2] = Grid_Project( grids[1][1], dimInds ) + grids[1][0]->nPoints;

		dimInds[0]++;
		incEls[e_i][3] = Grid_Project( grids[1][1], dimInds ) + grids[1][0]->nPoints;
		dimInds[0]--;

		dimInds[2]++;
		incEls[e_i][4] = Grid_Project( grids[1][0], dimInds );

		dimInds[1]++;
		incEls[e_i][5] = Grid_Project( grids[1][0], dimInds );
		dimInds[1]--;

		incEls[e_i][6] = Grid_Project( grids[1][1], dimInds ) + grids[1][0]->nPoints;

		dimInds[0]++;
		incEls[e_i][7] = Grid_Project( grids[1][1], dimInds ) + grids[1][0]->nPoints;
		dimInds[0]--;
		dimInds[2]--;

		incEls[e_i][8] = Grid_Project( grids[1][2], dimInds ) + grids[1][0]->nPoints + grids[1][1]->nPoints;

		dimInds[0]++;
		incEls[e_i][9] = Grid_Project( grids[1][2], dimInds ) + grids[1][0]->nPoints + grids[1][1]->nPoints;
		dimInds[0]--;

		dimInds[1]++;
		incEls[e_i][10] = Grid_Project( grids[1][2], dimInds ) + grids[1][0]->nPoints + grids[1][1]->nPoints;

		dimInds[0]++;
		incEls[e_i][11] = Grid_Project( grids[1][2], dimInds ) + grids[1][0]->nPoints + grids[1][1]->nPoints;
		dimInds[0]--;
		dimInds[1]--;
	}

	CartesianGenerator_MapToDomain( self, topo->domains[MT_EDGE], topo->domains[MT_VOLUME]->nDomains, 
					nIncEls, incEls );
	MeshTopology_SetIncidence( topo, MT_VOLUME, MT_EDGE, nIncEls, incEls );
	FreeArray( nIncEls );
	FreeArray( incEls );
	FreeArray( dimInds );

	MPI_Barrier( self->mpiComm );
	Journal_Printf( stream, "... done.\n" );
	Stream_UnIndent( stream );
*/
}

void _CartesianGenerator_GenVolumeFaceInc( void* meshGenerator, MeshTopology* topo, Grid*** grids ) {
/*
	CartesianGenerator*	self = (CartesianGenerator*)meshGenerator;
	Stream*			stream = Journal_Register( Info_Type, self->type );
	unsigned*		nIncEls;
	unsigned**		incEls;
	unsigned*		dimInds;
	unsigned		e_i;

	assert( self && Stg_CheckType( self, CartesianGenerator ) );
	assert( topo && Stg_CheckType( topo, MeshTopology ) );
	assert( grids );
	assert( topo->nDims >= 3 );

	Journal_Printf( stream, "Generating volume-face incidence...\n" );
	Stream_Indent( stream );

	nIncEls = Memory_Alloc_Array_Unnamed( unsigned, topo->domains[topo->nDims]->nDomains );
	incEls = Memory_Alloc_2DArray_Unnamed( unsigned, topo->domains[topo->nDims]->nDomains, 6 );
	dimInds = Memory_Alloc_Array_Unnamed( unsigned, topo->nDims );
	for( e_i = 0; e_i < topo->domains[MT_VOLUME]->nDomains; e_i++ ) {
		unsigned	gInd = Sync_DomainToGlobal( topo->domains[MT_VOLUME], e_i );

		nIncEls[e_i] = 6;
		Grid_Lift( grids[topo->nDims][0], gInd, dimInds );

		incEls[e_i][0] = Grid_Project( grids[2][0], dimInds );

		dimInds[2]++;
		incEls[e_i][1] = Grid_Project( grids[2][0], dimInds );
		dimInds[2]--;

		incEls[e_i][2] = Grid_Project( grids[2][1], dimInds ) + grids[2][0]->nPoints;

		dimInds[1]++;
		incEls[e_i][3] = Grid_Project( grids[2][1], dimInds ) + grids[2][0]->nPoints;
		dimInds[1]--;

		incEls[e_i][4] = Grid_Project( grids[2][2], dimInds ) + grids[2][0]->nPoints + grids[2][1]->nPoints;

		dimInds[0]++;
		incEls[e_i][5] = Grid_Project( grids[2][2], dimInds ) + grids[2][0]->nPoints + grids[2][1]->nPoints;
		dimInds[0]--;
	}

	CartesianGenerator_MapToDomain( self, topo->domains[MT_FACE], topo->domains[MT_VOLUME]->nDomains, 
					nIncEls, incEls );
	MeshTopology_SetIncidence( topo, MT_VOLUME, MT_FACE, nIncEls, incEls );
	FreeArray( nIncEls );
	FreeArray( incEls );
	FreeArray( dimInds );

	MPI_Barrier( self->mpiComm );
	Journal_Printf( stream, "... done.\n" );
	Stream_UnIndent( stream );
*/
}

void _CartesianGenerator_GenFaceVertexInc( void* meshGenerator, MeshTopology* topo, Grid*** grids ) {
/*
	CartesianGenerator*	self = (CartesianGenerator*)meshGenerator;
	Stream*			stream = Journal_Register( Info_Type, self->type );
	unsigned*		nIncEls;
	unsigned**		incEls;
	unsigned*		dimInds;
	unsigned		e_i;

	assert( self && Stg_CheckType( self, CartesianGenerator ) );
	assert( topo && Stg_CheckType( topo, MeshTopology ) );
	assert( grids );
	assert( topo->nDims >= 2 );

	Journal_Printf( stream, "Generating face-vertex incidence...\n" );
	Stream_Indent( stream );

	nIncEls = Memory_Alloc_Array_Unnamed( unsigned, topo->domains[MT_FACE]->nDomains );
	incEls = Memory_Alloc_2DArray_Unnamed( unsigned, topo->domains[MT_FACE]->nDomains, 4 );
	dimInds = Memory_Alloc_Array_Unnamed( unsigned, topo->nDims );
	for( e_i = 0; e_i < topo->domains[MT_FACE]->nDomains; e_i++ ) {
		unsigned	gInd = Sync_DomainToGlobal( topo->domains[MT_FACE], e_i );

		nIncEls[e_i] = 4;

		if( gInd < grids[2][0]->nPoints ) {
			Grid_Lift( grids[2][0], gInd, dimInds );

			incEls[e_i][0] = Grid_Project( grids[0][0], dimInds );

			dimInds[0]++;
			incEls[e_i][1] = Grid_Project( grids[0][0], dimInds );
			dimInds[0]--;

			dimInds[1]++;
			incEls[e_i][2] = Grid_Project( grids[0][0], dimInds );

			dimInds[0]++;
			incEls[e_i][3] = Grid_Project( grids[0][0], dimInds );
			dimInds[0]--;
			dimInds[1]--;
		}
		else if( gInd < grids[2][0]->nPoints + grids[2][1]->nPoints ) {
			assert( topo->nDims >= 3 );

			Grid_Lift( grids[2][1], gInd - grids[2][0]->nPoints, dimInds );

			incEls[e_i][0] = Grid_Project( grids[0][0], dimInds );

			dimInds[0]++;
			incEls[e_i][1] = Grid_Project( grids[0][0], dimInds );
			dimInds[0]--;

			dimInds[2]++;
			incEls[e_i][2] = Grid_Project( grids[0][0], dimInds );

			dimInds[0]++;
			incEls[e_i][3] = Grid_Project( grids[0][0], dimInds );
			dimInds[0]--;
			dimInds[2]--;
		}
		else {
			assert( gInd < grids[2][0]->nPoints + grids[2][1]->nPoints + grids[2][2]->nPoints );
			assert( topo->nDims >= 3 );

			Grid_Lift( grids[2][2], gInd - grids[2][0]->nPoints - grids[2][1]->nPoints, dimInds );

			incEls[e_i][0] = Grid_Project( grids[0][0], dimInds );

			dimInds[1]++;
			incEls[e_i][1] = Grid_Project( grids[0][0], dimInds );
			dimInds[1]--;

			dimInds[2]++;
			incEls[e_i][2] = Grid_Project( grids[0][0], dimInds );

			dimInds[1]++;
			incEls[e_i][3] = Grid_Project( grids[0][0], dimInds );
			dimInds[1]--;
			dimInds[2]--;
		}
	}

	CartesianGenerator_MapToDomain( self, topo->domains[MT_VERTEX], topo->domains[MT_FACE]->nDomains, 
					nIncEls, incEls );
	MeshTopology_SetIncidence( topo, MT_FACE, MT_VERTEX, nIncEls, incEls );
	FreeArray( nIncEls );
	FreeArray( incEls );
	FreeArray( dimInds );

	MPI_Barrier( self->mpiComm );
	Journal_Printf( stream, "... done.\n" );
	Stream_UnIndent( stream );
*/
}

void _CartesianGenerator_GenFaceEdgeInc( void* meshGenerator, MeshTopology* topo, Grid*** grids ) {
/*
	CartesianGenerator*	self = (CartesianGenerator*)meshGenerator;
	Stream*			stream = Journal_Register( Info_Type, self->type );
	unsigned*		nIncEls;
	unsigned**		incEls;
	unsigned*		dimInds;
	unsigned		e_i;

	assert( self && Stg_CheckType( self, CartesianGenerator ) );
	assert( topo && Stg_CheckType( topo, MeshTopology ) );
	assert( grids );
	assert( topo->nDims >= 2 );

	Journal_Printf( stream, "Generating face-edge incidence...\n" );
	Stream_Indent( stream );

	nIncEls = Memory_Alloc_Array_Unnamed( unsigned, topo->domains[MT_FACE]->nDomains );
	incEls = Memory_Alloc_2DArray_Unnamed( unsigned, topo->domains[MT_FACE]->nDomains, 4 );
	dimInds = Memory_Alloc_Array_Unnamed( unsigned, topo->nDims );
	for( e_i = 0; e_i < topo->domains[MT_FACE]->nDomains; e_i++ ) {
		unsigned	gInd = Sync_DomainToGlobal( topo->domains[MT_FACE], e_i );

		nIncEls[e_i] = 4;

		if( gInd < grids[2][0]->nPoints ) {
			Grid_Lift( grids[2][0], gInd, dimInds );

			incEls[e_i][0] = Grid_Project( grids[1][0], dimInds );

			dimInds[1]++;
			incEls[e_i][1] = Grid_Project( grids[1][0], dimInds );
			dimInds[1]--;

			incEls[e_i][2] = Grid_Project( grids[1][1], dimInds ) + grids[1][0]->nPoints;

			dimInds[0]++;
			incEls[e_i][3] = Grid_Project( grids[1][1], dimInds ) + grids[1][0]->nPoints;
			dimInds[0]--;
		}
		else if( gInd < grids[2][0]->nPoints + grids[2][1]->nPoints ) {
			assert( topo->nDims >= 3 );

			Grid_Lift( grids[2][1], gInd - grids[2][0]->nPoints, dimInds );

			incEls[e_i][0] = Grid_Project( grids[1][0], dimInds );

			dimInds[2]++;
			incEls[e_i][1] = Grid_Project( grids[1][0], dimInds );
			dimInds[2]--;

			incEls[e_i][2] = Grid_Project( grids[1][2], dimInds ) + grids[1][0]->nPoints + grids[1][1]->nPoints;

			dimInds[0]++;
			incEls[e_i][3] = Grid_Project( grids[1][2], dimInds ) + grids[1][0]->nPoints + grids[1][1]->nPoints;
			dimInds[0]--;
		}
		else {
			assert( gInd < grids[2][0]->nPoints + grids[2][1]->nPoints + grids[2][2]->nPoints );
			assert( topo->nDims >= 3 );

			Grid_Lift( grids[2][2], gInd - grids[2][0]->nPoints - grids[2][1]->nPoints, dimInds );

			incEls[e_i][0] = Grid_Project( grids[1][1], dimInds ) + grids[1][0]->nPoints;

			dimInds[2]++;
			incEls[e_i][1] = Grid_Project( grids[1][1], dimInds ) + grids[1][0]->nPoints;
			dimInds[2]--;

			incEls[e_i][2] = Grid_Project( grids[1][2], dimInds ) + grids[1][0]->nPoints + grids[1][1]->nPoints;

			dimInds[1]++;
			incEls[e_i][3] = Grid_Project( grids[1][2], dimInds ) + grids[1][0]->nPoints + grids[1][1]->nPoints;
			dimInds[1]--;
		}
	}

	CartesianGenerator_MapToDomain( self, topo->domains[MT_EDGE], topo->domains[MT_FACE]->nDomains, 
					nIncEls, incEls );
	MeshTopology_SetIncidence( topo, MT_FACE, MT_EDGE, nIncEls, incEls );
	FreeArray( nIncEls );
	FreeArray( incEls );
	FreeArray( dimInds );

	MPI_Barrier( self->mpiComm );
	Journal_Printf( stream, "... done.\n" );
	Stream_UnIndent( stream );
*/
}

void _CartesianGenerator_GenEdgeVertexInc( void* meshGenerator, MeshTopology* topo, Grid*** grids ) {
/*
	CartesianGenerator*	self = (CartesianGenerator*)meshGenerator;
	Stream*			stream = Journal_Register( Info_Type, self->type );
	unsigned*		nIncEls;
	unsigned**		incEls;
	unsigned*		dimInds;
	unsigned		e_i;

	assert( self && Stg_CheckType( self, CartesianGenerator ) );
	assert( topo && Stg_CheckType( topo, MeshTopology ) );
	assert( grids );

	Journal_Printf( stream, "Generating edge-vertex incidence...\n" );
	Stream_Indent( stream );

	nIncEls = Memory_Alloc_Array_Unnamed( unsigned, topo->domains[MT_EDGE]->nDomains );
	incEls = Memory_Alloc_2DArray_Unnamed( unsigned, topo->domains[MT_EDGE]->nDomains, 2 );
	dimInds = Memory_Alloc_Array_Unnamed( unsigned, topo->nDims );
	for( e_i = 0; e_i < topo->domains[MT_EDGE]->nDomains; e_i++ ) {
		unsigned	gInd = Sync_DomainToGlobal( topo->domains[MT_EDGE], e_i );

		nIncEls[e_i] = 2;

		if( gInd < grids[1][0]->nPoints ) {
			Grid_Lift( grids[1][0], gInd, dimInds );

			incEls[e_i][0] = Grid_Project( grids[0][0], dimInds );

			dimInds[0]++;
			incEls[e_i][1] = Grid_Project( grids[0][0], dimInds );
			dimInds[0]--;
		}
		else if( gInd < grids[1][0]->nPoints + grids[1][1]->nPoints ) {
			assert( topo->nDims >= 2 );

			Grid_Lift( grids[1][1], gInd - grids[1][0]->nPoints, dimInds );

			incEls[e_i][0] = Grid_Project( grids[0][0], dimInds );

			dimInds[1]++;
			incEls[e_i][1] = Grid_Project( grids[0][0], dimInds );
			dimInds[1]--;
		}
		else {
			assert( gInd < grids[1][0]->nPoints + grids[1][1]->nPoints + grids[1][2]->nPoints );
			assert( topo->nDims >= 3 );

			Grid_Lift( grids[1][2], gInd - grids[1][0]->nPoints - grids[1][1]->nPoints, dimInds );

			incEls[e_i][0] = Grid_Project( grids[0][0], dimInds );

			dimInds[2]++;
			incEls[e_i][1] = Grid_Project( grids[0][0], dimInds );
			dimInds[2]--;
		}
	}

	CartesianGenerator_MapToDomain( self, topo->domains[MT_VERTEX], topo->domains[MT_EDGE]->nDomains, 
					nIncEls, incEls );
	MeshTopology_SetIncidence( topo, MT_EDGE, MT_VERTEX, nIncEls, incEls );
	FreeArray( nIncEls );
	FreeArray( incEls );
	FreeArray( dimInds );

	MPI_Barrier( self->mpiComm );
	Journal_Printf( stream, "... done.\n" );
	Stream_UnIndent( stream );
*/
}

void _CartesianGenerator_GenElementTypes( void* meshGenerator, Mesh* mesh ) {
	CartesianGenerator*	self = (CartesianGenerator*)meshGenerator;
	Stream*			stream;
	unsigned		nDomainEls;
	unsigned		e_i;

	assert( self && Stg_CheckType( self, CartesianGenerator ) );

	stream = Journal_Register( Info_Type, self->type );
	Journal_Printf( stream, "Generating element types...\n" );
	Stream_Indent( stream );

	mesh->nElTypes = 1;
	mesh->elTypes = Memory_Alloc_Array( Mesh_ElementType*, mesh->nElTypes, "Mesh::elTypes" );
	mesh->elTypes[0] = (Mesh_ElementType*)Mesh_HexType_New();
	Mesh_ElementType_SetMesh( mesh->elTypes[0], mesh );
	nDomainEls = Mesh_GetDomainSize( mesh, Mesh_GetDimSize( mesh ) );
	mesh->elTypeMap = Memory_Alloc_Array( unsigned, nDomainEls, "Mesh::elTypeMap" );
	for( e_i = 0; e_i < nDomainEls; e_i++ )
		mesh->elTypeMap[e_i] = 0;

	if( self->regular )
		Mesh_SetAlgorithms( mesh, Mesh_RegularAlgorithms_New( "" ) );

	MPI_Barrier( self->mpiComm );
	Journal_Printf( stream, "... element types are '%s',\n", mesh->elTypes[0]->type );
	Journal_Printf( stream, "... mesh algorithm type is '%s',\n", mesh->algorithms->type );
	Journal_Printf( stream, "... done.\n" );
	Stream_UnIndent( stream );
}


/*--------------------------------------------------------------------------------------------------------------------------
** Public Functions
*/

void CartesianGenerator_SetGeometryParams( void* meshGenerator, double* min, double* max ) {
	CartesianGenerator*	self = (CartesianGenerator*)meshGenerator;

	/* Sanity check. */
	assert( self );
	assert( !self->elGrid->nDims || (min && max) );

	/* Kill everything we have, geometrically. */
	KillArray( self->crdMin );
	KillArray( self->crdMax );

	/* Set the parameters. */
	if( self->elGrid->nDims ) {
		self->crdMin = Memory_Alloc_Array( double, self->elGrid->nDims, "CartesianGenerator::min" );
		self->crdMax = Memory_Alloc_Array( double, self->elGrid->nDims, "CartesianGenerator::max" );
		memcpy( self->crdMin, min, self->elGrid->nDims * sizeof(double) );
		memcpy( self->crdMax, max, self->elGrid->nDims * sizeof(double) );
	}
}

void CartesianGenerator_SetShadowDepth( void* meshGenerator, unsigned depth ) {
	CartesianGenerator*	self = (CartesianGenerator*)meshGenerator;

	/* Sanity check. */
	assert( self );

	self->shadowDepth = depth;
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Private Functions
*/

void CartesianGenerator_BuildDecomp( CartesianGenerator* self ) {
	unsigned	nProcs, rank;
	unsigned	nPos = 0;
	unsigned**	posNSubDomains = NULL;
	unsigned*	tmpSubDomains;
	double		bestRatio;
	unsigned	bestPos;
	unsigned	*myRankInds, *rankInds;
	unsigned	nNbrs, *nbrs;
	unsigned	p_i, d_i, r_i;

	/* Sanity check. */
	assert( self );

	/* Collect information regarding parallel configuration. */
	MPI_Comm_size( self->mpiComm, (int*)&nProcs );
	MPI_Comm_rank( self->mpiComm, (int*)&rank );

	/* Allocate for possible sub-domains. */
	tmpSubDomains = Memory_Alloc_Array( unsigned, self->elGrid->nDims, "" );

	/* Build a list of all acceptable decompositions. */
	CartesianGenerator_RecurseDecomps( self, 0, nProcs, 
					   tmpSubDomains, &nPos, &posNSubDomains );
	assert( nPos );

	/* Free unneeded memory. */
	FreeArray( tmpSubDomains );

	/* Now check for the best ratio. */
	bestRatio = HUGE_VAL;
	bestPos = -1;
	for( p_i = 0; p_i < nPos; p_i++ ) {
		double		curRatio = 0.0;
		unsigned	decompDims = 0;

		/* If decomposed in more dimensions than allowed, skip. */
		for( d_i = 0; d_i < self->elGrid->nDims; d_i++ )
			decompDims += (posNSubDomains[p_i][d_i] > 1) ? 1 : 0;
		if( self->maxDecompDims && decompDims > self->maxDecompDims )
			continue;

		/* Evaluate ratios. */
		for( d_i = 0; d_i < self->elGrid->nDims; d_i++ ) {
			double		nI = (double)self->elGrid->sizes[d_i] / (double)posNSubDomains[p_i][d_i];
			unsigned	d_j;

			for( d_j = d_i + 1; d_j < self->elGrid->nDims; d_j++ ) {
				double	nJ = (double)self->elGrid->sizes[d_j] / (double)posNSubDomains[p_i][d_j];

				curRatio += (nI > nJ) ? nI / nJ : nJ / nI;
			}
		}

		/* Better than best? */
		if( curRatio < bestRatio ) {
			bestRatio = curRatio;
			bestPos = p_i;
		}
	}
	assert( bestPos != -1 );

	/* Allocate for results. */
	self->origin = Memory_Alloc_Array( unsigned, self->elGrid->nDims, "CartesianGenerator::origin" );
	self->range = Memory_Alloc_Array( unsigned, self->elGrid->nDims, "CartesianGenerator::range" );

	/* Build a sub-domain grid. */
	self->procGrid = Grid_New();
	Grid_SetNumDims( self->procGrid, self->elGrid->nDims );
	Grid_SetSizes( self->procGrid, posNSubDomains[bestPos] );

	/* Free unneeded space. */
	FreeArray( posNSubDomains );

	/* Lift the rank to a parameterised offset. */
	Grid_Lift( self->procGrid, rank, self->origin );
	for( d_i = 0; d_i < self->elGrid->nDims; d_i++ ) {
		unsigned	base = self->elGrid->sizes[d_i] / self->procGrid->sizes[d_i];
		unsigned	mod = self->elGrid->sizes[d_i] % self->procGrid->sizes[d_i];
		unsigned	origin = self->origin[d_i];

		self->origin[d_i] *= base;
		self->range[d_i] = base;
		if( origin < mod ) {
			self->origin[d_i] += origin;
			self->range[d_i]++;
		}
		else
			self->origin[d_i] += mod;
	}

	/* Copy to the vertex equivalents. */
	self->vertOrigin = AllocArray( unsigned, self->elGrid->nDims );
	self->vertRange = AllocArray( unsigned, self->elGrid->nDims );
	for( d_i = 0; d_i < self->elGrid->nDims; d_i++ ) {
		self->vertOrigin[d_i] = self->origin[d_i];
		self->vertRange[d_i] = self->range[d_i] + 1;
	}

	/* Build the comm topology. */
	myRankInds = AllocArray( unsigned, Grid_GetNumDims( self->procGrid ) );
	rankInds = AllocArray( unsigned, Grid_GetNumDims( self->procGrid ) );
	Grid_Lift( self->procGrid, rank, myRankInds );
	nNbrs = 0;
	nbrs = NULL;
	for( r_i = 0; r_i < Grid_GetNumPoints( self->procGrid ); r_i++ ) {
		if( r_i == rank )
			continue;

		Grid_Lift( self->procGrid, r_i, rankInds );
		for( d_i = 0; d_i < Grid_GetNumDims( self->procGrid ); d_i++ ) {
			if( (myRankInds[d_i] > 0 && rankInds[d_i] < myRankInds[d_i] - 1) || 
			    (myRankInds[d_i] < self->procGrid->sizes[d_i] - 1 && rankInds[d_i] > myRankInds[d_i] + 1) )
			{
				break;
			}
		}
		if( d_i == Grid_GetNumDims( self->procGrid ) ) {
			nbrs = ReallocArray( nbrs, unsigned, nNbrs + 1 );
			nbrs[nNbrs++] = r_i;
		}
	}

	FreeArray( myRankInds );
	FreeArray( rankInds );

	self->comm = Comm_New();
	NewClass_AddRef( self->comm );
	Comm_SetMPIComm( self->comm, self->mpiComm );
	Comm_SetNeighbours( self->comm, nNbrs, nbrs );
	FreeArray( nbrs );
}

void CartesianGenerator_RecurseDecomps( CartesianGenerator* self, 
					unsigned dim, unsigned max, 
					unsigned* nSubDomains, 
					unsigned* nPos, unsigned*** posNSubDomains )
{
	unsigned	nProcs;
	unsigned	nSDs = 1;
	unsigned	d_i;

	MPI_Comm_size( self->mpiComm, (int*)&nProcs );

	/* If we're over the limit, return immediately. */
	for( d_i = 0; d_i < dim; d_i++ )
		nSDs *= nSubDomains[d_i];
	if( nSDs > nProcs )
		return;

	/* Where are we up to? */
	if( dim == self->elGrid->nDims ) {
		/* If this covers all processors, store it. */
		if( nSDs == nProcs ) {
			/* If we havn't already allocated do it now. */
			if( !*posNSubDomains ) {
				*nPos = 1;
				*posNSubDomains = Memory_Alloc_2DArray_Unnamed( unsigned, 2, self->elGrid->nDims );
			}
			else {
				/* Reallocate the arrays. */
				(*nPos)++;
				if( *nPos != 2 ) {
					*posNSubDomains = Memory_Realloc_2DArray( *posNSubDomains, unsigned, 
										  *nPos, self->elGrid->nDims );
				}
			}

			/* Store status. */
			memcpy( (*posNSubDomains)[(*nPos) - 1], nSubDomains, self->elGrid->nDims * sizeof(unsigned) );
		}
	}
	else {
		unsigned	p_i;

		/* Loop over all remaining */
		for( p_i = 0; p_i < max; p_i++ ) {
			/* Don't try and decompose more than this dimension allows. */
			if( p_i >= self->elGrid->sizes[dim] || 
			    (self->maxDecomp[dim] && p_i >= self->maxDecomp[dim]) )
			{
				break;
			}

			/* If we have a minimum decomp, skip until we reach it. */
			if( self->minDecomp[dim] && p_i < self->minDecomp[dim] - 1 )
				continue;

			/* Set the number of sub-domains. */
			nSubDomains[dim] = p_i + 1;

			/* Try this combination. */
			CartesianGenerator_RecurseDecomps( self, dim + 1, max - nSDs + 1, nSubDomains, 
							   nPos, posNSubDomains );
		}
	}
}

void CartesianGenerator_GenTopo( CartesianGenerator* self, MeshTopology* topo ) {
	Grid***		grids;
	Comm*		comm;
	unsigned	d_i, d_j;

	assert( self );
	assert( topo );

	MeshTopology_SetNumDims( topo, self->elGrid->nDims );

	/* Build additional grids for use in numbering. */
	grids = Memory_Alloc_2DArray_Unnamed( Grid*, topo->nTDims, topo->nTDims );
	for( d_i = 0; d_i < topo->nTDims; d_i++ )
		memset( grids[d_i], 0, topo->nTDims * sizeof(Grid*) );

	grids[topo->nDims][0] = self->elGrid;
	grids[0][0] = self->vertGrid;

	grids[1][0] = Grid_New();
	Grid_SetNumDims( grids[1][0], topo->nDims );
	for( d_i = 0; d_i < topo->nDims; d_i++ ) {
		if( d_i == 0 ) continue;
		grids[topo->nDims][0]->sizes[d_i]++;
	}
	Grid_SetSizes( grids[1][0], grids[topo->nDims][0]->sizes );
	for( d_i = 0; d_i < topo->nDims; d_i++ ) {
		if( d_i == 0 ) continue;
		grids[topo->nDims][0]->sizes[d_i]--;
	}

	if( topo->nDims >= 2 ) {
		grids[1][1] = Grid_New();
		Grid_SetNumDims( grids[1][1], topo->nDims );
		for( d_i = 0; d_i < topo->nDims; d_i++ ) {
			if( d_i == 1 ) continue;
			grids[topo->nDims][0]->sizes[d_i]++;
		}
		Grid_SetSizes( grids[1][1], grids[topo->nDims][0]->sizes );
		for( d_i = 0; d_i < topo->nDims; d_i++ ) {
			if( d_i == 1 ) continue;
			grids[topo->nDims][0]->sizes[d_i]--;
		}

		if( topo->nDims >= 3 ) {
			grids[1][2] = Grid_New();
			Grid_SetNumDims( grids[1][2], topo->nDims );
			for( d_i = 0; d_i < topo->nDims; d_i++ ) {
				if( d_i == 2 ) continue;
				grids[topo->nDims][0]->sizes[d_i]++;
			}
			Grid_SetSizes( grids[1][2], grids[topo->nDims][0]->sizes );
			for( d_i = 0; d_i < topo->nDims; d_i++ ) {
				if( d_i == 2 ) continue;
				grids[topo->nDims][0]->sizes[d_i]--;
			}

			grids[2][0] = Grid_New();
			Grid_SetNumDims( grids[2][0], topo->nDims );
			grids[topo->nDims][0]->sizes[2]++;
			Grid_SetSizes( grids[2][0], grids[topo->nDims][0]->sizes );
			grids[topo->nDims][0]->sizes[2]--;

			grids[2][1] = Grid_New();
			Grid_SetNumDims( grids[2][1], topo->nDims );
			grids[topo->nDims][0]->sizes[1]++;
			Grid_SetSizes( grids[2][1], grids[topo->nDims][0]->sizes );
			grids[topo->nDims][0]->sizes[1]--;

			grids[2][2] = Grid_New();
			Grid_SetNumDims( grids[2][2], topo->nDims );
			grids[topo->nDims][0]->sizes[0]++;
			Grid_SetSizes( grids[2][2], grids[topo->nDims][0]->sizes );
			grids[topo->nDims][0]->sizes[0]--;
		}
	}

	/* Generate topological elements. */
	if( self->enabledDims[0] )
		CartesianGenerator_GenVertices( self, topo, grids );
	if( self->enabledDims[self->nDims] )
		CartesianGenerator_GenElements( self, topo, grids );
	if( topo->nDims >= 2 ) {
		if( self->enabledDims[1] )
			CartesianGenerator_GenEdges( self, topo, grids );
		if( topo->nDims >= 3 ) {
			if( self->enabledDims[2] )
				CartesianGenerator_GenFaces( self, topo, grids );
		}
	}

	/* Generate topological incidence. */
	if( self->enabledInc[self->nDims][0] )
		CartesianGenerator_GenElementVertexInc( self, topo, grids );
	if( topo->nDims >= 2 ) {
		if( self->enabledInc[2][1] )
			CartesianGenerator_GenFaceEdgeInc( self, topo, grids );
		if( self->enabledInc[1][0] )
			CartesianGenerator_GenEdgeVertexInc( self, topo, grids );
		if( topo->nDims >= 3 ) {
			if( self->enabledInc[2][0] )
				CartesianGenerator_GenFaceVertexInc( self, topo, grids );
			if( self->enabledInc[3][1] )
				CartesianGenerator_GenVolumeEdgeInc( self, topo, grids );
			if( self->enabledInc[3][2] )
				CartesianGenerator_GenVolumeFaceInc( self, topo, grids );
		}
	}

	/* Set the shadow depth and correct incidence. */
	comm = MeshTopology_GetComm( topo );
	if( self->shadowDepth && Comm_GetNumNeighbours( comm ) > 0 ) {
		/* Build enough incidence to set shadow depth. */
		MeshTopology_InvertIncidence( topo, MT_VERTEX, topo->nDims );
		MeshTopology_ExpandIncidence( topo, topo->nDims );

		/* Set the shadow depth. */
		MeshTopology_SetShadowDepth( topo, self->shadowDepth );

		/* Kill up relations and neighbours. */
		MeshTopology_RemoveIncidence( topo, topo->nDims, topo->nDims );
		MeshTopology_RemoveIncidence( topo, 0, topo->nDims );
	}

	/* Complete all required relations. */
	for( d_i = 0; d_i < self->nDims; d_i++ ) {
		for( d_j = d_i + 1; d_j <= self->nDims; d_j++ ) {
			if( !self->enabledInc[d_i][d_j] )
				continue;

			MeshTopology_InvertIncidence( topo, d_i, d_j );
		}
	}
	for( d_i = 0; d_i <= self->nDims; d_i++ ) {
		if( self->enabledInc[d_i][d_i] ) {
			if( d_i == 0 )
				CartesianGenerator_CompleteVertexNeighbours( self, topo, grids );
			else {
				MeshTopology_ExpandIncidence( topo, d_i );
			}
		}
	}

	/* Free allocated grids. */
	grids[topo->nDims][0] = NULL;
	for( d_i = 1; d_i < topo->nDims; d_i++ ) {
		unsigned	d_j;

		for( d_j = 0; d_j < topo->nTDims; d_j++ )
			FreeObject( grids[d_i][d_j] );
	}
	FreeArray( grids );
}

void CartesianGenerator_GenEdges2D( CartesianGenerator* self, MeshTopology* topo, Grid*** grids ) {
/*
	Grid*		globalGrid;
	unsigned	nGlobalEls;
	Grid*		grid;
	unsigned	nEls[2];
	unsigned*	els;
	unsigned*	dimInds;
	unsigned	d_i, e_i;

	assert( self );
	assert( topo );
	assert( grids );
	assert( self->elGrid->nDims == 2 );

	globalGrid = Grid_New();
	Grid_SetNumDims( globalGrid, self->elGrid->nDims );
	self->elGrid->sizes[1]++;
	Grid_SetSizes( globalGrid, self->elGrid->sizes );
	self->elGrid->sizes[1]--;

	grid = Grid_New();
	Grid_SetNumDims( grid, self->elGrid->nDims );
	self->range[1]++;
	Grid_SetSizes( grid, self->range );
	self->range[1]--;

	nEls[0] = grid->sizes[0];
	nGlobalEls = globalGrid->sizes[0];
	for( d_i = 1; d_i < grid->nDims; d_i++ ) {
		nEls[0] *= grid->sizes[d_i];
		nGlobalEls *= globalGrid->sizes[d_i];
	}
	els = Memory_Alloc_Array_Unnamed( unsigned, nEls[0] );

	dimInds = Memory_Alloc_Array_Unnamed( unsigned, self->elGrid->nDims );

	for( e_i = 0; e_i < nEls[0]; e_i++ ) {
		Grid_Lift( grid, e_i, dimInds );
		for( d_i = 0; d_i < grid->nDims; d_i++ )
			dimInds[d_i] += self->origin[d_i];
		els[e_i] = Grid_Project( globalGrid, dimInds );
	}

	self->elGrid->sizes[0]++;
	Grid_SetSizes( globalGrid, self->elGrid->sizes );
	self->elGrid->sizes[0]--;

	self->range[0]++;
	Grid_SetSizes( grid, self->range );
	self->range[0]--;

	nEls[1] = grid->sizes[0];
	for( d_i = 1; d_i < grid->nDims; d_i++ )
		nEls[1] *= grid->sizes[d_i];
	els = Memory_Realloc_Array( els, unsigned, nEls[0] + nEls[1] );

	for( e_i = 0; e_i < nEls[1]; e_i++ ) {
		Grid_Lift( grid, e_i, dimInds );
		for( d_i = 0; d_i < grid->nDims; d_i++ )
			dimInds[d_i] += self->origin[d_i];
		els[nEls[0] + e_i] = nGlobalEls + Grid_Project( globalGrid, dimInds );
	}

	MeshTopology_SetElements( topo, MT_EDGE, nEls[0] + nEls[1], els );

	FreeArray( dimInds );
	FreeArray( els );
	FreeObject( grid );
	FreeObject( globalGrid );
*/
}

void CartesianGenerator_GenEdges3D( CartesianGenerator* self, MeshTopology* topo, Grid*** grids ) {
/*
	Grid*		globalGrid;
	unsigned	nGlobalEls[2];
	Grid*		grid;
	unsigned	nEls[3];
	unsigned*	els;
	unsigned*	dimInds;
	unsigned	d_i, e_i;

	assert( self );
	assert( topo );
	assert( grids );
	assert( self->elGrid->nDims == 3 );

	globalGrid = Grid_New();
	Grid_SetNumDims( globalGrid, self->elGrid->nDims );
	self->elGrid->sizes[1]++;
	self->elGrid->sizes[2]++;
	Grid_SetSizes( globalGrid, self->elGrid->sizes );
	self->elGrid->sizes[1]--;
	self->elGrid->sizes[2]--;

	grid = Grid_New();
	Grid_SetNumDims( grid, self->elGrid->nDims );
	self->range[1]++;
	self->range[2]++;
	Grid_SetSizes( grid, self->range );
	self->range[1]--;
	self->range[2]--;

	nEls[0] = grid->sizes[0];
	nGlobalEls[0] = globalGrid->sizes[0];
	for( d_i = 1; d_i < grid->nDims; d_i++ ) {
		nEls[0] *= grid->sizes[d_i];
		nGlobalEls[0] *= globalGrid->sizes[d_i];
	}
	els = Memory_Alloc_Array_Unnamed( unsigned, nEls[0] );

	dimInds = Memory_Alloc_Array_Unnamed( unsigned, self->elGrid->nDims );

	for( e_i = 0; e_i < nEls[0]; e_i++ ) {
		Grid_Lift( grid, e_i, dimInds );
		for( d_i = 0; d_i < grid->nDims; d_i++ )
			dimInds[d_i] += self->origin[d_i];
		els[e_i] = Grid_Project( globalGrid, dimInds );
	}

	self->elGrid->sizes[0]++;
	self->elGrid->sizes[2]++;
	Grid_SetSizes( globalGrid, self->elGrid->sizes );
	self->elGrid->sizes[0]--;
	self->elGrid->sizes[2]--;

	self->range[0]++;
	self->range[2]++;
	Grid_SetSizes( grid, self->range );
	self->range[0]--;
	self->range[2]--;

	nEls[1] = grid->sizes[0];
	nGlobalEls[1] = globalGrid->sizes[0];
	for( d_i = 1; d_i < grid->nDims; d_i++ ) {
		nEls[1] *= grid->sizes[d_i];
		nGlobalEls[1] *= globalGrid->sizes[d_i];
	}
	els = Memory_Realloc_Array( els, unsigned, nEls[0] + nEls[1] );

	for( e_i = 0; e_i < nEls[1]; e_i++ ) {
		Grid_Lift( grid, e_i, dimInds );
		for( d_i = 0; d_i < grid->nDims; d_i++ )
			dimInds[d_i] += self->origin[d_i];
		els[nEls[0] + e_i] = nGlobalEls[0] + Grid_Project( globalGrid, dimInds );
	}

	self->elGrid->sizes[0]++;
	self->elGrid->sizes[1]++;
	Grid_SetSizes( globalGrid, self->elGrid->sizes );
	self->elGrid->sizes[0]--;
	self->elGrid->sizes[1]--;

	self->range[0]++;
	self->range[1]++;
	Grid_SetSizes( grid, self->range );
	self->range[0]--;
	self->range[1]--;

	nEls[2] = grid->sizes[0];
	for( d_i = 1; d_i < grid->nDims; d_i++ )
		nEls[2] *= grid->sizes[d_i];
	els = Memory_Realloc_Array( els, unsigned, nEls[0] + nEls[1] + nEls[2] );

	for( e_i = 0; e_i < nEls[2]; e_i++ ) {
		Grid_Lift( grid, e_i, dimInds );
		for( d_i = 0; d_i < grid->nDims; d_i++ )
			dimInds[d_i] += self->origin[d_i];
		els[nEls[0] + nEls[1] + e_i] = nGlobalEls[0] + nGlobalEls[1] + Grid_Project( globalGrid, dimInds );
	}

	MeshTopology_SetElements( topo, MT_EDGE, nEls[0] + nEls[1] + nEls[2], els );

	FreeArray( dimInds );
	FreeArray( els );
	FreeObject( grid );
	FreeObject( globalGrid );
*/
}

void CartesianGenerator_CompleteVertexNeighbours( CartesianGenerator* self, MeshTopology* topo, Grid*** grids ) {
	Stream*		stream = Journal_Register( Info_Type, self->type );
	Sync*		sync;
	unsigned	nDims;
	unsigned	nVerts;
	unsigned*	inds;
	unsigned	*nNbrs, **nbrs;
	unsigned	domain, global;
	unsigned	v_i;

	assert( self );
	assert( topo );
	assert( grids );

	Journal_Printf( stream, "Generating vertex neighbours...\n" );
	Stream_Indent( stream );

	nDims = topo->nDims;
	sync = MeshTopology_GetDomain( topo, MT_VERTEX );
	nVerts = Sync_GetNumDomains( sync );
	inds = AllocArray( unsigned, nDims );
	nNbrs = AllocArray( unsigned, nVerts );
	nbrs = AllocArray2D( unsigned, nVerts, (nDims == 3) ? 6 : (nDims == 2) ? 4 : 2 );
	for( v_i = 0; v_i < nVerts; v_i++ ) {
		nNbrs[v_i] = 0;
		global = Sync_DomainToGlobal( sync, v_i );
		Grid_Lift( grids[0][0], global, inds );

		if( inds[0] > 0 ) {
			inds[0]--;
			domain = Grid_Project( grids[0][0], inds );
			if( Sync_TryGlobalToDomain( sync, domain, &domain ) ) {
				nbrs[v_i][nNbrs[v_i]] = domain;
				nNbrs[v_i]++;
			}
			inds[0]++;
		}

		if( inds[0] < grids[0][0]->sizes[0] - 1 ) {
			inds[0]++;
			domain = Grid_Project( grids[0][0], inds );
			if( Sync_TryGlobalToDomain( sync, domain, &domain ) ) {
				nbrs[v_i][nNbrs[v_i]] = domain;
				nNbrs[v_i]++;
			}
			inds[0]--;
		}

		if( nDims >= 2 ) {
			if( inds[1] > 0 ) {
				inds[1]--;
				domain = Grid_Project( grids[0][0], inds );
				if( Sync_TryGlobalToDomain( sync, domain, &domain ) ) {
					nbrs[v_i][nNbrs[v_i]] = domain;
					nNbrs[v_i]++;
				}
				inds[1]++;
			}

			if( inds[1] < grids[0][0]->sizes[1] - 1 ) {
				inds[1]++;
				domain = Grid_Project( grids[0][0], inds );
				if( Sync_TryGlobalToDomain( sync, domain, &domain ) ) {
					nbrs[v_i][nNbrs[v_i]] = domain;
					nNbrs[v_i]++;
				}
				inds[1]--;
			}

			if( nDims == 3 ) {
				if( inds[2] > 0 ) {
					inds[2]--;
					domain = Grid_Project( grids[0][0], inds );
					if( Sync_TryGlobalToDomain( sync, domain, &domain ) ) {
						nbrs[v_i][nNbrs[v_i]] = domain;
						nNbrs[v_i]++;
					}
					inds[2]++;
				}

				if( inds[2] < grids[0][0]->sizes[2] - 1 ) {
					inds[2]++;
					domain = Grid_Project( grids[0][0], inds );
					if( Sync_TryGlobalToDomain( sync, domain, &domain ) ) {
						nbrs[v_i][nNbrs[v_i]] = domain;
						nNbrs[v_i]++;
					}
					inds[2]--;
				}
			}
		}

		MeshTopology_SetIncidence( topo, MT_VERTEX, v_i, MT_VERTEX, nNbrs[v_i], nbrs[v_i] );
	}

	FreeArray( nNbrs );
	FreeArray( nbrs );
	FreeArray( inds );

	MPI_Barrier( self->mpiComm );
	Journal_Printf( stream, "... done.\n" );
	Stream_UnIndent( stream );
}

void CartesianGenerator_MapToDomain( CartesianGenerator* self, Sync* sync, 
				     unsigned nIncEls, unsigned* incEls )
{
	unsigned	inc_i;

	assert( self );
	assert( sync );
	assert( nIncEls );
	assert( incEls );

	for( inc_i = 0; inc_i < nIncEls; inc_i++ )
		incEls[inc_i] = Sync_GlobalToDomain( sync, incEls[inc_i] );
}

void CartesianGenerator_GenGeom( CartesianGenerator* self, Mesh* mesh ) {
	Stream*		stream = Journal_Register( Info_Type, self->type );
	Sync*		sync;
	Grid*		grid;
	unsigned*	inds;
	double*		steps;
	unsigned	n_i, d_i;

	assert( self );
	assert( mesh );

	Journal_Printf( stream, "Generating geometry...\n" );
	Stream_Indent( stream );

	/* Build grid and space for indices. */
	grid = self->vertGrid;
	inds = Memory_Alloc_Array_Unnamed( unsigned, mesh->topo->nDims );

	/* Calculate steps. */
	steps = Memory_Alloc_Array_Unnamed( double, mesh->topo->nDims );
	for( d_i = 0; d_i < mesh->topo->nDims; d_i++ )
		steps[d_i] = self->crdMax[d_i] - self->crdMin[d_i];

	/* Allocate for coordinates. */
	sync = (Sync*)MeshTopology_GetDomain( mesh->topo, 0 );
	mesh->verts = AllocNamedArray2D( double, Sync_GetNumDomains( sync ), 
					 mesh->topo->nDims, 
					 "Mesh::verts" );

	/* Loop over domain nodes. */
	for( n_i = 0; n_i < Sync_GetNumDomains( sync ); n_i++ ) {
		double*		vert;
		unsigned	gNode;

		gNode = Sync_DomainToGlobal( sync, n_i );
		Grid_Lift( grid, gNode, inds );
		vert = Mesh_GetVertex( mesh, n_i );

		/* Calculate coordinate. */
		for( d_i = 0; d_i < mesh->topo->nDims; d_i++ ) {
			vert[d_i] = self->crdMin[d_i] + 
				((double)inds[d_i] / (double)(grid->sizes[d_i] - 1)) * steps[d_i];
		}
	}

	/* Free resources. */
	FreeArray( inds );
	FreeArray( steps );

	MPI_Barrier( self->mpiComm );
	Journal_Printf( stream, "... done.\n" );
	Stream_UnIndent( stream );
}

void CartesianGenerator_Destruct( CartesianGenerator* self ) {
	assert( self );

	CartesianGenerator_DestructTopology( self );
	CartesianGenerator_DestructGeometry( self );
	KillArray( self->minDecomp );
	KillArray( self->maxDecomp );
	NewClass_RemoveRef( self->comm );
}

void CartesianGenerator_DestructTopology( CartesianGenerator* self ) {
	assert( self );

	self->maxDecompDims = 0;
	memset( self->minDecomp, 0, self->nDims * sizeof(unsigned) );
	memset( self->maxDecomp, 0, self->nDims * sizeof(unsigned) );
	KillObject( self->vertGrid );
	KillObject( self->elGrid );
	KillObject( self->procGrid );
	KillArray( self->origin );
	KillArray( self->range );
	KillArray( self->vertOrigin );
	KillArray( self->vertRange );
}

void CartesianGenerator_DestructGeometry( CartesianGenerator* self ) {
	assert( self );

	FreeArray( self->crdMin );
	FreeArray( self->crdMax );
}
