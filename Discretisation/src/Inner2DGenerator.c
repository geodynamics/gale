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
** $Id: Inner2DGenerator.c 3584 2006-05-16 11:11:07Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <mpi.h>

#include "StGermain/StGermain.h"
#include "Discretisation.h"


/* Textual name of this class */
const Type Inner2DGenerator_Type = "Inner2DGenerator";


/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/

Inner2DGenerator* Inner2DGenerator_New( Name name ) {
	return _Inner2DGenerator_New( sizeof(Inner2DGenerator), 
				 Inner2DGenerator_Type, 
				 _Inner2DGenerator_Delete, 
				 _Inner2DGenerator_Print, 
				 NULL, 
				 (void* (*)(Name))_Inner2DGenerator_New, 
				 _Inner2DGenerator_Construct, 
				 _Inner2DGenerator_Build, 
				 _Inner2DGenerator_Initialise, 
				 _Inner2DGenerator_Execute, 
				 _Inner2DGenerator_Destroy, 
				 name, 
				 NON_GLOBAL, 
				 _MeshGenerator_SetDimSize, 
				 Inner2DGenerator_Generate );
}

Inner2DGenerator* _Inner2DGenerator_New( Inner2DGENERATOR_DEFARGS ) {
	Inner2DGenerator*	self;

	/* Allocate memory */
	assert( sizeOfSelf >= sizeof(Inner2DGenerator) );
	self = (Inner2DGenerator*)_MeshGenerator_New( MESHGENERATOR_PASSARGS );

	/* Virtual info */

	/* Inner2DGenerator info */
	_Inner2DGenerator_Init( self );

	return self;
}

void _Inner2DGenerator_Init( Inner2DGenerator* self ) {
	assert( self && Stg_CheckType( self, Inner2DGenerator ) );

	self->elMesh = NULL;
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _Inner2DGenerator_Delete( void* generator ) {
	Inner2DGenerator*	self = (Inner2DGenerator*)generator;

	/* Delete the parent. */
	_MeshGenerator_Delete( self );
}

void _Inner2DGenerator_Print( void* generator, Stream* stream ) {
	Inner2DGenerator*	self = (Inner2DGenerator*)generator;
	
	/* Set the Journal for printing informations */
	Stream* generatorStream;
	generatorStream = Journal_Register( InfoStream_Type, "Inner2DGeneratorStream" );

	/* Print parent */
	Journal_Printf( stream, "Inner2DGenerator (ptr): (%p)\n", self );
	_MeshGenerator_Print( self, stream );
}

void _Inner2DGenerator_Construct( void* generator, Stg_ComponentFactory* cf, void* data ) {
	Inner2DGenerator*	self = (Inner2DGenerator*)generator;
	Mesh*		elMesh;

	assert( self );
	assert( cf );

	_MeshGenerator_Construct( self, cf, data );

	elMesh = Stg_ComponentFactory_ConstructByKey( cf, self->name, "elementMesh", Mesh, True, data );
	Inner2DGenerator_SetElementMesh( self, elMesh );
}

void _Inner2DGenerator_Build( void* generator, void* data ) {
	_MeshGenerator_Build( generator, data );
}

void _Inner2DGenerator_Initialise( void* generator, void* data ) {
	_MeshGenerator_Initialise( generator, data );
}

void _Inner2DGenerator_Execute( void* generator, void* data ) {
}

void _Inner2DGenerator_Destroy( void* generator, void* data ) {
}

void Inner2DGenerator_Generate( void* generator, void* _mesh ) {
	Inner2DGenerator*	self = (Inner2DGenerator*)generator;
	FeMesh*		mesh = (FeMesh*)_mesh;
	Grid**		grid;
	Grid*		elGrid;

	assert( self && Stg_CheckType( self, Inner2DGenerator ) );
	assert( mesh && Stg_CheckType( mesh, FeMesh ) );

	Inner2DGenerator_BuildTopology( self, mesh );
	Inner2DGenerator_BuildGeometry( self, mesh );
	Inner2DGenerator_BuildElementTypes( self, mesh );

	elGrid = *(Grid**)ExtensionManager_Get( self->elMesh->info, self->elMesh, 
					       ExtensionManager_GetHandle( self->elMesh->info, "elementGrid" ) );
	ExtensionManager_Add( mesh->info, "elementGrid", sizeof(Grid*) );
	grid = (Grid**)ExtensionManager_Get( mesh->info, mesh, 
					     ExtensionManager_GetHandle( mesh->info, "elementGrid" ) );
	*grid = Grid_New();
	Grid_SetNumDims( *grid, elGrid->nDims );
	Grid_SetSizes( *grid, elGrid->sizes );
}


/*--------------------------------------------------------------------------------------------------------------------------
** Public Functions
*/

void Inner2DGenerator_SetElementMesh( void* generator, void* mesh ) {
	Inner2DGenerator*	self = (Inner2DGenerator*)generator;

	assert( self && Stg_CheckType( self, Inner2DGenerator ) );

	self->elMesh = mesh;
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Private Functions
*/

void Inner2DGenerator_BuildTopology( Inner2DGenerator* self, FeMesh* mesh ) {
	Mesh*		elMesh;
	MeshTopology	*topo, *elTopo;
	unsigned	nDims;
	unsigned	nIncEls, *incEls;
	unsigned	nDomainEls;
	Decomp		*elDecomp, *nodeDecomp;
	Sync		*elSync, *nodeSync;
	int		nLocals, *locals;
	int		nRemotes, *remotes;
	unsigned	global;
	unsigned	e_i, l_i, r_i;

	assert( self );
	assert( mesh );

	elMesh = self->elMesh;
	nDims = Mesh_GetDimSize( elMesh );
	elTopo = Mesh_GetTopology( elMesh );
	elSync = Mesh_GetSync( elMesh, nDims );

	topo = Mesh_GetTopology( mesh );
	MeshTopology_SetComm( topo, MeshTopology_GetComm( elTopo ) );
	MeshTopology_SetNumDims( topo, nDims );
	MeshTopology_SetDomain( topo, nDims, elSync );

	/* Need to redefine the nodes, nDims + 1 per parent element. */
	elDecomp = (Decomp*)Sync_GetDecomp( elSync );
	nodeDecomp = Decomp_New();
	nLocals = Decomp_GetNumLocals( elDecomp ) * 3;
	locals = MemArray( int, nLocals, Inner2DGenerator_Type );
	for( l_i = 0; l_i < Decomp_GetNumLocals( elDecomp ); l_i++ ) {
		global = Decomp_LocalToGlobal( elDecomp, l_i );
		locals[l_i * 3 + 0] = global * 3;
		locals[l_i * 3 + 1] = global * 3 + 1;
		locals[l_i * 3 + 2] = global * 3 + 2;
	}
	Decomp_SetLocals( nodeDecomp, nLocals, locals );
	MemFree( locals );

	nodeSync = Sync_New();
	Sync_SetComm( nodeSync, Sync_GetComm( elSync ) );
	Sync_SetDecomp( nodeSync, nodeDecomp );
	nRemotes = Sync_GetNumRemotes( elSync ) * 3;
	remotes = MemArray( int, nRemotes, Inner2DGenerator_Type );
	for( r_i = 0; r_i < Sync_GetNumRemotes( elSync ); r_i++ ) {
		global = Sync_RemoteToGlobal( elSync, r_i );
		remotes[r_i * 3 + 0] = global * 3;
		remotes[r_i * 3 + 1] = global * 3 + 1;
		remotes[r_i * 3 + 2] = global * 3 + 2;
	}
	Sync_SetRemotes( nodeSync, nRemotes, remotes );
	MemFree( remotes );

	MeshTopology_SetDomain( topo, 0, nodeSync );

	/* Same shadow depth. */
	topo->shadDepth = elTopo->shadDepth;

	/* Build the incidence. */
	nDomainEls = Mesh_GetDomainSize( elMesh, nDims );
	nIncEls = 3;
	incEls = MemArray( unsigned, 3, Inner2DGenerator_Type );
	for( e_i = 0; e_i < nDomainEls; e_i++ ) {
		incEls[0] = e_i * 3;
		incEls[1] = e_i * 3 + 1;
		incEls[2] = e_i * 3 + 2;
		MeshTopology_SetIncidence( topo, nDims, e_i, 0, nIncEls, (int*)incEls );
	}
	FreeArray( incEls );

	MeshTopology_InvertIncidence( topo, MT_VERTEX, nDims );
}

void Inner2DGenerator_BuildGeometry( Inner2DGenerator* self, FeMesh* mesh ) {
	Mesh*		elMesh;
	double		localCrds[3][2] = {{-0.5, -0.5}, 
					   {0.5, -0.5}, 
					   {0, 0.5}};
	double		globalCrd[2];
	double		*vert;
	unsigned	nDims;
	unsigned	nDomainEls;
	unsigned	e_i;

	assert( self );
	assert( mesh );

	elMesh = self->elMesh;
	nDims = Mesh_GetDimSize( elMesh );
	nDomainEls = Mesh_GetDomainSize( elMesh, nDims );
	mesh->verts = AllocArray2D( double, nDomainEls, nDims );
	for( e_i = 0; e_i < nDomainEls; e_i++ ) {
		FeMesh_CoordLocalToGlobal( elMesh, e_i, localCrds[0], globalCrd );
		vert = Mesh_GetVertex( mesh, e_i );
		memcpy( vert, globalCrd, nDims * sizeof(double) );
	}
}

void Inner2DGenerator_BuildElementTypes( Inner2DGenerator* self, FeMesh* mesh ) {
	unsigned		nDomainEls;
	Mesh_Algorithms*	algs;
	unsigned		e_i;

	assert( self );
	assert( mesh );

	mesh->nElTypes = 1;
	mesh->elTypes = AllocNamedArray( Mesh_ElementType*, mesh->nElTypes, "Mesh::elTypes" );
	mesh->elTypes[0] = (Mesh_ElementType*)Mesh_CentroidType_New();
	Mesh_ElementType_SetMesh( mesh->elTypes[0], mesh );
	Mesh_CentroidType_SetElementMesh( mesh->elTypes[0], self->elMesh );
	nDomainEls = Mesh_GetDomainSize( mesh, Mesh_GetDimSize( mesh ) );
	mesh->elTypeMap = AllocNamedArray( unsigned, nDomainEls, "Mesh::elTypeMap" );
	for( e_i = 0; e_i < nDomainEls; e_i++ )
		mesh->elTypeMap[e_i] = 0;

	algs = Mesh_CentroidAlgorithms_New( "" );
	Mesh_CentroidAlgorithms_SetElementMesh( algs, self->elMesh );
	Mesh_SetAlgorithms( mesh, algs );
}
