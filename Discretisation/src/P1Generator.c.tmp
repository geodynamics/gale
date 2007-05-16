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
** $Id: P1Generator.c 3584 2006-05-16 11:11:07Z PatrickSunter $
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
const Type P1Generator_Type = "P1Generator";


/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/

P1Generator* P1Generator_New( Name name ) {
	return _P1Generator_New( sizeof(P1Generator), 
				 P1Generator_Type, 
				 _P1Generator_Delete, 
				 _P1Generator_Print, 
				 NULL, 
				 (void* (*)(Name))_P1Generator_New, 
				 _P1Generator_Construct, 
				 _P1Generator_Build, 
				 _P1Generator_Initialise, 
				 _P1Generator_Execute, 
				 _P1Generator_Destroy, 
				 name, 
				 NON_GLOBAL, 
				 _MeshGenerator_SetDimSize, 
				 P1Generator_Generate );
}

P1Generator* _P1Generator_New( P1GENERATOR_DEFARGS ) {
	P1Generator*	self;

	/* Allocate memory */
	assert( sizeOfSelf >= sizeof(P1Generator) );
	self = (P1Generator*)_MeshGenerator_New( MESHGENERATOR_PASSARGS );

	/* Virtual info */

	/* P1Generator info */
	_P1Generator_Init( self );

	return self;
}

void _P1Generator_Init( P1Generator* self ) {
	assert( self && Stg_CheckType( self, P1Generator ) );

	self->elMesh = NULL;
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _P1Generator_Delete( void* generator ) {
	P1Generator*	self = (P1Generator*)generator;

	/* Delete the parent. */
	_MeshGenerator_Delete( self );
}

void _P1Generator_Print( void* generator, Stream* stream ) {
	P1Generator*	self = (P1Generator*)generator;
	
	/* Set the Journal for printing informations */
	Stream* generatorStream;
	generatorStream = Journal_Register( InfoStream_Type, "P1GeneratorStream" );

	/* Print parent */
	Journal_Printf( stream, "P1Generator (ptr): (%p)\n", self );
	_MeshGenerator_Print( self, stream );
}

void _P1Generator_Construct( void* generator, Stg_ComponentFactory* cf, void* data ) {
	P1Generator*	self = (P1Generator*)generator;
	FeMesh*		elMesh;

	assert( self );
	assert( cf );

	_MeshGenerator_Construct( self, cf, data );

	elMesh = Stg_ComponentFactory_ConstructByKey( cf, self->name, "elementMesh", FeMesh, True, data );
	P1Generator_SetElementMesh( self, elMesh );
}

void _P1Generator_Build( void* generator, void* data ) {
	_MeshGenerator_Build( generator, data );
}

void _P1Generator_Initialise( void* generator, void* data ) {
	_MeshGenerator_Initialise( generator, data );
}

void _P1Generator_Execute( void* generator, void* data ) {
}

void _P1Generator_Destroy( void* generator, void* data ) {
}

void P1Generator_Generate( void* generator, void* _mesh ) {
	P1Generator*	self = (P1Generator*)generator;
	FeMesh*		mesh = (FeMesh*)_mesh;

	assert( self && Stg_CheckType( self, P1Generator ) );
	assert( mesh && Stg_CheckType( mesh, FeMesh ) );

	P1Generator_BuildTopology( self, mesh );
	P1Generator_BuildGeometry( self, mesh );
	P1Generator_BuildElementTypes( self, mesh );
}


/*--------------------------------------------------------------------------------------------------------------------------
** Public Functions
*/

void P1Generator_SetElementMesh( void* generator, void* mesh ) {
	P1Generator*	self = (P1Generator*)generator;

	assert( self && Stg_CheckType( self, P1Generator ) );

	self->elMesh = mesh;
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Private Functions
*/

void P1Generator_BuildTopology( P1Generator* self, FeMesh* mesh ) {
	FeMesh*		elMesh;
	MeshTopology	*topo, *elTopo;
	Decomp_Sync*	elSync;
	unsigned	nDims;
	unsigned	*nIncEls, **incEls;
	unsigned	nLocalEls;
	Decomp*		nodeDecomp;
	Decomp_Sync*	nodeSync;
	unsigned	nNodes, *nodes;
	unsigned	globalEl;
	unsigned	e_i;

	assert( self );
	assert( mesh );

	elMesh = self->elMesh;
	nDims = Mesh_GetDimSize( elMesh );
	elTopo = Mesh_GetTopology( elMesh );
	elSync = MeshTopology_GetSync( elTopo, nDims );

	topo = Mesh_GetTopology( mesh );
	MeshTopology_SetDimSize( topo, nDims );
	MeshTopology_SetSync( topo, nDims, elSync );
	topo->shadowDepth = elTopo->shadowDepth;

	/* Generate the vertices. */
	nLocalEls = FeMesh_GetElementLocalSize( elMesh );
	nodeDecomp = Decomp_New();
	Decomp_SetComm( nodeDecomp, Decomp_GetComm( Decomp_Sync_GetDecomp( elSync ) ) );
	nNodes = nLocalEls * 3;
	nodes = AllocArray( unsigned, nNodes );
	for( e_i = 0; e_i < nLocalEls; e_i++ ) {
		globalEl = FeMesh_ElementDomainToGlobal( mesh, e_i );
		nodes[3 * e_i + 0] = 3 * globalEl + 0;
		nodes[3 * e_i + 1] = 3 * globalEl + 1;
		nodes[3 * e_i + 2] = 3 * globalEl + 2;
	}
	Decomp_SetLocals( nodeDecomp, nNodes, nodes );
	FreeArray( nodes );
	nodeSync = Decomp_Sync_New();
	Decomp_Sync_SetDecomp( nodeSync, nodeDecomp );
	MeshTopology_SetSync( topo, MT_VERTEX, nodeSync );

	/* Generate the incidence. */
	nIncEls = AllocArray( unsigned, nLocalEls );
	incEls = AllocArray2D( unsigned, nLocalEls, 3 );
	for( e_i = 0; e_i < nLocalEls; e_i++ ) {
		nIncEls[e_i] = 3;
		incEls[e_i][0] = e_i * 3 + 0;
		incEls[e_i][1] = e_i * 3 + 1;
		incEls[e_i][2] = e_i * 3 + 2;
	}
	MeshTopology_SetIncidence( topo, nDims, MT_VERTEX, nIncEls, incEls );
	FreeArray( nIncEls );
	FreeArray( incEls );

	MeshTopology_Invert( topo, MT_VERTEX, nDims );
}

void P1Generator_BuildGeometry( P1Generator* self, FeMesh* mesh ) {
	FeMesh*			elMesh;
	double			localCrds[4][3] = {{-0.5, -0.5, -0.5}, {0.5, -0.5, -0.5}, 
						   {0.0, 0.707106, 0.0}, 
						   {0.0, 0.0, 0.707106}};
	double			*globalCrd, *vert;
	unsigned		nDims;
	unsigned		nLocalEls;
	unsigned		e_i, c_i;

	assert( self );
	assert( mesh );

	elMesh = self->elMesh;
	nDims = Mesh_GetDimSize( elMesh );
	nLocalEls = FeMesh_GetElementDomainSize( elMesh );
	mesh->verts = AllocArray2D( double, nLocalEls * 3, nDims );
	globalCrd = AllocArray( double, nDims );
	for( e_i = 0; e_i < nLocalEls; e_i++ ) {
		for( c_i = 0; c_i < 3; c_i++ ) {
			FeMesh_CoordLocalToGlobal( elMesh, e_i, localCrds[c_i], globalCrd );
			vert = Mesh_GetVertex( mesh, 3 * e_i + c_i );
			memcpy( vert, globalCrd, nDims * sizeof(double) );
		}
	}
	FreeArray( globalCrd );
}

void P1Generator_BuildElementTypes( P1Generator* self, FeMesh* mesh ) {
	unsigned		nLocalEls;
	Mesh_Algorithms*	algs;
	unsigned		e_i;

	assert( self );
	assert( mesh );

	mesh->nElTypes = 1;
	mesh->elTypes = AllocArray( Mesh_ElementType*, mesh->nElTypes );
	mesh->elTypes[0] = (Mesh_ElementType*)Mesh_CentroidType_New();
	Mesh_ElementType_SetMesh( mesh->elTypes[0], mesh );
	Mesh_CentroidType_SetElementMesh( mesh->elTypes[0], self->elMesh );
	nLocalEls = FeMesh_GetElementLocalSize( mesh );
	mesh->elTypeMap = AllocArray( unsigned, nLocalEls );
	for( e_i = 0; e_i < nLocalEls; e_i++ )
		mesh->elTypeMap[e_i] = 0;

	algs = (Mesh_Algorithms*)Mesh_CentroidAlgorithms_New( "" );
	Mesh_CentroidAlgorithms_SetElementMesh( algs, self->elMesh );
	Mesh_SetAlgorithms( mesh, algs );
}
