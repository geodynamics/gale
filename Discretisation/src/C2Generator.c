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
** $Id: C2Generator.c 3584 2006-05-16 11:11:07Z PatrickSunter $
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
const Type C2Generator_Type = "C2Generator";


/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/

C2Generator* C2Generator_New( Name name ) {
	return _C2Generator_New( sizeof(C2Generator), 
				 C2Generator_Type, 
				 _C2Generator_Delete, 
				 _C2Generator_Print, 
				 NULL, 
				 (void* (*)(Name))_C2Generator_New, 
				 _C2Generator_Construct, 
				 _C2Generator_Build, 
				 _C2Generator_Initialise, 
				 _C2Generator_Execute, 
				 _C2Generator_Destroy, 
				 name, 
				 NON_GLOBAL, 
				 CartesianGenerator_SetDimSize, 
				 CartesianGenerator_Generate, 
				 C2Generator_SetTopologyParams, 
				 _CartesianGenerator_GenElements, 
				 _CartesianGenerator_GenFaces, 
				 _CartesianGenerator_GenEdges, 
				 _CartesianGenerator_GenVertices, 
				 C2Generator_GenElementVertexInc, 
				 _CartesianGenerator_GenVolumeEdgeInc, 
				 _CartesianGenerator_GenVolumeFaceInc, 
				 C2Generator_GenFaceVertexInc, 
				 _CartesianGenerator_GenFaceEdgeInc, 
				 C2Generator_GenEdgeVertexInc, 
				 C2Generator_GenElementTypes );
}

C2Generator* _C2Generator_New( C2GENERATOR_DEFARGS ) {
	C2Generator*	self;

	/* Allocate memory */
	assert( sizeOfSelf >= sizeof(C2Generator) );
	self = (C2Generator*)_CartesianGenerator_New( CARTESIANGENERATOR_PASSARGS );

	/* Virtual info */

	/* C2Generator info */
	_C2Generator_Init( self );

	return self;
}

void _C2Generator_Init( C2Generator* self ) {
	assert( self );
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _C2Generator_Delete( void* meshGenerator ) {
	C2Generator*	self = (C2Generator*)meshGenerator;

	/* Delete the parent. */
	_CartesianGenerator_Delete( self );
}

void _C2Generator_Print( void* meshGenerator, Stream* stream ) {
	C2Generator*	self = (C2Generator*)meshGenerator;
	
	/* Set the Journal for printing informations */
	Stream* meshGeneratorStream;
	meshGeneratorStream = Journal_Register( InfoStream_Type, "C2GeneratorStream" );

	assert( self );

	/* Print parent */
	Journal_Printf( stream, "C2Generator (ptr): (%p)\n", self );
	_CartesianGenerator_Print( self, stream );
}

void _C2Generator_Construct( void* meshGenerator, Stg_ComponentFactory* cf, void* data ) {
	_CartesianGenerator_Construct( meshGenerator, cf, data );
}

void _C2Generator_Build( void* meshGenerator, void* data ) {
	_CartesianGenerator_Build( meshGenerator, data );
}

void _C2Generator_Initialise( void* meshGenerator, void* data ) {
}

void _C2Generator_Execute( void* meshGenerator, void* data ) {
}

void _C2Generator_Destroy( void* meshGenerator, void* data ) {
}

void C2Generator_SetTopologyParams( void* meshGenerator, unsigned* sizes, 
				    unsigned maxDecompDims, unsigned* minDecomp, unsigned* maxDecomp )
{
	C2Generator*	self = (C2Generator*)meshGenerator;
	unsigned*	vertSizes;
	unsigned	d_i;

	assert( self );

	_CartesianGenerator_SetTopologyParams( self, sizes, 
					       maxDecompDims, minDecomp, maxDecomp );
	vertSizes = AllocArray( unsigned, self->vertGrid->nDims );
	for( d_i = 0; d_i < self->vertGrid->nDims; d_i++ ) {
		vertSizes[d_i] = self->vertGrid->sizes[d_i] * 2 - 1;
		self->vertOrigin[d_i] *= 2;
		self->vertRange[d_i] = self->vertRange[d_i] * 2 - 1;
	}
	Grid_SetSizes( self->vertGrid, vertSizes );
	FreeArray( vertSizes );
}

void C2Generator_GenElementVertexInc( void* meshGenerator, MeshTopology* topo, Grid*** grids ) {
	C2Generator*	self = (C2Generator*)meshGenerator;
	Stream*		stream = Journal_Register( Info_Type, self->type );
	unsigned*	incEls;
	unsigned*	dimInds;
	unsigned	vertsPerEl;
	unsigned	nDims;
	unsigned	e_i, d_i;
	int nDomainEls;

	assert( self );
	assert( topo );
	assert( grids );

	Journal_Printf( stream, "Generating element-vertex incidence...\n" );
	Stream_Indent( stream );

	vertsPerEl = (topo->nDims == 1) ? 3 : (topo->nDims == 2) ? 9 : 27;

	nDims = topo->nDims;
	nDomainEls = Sync_GetNumDomains( MeshTopology_GetDomain( topo, nDims ) );
	incEls = Memory_Alloc_Array_Unnamed( unsigned, vertsPerEl );
	dimInds = Memory_Alloc_Array_Unnamed( unsigned, topo->nDims );
	for( e_i = 0; e_i < nDomainEls; e_i++ ) {
		unsigned	gInd = Sync_DomainToGlobal( MeshTopology_GetDomain( topo, nDims ), e_i );
		unsigned	curNode = 0;

		Grid_Lift( grids[topo->nDims][0], gInd, dimInds );
		for( d_i = 0; d_i < nDims; d_i++ )
			dimInds[d_i] *= 2;

		incEls[curNode++] = Grid_Project( grids[0][0], dimInds );
		dimInds[0]++;
		incEls[curNode++] = Grid_Project( grids[0][0], dimInds );
		dimInds[0]++;
		incEls[curNode++] = Grid_Project( grids[0][0], dimInds );
		dimInds[0] -= 2;

		if( topo->nDims >= 2 ) {
			dimInds[1]++;
			incEls[curNode++] = Grid_Project( grids[0][0], dimInds );
			dimInds[0]++;
			incEls[curNode++] = Grid_Project( grids[0][0], dimInds );
			dimInds[0]++;
			incEls[curNode++] = Grid_Project( grids[0][0], dimInds );
			dimInds[0] -= 2;
			dimInds[1]++;
			incEls[curNode++] = Grid_Project( grids[0][0], dimInds );
			dimInds[0]++;
			incEls[curNode++] = Grid_Project( grids[0][0], dimInds );
			dimInds[0]++;
			incEls[curNode++] = Grid_Project( grids[0][0], dimInds );
			dimInds[0] -= 2;
			dimInds[1] -= 2;

			if( topo->nDims >= 3 ) {
				dimInds[2]++;
				incEls[curNode++] = Grid_Project( grids[0][0], dimInds );
				dimInds[0]++;
				incEls[curNode++] = Grid_Project( grids[0][0], dimInds );
				dimInds[0]++;
				incEls[curNode++] = Grid_Project( grids[0][0], dimInds );
				dimInds[0] -= 2;
				dimInds[1]++;
				incEls[curNode++] = Grid_Project( grids[0][0], dimInds );
				dimInds[0]++;
				incEls[curNode++] = Grid_Project( grids[0][0], dimInds );
				dimInds[0]++;
				incEls[curNode++] = Grid_Project( grids[0][0], dimInds );
				dimInds[0] -= 2;
				dimInds[1]++;
				incEls[curNode++] = Grid_Project( grids[0][0], dimInds );
				dimInds[0]++;
				incEls[curNode++] = Grid_Project( grids[0][0], dimInds );
				dimInds[0]++;
				incEls[curNode++] = Grid_Project( grids[0][0], dimInds );
				dimInds[0] -= 2;
				dimInds[1] -= 2;
				dimInds[2]++;
				incEls[curNode++] = Grid_Project( grids[0][0], dimInds );
				dimInds[0]++;
				incEls[curNode++] = Grid_Project( grids[0][0], dimInds );
				dimInds[0]++;
				incEls[curNode++] = Grid_Project( grids[0][0], dimInds );
				dimInds[0] -= 2;
				dimInds[1]++;
				incEls[curNode++] = Grid_Project( grids[0][0], dimInds );
				dimInds[0]++;
				incEls[curNode++] = Grid_Project( grids[0][0], dimInds );
				dimInds[0]++;
				incEls[curNode++] = Grid_Project( grids[0][0], dimInds );
				dimInds[0] -= 2;
				dimInds[1]++;
				incEls[curNode++] = Grid_Project( grids[0][0], dimInds );
				dimInds[0]++;
				incEls[curNode++] = Grid_Project( grids[0][0], dimInds );
				dimInds[0]++;
				incEls[curNode++] = Grid_Project( grids[0][0], dimInds );
				dimInds[0] -= 2;
				dimInds[1] -= 2;
				dimInds[2] -= 2;
			}
		}
		CartesianGenerator_MapToDomain( (CartesianGenerator*)self, MeshTopology_GetDomain( topo, 0), 
						vertsPerEl, incEls );
		MeshTopology_SetIncidence( topo, topo->nDims, e_i, MT_VERTEX, vertsPerEl, incEls );
	}

	FreeArray( incEls );
	FreeArray( dimInds );

	MPI_Barrier( self->mpiComm );
	Journal_Printf( stream, "... done.\n" );
	Stream_UnIndent( stream );
}

void C2Generator_GenFaceVertexInc( void* meshGenerator, MeshTopology* topo, Grid*** grids ) {
	abort();
}

void C2Generator_GenEdgeVertexInc( void* meshGenerator, MeshTopology* topo, Grid*** grids ) {
	abort();
}

void C2Generator_GenElementTypes( void* meshGenerator, Mesh* mesh ) {
	C2Generator*	self = (C2Generator*)meshGenerator;
	Stream*		stream;
	unsigned	nDomainEls;
	unsigned	vertMap[8] = {0, 2, 6, 8, 18, 20, 24, 26};
	unsigned	e_i;

	assert( self );

	stream = Journal_Register( Info_Type, self->type );
	Journal_Printf( stream, "Generating element types...\n" );
	Stream_Indent( stream );

	mesh->nElTypes = 1;
	mesh->elTypes = AllocArray( Mesh_ElementType*, mesh->nElTypes );
	mesh->elTypes[0] = (Mesh_ElementType*)Mesh_HexType_New();
	Mesh_ElementType_SetMesh( mesh->elTypes[0], mesh );
	Mesh_HexType_SetVertexMap( mesh->elTypes[0], vertMap );
	nDomainEls = Mesh_GetDomainSize( mesh, Mesh_GetDimSize( mesh ) );
	mesh->elTypeMap = AllocArray( unsigned, nDomainEls );
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


/*----------------------------------------------------------------------------------------------------------------------------------
** Private Functions
*/
