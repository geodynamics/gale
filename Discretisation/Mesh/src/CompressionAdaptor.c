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
** $Id: CompressionAdaptor.c 3584 2006-05-16 11:11:07Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
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
#include "Decomp_Sync.h"
#include "MeshTopology.h"
#include "MeshClass.h"
#include "MeshGenerator.h"
#include "MeshAdaptor.h"
#include "CompressionAdaptor.h"


typedef double (CompressionAdaptor_DeformFunc)( CompressionAdaptor* self, Mesh* mesh, 
					    unsigned* globalSize, unsigned vertex, unsigned* vertexInds );


/* Textual name of this class */
const Type CompressionAdaptor_Type = "CompressionAdaptor";


/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/

CompressionAdaptor* CompressionAdaptor_New( Name name ) {
	return _CompressionAdaptor_New( sizeof(CompressionAdaptor), 
				    CompressionAdaptor_Type, 
				    _CompressionAdaptor_Delete, 
				    _CompressionAdaptor_Print, 
				    NULL, 
				    (void* (*)(Name))_CompressionAdaptor_New, 
				    _CompressionAdaptor_Construct, 
				    _CompressionAdaptor_Build, 
				    _CompressionAdaptor_Initialise, 
				    _CompressionAdaptor_Execute, 
				    _CompressionAdaptor_Destroy, 
				    name, 
				    NON_GLOBAL, 
				    _MeshGenerator_SetDimSize, 
				    CompressionAdaptor_Generate );
}

CompressionAdaptor* _CompressionAdaptor_New( COMPRESSIONADAPTOR_DEFARGS ) {
	CompressionAdaptor* self;
	
	/* Allocate memory */
	assert( sizeOfSelf >= sizeof(CompressionAdaptor) );
	self = (CompressionAdaptor*)_MeshGenerator_New( MESHADAPTOR_PASSARGS );

	/* Virtual info */

	/* CompressionAdaptor info */
	_CompressionAdaptor_Init( self );

	return self;
}

void _CompressionAdaptor_Init( CompressionAdaptor* self ) {
	self->grad = 0.0;
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _CompressionAdaptor_Delete( void* adaptor ) {
	CompressionAdaptor*	self = (CompressionAdaptor*)adaptor;

	/* Delete the parent. */
	_MeshGenerator_Delete( self );
}

void _CompressionAdaptor_Print( void* adaptor, Stream* stream ) {
	CompressionAdaptor*	self = (CompressionAdaptor*)adaptor;
	
	/* Set the Journal for printing informations */
	Stream* adaptorStream;
	adaptorStream = Journal_Register( InfoStream_Type, "CompressionAdaptorStream" );

	/* Print parent */
	Journal_Printf( stream, "CompressionAdaptor (ptr): (%p)\n", self );
	_MeshGenerator_Print( self, stream );
}

void _CompressionAdaptor_Construct( void* adaptor, Stg_ComponentFactory* cf, void* data ) {
	CompressionAdaptor*	self = (CompressionAdaptor*)adaptor;

	assert( self );
	assert( cf );

	/* Call parent construct. */
	_MeshAdaptor_Construct( self, cf, data );

	self->grad = Stg_ComponentFactory_GetDouble( cf, self->name, "initialSeparation", 0.0 );
}

void _CompressionAdaptor_Build( void* adaptor, void* data ) {
	_MeshAdaptor_Build( adaptor, data );
}

void _CompressionAdaptor_Initialise( void* adaptor, void* data ) {
	_MeshAdaptor_Initialise( adaptor, data );
}

void _CompressionAdaptor_Execute( void* adaptor, void* data ) {
}

void _CompressionAdaptor_Destroy( void* adaptor, void* data ) {
}

/*--------------------------------------------------------------------------------------------------------------------------
** Public Functions
*/

void CompressionAdaptor_Generate( void* adaptor, void* _mesh ) {
	CompressionAdaptor*		self = (CompressionAdaptor*)adaptor;
	Mesh*				mesh = (Mesh*)_mesh;
	Grid				*grid;
	unsigned*			inds;
	double				min;
	unsigned			gNode;
	double				x;
	unsigned			d_i, n_i;

	/* Build base mesh, which is assumed to be cartesian. */
	MeshGenerator_Generate( self->generator, mesh );

	/* If we're not 2D or 3D, forget about it. */
	if( mesh->topo->nDims != 2 && mesh->topo->nDims != 3 )
		return;

	/* Extract the cartesian information. */
	grid = *(Grid**)ExtensionManager_Get( mesh->info, mesh, 
					      ExtensionManager_GetHandle( mesh->info, "vertexGrid" ) );
	inds = AllocArray( unsigned, Mesh_GetDimSize( mesh ) );
	for( d_i = 0; d_i < Mesh_GetDimSize( mesh ); d_i++ )
		inds[d_i] = 0;
	gNode = Grid_Project( grid, inds );
	insist( Mesh_GlobalToDomain( mesh, MT_VERTEX, gNode, &gNode ) );
	min = mesh->verts[gNode][1];

	/* Loop over domain nodes. */
	for( n_i = 0; n_i < MeshTopology_GetDomainSize( mesh->topo, MT_VERTEX ); n_i++ ) {
		gNode = MeshTopology_DomainToGlobal( mesh->topo, MT_VERTEX, n_i );
		Grid_Lift( grid, gNode, inds );

		/* Deform this node. */
		x = mesh->verts[n_i][1] - min;
		mesh->verts[n_i][1] = 0.5 * self->grad * x * x + x + min;
	}

	/* Free resources. */
	FreeArray( inds );
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Private Functions
*/
