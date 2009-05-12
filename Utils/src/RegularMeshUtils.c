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
** $Id: RegularMeshUtils.c 4184 2007-09-25 07:54:17Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <mpi.h>

#include <StGermain/StGermain.h>
#include <StgDomain/Geometry/Geometry.h>
#include <StgDomain/Shape/Shape.h>
#include <StgDomain/Mesh/Mesh.h>

#include "types.h"
#include "RegularMeshUtils.h"


Index RegularMeshUtils_ascendingIJK_ToHughesNodeNumberMap[8] = { 0, 1, 3, 2, 4, 5, 7, 6 };


/*----------------------------------------------------------------------------------------------------------------------------------
** Mapping functions
*/

void RegularMeshUtils_Node_1DTo3D( void* _mesh, unsigned global, unsigned* inds ) {
	Mesh*	mesh = (Mesh*)_mesh;
	Grid**	grid;

	assert( mesh );
	assert( global < Mesh_GetGlobalSize( mesh, MT_VERTEX ) );
	assert( inds );

	grid = (Grid**)ExtensionManager_Get( mesh->info, mesh, 
					     ExtensionManager_GetHandle( mesh->info, "vertexGrid" ) );
	Grid_Lift( *grid, global, inds );
}

unsigned RegularMeshUtils_Node_3DTo1D( void* _mesh, unsigned* inds ) {
	Mesh*		mesh = (Mesh*)_mesh;
	Grid**		grid;

	assert( mesh );
	assert( inds );

	grid = (Grid**)ExtensionManager_Get( mesh->info, mesh, 
					     ExtensionManager_GetHandle( mesh->info, "vertexGrid" ) );

	return Grid_Project( *grid, inds );
}

void RegularMeshUtils_Element_1DTo3D( void* _mesh, unsigned global, unsigned* inds ) {
	Mesh*	mesh = (Mesh*)_mesh;
	Grid**	grid;

	assert( mesh );
	assert( global < Mesh_GetGlobalSize( mesh, MT_VERTEX ) );
	assert( inds );

	grid = (Grid**)ExtensionManager_Get( mesh->info, mesh, 
					     ExtensionManager_GetHandle( mesh->info, "elementGrid" ) );
	Grid_Lift( *grid, global, inds );
}

unsigned RegularMeshUtils_Element_3DTo1D( void* _mesh, unsigned* inds ) {
	Mesh*		mesh = (Mesh*)_mesh;
	Grid**		grid;

	assert( mesh );
	assert( inds );

	grid = (Grid**)ExtensionManager_Get( mesh->info, mesh, 
					     ExtensionManager_GetHandle( mesh->info, "elementGrid" ) );

	return Grid_Project( *grid, inds );
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Set functions
*/

IndexSet* RegularMeshUtils_CreateGlobalTopSet( void* _mesh ) {
	Mesh*		mesh = (Mesh*)_mesh;
	Grid*		grid;
	unsigned	nNodes;
	IndexSet*	set;
	IJK		ijk;
	unsigned	n_i;

	assert( mesh );
	assert( Mesh_GetDimSize( mesh ) >= 2 );

	grid = *(Grid**)ExtensionManager_Get( mesh->info, mesh, 
					      ExtensionManager_GetHandle( mesh->info, "vertexGrid" ) );

	nNodes = Mesh_GetDomainSize( mesh, MT_VERTEX );
	set = IndexSet_New( nNodes );

	for( n_i = 0; n_i < nNodes; n_i++ ) {
		RegularMeshUtils_Node_1DTo3D( mesh, Mesh_DomainToGlobal( mesh, MT_VERTEX, n_i ), ijk );
		if( ijk[1] == grid->sizes[1] - 1 )
			IndexSet_Add( set, n_i );
	}

	return set;
}

IndexSet* RegularMeshUtils_CreateGlobalBottomSet( void* _mesh ) {
	Mesh*		mesh = (Mesh*)_mesh;
	unsigned	nNodes;
	IndexSet*	set;
	IJK		ijk;
	unsigned	n_i;

	assert( mesh );
	assert( Mesh_GetDimSize( mesh ) >= 2 );

	nNodes = Mesh_GetDomainSize( mesh, MT_VERTEX );
	set = IndexSet_New( nNodes );

	for( n_i = 0; n_i < nNodes; n_i++ ) {
		RegularMeshUtils_Node_1DTo3D( mesh, Mesh_DomainToGlobal( mesh, MT_VERTEX, n_i ), ijk );
		if( ijk[1] == 0 )
			IndexSet_Add( set, n_i );
	}

	return set;
}

IndexSet* RegularMeshUtils_CreateGlobalLeftSet( void* _mesh ) {
	Mesh*		mesh = (Mesh*)_mesh;
	unsigned	nNodes;
	IndexSet*	set;
	IJK		ijk;
	unsigned	n_i;

	assert( mesh );

	nNodes = Mesh_GetDomainSize( mesh, MT_VERTEX );
	set = IndexSet_New( nNodes );

	for( n_i = 0; n_i < nNodes; n_i++ ) {
		RegularMeshUtils_Node_1DTo3D( mesh, Mesh_DomainToGlobal( mesh, MT_VERTEX, n_i ), ijk );
		if( ijk[0] == 0 )
			IndexSet_Add( set, n_i );
	}

	return set;
}

IndexSet* RegularMeshUtils_CreateGlobalRightSet( void* _mesh ) {
	Mesh*		mesh = (Mesh*)_mesh;
	Grid*		grid;
	unsigned	nNodes;
	IndexSet*	set;
	IJK		ijk;
	unsigned	n_i;

	assert( mesh );

	grid = *(Grid**)ExtensionManager_Get( mesh->info, mesh, 
					      ExtensionManager_GetHandle( mesh->info, "vertexGrid" ) );

	nNodes = Mesh_GetDomainSize( mesh, MT_VERTEX );
	set = IndexSet_New( nNodes );

	for( n_i = 0; n_i < nNodes; n_i++ ) {
		RegularMeshUtils_Node_1DTo3D( mesh, Mesh_DomainToGlobal( mesh, MT_VERTEX, n_i ), ijk );
		if( ijk[0] == grid->sizes[0] - 1 )
			IndexSet_Add( set, n_i );
	}

	return set;
}

IndexSet* RegularMeshUtils_CreateGlobalFrontSet( void* _mesh ) {
	Mesh*		mesh = (Mesh*)_mesh;
	Grid*		grid;
	unsigned	nNodes;
	IndexSet*	set;
	IJK		ijk;
	unsigned	n_i;

	assert( mesh );
	assert( Mesh_GetDimSize( mesh ) >= 3 );

	grid = *(Grid**)ExtensionManager_Get( mesh->info, mesh, 
					      ExtensionManager_GetHandle( mesh->info, "vertexGrid" ) );

	nNodes = Mesh_GetDomainSize( mesh, MT_VERTEX );
	set = IndexSet_New( nNodes );

	for( n_i = 0; n_i < nNodes; n_i++ ) {
		RegularMeshUtils_Node_1DTo3D( mesh, Mesh_DomainToGlobal( mesh, MT_VERTEX, n_i ), ijk );
		if( ijk[2] == grid->sizes[2] - 1 )
			IndexSet_Add( set, n_i );
	}

	return set;
}

IndexSet* RegularMeshUtils_CreateGlobalBackSet( void* _mesh ) {
	Mesh*		mesh = (Mesh*)_mesh;
	unsigned	nNodes;
	IndexSet*	set;
	IJK		ijk;
	unsigned	n_i;

	assert( mesh );
	assert( Mesh_GetDimSize( mesh ) >= 3 );

	nNodes = Mesh_GetDomainSize( mesh, MT_VERTEX );
	set = IndexSet_New( nNodes );

	for( n_i = 0; n_i < nNodes; n_i++ ) {
		RegularMeshUtils_Node_1DTo3D( mesh, Mesh_DomainToGlobal( mesh, MT_VERTEX, n_i ), ijk );
		if( ijk[2] == 0 )
			IndexSet_Add( set, n_i );
	}

	return set;
}

IndexSet* RegularMeshUtils_CreateGlobalInnerTopSet( void* _mesh ) {
	Mesh*		mesh = (Mesh*)_mesh;
	Grid*		grid;
	unsigned	nNodes;
	IndexSet*	set;
	IJK		ijk;
	unsigned	n_i;
        int nDims;

	assert( mesh );
	assert( Mesh_GetDimSize( mesh ) >= 2 );

	grid = *(Grid**)ExtensionManager_Get( mesh->info, mesh, 
					      ExtensionManager_GetHandle( mesh->info, "vertexGrid" ) );

        nDims = Mesh_GetDimSize( mesh );
	nNodes = Mesh_GetDomainSize( mesh, MT_VERTEX );
	set = IndexSet_New( nNodes );

	for( n_i = 0; n_i < nNodes; n_i++ ) {
		RegularMeshUtils_Node_1DTo3D( mesh, Mesh_DomainToGlobal( mesh, MT_VERTEX, n_i ), ijk );
		if( ijk[1] == grid->sizes[1] - 1 && 
		    (ijk[0] != grid->sizes[0] - 1 || (nDims == 3 && ijk[2] != grid->sizes[2] - 1)) && 
		    (ijk[0] != 0 || (nDims == 3 && ijk[2] != grid->sizes[2] - 1)) && 
		    (ijk[0] != grid->sizes[0] - 1 || (nDims == 3 && ijk[2] != 0)) && 
		    (ijk[0] != 0 || (nDims == 3 && ijk[2] != 0)) )
		{
			IndexSet_Add( set, n_i );
		}
	}

	return set;
}

IndexSet* RegularMeshUtils_CreateGlobalInnerBottomSet( void* _mesh ) {
	Mesh*		mesh = (Mesh*)_mesh;
	Grid*		grid;
	unsigned	nNodes;
	IndexSet*	set;
	IJK		ijk;
	unsigned	n_i;
        int nDims;

	assert( mesh );
	assert( Mesh_GetDimSize( mesh ) >= 2 );

	grid = *(Grid**)ExtensionManager_Get( mesh->info, mesh, 
					      ExtensionManager_GetHandle( mesh->info, "vertexGrid" ) );

        nDims = Mesh_GetDimSize( mesh );
	nNodes = Mesh_GetDomainSize( mesh, MT_VERTEX );
	set = IndexSet_New( nNodes );

	for( n_i = 0; n_i < nNodes; n_i++ ) {
		RegularMeshUtils_Node_1DTo3D( mesh, Mesh_DomainToGlobal( mesh, MT_VERTEX, n_i ), ijk );
		if(ijk[1] == 0 &&
		   (ijk[0] != grid->sizes[0] - 1 && ijk[0] != 0) &&
		   (nDims != 3 || (ijk[2] != grid->sizes[2] - 1 && ijk[2] != 0)))
		{
		    IndexSet_Add( set, n_i );
		}
/* 		if( ijk[1] == 0 &&  */
/* 		    (ijk[0] != grid->sizes[0] - 1 || (nDims == 3 && ijk[2] != grid->sizes[2] - 1)) &&  */
/* 		    (ijk[0] != 0 || (nDims == 3 && ijk[2] != grid->sizes[2] - 1)) &&  */
/* 		    (ijk[0] != grid->sizes[0] - 1 || (nDims == 3 && ijk[2] != 0)) &&  */
/* 		    (ijk[0] != 0 || (nDims == 3 && ijk[2] != 0)) ) */
/* 		{ */
/* 			IndexSet_Add( set, n_i ); */
/* 		} */
	}

	return set;
}

IndexSet* RegularMeshUtils_CreateGlobalInnerLeftSet( void* _mesh ) {
	Mesh*		mesh = (Mesh*)_mesh;
	Grid*		grid;
	unsigned	nNodes;
	IndexSet*	set;
	IJK		ijk;
	unsigned	n_i;
        int nDims;

	assert( mesh );
	assert( Mesh_GetDimSize( mesh ) >= 2 );

	grid = *(Grid**)ExtensionManager_Get( mesh->info, mesh, 
					      ExtensionManager_GetHandle( mesh->info, "vertexGrid" ) );

        nDims = Mesh_GetDimSize( mesh );
	nNodes = Mesh_GetDomainSize( mesh, MT_VERTEX );
	set = IndexSet_New( nNodes );

	for( n_i = 0; n_i < nNodes; n_i++ ) {
		RegularMeshUtils_Node_1DTo3D( mesh, Mesh_DomainToGlobal( mesh, MT_VERTEX, n_i ), ijk );
		if( ijk[0] == 0 && 
		    (ijk[1] != grid->sizes[1] - 1 || (nDims == 3 && ijk[2] != grid->sizes[2] - 1)) && 
		    (ijk[1] != 0 || (nDims == 3 && ijk[2] != grid->sizes[2] - 1)) && 
		    (ijk[1] != grid->sizes[1] - 1 || (nDims == 3 && ijk[2] != 0)) && 
		    (ijk[1] != 0 || (nDims == 3 && ijk[2] != 0)) )
		{
			IndexSet_Add( set, n_i );
		}
	}

	return set;
}

IndexSet* RegularMeshUtils_CreateGlobalInnerRightSet( void* _mesh ) {
	Mesh*		mesh = (Mesh*)_mesh;
	Grid*		grid;
	unsigned	nNodes;
	IndexSet*	set;
	IJK		ijk;
	unsigned	n_i;
        int nDims;

	assert( mesh );
	assert( Mesh_GetDimSize( mesh ) >= 2 );

	grid = *(Grid**)ExtensionManager_Get( mesh->info, mesh, 
					      ExtensionManager_GetHandle( mesh->info, "vertexGrid" ) );

        nDims = Mesh_GetDimSize( mesh );
	nNodes = Mesh_GetDomainSize( mesh, MT_VERTEX );
	set = IndexSet_New( nNodes );

	for( n_i = 0; n_i < nNodes; n_i++ ) {
		RegularMeshUtils_Node_1DTo3D( mesh, Mesh_DomainToGlobal( mesh, MT_VERTEX, n_i ), ijk );
		if( ijk[0] == grid->sizes[0] - 1 && 
		    (ijk[1] != grid->sizes[1] - 1 || (nDims == 3 && ijk[2] != grid->sizes[2] - 1)) && 
		    (ijk[1] != 0 || (nDims == 3 && ijk[2] != grid->sizes[2] - 1)) && 
		    (ijk[1] != grid->sizes[1] - 1 || (nDims == 3 && ijk[2] != 0)) && 
		    (ijk[1] != 0 || (nDims == 3 && ijk[2] != 0)) )
		{
			IndexSet_Add( set, n_i );
		}
	}

	return set;
}

IndexSet* RegularMeshUtils_CreateGlobalInnerFrontSet( void* _mesh ) {
	Mesh*		mesh = (Mesh*)_mesh;
	Grid*		grid;
	unsigned	nNodes;
	IndexSet*	set;
	IJK		ijk;
	unsigned	n_i;

	assert( mesh );
	assert( Mesh_GetDimSize( mesh ) >= 2 );

	grid = *(Grid**)ExtensionManager_Get( mesh->info, mesh, 
					      ExtensionManager_GetHandle( mesh->info, "vertexGrid" ) );

	nNodes = Mesh_GetDomainSize( mesh, MT_VERTEX );
	set = IndexSet_New( nNodes );

	for( n_i = 0; n_i < nNodes; n_i++ ) {
		RegularMeshUtils_Node_1DTo3D( mesh, Mesh_DomainToGlobal( mesh, MT_VERTEX, n_i ), ijk );
		if( ijk[0] == grid->sizes[2] - 1 && 
		    (ijk[0] != grid->sizes[0] - 1 || ijk[1] != grid->sizes[1] - 1 ) && 
		    (ijk[0] != 0 || ijk[1] != grid->sizes[1] - 1 ) && 
		    (ijk[0] != grid->sizes[0] - 1 || ijk[1] != 0 ) && 
		    (ijk[0] != 0 || ijk[1] != 0 ) )
		{
			IndexSet_Add( set, n_i );
		}
	}

	return set;
}

IndexSet* RegularMeshUtils_CreateGlobalInnerBackSet( void* _mesh ) {
	Mesh*		mesh = (Mesh*)_mesh;
	Grid*		grid;
	unsigned	nNodes;
	IndexSet*	set;
	IJK		ijk;
	unsigned	n_i;

	assert( mesh );
	assert( Mesh_GetDimSize( mesh ) >= 2 );

	grid = *(Grid**)ExtensionManager_Get( mesh->info, mesh, 
					      ExtensionManager_GetHandle( mesh->info, "vertexGrid" ) );

	nNodes = Mesh_GetDomainSize( mesh, MT_VERTEX );
	set = IndexSet_New( nNodes );

	for( n_i = 0; n_i < nNodes; n_i++ ) {
		RegularMeshUtils_Node_1DTo3D( mesh, Mesh_DomainToGlobal( mesh, MT_VERTEX, n_i ), ijk );
		if( ijk[0] == 0 && 
		    (ijk[0] != grid->sizes[0] - 1 || ijk[1] != grid->sizes[1] - 1 ) && 
		    (ijk[0] != 0 || ijk[1] != grid->sizes[1] - 1 ) && 
		    (ijk[0] != grid->sizes[0] - 1 || ijk[1] != 0 ) && 
		    (ijk[0] != 0 || ijk[1] != 0 ) )
		{
			IndexSet_Add( set, n_i );
		}
	}

	return set;
}

IndexSet* RegularMeshUtils_CreateGlobalBottomLeftFrontSet( void* _mesh ) {
	Mesh*		mesh = (Mesh*)_mesh;
	Grid*		grid;
	unsigned	nNodes;
	IndexSet*	set;
	IJK		ijk;
	unsigned	n_i;

	assert( mesh );
	assert( Mesh_GetDimSize( mesh ) >= 2 );

	grid = *(Grid**)ExtensionManager_Get( mesh->info, mesh, 
					      ExtensionManager_GetHandle( mesh->info, "vertexGrid" ) );

	nNodes = Mesh_GetDomainSize( mesh, MT_VERTEX );
	set = IndexSet_New( nNodes );

	for( n_i = 0; n_i < nNodes; n_i++ ) {
		RegularMeshUtils_Node_1DTo3D( mesh, Mesh_DomainToGlobal( mesh, MT_VERTEX, n_i ), ijk );
		if( ijk[0] == 0 && 
		    ijk[1] == 0 && 
		    ijk[2] == grid->sizes[2] - 1 )
		{
			IndexSet_Add( set, n_i );
		}
	}

	return set;
}

IndexSet* RegularMeshUtils_CreateGlobalBottomRightFrontSet( void* _mesh ) {
	Mesh*		mesh = (Mesh*)_mesh;
	Grid*		grid;
	unsigned	nNodes;
	IndexSet*	set;
	IJK		ijk;
	unsigned	n_i;

	assert( mesh );
	assert( Mesh_GetDimSize( mesh ) >= 2 );

	grid = *(Grid**)ExtensionManager_Get( mesh->info, mesh, 
					      ExtensionManager_GetHandle( mesh->info, "vertexGrid" ) );

	nNodes = Mesh_GetDomainSize( mesh, MT_VERTEX );
	set = IndexSet_New( nNodes );

	for( n_i = 0; n_i < nNodes; n_i++ ) {
		RegularMeshUtils_Node_1DTo3D( mesh, Mesh_DomainToGlobal( mesh, MT_VERTEX, n_i ), ijk );
		if( ijk[0] == grid->sizes[0] - 1 && 
		    ijk[1] == 0 && 
		    ijk[2] == grid->sizes[2] - 1 )
		{
			IndexSet_Add( set, n_i );
		}
	}

	return set;
}

IndexSet* RegularMeshUtils_CreateGlobalTopLeftFrontSet( void* _mesh ) {
	Mesh*		mesh = (Mesh*)_mesh;
	Grid*		grid;
	unsigned	nNodes;
	IndexSet*	set;
	IJK		ijk;
	unsigned	n_i;
	int nDims;

	assert( mesh );
	assert( Mesh_GetDimSize( mesh ) >= 2 );

	grid = *(Grid**)ExtensionManager_Get( mesh->info, mesh, 
					      ExtensionManager_GetHandle( mesh->info, "vertexGrid" ) );

        nDims = Mesh_GetDimSize( mesh );
	nNodes = Mesh_GetDomainSize( mesh, MT_VERTEX );
	set = IndexSet_New( nNodes );

	for( n_i = 0; n_i < nNodes; n_i++ ) {
		RegularMeshUtils_Node_1DTo3D( mesh, Mesh_DomainToGlobal( mesh, MT_VERTEX, n_i ), ijk );
		if( ijk[0] == 0 && 
		    ijk[1] == grid->sizes[1] - 1 && 
		    (nDims != 3 || ijk[2] == grid->sizes[2] - 1) )
		{
			IndexSet_Add( set, n_i );
		}
	}

	return set;
}

IndexSet* RegularMeshUtils_CreateGlobalTopRightFrontSet( void* _mesh ) {
	Mesh*		mesh = (Mesh*)_mesh;
	Grid*		grid;
	unsigned	nNodes;
	IndexSet*	set;
	IJK		ijk;
	unsigned	n_i;
	int nDims;

	assert( mesh );
	assert( Mesh_GetDimSize( mesh ) >= 2 );

	grid = *(Grid**)ExtensionManager_Get( mesh->info, mesh, 
					      ExtensionManager_GetHandle( mesh->info, "vertexGrid" ) );

        nDims = Mesh_GetDimSize( mesh );
	nNodes = Mesh_GetDomainSize( mesh, MT_VERTEX );
	set = IndexSet_New( nNodes );

	for( n_i = 0; n_i < nNodes; n_i++ ) {
		RegularMeshUtils_Node_1DTo3D( mesh, Mesh_DomainToGlobal( mesh, MT_VERTEX, n_i ), ijk );
		if( ijk[0] == grid->sizes[0] - 1 && 
		    ijk[1] == grid->sizes[1] - 1 && 
		    (nDims != 3 || ijk[2] == grid->sizes[2] - 1) )
		{
			IndexSet_Add( set, n_i );
		}
	}

	return set;
}

IndexSet* RegularMeshUtils_CreateGlobalBottomLeftBackSet( void* _mesh ) {
	Mesh*		mesh = (Mesh*)_mesh;
	Grid*		grid;
	unsigned	nNodes;
	IndexSet*	set;
	int nDims;
	IJK		ijk;
	unsigned	n_i;

	assert( mesh );
	assert( Mesh_GetDimSize( mesh ) >= 2 );

	grid = *(Grid**)ExtensionManager_Get( mesh->info, mesh, 
					      ExtensionManager_GetHandle( mesh->info, "vertexGrid" ) );

        nDims = Mesh_GetDimSize( mesh );
	nNodes = Mesh_GetDomainSize( mesh, MT_VERTEX );
	set = IndexSet_New( nNodes );

	for( n_i = 0; n_i < nNodes; n_i++ ) {
		RegularMeshUtils_Node_1DTo3D( mesh, Mesh_DomainToGlobal( mesh, MT_VERTEX, n_i ), ijk );
		if( ijk[0] == 0 && 
		    ijk[1] == 0 && 
		    (nDims != 3 || ijk[2] == 0))
		{
			IndexSet_Add( set, n_i );
		}
	}

	return set;
}

IndexSet* RegularMeshUtils_CreateGlobalBottomRightBackSet( void* _mesh ) {
	Mesh*		mesh = (Mesh*)_mesh;
	Grid*		grid;
	unsigned	nNodes;
	int nDims;
	IndexSet*	set;
	IJK		ijk;
	unsigned	n_i;

	assert( mesh );
	assert( Mesh_GetDimSize( mesh ) >= 2 );

	grid = *(Grid**)ExtensionManager_Get( mesh->info, mesh, 
					      ExtensionManager_GetHandle( mesh->info, "vertexGrid" ) );

        nDims = Mesh_GetDimSize( mesh );
	nNodes = Mesh_GetDomainSize( mesh, MT_VERTEX );
	set = IndexSet_New( nNodes );

	for( n_i = 0; n_i < nNodes; n_i++ ) {
		RegularMeshUtils_Node_1DTo3D( mesh, Mesh_DomainToGlobal( mesh, MT_VERTEX, n_i ), ijk );
		if( ijk[0] == grid->sizes[0] - 1 && 
		    ijk[1] == 0 && 
		    (nDims != 3 || ijk[2] == 0) )
		{
			IndexSet_Add( set, n_i );
		}
	}

	return set;
}

IndexSet* RegularMeshUtils_CreateGlobalTopLeftBackSet( void* _mesh ) {
	Mesh*		mesh = (Mesh*)_mesh;
	Grid*		grid;
	unsigned	nNodes;
	IndexSet*	set;
	IJK		ijk;
	unsigned	n_i;
	int nDims;

	assert( mesh );
	assert( Mesh_GetDimSize( mesh ) >= 2 );

	grid = *(Grid**)ExtensionManager_Get( mesh->info, mesh, 
					      ExtensionManager_GetHandle( mesh->info, "vertexGrid" ) );

        nDims = Mesh_GetDimSize( mesh );
	nNodes = Mesh_GetDomainSize( mesh, MT_VERTEX );
	set = IndexSet_New( nNodes );

	for( n_i = 0; n_i < nNodes; n_i++ ) {
		RegularMeshUtils_Node_1DTo3D( mesh, Mesh_DomainToGlobal( mesh, MT_VERTEX, n_i ), ijk );
		if( ijk[0] == 0 && 
		    ijk[1] == grid->sizes[1] - 1 && 
		    (nDims != 3 || ijk[2] == 0) )
		{
			IndexSet_Add( set, n_i );
		}
	}

	return set;
}

IndexSet* RegularMeshUtils_CreateGlobalTopRightBackSet( void* _mesh ) {
	Mesh*		mesh = (Mesh*)_mesh;
	Grid*		grid;
	unsigned	nNodes;
	IndexSet*	set;
	IJK		ijk;
	unsigned	n_i;
	int nDims;

	assert( mesh );
	assert( Mesh_GetDimSize( mesh ) >= 2 );

	grid = *(Grid**)ExtensionManager_Get( mesh->info, mesh, 
					      ExtensionManager_GetHandle( mesh->info, "vertexGrid" ) );

        nDims = Mesh_GetDimSize( mesh );
	nNodes = Mesh_GetDomainSize( mesh, MT_VERTEX );
	set = IndexSet_New( nNodes );

	for( n_i = 0; n_i < nNodes; n_i++ ) {
		RegularMeshUtils_Node_1DTo3D( mesh, Mesh_DomainToGlobal( mesh, MT_VERTEX, n_i ), ijk );
		if( ijk[0] == grid->sizes[0] - 1 && 
		    ijk[1] == grid->sizes[1] - 1 && 
		    (nDims != 3 || ijk[2] == 0) )
		{
			IndexSet_Add( set, n_i );
		}
	}

	return set;
}

IndexSet* RegularMeshUtils_CreateLocalInGlobalTopSet( void* _mesh ) {
	Mesh*		mesh = (Mesh*)_mesh;
	Grid*		grid;
	unsigned	nNodes;
	IndexSet*	set;
	IJK		ijk;
	unsigned	n_i;

	assert( mesh );
	assert( Mesh_GetDimSize( mesh ) >= 2 );

	grid = *(Grid**)ExtensionManager_Get( mesh->info, mesh, 
					      ExtensionManager_GetHandle( mesh->info, "vertexGrid" ) );

	nNodes = Mesh_GetLocalSize( mesh, MT_VERTEX );
	set = IndexSet_New( nNodes );

	for( n_i = 0; n_i < nNodes; n_i++ ) {
		RegularMeshUtils_Node_1DTo3D( mesh, Mesh_DomainToGlobal( mesh, MT_VERTEX, n_i ), ijk );
		if( ijk[1] == grid->sizes[1] - 1 )
			IndexSet_Add( set, n_i );
	}

	return set;
}

IndexSet* RegularMeshUtils_CreateLocalInGlobalBottomSet( void* _mesh ) {
	Mesh*		mesh = (Mesh*)_mesh;
	unsigned	nNodes;
	IndexSet*	set;
	IJK		ijk;
	unsigned	n_i;

	assert( mesh );
	assert( Mesh_GetDimSize( mesh ) >= 2 );

	nNodes = Mesh_GetLocalSize( mesh, MT_VERTEX );
	set = IndexSet_New( nNodes );

	for( n_i = 0; n_i < nNodes; n_i++ ) {
		RegularMeshUtils_Node_1DTo3D( mesh, Mesh_DomainToGlobal( mesh, MT_VERTEX, n_i ), ijk );
		if( ijk[1] == 0 )
			IndexSet_Add( set, n_i );
	}

	return set;
}

IndexSet* RegularMeshUtils_CreateLocalInGlobalLeftSet( void* _mesh ) {
	Mesh*		mesh = (Mesh*)_mesh;
	unsigned	nNodes;
	IndexSet*	set;
	IJK		ijk;
	unsigned	n_i;

	assert( mesh );

	nNodes = Mesh_GetLocalSize( mesh, MT_VERTEX );
	set = IndexSet_New( nNodes );

	for( n_i = 0; n_i < nNodes; n_i++ ) {
		RegularMeshUtils_Node_1DTo3D( mesh, Mesh_DomainToGlobal( mesh, MT_VERTEX, n_i ), ijk );
		if( ijk[0] == 0 )
			IndexSet_Add( set, n_i );
	}

	return set;
}

IndexSet* RegularMeshUtils_CreateLocalInGlobalRightSet( void* _mesh ) {
	Mesh*		mesh = (Mesh*)_mesh;
	Grid*		grid;
	unsigned	nNodes;
	IndexSet*	set;
	IJK		ijk;
	unsigned	n_i;

	assert( mesh );

	grid = *(Grid**)ExtensionManager_Get( mesh->info, mesh, 
					      ExtensionManager_GetHandle( mesh->info, "vertexGrid" ) );

	nNodes = Mesh_GetLocalSize( mesh, MT_VERTEX );
	set = IndexSet_New( nNodes );

	for( n_i = 0; n_i < nNodes; n_i++ ) {
		RegularMeshUtils_Node_1DTo3D( mesh, Mesh_DomainToGlobal( mesh, MT_VERTEX, n_i ), ijk );
		if( ijk[0] == grid->sizes[0] - 1 )
			IndexSet_Add( set, n_i );
	}

	return set;
}

IndexSet* RegularMeshUtils_CreateLocalInGlobalFrontSet( void* _mesh ) {
	Mesh*		mesh = (Mesh*)_mesh;
	Grid*		grid;
	unsigned	nNodes;
	IndexSet*	set;
	IJK		ijk;
	unsigned	n_i;

	assert( mesh );
	assert( Mesh_GetDimSize( mesh ) >= 3 );

	grid = *(Grid**)ExtensionManager_Get( mesh->info, mesh, 
					      ExtensionManager_GetHandle( mesh->info, "vertexGrid" ) );

	nNodes = Mesh_GetLocalSize( mesh, MT_VERTEX );
	set = IndexSet_New( nNodes );

	for( n_i = 0; n_i < nNodes; n_i++ ) {
		RegularMeshUtils_Node_1DTo3D( mesh, Mesh_DomainToGlobal( mesh, MT_VERTEX, n_i ), ijk );
		if( ijk[2] == grid->sizes[2] - 1 )
			IndexSet_Add( set, n_i );
	}

	return set;
}

IndexSet* RegularMeshUtils_CreateLocalInGlobalBackSet( void* _mesh ) {
	Mesh*		mesh = (Mesh*)_mesh;
	unsigned	nNodes;
	IndexSet*	set;
	IJK		ijk;
	unsigned	n_i;

	assert( mesh );
	assert( Mesh_GetDimSize( mesh ) >= 3 );

	nNodes = Mesh_GetLocalSize( mesh, MT_VERTEX );
	set = IndexSet_New( nNodes );

	for( n_i = 0; n_i < nNodes; n_i++ ) {
		RegularMeshUtils_Node_1DTo3D( mesh, Mesh_DomainToGlobal( mesh, MT_VERTEX, n_i ), ijk );
		if( ijk[2] == 0 )
			IndexSet_Add( set, n_i );
	}

	return set;
}

IndexSet* RegularMeshUtils_CreateContactTopSet( void* _mesh, int lowDepth, int uppDepth ) {
   Mesh* mesh = (Mesh*)_mesh;
   Grid* grid;
   int nNodes;
   IndexSet* set;
   int ijk[2], left, right, bottom, top;
   int ii;

   assert( mesh );
   assert( Mesh_GetDimSize( mesh ) == 2 );

   grid = *Mesh_GetExtension( mesh, Grid**, "vertexGrid" );
   nNodes = Mesh_GetDomainSize( mesh, 0 );
   set = IndexSet_New( nNodes );

/*
   left = depth + 1;
   right = grid->sizes[0] - (depth + 1) - 1;
   bottom = 1;
   top = depth;
*/
   left = lowDepth;
   right = grid->sizes[0] - 1 - uppDepth;
   bottom = grid->sizes[1] - 1;
   top = grid->sizes[1] - 1;
   for( ii = 0; ii < nNodes; ii++ ) {
      Grid_Lift( grid, Mesh_DomainToGlobal( mesh, 0, ii ), ijk );
      if( ijk[0] >= left && ijk[0] <= right && ijk[1] >= bottom && ijk[1] <= top )
	 IndexSet_Add( set, ii );
   }

   return set;
}

IndexSet* RegularMeshUtils_CreateContactBottomSet( void* _mesh, int lowDepth, int uppDepth, int inDepth ) {
   Mesh* mesh = (Mesh*)_mesh;
   Grid* grid;
   int nNodes;
   IndexSet* set;
   int ijk[2], left, right, bottom, top;
   int ii;

   assert( mesh );
   assert( Mesh_GetDimSize( mesh ) == 2 );

   grid = *Mesh_GetExtension( mesh, Grid**, "vertexGrid" );
   nNodes = Mesh_GetDomainSize( mesh, 0 );
   set = IndexSet_New( nNodes );

/*
   left = depth + 1;
   right = grid->sizes[0] - (depth + 1) - 1;
   bottom = 1;
   top = depth;
*/
   left = lowDepth;
   right = grid->sizes[0] - 1 - uppDepth;
   bottom = 0;
   top = inDepth;
   for( ii = 0; ii < nNodes; ii++ ) {
      Grid_Lift( grid, Mesh_DomainToGlobal( mesh, 0, ii ), ijk );
      if( ijk[0] >= left && ijk[0] <= right && ijk[1] >= bottom && ijk[1] <= top ) {
	 IndexSet_Add( set, ii );
      }
   }

   return set;
}

IndexSet* RegularMeshUtils_CreateContactLeftSet( void* _mesh, int lowDepth, int uppDepth ) {
   Mesh* mesh = (Mesh*)_mesh;
   Grid* grid;
   int nNodes;
   IndexSet* set;
   int ijk[2], left, right, bottom, top;
   int ii;

   assert( mesh );
   assert( Mesh_GetDimSize( mesh ) == 2 );

   grid = *Mesh_GetExtension( mesh, Grid**, "vertexGrid" );
   nNodes = Mesh_GetDomainSize( mesh, 0 );
   set = IndexSet_New( nNodes );

/*
   left = 1;
   right = depth;
   bottom = depth + 1;
   top = grid->sizes[1] - 2;
*/
   left = 0;
   right = 0;
   bottom = lowDepth;
   top = grid->sizes[1] - 1 - uppDepth;
   for( ii = 0; ii < nNodes; ii++ ) {
      Grid_Lift( grid, Mesh_DomainToGlobal( mesh, 0, ii ), ijk );
      if( ijk[0] >= left && ijk[0] <= right && ijk[1] >= bottom && ijk[1] <= top )
	 IndexSet_Add( set, ii );
   }

   return set;
}

IndexSet* RegularMeshUtils_CreateContactRightSet( void* _mesh, int lowDepth, int uppDepth, int inDepth ) {
   Mesh* mesh = (Mesh*)_mesh;
   Grid* grid;
   int nNodes;
   IndexSet* set;
   int ijk[2], left, right, bottom, top;
   int ii;

   assert( mesh );
   assert( Mesh_GetDimSize( mesh ) == 2 );

   grid = *Mesh_GetExtension( mesh, Grid**, "vertexGrid" );
   nNodes = Mesh_GetDomainSize( mesh, 0 );
   set = IndexSet_New( nNodes );

/*
   left = grid->sizes[0] - depth - 1;
   right = grid->sizes[0] - 2;
   bottom = depth + 1;
   top = grid->sizes[1] - 2;
*/
   left = grid->sizes[0] - 1 - inDepth;
   right = grid->sizes[0] - 1;
   bottom = lowDepth;
   top = grid->sizes[1] - 1- uppDepth;
   for( ii = 0; ii < nNodes; ii++ ) {
      Grid_Lift( grid, Mesh_DomainToGlobal( mesh, 0, ii ), ijk );
      if( ijk[0] >= left && ijk[0] <= right && ijk[1] >= bottom && ijk[1] <= top )
	 IndexSet_Add( set, ii );
   }

   return set;
}

Node_DomainIndex RegularMeshUtils_GetDiagOppositeAcrossElementNodeIndex( void* _mesh, 
									 Element_DomainIndex refElement_dI, 
									 Node_DomainIndex refNode_dI )
{
	Mesh*              mesh = (Mesh*)_mesh;
	const Node_Index   oppositeNodesMap2D[] = { 3, 2, 1, 0 };
	Node_Index         oppositeNodesMap3D[] = { 7, 6, 5, 4, 3, 2, 1, 0 };
	Node_DomainIndex*  currElementNodes = NULL;
	Node_Index         currElementNodeCount = 0;
	Node_Index         refNode_eI = 0;
	Node_DomainIndex   oppositeNode_dI = 0;
	Node_Index         oppositeNode_eI = 0;
	Stream*            errorStrm = Journal_Register( Error_Type, "RegularMeshUtils" );
	IArray*		   inc;

	Journal_Firewall( Mesh_GetElementType( mesh, refElement_dI )->type == Mesh_HexType_Type, errorStrm, 
			  "Error (%s:%s:%d):\n\tIncorrect element type (%s); require %s.\n", 
			  __func__, __FILE__, __LINE__, Mesh_GetElementType( mesh, refElement_dI )->type, 
			  Mesh_HexType_Type );

#if 0
	Journal_Firewall( CornerNL_Type == mesh->layout->nodeLayout->type , errorStr,
		"Error- in %s: Given mesh has node layout of type \"%s\", different to "
		"required type \"%s\".\n", __func__, mesh->layout->nodeLayout->type, CornerNL_Type );
#endif

	inc = IArray_New();
	Mesh_GetIncidence( mesh, Mesh_GetDimSize( mesh ), refElement_dI, MT_VERTEX, 
			   inc );
	currElementNodeCount = IArray_GetSize( inc );
	currElementNodes = IArray_GetPtr( inc );

	/* Find index of reference node within reference element */
	for( refNode_eI = 0; refNode_eI < currElementNodeCount; refNode_eI++ ) {
		if ( refNode_dI == currElementNodes[refNode_eI] )
			break;
	}
	Journal_Firewall( refNode_eI < currElementNodeCount, errorStrm,
		"Error - in %s(): Reference node %d (domain) not found within reference element %d (domain).\n",
		__func__, refNode_dI, refElement_dI );

	/* Use mapping table to get diagonally opposite node, then convert back to domain index */
	
	if ( Mesh_GetDimSize( mesh ) == 2 ) {
		oppositeNode_eI = oppositeNodesMap2D[refNode_eI];
	}
	else {
		oppositeNode_eI = oppositeNodesMap3D[refNode_eI];
	}

	oppositeNode_dI = currElementNodes[oppositeNode_eI];
	NewClass_Delete( inc );
	return oppositeNode_dI;
}
