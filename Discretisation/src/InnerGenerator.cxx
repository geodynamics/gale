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
** $Id: InnerGenerator.c 3584 2006-05-16 11:11:07Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <mpi.h>

#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include "Discretisation.h"


/* Textual name of this class */
const Type InnerGenerator_Type = "InnerGenerator";


/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/

InnerGenerator* InnerGenerator_New( Name name, AbstractContext* context ) {
	/* Variables set in this function */
	SizeT                                              _sizeOfSelf = sizeof(InnerGenerator);
	Type                                                      type = InnerGenerator_Type;
	Stg_Class_DeleteFunction*                              _delete = _InnerGenerator_Delete;
	Stg_Class_PrintFunction*                                _print = _InnerGenerator_Print;
	Stg_Class_CopyFunction*                                  _copy = NULL;
	Stg_Component_DefaultConstructorFunction*  _defaultConstructor = (void* (*)(Name))_InnerGenerator_New;
	Stg_Component_ConstructFunction*                    _construct = _InnerGenerator_AssignFromXML;
	Stg_Component_BuildFunction*                            _build = _InnerGenerator_Build;
	Stg_Component_InitialiseFunction*                  _initialise = _InnerGenerator_Initialise;
	Stg_Component_ExecuteFunction*                        _execute = _InnerGenerator_Execute;
	Stg_Component_DestroyFunction*                        _destroy = NULL;
	AllocationType                              nameAllocationType = NON_GLOBAL;
	MeshGenerator_SetDimSizeFunc*                   setDimSizeFunc = _MeshGenerator_SetDimSize;
	MeshGenerator_GenerateFunc*                       generateFunc = (MeshGenerator_GenerateFunc*)InnerGenerator_Generate;

	InnerGenerator* self = _InnerGenerator_New(  INNERGENERATOR_PASSARGS  );

   _MeshGenerator_Init( (MeshGenerator*)self, context );
   _InnerGenerator_Init( self );

   return self;
}

InnerGenerator* _InnerGenerator_New(  INNERGENERATOR_DEFARGS  ) {
	InnerGenerator*	self;

	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(InnerGenerator) );
	self = (InnerGenerator*)_MeshGenerator_New(  MESHGENERATOR_PASSARGS  );

	return self;
}

void _InnerGenerator_Init( InnerGenerator* self ) {
	assert( self && Stg_CheckType( self, InnerGenerator ) );

	self->elMesh = NULL;
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _InnerGenerator_Delete( void* generator ) {
	InnerGenerator*	self = (InnerGenerator*)generator;

	/* Delete the parent. */
	_MeshGenerator_Delete( self );
}

void _InnerGenerator_Print( void* generator, Stream* stream ) {
	InnerGenerator*	self = (InnerGenerator*)generator;
	
	/* Print parent */
	Journal_Printf( stream, "InnerGenerator (ptr): (%p)\n", self );
	_MeshGenerator_Print( self, stream );
}

void _InnerGenerator_AssignFromXML( void* generator, Stg_ComponentFactory* cf, void* data ) {
	InnerGenerator*	self = (InnerGenerator*)generator;
	Mesh*		elMesh;

	assert( self );
	assert( cf );

	_MeshGenerator_AssignFromXML( self, cf, data );

	elMesh = Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"elementMesh", Mesh, True, data  );
	InnerGenerator_SetElementMesh( self, elMesh );
}

void _InnerGenerator_Build( void* generator, void* data ) {
	_MeshGenerator_Build( generator, data );
}

void _InnerGenerator_Initialise( void* generator, void* data ) {
	_MeshGenerator_Initialise( generator, data );
}

void _InnerGenerator_Execute( void* generator, void* data ) {
}

void _InnerGenerator_Destroy( void* generator, void* data ) {
	InnerGenerator*	self = (InnerGenerator*)generator;

   Stg_Component_Destroy( self->elMesh, data, False );
   _MeshGenerator_Destroy( self, data );
}

void InnerGenerator_Generate( void* generator, void* _mesh ) {
	InnerGenerator*	self = (InnerGenerator*)generator;
	FeMesh*		mesh = (FeMesh*)_mesh;
	Grid**		grid;
	Grid*		elGrid;

	assert( self && Stg_CheckType( self, InnerGenerator ) );
	assert( mesh && Stg_CheckType( mesh, FeMesh ) );

	InnerGenerator_BuildTopology( self, mesh );
	InnerGenerator_BuildGeometry( self, mesh );
	InnerGenerator_BuildElementTypes( self, mesh );

	elGrid = *(Grid**)ExtensionManager_Get( self->elMesh->info, self->elMesh, 
					       ExtensionManager_GetHandle( self->elMesh->info, (Name)"elementGrid" )  );
	ExtensionManager_Add( mesh->info, (Name)"elementGrid", sizeof(Grid*) );
	grid = (Grid** )ExtensionManager_Get( mesh->info, mesh, 
					     ExtensionManager_GetHandle( mesh->info, (Name)"elementGrid" ) );
	*grid = Grid_New( );
	Grid_SetNumDims( *grid, elGrid->nDims );
	Grid_SetSizes( *grid, elGrid->sizes );
}


/*--------------------------------------------------------------------------------------------------------------------------
** Public Functions
*/

void InnerGenerator_SetElementMesh( void* generator, void* mesh ) {
	InnerGenerator*	self = (InnerGenerator*)generator;

	assert( self && Stg_CheckType( self, InnerGenerator ) );

	self->elMesh = (Mesh*)mesh;
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Private Functions
*/

void InnerGenerator_BuildTopology( InnerGenerator* self, FeMesh* mesh ) {
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
	unsigned	e_i;
        int             l_i, r_i;

	assert( self );
	assert( mesh );

	elMesh = self->elMesh;
	nDims = Mesh_GetDimSize( elMesh );
	elTopo = Mesh_GetTopology( elMesh );
	elSync = Mesh_GetSync( elMesh, (MeshTopology_Dim)nDims );

	topo = Mesh_GetTopology( mesh );
	MeshTopology_SetComm( topo, MeshTopology_GetComm( elTopo ) );
	MeshTopology_SetNumDims( topo, nDims );
	IGraph_SetDomain( topo, nDims, elSync );

	/* Need to redefine the nodes, nDims + 1 per parent element. */
	elDecomp = (Decomp*)Sync_GetDecomp( elSync );
	nodeDecomp = Decomp_New();

	if( nDims == 2 ) {
		nLocals = Decomp_GetNumLocals( elDecomp ) * 3;
		locals = MemArray( int, nLocals, InnerGenerator_Type );
		for( l_i = 0; l_i < Decomp_GetNumLocals( elDecomp ); l_i++ ) {
			global = Decomp_LocalToGlobal( elDecomp, l_i );
			locals[l_i * 3 + 0] = global * 3;
			locals[l_i * 3 + 1] = global * 3 + 1;
			locals[l_i * 3 + 2] = global * 3 + 2;
		}
	}
	else if( nDims == 3 ) {
		nLocals = Decomp_GetNumLocals( elDecomp ) * 4;
		locals = MemArray( int, nLocals, InnerGenerator_Type );
		for( l_i = 0; l_i < Decomp_GetNumLocals( elDecomp ); l_i++ ) {
			global = Decomp_LocalToGlobal( elDecomp, l_i );
			locals[l_i * 4 + 0] = global * 4;
			locals[l_i * 4 + 1] = global * 4 + 1;
			locals[l_i * 4 + 2] = global * 4 + 2;
			locals[l_i * 4 + 3] = global * 4 + 3;
		}
	}
	Decomp_SetLocals( nodeDecomp, nLocals, locals );
	MemFree( locals );

	nodeSync = Sync_New();
	Sync_SetComm( nodeSync, Sync_GetComm( elSync ) );
	Sync_SetDecomp( nodeSync, nodeDecomp );

	if( nDims == 2 ) {
		nRemotes = Sync_GetNumRemotes( elSync ) * 3;
		remotes = MemArray( int, nRemotes, InnerGenerator_Type );
		for( r_i = 0; r_i < Sync_GetNumRemotes( elSync ); r_i++ ) {
			global = Sync_RemoteToGlobal( elSync, r_i );
			remotes[r_i * 3 + 0] = global * 3;
			remotes[r_i * 3 + 1] = global * 3 + 1;
			remotes[r_i * 3 + 2] = global * 3 + 2;
		}
	}
	else if( nDims == 3 ) {
		nRemotes = Sync_GetNumRemotes( elSync ) * 4;
		remotes = MemArray( int, nRemotes, InnerGenerator_Type );
		for( r_i = 0; r_i < Sync_GetNumRemotes( elSync ); r_i++ ) {
			global = Sync_RemoteToGlobal( elSync, r_i );
			remotes[r_i * 4 + 0] = global * 4;
			remotes[r_i * 4 + 1] = global * 4 + 1;
			remotes[r_i * 4 + 2] = global * 4 + 2;
			remotes[r_i * 4 + 3] = global * 4 + 3;
		}
	}
	Sync_SetRemotes( nodeSync, nRemotes, remotes );
	MemFree( remotes );

	IGraph_SetDomain( topo, 0, nodeSync );

	/* Same shadow depth. */
	topo->shadDepth = elTopo->shadDepth;

	/* Build the incidence. */
	nDomainEls = Mesh_GetDomainSize( elMesh, (MeshTopology_Dim)nDims );
	if( nDims == 2 ) {
		nIncEls = 3;
		incEls = MemArray( unsigned, 3, InnerGenerator_Type );
		for( e_i = 0; e_i < nDomainEls; e_i++ ) {
			incEls[0] = e_i * 3;
			incEls[1] = e_i * 3 + 1;
			incEls[2] = e_i * 3 + 2;
			IGraph_SetIncidence( topo, nDims, e_i, 0, nIncEls, (int*)incEls );
		}
	}
	else if( nDims == 3 ) {
		nIncEls = 4;
		incEls = MemArray( unsigned, 4, InnerGenerator_Type );
		for( e_i = 0; e_i < nDomainEls; e_i++ ) {
			incEls[0] = e_i * 4;
			incEls[1] = e_i * 4 + 1;
			incEls[2] = e_i * 4 + 2;
			incEls[3] = e_i * 4 + 3;
			IGraph_SetIncidence( topo, nDims, e_i, 0, nIncEls, (int*)incEls );
		}
	}
	FreeArray( incEls );

	IGraph_InvertIncidence( topo, MT_VERTEX, nDims );
}

void InnerGenerator_SetCoordinates( InnerGenerator* self, FeMesh* mesh );

void InnerGenerator_BuildGeometry( InnerGenerator* self, FeMesh* mesh ) {
  assert( self );
  assert( mesh );

  Mesh *elMesh = self->elMesh;
  unsigned nDims = Mesh_GetDimSize( elMesh );
  unsigned nDomainEls = Mesh_GetDomainSize( elMesh, (MeshTopology_Dim)nDims );

  if( nDims == 2 ) {
    mesh->verts = AllocArray2D( double, nDomainEls * 3, nDims );
  }
  else if( nDims == 3 ) {
    mesh->verts = AllocArray2D( double, nDomainEls * 4, nDims );
  }
  InnerGenerator_SetCoordinates(self, mesh);
}

void InnerGenerator_SetCoordinates( InnerGenerator* self, FeMesh* mesh ) {
  assert( self );
  assert( mesh );

  Mesh *elMesh = self->elMesh;
  unsigned nDims = Mesh_GetDimSize( elMesh );
  unsigned nDomainEls = Mesh_GetDomainSize( elMesh, (MeshTopology_Dim)nDims );
  double *vert;

  if( nDims == 2 ) {
    double localCrds[3][2] = {{-0.5, -0.5}, 
                              {0.5, -0.5}, 
                              {0, 0.5}};
    double globalCrd[2];
    for(unsigned e_i = 0; e_i < nDomainEls; e_i++ ) {
      unsigned elInd = e_i * 3;

      FeMesh_CoordLocalToGlobal( elMesh, e_i, localCrds[0], globalCrd );
      vert = Mesh_GetVertex( mesh, elInd );
      memcpy( vert, globalCrd, nDims * sizeof(double) );

      FeMesh_CoordLocalToGlobal( elMesh, e_i, localCrds[1], globalCrd );
      vert = Mesh_GetVertex( mesh, elInd + 1 );
      memcpy( vert, globalCrd, nDims * sizeof(double) );

      FeMesh_CoordLocalToGlobal( elMesh, e_i, localCrds[2], globalCrd );
      vert = Mesh_GetVertex( mesh, elInd + 2 );
      memcpy( vert, globalCrd, nDims * sizeof(double) );
    }
  }

  else if( nDims == 3 ) {
    double localCrds3D[4][3] = { {-0.5, -0.5, -0.5},
                                 {0.5, -0.5, -0.5},
                                 {0.0, 0.5, -0.5},
                                 {0.0, 0.0, 0.5} };
    double globalCrd3D[3];

    for(unsigned e_i = 0; e_i < nDomainEls; e_i++ ) {
      unsigned elInd = e_i * 4;

      FeMesh_CoordLocalToGlobal( elMesh, e_i, localCrds3D[0], globalCrd3D );
      vert = Mesh_GetVertex( mesh, elInd );
      memcpy( vert, globalCrd3D, nDims * sizeof(double) );

      FeMesh_CoordLocalToGlobal( elMesh, e_i, localCrds3D[1], globalCrd3D );
      vert = Mesh_GetVertex( mesh, elInd + 1 );
      memcpy( vert, globalCrd3D, nDims * sizeof(double) );

      FeMesh_CoordLocalToGlobal( elMesh, e_i, localCrds3D[2], globalCrd3D );
      vert = Mesh_GetVertex( mesh, elInd + 2 );
      memcpy( vert, globalCrd3D, nDims * sizeof(double) );

      FeMesh_CoordLocalToGlobal( elMesh, e_i, localCrds3D[3], globalCrd3D );
      vert = Mesh_GetVertex( mesh, elInd + 3 );
      memcpy( vert, globalCrd3D, nDims * sizeof(double) );
    }
  }

}


void InnerGenerator_BuildElementTypes( InnerGenerator* self, FeMesh* mesh ) {
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

	algs = (Mesh_Algorithms*)Mesh_CentroidAlgorithms_New( "", NULL );
	Mesh_CentroidAlgorithms_SetElementMesh( algs, self->elMesh );
	Mesh_SetAlgorithms( mesh, algs );
}


