/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003, Victorian Partnership for Advanced Computing (VPAC) Ltd, 
** 110 Victoria Street, Melbourne, 3053, Australia.
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
** $Id: RegularRemesher.c 3952 2007-01-09 06:24:06Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdlib.h>
#include <string.h>
#include <StGermain/StGermain.h>
#include <StgDomain/Geometry/Geometry.h>
#include <StgDomain/Shape/Shape.h>
#include <StgDomain/Mesh/Mesh.h>

#include "types.h"
#include "NewRemesher.h"
#include "RegularRemesher.h"
#include "StGermain/Base/Foundation/ClassDef.h"


void _RegularRemesher_Init( void* _self ) {
   RegularRemesher* self = Class_Cast( _self, RegularRemesher );

   _NewRemesher_Init( self );
   self->remeshDims = ISet_New();
   self->staticWalls = Class_Array( self, ISet*, 2 );
   self->staticWalls[0] = ISet_New();
   self->staticWalls[1] = ISet_New();
   ISet_SetMaxSize( self->remeshDims, 10 );
   ISet_SetMaxSize( self->staticWalls[0], 10 );
   ISet_SetMaxSize( self->staticWalls[1], 10 );
   self->syncs = NULL;
   self->crds = NULL;
   self->nWallVerts = NULL;
   self->wallVerts = NULL;
   self->wallCrds = NULL;
}

void _RegularRemesher_Copy( void* _self, const void* _op ) {
   RegularRemesher* self = Class_Cast( _self, RegularRemesher );
   const RegularRemesher* op = Class_Cast( _op, RegularRemesher );

   _NewRemesher_Copy( self, op );
   /* TODO */
}

void _RegularRemesher_Destruct( void* _self ) {
   RegularRemesher* self = Class_Cast( _self, RegularRemesher );
   int nDims;
   int d_i;

   if( self->syncs ) {
      nDims = self->nDims;
      for( d_i = 0; d_i < nDims; d_i++ ) {
   	 NewClass_RemoveRef( self->syncs[d_i] );
	 Class_Free( self, self->crds[d_i] );
	 Class_Free( self, self->nWallVerts[d_i] );
	 Class_Free( self, self->wallVerts[d_i][0] );
	 Class_Free( self, self->wallVerts[d_i][1] );
	 Class_Free( self, self->wallVerts[d_i] );
	 Class_Free( self, self->wallCrds[d_i][0] );
	 Class_Free( self, self->wallCrds[d_i][1] );
	 Class_Free( self, self->wallCrds[d_i] );
      }
      Class_Free( self, self->syncs );
      Class_Free( self, self->crds );
      Class_Free( self, self->nWallVerts );
      Class_Free( self, self->wallVerts );
      Class_Free( self, self->wallCrds );
   }
   NewClass_Delete( self->remeshDims );
   NewClass_Delete( self->staticWalls[0] );
   NewClass_Delete( self->staticWalls[1] );
   Class_Free( self, self->staticWalls );
   _NewClass_Destruct( self );
}

void _RegularRemesher_Print( const void* _self, Stream* stream ) {
   RegularRemesher* self = Class_Cast( _self, RegularRemesher );

   _NewRemesher_Print( self, stream );
}

void _RegularRemesher_Remesh( void* _self ) {
   RegularRemesher* self = Class_Cast( _self, RegularRemesher );
   Mesh* mesh;
   Grid *vGrid;
   Sync* meshSync;
   int nDims, nVerts;
   int center, ind, *inds;
   double leftCrd, rightCrd;
   int d_i, v_i, w_i;

   assert( self->mesh );

   mesh = self->mesh;
   meshSync = Mesh_GetSync( mesh, 0 );
   nDims = Mesh_GetDimSize( mesh );
   nVerts = Mesh_GetLocalSize( mesh, 0 );
   vGrid = *Mesh_GetExtension( mesh, Grid**, "vertexGrid" );
   inds = Class_Array( self, int, nDims );

   for( d_i = 0; d_i < nDims; d_i++ ) {
      for( w_i = 0; w_i < 2; w_i++ ) {
	 if( !ISet_Has( self->staticWalls[w_i], d_i ) )
	    continue;
	 for( v_i = 0; v_i < self->nWallVerts[d_i][w_i]; v_i++ ) {
	    mesh->verts[self->wallVerts[d_i][w_i][v_i]][d_i] = 
	       self->wallCrds[d_i][w_i][v_i];
	 }
      }
   }

   for( d_i = 0; d_i < nDims; d_i++ ) {
      if( !ISet_Has( self->remeshDims, d_i ) )
	 continue;

      Sync_SyncArray( self->syncs[d_i], 
		      mesh->verts[0] + d_i, nDims * sizeof(double), 
		      self->crds[d_i], sizeof(double), 
		      sizeof(double) );

      for( v_i = 0; v_i < nVerts; v_i++ ) {
	 ind = Sync_DomainToGlobal( meshSync, v_i );
	 Grid_Lift( vGrid, ind, (unsigned*)inds );
	 center = inds[d_i];
	 if( center == 0 || center == vGrid->sizes[d_i] - 1 )
	    continue;

	 inds[d_i] = 0;
	 ind = Grid_Project( vGrid, (unsigned*)inds );
	 ind = Sync_GlobalToDomain( self->syncs[d_i], ind );
	 if( ind >= nVerts )
	    leftCrd = self->crds[d_i][ind - nVerts];
	 else
	    leftCrd = Mesh_GetVertex( mesh, ind)[d_i];

	 inds[d_i] = vGrid->sizes[d_i] - 1;
	 ind = Grid_Project( vGrid, (unsigned*)inds );
	 ind = Sync_GlobalToDomain( self->syncs[d_i], ind );
	 if( ind >= nVerts )
	    rightCrd = self->crds[d_i][ind - nVerts];
	 else
	    rightCrd = Mesh_GetVertex( mesh, ind)[d_i];

	 mesh->verts[v_i][d_i] = leftCrd + 
	    (double)center * (rightCrd - leftCrd) / 
	    (double)(vGrid->sizes[d_i] - 1);
      }
   }

   Class_Free( self, inds );
   Mesh_Sync( mesh );
   Mesh_DeformationUpdate( mesh );
}

void RegularRemesher_Build( void* _self ) {
   RegularRemesher* self = Class_Cast( _self, RegularRemesher );
   Mesh* mesh;
   Sync* meshSync;
   const Decomp* meshDecomp;
   ISet* wallSet;
   Grid* vGrid;
   int *inds, ind;
   int nDims, nVerts;
   int nRems, *rems;
   int d_i, w_i, v_i;

   mesh = self->mesh;
   Stg_Component_Build( mesh, NULL, False );
   meshSync = Mesh_GetSync( mesh, 0 );
   meshDecomp = Sync_GetDecomp( meshSync );
   vGrid = *Mesh_GetExtension( mesh, Grid**, "vertexGrid" );
   nDims = Mesh_GetDimSize( mesh );
   self->nDims = nDims;
   nVerts = Mesh_GetLocalSize( mesh, 0 );
   inds = Class_Array( self, int, nDims );

   if( !self->syncs ) {
      self->syncs = Class_Array( self, Sync*, nDims );
      self->crds = Class_Array( self, double*, nDims );
      self->nWallVerts = Class_Array( self, int*, nDims );
      self->wallVerts = Class_Array( self, int**, nDims );
      self->wallCrds = Class_Array( self, double**, nDims );
      for( d_i = 0; d_i < nDims; d_i++ ) {
	 self->nWallVerts[d_i] = Class_Array( self, int, 2 );
	 self->wallVerts[d_i] = Class_Array( self, int*, 2 );
	 self->wallCrds[d_i] = Class_Array( self, double*, 2 );
	 memset( self->nWallVerts[d_i], 0, 2 * sizeof(int) );
	 memset( self->wallVerts[d_i], 0, 2 * sizeof(int*) );
	 memset( self->wallCrds[d_i], 0, 2 * sizeof(double*) );
      }
      memset( self->syncs, 0, nDims * sizeof(Sync*) );
      memset( self->crds, 0, nDims * sizeof(double*) );
   }

   wallSet = ISet_New();
   ISet_SetMaxSize( wallSet, nVerts );
   for( d_i = 0; d_i < nDims; d_i++ ) {
      for( w_i = 0; w_i < 2; w_i++ ) {
	 Class_Free( self, self->wallVerts[d_i][w_i] );
	 Class_Free( self, self->wallCrds[d_i][w_i] );
	 if( !ISet_Has( self->staticWalls[w_i], d_i ) ) {
	    self->nWallVerts[d_i][w_i] = 0;
	    self->wallVerts[d_i][w_i] = NULL;
	    self->wallCrds[d_i][w_i] = NULL;
	    continue;
	 }

	 for( v_i = 0; v_i < nVerts; v_i++ ) {
	    ind = Sync_DomainToGlobal( meshSync, v_i );
	    Grid_Lift( vGrid, ind, (unsigned*)inds );
	    if( (w_i == 0 && inds[d_i] == 0) || 
		(w_i == 1 && inds[d_i] == vGrid->sizes[d_i] - 1) )
	    {
	       ISet_Insert( wallSet, v_i );
	    }
	 }

	 self->nWallVerts[d_i][w_i] = ISet_GetSize( wallSet );
	 self->wallVerts[d_i][w_i] = Class_Array( 
	    self, int, self->nWallVerts[d_i][w_i] );
	 self->wallCrds[d_i][w_i] = Class_Array( 
	    self, double, self->nWallVerts[d_i][w_i] );
	 ISet_GetArray( wallSet, self->wallVerts[d_i][w_i] );
	 ISet_Clear( wallSet );
	 for( v_i = 0; v_i < self->nWallVerts[d_i][w_i]; v_i++ ) {
	    self->wallCrds[d_i][w_i][v_i] = 
	       mesh->verts[self->wallVerts[d_i][w_i][v_i]][d_i];
	 }
      }
   }

   for( d_i = 0; d_i < nDims; d_i++ ) {
      NewClass_RemoveRef( self->syncs[d_i] );
      Class_Free( self, self->crds[d_i] );
      if( !ISet_Has( self->remeshDims, d_i ) ) {
	 self->syncs[d_i] = NULL;
	 self->crds[d_i] = NULL;
	 continue;
      }

      for( v_i = 0; v_i < nVerts; v_i++ ) {
	 ind = Sync_DomainToGlobal( meshSync, v_i );
	 Grid_Lift( vGrid, ind, (unsigned*)inds );
	 inds[d_i] = 0;
	 ind = Grid_Project( vGrid, (unsigned*)inds );
	 if( !Sync_TryGlobalToDomain( meshSync, ind, &ind ) )
	    ISet_TryInsert( wallSet, ind );

	 inds[d_i] = vGrid->sizes[d_i] - 1;
	 ind = Grid_Project( vGrid, (unsigned*)inds );
	 if( !Sync_TryGlobalToDomain( meshSync, ind, &ind ) )
	    ISet_TryInsert( wallSet, ind );
      }

      nRems = ISet_GetSize( wallSet );
      rems = Class_Array( self, int, nRems );
      ISet_GetArray( wallSet, rems );
      ISet_Clear( wallSet );
      self->syncs[d_i] = Sync_New();
      NewClass_AddRef( self->syncs[d_i] );
      Sync_SetDecomp( self->syncs[d_i], meshDecomp );
      Sync_FindRemotes( self->syncs[d_i], nRems, rems );
      Class_Free( self, rems );
      self->crds[d_i] = Class_Array( self, double, nRems );
   }

   Class_Free( self, inds );
   NewClass_Delete( wallSet );
}

void RegularRemesher_SetRemeshState( void* _self, int dim, Bool state ) {
   RegularRemesher* self = Class_Cast( _self, RegularRemesher );

   if( state )
      ISet_TryInsert( self->remeshDims, dim );
   else if( ISet_Has( self->remeshDims, dim ) )
      ISet_Remove( self->remeshDims, dim );
}

void RegularRemesher_SetStaticWall( void* _self, int dim, int wall, Bool state ) {
   RegularRemesher* self = Class_Cast( _self, RegularRemesher );

   if( state )
      ISet_TryInsert( self->staticWalls[wall], dim );
   else if( ISet_Has( self->staticWalls[wall], dim ) )
      ISet_Remove( self->staticWalls[wall], dim );
}
