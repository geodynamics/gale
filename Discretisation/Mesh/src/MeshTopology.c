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
** $Id: MeshTopology.c 3952 2007-01-09 06:24:06Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdlib.h>
#include <string.h>
#include "StGermain/Base/Base.h"
#include "types.h"
#include "Decomp.h"
#include "Sync.h"
#include "MeshTopology.h"
#include "StGermain/Base/Foundation/ClassDef.h"


void _MeshTopology_Construct( void* _self ) {
   MeshTopology* self = (MeshTopology*)_self;

   _NewClass_Construct( self );
   self->nDims = 0;
   self->nTDims = 0;
   self->shadDepth = 0;
   self->comm = NULL;
   self->locals = NULL;
   self->remotes = NULL;
   self->nGhosts = 0;
   self->ghosts = NULL;
   self->nIncEls = NULL;
   self->incEls = NULL;
}

void _MeshTopology_Destruct( void* _self ) {
   MeshTopology* self = (MeshTopology*)_self;

   MeshTopology_Clear( self );
   _NewClass_Destruct( self );
}

void _MeshTopology_Copy( void* _self, const void* op ) {
   /*MeshTopology* self = (MeshTopology*)_self;*/

   assert( 0 );
   /* TODO */
}

SizeT _MeshTopology_CalcMem( const void* _self, PtrMap* ptrs ) {
   MeshTopology* self = (MeshTopology*)_self;
   SizeT mem;
   int d_i;

   if( PtrMap_Find( ptrs, (void*)self ) )
      return 0;
   mem = _NewClass_CalcMem( self, ptrs );
   for( d_i = 0; d_i < self->nTDims; d_i++ ) {
      if( self->locals )
	 mem += NewClass_CalcMem( self->locals[d_i], ptrs );
      if( self->remotes )
	 mem += NewClass_CalcMem( self->remotes[d_i], ptrs );
   }
   return mem;
}

void MeshTopology_SetNumDims( void* _self, int nDims ) {
   MeshTopology* self = (MeshTopology*)_self;
   int d_i;

   MeshTopology_ClearDims( self );
   self->nDims = nDims;
   if( self->nDims > 0 )
      self->nTDims = nDims + 1;
   self->locals = Class_Array( self, Decomp*, self->nTDims );
   self->remotes = Class_Array( self, Sync*, self->nTDims );
   self->nIncEls = Class_Array2D( self, int*, self->nTDims, self->nTDims );
   self->incEls = Class_Array2D( self, int**, self->nTDims, self->nTDims );
   for( d_i = 0; d_i < self->nTDims; d_i++ ) {
      self->locals[d_i] = Decomp_New();
      NewClass_AddRef( self->locals[d_i] );
      self->remotes[d_i] = Sync_New();
      NewClass_AddRef( self->remotes[d_i] );
      Decomp_SetComm( self->locals[d_i], self->comm );
      Sync_SetDecomp( self->remotes[d_i], self->locals[d_i] );
      memset( self->nIncEls[d_i], 0, sizeof(int**) * self->nTDims );
      memset( self->incEls[d_i], 0, sizeof(int***) * self->nTDims );
   }
}

void MeshTopology_SetComm( void* _self, const Comm* comm ) {
   MeshTopology* self = (MeshTopology*)_self;
   const Comm* curComm;
   int d_i;

   assert( self );
   self->comm = (Comm*)comm;
   if( comm ) {
      NewClass_AddRef( (Comm*)comm );
      for( d_i = 0; d_i < self->nTDims; d_i++ ) {
	 curComm = Decomp_GetComm( self->locals[d_i] );
	 if( curComm != comm )
	    Decomp_SetComm( self->locals[d_i], self->comm );
      }
   }
}

void MeshTopology_SetDomain( void* _self, int dim, Sync* sync ) {
   MeshTopology* self = (MeshTopology*)_self;

   assert( self && dim < self->nTDims );
   NewClass_RemoveRef( self->locals[dim] );
   NewClass_RemoveRef( self->remotes[dim] );
   self->remotes[dim] = sync;
   NewClass_AddRef( sync );
   if( sync ) {
      self->locals[dim] = (Decomp*)Sync_GetDecomp( sync );
      NewClass_AddRef( self->locals[dim] );
   }
}

void MeshTopology_SetLocalElements( void* _self, int dim, int nEls, 
				    const int* globals )
{
   MeshTopology* self = (MeshTopology*)_self;
   int d_i, e_i;

   assert( self );
   assert( dim < self->nTDims );
   assert( !nEls || globals );
   for( d_i = 0; d_i < self->nTDims; d_i++ ) {
      for( e_i = 0; e_i < Sync_GetNumDomains( self->remotes[dim] ); d_i++ )
	 Class_Free( self, self->incEls[dim][d_i][e_i] );
      Class_Free( self, self->incEls[dim][d_i] );
      Class_Free( self, self->nIncEls[dim][d_i] );
      self->incEls[dim][d_i] = NULL;
      self->nIncEls[dim][d_i] = NULL;
   }
   Decomp_SetLocals( self->locals[dim], nEls, globals );
   Sync_SetDecomp( self->remotes[dim], self->locals[dim] );
}

void MeshTopology_AddLocalElements( void* _self, int dim, int nEls, 
				    const int* globals )
{
   MeshTopology* self = (MeshTopology*)_self;

   assert( 0 );
   /* TODO */

   assert( self );
   assert( dim < self->nTDims );
   assert( !nEls || globals );
   Decomp_AddLocals( self->locals[dim], nEls, globals );
   /* TODO: Expand this dimensions incidence if some incidence already set. */
}

void MeshTopology_RemoveLocalElements( void* _self, int dim, int nEls, 
				       const int* globals, IMap* map )
{
   /*MeshTopology* self = (MeshTopology*)_self;*/

   assert( 0 );
   /* TODO: Method body goes here */
}

void MeshTopology_SetRemoteElements( void* _self, int dim, int nEls, 
				     const int* globals )
{
   MeshTopology* self = (MeshTopology*)_self;
   int nDoms;
   int d_i, e_i;

   assert( self );
   assert( dim < self->nTDims );
   assert( !nEls || globals );
   for( d_i = 0; d_i < self->nTDims; d_i++ ) {
      if( self->nIncEls[dim][d_i] ) {
	 for( e_i = Decomp_GetNumLocals( self->locals[dim] );
	      e_i < Sync_GetNumDomains( self->remotes[dim] );
	      e_i++ )
	 {
	    Class_Free( self, self->incEls[dim][d_i][e_i] );
	 }
      }
   }

   Sync_SetRemotes( self->remotes[dim], nEls, globals );

   for( d_i = 0; d_i < self->nTDims; d_i++ ) {
      if( self->nIncEls[dim][d_i] ) {
	 nDoms = Sync_GetNumDomains( self->remotes[dim] );
	 self->nIncEls[dim][d_i] = Class_Rearray( self, self->incEls[dim][d_i], 
						  int, nDoms );
	 self->incEls[dim][d_i] = Class_Rearray( self, self->incEls[dim][d_i], 
						 int*, nDoms );
	 for( e_i = Decomp_GetNumLocals( self->locals[dim] );
	      e_i < nDoms; 
	      e_i++ )
	 {
	    self->nIncEls[dim][d_i][e_i] = 0;
	    self->incEls[dim][d_i][e_i] = NULL;
	 }
      }
   }
}

void MeshTopology_AddRemoteElements( void* _self, int dim, int nEls, 
				     const int* globals )
{
   MeshTopology* self = (MeshTopology*)_self;
   int nOldDoms, nDoms;
   int d_i, e_i;

   assert( self );
   assert( dim < self->nTDims );
   assert( !nEls || globals );
   nOldDoms = Sync_GetNumDomains( self->remotes[dim] );
   Sync_AddRemotes( self->remotes[dim], nEls, globals );
   for( d_i = 0; d_i < self->nTDims; d_i++ ) {
      if( self->nIncEls[dim][d_i] ) {
	 nDoms = Sync_GetNumDomains( self->remotes[dim] );;
	 self->nIncEls[dim][d_i] = Class_Rearray( self, self->incEls[dim][d_i], 
						  int, nDoms );
	 self->incEls[dim][d_i] = Class_Rearray( self, self->incEls[dim][d_i], 
						 int*, nDoms );
	 for( e_i = nOldDoms; e_i < nDoms; e_i++ ) {
	    self->nIncEls[dim][d_i][e_i] = 0;
	    self->incEls[dim][d_i][e_i] = NULL;
	 }
      }
   }
}

void MeshTopology_RemoveRemoteElements( void* _self, int dim, int nEls, 
					const int* globals, IMap* map )
{
   /*MeshTopology* self = (MeshTopology*)_self;*/

   assert( 0 );
   /* TODO: Method body goes here */
}

void MeshTopology_SetIncidence( void* _self, int fromDim, int fromEl, 
				int toDim, int nIncEls, const int* incEls  )
{
   MeshTopology* self = (MeshTopology*)_self;
   int nDoms;

   assert( self );
   assert( fromDim < self->nTDims && toDim < self->nTDims );
   assert( self->locals[fromDim] );
   if( !self->nIncEls[fromDim][toDim] ) {
      nDoms = Sync_GetNumDomains( self->remotes[fromDim] );
      self->nIncEls[fromDim][toDim] = Class_Array( self, int, nDoms );
      self->incEls[fromDim][toDim] = Class_Array( self, int*, nDoms );
   }

   self->nIncEls[fromDim][toDim][fromEl] = nIncEls;
   self->incEls[fromDim][toDim][fromEl] = Class_Array( self, int, nIncEls );
   memcpy( self->incEls[fromDim][toDim][fromEl], incEls, nIncEls * sizeof(int) );
}

void MeshTopology_InvertIncidence( void* _self, int fromDim, int toDim ) {
   MeshTopology* self = (MeshTopology*)_self;
   int fromSize, toSize;
   int *nInvIncEls, **invIncEls;
   int *nIncEls, **incEls;
   int elInd;
   int e_i, inc_i;

   assert( self );
   fromSize = Sync_GetNumDomains( self->remotes[fromDim] );
   toSize = Sync_GetNumDomains( self->remotes[toDim] );
   nInvIncEls = self->nIncEls[toDim][fromDim];
   invIncEls = self->incEls[toDim][fromDim];
   nIncEls = Class_Array( self, int, fromSize );
   incEls = Class_Array( self, int*, fromSize );
   memset( nIncEls, 0, fromSize * sizeof(int) );
   for( e_i = 0; e_i < toSize; e_i++ ) {
      for( inc_i = 0; inc_i < nInvIncEls[e_i]; inc_i++ )
	 nIncEls[invIncEls[e_i][inc_i]]++;
   }

   for( e_i = 0; e_i < fromSize; e_i++ )
      incEls[e_i] = Class_Array( self, int, nIncEls[e_i] );
   memset( nIncEls, 0, fromSize * sizeof(unsigned) );
   for( e_i = 0; e_i < toSize; e_i++ ) {
      for( inc_i = 0; inc_i < nInvIncEls[e_i]; inc_i++ )
	 elInd = invIncEls[e_i][inc_i];
	 incEls[elInd][nIncEls[elInd]++] = e_i;
   }

   for( e_i = 0; e_i < fromSize; e_i++ ) {
      MeshTopology_SetIncidence( self, fromDim, e_i, toDim, nIncEls[e_i], incEls[e_i] );
      Class_Free( self, incEls[e_i] );
   }
   FreeArray( nIncEls );
   FreeArray( incEls );
}

void MeshTopology_SetGhostVertices( void* _self, int nGhosts, const int* globals ) {
   MeshTopology* self = (MeshTopology*)_self;

   assert( self && (!nGhosts || globals) );
   self->ghosts = Class_Rearray( self, self->ghosts, int, nGhosts );
   memcpy( self->ghosts, globals, nGhosts * sizeof(int) );
   self->nGhosts = nGhosts;
}

void MeshTopology_SetShadowDepth( void* _self, int depth ) {
/*
  MeshTopology* self = (MeshTopology*)_self;
  int nLocals;
  const int* locals;

  assert( self && depth >= 0 );
  assert( self->comm );

  nNbrs = Comm_GetNumNeighbours( self->comm );
  nNbrGhosts = Class_Array( self, int, nNbrs );
  nbrGhosts = Class_Array( self, int*, nNbrs );
  Comm_AllgatherInit( self->comm, self->nGhosts, nNbrGhosts, sizeof(int) );
  for( n_i = 0; n_i < nNbrs; n_i++ )
  nbrGhosts[n_i] = Class_Array( self, int, nNbrGhosts[n_i] );
  Comm_AllgatherBegin( self->comm, self->ghosts, nbrGhosts );
  Comm_AllgatherEnd( self->comm );

  ISet_Init( locSet );
  ISet_Init( ghostSet );
  ISet_Init( isect );
  Decomp_GetLocals( self->locals[0], &nLocals, &locals );
  ISet_UseArray( locSet, nLocals, locals );
  for( n_i = 0; n_i < nNbrs; n_i++ ) {
  ISet_UseArray( ghostSet, nNbrGhosts[n_i], nbrGhosts[n_i] );
  Class_Free( self, nbrGhosts[n_i] );
  ISet_Isect( locSet, ghostSet, isect );
  ISet_Clear( ghostSet );
  nNbrGhosts[n_i] = ISet_GetSize( isect );
  nbrGhosts[n_i] = Class_Array( self, int, nNbrGhosts[n_i] );
  ISet_GetArray( isect, nNbrGhosts + n_i, nbrGhosts[n_i] );
  ISet_Clear( isect );
  }
  ISet_Destruct( locSet );
  ISet_Destruct( ghostSet );

  for( n_i = 0; n_i < nNbrs; n_i++ ) {
  for( v_i = 0; v_i < nNbrGhosts[n_i]; v_i++ ) {
  insist( Decomp_GlobalToLocal( self->locals[0], nbrGhosts[n_i][v_i], &loc ) );
  for( inc_i = 0; inc_i < self->nIncEls[0][nDims][loc]; inc_i++ )
  ISet_TryInsert( isect, self->incEls[0][nDims][loc][inc_i] );
  }
  nNbrGhosts[n_i] = ISet_GetSize( isect );
  nbrGhosts[n_i] = Class_Rearray( self, nbrGhosts[n_i], int, nNbrGhosts[n_i] );
  ISet_GetArray( isect, nNbrGhosts + n_i, nbrGhosts[n_i] );
  for( s_i = 0; s_i < depth - 1; s_i++ )
  MeshTopology_ExpandShadows( self, nNbrGhosts + n_i, nbrGhosts + n_i, isect );
  ISet_Clear( isect );
  }

  Comm_AlltoallInit( self->comm, nNbr
*/
}

void MeshTopology_Clear( void* self ) {
   MeshTopology_ClearDims( self );
   if( ((MeshTopology*)self)->comm )
      NewClass_RemoveRef( ((MeshTopology*)self)->comm );
   ((MeshTopology*)self)->comm = NULL;
}

void MeshTopology_ClearDims( void* _self ) {
   MeshTopology* self = (MeshTopology*)_self;
   int d_i;

   MeshTopology_ClearElements( self );
   for( d_i = 0; d_i < self->nTDims; d_i++ ) {
      if( self->locals && self->locals[d_i] )
	 NewClass_RemoveRef( self->locals[d_i] );
      if( self->remotes && self->remotes[d_i] )
	 NewClass_RemoveRef( self->remotes[d_i] );
   }
   Class_Free( self, self->locals );
   Class_Free( self, self->remotes );
   Class_Free( self, self->nIncEls );
   Class_Free( self, self->incEls );

   self->nDims = 0;
   self->nTDims = 0;
   self->shadDepth = 0;
   self->locals = NULL;
   self->remotes = NULL;
   self->nIncEls = NULL;
   self->incEls = NULL;
}

void MeshTopology_ClearElements( void* _self ) {
   MeshTopology* self = (MeshTopology*)_self;
   int d_i;

   MeshTopology_ClearIncidence( self );
   for( d_i = 0; d_i < self->nTDims; d_i++ ) {
      if( self->locals )
	 Decomp_ClearLocals( self->locals[d_i] );
      if( self->remotes && self->remotes[d_i] )
	 Sync_ClearRemotes( self->remotes[d_i] );
   }
}

void MeshTopology_ClearIncidence( void* _self ) {
   MeshTopology* self = (MeshTopology*)_self;
   int d_i, d_j, e_i;

   assert( self );
   for( d_i = 0; d_i < self->nTDims; d_i++ ) {
      for( d_j = 0; d_j < self->nTDims; d_j++ ) {
	 if( self->nIncEls[d_i][d_j] ) {
	    for( e_i = 0; e_i < Decomp_GetNumLocals( self->locals[d_i] ); e_i++ )
	       Class_Free( self, self->incEls[d_i][d_j][e_i] );
	    Class_Free( self, self->incEls[d_i][d_j] );
	    Class_Free( self, self->nIncEls[d_i][d_j] );
	    self->nIncEls[d_i][d_j] = NULL;
	    self->incEls[d_i][d_j] = NULL;
	 }
      }
   }
}

int MeshTopology_GetNumDims( const void* self ) {
   assert( self );
   return ((MeshTopology*)self)->nDims;
}

const Comm* MeshTopology_GetComm( const void* self ) {
   assert( self );
   return ((MeshTopology*)self)->comm;
}

Bool MeshTopology_HasDomain( const void* self, int dim ) {
   assert( self && dim < ((MeshTopology*)self)->nTDims );
   return Sync_GetNumDomains( ((MeshTopology*)self)->remotes[dim] ) ? 
      True : False;
}

const Sync* MeshTopology_GetDomain( const void* self, int dim ) {
   assert( self && dim < ((MeshTopology*)self)->nTDims );
   return ((MeshTopology*)self)->remotes[dim];
}

Bool MeshTopology_HasIncidence( const void* self, int fromDim, int toDim ) {
   assert( self );
   assert( fromDim < ((MeshTopology*)self)->nTDims );
   assert( toDim < ((MeshTopology*)self)->nTDims );
   return ((MeshTopology*)self)->nIncEls[fromDim][toDim] ? True : False;
}

int MeshTopology_GetIncidenceSize( const void* self, int fromDim, int fromEl, 
				   int toDim )
{
   assert( self );
   assert( fromDim < ((MeshTopology*)self)->nTDims );
   assert( toDim < ((MeshTopology*)self)->nTDims );
   assert( fromEl < Sync_GetNumDomains( ((MeshTopology*)self)->remotes[fromDim] ) );
   return ((MeshTopology*)self)->nIncEls[fromDim][toDim][fromEl];
}

void MeshTopology_GetIncidence( const void* self, int fromDim, int fromEl, 
				int toDim, int* nIncEls, const int** incEls )
{
   assert( self );
   assert( fromDim < ((MeshTopology*)self)->nTDims );
   assert( toDim < ((MeshTopology*)self)->nTDims );
   assert( fromEl < Sync_GetNumDomains( ((MeshTopology*)self)->remotes[fromDim] ) );
   *nIncEls = ((MeshTopology*)self)->nIncEls[fromDim][toDim][fromEl];
   *incEls = ((MeshTopology*)self)->incEls[fromDim][toDim][fromEl];
}
