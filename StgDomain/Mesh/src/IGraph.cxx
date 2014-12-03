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
** $Id: IGraph.c 3952 2007-01-09 06:24:06Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdlib.h>
#include <string.h>
#include "StGermain/StGermain.h"
#include "types.h"
#include "Decomp.h"
#include "Sync.h"
#include "MeshTopology.h"
#include "IGraph.h"
#include "StGermain/Base/Foundation/ClassDef.h"


void IGraph_PickleIncidenceInit( IGraph* self, int dim, int nEls, int* els, int* nBytes );
void IGraph_PickleIncidence( IGraph* self, int dim, int nEls, int* els, stgByte* bytes );
void IGraph_UnpickleIncidence( IGraph* self, int dim, int nBytes, stgByte* bytes );
int IGraph_Cmp( const void* l, const void* r );


void _IGraph_Init( void* _self ) {
   IGraph* self = (IGraph*)_self;

   _MeshTopology_Init( self );
   self->nDims = (MeshTopology_Dim)0;
   self->nTDims = 0;
   self->shadDepth = 0;
   self->comm = NULL;
   self->locals = NULL;
   self->remotes = NULL;
   self->nBndEls = NULL;
   self->bndEls = NULL;
   self->nIncEls = NULL;
   self->incEls = NULL;
}

void _IGraph_Destruct( void* _self ) {
   IGraph* self = (IGraph*)_self;

   IGraph_Clear( self );
   _MeshTopology_Destruct( self );
}

void _IGraph_Copy( void* _self, const void* op ) {
   /*IGraph* self = (IGraph*)_self;*/

   assert( 0 );
   /* TODO */
}

SizeT _IGraph_CalcMem( const void* _self, PtrMap* ptrs ) {
   IGraph* self = (IGraph*)_self;
   SizeT mem;
   int d_i;

   if( PtrMap_Find( ptrs, (void*)self ) )
      return 0;
   mem = _MeshTopology_CalcMem( self, ptrs );

	for( d_i = 0; d_i < self->nTDims; d_i++ ) {
      if( self->locals )
			mem += NewClass_CalcMem( self->locals[d_i], ptrs );
		if( self->remotes )
			mem += NewClass_CalcMem( self->remotes[d_i], ptrs );
   }
   return mem;
}

void _IGraph_SetNumDims( void* _self, int nDims ) {
   IGraph* self = (IGraph*)_self;
   int d_i;

   IGraph_ClearDims( self );
   _MeshTopology_SetNumDims( self, nDims );

   self->locals = Class_Array( self, Decomp*, self->nTDims );
   self->remotes = Class_Array( self, Sync*, self->nTDims );
   self->nIncEls = Class_Array2D( self, int*, self->nTDims, self->nTDims );
   self->incEls = Class_Array2D( self, int**, self->nTDims, self->nTDims );

   for( d_i = 0; d_i < self->nTDims; d_i++ ) {
		self->locals[d_i] = Decomp_New();
      NewClass_AddRef( self->locals[d_i] );
      self->remotes[d_i] = Sync_New();
		NewClass_AddRef( self->remotes[d_i] );

		if( self->comm )
			Decomp_SetMPIComm( self->locals[d_i], Comm_GetMPIComm( self->comm ) );
		Sync_SetDecomp( self->remotes[d_i], self->locals[d_i] );

      if( self->comm )
			Sync_SetComm( self->remotes[d_i], self->comm );

		memset( self->nIncEls[d_i], 0, sizeof(int**) * self->nTDims );
		memset( self->incEls[d_i], 0, sizeof(int***) * self->nTDims );
   }
}

void IGraph_SetComm( void* _self, const Comm* comm ) {
   IGraph* self = (IGraph*)_self;
   const Comm* curComm;
   int d_i;

   assert( self );

   IGraph_ClearElements( self );
   _MeshTopology_SetComm( self, comm );

   if( comm ) {
      for( d_i = 0; d_i < self->nTDims; d_i++ ) {
			curComm = Sync_GetComm( self->remotes[d_i] );

			if( curComm != comm ) {
				Sync_SetComm( self->remotes[d_i], self->comm );
				Decomp_SetMPIComm( self->locals[d_i], Comm_GetMPIComm( self->comm ) );
			}
      }
   }
}

void IGraph_SetDomain( void* _self, int dim, Sync* sync ) {
   IGraph* self = (IGraph*)_self;

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

void IGraph_SetElements( void* _self, int dim, int nEls, const int* globals ) {
   IGraph* self = (IGraph*)_self;
   int rank;
   int nNbrs;
   const int *nbrs;
   int nSubEls;
   const int *subEls;
   int rem, netRem;
   int *nNbrEls, **nbrEls;
   IArray** isects;
   ISet localsObj, *locals = &localsObj;
   ISet remotesObj, *remotes = &remotesObj;
   int nCurEls, *curEls;
   MPI_Comm mpiComm;
   int n_i, e_i;

   assert( self && dim < self->nTDims );
   assert( !nEls || globals );
   assert( self->comm );

   Comm_GetNeighbours( self->comm, &nNbrs, &nbrs );
   if( !nNbrs ) {
      IGraph_SetLocalElements( self, dim, nEls, globals );
      return;
   }

   ISet_Construct( locals );
   mpiComm = Comm_GetMPIComm( self->comm );
   insist( MPI_Comm_rank( mpiComm, &rank ), == MPI_SUCCESS );
   isects = Class_Array( self, IArray*, nNbrs );
   for( n_i = 0; n_i < nNbrs; n_i++ )
      isects[n_i] = IArray_New();

   ISet_UseArray( locals, nEls, globals );
   nSubEls = (nEls < 1000) ? nEls : 1000;
   rem = nEls;
   subEls = globals;
   nNbrEls = Class_Array( self, int, nNbrs );
   nbrEls = Class_Array( self, int*, nNbrs );

	do {
		Comm_AllgatherInit( self->comm, nSubEls, nNbrEls, sizeof(int) );

      for( n_i = 0; n_i < nNbrs; n_i++ )
			nbrEls[n_i] = Class_Array( self, int, nNbrEls[n_i] );

		Comm_AllgatherBegin( self->comm, subEls, (void**)nbrEls );
      Comm_AllgatherEnd( self->comm );

		for( n_i = 0; n_i < nNbrs; n_i++ ) {
			for( e_i = 0; e_i < nNbrEls[n_i]; e_i++ ) {
				if( ISet_Has( locals, nbrEls[n_i][e_i] ) )
					IArray_Append( isects[n_i], nbrEls[n_i][e_i] );
			}
			Class_Free( self, nbrEls[n_i] );
		}

      subEls += nSubEls;
      rem -= nSubEls;
      nSubEls = (rem < 1000) ? rem : 1000;
      insist( MPI_Allreduce( &rem, &netRem, 1, MPI_INT, MPI_SUM, mpiComm ), == MPI_SUCCESS );
   } while( netRem );
   Class_Free( self, nNbrEls );
   Class_Free( self, nbrEls );

   ISet_Construct( remotes );
   ISet_SetMaxSize( remotes, nEls );
   for( n_i = 0; n_i < nNbrs; n_i++ ) {
      IArray_GetArray( isects[n_i], &nCurEls, (const int**)&curEls );
      if( nbrs[n_i] < rank ) {
			for( e_i = 0; e_i < nCurEls; e_i++ ) {
				ISet_TryRemove( locals, curEls[e_i] );
				ISet_TryInsert( remotes, curEls[e_i] );
			}
      }
      NewClass_Delete( isects[n_i] );
   }
   Class_Free( self, isects );

   nCurEls = ISet_GetSize( locals );
   curEls = Class_Array( self, int, nCurEls );
   ISet_GetArray( locals, curEls );
   ISet_Destruct( locals );
   qsort( curEls, nCurEls, sizeof(int), IGraph_Cmp );
   IGraph_SetLocalElements( self, dim, nCurEls, curEls );
   Class_Free( self, curEls );
   nCurEls = ISet_GetSize( remotes );
   curEls = Class_Array( self, int, nCurEls );
   ISet_GetArray( remotes, curEls );
   ISet_Destruct( remotes );
   qsort( curEls, nCurEls, sizeof(int), IGraph_Cmp );
   IGraph_SetRemoteElements( self, dim, nCurEls, curEls );
   Class_Free( self, curEls );
}

void IGraph_SetLocalElements( void* _self, int dim, int nEls, const int* globals ) {
   IGraph* self = (IGraph*)_self;
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

void IGraph_AddLocalElements( void* _self, int dim, int nEls, const int* globals ) {
   IGraph* self = (IGraph*)_self;

   assert( 0 );
   /* TODO */

   assert( self );
   assert( dim < self->nTDims );
   assert( !nEls || globals );
   Decomp_AddLocals( self->locals[dim], nEls, globals );
   /* TODO: Expand this dimensions incidence if some incidence already set. */
}

void IGraph_RemoveLocalElements( void* _self, int dim, int nEls, const int* globals, IMap* map ) {
   /*IGraph* self = (IGraph*)_self;*/

   assert( 0 );
   /* TODO: Method body goes here */
}

void IGraph_SetRemoteElements( void* _self, int dim, int nEls, const int* globals ) {
   IGraph* self = (IGraph*)_self;
   int nDoms;
   int d_i, e_i;

   assert( self );
   assert( dim < self->nTDims );
   assert( !nEls || globals );
   for( d_i = 0; d_i < self->nTDims; d_i++ ) {
		if( self->nIncEls[dim][d_i] ) {
			for( e_i = Decomp_GetNumLocals( self->locals[dim] ); e_i < Sync_GetNumDomains( self->remotes[dim] ); e_i++ ) {
				Class_Free( self, self->incEls[dim][d_i][e_i] );
			}
      }
   }

   Sync_SetRemotes( self->remotes[dim], nEls, globals );

   for( d_i = 0; d_i < self->nTDims; d_i++ ) {
      if( self->nIncEls[dim][d_i] ) {
			nDoms = Sync_GetNumDomains( self->remotes[dim] );
			self->nIncEls[dim][d_i] = Class_Rearray( self, self->incEls[dim][d_i], int, nDoms );
			self->incEls[dim][d_i] = Class_Rearray( self, self->incEls[dim][d_i], int*, nDoms );

			for( e_i = Decomp_GetNumLocals( self->locals[dim] ); e_i < nDoms; e_i++ ) {
				self->nIncEls[dim][d_i][e_i] = 0;
				self->incEls[dim][d_i][e_i] = NULL;
			}
      }
   }
}

void IGraph_AddRemoteElements( void* _self, int dim, int nEls, const int* globals ) {
   IGraph* self = (IGraph*)_self;
   int nOldDoms, nDoms;
   int d_i, e_i;

   assert( self );
   assert( dim < self->nTDims );
   assert( !nEls || globals );
   nOldDoms = Sync_GetNumDomains( self->remotes[dim] );
   Sync_AddRemotes( self->remotes[dim], nEls, globals );
   for( d_i = 0; d_i < self->nTDims; d_i++ ) {
      if( self->nIncEls[dim][d_i] ) {
			nDoms = Sync_GetNumDomains( self->remotes[dim] );
			self->nIncEls[dim][d_i] = Class_Rearray( self, self->nIncEls[dim][d_i], int, nDoms );
			self->incEls[dim][d_i] = Class_Rearray( self, self->incEls[dim][d_i], int*, nDoms );

			for( e_i = nOldDoms; e_i < nDoms; e_i++ ) {
				self->nIncEls[dim][d_i][e_i] = 0;
				self->incEls[dim][d_i][e_i] = NULL;
			}
      }
   }
}

void IGraph_RemoveRemoteElements( void* _self, int dim, int nEls, const int* globals, IMap* map ) {
   /*IGraph* self = (IGraph*)_self;*/

   assert( 0 );
   /* TODO: Method body goes here */
}

void IGraph_SetBoundaryElements( void* _self, int dim, int nEls, const int* els ) {
   IGraph* self = Class_Cast( _self, IGraph );

   assert( !nEls || els );
   assert( dim < self->nTDims );

   if( !self->nBndEls ) {
      self->nBndEls = Class_Array( self, int, self->nTDims );
      memset( self->nBndEls, 0, sizeof(int) * self->nTDims );
   }
   if( !self->bndEls ) {
      self->bndEls = Class_Array( self, int*, self->nTDims );
      memset( self->bndEls, 0, sizeof(int*) * self->nTDims );
   }
   self->nBndEls[dim] = nEls;
   self->bndEls[dim] = Class_Rearray( self, self->bndEls[dim], int, nEls );
   memcpy( self->bndEls[dim], els, sizeof(int) * nEls );
}

void IGraph_SetIncidence( void* _self, int fromDim, int fromEl, int toDim, int nIncEls, const int* incEls  ) {
   IGraph* self = (IGraph*)_self;
   int nDoms;

   assert( self );
   assert( fromDim < self->nTDims && toDim < self->nTDims );
   assert( self->locals[fromDim] );
   if( !self->nIncEls[fromDim][toDim] ) {
      nDoms = Sync_GetNumDomains( self->remotes[fromDim] );
      self->nIncEls[fromDim][toDim] = Class_Array( self, int, nDoms );
      self->incEls[fromDim][toDim] = Class_Array( self, int*, nDoms );
      memset( self->incEls[fromDim][toDim], 0, sizeof(int*) * nDoms );
   }

   self->nIncEls[fromDim][toDim][fromEl] = nIncEls;
   self->incEls[fromDim][toDim][fromEl] = Class_Rearray( self, self->incEls[fromDim][toDim][fromEl], int, nIncEls );
   memcpy( self->incEls[fromDim][toDim][fromEl], incEls, nIncEls * sizeof(int) );
}

void IGraph_RemoveIncidence( void* _self, int fromDim, int toDim ) {
   IGraph* self = (IGraph*)_self;
   int nEls;
   int e_i;

   assert( self );
   assert( fromDim < self->nTDims );
   assert( toDim < self->nTDims );

   nEls = Sync_GetNumDomains( self->remotes[fromDim] );
   for( e_i = 0; e_i < nEls; e_i++ )
      Class_Free( self, self->incEls[fromDim][toDim][e_i] );
   Class_Free( self, self->incEls[fromDim][toDim] );
   Class_Free( self, self->nIncEls[fromDim][toDim] );
   self->incEls[fromDim][toDim] = NULL;
   self->nIncEls[fromDim][toDim] = NULL;
}

void IGraph_InvertIncidence( void* _self, int fromDim, int toDim ) {
   IGraph* self = (IGraph*)_self;
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
   memset( nIncEls, 0, fromSize * sizeof(int) );
   for( e_i = 0; e_i < toSize; e_i++ ) {
      for( inc_i = 0; inc_i < nInvIncEls[e_i]; inc_i++ )
			nIncEls[invIncEls[e_i][inc_i]]++;
   }

   incEls = Class_Array( self, int*, fromSize );
   for( e_i = 0; e_i < fromSize; e_i++ )
      incEls[e_i] = Class_Array( self, int, nIncEls[e_i] );
   memset( nIncEls, 0, fromSize * sizeof(unsigned) );
   for( e_i = 0; e_i < toSize; e_i++ ) {
		for( inc_i = 0; inc_i < nInvIncEls[e_i]; inc_i++ ) {
			elInd = invIncEls[e_i][inc_i];
			incEls[elInd][nIncEls[elInd]++] = e_i;
      }
   }

   for( e_i = 0; e_i < fromSize; e_i++ ) {
      IGraph_SetIncidence( self, fromDim, e_i, toDim, nIncEls[e_i], incEls[e_i] );
      Class_Free( self, incEls[e_i] );
   }
   Class_Free( self, nIncEls );
   Class_Free( self, incEls );
}

void IGraph_ExpandIncidence( void* _self, int dim ) {
   IGraph* self = (IGraph*)_self;
   ISet nbrSetObj, *nbrSet = &nbrSetObj;
   int nEls;
   int nCurNbrs, maxNbrs;
   int nIncEls, *incEls;
   int nUpEls, *upEls;
   int e_i, inc_i, inc_j;

   assert( self );
   assert( dim < self->nTDims );
   assert( dim > 0 );
   assert( self->nIncEls[dim][0] );

   nEls = Sync_GetNumDomains( self->remotes[dim] );
   maxNbrs = 0;
   for( e_i = 0; e_i < nEls; e_i++ ) {
      nCurNbrs = 0;
      nIncEls = self->nIncEls[dim][0][e_i];
      incEls = self->incEls[dim][0][e_i];
      for( inc_i = 0; inc_i < nIncEls; inc_i++ )
			nCurNbrs += self->nIncEls[0][dim][incEls[inc_i]];
			if( nCurNbrs > maxNbrs )
				maxNbrs = nCurNbrs;
	}

   ISet_Construct( nbrSet );
   ISet_SetMaxSize( nbrSet, maxNbrs );
   if( !self->nIncEls[dim][dim] ) {
      self->nIncEls[dim][dim] = Class_Array( self, int, nEls );
      memset( self->nIncEls[dim][dim], 0, nEls * sizeof(int) );
   }
   if( !self->incEls[dim][dim] ) {
      self->incEls[dim][dim] = Class_Array( self, int*, nEls );
      memset( self->incEls[dim][dim], 0, nEls * sizeof(int*) );
   }
   for( e_i = 0; e_i < nEls; e_i++ ) {
      nIncEls = self->nIncEls[dim][0][e_i];
      incEls = self->incEls[dim][0][e_i];
      for( inc_i = 0; inc_i < nIncEls; inc_i++ ) {
			nUpEls = self->nIncEls[0][dim][incEls[inc_i]];
			upEls = self->incEls[0][dim][incEls[inc_i]];
			for( inc_j = 0; inc_j < nUpEls; inc_j++ ) {
				if( upEls[inc_j] == e_i )
					continue;
				ISet_TryInsert( nbrSet, upEls[inc_j] );
			}
      }
      self->nIncEls[dim][dim][e_i] = ISet_GetSize( nbrSet );
      self->incEls[dim][dim][e_i] = Class_Rearray( self, self->incEls[dim][dim][e_i], int, self->nIncEls[dim][dim][e_i] );
      ISet_GetArray( nbrSet, self->incEls[dim][dim][e_i] );
      ISet_Clear( nbrSet );
   }
   ISet_Destruct( nbrSet );
}

void _IGraph_SetShadowDepth( void* _self, int depth ) {
   IGraph* self = (IGraph*)_self;
   int nNbrs, nDims;
   Sync* vertSync;
   ISet ghostSetObj, *ghostSet = &ghostSetObj;
   ISet mySetObj, *mySet = &mySetObj;
   IArray* isects;
   int nGhosts, *ghosts;
   int nLocals, nRemotes;
   const int* locals, *remotes;
   int *nNbrGhosts, **nbrGhosts;
   int *nBytes, *nRecvBytes;
   stgByte **bytes, **recvBytes;
   int nIncEls, *incEls;
   int **nLowEls, ***lowEls;
   int *nShdEls, **shdEls;
   int el, dom;
   int n_i, s_i, l_i, inc_i;
   int v_i, g_i, e_i, d_i;

   assert( self && depth >= 0 );
   assert( self->comm && self->nTDims );

   _MeshTopology_SetShadowDepth( self, depth );

   /* Build ghost set. */
   nDims = self->nDims;
   nNbrs = Comm_GetNumNeighbours( self->comm );
   vertSync = self->remotes[0];
   nGhosts = 0;
   for( n_i = 0; n_i < nNbrs; n_i++ )
      nGhosts += vertSync->nSnks[n_i] + vertSync->nSrcs[n_i];
   ISet_Construct( ghostSet );
   ISet_SetMaxSize( ghostSet, nGhosts );
   for( n_i = 0; n_i < nNbrs; n_i++ ) {
      for( s_i = 0; s_i < vertSync->nSnks[n_i]; s_i++ ) {
			el = Decomp_LocalToGlobal( self->locals[0], vertSync->snks[n_i][s_i] );
			ISet_TryInsert( ghostSet, el );
      }
      for( s_i = 0; s_i < vertSync->nSrcs[n_i]; s_i++ ) {
			el = Sync_RemoteToGlobal( vertSync, vertSync->srcs[n_i][s_i] );
			ISet_TryInsert( ghostSet, el );
      }
   }
   nGhosts = ISet_GetSize( ghostSet );
   ghosts = Class_Array( self, int, nGhosts );
   ISet_GetArray( ghostSet, ghosts );
   ISet_Destruct( ghostSet );

   /* Gather neighbouring ghost sets. */
   nNbrGhosts = Class_Array( self, int, nNbrs );
   nbrGhosts = Class_Array( self, int*, nNbrs );
   Comm_AllgatherInit( self->comm, nGhosts, nNbrGhosts, sizeof(int) );
   for( n_i = 0; n_i < nNbrs; n_i++ )
      nbrGhosts[n_i] = Class_Array( self, int, nNbrGhosts[n_i] );
   Comm_AllgatherBegin( self->comm, ghosts, (void**)nbrGhosts );
   Comm_AllgatherEnd( self->comm );
   Class_Free( self, ghosts );

   /* Build intersections. */
   ISet_Construct( mySet );
   ISet_SetMaxSize( mySet, Sync_GetNumDomains( vertSync ) );
   Decomp_GetLocals( self->locals[0], &nLocals, &locals );
   for( l_i = 0; l_i < nLocals; l_i++ )
      ISet_Insert( mySet, locals[l_i] );
   Sync_GetRemotes( self->remotes[0], &nRemotes, &remotes );
   for( l_i = 0; l_i < nRemotes; l_i++ )
      ISet_Insert( mySet, remotes[l_i] );
   isects = Class_Array( self, IArray, nNbrs );
   for( n_i = 0; n_i < nNbrs; n_i++ ) {
      IArray_Construct( isects + n_i );
      for( g_i = 0; g_i < nNbrGhosts[n_i]; g_i++ ) {
			if( ISet_Has( mySet, nbrGhosts[n_i][g_i] ) )
				IArray_Append( isects + n_i, nbrGhosts[n_i][g_i] );
		}
      Class_Free( self, nbrGhosts[n_i] );
   }

   /* Convert vertices to shadowed elements. */
   ISet_Clear( mySet );
   ISet_SetMaxSize( mySet, Decomp_GetNumLocals( self->locals[self->nDims] ) );
   for( n_i = 0; n_i < nNbrs; n_i++ ) {
      IArray_GetArray( isects + n_i, nNbrGhosts + n_i, (const int**)(nbrGhosts + n_i) );

      for( v_i = 0; v_i < nNbrGhosts[n_i]; v_i++ ) {
			dom = Sync_GlobalToDomain( vertSync, nbrGhosts[n_i][v_i] );

			for( inc_i = 0; inc_i < self->nIncEls[0][nDims][dom]; inc_i++ ) {
				el = self->incEls[0][nDims][dom][inc_i];

				if( el < Decomp_GetNumLocals( self->locals[nDims] ) ) {
					el = Decomp_LocalToGlobal( self->locals[nDims], el );
					ISet_TryInsert( mySet, el );
				}
			}
      }
      IArray_Destruct( isects + n_i );
      nNbrGhosts[n_i] = ISet_GetSize( mySet );
      nbrGhosts[n_i] = Class_Array( self, int, nNbrGhosts[n_i] );
      ISet_GetArray( mySet, nbrGhosts[n_i] );
      ISet_Clear( mySet );
   }
   Class_Free( self, isects );
   ISet_SetMaxSize( mySet, 0 );

   /* Transfer elements. */
   nShdEls = Class_Array( self, int, nNbrs );
   shdEls = Class_Array( self, int*, nNbrs );
   Comm_AlltoallInit( self->comm, nNbrGhosts, nShdEls, sizeof(int) );
   for( n_i = 0; n_i < nNbrs; n_i++ )
      shdEls[n_i] = Class_Array( self, int, nShdEls[n_i] );
   Comm_AlltoallBegin( self->comm, (const void**)nbrGhosts, (void**)shdEls );
   Comm_AlltoallEnd( self->comm );
   dom = 0;
   for( n_i = 0; n_i < nNbrs; n_i++ )
      dom += nShdEls[n_i];
   ISet_SetMaxSize( mySet, dom );
   for( n_i = 0; n_i < nNbrs; n_i++ ) {
      for( e_i = 0; e_i < nShdEls[n_i]; e_i++ )
			ISet_Insert( mySet, shdEls[n_i][e_i] );
			Class_Free( self, shdEls[n_i] );
   }
   nShdEls[0] = ISet_GetSize( mySet );
   shdEls[0] = Class_Array( self, int, nShdEls[0] );
   ISet_GetArray( mySet, shdEls[0] );
   ISet_Clear( mySet );
   ISet_SetMaxSize( mySet, 0 );
   qsort( shdEls[0], nShdEls[0], sizeof(int), IGraph_Cmp );
   IGraph_AddRemoteElements( self, nDims, nShdEls[0], shdEls[0] );
   Class_Free( self, shdEls[0] );

   /* Transfer lower level shadowed elements. */
   for( n_i = 0; n_i < nNbrs; n_i++ ) {
      for( e_i = 0; e_i < nNbrGhosts[n_i]; e_i++ ) {
			nbrGhosts[n_i][e_i] = Decomp_GlobalToLocal( self->locals[nDims], nbrGhosts[n_i][e_i] );
      }
   }
   nLowEls = Class_Array( self, int*, self->nTDims );
   lowEls = Class_Array( self, int**, self->nTDims );
	for( d_i = nDims - 1; d_i >= 0; d_i-- ) {
      if( !self->nIncEls[nDims][d_i] ) {
			nLowEls[d_i] = NULL;
			lowEls[d_i] = NULL;
			continue;
		}
      nLowEls[d_i] = Class_Array( self, int, nNbrs );
      lowEls[d_i] = Class_Array( self, int*, nNbrs );
      ISet_SetMaxSize( mySet, Decomp_GetNumLocals( self->locals[d_i] ) );

      for( n_i = 0; n_i < nNbrs; n_i++ ) {
			for( e_i = 0; e_i < nNbrGhosts[n_i]; e_i++ ) {
				nIncEls = self->nIncEls[nDims][d_i][nbrGhosts[n_i][e_i]];
				incEls = self->incEls[nDims][d_i][nbrGhosts[n_i][e_i]];

				for( inc_i = 0; inc_i < nIncEls; inc_i++ ) {
					if( incEls[inc_i] >= Decomp_GetNumLocals( self->locals[d_i] ) )
						continue;
					el = Decomp_LocalToGlobal( self->locals[d_i], incEls[inc_i] );
					ISet_TryInsert( mySet, el );
				}
			}
			nLowEls[d_i][n_i] = ISet_GetSize( mySet );
			lowEls[d_i][n_i] = Class_Array( self, int, nLowEls[d_i][n_i] );
			ISet_GetArray( mySet, lowEls[d_i][n_i] );
			ISet_Clear( mySet );
		}

      Comm_AlltoallInit( self->comm, nLowEls[d_i], nShdEls, sizeof(int) );
		for( n_i = 0; n_i < nNbrs; n_i++ )
			shdEls[n_i] = Class_Array( self, int, nShdEls[n_i] );
      Comm_AlltoallBegin( self->comm, (const void**)lowEls[d_i], (void**)shdEls );
      Comm_AlltoallEnd( self->comm );

      dom = 0;
      for( n_i = 0; n_i < nNbrs; n_i++ )
			dom += nShdEls[n_i];
      ISet_SetMaxSize( mySet, dom );
		for( n_i = 0; n_i < nNbrs; n_i++ ) {
			for( s_i = 0; s_i < nShdEls[n_i]; s_i++ ) {
				if( !Sync_TryGlobalToDomain( self->remotes[d_i], shdEls[n_i][s_i], &el ) ) {
					ISet_Insert( mySet, shdEls[n_i][s_i] );
				}
			}
			Class_Free( self, shdEls[n_i] );
      }
      nShdEls[0] = ISet_GetSize( mySet );
      shdEls[0] = Class_Array( self, int, nShdEls[0] );
      ISet_GetArray( mySet, shdEls[0] );
      ISet_Clear( mySet );
      qsort( shdEls[0], nShdEls[0], sizeof(int), IGraph_Cmp );
      IGraph_AddRemoteElements( self, d_i, nShdEls[0], shdEls[0] );
      Class_Free( self, shdEls[0] );
   }
   ISet_Destruct( mySet );
   Class_Free( self, shdEls );
   Class_Free( self, nShdEls );

   /* Transfer shadowed incidence. */
   nBytes = Class_Array( self, int, nNbrs );
   bytes = Class_Array( self, stgByte*, nNbrs );
   nRecvBytes = Class_Array( self, int, nNbrs );
   recvBytes = Class_Array( self, stgByte*, nNbrs );
   nLowEls[nDims] = nNbrGhosts;
   lowEls[nDims] = nbrGhosts;
   for( d_i = 0; d_i < nDims; d_i++ ) {
		if( !nLowEls[d_i] )
			continue;
		for( n_i = 0; n_i < nNbrs; n_i++ ) {
			for( e_i = 0; e_i < nLowEls[d_i][n_i]; e_i++ ) {
				lowEls[d_i][n_i][e_i] = Decomp_GlobalToLocal( self->locals[d_i], lowEls[d_i][n_i][e_i] );
			}
		}
   }
   for( d_i = nDims; d_i >= 0; d_i-- ) {
		if( !nLowEls[d_i] )
			continue;
		if( d_i == 0 ) {
			for( n_i = 0; n_i < nNbrs; n_i++ )
				Class_Free( self, lowEls[0][n_i] );
			Class_Free( self, lowEls[0] );
			Class_Free( self, nLowEls[0] );
			continue;
		}

		for( n_i = 0; n_i < nNbrs; n_i++ ) {
			IGraph_PickleIncidenceInit( self, d_i, nLowEls[d_i][n_i], lowEls[d_i][n_i], nBytes + n_i );
			bytes[n_i] = Class_Array( self, stgByte, nBytes[n_i] );
			IGraph_PickleIncidence( self, d_i, nLowEls[d_i][n_i], lowEls[d_i][n_i], bytes[n_i] );
			Class_Free( self, lowEls[d_i][n_i] );
		}
      Class_Free( self, nLowEls[d_i] );
      Class_Free( self, lowEls[d_i] );

      Comm_AlltoallInit( self->comm, nBytes, nRecvBytes, sizeof(stgByte) );
		for( n_i = 0; n_i < nNbrs; n_i++ )
			recvBytes[n_i] = Class_Array( self, stgByte, nRecvBytes[n_i] );
		Comm_AlltoallBegin( self->comm, (const void**)bytes, (void**)recvBytes );
		Comm_AlltoallEnd( self->comm );
      for( n_i = 0; n_i < nNbrs; n_i++ )
			Class_Free( self, bytes[n_i] );

      for( n_i = 0; n_i < nNbrs; n_i++ ) {
			IGraph_UnpickleIncidence( self, d_i, nRecvBytes[n_i], recvBytes[n_i] );
			Class_Free( self, recvBytes[n_i] );
      }
   }
   Class_Free( self, nBytes );
   Class_Free( self, bytes );
   Class_Free( self, lowEls );
   Class_Free( self, nLowEls );
   Class_Free( self, recvBytes );
   Class_Free( self, nRecvBytes );
}

void IGraph_Clear( void* self ) {
   IGraph_ClearDims( self );
   if( ((IGraph*)self)->comm )
      NewClass_RemoveRef( ((IGraph*)self)->comm );
   ((IGraph*)self)->comm = NULL;
	((NewClass*)self)->curAllocd = 0;
}

void IGraph_ClearDims( void* _self ) {
   IGraph* self = (IGraph*)_self;
   int d_i;

   IGraph_ClearElements( self );
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

   self->nDims = (MeshTopology_Dim)0;
   self->nTDims = 0;
   self->shadDepth = 0;
   self->locals = NULL;
   self->remotes = NULL;
   self->nIncEls = NULL;
   self->incEls = NULL;
}

void IGraph_ClearElements( void* _self ) {
   IGraph* self = (IGraph*)_self;
   int d_i;

   IGraph_ClearIncidence( self );
   for( d_i = 0; d_i < self->nTDims; d_i++ ) {
      if( self->locals )
			Decomp_ClearLocals( self->locals[d_i] );
      if( self->remotes && self->remotes[d_i] )
			Sync_ClearRemotes( self->remotes[d_i] );
      if( self->bndEls )
			Class_Free( self, self->bndEls[d_i] );
   }
   Class_Free( self, self->bndEls );
   Class_Free( self, self->nBndEls );
   self->bndEls = NULL;
   self->nBndEls = NULL;
}

void IGraph_ClearIncidence( void* _self ) {
   IGraph* self = (IGraph*)_self;
   int d_i, d_j, e_i;

   assert( self );
   for( d_i = 0; d_i < self->nTDims; d_i++ ) {
		for( d_j = 0; d_j < self->nTDims; d_j++ ) {
			if( self->nIncEls[d_i][d_j] ) {
				for( e_i = 0; e_i < Sync_GetNumDomains( self->remotes[d_i] ); e_i++ )
					Class_Free( self, self->incEls[d_i][d_j][e_i] );
				Class_Free( self, self->incEls[d_i][d_j] );
				Class_Free( self, self->nIncEls[d_i][d_j] );
				self->nIncEls[d_i][d_j] = NULL;
				self->incEls[d_i][d_j] = NULL;
			}
      }
   }
}

int IGraph_GetNumDims( const void* self ) {
   assert( self );
   return ((IGraph*)self)->nDims;
}

const Comm* IGraph_GetComm( const void* self ) {
   assert( self );
   return ((IGraph*)self)->comm;
}

Bool IGraph_HasDomain( const void* self, int dim ) {
   assert( self && dim < ((IGraph*)self)->nTDims );
   return Sync_GetNumDomains( ((IGraph*)self)->remotes[dim] ) ? 
      True : False;
}

const Sync* IGraph_GetDomain( const void* self, int dim ) {
   assert( self && dim < ((IGraph*)self)->nTDims );
   return ((IGraph*)self)->remotes[dim];
}

void IGraph_GetBoundaryElements( const void* self, int dim, int* nEls, const int** els ) {
   assert( Class_IsSuper( self, IGraph ) );
   assert( dim < ((IGraph*)self)->nTDims );
   assert( nEls );
   assert( els );

   *nEls = ((IGraph*)self)->nBndEls ? ((IGraph*)self)->nBndEls[dim] : 0;
   *els = ((IGraph*)self)->bndEls ? ((IGraph*)self)->bndEls[dim] : NULL;
}

Bool IGraph_HasIncidence( const void* self, int fromDim, int toDim ) {
   assert( self );
   assert( fromDim < ((IGraph*)self)->nTDims );
   assert( toDim < ((IGraph*)self)->nTDims );
   return ((IGraph*)self)->nIncEls[fromDim][toDim] ? True : False;
}

int IGraph_GetIncidenceSize( const void* self, int fromDim, int fromEl, int toDim ) {
   assert( self );
   assert( fromDim < ((IGraph*)self)->nTDims );
   assert( toDim < ((IGraph*)self)->nTDims );
   assert( fromEl < Sync_GetNumDomains( ((IGraph*)self)->remotes[fromDim] ) );
   return ((IGraph*)self)->nIncEls[fromDim][toDim][fromEl];
}

void _IGraph_GetIncidence( void* self, int fromDim, int fromEl, int toDim, IArray* inc ) {
   assert( self );
   assert( fromDim < ((IGraph*)self)->nTDims );
   assert( toDim < ((IGraph*)self)->nTDims );
   assert( fromEl < Sync_GetNumDomains( ((IGraph*)self)->remotes[fromDim] ) );
   assert( inc );

   IArray_SoftResize( inc, ((IGraph*)self)->nIncEls[fromDim][toDim][fromEl] );
   memcpy( inc->ptr, ((IGraph*)self)->incEls[fromDim][toDim][fromEl], IArray_GetSize( inc ) * sizeof(int) );
}

void IGraph_PrintIncidence( const void* _self, int fromDim, int toDim ) {
   IGraph* self = (IGraph*)_self;
   int nEls, global;
   int nIncEls, *incEls;
   int e_i, inc_i;

   assert( self );
   assert( toDim < self->nTDims );
   assert( fromDim < self->nTDims );

   nEls = Sync_GetNumDomains( self->remotes[fromDim] );
   printf( "Printing incidence for %d elements:\n", nEls );
   for( e_i = 0; e_i < nEls; e_i++ ) {
      global = Sync_DomainToGlobal( self->remotes[fromDim], e_i );
      nIncEls = self->nIncEls[fromDim][toDim][e_i];
      incEls = self->incEls[fromDim][toDim][e_i];
      printf( "   %d, %d incident elements:\n", global, nIncEls );
      for( inc_i = 0; inc_i < nIncEls; inc_i++ ) {
			printf( "      %d\n", incEls[inc_i] );
      }
   }
}

void IGraph_PickleIncidenceInit( IGraph* self, int dim, int nEls, int* els, int* nBytes ) {
   int size;
   int d_i, e_i;

   assert( self );
   assert( dim < self->nTDims );
   assert( !nEls || els );
   assert( nBytes );

   size = 1;
   for( e_i = 0; e_i < nEls; e_i++ ) {
      size += 1;
      for( d_i = 0; d_i < dim; d_i++ ) {
			size += 1;
			if( self->nIncEls[dim][d_i] )
				size += self->nIncEls[dim][d_i][els[e_i]];
		}
	}
   *nBytes = size * sizeof(int);
}

void IGraph_PickleIncidence( IGraph* self, int dim, int nEls, int* els, stgByte* bytes ) {
   Sync* sync;
   int curEntry, *entries;
   int nIncEls, *incEls;
   int inc_i, d_i, e_i;

   assert( self );
   assert( dim < self->nTDims );
   assert( !nEls || els );
   assert( bytes );

   sync = self->remotes[dim];
   entries = (int*)bytes;
   entries[0] = nEls;
   curEntry = 1;
   for( e_i = 0; e_i < nEls; e_i++ ) {
      entries[curEntry++] = Sync_DomainToGlobal( sync, els[e_i] );
      for( d_i = 0; d_i < dim; d_i++ ) {
			if( !self->nIncEls[dim][d_i] ) {
				entries[curEntry++] = 0;
				continue;
			}
			nIncEls = self->nIncEls[dim][d_i][els[e_i]];
			entries[curEntry++] = nIncEls;
			incEls = self->incEls[dim][d_i][els[e_i]];
			for( inc_i = 0; inc_i < nIncEls; inc_i++ ) {
				entries[curEntry++] = Sync_DomainToGlobal( self->remotes[d_i], incEls[inc_i] );
			}
      }
   }
}

void IGraph_UnpickleIncidence( IGraph* self, int dim, int nBytes, stgByte* bytes ) {
   Sync* sync;
   int nEls, el;
   int nIncEls, **incEls;
   int curEntry, *entries;
   int inc_i, e_i, d_i;

   assert( self );
   assert( dim < self->nTDims );
   assert( nBytes && bytes );

   sync = self->remotes[dim];
   entries = (int*)bytes;
   nEls = entries[0];
   curEntry = 1;
   for( e_i = 0; e_i < nEls; e_i++ ) {
      el = Sync_GlobalToDomain( sync, entries[curEntry++] );
      for( d_i = 0; d_i < dim; d_i++ ) {
			nIncEls = entries[curEntry++];
			if( !self->nIncEls[dim][d_i] ) {
				if( !nIncEls )
					continue;
				self->nIncEls[dim][d_i] = Class_Array( self, int, Sync_GetNumDomains( sync ) );
				memset( self->nIncEls[dim][d_i], 0, sizeof(int) * Sync_GetNumDomains( sync ) );
			}
			self->nIncEls[dim][d_i][el] = nIncEls;
			if( !nIncEls ) {
				if( self->incEls[dim][d_i] )
					Class_Free( self, self->incEls[dim][d_i][el] );
				continue;
			}
			if( !self->incEls[dim][d_i] ) {
				self->incEls[dim][d_i] = Class_Array( self, int*, Sync_GetNumDomains( sync ) );
				memset( self->incEls[dim][d_i], 0, 
				Sync_GetNumDomains( sync ) * sizeof(int*) );
			}
			incEls = self->incEls[dim][d_i];
			incEls[el] = Class_Rearray( self, incEls[el], int, nIncEls );
			for( inc_i = 0; inc_i < nIncEls; inc_i++ ) {
				incEls[el][inc_i] = Sync_GlobalToDomain( self->remotes[d_i], entries[curEntry++] );
			}
      }
   }
}

int IGraph_Cmp( const void* l, const void* r ) {
   assert( *(int*)l != *(int*)r );
   return (*(int*)l < *(int*)r) ? -1 : 1;
}
