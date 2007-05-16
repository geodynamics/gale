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
** $Id: Sync.c 3952 2007-01-09 06:24:06Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdlib.h>
#include <string.h>
#include "StGermain/Base/Base.h"
#include "types.h"
#include "Decomp.h"
#include "Sync.h"
#include "StGermain/Base/Foundation/ClassDef.h"


void Sync_UpdateTables( Sync* self );
void Sync_UpdateShared( Sync* self );
void Sync_UpdateOwners( Sync* self );
void Sync_ClearTables( Sync* self );


void _Sync_Construct( void* _self ) {
   Sync* self = (Sync*)_self;

   _NewClass_Construct( self );
   self->nDomains = 0;
   self->decomp = NULL;
   self->remotes = &self->remotesObj;
   IArray_Init( self->remotes );
   self->owners = NULL;
   self->nShared = 0;
   self->shared = NULL;
   self->nSharers = NULL;
   self->sharers = NULL;
   self->gr = &self->grObj;
   IMap_Init( self->gr );
   self->ls = &self->lsObj;
   IMap_Init( self->ls );
   self->nSrcs = NULL;
   self->srcs = NULL;
   self->nSnks = NULL;
   self->snks = NULL;
}

void _Sync_Destruct( void* _self ) {
   Sync* self = (Sync*)_self;

   Sync_Clear( self );
   IArray_Destruct( self->remotes );
   IMap_Destruct( self->gr );
   IMap_Destruct( self->ls );
   _NewClass_Destruct( self );
}

void _Sync_Copy( void* _self, const void* op ) {
   /*Sync* self = (Sync*)_self;*/

   /* TODO: Method body goes here */
   assert( 0 );
}

SizeT _Sync_CalcMem( const void* _self, PtrMap* ptrs ) {
   const Sync* self = (const Sync*)_self;
   SizeT mem;

   if( PtrMap_Find( ptrs, (void*)self ) )
      return 0;
   mem = _NewClass_CalcMem( self, ptrs );
   if( self->decomp )
      mem += NewClass_CalcMem( self->decomp, ptrs );
   mem += NewClass_CalcMem( self->remotes, ptrs );
   mem += NewClass_CalcMem( self->gr, ptrs );
   mem += NewClass_CalcMem( self->ls, ptrs );
   return mem;
}

void Sync_SetDecomp( void* _self, const Decomp* decomp ) {
   Sync* self = (Sync*)_self;

   Sync_Clear( self );
   self->decomp = (Decomp*)decomp;
   if( self->decomp ) {
      NewClass_AddRef( self->decomp );
      self->nDomains = Decomp_GetNumLocals( self->decomp );
      Sync_UpdateTables( self );
      Sync_UpdateShared( self );
      Sync_UpdateOwners( self );
   }
}

void Sync_SetRemotes( void* _self, int nRemotes, const int* remotes ) {
   Sync* self = (Sync*)_self;

   Sync_ClearRemotes( self );
   Sync_AddRemotes( self, nRemotes, remotes );
}

void Sync_AddRemotes( void* _self, int nRemotes, const int* remotes ) {
   Sync* self = (Sync*)_self;
   int r_i;

   assert( self && self->decomp );
   IMap_SetMaxItems( self->gr, IArray_GetSize( self->remotes ) + nRemotes );
   IArray_Add( self->remotes, nRemotes, remotes );
   for( r_i = 0; r_i < nRemotes; r_i++ )
      IMap_Insert( self->gr, remotes[r_i], r_i );
   Sync_UpdateTables( self );
   Sync_UpdateShared( self );
   Sync_UpdateOwners( self );
   self->nDomains += nRemotes;
}

void Sync_RemoveRemotes( void* _self, int nRemotes, const int* remotes, IMap* map ) {
   Sync* self = (Sync*)_self;

   assert( self && self->decomp );
   IArray_Remove( self->remotes, nRemotes, remotes, map );
   Sync_UpdateTables( self );
   Sync_UpdateShared( self );
   Sync_UpdateOwners( self );
   self->nDomains -= nRemotes;
}

void Sync_Clear( void* _self ) {
   Sync* self = (Sync*)_self;

   Sync_ClearRemotes( self );
   if( self->decomp )
      NewClass_RemoveRef( self->decomp );
   self->decomp = NULL;
   self->nDomains = 0;
}

void Sync_ClearRemotes( void* _self ) {
   Sync* self = (Sync*)_self;
   int s_i;

   Sync_ClearTables( self );
   IArray_Clear( self->remotes );
   for( s_i = 0; s_i < self->nShared; s_i++ )
      Class_Free( self, self->sharers[s_i] );
   Class_Free( self, self->shared );
   Class_Free( self, self->nSharers );
   Class_Free( self, self->sharers );
   Class_Free( self, self->owners );
   IMap_Clear( self->gr );
   IMap_Clear( self->ls );
   if( self->decomp )
      self->nDomains = Decomp_GetNumLocals( self->decomp );
   else
      self->nDomains = 0;
   self->nShared = 0;
   self->shared = NULL;
   self->nSharers = NULL;
   self->sharers = NULL;
   self->owners = NULL;
}

const Decomp* Sync_GetDecomp( const void* self ) {
   assert( self );
   return ((Sync*)self)->decomp;
}

int Sync_GetNumRemotes( const void* self ) {
   assert( self );
   return IArray_GetSize( ((Sync*)self)->remotes );
}

int Sync_GetNumDomains( const void* self ) {
   assert( self );
   return ((Sync*)self)->nDomains;
}

int Sync_GetNumShared( const void* self ) {
   assert( self );
   return ((Sync*)self)->nShared;
}

int Sync_GetNumSharers( const void* self, int shared ) {
   assert( self && shared < ((Sync*)self)->nShared );
   assert( ((Sync*)self)->nSharers );
   return ((Sync*)self)->nSharers[shared];
}

void Sync_GetRemotes( const void* self, int* nRemotes, const int** remotes ) {
   assert( self );
   *nRemotes = IArray_GetSize( ((Sync*)self)->remotes );
   *remotes = IArray_GetPtr( ((Sync*)self)->remotes );
}

int Sync_GetOwner( const void* self, int remote ) {
   assert( self && remote < IArray_GetSize( ((Sync*)self)->remotes ) );
   assert( ((Sync*)self)->owners );
   return ((Sync*)self)->owners[remote];
}

void Sync_GetShared( const void* self, int* nShared, const int** shared ) {
   assert( self );
   assert( !((Sync*)self)->nShared || shared );
   if( nShared )
      *nShared = ((Sync*)self)->nShared;
   *shared = ((Sync*)self)->shared;
}

void Sync_GetSharers( const void* self, int shared, int* nSharers, const int** sharers ) {
   assert( self && shared < ((Sync*)self)->nShared );
   assert( ((Sync*)self)->nSharers && sharers );
   if( nSharers )
      *nSharers = ((Sync*)self)->nSharers[shared];
   *sharers = ((Sync*)self)->sharers[shared];
}

int Sync_RemoteToGlobal( const void* self, int remote ) {
   assert( self && remote < IArray_GetSize( ((Sync*)self)->remotes ) );
   return IArray_GetPtr( ((Sync*)self)->remotes )[remote];
}

Bool Sync_GlobalToRemote( const void* self, int global, int* remote ) {
   assert( self && global < Decomp_GetNumGlobals( ((Sync*)self)->decomp ) );
   return IMap_Try( ((Sync*)self)->gr, global, remote );
}

int Sync_DomainToGlobal( const void* self, int domain ) {
   assert( self );
   assert( domain < Decomp_GetNumLocals( ((Sync*)self)->decomp ) + 
	   IArray_GetSize( ((Sync*)self)->remotes ) );
   if( domain < Decomp_GetNumLocals( ((Sync*)self)->decomp ) )
      return Decomp_LocalToGlobal( ((Sync*)self)->decomp, domain );
   else {
      return IArray_GetPtr( ((Sync*)self)->remotes )
	 [domain - Decomp_GetNumLocals( ((Sync*)self)->decomp )];
   }
}

Bool Sync_GlobalToDomain( const void* self, int global, int* domain ) {
   assert( self );
   assert( global < Decomp_GetNumGlobals( ((Sync*)self)->decomp ) );
   if( !Decomp_GlobalToLocal( ((Sync*)self)->decomp, global, domain ) ) {
      if( IMap_Try( ((Sync*)self)->gr, global, domain ) ) {
	 *domain += Decomp_GetNumLocals( ((Sync*)self)->decomp );
	 return True;
      }
   }
   else
      return True;
}

int Sync_SharedToLocal( const void* self, int shared ) {
   assert( self && shared < ((Sync*)self)->nShared );
   return ((Sync*)self)->shared[shared];
}

Bool Sync_LocalToShared( const void* self, int local, int* shared ) {
   assert( self && local < Decomp_GetNumLocals( ((Sync*)self)->decomp ) );
   assert( shared );
   return IMap_Try( ((Sync*)self)->ls, local, shared );
}

void Sync_SyncArray( const void* _self, 
		     const void* local, size_t localStride, 
		     const void* remote, size_t remoteStride, 
		     size_t itmSize )
{
   Sync* self = (Sync*)_self;
   const Comm* comm;
   int nNbrs;
   int* nSrcs;
   StgByte **srcs, **snks;
   int n_i, s_i;

   assert( self );
   comm = Decomp_GetComm( self->decomp );
   nNbrs = Comm_GetNumNeighbours( comm );
   snks = Class_Array( self, StgByte*, nNbrs );
   for( n_i = 0; n_i < nNbrs; n_i++ ) {
      snks[n_i] = Class_Array( self, StgByte, self->nSnks[n_i] * itmSize );
      for( s_i = 0; s_i < self->nSnks[n_i]; s_i++ ) {
	 memcpy( snks[n_i] + s_i * itmSize, 
		 (StgByte*)local + self->snks[n_i][s_i] * localStride, 
		 itmSize );
      }
   }

   nSrcs = Class_Array( self, int, nNbrs );
   Comm_AlltoallInit( comm, self->nSnks, nSrcs, itmSize );
   srcs = Class_Array( self, StgByte*, nNbrs );
   for( n_i = 0; n_i < nNbrs; n_i++ )
      srcs[n_i] = Class_Array( self, StgByte, nSrcs[n_i] * itmSize );
   Comm_AlltoallBegin( comm, (const void**)snks, (void**)srcs );
   Comm_AlltoallEnd( comm );
   Class_Free( self, nSrcs );
   for( n_i = 0; n_i < nNbrs; n_i++ )
      Class_Free( self, snks[n_i] );
   Class_Free( self, snks );

   for( n_i = 0; n_i < nNbrs; n_i++ ) {
      for( s_i = 0; s_i < self->nSrcs[n_i]; s_i++ ) {
	 memcpy( (StgByte*)remote + self->srcs[n_i][s_i] * remoteStride, 
		 srcs[n_i] + s_i * itmSize, 
		 itmSize );
      }
      Class_Free( self, srcs[n_i] );
   }
   Class_Free( self, srcs );
}

void Sync_UpdateTables( Sync* self ) {
   const Comm* comm;
   int nNbrs;
   ISet theirRemsObj, *theirRems = &theirRemsObj;
   ISet myLocsObj, *myLocs = &myLocsObj;
   ISet isectObj, *isect = &isectObj;
   int nLocals, nRems;
   const int *locals, *rems;
   int* nAllRems, **allRems;
   int* nFnds, **fnds;
   int n_i, s_i;

   Sync_ClearTables( self );
   comm = Decomp_GetComm( self->decomp );
   nNbrs = Comm_GetNumNeighbours( comm );
   ISet_Init( theirRems );
   ISet_Init( myLocs );
   ISet_Init( isect );
   Decomp_GetLocals( self->decomp, &nLocals, &locals );
   ISet_UseArray( myLocs, nLocals, locals );
   nRems = IArray_GetSize( self->remotes );
   rems = IArray_GetPtr( self->remotes );

   nAllRems = Class_Array( self, int, nNbrs );
   allRems = Class_Array( self, int*, nNbrs );
   Comm_AllgatherInit( comm, nRems, nAllRems, sizeof(int) );
   for( n_i = 0; n_i < nNbrs; n_i++ )
      allRems[n_i] = Class_Array( self, int, nAllRems[n_i] );
   Comm_AllgatherBegin( comm, rems, (void**)allRems );
   Comm_AllgatherEnd( comm );

   for( n_i = 0; n_i < nNbrs; n_i++ ) {
      ISet_UseArray( theirRems, nAllRems[n_i], allRems[n_i] );
      ISet_Isect( myLocs, theirRems, isect );
      ISet_Clear( theirRems );
      allRems[n_i] = Class_Rearray( self, allRems[n_i], int, 
				    ISet_GetNumItems( isect ) );
      ISet_GetArray( isect, nAllRems + n_i, allRems[n_i] );
   }

   nFnds = Class_Array( self, int, nNbrs );
   fnds = Class_Array( self, int*, nNbrs );
   Comm_AlltoallInit( comm, nAllRems, nFnds, sizeof(int) );
   for( n_i = 0; n_i < nNbrs; n_i++ )
      fnds[n_i] = Class_Array( self, int, nFnds[n_i] );
   Comm_AlltoallBegin( comm, (const void**)allRems, (void**)fnds );
   Comm_AlltoallEnd( comm );

   ISet_Destruct( theirRems );
   ISet_Destruct( myLocs );
   ISet_Destruct( isect );

   self->nSrcs = nFnds;
   self->srcs = fnds;
   self->nSnks = nAllRems;
   self->snks = allRems;

   for( n_i = 0; n_i < nNbrs; n_i++ ) {
      for( s_i = 0; s_i < self->nSrcs[n_i]; s_i++ ) {
	 insist( Sync_GlobalToRemote( self, self->srcs[n_i][s_i], 
				      self->srcs[n_i] + s_i ), == True );
      }
      for( s_i = 0; s_i < self->nSnks[n_i]; s_i++ ) {
	 insist( Decomp_GlobalToLocal( self->decomp, self->snks[n_i][s_i], 
				       self->snks[n_i] + s_i ), == True );
      }
   }
}

void Sync_UpdateShared( Sync* self ) {
   ISet sharedSetObj, *sharedSet = &sharedSetObj;
   const Comm* comm;
   int shared, nNbrs;
   int n_i, s_i;

   assert( self );
   for( s_i = 0; s_i < self->nShared; s_i++ )
      Class_Free( self, self->sharers[s_i] );

   comm = Decomp_GetComm( self->decomp );
   nNbrs = Comm_GetNumNeighbours( comm );
   ISet_Init( sharedSet );
   for( n_i = 0; n_i < nNbrs; n_i++ ) {
      for( s_i = 0; s_i < self->nSnks[n_i]; s_i++ )
	 ISet_TryInsert( sharedSet, self->snks[n_i][s_i] );
   }
   self->nShared = ISet_GetNumItems( sharedSet );
   self->shared = Class_Rearray( self, self->shared, int, self->nShared );
   ISet_GetArray( sharedSet, NULL, self->shared );
   ISet_Destruct( sharedSet );

   IMap_Clear( self->ls );
   IMap_SetMaxItems( self->ls, self->nShared );
   for( s_i = 0; s_i < self->nShared; s_i++ )
      IMap_Insert( self->ls, self->shared[s_i], s_i );

   self->nSharers = Class_Rearray( self, self->nSharers, int, self->nShared );
   memset( self->nSharers, 0, self->nShared * sizeof(int) );
   for( n_i = 0; n_i < nNbrs; n_i++ ) {
      for( s_i = 0; s_i < self->nSnks[n_i]; s_i++ ) {
	 shared = IMap_Map( self->ls, self->snks[n_i][s_i] );
	 self->nSharers[shared]++;
      }
   }
   self->sharers = Class_Rearray( self, self->sharers, int*, self->nShared );
   for( s_i = 0; s_i < self->nShared; s_i++ )
      self->sharers[s_i] = Class_Array( self, int, self->nSharers[s_i] );
   memset( self->nSharers, 0, self->nShared * sizeof(int) );
   for( n_i = 0; n_i < nNbrs; n_i++ ) {
      for( s_i = 0; s_i < self->nSnks[n_i]; s_i++ ) {
	 shared = IMap_Map( self->ls, self->snks[n_i][s_i] );
	 self->sharers[shared][self->nSharers[shared]++] = n_i;
      }
   }
}

void Sync_UpdateOwners( Sync* self ) {
   const Comm* comm;
   int nRemotes, nNbrs;
   int n_i, s_i;

   assert( self );
   comm = Decomp_GetComm( self->decomp );
   nNbrs = Comm_GetNumNeighbours( comm );
   nRemotes = IArray_GetSize( self->remotes );
   self->owners = Class_Rearray( self, self->owners, int, nRemotes );
   for( n_i = 0; n_i < nNbrs; n_i++ ) {
      for( s_i = 0; s_i < self->nSrcs[n_i]; s_i++ )
	 self->owners[self->srcs[n_i][s_i]] = n_i;
   }
}

void Sync_ClearTables( Sync* self ) {
   const Comm* comm;
   int n_i;

   assert( self );
   if( self->decomp ) {
      comm = Decomp_GetComm( self->decomp );
      if( comm ) {
	 for( n_i = 0; n_i < Comm_GetNumNeighbours( comm ); n_i++ ) {
	    if( self->srcs )
	       Class_Free( self, self->srcs[n_i] );
	    if( self->snks )
	       Class_Free( self, self->snks[n_i] );
	 }
      }
   }
   Class_Free( self, self->nSrcs );
   Class_Free( self, self->srcs );
   Class_Free( self, self->nSnks );
   Class_Free( self, self->snks );
   self->nSrcs = NULL;
   self->srcs = NULL;
   self->nSnks = NULL;
   self->snks = NULL;
}
