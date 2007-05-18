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
** $Id: Decomp.c 3952 2007-01-09 06:24:06Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdlib.h>
#include <mpi.h>
#include "StGermain/Base/Base.h"
#include "types.h"
#include "Decomp.h"
#include "StGermain/Base/Foundation/ClassDef.h"


void Decomp_Update( Decomp* self );


void _Decomp_Construct( void* _self ) {
   Decomp* self = (Decomp*)_self;

   _NewClass_Construct( self );
   self->mpiComm = MPI_COMM_WORLD;
   self->nGlobals = 0;
   self->locals = &self->localsObj;
   IArray_Init( self->locals );
   self->inv = &self->invObj;
   IMap_Init( self->inv );
}

void _Decomp_Destruct( void* _self ) {
   Decomp* self = (Decomp*)_self;

   Decomp_Clear( self );
   IArray_Destruct( self->locals );
   IMap_Destruct( self->inv );
   _NewClass_Destruct( self );
}

void _Decomp_Copy( void* _self, const void* _op ) {
   Decomp* self = (Decomp*)_self;
   const Decomp* op = (const Decomp*)_op;

   assert( self && op );
   self->mpiComm = op->mpiComm;
   self->nGlobals = op->nGlobals;
   IArray_Copy( self->locals, op->locals );
   IMap_Copy( self->inv, op->inv );
}

SizeT _Decomp_CalcMem( const void* _self, PtrMap* ptrs ) {
   const Decomp* self = (const Decomp*)_self;
   SizeT mem;

   if( PtrMap_Find( ptrs, (void*)self ) )
      return 0;
   mem = _NewClass_CalcMem( self, ptrs );
   mem += NewClass_CalcMem( self->locals, ptrs );
   mem += NewClass_CalcMem( self->inv, ptrs );
   return mem;
}

void Decomp_SetMPIComm( void* _self, MPI_Comm mpiComm ) {
   Decomp* self = (Decomp*)_self;

   assert( self );
   self->mpiComm = mpiComm;
}

void Decomp_SetLocals( void* _self, int nLocals, const int* locals ) {
   Decomp* self = (Decomp*)_self;
   int l_i;

   assert( self && (!nLocals || locals) );
   IArray_Set( self->locals, nLocals, locals );
   Decomp_Update( self );
   IMap_Clear( self->inv );
   IMap_SetMaxSize( self->inv, nLocals );
   for( l_i = 0; l_i < nLocals; l_i++ )
      IMap_Insert( self->inv, locals[l_i], l_i );
}

void Decomp_AddLocals( void* _self, int nLocals, const int* locals ) {
   Decomp* self = (Decomp*)_self;
   int l_i;

   assert( self && (!nLocals || locals) );
   IMap_SetMaxSize( self->inv, IArray_GetSize( self->locals ) + nLocals );
   IArray_Add( self->locals, nLocals, locals );
   Decomp_Update( self );
   for( l_i = 0; l_i < nLocals; l_i++ )
      IMap_Insert( self->inv, locals[l_i], l_i );
}

void Decomp_RemoveLocals( void* _self, int nLocals, const int* locals, IMap* map ) {
   Decomp* self = (Decomp*)_self;
   int oldSize;
   int l_i;

   assert( self && (!nLocals || locals) );
   oldSize = IArray_GetSize( self->locals );
   IArray_Remove( self->locals, nLocals, locals, map );
   Decomp_Update( self );
   for( l_i = 0; l_i < nLocals; l_i++ )
      IMap_Remove( self->inv, locals[l_i] );
   IMap_SetMaxSize( self->inv, oldSize - nLocals );
}

void Decomp_Clear( void* _self ) {
   Decomp* self = (Decomp*)_self;

   Decomp_ClearLocals( self );
   self->mpiComm = MPI_COMM_WORLD;
   self->nGlobals = 0;
}

void Decomp_ClearLocals( void* _self ) {
   Decomp* self = (Decomp*)_self;

   assert( self );
   IArray_Clear( self->locals );
   IMap_Clear( self->inv );
}

MPI_Comm Decomp_GetComm( const void* self ) {
   assert( self );
   return ((Decomp*)self)->mpiComm;
}

int Decomp_GetNumGlobals( const void* self ) {
   assert( self );
   return ((Decomp*)self)->nGlobals;
}

int Decomp_GetNumLocals( const void* self ) {
   assert( self );
   return IArray_GetSize( ((Decomp*)self)->locals );
}

void Decomp_GetLocals( const void* self, int* nLocals, const int** locals ) {
   assert( self );
   *nLocals = IArray_GetSize( ((Decomp*)self)->locals );
   *locals = IArray_GetPtr( ((Decomp*)self)->locals );
}

int Decomp_LocalToGlobal( const void* self, int local ) {
   assert( self && local < IArray_GetSize( ((Decomp*)self)->locals ) );
   return IArray_GetPtr( ((Decomp*)self)->locals )[local];
}

int Decomp_GlobalToLocal( const void* self, int global ) {
   assert( self && global < ((Decomp*)self)->nGlobals );
   return IMap_Map( ((Decomp*)self)->inv, global );
}

Bool Decomp_TryGlobalToLocal( const void* self, int global, int* local ) {
   assert( self && global < ((Decomp*)self)->nGlobals );
   return IMap_TryMap( ((Decomp*)self)->inv, global, local );
}

void Decomp_Update( Decomp* self ) {
   int nLocals;

   assert( self );
   if( self->mpiComm ) {
      nLocals = IArray_GetSize( self->locals );
      insist( MPI_Allreduce( &nLocals, &self->nGlobals, 1, MPI_INT, MPI_SUM, 
			     self->mpiComm ), == MPI_SUCCESS );
   }
   else
      self->nGlobals = 0;
}
