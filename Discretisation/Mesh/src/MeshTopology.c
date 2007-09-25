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


void _MeshTopology_Init( void* _self ) {
   MeshTopology* self = (MeshTopology*)_self;

   _NewClass_Init( self );
   self->nDims = 0;
   self->nTDims = 0;
   self->shadDepth = 0;
   self->comm = NULL;
}

void _MeshTopology_Destruct( void* _self ) {
   MeshTopology* self = (MeshTopology*)_self;

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

   if( PtrMap_Find( ptrs, (void*)self ) )
      return 0;
   mem = _NewClass_CalcMem( self, ptrs );
   return mem;
}

void _MeshTopology_SetNumDims( void* _self, int nDims ) {
   MeshTopology* self = (MeshTopology*)_self;

   self->nDims = nDims;
   if( self->nDims > 0 )
      self->nTDims = nDims + 1;
}

void _MeshTopology_SetComm( void* _self, const Comm* comm ) {
   MeshTopology* self = (MeshTopology*)_self;

   NewClass_RemoveRef( self->comm );
   self->comm = (Comm*)comm;
   if( comm )
      NewClass_AddRef( (Comm*)comm );
}

void _MeshTopology_SetShadowDepth( void* _self, int depth ) {
   MeshTopology* self = (MeshTopology*)_self;

   assert( self && depth >= 0 );
   assert( self->comm && self->nTDims );

   self->shadDepth = depth;
}

int MeshTopology_GetNumDims( const void* self ) {
   assert( self );
   return ((MeshTopology*)self)->nDims;
}

const Comm* MeshTopology_GetComm( const void* self ) {
   assert( self );
   return ((MeshTopology*)self)->comm;
}
