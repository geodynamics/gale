/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003, Victorian Partnership for Advanced Computing (VPAC) Ltd, 
110 Victoria Street, Melbourne, 3053, Australia.
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
** $Id: IArray.c 3952 2007-01-09 06:24:06Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "StGermain/Base/Foundation/Foundation.h"
#include "types.h"
#include "ISet.h"
#include "IMap.h"
#include "IArray.h"
#include "StGermain/Base/Foundation/ClassDef.h"


int IArray_Cmp( const void* l, const void* r );


void _IArray_Construct( void* _self ) {
   IArray* self = (IArray*)_self;

   _NewClass_Construct( self );
   self->delta = 100;
   self->maxSize = 0;
   self->size = 0;
   self->ptr = NULL;
}

void _IArray_Destruct( void* self ) {
   IArray_Clear( self );
   _NewClass_Destruct( self );
}

void _IArray_Copy( void* _self, const void* _op ) {
   IArray* self = (IArray*)_self;
   const IArray* op = (const IArray*)_op;

   assert( self );
   self->delta = op->delta;
   self->maxSize = op->maxSize;
   self->size = op->size;
   self->ptr = Class_Array( self, int, self->size );
   memcpy( self->ptr, op->ptr, sizeof(int) * self->size );
}

void IArray_SetDelta( void* _self, int delta ) {
   IArray* self = (IArray*)_self;

   assert( self );
   self->delta = delta;
}

void IArray_Resize( void* _self, int size ) {
   IArray* self = (IArray*)_self;

   assert( self && self->delta );
   self->maxSize = size / self->delta + ((size % self->delta) ? 1 : 0);
   self->maxSize *= self->delta;
   self->size = size;
   self->ptr = Class_Rearray( self, self->ptr, int, self->maxSize );
}

void IArray_Set( void* _self, int nItms, const int* itms ) {
   IArray* self = (IArray*)_self;

   assert( self && (!nItms || itms) );
   IArray_Resize( self, nItms );
   memcpy( self->ptr, itms, nItms * sizeof(int) );
}

void IArray_Add( void* _self, int nItms, const int* itms ) {
   IArray* self = (IArray*)_self;
   int oldSize;

   assert( self && (!nItms || itms) );
   oldSize = self->size;
   IArray_Resize( self, oldSize + nItms );
   memcpy( self->ptr + oldSize, itms, nItms * sizeof(int) );
}

void IArray_Remove( void* _self, int nItms, const int* locals, IMap* map ) {
   IArray* self = (IArray*)_self;
   ISet toRemObj, *toRem = &toRemObj;
   int* ord, pos;
   int i_i;

   assert( self );
   ISet_Init( toRem );
   ISet_UseArray( toRem, nItms, locals );
   ord = Class_Array( self, int, ISet_GetSize( toRem ) );
   memcpy( ord, locals, nItms * sizeof(int) );
   qsort( ord, nItms, sizeof(int), IArray_Cmp );
   IMap_Clear( map );
   IMap_SetMaxSize( map, nItms );
   for( i_i = 0, pos = self->size - 1; 
	i_i < nItms && pos > ord[i_i];
	i_i++, pos-- )
   {
      while( ISet_Has( toRem, pos ) && pos > ord[i_i] )
	 pos--;
      if( pos <= ord[i_i] )
	 break;
      self->ptr[ord[i_i]] = self->ptr[pos];
      IMap_Insert( map, pos, ord[i_i] );
   }
   if( IMap_GetSize( map ) < nItms )
      IMap_SetMaxSize( map, IMap_GetSize( map ) );
   ISet_Destruct( toRem );
   Class_Free( self, ord );

   IArray_Resize( self, self->size - nItms );
}

void IArray_Append( void* _self, int itm ) {
   IArray* self = (IArray*)_self;

   assert( self );
   if( self->size == self->maxSize )
      IArray_Resize( self, self->size + 1 );
   else
      self->size++;
   self->ptr[self->size - 1] = itm;
}

void IArray_Clear( void* self ) {
   IArray_Resize( self, 0 );
}

int IArray_GetSize( const void* self ) {
   assert( self );
   return ((IArray*)self)->size;
}

const int* IArray_GetPtr( const void* self ) {
   assert( self );
   return ((IArray*)self)->ptr;
}

int IArray_Cmp( const void* l, const void* r ) {
   assert( *(int*)l != *(int*)r );
   return (*(int*)l > *(int*)r) ? 1 : -1;
}
