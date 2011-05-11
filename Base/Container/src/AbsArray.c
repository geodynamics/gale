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
** $Id: AbsArray.c 3952 2007-01-09 06:24:06Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <StGermain/Base/Foundation/Foundation.h>

#include "types.h"
#include "AbsArray.h"
#include "StGermain/Base/Foundation/ClassDef.h"


void _AbsArray_Init( void* _self ) {
   AbsArray* self = Class_Cast( _self, AbsArray );

   _NewClass_Init( self );
   self->itmSize = 0;
   self->delta = 100;
   self->maxSize = 0;
   self->size = 0;
   self->ptr = NULL;
}

void _AbsArray_Destruct( void* self ) {
   AbsArray_Clear( self );
   _NewClass_Destruct( self );
}

void _AbsArray_Copy( void* _self, const void* _op ) {
   AbsArray* self = Class_Cast( _self, AbsArray );
   const AbsArray* op = Class_ConstCast( _op, AbsArray );

   self->itmSize = op->itmSize;
   self->delta = op->delta;
   self->maxSize = op->maxSize;
   self->size = op->size;
   self->ptr = Class_Array( self, stgByte, self->size * self->itmSize );
   memcpy( self->ptr, op->ptr, self->itmSize * self->size );
}

void AbsArray_SetItemSize( void* _self, int itmSize ) {
   AbsArray* self = Class_Cast( _self, AbsArray );

   AbsArray_Clear( self );
   self->itmSize = itmSize;
}

void AbsArray_SetDelta( void* _self, int delta ) {
   AbsArray* self = Class_Cast( _self, AbsArray );

   self->delta = delta;
}

void AbsArray_Resize( void* _self, int size ) {
   AbsArray* self = Class_Cast( _self, AbsArray );

   assert( self->delta );
   self->maxSize = size / self->delta + ((size % self->delta) ? 1 : 0);
   self->maxSize *= self->delta;
   self->size = size;
   self->ptr = Class_Rearray( self, self->ptr, stgByte, self->maxSize * self->itmSize );
}

void AbsArray_SoftResize( void* _self, int size ) {
   AbsArray* self = (AbsArray*)_self;

   assert( self && self->delta );
   if( size > self->maxSize )
      AbsArray_Resize( self, size );
   else
      self->size = size;
}

void AbsArray_Offset( void* _self, int delta ) {
   AbsArray* self = Class_Cast( _self, AbsArray );

   if ( delta > 0 ) {
      if (( self->size + delta ) >= self->maxSize ) {
	 self->maxSize += self->delta;
	 self->ptr = Class_Rearray( self, self->ptr, stgByte, self->maxSize * self->itmSize );
      }
      self->size += delta;
   }
   else if( delta < 0 ) {
      assert( self->size >= -delta );
      if( self->size == self->maxSize - 2 * self->delta )
	 self->ptr = Class_Rearray( self, self->ptr, stgByte, (self->maxSize - self->delta ) * self->itmSize );
   }
}

void AbsArray_Push( void* _self, void* itm ) {
   AbsArray* self = Class_Cast( _self, AbsArray );

   if( self->size == self->maxSize )
      AbsArray_Resize( self, self->size + 1 );
   else
      self->size++;
   memcpy( self->ptr + self->size - 1, itm, self->itmSize );
}

void AbsArray_Set( void* _self, void* itm ) {
   AbsArray* self = Class_Cast( _self, AbsArray );
   int i_i;

   for( i_i = 0; i_i < self->size; i_i++ )
      memcpy( self->ptr + i_i, itm, self->itmSize );
}

void AbsArray_Clear( void* _self ) {
   AbsArray* self = Class_Cast( _self, AbsArray );

   AbsArray_Resize( self, 0 );
}

int AbsArray_GetSize( const void* _self ) {
   AbsArray* self = Class_Cast( _self, AbsArray );

   return self->size;
}

const int* AbsArray_GetPtr( const void* _self ) {
   AbsArray* self = Class_Cast( _self, AbsArray );

   return (const int*)(self->ptr);
}


