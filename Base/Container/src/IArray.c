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


void _IArray_Construct( void* _self ) {
   IArray* self = (IArray*)_self;

   _NewClass_Construct( self );
   self->size = 0;
   self->ptr = NULL;
}

void _IArray_Destruct( void* _self ) {
   IArray* self = (IArray*)_self;

   IArray_Clear( self );
   _NewClass_Destruct( self );
}

void _IArray_Copy( void* _self, const void* _op ) {
   IArray* self = (IArray*)_self;
   const IArray* op = (const IArray*)_op;

   assert( self );
   self->size = op->size;
   self->ptr = Class_Array( self, int, self->size );
   memcpy( self->ptr, op->ptr, sizeof(int) * self->size );
}

void IArray_Set( void* _self, int nItms, const int* itms ) {
   IArray* self = (IArray*)_self;

   assert( self && (!nItms || itms) );
   self->size = nItms;
   self->ptr = Class_Rearray( self, self->ptr, int, nItms );
   memcpy( self->ptr, itms, nItms * sizeof(int) );
}

void IArray_Add( void* _self, int nItms, const int* itms ) {
   IArray* self = (IArray*)_self;

   assert( self && (!nItms || itms) );
   self->ptr = Class_Rearray( self, self->ptr, int, self->size + nItms );
   memcpy( self->ptr + self->size, itms, nItms * sizeof(int) );
   self->size += nItms;
}

void IArray_Remove( void* _self, int nItms, const int* locals, IMap* map ) {
   IArray* self = (IArray*)_self;
   ISet toRemObj, *toRem = &toRemObj;
   int* ord;
   int pos;
   int i_i;

   assert( self );

   ISet_Init( toRem );
   ISet_UseArray( toRem, nItms, locals );
   ord = Class_Array( self, int, ISet_GetNumItems( toRem ) );
   ISet_GetArray( toRem, NULL, ord );
   IMap_Clear( map );
   IMap_SetMaxItems( map, nItms );
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
   if( IMap_GetNumItems( map ) < nItms )
      IMap_SetMaxItems( map, IMap_GetNumItems( map ) );
   ISet_Destruct( toRem );
   Class_Free( self, ord );

   self->size -= nItms;
   self->ptr = Class_Rearray( self, self->ptr, int, self->size );
}

void IArray_Clear( void* _self ) {
   IArray* self = (IArray*)_self;

   assert( self );
   Class_Free( self, self->ptr );
   self->size = 0;
   self->ptr = NULL;
}

int IArray_GetSize( const void* self ) {
   assert( self );
   return ((IArray*)self)->size;
}

const int* IArray_GetPtr( const void* self ) {
   assert( self );
   return ((IArray*)self)->ptr;
}

