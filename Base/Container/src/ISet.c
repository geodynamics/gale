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
** $Id: ISet.c 3952 2007-01-09 06:24:06Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "StGermain/Base/Foundation/Foundation.h"
#include "types.h"
#include "IArray.h"
#include "Iter.h"
#include "ISetIter.h"
#include "ISet.h"
#include "StGermain/Base/Foundation/ClassDef.h"


const double ISet_TableFactor = 1.18;


void _ISet_Construct( void* _self ) {
   ISet* self = (ISet*)_self;

   _NewClass_Construct( self );
   self->maxSize = 0;
   self->curSize = 0;
   self->tblSize = 0;
   self->tbl = NULL;
   self->used = NULL;
}

void _ISet_Destruct( void* self ) {
   ISet_Clear( self );
   Class_Free( self, ((ISet*)self)->tbl );
   Class_Free( self, ((ISet*)self)->used );
   _NewClass_Destruct( self );
}

void _ISet_Copy( void* _self, const void* _op ) {
   ISet* self = (ISet*)_self;
   const ISet* op = (const ISet*)_op;
   ISetIter iterObj, *iter = &iterObj;

   assert( self && op );
   ISet_Clear( self );
   ISet_SetMaxSize( self, op->maxSize );
   ISetIter_Init( iter );
   for( ISet_First( op, iter ); Iter_IsValid( iter ); ISetIter_Next( iter ) )
      ISet_Insert( self, ISetIter_GetKey( iter ) );
   NewClass_Destruct( iter );
}

void ISet_UseArray( void* _self, int size, const int* array ) {
   ISet* self = (ISet*)_self;
   int i_i;

   ISet_Clear( self );
   ISet_SetMaxSize( self, size );
   for( i_i = 0; i_i < size; i_i++ )
      ISet_Insert( self, array[i_i] );
}

void ISet_SetMaxSize( void* _self, int maxSize ) {
   ISet* self = (ISet*)_self;
   int nOldItms, *keys;
   ISetItem* itm;
   int i_i;

   assert( self );
   nOldItms = self->curSize;
   keys = Class_Array( self, int, self->curSize );
   ISet_GetArray( self, keys );

   ISet_Clear( self );
   self->maxSize = maxSize;
   self->curSize = 0;
   self->tblSize = (int)((double)maxSize * ISet_TableFactor);
   self->tblSize += (self->tblSize % 2) ? 0 : 1;
   self->tbl = Class_Rearray( self, self->tbl, ISetItem, self->tblSize );
   for( i_i = 0; i_i < self->tblSize; i_i++ ) {
      itm = self->tbl + i_i;
      itm->key = 0;
      itm->next = NULL;
   }
   self->used = Class_Rearray( self, self->used, Bool, self->tblSize );
   memset( self->used, 0, self->tblSize* sizeof(Bool) );

   for( i_i = 0; i_i < nOldItms; i_i++ )
      ISet_Insert( self, keys[i_i] );
   Class_Free( self, keys );
}

void ISet_Insert( void* _self, int key ) {
   ISet* self = (ISet*)_self;
   ISetItem *itm, *cur;
   int ind;

   assert( self );
   ind = ISet_Hash( self, key );
   assert( ind < self->tblSize );
   itm = self->tbl + ind;
   if( !self->used[ind] ) {
      itm->key = key;
      itm->next = NULL;
      self->used[ind] = True;
   }
   else {
#ifndef NDEBUG
      cur = itm;
      do {
	 if( cur->key == key ) {
	    assert( cur->key != key );
	 }
	 cur = cur->next;
      } while( cur );
#endif
      cur = itm->next;
      itm->next = Class_Alloc( self, ISetItem );
      itm->next->key = key;
      itm->next->next = cur;
   }
   insist( ++self->curSize, <= self->maxSize );
}

Bool ISet_TryInsert( void* _self, int key ) {
   ISet* self = (ISet*)_self;
   ISetItem *itm, *cur;
   int ind;

   assert( self );
   ind = ISet_Hash( self, key );
   assert( ind < self->tblSize );
   itm = self->tbl + ind;
   if( !self->used[ind] ) {
      itm->key = key;
      itm->next = NULL;
      self->used[ind] = True;
   }
   else {
      cur = itm;
      do {
	 if( cur->key == key )
	    return False;
	 cur = cur->next;
      } while( cur );
      cur = itm->next;
      itm->next = Class_Alloc( self, ISetItem );
      itm->next->key = key;
      itm->next->next = cur;
   }
   insist( ++self->curSize, <= self->maxSize );
   return True;
}

void ISet_Remove( void* _self, int key ) {
   ISet* self = (ISet*)_self;
   ISetItem *itm, *prev, *toDel;
   int ind;

   assert( self );
   ind = ISet_Hash( self, key );
   assert( ind < self->tblSize );
   assert( self->used[ind] );
   itm = self->tbl + ind;
   if( itm->key == key ) {
      toDel = itm->next;
      if( toDel ) {
	 itm->key = toDel->key;
	 itm->next = toDel->next;
      }
      else
	 self->used[ind] = False;
   }
   else {
      prev = itm;
      toDel = itm->next;
      while( toDel ) {
	 if( toDel->key == key ) {
	    prev->next = toDel->next;
	    break;
	 }
	 prev = toDel;
	 toDel = toDel->next;
      }
      assert( toDel );
   }
   if( toDel )
      Class_Free( self, toDel );
   self->curSize--;
}

Bool ISet_TryRemove( void* _self, int key ) {
   ISet* self = (ISet*)_self;
   ISetItem *itm, *prev, *toDel;
   int ind;

   assert( self );
   ind = ISet_Hash( self, key );
   assert( ind < self->tblSize );
   if( !self->used[ind] )
      return False;
   itm = self->tbl + ind;
   if( itm->key == key ) {
      toDel = itm->next;
      if( toDel ) {
	 itm->key = toDel->key;
	 itm->next = toDel->next;
      }
      else
	 self->used[ind] = False;
   }
   else {
      prev = itm;
      toDel = itm->next;
      while( toDel ) {
	 if( toDel->key == key ) {
	    prev->next = toDel->next;
	    break;
	 }
	 prev = toDel;
	 toDel = toDel->next;
      }
      if( !toDel )
	 return False;
   }
   if( toDel )
      Class_Free( self, toDel );
   self->curSize--;
   return True;
}

void ISet_Clear( void* _self ) {
   ISet* self = (ISet*)_self;
   ISetItem *itm, *cur, *nxt;
   int i_i;

   assert( self );
   for( i_i = 0; i_i < self->tblSize; i_i++ ) {
      self->used[i_i] = False;
      itm = self->tbl + i_i;
      cur = itm->next;
      while( cur ) {
	 nxt = cur->next;
	 Class_Free( self, cur );
	 cur = nxt;
      }
      itm->next = NULL;
   }
   self->curSize = 0;
}

void ISet_Union( void* _self, const void* _op ) {
   ISet* self = (ISet*)_self;
   const ISet* op = (const ISet*)_op;
   ISetIter iterObj, *iter = &iterObj;
   IArray arrayObj, *array = &arrayObj;
   int arrSize, key;
   const int* arrPtr;
   int i_i;

   assert( self && op );
   ISetIter_Init( iter );
   IArray_Init( array );
   for( ISet_First( op, iter ); Iter_IsValid( iter ); ISetIter_Next( iter ) ) {
      key = ISetIter_GetKey( iter );
      if( !ISet_Has( self, key ) )
	 IArray_Append( array, key );
   }
   NewClass_Destruct( iter );

   arrSize = IArray_GetSize( array );
   arrPtr = IArray_GetPtr( array );
   ISet_SetMaxSize( self, self->curSize + arrSize );
   for( i_i = 0; i_i < arrSize; i_i++ )
      ISet_Insert( self, arrPtr[i_i] );
   IArray_Destruct( array );
}

void ISet_Isect( void* _self, const void* _op ) {
   ISet* self = (ISet*)_self;
   const ISet* op = (const ISet*)_op;
   ISetIter iterObj, *iter = &iterObj;
   IArray arrayObj, *array = &arrayObj;
   int arrSize, key;
   const int* arrPtr;
   int i_i;

   assert( self && op );
   ISetIter_Init( iter );
   IArray_Init( array );
   for( ISet_First( op, iter ); Iter_IsValid( iter ); ISetIter_Next( iter ) ) {
      key = ISetIter_GetKey( iter );
      if( ISet_Has( self, key ) )
	 IArray_Append( array, key );
   }
   NewClass_Destruct( iter );

   arrSize = IArray_GetSize( array );
   arrPtr = IArray_GetPtr( array );
   ISet_Clear( self );
   ISet_SetMaxSize( self, arrSize );
   for( i_i = 0; i_i < arrSize; i_i++ )
      ISet_Insert( self, arrPtr[i_i] );
   IArray_Destruct( array );
}

void ISet_Subtr( void* _self, const void* _op ) {
   ISet* self = (ISet*)_self;
   const ISet* op = (const ISet*)_op;
   ISetIter iterObj, *iter = &iterObj;
   IArray arrayObj, *array = &arrayObj;
   int arrSize, key;
   const int* arrPtr;
   int i_i;

   assert( self && op );
   ISetIter_Init( iter );
   IArray_Init( array );
   for( ISet_First( self, iter ); Iter_IsValid( iter ); ISetIter_Next( iter ) ) {
      key = ISetIter_GetKey( iter );
      if( !ISet_Has( op, key ) )
	 IArray_Append( array, key );
   }
   NewClass_Destruct( iter );

   arrSize = IArray_GetSize( array );
   arrPtr = IArray_GetPtr( array );
   ISet_Clear( self );
   ISet_SetMaxSize( self, arrSize );
   for( i_i = 0; i_i < arrSize; i_i++ )
      ISet_Insert( self, arrPtr[i_i] );
   IArray_Destruct( array );
}

int ISet_GetMaxSize( const void* self ) {
   assert( self );
   return ((ISet*)self)->maxSize;
}

int ISet_GetSize( const void* self ) {
   assert( self );
   return ((ISet*)self)->curSize;
}

void ISet_GetArray( const void* _self, int* keys ) {
   const ISet* self = (const ISet*)_self;
   ISetIter iterObj, *iter = &iterObj;
   int i_i;

   assert( self );
   ISetIter_Init( iter );
   for( i_i = 0, ISet_First( self, iter );
	Iter_IsValid( iter );
	i_i++, ISetIter_Next( iter ) )
   {
      keys[i_i] = ISetIter_GetKey( iter );
   }
   NewClass_Destruct( iter );
}

Bool ISet_Has( const void* _self, int key ) {
   const ISet* self = (const ISet*)_self;
   ISetItem* itm;
   int ind;

   assert( self );
   ind = ISet_Hash( self, key );
   assert( ind < self->tblSize );
   if( !self->used[ind] )
      return False;
   itm = self->tbl + ind;
   if( itm->key != key ) {
      while( (itm = itm->next) ) {
	 if( itm->key == key )
	    break;
      }
   }
   return itm ? True : False;
}

int ISet_Hash( const void* self, int key ) {
   return key % ((ISet*)self)->tblSize;
}

void ISet_First( const void* _self, ISetIter* iter ) {
   const ISet* self = (ISet*)_self;
   int i_i;

   assert( self && iter );
   for( i_i = 0; i_i < self->tblSize; i_i++ ) {
      if( self->used[i_i] ) {
	 iter->iset = (ISet*)self;
	 iter->tblInd = i_i;
	 iter->cur = self->tbl + i_i;
	 iter->valid = True;
	 return;
      }
   }
   iter->valid = False;
}
