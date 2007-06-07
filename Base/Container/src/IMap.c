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
** $Id: IMap.c 3952 2007-01-09 06:24:06Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "StGermain/Base/Foundation/Foundation.h"
#include "types.h"
#include "Iter.h"
#include "IMapIter.h"
#include "IMap.h"
#include "StGermain/Base/Foundation/ClassDef.h"


const double IMap_TableFactor = 1.18;


void _IMap_Init( void* _self ) {
   IMap* self = (IMap*)_self;

   _NewClass_Init( self );
   self->maxSize = 0;
   self->curSize = 0;
   self->tblSize = 0;
   self->tbl = NULL;
   self->used = NULL;
   IMap_SetMaxSize( self, 0 );
}

void _IMap_Destruct( void* self ) {
   IMap_Clear( self );
   Class_Free( self, ((IMap*)self)->tbl );
   Class_Free( self, ((IMap*)self)->used );
   _NewClass_Destruct( self );
}

void _IMap_Copy( void* _self, const void* _op ) {
   IMap* self = (IMap*)_self;
   const IMap* op = (const IMap*)_op;
   IMapIter iter;

   assert( self && op );
   IMap_Clear( self );
   IMap_SetMaxSize( self, op->maxSize );
   IMapIter_Construct( &iter );
   for( IMap_First( op, &iter ); Iter_IsValid( &iter ); IMapIter_Next( &iter ) )
      IMap_Insert( self, IMapIter_GetKey( &iter ), IMapIter_GetValue( &iter ) );
   NewClass_Destruct( &iter );
}

void IMap_SetMaxSize( void* _self, int maxSize ) {
   IMap* self = (IMap*)_self;
   int nOldItms, *keys, *vals;
   IMapIter iterObj, *iter = &iterObj;
   IMapItem* itm;
   int i_i;

   assert( self );
   nOldItms = self->curSize;
   keys = Class_Array( self, int, self->curSize );
   vals = Class_Array( self, int, self->curSize );
   IMapIter_Construct( iter );
   for( i_i = 0, IMap_First( self, iter );
	Iter_IsValid( iter );
	i_i++, IMapIter_Next( iter ) )
   {
      keys[i_i] = IMapIter_GetKey( iter );
      vals[i_i] = IMapIter_GetValue( iter );
   }
   NewClass_Destruct( iter );

   IMap_Clear( self );
   self->maxSize = maxSize;
   self->curSize = 0;
   self->tblSize = (int)((double)maxSize * IMap_TableFactor);
   self->tblSize += (self->tblSize % 2) ? 0 : 1;
   self->tbl = Class_Rearray( self, self->tbl, IMapItem, self->tblSize );
   for( i_i = 0; i_i < self->tblSize; i_i++ ) {
      itm = self->tbl + i_i;
      itm->key = 0;
      itm->val = 0;
      itm->next = NULL;
   }
   self->used = Class_Rearray( self, self->used, Bool, self->tblSize );
   memset( self->used, 0, self->tblSize* sizeof(Bool) );

   for( i_i = 0; i_i < nOldItms; i_i++ )
      IMap_Insert( self, keys[i_i], vals[i_i] );
   Class_Free( self, keys );
   Class_Free( self, vals );
}

void IMap_Insert( void* _self, int key, int val ) {
   IMap* self = (IMap*)_self;
   IMapItem *itm, *cur;
   int ind;

   assert( self );
   ind = IMap_Hash( self, key );
   assert( ind < self->tblSize );
   itm = self->tbl + ind;
   if( !self->used[ind] ) {
      itm->key = key;
      itm->val = val;
      itm->next = NULL;
      self->used[ind] = True;
   }
   else {
#ifndef NDEBUG
      cur = itm;
      do {
	 assert( cur->key != key );
	 cur = cur->next;
      } while( cur );
#endif
      cur = itm->next;
      itm->next = Class_Alloc( self, IMapItem );
      itm->next->key = key;
      itm->next->val = val;
      itm->next->next = cur;
   }
   insist( ++self->curSize, <= self->maxSize );
}

void IMap_SetValue( void* _self, int key, int val ) {
   IMap* self = (IMap*)_self;
   IMapItem *itm;
   int ind;

   assert( self );
   ind = IMap_Hash( self, key );
   assert( ind < self->tblSize );
   assert( self->used[ind] );
   itm = self->tbl + ind;
   do {
      if( itm->key == key )
	 break;
      itm = itm->next;
   } while( itm );
   assert( itm );
   itm->val = val;
}

void IMap_Remove( void* _self, int key ) {
   IMap* self = (IMap*)_self;
   IMapItem *itm, *prev, *toDel;
   int ind;

   assert( self );
   ind = IMap_Hash( self, key );
   assert( ind < self->tblSize );
   assert( self->used[ind] );
   itm = self->tbl + ind;
   if( itm->key == key ) {
      toDel = itm->next;
      if( toDel ) {
	 itm->key = toDel->key;
	 itm->val = toDel->val;
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

void IMap_Clear( void* _self ) {
   IMap* self = (IMap*)_self;
   IMapItem *itm, *cur, *nxt;
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

int IMap_GetSize( const void* self ) {
   assert( self );
   return ((IMap*)self)->curSize;
}

int IMap_Map( const void* _self, int key ) {
   const IMap* self = (const IMap*)_self;
   IMapItem* itm;
   int ind;

   assert( self );
   ind = IMap_Hash( self, key );
   assert( ind < self->tblSize );
   assert( self->used[ind] );
   itm = self->tbl + ind;
   do {
      if( itm->key == key )
	 break;
      itm = itm->next;
   } while( itm );
   assert( itm );
   return itm->val;
}

Bool IMap_TryMap( const void* _self, int key, int* val ) {
   const IMap* self = (const IMap*)_self;
   IMapItem* itm;
   int ind;

   assert( self && val );
   ind = IMap_Hash( self, key );
   assert( ind < self->tblSize );
   if( !self->used[ind] )
      return False;
   itm = self->tbl + ind;
   do {
      if( itm->key == key )
	 break;
      itm = itm->next;
   } while( itm );
   if( !itm )
      return False;
   *val = itm->val;
   return True;
}

Bool IMap_Has( const void* _self, int key ) {
   const IMap* self = (const IMap*)_self;
   IMapItem* itm;
   int ind;

   assert( self );
   ind = IMap_Hash( self, key );
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

int IMap_Hash( const void* self, int key ) {
   return key % ((IMap*)self)->tblSize;
}

void IMap_First( const void* _self, IMapIter* iter ) {
   const IMap* self = (IMap*)_self;
   int i_i;

   assert( self && iter );
   for( i_i = 0; i_i < self->tblSize; i_i++ ) {
      if( self->used[i_i] ) {
	 iter->imap = (IMap*)self;
	 iter->tblInd = i_i;
	 iter->cur = self->tbl + i_i;
	 iter->valid = True;
	 return;
      }
   }
   iter->valid = False;
}
