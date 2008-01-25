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
#include "Iter.h"
#include "ISetIter.h"
#include "ISet.h"
#include "StGermain/Base/Foundation/ClassDef.h"


ISetItem* ISet_CopyItem( ISet* self, ISetItem* itm );
void ISet_BisectArray( ISet* self, int nItms, const int* itms );
void ISet_DelItem( ISet* self, ISetItem* itm );
int ISet_Flatten( ISetItem* itm, int* keys, int ind );
void ISet_UnionItem( ISetItem* itm, ISet* dst );
void ISet_IsectItem( ISetItem* itm, const ISet* set, int* nKeys, int* keys );
void ISet_SubtrItem( ISetItem* itm, const ISet* set, int* nKeys, int* keys );
int ISet_CmpInt( const void* l, const void* r );


void _ISet_Construct( void* _self ) {
   ISet* self = (ISet*)_self;

   _NewClass_Construct( self );
   self->nItms = 0;
   self->depth = 0;
   self->root = NULL;
}

void _ISet_Destruct( void* self ) {
   ISet_Clear( self );
   _NewClass_Destruct( self );
}

void _ISet_Copy( void* _self, const void* _op ) {
   ISet* self = (ISet*)_self;
   const ISet* op = (const ISet*)_op;

   assert( op );

   if( self->root )
      ISet_DelItem( self, self->root );
   self->root = ISet_CopyItem( self, op->root );
   self->nItms = op->nItms;
   self->depth = op->depth;
}

void ISet_UseArray( void* self, int nItms, const int* itms ) {
   int* ord;

   ISet_Clear( self );
   ord = Class_Array( self, int, nItms );
   memcpy( ord, itms, sizeof(int) * nItms );
   qsort( ord, nItms, sizeof(int), ISet_CmpInt );
   ISet_BisectArray( (ISet*)self, nItms, itms );
   Class_Free( self, ord );
}

void ISet_Insert( void* _self, int key ) {
   ISet* self = (ISet*)_self;
   ISetItem** dst;
   int depth;

   assert( self );

   dst = &self->root;
   depth = 0;
   while( *dst ) {
      assert( key != (*dst)->key );
      if( key < (*dst)->key )
	 dst = &(*dst)->left;
      else
	 dst = &(*dst)->right;
      depth++;
   }

   *dst = Class_Alloc( self, ISetItem );
   (*dst)->key = key;
   (*dst)->left = NULL;
   (*dst)->right = NULL;

   self->nItms++;
   if( depth > self->depth )
      self->depth = depth;
}

void ISet_TryInsert( void* _self, int key ) {
   ISet* self = (ISet*)_self;
   ISetItem** dst;
   int depth;

   assert( self );

   dst = &self->root;
   depth = 0;
   while( *dst ) {
      if( key == (*dst)->key )
	 return;
      else if( key < (*dst)->key )
	 dst = &(*dst)->left;
      else
	 dst = &(*dst)->right;
      depth++;
   }

   *dst = Class_Alloc( self, ISetItem );
   (*dst)->key = key;
   (*dst)->left = NULL;
   (*dst)->right = NULL;

   self->nItms++;
   if( depth > self->depth )
      self->depth = depth;
}

void ISet_Clear( void* _self ) {
   ISet* self = (ISet*)_self;

   assert( self );

   if( self->root ) {
      ISet_DelItem( self, self->root );
      self->root = NULL;
      self->nItms = 0;
      self->depth = 0;
   }
}

void ISet_Balance( void* self ) {
   int nItms, *itms;

   itms = Class_Array( self, int, ((ISet*)self)->nItms );
   ISet_GetArray( self, &nItms, itms );
   ISet_Clear( self );
   ISet_BisectArray( (ISet*)self, nItms, itms );
   Class_Free( self, itms );
}

int ISet_GetNumItems( const void* self ) {
   assert( self );
   return ((ISet*)self)->nItms;
}

void ISet_GetArray( const void* self, int* nItms, int* itms ) {
   assert( self && (!((ISet*)self)->nItms || itms) );
   if( nItms )
      *nItms = ((ISet*)self)->nItms;
   ISet_Flatten( ((ISet*)self)->root, itms, 0 );
}

Bool ISet_Has( const void* _self, int key ) {
   const ISet* self = (const ISet*)_self;
   ISetItem* itm;

   assert( self );

   itm = self->root;
   while( itm ) {
      if( key == itm->key )
	 return True;
      else if( key < itm->key )
	 itm = itm->left;
      else
	 itm = itm->right;
   }
   return False;
}

void ISet_Union( const void* _self, const void* _op, void* _dst ) {
   const ISet* self = (const ISet*)_self;
   const ISet* op = (const ISet*)_op;
   ISet* dst = (ISet*)_dst;
   const ISet *big, *small;

   assert( self && op && dst );

   if( self->nItms > op->nItms ) {
      big = self;
      small = op;
   }
   else {
      big = op;
      small = self;
   }

   ISet_Copy( dst, big );
   ISet_UnionItem( small->root, dst );
}

void ISet_Isect( const void* _self, const void* _op, void* _dst ) {
   const ISet* self = (const ISet*)_self;
   const ISet* op = (const ISet*)_op;
   ISet* dst = (ISet*)_dst;
   const ISet *big, *small;
   int nKeys, *keys;

   assert( self && op && dst );

   if( self->nItms > op->nItms ) {
      big = self;
      small = op;
   }
   else {
      big = op;
      small = self;
   }

   nKeys = 0;
   keys = Class_Array( self, int, small->nItms );
   ISet_IsectItem( small->root, big, &nKeys, keys );
   ISet_UseArray( dst, nKeys, keys );
   Class_Free( self, keys );
}

void ISet_Subtr( const void* _self, const void* _op, void* _dst ) {
   const ISet* self = (const ISet*)_self;
   const ISet* op = (const ISet*)_op;
   ISet* dst = (ISet*)_dst;
   int nKeys, *keys;;

   assert( self && op && dst );

   nKeys = 0;
   keys = Class_Array( self, int, self->nItms );
   ISet_SubtrItem( self->root, op, &nKeys, keys );
   ISet_UseArray( dst, nKeys, keys );
   Class_Free( self, keys );
}

void ISet_First( const void* _self, ISetIter* iter ) {
   const ISet* self = (ISet*)_self;

   assert( self && iter );

   iter->iset = (ISet*)self;
   iter->depth = 0;
   iter->stack = Class_Rearray( self, iter->stack, ISetItem*, self->depth );
   iter->cur = self->root;
   while( iter->cur->left ) {
      iter->stack[iter->depth++] = iter->cur;
      iter->cur = iter->cur->left;
   }
   iter->valid = True;
}

ISetItem* ISet_CopyItem( ISet* self, ISetItem* itm ) {
   ISetItem* newItm;

   if( !itm )
      return NULL;
   newItm = Class_Alloc( self, ISetItem );
   newItm->key = itm->key;
   newItm->left = ISet_CopyItem( self, itm->left );
   newItm->right = ISet_CopyItem( self, itm->right );
   return newItm;
}

void ISet_BisectArray( ISet* self, int nItms, const int* itms ) {
   int mid;

   assert( self && (!nItms || itms) );

   if( !nItms )
      return;
   mid = nItms / 2;
   ISet_Insert( self, itms[mid] );
   ISet_BisectArray( self, mid, itms );
   ISet_BisectArray( self, nItms - mid - 1, itms + mid + 1 );
}

void ISet_DelItem( ISet* self, ISetItem* itm ) {
   assert( itm );

   if( itm->left )
      ISet_DelItem( self, itm->left );
   if( itm->right )
      ISet_DelItem( self, itm->right );
   Class_Free( self, itm );
}

int ISet_Flatten( ISetItem* itm, int* keys, int ind ) {
   if( itm ) {
      assert( keys );
      ind = ISet_Flatten( itm->left, keys, ind );
      keys[ind++] = itm->key;
      ind = ISet_Flatten( itm->right, keys, ind );
   }
   return ind;
}

void ISet_UnionItem( ISetItem* itm, ISet* dst ) {
   if( !itm ) return;
   ISet_TryInsert( dst, itm->key );
   ISet_UnionItem( itm->left, dst );
   ISet_UnionItem( itm->right, dst );
}

void ISet_IsectItem( ISetItem* itm, const ISet* set, int* nKeys, int* keys ) {
   if( !itm ) return;
   if( ISet_Has( set, itm->key ) )
      keys[(*nKeys)++] = itm->key;
   ISet_IsectItem( itm->left, set, nKeys, keys );
   ISet_IsectItem( itm->right, set, nKeys, keys );
}

void ISet_SubtrItem( ISetItem* itm, const ISet* set, int* nKeys, int* keys ) {
   if( !itm ) return;
   if( !ISet_Has( set, itm->key ) )
      keys[(*nKeys)++] = itm->key;
   ISet_SubtrItem( itm->left, set, nKeys, keys );
   ISet_SubtrItem( itm->right, set, nKeys, keys );
}

int ISet_CmpInt( const void* l, const void* r ) {
   assert( *(int*)l != *(int*)r );
   return (*(int*)l < *(int*)r) ? -1 : 1;
}
