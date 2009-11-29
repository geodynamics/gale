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
** $Id: STreeMap.c 3952 2007-01-09 06:24:06Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <StGermain/Base/Foundation/Foundation.h>

#include "types.h"
#include "STree.h"
#include "STreeMap.h"
#include "StGermain/Base/Foundation/ClassDef.h"


void _STreeMap_Init( void* _self ) {
   STreeMap* self = Class_Cast( _self, STreeMap );

   _STree_Init( self );
   self->keySize = 0;
   self->valSize = 0;
}

void _STreeMap_Destruct( void* _self ) {
   STreeMap* self = Class_Cast( _self, STreeMap );

   _STree_Destruct( self );
}

void STreeMap_SetItemSize( void* _self, int keySize, int valSize ) {
   STreeMap* self = Class_Cast( _self, STreeMap );

   STree_SetItemSize( self, keySize + valSize );
   self->keySize = keySize;
   self->valSize = valSize;
}

void STreeMap_Insert( void* _self, const void* key, const void* val ) {
   STreeMap* self = Class_Cast( _self, STreeMap );
   stgByte* tmp;

   tmp = (stgByte*)malloc( self->keySize + self->valSize );
   assert( key );
   memcpy( tmp, key, self->keySize );
   if( val )
      memcpy( tmp + self->keySize, val, self->valSize );
   else
      memset( tmp + self->keySize, 0, self->valSize );
   STree_Insert( self, tmp );
   if( tmp )
     free( tmp );
}

void STreeMap_Remove( void* _self, const void* key ) {
   STreeMap* self = Class_Cast( _self, STreeMap );

   STree_Remove( self, key );
}

void* STreeMap_Map( const void* _self, const void* key ) {
   STreeMap* self = Class_Cast( _self, STreeMap );
   const STreeNode* node;
   int res;

   assert( self->root );
   node = self->root;
   while( node && (res = self->cmp( key, node->data )) )
     node = (res < 0) ? node->left : node->right;
   assert( node );
   return node->data + self->keySize;
}

Bool STreeMap_HasKey( const void* _self, const void* key ) {
   STreeMap* self = Class_Cast( _self, STreeMap );

   return STree_Has( self, key );
}


