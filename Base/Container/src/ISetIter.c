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
** $Id: ISetIter.c 3952 2007-01-09 06:24:06Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "StGermain/Base/Foundation/Foundation.h"
#include "types.h"
#include "Iter.h"
#include "ISet.h"
#include "ISetIter.h"
#include "StGermain/Base/Foundation/ClassDef.h"


void _ISetIter_Construct( void* _self ) {
   ISetIter* self = (ISetIter*)_self;

   _Iter_Construct( self );
   self->iset = NULL;
   self->depth = 0;
   self->stack = NULL;
   self->cur = NULL;
}

void _ISetIter_Destruct( void* self ) {
   assert( self );
   Class_Free( self, ((ISetIter*)self)->stack );
   _NewClass_Destruct( self );
}

void _ISetIter_Copy( void* _self, const void* _op ) {
   ISetIter* self = (ISetIter*)_self;
   const ISetIter* op = (const ISetIter*)_op;

   _Iter_Copy( self, op );
   self->iset = op->iset;
   self->depth = op->depth;
   if( op->stack && op->iset ) {
      self->stack = Class_Array( self, ISetItem*, op->depth );
      memcpy( self->stack, op->stack, sizeof(ISetItem*) * op->depth );
   }
   else
      self->stack = NULL;
   self->cur = op->cur;
}

void _ISetIter_Next( void* _self ) {
   ISetIter* self = (ISetIter*)_self;

   assert( self );
   assert( self->stack && self->cur && self->valid );

   if( self->cur->right ) {
      self->cur = self->cur->right;
      while( self->cur->left ) {
	 self->stack[self->depth++] = self->cur;
	 self->cur = self->cur->left;
      }
   }
   else {
      if( self->depth )
	 self->cur = self->stack[--self->depth];
      else
	 self->valid = False;
   }
}

int ISetIter_GetKey( const void* self ) {
   assert( self );
   assert( ((ISetIter*)self)->valid && ((ISetIter*)self)->cur );
   return ((ISetIter*)self)->cur->key;
}
