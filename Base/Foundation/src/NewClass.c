/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street, Melbourne, 3053, Australia.
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
** $Id: Class.c 3952 2007-01-09 06:24:06Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "types.h"
#include "forwardDecl.h"
#include "MemoryTag.h"
#include "Memory.h"
#include "Class.h"
#include "NewClass.h"
#include "ClassDef.h"


void _NewClass_Construct( void* self ) {
   assert( self );
   ((NewClass*)self)->nRefsToSelf = 0;
   ((NewClass*)self)->curAllocd = 0;
}

void _NewClass_Destruct( void* self ) {
   assert( self );
   assert( ((NewClass*)self)->curAllocd == 0 );
}

void _NewClass_Copy( void* self, const void* op ) {
   assert( self && op
 );
}

void _NewClass_Print( void* self, Stream* stream ) {
}

SizeT _NewClass_CalcMem( const void* _self, struct PtrMap* ptrs ) {
   const NewClass* self = (const NewClass*)_self;

   assert( self && ptrs );
   if( PtrMap_Find( ptrs, (void*)self ) )
      return 0;
   PtrMap_Append( ptrs, (void*)self, (void*)self );
   return self->curAllocd;
}

void NewClass_Delete( void* self ) {
   if( !self ) return;
   NewClass_Destruct( self );
   MemFree( self );
}

void NewClass_AddRef( void* _self ) {
   NewClass* self = (NewClass*)_self;

   if( !self ) return;
   assert( self->nRefsToSelf >= 0 );
   self->nRefsToSelf++;
}

void NewClass_RemoveRef( void* _self ) {
   NewClass* self = (NewClass*)_self;

   if( !self ) return;
   assert( self->nRefsToSelf > 0 );
   if( !(--self->nRefsToSelf) )
      NewClass_Delete( self );
}

void* NewClass_Dup( const void* _self ) {
   const NewClass* self = (const NewClass*)_self;
   void* dup;

   assert( self );

   dup = self->newFunc();
   NewClass_Copy( dup, self );
   return dup;
}

Type NewClass_GetType( const void* _self ) {
   const NewClass* self = (const NewClass*)_self;

   assert( self );
   return self->type;
}

SizeT NewClass_GetMemUsage( const void* _self ) {
   const NewClass* self = (const NewClass*)_self;
   struct PtrMap* ptrs;
   SizeT mem;

   assert( self );
   ptrs = PtrMap_New( 10 );
   mem = NewClass_CalcMem( self, ptrs );
   Stg_Class_Delete( ptrs );
   return mem;
}
