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
** $Id: NewObject.c 3952 2007-01-09 06:24:06Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdlib.h>
#include <string.h>
#include "types.h"
#include "debug.h"
#include "MemoryTag.h"
#include "Memory.h"
#include "NewClass.h"
#include "NewObject.h"
#include "StGermain/Base/Foundation/ClassDef.h"


void _NewObject_Init( void* _self ) {
   NewObject* self = (NewObject*)_self;

   _NewClass_Init( self );
   self->name = NULL;
}

void _NewObject_Destruct( void* _self ) {
   NewObject* self = (NewObject*)_self;
   assert( self );

   Class_Free( self, self->name );
   _NewClass_Destruct( self );
}

void _NewObject_Copy( void* self, const void* op ) {
   _NewClass_Copy( self, op );
   NewObject_SetName( self, ((NewObject*)op)->name );
}

void NewObject_SetName( void* _self, const char* name ) {
   NewObject* self = (NewObject*)_self;
   int len;
   assert( self );

   len = name ? strlen( name ) + 1 : 0;
   self->name = Class_Rearray( self, self->name, char, len );
   if( name )
      strcpy( self->name, name );
}

const char* NewObject_GetName( void* self ) {
   assert( self );

   return ((NewObject*)self)->name;
}



