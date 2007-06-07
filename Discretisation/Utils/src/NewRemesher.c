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
** $Id: NewRemesher.c 3952 2007-01-09 06:24:06Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdlib.h>
#include <StGermain/Base/Base.h>
#include <StGermain/Discretisation/Mesh/Mesh.h>
#include "types.h"
#include "NewRemesher.h"
#include "StGermain/Base/Foundation/ClassDef.h"


void _NewRemesher_Init( void* _self ) {
   NewRemesher* self = Class_Cast( _self, NewRemesher );

   _NewClass_Init( self );
   self->mesh = NULL;
}

void _NewRemesher_Copy( void* _self, const void* _op ) {
   NewRemesher* self = Class_Cast( _self, NewRemesher );
   const NewRemesher* op = Class_Cast( _op, NewRemesher );

   _NewClass_Copy( self, op );
   self->mesh = op->mesh;
}

void _NewRemesher_Print( const void* _self, Stream* stream ) {
   NewRemesher* self = Class_Cast( _self, NewRemesher );

   _NewClass_Print( self, stream );
}

void NewRemesher_SetMesh( void* _self, void* mesh ) {
   NewRemesher* self = Class_Cast( _self, NewRemesher );

   /*assert( Class_IsSuper( mesh, Mesh ) );*/

   self->mesh = (Mesh*)mesh;
}

Mesh* NewRemesher_GetMesh( const void* self ) {
   assert( Class_IsSuper( self, NewRemesher ) );
   return ((NewRemesher*)self)->mesh;
}
