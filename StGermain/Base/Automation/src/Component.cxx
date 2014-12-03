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
** $Id: Component.c 3952 2007-01-09 06:24:06Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdlib.h>
#include "StGermain/Base/Foundation/Foundation.h"
#include "StGermain/Base/IO/IO.h"
#include "types.h"
#include "Stg_ComponentFactory.h"
#include "Component.h"
#include "StGermain/Base/Foundation/ClassDef.h"


void _stgComponent_Init( void* _self ) {
   stgComponent* self = Class_Cast( _self, stgComponent );

   _NewObject_Init( self );
   self->isInstantiated = False;
   self->isBuilt = False;
   self->isInitialised = False;
   self->hasExecuted = False;
}

void _stgComponent_Copy( void* _self, const void* _op ) {
   stgComponent* self = Class_Cast( _self, stgComponent );
   stgComponent* op = Class_Cast( _self, stgComponent );

   _NewObject_Copy( self, op );
   self->isInstantiated = op->isInstantiated;
   self->isBuilt = op->isBuilt;
   self->isInitialised = op->isInitialised;
   self->hasExecuted = op->hasExecuted;
}

void stgComponent_Instantiate( void* _self, Stg_ComponentFactory* cf, 
				void* data, Bool force )
{
   stgComponent* self = Class_Cast( _self, stgComponent );

   if( force || !self->isInstantiated ) {
      self->isInstantiated = True;
      stgComponent_InstantiateSelf( self, cf, data );
   }
}

void stgComponent_Build( void* _self, void* data, Bool force ) {
   stgComponent* self = Class_Cast( _self, stgComponent );

   if( force || !self->isBuilt ) {
      self->isBuilt = True;
      stgComponent_BuildSelf( self, data );
   }
}

void stgComponent_Initialise( void* _self, void* data, Bool force ) {
   stgComponent* self = Class_Cast( _self, stgComponent );

   if( force || !self->isInitialised ) {
      self->isInitialised = True;
      stgComponent_InitialiseSelf( self, data );
   }
}

void stgComponent_Execute( void* _self, void* data, Bool force ) {
   stgComponent* self = Class_Cast( _self, stgComponent );

   if( force || !self->hasExecuted ) {
      self->hasExecuted = True;
      stgComponent_ExecuteSelf( self, data );
   }
}

Bool stgComponent_IsInstantiated( const void* self ) {
   assert( Class_IsSuper( self, stgComponent ) );
   return ((stgComponent*)self)->isInstantiated;
}

Bool Component_IsBuilt( const void* self ) {
   assert( Class_IsSuper( self, stgComponent ) );
   return ((stgComponent*)self)->isBuilt;
}

Bool Component_IsInitialised( const void* self ) {
   assert( Class_IsSuper( self, stgComponent ) );
   return ((stgComponent*)self)->isInitialised;
}

Bool Component_HasExecuted( void* self ) {
   assert( Class_IsSuper( self, stgComponent ) );
   return ((stgComponent*)self)->hasExecuted;
}


