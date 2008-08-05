/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
** $Id: Remesher.c 2192 2004-10-15 02:45:38Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <mpi.h>

#include <StGermain/StGermain.h>
#include <StgDomain/Geometry/Geometry.h>
#include <StgDomain/Shape/Shape.h>
#include <StgDomain/Mesh/Mesh.h>

#include "types.h"
#include "Remesher.h"


/* Textual name of this class */
const Type Remesher_Type = "Remesher";


/*
** Constructors */

Remesher* _Remesher_New( REMESHER_ARGS ) {
   Remesher*	self;

   /* Allocate memory. */
   self = (Remesher*)_Stg_Component_New( STG_COMPONENT_PASSARGS );

   /* Virtual functions. */
   self->remeshFunc = remeshFunc;

   /* Remesher info */
   _Remesher_Init( self );

   return self;
}

void _Remesher_Init( Remesher* self ) {
   /* Remesher info */
   self->mesh = NULL;
}


/*
** Virtual functions */

void _Remesher_Delete( void* remesher ) {
   Remesher*	self = (Remesher*)remesher;

   /* Delete the class itself */

   /* Delete parent */
   _Stg_Component_Delete( self );
}

void _Remesher_Print( void* remesher, Stream* stream ) {
   Remesher*	self = (Remesher*)remesher;
   Stream*		myStream;
	
   /* Set the Journal for printing informations */
   myStream = Journal_Register( InfoStream_Type, "RemesherStream" );

   /* Print parent */
   _Stg_Component_Print( self, stream );

   /* General info */
   Journal_Printf( myStream, "Remesher (ptr): (%p)\n", self );

   /* Virtual info */

   /* Remesher info */
}

Remesher* _Remesher_DefaultNew( Name name ) {
   return _Remesher_New( sizeof(Remesher), Remesher_Type, _Remesher_Delete,
                         _Remesher_Print, NULL, (void*(*)(Name))_Remesher_DefaultNew,
                         _Remesher_Construct, _Remesher_Build, _Remesher_Initialise,
                         _Remesher_Execute, _Remesher_Destroy, name, False, NULL );
}

void _Remesher_Construct( void* remesher, Stg_ComponentFactory* cf, void* data ) {
   Remesher*	self = (Remesher*)remesher;
   char*		meshName;

   assert( self );
   assert( cf );
   assert( cf->componentDict );

   self->mesh = Stg_ComponentFactory_ConstructByKey( cf, self->name, "mesh", Mesh, True, data );
}

void _Remesher_Build( void* remesher, void* data ) {
   Remesher*	self = (Remesher*)remesher;

   assert( self );
   assert( self->mesh );

   /* Build the mesh. */
   Stg_Component_Build( self->mesh, NULL, False );
}

void _Remesher_Initialise( void* remesher, void* data ) {
   Remesher*	self = (Remesher*)remesher;

   assert( self );
}

void _Remesher_Execute( void* remesher, void* data ) {
   Remesher*	self = (Remesher*)remesher;

   assert( self );
}

void _Remesher_Destroy( void* remesher, void* data ) {
   Remesher*	self = (Remesher*)remesher;

   assert( self );

   /* TODO: If delete deletes, what does destroy do? */
}


/*
** Public Functions */
