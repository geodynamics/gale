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
** $Id: RegularRemesherCmpt.c 2192 2004-10-15 02:45:38Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifdef HAVE_PETSC

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
#include "NewRemesher.h"
#include "RegularRemesher.h"
#include "RegularRemesherCmpt.h"


/* Textual name of this class */
const Type RegularRemesherCmpt_Type = "RegularRemesherCmpt";


/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/

#define REMESHER_DEFARGS				\
   sizeof(RegularRemesherCmpt),				\
      RegularRemesherCmpt_Type,				\
      _RegularRemesherCmpt_Delete,			\
      _RegularRemesherCmpt_Print,			\
      NULL,						\
      (void*(*)(Name))_RegularRemesherCmpt_DefaultNew,	\
      _RegularRemesherCmpt_AssignFromXML,			\
      _RegularRemesherCmpt_Build,			\
      _RegularRemesherCmpt_Initialise,			\
      _RegularRemesherCmpt_Execute,			\
      _RegularRemesherCmpt_Destroy,			\
      name,						\
      False,						\
      NULL


RegularRemesherCmpt* RegularRemesherCmpt_New( Name name ) {
   return _RegularRemesherCmpt_New( REMESHER_DEFARGS );
}


RegularRemesherCmpt* _RegularRemesherCmpt_New( REMESHER_ARGS ) {
   RegularRemesherCmpt*	self;

   /* Allocate memory. */
   self = (RegularRemesherCmpt*)_Remesher_New( REMESHER_PASSARGS );

   /* RegularRemesherCmpt info */
   _RegularRemesherCmpt_Init( self );

   return self;
}


void RegularRemesherCmpt_Init( RegularRemesherCmpt* self ) {
   assert( 0 ); /* TODO */
#if 0
   /* General info */
   self->type = RegularRemesherCmpt_Type;
   self->_sizeOfSelf = sizeof(RegularRemesherCmpt);
   self->_deleteSelf = False;
	
   /* Virtual info */
   self->_delete = _RegularRemesherCmpt_Delete;
   self->_print = _RegularRemesherCmpt_Print;
   self->_copy = NULL;
   _Stg_Class_Init( (Stg_Class*)self );
	
   /* RegularRemesherCmpt info */
   _RegularRemesherCmpt_Init( self );
#endif
}


void _RegularRemesherCmpt_Init( RegularRemesherCmpt* self ) {
   /* RegularRemesherCmpt info */
   self->regRemesh = RegularRemesher_New();
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _RegularRemesherCmpt_Delete( void* remesher ) {
   RegularRemesherCmpt*	self = (RegularRemesherCmpt*)remesher;

   /* Delete the class itself */
   NewClass_Delete( self->regRemesh );

   /* Delete parent */
   _Stg_Component_Delete( remesher );
}


void _RegularRemesherCmpt_Print( void* remesher, Stream* stream ) {
   RegularRemesherCmpt*	self = (RegularRemesherCmpt*)remesher;
   Stream*		myStream;
	
   /* Set the Journal for printing informations */
   myStream = Journal_Register( InfoStream_Type, "RegularRemesherCmptStream" );

   /* Print parent */
   _Stg_Component_Print( self, stream );

   /* General info */
   Journal_Printf( myStream, "RegularRemesherCmpt (ptr): (%p)\n", self );

   /* Virtual info */

   /* RegularRemesherCmpt info */
}


RegularRemesherCmpt* _RegularRemesherCmpt_DefaultNew( Name name ) {
   return _RegularRemesherCmpt_New( REMESHER_DEFARGS );
}


void _RegularRemesherCmpt_AssignFromXML( void* remesher, Stg_ComponentFactory* cf, void* data ) {
   RegularRemesherCmpt* self = (RegularRemesherCmpt*)remesher;
   Dictionary* dict;
   Dictionary_Entry_Value* list;
   Mesh* mesh;
   int nItms, dim, wall;
   int i_i;

   assert( self );
   assert( cf );

   _RegularRemesherCmpt_Init( self );

   mesh = Stg_ComponentFactory_ConstructByKey( cf, self->name, "mesh", Mesh, True, data );
   NewRemesher_SetMesh( self->regRemesh, mesh );

   self->regRemesh->contactDepth = Stg_ComponentFactory_GetInt( cf, self->name, "contactDepth", 0 );
   self->regRemesh->contactSize = Stg_ComponentFactory_GetDouble( cf, self->name, "contactSize", 0.0 );
   self->regRemesh->diffuseCorners = Stg_ComponentFactory_GetBool(
      cf, self->name, "diffuseCorners", False );
   self->regRemesh->diffuseSurface = Stg_ComponentFactory_GetBool(
      cf, self->name, "diffuseSurface", False );
   self->regRemesh->diffusionCoef = Stg_ComponentFactory_GetDouble(
      cf, self->name, "diffusionCoef", 1.0 );
   self->regRemesh->ctx = Stg_ComponentFactory_ConstructByName(
      cf, "context", AbstractContext, True, data );

   dict = Dictionary_Entry_Value_AsDictionary( Dictionary_Get( cf->componentDict, self->name ) );
   list = Dictionary_Get( dict, "remeshDims" );
   if( list ) {
      nItms = Dictionary_Entry_Value_GetCount( list );
      for( i_i = 0; i_i < nItms; i_i++ ) {
	 dim = Dictionary_Entry_Value_AsInt( 
	    Dictionary_Entry_Value_GetElement( list, i_i ) );
	 RegularRemesher_SetRemeshState( self->regRemesh, dim, True );
      }
   }

   list = Dictionary_Get( dict, "staticWalls" );
   if( list ) {
      nItms = Dictionary_Entry_Value_GetCount( list );
      assert( nItms % 2 == 0 );
      for( i_i = 0; i_i < nItms; i_i += 2 ) {
	 dim = Dictionary_Entry_Value_AsInt( 
	    Dictionary_Entry_Value_GetElement( list, i_i ) );
	 wall = Dictionary_Entry_Value_AsInt( 
	    Dictionary_Entry_Value_GetElement( list, i_i + 1 ) );
	 RegularRemesher_SetStaticWall( self->regRemesh, dim, wall, True );
      }
   }
}


void _RegularRemesherCmpt_Build( void* remesher, void* data ) {
   RegularRemesherCmpt*	self = (RegularRemesherCmpt*)remesher;

   assert( self );

   RegularRemesher_Build( self->regRemesh );
}


void _RegularRemesherCmpt_Initialise( void* remesher, void* data ) {
   RegularRemesherCmpt*	self = (RegularRemesherCmpt*)remesher;

   assert( self );
}


void _RegularRemesherCmpt_Execute( void* remesher, void* data ) {
   RegularRemesherCmpt*	self = (RegularRemesherCmpt*)remesher;

   assert( self );

   RegularRemesher_Remesh( self->regRemesh );
}


void _RegularRemesherCmpt_Destroy( void* remesher, void* data ) {
   RegularRemesherCmpt*	self = (RegularRemesherCmpt*)remesher;

   assert( self );
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Public Functions
*/


/*----------------------------------------------------------------------------------------------------------------------------------
** Private Functions
*/

#endif
