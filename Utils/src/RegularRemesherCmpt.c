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

RegularRemesherCmpt* RegularRemesherCmpt_New( Name name, AbstractContext* context, Mesh* mesh, RegularRemesher* regRemesh ) {
   RegularRemesherCmpt*	self = _RegularRemesherCmpt_DefaultNew( name );

	self->isConstructed = True;
	_Remesher_Init( self, context, mesh );
   _RegularRemesherCmpt_Init( self, regRemesh );

	return self;
}

RegularRemesherCmpt* _RegularRemesherCmpt_DefaultNew( Name name ) {
   return _RegularRemesherCmpt_New(
		sizeof(RegularRemesherCmpt),
		RegularRemesherCmpt_Type,
		_RegularRemesherCmpt_Delete,
		_RegularRemesherCmpt_Print,
		NULL,
		(void*(*)(Name))_RegularRemesherCmpt_DefaultNew,
		_RegularRemesherCmpt_AssignFromXML,
		_RegularRemesherCmpt_Build,
		_RegularRemesherCmpt_Initialise,
		_RegularRemesherCmpt_Execute,
		_RegularRemesherCmpt_Destroy,
      name,
		NON_GLOBAL,
		NULL );
}

RegularRemesherCmpt* _RegularRemesherCmpt_New( REGULARREMESHERCMPT_DEFARGS ) {
   RegularRemesherCmpt*	self;

   /* Allocate memory. */
   self = (RegularRemesherCmpt*)_Remesher_New( REMESHER_PASSARGS );

   /* RegularRemesherCmpt info */

   return self;
}

void _RegularRemesherCmpt_Init( void* remesher, RegularRemesher* regRemesh ) {
   RegularRemesherCmpt*	self = (RegularRemesherCmpt*)remesher;

   /* RegularRemesherCmpt info */
   self->regRemesh = regRemesh;
}

/*----------------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _RegularRemesherCmpt_Delete( void* remesher ) {
   RegularRemesherCmpt*	self = (RegularRemesherCmpt*)remesher;

   /* Delete parent */
   _Stg_Component_Delete( self );
}

void _RegularRemesherCmpt_Print( void* remesher, Stream* stream ) {
   RegularRemesherCmpt*	self = (RegularRemesherCmpt*)remesher;
   Stream*					myStream;
	
   /* Set the Journal for printing informations */
   myStream = Journal_Register( InfoStream_Type, "RegularRemesherCmptStream" );

   /* Print parent */
   _Stg_Component_Print( self, stream );

   /* General info */
   Journal_Printf( myStream, "RegularRemesherCmpt (ptr): (%p)\n", self );

   /* Virtual info */

   /* RegularRemesherCmpt info */
}

void _RegularRemesherCmpt_AssignFromXML( void* remesher, Stg_ComponentFactory* cf, void* data ) {
   RegularRemesherCmpt*		self = (RegularRemesherCmpt*)remesher;
   Dictionary*					dict;
   Dictionary_Entry_Value*	list;
   Mesh*							mesh;
   int							nItms, dim, wall;
   int							i_i;
	RegularRemesher*			regRemesh;

   assert( self );
   assert( cf );

	regRemesh = RegularRemesher_New();
   mesh = Stg_ComponentFactory_ConstructByKey( cf, self->name, "mesh", Mesh, True, data );
   NewRemesher_SetMesh( regRemesh, mesh );

   regRemesh->contactDepth = Stg_ComponentFactory_GetInt( cf, self->name, "contactDepth", 0 );
   regRemesh->contactSize = Stg_ComponentFactory_GetDouble( cf, self->name, "contactSize", 0.0 );
   regRemesh->diffuseCorners = Stg_ComponentFactory_GetBool( cf, self->name, "diffuseCorners", False );
   regRemesh->diffuseSurface = Stg_ComponentFactory_GetBool( cf, self->name, "diffuseSurface", False );
   regRemesh->diffusionCoef = Stg_ComponentFactory_GetDouble( cf, self->name, "diffusionCoef", 1.0 );
   regRemesh->ctx = Stg_ComponentFactory_ConstructByName( cf, "context", AbstractContext, True, data );

   dict = Dictionary_Entry_Value_AsDictionary( Dictionary_Get( cf->componentDict, self->name ) );
   list = Dictionary_Get( dict, "remeshDims" );
   if( list ) {
		nItms = Dictionary_Entry_Value_GetCount( list );
		for( i_i = 0; i_i < nItms; i_i++ ) {
			dim = Dictionary_Entry_Value_AsInt( 
			Dictionary_Entry_Value_GetElement( list, i_i ) );
			RegularRemesher_SetRemeshState( regRemesh, dim, True );
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
			RegularRemesher_SetStaticWall( regRemesh, dim, wall, True );
      }
   }
   _RegularRemesherCmpt_Init( self, regRemesh );
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

   /* Delete the class itself */
   NewClass_Delete( self->regRemesh );
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Public Functions
*/


/*----------------------------------------------------------------------------------------------------------------------------------
** Private Functions
*/

#endif
