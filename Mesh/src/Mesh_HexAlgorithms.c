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
** $Id: Mesh_HexAlgorithms.c 3584 2006-05-16 11:11:07Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <mpi.h>

#include <StGermain/StGermain.h>
#include <StgDomain/Geometry/Geometry.h>

#include "Mesh.h"


/* Textual name of this class */
const Type Mesh_HexAlgorithms_Type = "Mesh_HexAlgorithms";

/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/

Mesh_HexAlgorithms* Mesh_HexAlgorithms_New( Name name, AbstractContext* context) {
   Mesh_HexAlgorithms* self = _Mesh_HexAlgorithms_New( sizeof(Mesh_HexAlgorithms), 
					Mesh_HexAlgorithms_Type, 
					_Mesh_HexAlgorithms_Delete, 
					_Mesh_HexAlgorithms_Print, 
					NULL, 
					(void* (*)(Name))_Mesh_HexAlgorithms_New, 
					_Mesh_HexAlgorithms_AssignFromXML, 
					_Mesh_HexAlgorithms_Build, 
					_Mesh_HexAlgorithms_Initialise, 
					_Mesh_HexAlgorithms_Execute, 
					_Mesh_HexAlgorithms_Destroy, 
					name, 
					NON_GLOBAL, 
					_Mesh_Algorithms_SetMesh, 
					_Mesh_Algorithms_Update, 
					_Mesh_Algorithms_NearestVertex, 
					_Mesh_Algorithms_Search, 
					_Mesh_Algorithms_SearchElements, 
					_Mesh_Algorithms_GetMinimumSeparation, 
					_Mesh_Algorithms_GetLocalCoordRange, 
					_Mesh_Algorithms_GetDomainCoordRange, 
					_Mesh_Algorithms_GetGlobalCoordRange );
	/* Mesh_HexAlgorithms info */
	_Mesh_Algorithms_Init( self, context );
	_Mesh_HexAlgorithms_Init( self );

   return self;

}

Mesh_HexAlgorithms* _Mesh_HexAlgorithms_New( MESH_HEXALGORITHMS_DEFARGS ) {
	Mesh_HexAlgorithms* self;
	
	/* Allocate memory */
	assert( sizeOfSelf >= sizeof(Mesh_HexAlgorithms) );
	self = (Mesh_HexAlgorithms*)_Mesh_Algorithms_New( MESH_ALGORITHMS_PASSARGS );

	return self;
}

void _Mesh_HexAlgorithms_Init( void* hexAlgorithms ) {
	Mesh_HexAlgorithms*	self = (Mesh_HexAlgorithms*)hexAlgorithms;
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _Mesh_HexAlgorithms_Delete( void* hexAlgorithms ) {
	Mesh_HexAlgorithms*	self = (Mesh_HexAlgorithms*)hexAlgorithms;

	/* Delete the parent. */
	_Mesh_Algorithms_Delete( self );
}

void _Mesh_HexAlgorithms_Print( void* hexAlgorithms, Stream* stream ) {
	Mesh_HexAlgorithms*	self = (Mesh_HexAlgorithms*)hexAlgorithms;
	
	/* Set the Journal for printing informations */
	Stream* hexAlgorithmsStream;
	hexAlgorithmsStream = Journal_Register( InfoStream_Type, "Mesh_HexAlgorithmsStream" );

	/* Print parent */
	Journal_Printf( stream, "Mesh_HexAlgorithms (ptr): (%p)\n", self );
	_Mesh_Algorithms_Print( self, stream );
}

void _Mesh_HexAlgorithms_AssignFromXML( void* hexAlgorithms, Stg_ComponentFactory* cf, void* data ) {

   Mesh_HexAlgorithms*	self = (Mesh_HexAlgorithms*)hexAlgorithms;
   _Mesh_Algorithms_AssignFromXML( self, cf, data );
   _Mesh_HexAlgorithms_Init( self );
}

void _Mesh_HexAlgorithms_Build( void* hexAlgorithms, void* data ) {
	Mesh_HexAlgorithms*	self = (Mesh_HexAlgorithms*)hexAlgorithms;
	_Mesh_Algorithms_Build( self, data );
}

void _Mesh_HexAlgorithms_Initialise( void* hexAlgorithms, void* data ) {
	Mesh_HexAlgorithms*	self = (Mesh_HexAlgorithms*)hexAlgorithms;
	_Mesh_Algorithms_Initialise( self, data );
}

void _Mesh_HexAlgorithms_Execute( void* hexAlgorithms, void* data ) {
	Mesh_HexAlgorithms*	self = (Mesh_HexAlgorithms*)hexAlgorithms;
	_Mesh_Algorithms_Execute( self, data );
}

void _Mesh_HexAlgorithms_Destroy( void* hexAlgorithms, void* data ) {
	Mesh_HexAlgorithms*	self = (Mesh_HexAlgorithms*)hexAlgorithms;
	_Mesh_Algorithms_Destroy( self, data );
}


/*--------------------------------------------------------------------------------------------------------------------------
** Public Functions
*/


/*----------------------------------------------------------------------------------------------------------------------------------
** Private Functions
*/
