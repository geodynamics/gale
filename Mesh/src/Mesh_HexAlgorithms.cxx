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
	/* Variables set in this function */
	SizeT                                                   _sizeOfSelf = sizeof(Mesh_HexAlgorithms);
	Type                                                           type = Mesh_HexAlgorithms_Type;
	Stg_Class_DeleteFunction*                                   _delete = _Mesh_HexAlgorithms_Delete;
	Stg_Class_PrintFunction*                                     _print = _Mesh_HexAlgorithms_Print;
	Stg_Class_CopyFunction*                                       _copy = NULL;
	Stg_Component_DefaultConstructorFunction*       _defaultConstructor = (void* (*)(Name))_Mesh_HexAlgorithms_New;
	Stg_Component_ConstructFunction*                         _construct = _Mesh_HexAlgorithms_AssignFromXML;
	Stg_Component_BuildFunction*                                 _build = _Mesh_HexAlgorithms_Build;
	Stg_Component_InitialiseFunction*                       _initialise = _Mesh_HexAlgorithms_Initialise;
	Stg_Component_ExecuteFunction*                             _execute = _Mesh_HexAlgorithms_Execute;
	Stg_Component_DestroyFunction*                             _destroy = _Mesh_HexAlgorithms_Destroy;
	AllocationType                                   nameAllocationType = NON_GLOBAL;
	Mesh_Algorithms_SetMeshFunc*                            setMeshFunc = _Mesh_Algorithms_SetMesh;
	Mesh_Algorithms_UpdateFunc*                              updateFunc = _Mesh_Algorithms_Update;
	Mesh_Algorithms_NearestVertexFunc*                nearestVertexFunc = _Mesh_Algorithms_NearestVertex;
	Mesh_Algorithms_SearchFunc*                              searchFunc = _Mesh_Algorithms_Search;
	Mesh_Algorithms_SearchElementsFunc*              searchElementsFunc = _Mesh_Algorithms_SearchElements;
	Mesh_Algorithms_GetMinimumSeparationFunc*  getMinimumSeparationFunc = _Mesh_Algorithms_GetMinimumSeparation;
	Mesh_Algorithms_GetLocalCoordRangeFunc*      getLocalCoordRangeFunc = _Mesh_Algorithms_GetLocalCoordRange;
	Mesh_Algorithms_GetDomainCoordRangeFunc*    getDomainCoordRangeFunc = _Mesh_Algorithms_GetDomainCoordRange;
	Mesh_Algorithms_GetGlobalCoordRangeFunc*    getGlobalCoordRangeFunc = _Mesh_Algorithms_GetGlobalCoordRange;

   Mesh_HexAlgorithms* self = _Mesh_HexAlgorithms_New(  MESH_HEXALGORITHMS_PASSARGS  );
	/* Mesh_HexAlgorithms info */
	_Mesh_Algorithms_Init( (Mesh_Algorithms*)self, context );
	_Mesh_HexAlgorithms_Init( self );

   return self;

}

Mesh_HexAlgorithms* _Mesh_HexAlgorithms_New(  MESH_HEXALGORITHMS_DEFARGS  ) {
	Mesh_HexAlgorithms* self;
	
	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(Mesh_HexAlgorithms) );
	self = (Mesh_HexAlgorithms*)_Mesh_Algorithms_New(  MESH_ALGORITHMS_PASSARGS  );

	return self;
}

void _Mesh_HexAlgorithms_Init( void* hexAlgorithms ) {
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


