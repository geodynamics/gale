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
** $Id: FeMesh_ElementType.c 3584 2006-05-16 11:11:07Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <mpi.h>

#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>

#include "Discretisation.h"


/* Textual name of this class */
const Type FeMesh_ElementType_Type = "FeMesh_ElementType";


/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/

FeMesh_ElementType* FeMesh_ElementType_New( Name name ) {
	/* Variables set in this function */
	SizeT                                                    _sizeOfSelf = sizeof(FeMesh_ElementType);
	Type                                                            type = FeMesh_ElementType_Type;
	Stg_Class_DeleteFunction*                                    _delete = _FeMesh_ElementType_Delete;
	Stg_Class_PrintFunction*                                      _print = _FeMesh_ElementType_Print;
	Stg_Class_CopyFunction*                                        _copy = NULL;
	Mesh_ElementType_UpdateFunc*                              updateFunc = Mesh_HexType_Update;
	Mesh_ElementType_ElementHasPointFunc*            elementHasPointFunc = FeMesh_ElementType_ElementHasPoint;
	Mesh_ElementType_GetMinimumSeparationFunc*  getMinimumSeparationFunc = Mesh_HexType_GetMinimumSeparation;
	Mesh_ElementType_GetCentroidFunc*                    getCentroidFunc = _Mesh_ElementType_GetCentroid;

	return _FeMesh_ElementType_New(  FEMESH_ELEMENTTYPE_PASSARGS  );
}

FeMesh_ElementType* _FeMesh_ElementType_New(  FEMESH_ELEMENTTYPE_DEFARGS  ) {
	FeMesh_ElementType* self;

	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(FeMesh_ElementType) );
	self = (FeMesh_ElementType*)_Mesh_HexType_New(  MESH_HEXTYPE_PASSARGS  );

	/* Virtual info */

	/* FeMesh_ElementType info */
	_FeMesh_ElementType_Init( self );

	return self;
}

void _FeMesh_ElementType_Init( FeMesh_ElementType* self ) {
	assert( self && Stg_CheckType( self, FeMesh_ElementType ) );

	_Mesh_HexType_Init( (Mesh_HexType*)self );
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _FeMesh_ElementType_Delete( void* elementType ) {
	FeMesh_ElementType*	self = (FeMesh_ElementType*)elementType;

	/* Delete the parent. */
	_Mesh_HexType_Delete( self );
}

void _FeMesh_ElementType_Print( void* elementType, Stream* stream ) {
	FeMesh_ElementType*	self = (FeMesh_ElementType*)elementType;
	Stream*			elementTypeStream;

	elementTypeStream = Journal_Register( InfoStream_Type, (Name)"FeMesh_ElementTypeStream"  );

	/* Print parent */
	Journal_Printf( stream, "FeMesh_ElementType (ptr): (%p)\n", self );
	_Mesh_HexType_Print( self, stream );
}

Bool FeMesh_ElementType_ElementHasPoint( void* hexType, unsigned elInd, double* point, 
					 MeshTopology_Dim* dim, unsigned* ind )
{
	FeMesh_ElementType* self = (FeMesh_ElementType*)hexType;
	int nDims, ii;

	assert( self && Stg_CheckType( self, FeMesh_ElementType ) );
	assert( Mesh_GetDimSize( self->mesh ) <= 3 );

	FeMesh_CoordGlobalToLocal( self->mesh, elInd, point, self->local );
	nDims = Mesh_GetDimSize( self->mesh );
	for( ii = 0; ii < nDims; ii++ ) {
		if( self->local[ii] < -1.0 || self->local[ii] > 1.0 )
			return False;
	}

	*dim = nDims;
	*ind = elInd;
	return True;
}


/*--------------------------------------------------------------------------------------------------------------------------
** Public Functions
*/


/*----------------------------------------------------------------------------------------------------------------------------------
** Private Functions
*/



