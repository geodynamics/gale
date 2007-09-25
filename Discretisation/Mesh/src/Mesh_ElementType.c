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
** $Id: Mesh_ElementType.c 3584 2006-05-16 11:11:07Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <mpi.h>

#include "Base/Base.h"
#include "Mesh.h"


/* Textual name of this class */
const Type Mesh_ElementType_Type = "Mesh_ElementType";


/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/

Mesh_ElementType* _Mesh_ElementType_New( MESH_ELEMENTTYPE_DEFARGS ) {
	Mesh_ElementType*	self;

	/* Allocate memory */
	assert( sizeOfSelf >= sizeof(Mesh_ElementType) );
	self = (Mesh_ElementType*)_Stg_Class_New( STG_CLASS_PASSARGS );

	/* Virtual info */
	self->updateFunc = updateFunc;
	self->elementHasPointFunc = elementHasPointFunc;
	self->getMinimumSeparationFunc = getMinimumSeparationFunc;
	self->getCentroidFunc = getCentroidFunc;

	/* Mesh_ElementType info */
	_Mesh_ElementType_Init( self );

	return self;
}

void _Mesh_ElementType_Init( Mesh_ElementType* self ) {
	assert( self && Stg_CheckType( self, Mesh_ElementType ) );

	self->mesh = NULL;
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _Mesh_ElementType_Delete( void* elementType ) {
	Mesh_ElementType*	self = (Mesh_ElementType*)elementType;

	/* Delete the parent. */
	_Stg_Class_Delete( self );
}

void _Mesh_ElementType_Print( void* elementType, Stream* stream ) {
	Mesh_ElementType*	self = (Mesh_ElementType*)elementType;
	Stream*			elementTypeStream;

	elementTypeStream = Journal_Register( InfoStream_Type, "Mesh_ElementTypeStream" );

	/* Print parent */
	Journal_Printf( stream, "Mesh_ElementType (ptr): (%p)\n", self );
	_Stg_Class_Print( self, stream );
}

void _Mesh_ElementType_GetCentroid( void* elementType, unsigned element, double* centroid ) {
	Mesh_ElementType*	self = (Mesh_ElementType*)elementType;
	Mesh*			mesh;
	IArray*			inc;
	unsigned		nIncVerts;
	const int		*incVerts;
	unsigned		nDims;
	double			denom;
	unsigned		d_i, v_i;

	assert( self );

	mesh = self->mesh;
	nDims = Mesh_GetDimSize( mesh );
	inc = IArray_New();
	Mesh_GetIncidence( mesh, nDims, element, MT_VERTEX, inc );
	nIncVerts = (unsigned)IArray_GetSize( inc );
	incVerts = IArray_GetPtr( inc );

	assert( nIncVerts );
	denom = 1.0 / (double)nIncVerts;

	for( d_i = 0; d_i < nDims; d_i++ ) {
		centroid[d_i] = Mesh_GetVertex( mesh, incVerts[0] )[d_i];
		for( v_i = 1; v_i < nIncVerts; v_i++ )
			centroid[d_i] += Mesh_GetVertex( mesh, incVerts[v_i] )[d_i];
		centroid[d_i] *= denom;
	}

	NewClass_Delete( inc );
}


/*--------------------------------------------------------------------------------------------------------------------------
** Public Functions
*/

void Mesh_ElementType_SetMesh( void* elementType, void* mesh ) {
	Mesh_ElementType*	self = (Mesh_ElementType*)elementType;

	assert( self && Stg_CheckType( self, Mesh_ElementType ) );
	assert( !mesh || Stg_CheckType( mesh, Mesh ) );

	self->mesh = mesh;
	Mesh_ElementType_Update( self );
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Private Functions
*/
