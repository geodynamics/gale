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
** $Id: Mesh_CentroidType.c 3584 2006-05-16 11:11:07Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <mpi.h>

#include "Base/Base.h"
#include "Discretisation/Geometry/Geometry.h"

#include "types.h"
#include "shortcuts.h"
#include "Decomp.h"
#include "Sync.h"
#include "MeshTopology.h"
#include "Mesh_ElementType.h"
#include "MeshClass.h"
#include "Mesh_CentroidType.h"


/* Textual name of this class */
const Type Mesh_CentroidType_Type = "Mesh_CentroidType";


/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/

Mesh_CentroidType* Mesh_CentroidType_New( Name name ) {
	return _Mesh_CentroidType_New( sizeof(Mesh_CentroidType), 
				       Mesh_CentroidType_Type, 
				       _Mesh_CentroidType_Delete, 
				       _Mesh_CentroidType_Print, 
				       NULL, 
				       Mesh_CentroidType_Update, 
				       Mesh_CentroidType_ElementHasPoint, 
				       Mesh_CentroidType_GetMinimumSeparation, 
				       Mesh_CentroidType_GetCentroid );
}

Mesh_CentroidType* _Mesh_CentroidType_New( MESH_CENTROIDTYPE_DEFARGS ) {
	Mesh_CentroidType* self;

	/* Allocate memory */
	assert( sizeOfSelf >= sizeof(Mesh_CentroidType) );
	self = (Mesh_CentroidType*)_Mesh_ElementType_New( MESH_ELEMENTTYPE_PASSARGS );

	/* Virtual info */

	/* Mesh_CentroidType info */
	_Mesh_CentroidType_Init( self );

	return self;
}

void _Mesh_CentroidType_Init( Mesh_CentroidType* self ) {
	assert( self && Stg_CheckType( self, Mesh_CentroidType ) );

	self->elMesh = NULL;
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _Mesh_CentroidType_Delete( void* centroidType ) {
	Mesh_CentroidType*	self = (Mesh_CentroidType*)centroidType;

	/* Delete the parent. */
	_Mesh_ElementType_Delete( self );
}

void _Mesh_CentroidType_Print( void* centroidType, Stream* stream ) {
	Mesh_CentroidType*	self = (Mesh_CentroidType*)centroidType;
	Stream*			centroidTypeStream;

	centroidTypeStream = Journal_Register( InfoStream_Type, "Mesh_CentroidTypeStream" );

	/* Print parent */
	Journal_Printf( stream, "Mesh_CentroidType (ptr): (%p)\n", self );
	_Mesh_ElementType_Print( self, stream );
}

void Mesh_CentroidType_Update( void* centroidType ) {
	Mesh_CentroidType*	self = (Mesh_CentroidType*)centroidType;

	assert( self && Stg_CheckType( self, Mesh_CentroidType ) );
}

Bool Mesh_CentroidType_ElementHasPoint( void* centroidType, unsigned elInd, double* point, 
					MeshTopology_Dim* dim, unsigned* ind )
{
	Mesh_CentroidType*	self = (Mesh_CentroidType*)centroidType;

	assert( self && Stg_CheckType( self, Mesh_CentroidType ) );

	return Mesh_ElementHasPoint( self->elMesh, elInd, point, dim, ind );
}

double Mesh_CentroidType_GetMinimumSeparation( void* centroidType, unsigned elInd, double* perDim ) {
	Mesh_CentroidType*	self = (Mesh_CentroidType*)centroidType;
	Mesh_ElementType*	elType;

	assert( self && Stg_CheckType( self, Mesh_CentroidType ) );

	elType = Mesh_GetElementType( self->elMesh, elInd );

	return Mesh_ElementType_GetMinimumSeparation( elType, elInd, perDim );
}

void Mesh_CentroidType_GetCentroid( void* centroidType, unsigned element, double* centroid ) {
	Mesh_CentroidType*	self = (Mesh_CentroidType*)centroidType;
	unsigned		nInc, *inc;

	assert( self && Stg_CheckType( self, Mesh_CentroidType ) );

	Mesh_GetIncidence( self->mesh, Mesh_GetDimSize( self->mesh ), element, MT_VERTEX, &nInc, &inc );
	assert( nInc == 1 );
	memcpy( centroid, Mesh_GetVertex( self->mesh, inc[0] ), Mesh_GetDimSize( self->mesh ) * sizeof(unsigned) );
}


/*--------------------------------------------------------------------------------------------------------------------------
** Public Functions
*/

void Mesh_CentroidType_SetElementMesh( void* centroidType, void* mesh ) {
	Mesh_CentroidType*	self = (Mesh_CentroidType*)centroidType;

	assert( self && Stg_CheckType( self, Mesh_CentroidType ) );
	assert( !mesh || Stg_CheckType( mesh, Mesh ) );

	self->elMesh = mesh;
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Private Functions
*/
