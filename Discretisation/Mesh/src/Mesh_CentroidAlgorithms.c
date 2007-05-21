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
** $Id: Mesh_CentroidAlgorithms.c 3584 2006-05-16 11:11:07Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <mpi.h>

#include "Base/Base.h"
#include "Discretisation/Geometry/Geometry.h"

#include "Mesh.h"


/* Textual name of this class */
const Type Mesh_CentroidAlgorithms_Type = "Mesh_CentroidAlgorithms";

/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/

Mesh_CentroidAlgorithms* Mesh_CentroidAlgorithms_New( Name name ) {
	return _Mesh_CentroidAlgorithms_New( sizeof(Mesh_CentroidAlgorithms), 
					     Mesh_CentroidAlgorithms_Type, 
					     _Mesh_CentroidAlgorithms_Delete, 
					     _Mesh_CentroidAlgorithms_Print, 
					     NULL, 
					     (void* (*)(Name))_Mesh_CentroidAlgorithms_New, 
					     _Mesh_CentroidAlgorithms_Construct, 
					     _Mesh_CentroidAlgorithms_Build, 
					     _Mesh_CentroidAlgorithms_Initialise, 
					     _Mesh_CentroidAlgorithms_Execute, 
					     _Mesh_CentroidAlgorithms_Destroy, 
					     name, 
					     NON_GLOBAL, 
					     _Mesh_Algorithms_SetMesh, 
					     Mesh_CentroidAlgorithms_Update, 
					     Mesh_CentroidAlgorithms_NearestVertex, 
					     Mesh_CentroidAlgorithms_Search, 
					     Mesh_CentroidAlgorithms_SearchElements, 
					     _Mesh_Algorithms_GetMinimumSeparation, 
					     Mesh_CentroidAlgorithms_GetLocalCoordRange, 
					     Mesh_CentroidAlgorithms_GetDomainCoordRange, 
					     Mesh_CentroidAlgorithms_GetGlobalCoordRange );
}

Mesh_CentroidAlgorithms* _Mesh_CentroidAlgorithms_New( MESH_HEXALGORITHMS_DEFARGS ) {
	Mesh_CentroidAlgorithms* self;
	
	/* Allocate memory */
	assert( sizeOfSelf >= sizeof(Mesh_CentroidAlgorithms) );
	self = (Mesh_CentroidAlgorithms*)_Mesh_Algorithms_New( MESH_ALGORITHMS_PASSARGS );

	/* Virtual info */

	/* Mesh_CentroidAlgorithms info */
	_Mesh_CentroidAlgorithms_Init( self );

	return self;
}

void _Mesh_CentroidAlgorithms_Init( Mesh_CentroidAlgorithms* self ) {
	assert( self && Stg_CheckType( self, Mesh_CentroidAlgorithms ) );

	self->elMesh = NULL;
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _Mesh_CentroidAlgorithms_Delete( void* centroidAlgorithms ) {
	Mesh_CentroidAlgorithms*	self = (Mesh_CentroidAlgorithms*)centroidAlgorithms;

	/* Delete the parent. */
	_Mesh_Algorithms_Delete( self );
}

void _Mesh_CentroidAlgorithms_Print( void* centroidAlgorithms, Stream* stream ) {
	Mesh_CentroidAlgorithms*	self = (Mesh_CentroidAlgorithms*)centroidAlgorithms;
	
	/* Set the Journal for printing informations */
	Stream* centroidAlgorithmsStream;
	centroidAlgorithmsStream = Journal_Register( InfoStream_Type, "Mesh_CentroidAlgorithmsStream" );

	/* Print parent */
	Journal_Printf( stream, "Mesh_CentroidAlgorithms (ptr): (%p)\n", self );
	_Mesh_Algorithms_Print( self, stream );
}

void _Mesh_CentroidAlgorithms_Construct( void* centroidAlgorithms, Stg_ComponentFactory* cf, void* data ) {
}

void _Mesh_CentroidAlgorithms_Build( void* centroidAlgorithms, void* data ) {
}

void _Mesh_CentroidAlgorithms_Initialise( void* centroidAlgorithms, void* data ) {
}

void _Mesh_CentroidAlgorithms_Execute( void* centroidAlgorithms, void* data ) {
}

void _Mesh_CentroidAlgorithms_Destroy( void* centroidAlgorithms, void* data ) {
}

void Mesh_CentroidAlgorithms_Update( void* centroidAlgorithms ) {
}

unsigned Mesh_CentroidAlgorithms_NearestVertex( void* centroidAlgorithms, double* point ) {
	Mesh_CentroidAlgorithms*	self = (Mesh_CentroidAlgorithms*)centroidAlgorithms;
	unsigned			elInd;

	assert( self && Stg_CheckType( self, Mesh_CentroidAlgorithms ) );

	if( Mesh_SearchElements( self->elMesh, point, &elInd ) ) {
		unsigned	nInc, *inc;

		Mesh_GetIncidence( self->elMesh, Mesh_GetDimSize( self->mesh ), elInd, MT_VERTEX, 
				   &nInc, &inc );
		assert( nInc == 1 );
		return inc[0];
	}
	else
		return _Mesh_Algorithms_NearestVertex( self, point );
}

Bool Mesh_CentroidAlgorithms_Search( void* centroidAlgorithms, double* point, 
				     MeshTopology_Dim* dim, unsigned* ind )
{
	Mesh_CentroidAlgorithms*	self = (Mesh_CentroidAlgorithms*)centroidAlgorithms;

	assert( self && Stg_CheckType( self, Mesh_CentroidAlgorithms ) );

	return Mesh_Search( self->elMesh, point, dim, ind );
}

Bool Mesh_CentroidAlgorithms_SearchElements( void* centroidAlgorithms, double* point, 
					     unsigned* elInd )
{
	Mesh_CentroidAlgorithms*	self = (Mesh_CentroidAlgorithms*)centroidAlgorithms;

	assert( self && Stg_CheckType( self, Mesh_CentroidAlgorithms ) );

	return Mesh_SearchElements( self->elMesh, point, elInd );
}

void Mesh_CentroidAlgorithms_GetLocalCoordRange( void* algorithms, double* min, double* max ) {
	Mesh_CentroidAlgorithms* self = (Mesh_CentroidAlgorithms*)algorithms;

	assert( self );
	Mesh_GetLocalCoordRange( self->elMesh, min, max );
}

void Mesh_CentroidAlgorithms_GetDomainCoordRange( void* algorithms, double* min, double* max ) {
	Mesh_CentroidAlgorithms* self = (Mesh_CentroidAlgorithms*)algorithms;

	assert( self );
	Mesh_GetDomainCoordRange( self->elMesh, min, max );
}

void Mesh_CentroidAlgorithms_GetGlobalCoordRange( void* algorithms, double* min, double* max ) {
	Mesh_CentroidAlgorithms* self = (Mesh_CentroidAlgorithms*)algorithms;

	assert( self );
	Mesh_GetGlobalCoordRange( self->elMesh, min, max );
}


/*--------------------------------------------------------------------------------------------------------------------------
** Public Functions
*/

void Mesh_CentroidAlgorithms_SetElementMesh( void* centroidAlgorithms, void* mesh ) {
	Mesh_CentroidAlgorithms*	self = (Mesh_CentroidAlgorithms*)centroidAlgorithms;

	assert( self && Stg_CheckType( self, Mesh_CentroidAlgorithms ) );

	self->elMesh = mesh;
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Private Functions
*/
