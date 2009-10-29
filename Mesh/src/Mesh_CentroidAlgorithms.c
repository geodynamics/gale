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

#include <StGermain/StGermain.h>
#include <StgDomain/Geometry/Geometry.h>

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
					     _Mesh_CentroidAlgorithms_AssignFromXML, 
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
	self->incArray = IArray_New();
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _Mesh_CentroidAlgorithms_Delete( void* centroidAlgorithms ) {
	Mesh_CentroidAlgorithms*	self = (Mesh_CentroidAlgorithms*)centroidAlgorithms;

	/* Delete the parent. */
	_Mesh_Algorithms_Delete( self );
}

void _Mesh_CentroidAlgorithms_Print( void* _centroidAlgorithms, Stream* stream ) {
	Mesh_CentroidAlgorithms*	centroid_Algorithms= (Mesh_CentroidAlgorithms*)_centroidAlgorithms;
	
	/* Set the Journal for printing informations */
	Stream* centroidAlgorithmsStream;
	centroidAlgorithmsStream = Journal_Register( InfoStream_Type, "Mesh_CentroidAlgorithmsStream" );

	/* Print parent */
	Journal_Printf( stream, "Mesh_CentroidAlgorithms (ptr): (%p)\n", centroid_Algorithms );
	_Mesh_Algorithms_Print( centroid_Algorithms, stream );
}

void _Mesh_CentroidAlgorithms_AssignFromXML( void* _centroidAlgorithms, Stg_ComponentFactory* cf, void* data ) {
    
    Mesh_CentroidAlgorithms*      centroidAlgorithms = (Mesh_CentroidAlgorithms*)_centroidAlgorithms;
    
    _Mesh_Algorithms_AssignFromXML( centroidAlgorithms, cf, data );
}

void _Mesh_CentroidAlgorithms_Build( void* _centroidAlgorithms, void* data ) {
    
    Mesh_CentroidAlgorithms*      centroidAlgorithms = (Mesh_CentroidAlgorithms*)_centroidAlgorithms;
    
    _Mesh_Algorithms_Build( centroidAlgorithms, data );
    Stg_Component_Build( centroidAlgorithms->elMesh, data, False );  
}

void _Mesh_CentroidAlgorithms_Initialise( void* _centroidAlgorithms, void* data ) {
    Mesh_CentroidAlgorithms*      centroidAlgorithms = (Mesh_CentroidAlgorithms*)_centroidAlgorithms;
    
    _Mesh_Algorithms_Initialise( centroidAlgorithms, data );
    Stg_Component_Initialise( centroidAlgorithms->elMesh, data, False );  
}

void _Mesh_CentroidAlgorithms_Execute( void* _centroidAlgorithms, void* data ) {
    
    Mesh_CentroidAlgorithms*      centroidAlgorithms = (Mesh_CentroidAlgorithms*)_centroidAlgorithms;
    
    _Mesh_Algorithms_Initialise( centroidAlgorithms, data );
    Stg_Component_Initialise( centroidAlgorithms->elMesh, data, False );  
}

void _Mesh_CentroidAlgorithms_Destroy( void* _centroidAlgorithms, void* data ) {
    
    Mesh_CentroidAlgorithms*      centroidAlgorithms = (Mesh_CentroidAlgorithms*)_centroidAlgorithms;
    
    _Mesh_Algorithms_Destroy( centroidAlgorithms, data );
    Stg_Component_Destroy( centroidAlgorithms->elMesh, data, False );  
}

void Mesh_CentroidAlgorithms_Update( void* centroidAlgorithms ) {
}

#define Vec_Sep( nDims, v0, v1 )									\
	(((v0)[0] - (v1)[0]) * ((v0)[0] - (v1)[0]) +							\
	 (((nDims) >= 2) ? (((v0)[1] - (v1)[1]) * ((v0)[1] - (v1)[1]) + 				\
			    (((nDims) == 3) ? (((v0)[2] - (v1)[2]) * ((v0)[2] - (v1)[2])) : 0)) : 0))

unsigned Mesh_CentroidAlgorithms_NearestVertex( void* centroidAlgorithms, double* point ) {
	Mesh_CentroidAlgorithms*	self = (Mesh_CentroidAlgorithms*)centroidAlgorithms;
	unsigned			elInd;
	double dist, nearDist;
	unsigned near;
	unsigned nDims;
	double* vert;
	unsigned inc_i;

	assert( self );

	if( Mesh_SearchElements( self->elMesh, point, &elInd ) ) {
		unsigned	nInc, *inc;

		nDims = Mesh_GetDimSize( self->mesh );
		Mesh_GetIncidence( self->elMesh, Mesh_GetDimSize( self->mesh ), elInd, MT_VERTEX, 
				   self->incArray );
		nInc = IArray_GetSize( self->incArray );
		inc = (unsigned*)IArray_GetPtr( self->incArray );
		near = inc[0];
		vert = Mesh_GetVertex( self->mesh, inc[0] );
		nearDist = Vec_Sep( nDims, vert, point );
		for( inc_i = 1; inc_i < nInc; inc_i++ ) {
			vert = Mesh_GetVertex( self->mesh, inc[inc_i] );
			dist = Vec_Sep( nDims, vert, point );
			if( dist < nearDist ) {
				near = inc[inc_i];
				nearDist = dist;
			}
		}
		return near;
	}
	else
		return _Mesh_Algorithms_NearestVertex( self, point );
}

Bool Mesh_CentroidAlgorithms_Search( void* centroidAlgorithms, double* point, 
				     MeshTopology_Dim* dim, unsigned* ind )
{
	Mesh_CentroidAlgorithms*	self = (Mesh_CentroidAlgorithms*)centroidAlgorithms;

	assert( self );

	return Mesh_Search( self->elMesh, point, dim, ind );
}

Bool Mesh_CentroidAlgorithms_SearchElements( void* centroidAlgorithms, double* point, 
					     unsigned* elInd )
{
	Mesh_CentroidAlgorithms*	self = (Mesh_CentroidAlgorithms*)centroidAlgorithms;

	assert( self );

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
