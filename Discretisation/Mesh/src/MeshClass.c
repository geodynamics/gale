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
** $Id: MeshClass.c 4081 2007-04-27 06:20:07Z LukeHodkinson $
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
const Type Mesh_Type = "Mesh";


/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/

Mesh* Mesh_New( Name name ) {
	return _Mesh_New( sizeof(Mesh), 
			  Mesh_Type, 
			  _Mesh_Delete, 
			  _Mesh_Print, 
			  NULL, 
			  (void* (*)(Name))_Mesh_New, 
			  _Mesh_Construct, 
			  _Mesh_Build, 
			  _Mesh_Initialise, 
			  _Mesh_Execute, 
			  _Mesh_Destroy, 
			  name, 
			  NON_GLOBAL );
}

Mesh* _Mesh_New( MESH_DEFARGS ) {
	Mesh* self;
	
	/* Allocate memory */
	assert( sizeOfSelf >= sizeof(Mesh) );
	self = (Mesh*)_Stg_Component_New( STG_COMPONENT_PASSARGS );

	/* Virtual info */

	/* Mesh info */
	_Mesh_Init( self );

	return self;
}

void _Mesh_Init( Mesh* self ) {
	self->topo = MeshTopology_New( "" );
	self->verts = NULL;
	self->vertSyncArray = NULL;

	self->vars = List_New();
	List_SetItemSize( self->vars, sizeof(MeshVariable*) );

	self->minSep = 0.0;
	self->minAxialSep = NULL;
	self->minLocalCrd = NULL;
	self->maxLocalCrd = NULL;
	self->minDomainCrd = NULL;
	self->maxDomainCrd = NULL;
	self->minGlobalCrd = NULL;
	self->maxGlobalCrd = NULL;

	self->algorithms = Mesh_Algorithms_New( "" );
	Mesh_Algorithms_SetMesh( self->algorithms, self );
	self->nElTypes = 0;
	self->elTypes = NULL;
	self->elTypeMap = NULL;

	self->topoDataSizes = UIntMap_New();
	self->topoDataInfos = NULL;
	self->topoDatas = NULL;
	self->dataSyncArrays = NULL;
	self->info = ExtensionManager_New_OfExistingObject( "mesh_info", self );

	self->generator = NULL;
	self->emReg = NULL;
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _Mesh_Delete( void* mesh ) {
	Mesh*	self = (Mesh*)mesh;

	Mesh_Destruct( self );
	KillObject( self->algorithms );
	KillObject( self->topo );
	KillObject( self->info );

	/* Delete the parent. */
	_Stg_Component_Delete( self );
}

void _Mesh_Print( void* mesh, Stream* stream ) {
	Mesh*	self = (Mesh*)mesh;
	
	/* Set the Journal for printing informations */
	Stream* meshStream;
	meshStream = Journal_Register( InfoStream_Type, "MeshStream" );

	/* Print parent */
	Journal_Printf( stream, "Mesh (ptr): (%p)\n", self );
	_Stg_Component_Print( self, stream );
}

void _Mesh_Construct( void* mesh, Stg_ComponentFactory* cf, void* data ) {
}

void _Mesh_Build( void* mesh, void* data ) {
	Mesh*		self = (Mesh*)mesh;
	unsigned	nDims;
	unsigned	d_i;

	assert( self );

	if( !self->generator )
		return;

	MeshGenerator_Generate( self->generator, self );

	nDims = Mesh_GetDimSize( self );
	if( !nDims )
		return;

	self->vertSyncArray = Decomp_Sync_Array_New();
	Decomp_Sync_Array_SetSync( self->vertSyncArray, self->topo->domains[MT_VERTEX] );
	Decomp_Sync_Array_SetMemory( self->vertSyncArray, 
				     self->verts[0], self->verts[Mesh_GetLocalSize( self, MT_VERTEX )], 
				     nDims * sizeof(double), nDims * sizeof(double), 
				     nDims * sizeof(double) );

	self->topoDataInfos = Memory_Alloc_Array( ExtensionManager*, nDims, "mesh::topoDataInfos" );
	self->topoDatas = Memory_Alloc_Array( void*, nDims, "mesh::topoDatas" );
	self->dataSyncArrays = Memory_Alloc_Array( Decomp_Sync_Array*, nDims, "mesh::dataSyncArrays" );

	for( d_i = 0; d_i < nDims; d_i++ ) {
		char		name[20];
		unsigned	size;

		if( !UIntMap_Map( self->topoDataSizes, d_i, &size ) || !size ||
		    !Mesh_GetDomainSize( self, d_i ) )
		{
			self->topoDataInfos[d_i] = NULL;
			self->topoDatas[d_i] = NULL;
			self->dataSyncArrays[d_i] = NULL;
			continue;
		}

		sprintf( name, "topoData(%d)", d_i );
		self->topoDataInfos[d_i] = ExtensionManager_New_OfStruct( name, size );
		self->topoDatas[d_i] = (void*)ExtensionManager_Malloc( self->topoDataInfos[d_i], Mesh_GetDomainSize( self, d_i ) );

		self->dataSyncArrays[d_i] = Decomp_Sync_Array_New();
		Decomp_Sync_Array_SetSync( self->dataSyncArrays[d_i], self->topo->domains[MT_VERTEX] );
		/* TODO: Why setting verts as the memory to be synced?
		Decomp_Sync_Array_SetMemory( self->dataSyncArrays[d_i], 
					     self->verts[0], self->verts[Mesh_GetLocalSize( self, d_i )], 
					     size, size, 
					     size );
		*/
	}

	/*
	** Set up the geometric information.
	*/

	self->minAxialSep = Memory_Alloc_Array( double, nDims, "Mesh::minAxialSep" );
	self->minLocalCrd = Memory_Alloc_Array( double, nDims, "Mesh::minLocalCrd" );
	self->maxLocalCrd = Memory_Alloc_Array( double, nDims, "Mesh::maxLocalCrd" );
	self->minDomainCrd = Memory_Alloc_Array( double, nDims, "Mesh::minLocalCrd" );
	self->maxDomainCrd = Memory_Alloc_Array( double, nDims, "Mesh::maxLocalCrd" );
	self->minGlobalCrd = Memory_Alloc_Array( double, nDims, "Mesh::minGlobalCrd" );
	self->maxGlobalCrd = Memory_Alloc_Array( double, nDims, "Mesh::maxGlobalCrd" );
	Mesh_DeformationUpdate( self );
}

void _Mesh_Initialise( void* mesh, void* data ) {
}

void _Mesh_Execute( void* mesh, void* data ) {
}

void _Mesh_Destroy( void* mesh, void* data ) {
}


/*--------------------------------------------------------------------------------------------------------------------------
** Public Functions
*/

void Mesh_SetExtensionManagerRegister( void* mesh, void* extMgrReg ) {
	Mesh*	self = (Mesh*)mesh;

	assert( self );

	self->emReg = extMgrReg;
	if( extMgrReg )
		ExtensionManager_Register_Add( extMgrReg, self->info );
}

void Mesh_SetTopologyDataSize( void* mesh, MeshTopology_Dim dim, unsigned size ) {
	Mesh*	self = (Mesh*)mesh;

	assert( self );

	UIntMap_Insert( self->topoDataSizes, dim, size );
}

void Mesh_SetGenerator( void* mesh, void* generator ) {
	Mesh*	self = (Mesh*)mesh;

	assert( self );

	Mesh_Destruct( self );
	self->generator = generator;
}

void Mesh_SetAlgorithms( void* mesh, void* algorithms ) {
	Mesh*	self = (Mesh*)mesh;

	assert( self && Stg_CheckType( self, Mesh ) );

	FreeObject( self->algorithms );
	if( algorithms ) {
		assert( Stg_CheckType( algorithms, Mesh_Algorithms ) );
		self->algorithms = algorithms;
	}
	else
		self->algorithms = Mesh_Algorithms_New( "" );
	Mesh_Algorithms_SetMesh( self->algorithms, self );
}

unsigned Mesh_GetDimSize( void* mesh ) {
	Mesh*	self = (Mesh*)mesh;

	assert( self );
	assert( self->topo );

	return self->topo->nDims;
}

unsigned Mesh_GetGlobalSize( void* mesh, MeshTopology_Dim dim ) {
	Mesh*	self = (Mesh*)mesh;

	assert( self );
	assert( self->topo );

	return MeshTopology_GetGlobalSize( self->topo, dim );
}

unsigned Mesh_GetLocalSize( void* mesh, MeshTopology_Dim dim ) {
	Mesh*	self = (Mesh*)mesh;

	assert( self );
	assert( self->topo );

	return MeshTopology_GetLocalSize( self->topo, dim );
}

unsigned Mesh_GetRemoteSize( void* mesh, MeshTopology_Dim dim ) {
	Mesh*	self = (Mesh*)mesh;

	assert( self );
	assert( self->topo );

	return MeshTopology_GetRemoteSize( self->topo, dim );
}

unsigned Mesh_GetDomainSize( void* mesh, MeshTopology_Dim dim ) {
	Mesh*	self = (Mesh*)mesh;

	assert( self );
	assert( self->topo );

	return MeshTopology_GetDomainSize( self->topo, dim );
}

unsigned Mesh_GetSharedSize( void* mesh, MeshTopology_Dim dim ) {
	Mesh*	self = (Mesh*)mesh;

	assert( self );
	assert( self->topo );

	return MeshTopology_GetSharedSize( self->topo, dim );
}

MeshTopology* Mesh_GetTopology( void* mesh ) {
	Mesh*	self = (Mesh*)mesh;

	assert( self );

	return self->topo;
}

Decomp_Sync* Mesh_GetSync( void* mesh, MeshTopology_Dim dim ) {
	Mesh*	self = (Mesh*)mesh;

	assert( self );

	return MeshTopology_GetSync( self->topo, dim );
}

Bool Mesh_GlobalToDomain( void* mesh, MeshTopology_Dim dim, unsigned global, unsigned* domain ) {
	Mesh*	self = (Mesh*)mesh;

	assert( self );
	assert( self->topo );

	return MeshTopology_GlobalToDomain( self->topo, dim, global, domain );
}

unsigned Mesh_DomainToGlobal( void* mesh, MeshTopology_Dim dim, unsigned domain ) {
	Mesh*	self = (Mesh*)mesh;

	assert( self );
	assert( self->topo );

	return MeshTopology_DomainToGlobal( self->topo, dim, domain );
}

Bool Mesh_DomainToShared( void* mesh, MeshTopology_Dim dim, unsigned domain, unsigned* shared ) {
	Mesh*	self = (Mesh*)mesh;

	assert( self );
	assert( self->topo );

	return MeshTopology_DomainToShared( self->topo, dim, domain, shared );
}

unsigned Mesh_SharedToDomain( void* mesh, MeshTopology_Dim dim, unsigned shared ) {
	Mesh*	self = (Mesh*)mesh;

	assert( self );
	assert( self->topo );

	return MeshTopology_SharedToDomain( self->topo, dim, shared );
}

unsigned Mesh_GetOwner( void* mesh, MeshTopology_Dim dim, unsigned remote ) {
	Mesh*	self = (Mesh*)mesh;

	assert( self );

	return MeshTopology_GetOwner( self->topo, dim, remote );
}

void Mesh_GetSharers( void* mesh, MeshTopology_Dim dim, unsigned shared, 
		      unsigned* nSharers, unsigned** sharers )
{
	Mesh*	self = (Mesh*)mesh;

	assert( self );
	assert( self->topo );

	MeshTopology_GetSharers( self->topo, dim, shared, nSharers, sharers );
}

Bool Mesh_HasIncidence( void* mesh, MeshTopology_Dim fromDim, MeshTopology_Dim toDim ) {
	Mesh*	self = (Mesh*)mesh;

	assert( self );
	assert( self->topo );

	return MeshTopology_HasIncidence( self->topo, fromDim, toDim );
}

unsigned Mesh_GetIncidenceSize( void* mesh, MeshTopology_Dim fromDim, unsigned fromInd, 
				MeshTopology_Dim toDim )
{
	Mesh*	self = (Mesh*)mesh;

	assert( self );

	return MeshTopology_GetIncidenceSize( self->topo, fromDim, fromInd, toDim );
}

void Mesh_GetIncidence( void* mesh, MeshTopology_Dim fromDim, unsigned fromInd, MeshTopology_Dim toDim, 
			unsigned* nInc, unsigned** inc )
{
	Mesh*	self = (Mesh*)mesh;

	assert( self );
	assert( self->topo );

	MeshTopology_GetIncidence( self->topo, fromDim, fromInd, toDim, nInc, inc );
}

unsigned Mesh_NearestVertex( void* mesh, double* point ) {
	Mesh*	self = (Mesh*)mesh;

	assert( self );

	return Mesh_Algorithms_NearestVertex( self->algorithms, point );
}

Bool Mesh_Search( void* mesh, double* point, 
		  MeshTopology_Dim* dim, unsigned* ind )
{
	Mesh*	self = (Mesh*)mesh;

	assert( self && Stg_CheckType( self, Mesh ) );

	return Mesh_Algorithms_Search( self->algorithms, point, dim, ind );
}

Bool Mesh_SearchElements( void* mesh, double* point, unsigned* elInd ) {
	Mesh*	self = (Mesh*)mesh;

	assert( self && Stg_CheckType( self, Mesh ) );

	return Mesh_Algorithms_SearchElements( self->algorithms, point, elInd );
}

Bool Mesh_ElementHasPoint( void* mesh, unsigned element, double* point, 
			   MeshTopology_Dim* dim, unsigned* ind )
{
	Mesh*	self = (Mesh*)mesh;

	assert( self );
	assert( element < Mesh_GetDomainSize( self, Mesh_GetDimSize( self ) ) );
	assert( self->elTypeMap );
	assert( self->elTypeMap[element] < self->nElTypes );
	assert( self->elTypes[self->elTypeMap[element]] );

	return Mesh_ElementType_ElementHasPoint( self->elTypes[self->elTypeMap[element]], element, point, 
						 dim, ind );
}

Mesh_ElementType* Mesh_GetElementType( void* mesh, unsigned element ) {
	Mesh*	self = (Mesh*)mesh;

	assert( self );
	assert( element < Mesh_GetDomainSize( self, Mesh_GetDimSize( self ) ) );
	assert( self->elTypeMap );
	assert( self->elTypeMap[element] < self->nElTypes );
	assert( self->elTypes );

	return self->elTypes[self->elTypeMap[element]];
}

CommTopology* Mesh_GetCommTopology( void* mesh, MeshTopology_Dim dim ) {
	Mesh*	self = (Mesh*)mesh;

	assert( self );

	return MeshTopology_GetCommTopology( self->topo, dim );
}

double* Mesh_GetVertex( void* mesh, unsigned domain ) {
	Mesh*	self = (Mesh*)mesh;

	assert( self );
	assert( domain < Mesh_GetDomainSize( self, MT_VERTEX ) );
	assert( self->verts );

	return self->verts[domain];
}

void* Mesh_GetTopologyData( void* mesh, MeshTopology_Dim dim ) {
	Mesh*	self = (Mesh*)mesh;

	assert( self );
	assert( self->topoDatas );
	assert( dim < Mesh_GetDimSize( self ) );

	return self->topoDatas[dim];
}

void Mesh_GetMinimumSeparation( void* mesh, double* minSep, double* axial ) {
	Mesh*	self = (Mesh*)mesh;

	assert( self );
	assert( minSep );

	*minSep = self->minSep;
	if( axial )
		memcpy( axial, self->minAxialSep, Mesh_GetDimSize( self ) * sizeof(double) );
}

void Mesh_GetLocalCoordRange( void* mesh, double* min, double* max ) {
	Mesh*	self = (Mesh*)mesh;

	assert( self );
	assert( min );
	assert( max );

	memcpy( min, self->minLocalCrd, Mesh_GetDimSize( self ) * sizeof(double) );
	memcpy( max, self->maxLocalCrd, Mesh_GetDimSize( self ) * sizeof(double) );
}

void Mesh_GetDomainCoordRange( void* mesh, double* min, double* max ) {
	Mesh*	self = (Mesh*)mesh;

	assert( self );
	assert( min );
	assert( max );

	memcpy( min, self->minDomainCrd, Mesh_GetDimSize( self ) * sizeof(double) );
	memcpy( max, self->maxDomainCrd, Mesh_GetDimSize( self ) * sizeof(double) );
}

void Mesh_GetGlobalCoordRange( void* mesh, double* min, double* max ) {
	Mesh*	self = (Mesh*)mesh;

	assert( self );

	memcpy( min, self->minGlobalCrd, Mesh_GetDimSize( self ) * sizeof(double) );
	memcpy( max, self->maxGlobalCrd, Mesh_GetDimSize( self ) * sizeof(double) );
}

void Mesh_DeformationUpdate( void* mesh ) {
	Mesh*	self = (Mesh*)mesh;

	assert( self );

	self->minSep = Mesh_Algorithms_GetMinimumSeparation( self->algorithms, self->minAxialSep );
	Mesh_Algorithms_GetLocalCoordRange( self->algorithms, self->minLocalCrd, self->maxLocalCrd );
	Mesh_Algorithms_GetDomainCoordRange( self->algorithms, self->minDomainCrd, self->maxDomainCrd );
	Mesh_Algorithms_GetGlobalCoordRange( self->algorithms, self->minGlobalCrd, self->maxGlobalCrd );

	Mesh_Algorithms_Update( self->algorithms );
}

void Mesh_Sync( void* mesh ) {
	Mesh*	self = (Mesh*)mesh;

	assert( self );

	if( self->vertSyncArray )
		Decomp_Sync_Array_Sync( self->vertSyncArray );

	if( self->dataSyncArrays ) {
		unsigned	d_i;

		for( d_i = 0; d_i < Mesh_GetDimSize( self ); d_i++ ) {
			if( self->dataSyncArrays[d_i] )
				Decomp_Sync_Array_Sync( self->dataSyncArrays[d_i] );
		}
	}
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Private Functions
*/

void Mesh_Destruct( Mesh* self ) {
	unsigned	et_i, v_i;

	for( et_i = 0; et_i < self->nElTypes; et_i++ )
		FreeObject( self->elTypes[et_i] );
	KillArray( self->elTypes );
	KillArray( self->elTypeMap );
	self->nElTypes = 0;

	KillArray( self->verts );
	KillObject( self->vertSyncArray );

	self->generator = NULL;
	self->emReg = NULL;

	for( v_i = 0; v_i < List_GetSize( self->vars ); v_i++ ) {
		MeshVariable*	var;

		var = *(MeshVariable**)List_GetItem( self->vars, v_i );
		MeshVariable_SetMesh( var, NULL );
	}
	List_Clear( self->vars );
}
