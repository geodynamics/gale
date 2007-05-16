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
** $Id: ElementCellLayout.c 4102 2007-05-16 01:09:00Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include "Base/Base.h"
#include "Discretisation/Geometry/Geometry.h"
#include "Discretisation/Shape/Shape.h"
#include "Discretisation/Mesh/Mesh.h"

#include "types.h"
#include "shortcuts.h"
#include "ShadowInfo.h"
#include "CellLayout.h"
#include "ElementCellLayout.h"

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "StandardParticle.h"

const Type ElementCellLayout_Type = "ElementCellLayout";

void* _ElementCellLayout_DefaultNew( Name name ){
	return (void*) _ElementCellLayout_New( 
			sizeof(ElementCellLayout),
			ElementCellLayout_Type, 
			_ElementCellLayout_Delete,
			_ElementCellLayout_Print,
			_ElementCellLayout_Copy, 
			_ElementCellLayout_DefaultNew,
			_ElementCellLayout_Construct,
			_ElementCellLayout_Build,
			_ElementCellLayout_Initialise,
			_ElementCellLayout_Execute,
			_ElementCellLayout_Destroy,
			name,
			False,
			_ElementCellLayout_CellLocalCount,
			_ElementCellLayout_CellShadowCount,
			_ElementCellLayout_PointCount,
			_ElementCellLayout_InitialisePoints, 
			_ElementCellLayout_MapElementIdToCellId,
			_ElementCellLayout_IsInCell,
			_ElementCellLayout_CellOf,
			_ElementCellLayout_GetShadowInfo,
			NULL );
}

ElementCellLayout* ElementCellLayout_New( 
		Name name,
		void* mesh )
{ 
	return _ElementCellLayout_New( sizeof(ElementCellLayout), ElementCellLayout_Type, _ElementCellLayout_Delete,
		_ElementCellLayout_Print, _ElementCellLayout_Copy, _ElementCellLayout_DefaultNew,
		_ElementCellLayout_Construct, _ElementCellLayout_Build, _ElementCellLayout_Initialise,
		_ElementCellLayout_Execute, _ElementCellLayout_Destroy, name, True,
		_ElementCellLayout_CellLocalCount, _ElementCellLayout_CellShadowCount,
		_ElementCellLayout_PointCount, _ElementCellLayout_InitialisePoints, 
		_ElementCellLayout_MapElementIdToCellId, _ElementCellLayout_IsInCell,
		_ElementCellLayout_CellOf, _ElementCellLayout_GetShadowInfo, mesh );
}

void ElementCellLayout_Init( 
		ElementCellLayout* self, 
		Name name,
		void* mesh )
{ 
	/* General info */
	self->type = ElementCellLayout_Type;
	self->_sizeOfSelf = sizeof(ElementCellLayout);
	self->_deleteSelf = False;
	
	/* Virtual info */
	self->_delete = _ElementCellLayout_Delete;
	self->_print = _ElementCellLayout_Print;
	self->_copy = _ElementCellLayout_Copy;
	self->_defaultConstructor = _ElementCellLayout_DefaultNew;
	self->_construct = _ElementCellLayout_Construct;
	self->_build = _ElementCellLayout_Build;
	self->_initialise = _ElementCellLayout_Initialise;
	self->_execute = _ElementCellLayout_Execute;
	self->_destroy = _ElementCellLayout_Destroy;
	self->_cellLocalCount = _ElementCellLayout_CellLocalCount,
	self->_cellShadowCount = _ElementCellLayout_CellShadowCount,
	self->_pointCount = _ElementCellLayout_PointCount,
	self->_initialisePoints = _ElementCellLayout_InitialisePoints,
	self->_mapElementIdToCellId = _ElementCellLayout_MapElementIdToCellId,
	self->_isInCell = _ElementCellLayout_IsInCell,
	self->_cellOf = _ElementCellLayout_CellOf,
	self->_getShadowInfo = _ElementCellLayout_GetShadowInfo,

	_Stg_Object_Init( (Stg_Object*)self, name, NON_GLOBAL );
	_Stg_Component_Init( (Stg_Component*)self );
	_CellLayout_Init( (CellLayout*)self );
	
	/* ElementCellLayout info */
	_ElementCellLayout_Init( self, mesh );
}

ElementCellLayout* _ElementCellLayout_New( 
		SizeT					_sizeOfSelf,
		Type					type,
		Stg_Class_DeleteFunction*			_delete,
		Stg_Class_PrintFunction*			_print,
		Stg_Class_CopyFunction*			_copy, 
		Stg_Component_DefaultConstructorFunction*	_defaultConstructor,
		Stg_Component_ConstructFunction*			_construct,
		Stg_Component_BuildFunction*		_build,
		Stg_Component_InitialiseFunction*		_initialise,
		Stg_Component_ExecuteFunction*		_execute,
		Stg_Component_DestroyFunction*		_destroy,
		Name							name,
		Bool							initFlag,
		CellLayout_CellCountFunction*		_cellLocalCount,
		CellLayout_CellCountFunction*		_cellShadowCount,
		CellLayout_PointCountFunction*		_pointCount,
		CellLayout_InitialisePointsFunction*	_initialisePoints,
		CellLayout_MapElementIdToCellIdFunction*	_mapElementIdToCellId,		
		CellLayout_IsInCellFunction*		_isInCell, 
		CellLayout_CellOfFunction*		_cellOf,
		CellLayout_GetShadowInfoFunction*	_getShadowInfo,		
		void*					mesh ) 
{
	ElementCellLayout* self;
	
	/* Allocate memory */
	self = (ElementCellLayout*)_CellLayout_New( 
		_sizeOfSelf, 
		type,
		_delete,
		_print,
		_copy,
		_defaultConstructor,
		_construct,
		_build,
		_initialise,
		_execute,
		_destroy,
		name,
		initFlag,
		_cellLocalCount,
		_cellShadowCount,
		_pointCount,
		_initialisePoints,
		_mapElementIdToCellId,
		_isInCell,
		_cellOf,
		_getShadowInfo );
	
	/* General info */
	
	/* Virtual info */
	
	/* ElementCellLayout info */
	if( initFlag ){
		_ElementCellLayout_Init( self, mesh );
	}
	
	return self;
}

void _ElementCellLayout_Init( ElementCellLayout* self, void* mesh ) { 
	/* General and Virtual info should already be set */
	
	/* ElementCellInfo info */
	self->mesh = (Mesh*)mesh;
	self->isConstructed = True;
}


void _ElementCellLayout_Delete( void* elementCellLayout ) {
	ElementCellLayout* self = (ElementCellLayout*)elementCellLayout;
	
	/* Stg_Class_Delete parent class */
	_CellLayout_Delete( self );
}

void _ElementCellLayout_Print( void* elementCellLayout, Stream* stream ) {
	ElementCellLayout* self = (ElementCellLayout*)elementCellLayout;

	/* Set the Journal for printing informations */
	Stream* elementCellLayoutStream = stream;
	
	/* General info */
	Journal_Printf( elementCellLayoutStream, "ElementCellLayout (ptr): %p\n", self );
	
	/* Parent class info */
	_CellLayout_Print( self, elementCellLayoutStream );
	
	/* Virtual info */
	
	/* ElementCellLayout info */
	Journal_Printf( elementCellLayoutStream, "\tmesh (ptr): %p\n", self->mesh );
}


void* _ElementCellLayout_Copy( void* elementCellLayout, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	ElementCellLayout*	self = (ElementCellLayout*)elementCellLayout;
	ElementCellLayout*	newElementCellLayout;
	PtrMap*			map = ptrMap;
	Bool			ownMap = False;
	
	if( !map ) {
		map = PtrMap_New( 10 );
		ownMap = True;
	}
	
	newElementCellLayout = _CellLayout_Copy( self, dest, deep, nameExt, ptrMap );
	
	if( deep ) {
		newElementCellLayout->mesh = (Mesh*)Stg_Class_Copy( self->mesh, NULL, deep, nameExt, map );
	}
	else {
		newElementCellLayout->mesh = self->mesh;
	}
	
	if( ownMap ) {
		Stg_Class_Delete( map );
	}
	
	return (void*)newElementCellLayout;
}

void _ElementCellLayout_Construct( void* elementCellLayout, Stg_ComponentFactory *cf, void* data ){
	ElementCellLayout* self = (ElementCellLayout*)elementCellLayout;
	Mesh*              mesh;

	mesh =  Stg_ComponentFactory_ConstructByKey(  cf,  self->name,  "Mesh", Mesh,  True, data ) ;
	
	_CellLayout_Init( (CellLayout*)self );
	_ElementCellLayout_Init(self, mesh);
}
	
void _ElementCellLayout_Build( void *elementCellLayout, void *data ){
	ElementCellLayout*	self = (ElementCellLayout*)elementCellLayout;

	Build( self->mesh, NULL, False );

	if( !Mesh_HasIncidence( self->mesh, Mesh_GetDimSize( self->mesh ), MT_VERTEX ) ) {
		Stream* elementCellLayoutStream = Journal_Register( ErrorStream_Type, self->type );
		Journal_Printf( elementCellLayoutStream, "Warning: Mesh not configured to build element node table. "
			"Activating it now.\n" );
		abort();
	}

	ElementCellLayout_BuildShadowInfo( self );
}
	
void _ElementCellLayout_Initialise( void *elementCellLayout, void *data ){
	ElementCellLayout* self = (ElementCellLayout*)elementCellLayout;
	Initialise( self->mesh, data, False );
}
	
void _ElementCellLayout_Execute( void *elementCellLayout, void *data ){
	
}

void _ElementCellLayout_Destroy( void *elementCellLayout, void *data ){
	
}

Cell_Index _ElementCellLayout_CellLocalCount( void* elementCellLayout ) {
	ElementCellLayout* self = (ElementCellLayout*)elementCellLayout;
	return Mesh_GetLocalSize( self->mesh, Mesh_GetDimSize( self->mesh ) );
}

Cell_Index _ElementCellLayout_CellShadowCount( void* elementCellLayout ) {
	ElementCellLayout* self = (ElementCellLayout*)elementCellLayout;
	return Mesh_GetRemoteSize( self->mesh, Mesh_GetDimSize( self->mesh ) );
}

Cell_PointIndex _ElementCellLayout_PointCount( void* elementCellLayout, Cell_Index cellIndex ) {
	ElementCellLayout* self = (ElementCellLayout*)elementCellLayout;
	unsigned		nInc;
	unsigned*		inc;

	Mesh_GetIncidence( self->mesh, Mesh_GetDimSize( self->mesh ), cellIndex, MT_VERTEX, 
			   &nInc, &inc );
	return nInc;
}

void _ElementCellLayout_InitialisePoints( void* elementCellLayout, Cell_Index cellIndex, Cell_PointIndex pointCount, 
					  double*** points )
{
	ElementCellLayout* self = (ElementCellLayout*)elementCellLayout;
	Cell_PointIndex point_I;
	
	/* point to the mesh's node's coordinates */
	for( point_I = 0; point_I < pointCount; point_I++ ) {
		unsigned	nInc;
		unsigned*	inc;

	   	Mesh_GetIncidence( self->mesh, Mesh_GetDimSize( self->mesh ), cellIndex, MT_VERTEX, 
				   &nInc, &inc );
		points[point_I] = &self->mesh->verts[inc[point_I]];
	}
}


Cell_Index _ElementCellLayout_MapElementIdToCellId( void* elementCellLayout, unsigned element_dI ) {
	
	#ifdef CAUTIOUS
	{
		ElementCellLayout*      self = (ElementCellLayout*)elementCellLayout;
		Stream* errorStr = Journal_Register( Error_Type, self->type );
		Journal_Firewall( element_dI < Mesh_GetDomainSize( self->mesh, Mesh_GetDimSize( self->mesh ) ), errorStr, "Error - in %s(): User asked "
			"for cell corresponding to element %d, but the mesh that this cell layout is based on only "
			"has %d elements.\n", __func__, element_dI, Mesh_GetDomainSize( self->mesh, Mesh_GetDimSize( self->mesh ) ) );
	}	
	#endif
	
	return element_dI;
}


Bool _ElementCellLayout_IsInCell( void* elementCellLayout, Cell_Index cellIndex, void* _particle ) {
	ElementCellLayout*      self     = (ElementCellLayout*)elementCellLayout;
	GlobalParticle*	        particle = (GlobalParticle*)_particle;
	unsigned		elDim, elInd;

	return Mesh_ElementHasPoint( self->mesh, cellIndex, particle->coord, &elDim, &elInd );
}

Cell_Index _ElementCellLayout_CellOf( void* elementCellLayout, void* _particle ) {
	ElementCellLayout*      self     = (ElementCellLayout*)elementCellLayout;
	GlobalParticle*	        particle = (GlobalParticle*)_particle;
	unsigned		elInd;

	if( !Mesh_SearchElements( self->mesh, particle->coord, &elInd ) )
		elInd = Mesh_GetDomainSize( self->mesh, Mesh_GetDimSize( self->mesh ) );

	return elInd;
}


ShadowInfo* _ElementCellLayout_GetShadowInfo( void* elementCellLayout ) {
	ElementCellLayout*      self = (ElementCellLayout*)elementCellLayout;

	return &self->cellShadowInfo;
}


void ElementCellLayout_BuildShadowInfo( ElementCellLayout* self ) {
	unsigned	nDims;
	Comm*		comm;
	unsigned	nIncProcs;
	unsigned*	incProcs;
	unsigned	n_i;

	nDims = Mesh_GetDimSize( self->mesh );
	comm = Mesh_GetCommTopology( self->mesh, nDims );
	Comm_GetNeighbours( comm, &nIncProcs, &incProcs );

	/* Extract neighbouring proc information. */
	self->cellShadowInfo.procNbrInfo = Memory_Alloc_Unnamed( ProcNbrInfo );
	self->cellShadowInfo.procNbrInfo->procNbrCnt = nIncProcs;
	self->cellShadowInfo.procNbrInfo->procNbrTbl = AllocArray( unsigned, nIncProcs );
	memcpy( self->cellShadowInfo.procNbrInfo->procNbrTbl, incProcs, nIncProcs * sizeof(unsigned) );

	/* Count shadow info. */
	if( nIncProcs ) {
		self->cellShadowInfo.procShadowedCnt = AllocArray( unsigned, nIncProcs );
		memset( self->cellShadowInfo.procShadowedCnt, 0, nIncProcs * sizeof(unsigned) );
		self->cellShadowInfo.procShadowCnt = AllocArray( unsigned, nIncProcs );
		memset( self->cellShadowInfo.procShadowCnt, 0, nIncProcs * sizeof(unsigned) );
	}
	for( n_i = 0; n_i < Mesh_GetSharedSize( self->mesh, nDims ); n_i++ ) {
		unsigned	nSharers;
		unsigned*	sharers;
		unsigned	s_i;

		Mesh_GetSharers( self->mesh, nDims, n_i, 
				 &nSharers, &sharers );
		for( s_i = 0; s_i < nSharers; s_i++ )
			self->cellShadowInfo.procShadowedCnt[sharers[s_i]]++;
	}
	for( n_i = 0; n_i < Mesh_GetRemoteSize( self->mesh, nDims ); n_i++ ) {
		unsigned	owner;

		owner = Mesh_GetOwner( self->mesh, nDims, n_i );
		self->cellShadowInfo.procShadowCnt[owner]++;
	}

	/* Build shadow info indices. */
	if( nIncProcs ) {
		self->cellShadowInfo.procShadowedTbl = Memory_Alloc_2DComplex_Unnamed( unsigned, nIncProcs, 
										       self->cellShadowInfo.procShadowedCnt );
		self->cellShadowInfo.procShadowTbl = Memory_Alloc_2DComplex_Unnamed( unsigned, nIncProcs, 
										     self->cellShadowInfo.procShadowCnt );
		memset( self->cellShadowInfo.procShadowedCnt, 0, nIncProcs * sizeof(unsigned) );
		memset( self->cellShadowInfo.procShadowCnt, 0, nIncProcs * sizeof(unsigned) );
	}
	for( n_i = 0; n_i < Mesh_GetSharedSize( self->mesh, nDims ); n_i++ ) {
		unsigned	local;
		unsigned	curInd;
		unsigned	nSharers;
		unsigned*	sharers;
		unsigned	s_i;

		local = Mesh_SharedToLocal( self->mesh, nDims, n_i );

		Mesh_GetSharers( self->mesh, nDims, n_i, 
				 &nSharers, &sharers );
		for( s_i = 0; s_i < nSharers; s_i++ ) {
			curInd = self->cellShadowInfo.procShadowedCnt[sharers[s_i]]++;
			self->cellShadowInfo.procShadowedTbl[sharers[s_i]][curInd] = local;
		}
	}
	for( n_i = 0; n_i < Mesh_GetRemoteSize( self->mesh, nDims ); n_i++ ) {
		unsigned	domain;
		unsigned	curInd;
		unsigned	owner;

		domain = Mesh_GetLocalSize( self->mesh, nDims ) + n_i;
		owner = Mesh_GetOwner( self->mesh, nDims, n_i );
		curInd = self->cellShadowInfo.procShadowCnt[owner]++;
		self->cellShadowInfo.procShadowTbl[owner][curInd] = domain;
	}
}
