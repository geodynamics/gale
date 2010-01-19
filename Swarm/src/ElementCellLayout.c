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
** $Id: ElementCellLayout.c 4184 2007-09-25 07:54:17Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/Geometry/Geometry.h>
#include <StgDomain/Shape/Shape.h>
#include <StgDomain/Mesh/Mesh.h>

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

ElementCellLayout* ElementCellLayout_New( Name name, AbstractContext* context, void* mesh ) { 
	ElementCellLayout* self = _ElementCellLayout_DefaultNew( name );

	self->isConstructed = True;
	_CellLayout_Init( (CellLayout*)self, context );
	_ElementCellLayout_Init( self, mesh );

	return self;
}

void* _ElementCellLayout_DefaultNew( Name name ){
	/* Variables set in this function */
	SizeT                                                _sizeOfSelf = sizeof(ElementCellLayout);
	Type                                                        type = ElementCellLayout_Type;
	Stg_Class_DeleteFunction*                                _delete = _ElementCellLayout_Delete;
	Stg_Class_PrintFunction*                                  _print = _ElementCellLayout_Print;
	Stg_Class_CopyFunction*                                    _copy = _ElementCellLayout_Copy;
	Stg_Component_DefaultConstructorFunction*    _defaultConstructor = _ElementCellLayout_DefaultNew;
	Stg_Component_ConstructFunction*                      _construct = _ElementCellLayout_AssignFromXML;
	Stg_Component_BuildFunction*                              _build = _ElementCellLayout_Build;
	Stg_Component_InitialiseFunction*                    _initialise = _ElementCellLayout_Initialise;
	Stg_Component_ExecuteFunction*                          _execute = _ElementCellLayout_Execute;
	Stg_Component_DestroyFunction*                          _destroy = _ElementCellLayout_Destroy;
	AllocationType                                nameAllocationType = NON_GLOBAL;
	CellLayout_CellCountFunction*                    _cellLocalCount = _ElementCellLayout_CellLocalCount;
	CellLayout_CellCountFunction*                   _cellShadowCount = _ElementCellLayout_CellShadowCount;
	CellLayout_PointCountFunction*                       _pointCount = _ElementCellLayout_PointCount;
	CellLayout_InitialisePointsFunction*           _initialisePoints = _ElementCellLayout_InitialisePoints;
	CellLayout_MapElementIdToCellIdFunction*   _mapElementIdToCellId = _ElementCellLayout_MapElementIdToCellId;
	CellLayout_IsInCellFunction*                           _isInCell = _ElementCellLayout_IsInCell;
	CellLayout_CellOfFunction*                               _cellOf = _ElementCellLayout_CellOf;
	CellLayout_GetShadowInfoFunction*                 _getShadowInfo = _ElementCellLayout_GetShadowInfo;

	return (void*) _ElementCellLayout_New(  ELEMENTCELLLAYOUT_PASSARGS  );
}

ElementCellLayout* _ElementCellLayout_New(  ELEMENTCELLLAYOUT_DEFARGS  ) {
	ElementCellLayout* self;
	
	/* Allocate memory */
	self = (ElementCellLayout*)_CellLayout_New(  CELLLAYOUT_PASSARGS  );
	
	/* General info */
	
	/* Virtual info */
	
	/* ElementCellLayout info */
	
	return self;
}


void _ElementCellLayout_Init( ElementCellLayout* self, void* mesh ) { 
	/* General and Virtual info should already be set */
	
	/* ElementCellInfo info */
	self->mesh = (Mesh*)mesh;
	self->incArray = IArray_New();
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

void _ElementCellLayout_AssignFromXML( void* elementCellLayout, Stg_ComponentFactory *cf, void* data ){
	ElementCellLayout* self = (ElementCellLayout*)elementCellLayout;
	Mesh*              mesh;

	_CellLayout_AssignFromXML( self, cf, data );

	mesh =  Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"Mesh", Mesh, True, data  ) ;
	
	_ElementCellLayout_Init( self, mesh );
}
	
void _ElementCellLayout_Build( void *elementCellLayout, void *data ){
	ElementCellLayout*	self = (ElementCellLayout*)elementCellLayout;

	Stg_Component_Build( self->mesh, NULL, False );

	if( !Mesh_HasIncidence( self->mesh, Mesh_GetDimSize( self->mesh ), MT_VERTEX ) ) {
		Stream* elementCellLayoutStream = Journal_Register( ErrorStream_Type, (Name)self->type  );
		Journal_Printf( elementCellLayoutStream, "Warning: Mesh not configured to build element node table. "
			"Activating it now.\n" );
		abort();
	}

	ElementCellLayout_BuildShadowInfo( self );
}
	
void _ElementCellLayout_Initialise( void *elementCellLayout, void *data ){
	ElementCellLayout* self = (ElementCellLayout*)elementCellLayout;
	Stg_Component_Initialise( self->mesh, data, False );
}
	
void _ElementCellLayout_Execute( void *elementCellLayout, void *data ){
}

void _ElementCellLayout_Destroy( void *elementCellLayout, void *data ){
	ElementCellLayout* self = (ElementCellLayout*)elementCellLayout;
		
	ElementCellLayout_DestroyShadowInfo( self );
	NewClass_Delete( self->incArray );
	
	_CellLayout_Destroy( self, data );
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

	Mesh_GetIncidence( self->mesh, Mesh_GetDimSize( self->mesh ), cellIndex, MT_VERTEX, 
			   self->incArray );
	return IArray_GetSize( self->incArray );
}

void _ElementCellLayout_InitialisePoints( void* elementCellLayout, Cell_Index cellIndex, Cell_PointIndex pointCount, 
					  double*** points )
{
	ElementCellLayout* self = (ElementCellLayout*)elementCellLayout;
	Cell_PointIndex point_I;
	unsigned	nInc;
	unsigned*	inc;

	Mesh_GetIncidence( self->mesh, Mesh_GetDimSize( self->mesh ), cellIndex, MT_VERTEX, 
			   self->incArray );
	nInc = IArray_GetSize( self->incArray );
	inc = (unsigned*)IArray_GetPtr( self->incArray );
	
	/* point to the mesh's node's coordinates */
	for( point_I = 0; point_I < pointCount; point_I++ ) {
		points[point_I] = &self->mesh->verts[inc[point_I]];
	}
}


Cell_Index _ElementCellLayout_MapElementIdToCellId( void* elementCellLayout, unsigned element_dI ) {
	
	#ifdef CAUTIOUS
	{
		ElementCellLayout*      self = (ElementCellLayout*)elementCellLayout;
		Stream* errorStr = Journal_Register( Error_Type, (Name)self->type  );
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

void ElementCellLayout_DestroyShadowInfo( ElementCellLayout* self ) {
	unsigned	nIncProcs = self->cellShadowInfo.procNbrInfo->procNbrCnt;

	/* Extract neighbouring proc information. */
	Memory_Free( self->cellShadowInfo.procNbrInfo->procNbrTbl );
   if( nIncProcs ) {
      Memory_Free( self->cellShadowInfo.procShadowedCnt );
      Memory_Free( self->cellShadowInfo.procShadowCnt );
      Memory_Free( self->cellShadowInfo.procShadowedTbl );
      Memory_Free( self->cellShadowInfo.procShadowTbl );
   }
	Memory_Free( self->cellShadowInfo.procNbrInfo );
}

void ElementCellLayout_BuildShadowInfo( ElementCellLayout* self ) {
	unsigned	nDims;
	Comm*		comm;
	int	        nIncProcs;
	const int*      incProcs;
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
		int	nSharers;
		const int*	sharers;
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
		int	        nSharers;
		const int*	        sharers;
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


