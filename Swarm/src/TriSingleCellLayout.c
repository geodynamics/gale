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
** $Id: TriSingleCellLayout.c 4081 2007-04-27 06:20:07Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>

#include <StgDomain/Geometry/Geometry.h>
#include <StgDomain/Shape/Shape.h>
#include <StgDomain/Mesh/Mesh.h>

#include "types.h"
#include "shortcuts.h"
#include "CellLayout.h"
#include "TriSingleCellLayout.h"

#include <stdio.h>
#include <string.h>
#include <assert.h>

const Type TriSingleCellLayout_Type = "TriSingleCellLayout";


TriSingleCellLayout* TriSingleCellLayout_New( Name name, AbstractContext* context, int dim, Dictionary* dictionary ) { 
	TriSingleCellLayout* self = _TriSingleCellLayout_DefaultNew( name );

	self->isConstructed = True;
	_CellLayout_Init( (CellLayout*)self, context );
	_TriSingleCellLayout_Init( self, dictionary, dim );

	return self;
}

TriSingleCellLayout* _TriSingleCellLayout_DefaultNew( Name name ) {
	return (TriSingleCellLayout*)_TriSingleCellLayout_New(
		sizeof(TriSingleCellLayout),
		TriSingleCellLayout_Type,
		_TriSingleCellLayout_Delete,
		_TriSingleCellLayout_Print,
		_TriSingleCellLayout_Copy,
		(Stg_Component_DefaultConstructorFunction*)_TriSingleCellLayout_DefaultNew,
		_TriSingleCellLayout_AssignFromXML,
		_TriSingleCellLayout_Build,
		_TriSingleCellLayout_Initialise,
		_TriSingleCellLayout_Execute,
		_TriSingleCellLayout_Destroy,
		name,
		NON_GLOBAL,
		_TriSingleCellLayout_CellLocalCount,
		_TriSingleCellLayout_CellShadowCount,
		_TriSingleCellLayout_PointCount,
		_TriSingleCellLayout_InitialisePoints,
		_TriSingleCellLayout_MapElementIdToCellId,
		_TriSingleCellLayout_IsInCell,
		_TriSingleCellLayout_CellOf,
		_TriSingleCellLayout_GetShadowInfo,
		0,
		NULL );
}

TriSingleCellLayout* _TriSingleCellLayout_New( TRISINGLECELLLAYOUT_DEFARGS ) {
	TriSingleCellLayout* self;
	
	/* Allocate memory */
	self = (TriSingleCellLayout*)_CellLayout_New( CELLLAYOUT_PASSARGS );
	
	/* General info */
	self->dictionary = dictionary;
	
	/* Virtual info */
	
	/* TriSingleCellLayout info */
	
	return self;
}

void _TriSingleCellLayout_Init( TriSingleCellLayout* self, Dictionary* dictionary, int dim ) { 
	/* General and Virtual info should already be set */
	
	/* SingleCellInfo info */
	self->dictionary = dictionary;
	self->dim = dim;
}

void _TriSingleCellLayout_Delete( void* triSingleCellLayout ) {
	TriSingleCellLayout* self = (TriSingleCellLayout*)triSingleCellLayout;
	
	/* Stg_Class_Delete parent class */
	_CellLayout_Delete( self );
}

void _TriSingleCellLayout_Print( void* triSingleCellLayout, Stream* stream ) {
	TriSingleCellLayout* self = (TriSingleCellLayout*)triSingleCellLayout;
	
	/* Set the Journal for printing informations */
	Stream* triSingleCellLayoutStream = stream;
	
	/* General info */
	Journal_Printf( triSingleCellLayoutStream, "TriSingleCellLayout (ptr): %p\n", self ); 
	
	/* Parent class info */
	_CellLayout_Print( self, stream );
	
	/* Virtual info */
	
	/* TriSingleCellLayout info */
	Journal_Printf( triSingleCellLayoutStream, "self->dim: %u", self->dim );
}


void* _TriSingleCellLayout_Copy( void* triSingleCellLayout, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	TriSingleCellLayout*	self = (TriSingleCellLayout*)triSingleCellLayout;
	TriSingleCellLayout*	newTriSingleCellLayout;
	
	newTriSingleCellLayout = _CellLayout_Copy( self, dest, deep, nameExt, ptrMap );
	
	newTriSingleCellLayout->dictionary = self->dictionary;
	newTriSingleCellLayout->dim = self->dim;
	
	return (void*)newTriSingleCellLayout;
}
	
void _TriSingleCellLayout_AssignFromXML( void *triSingleCellLayout, Stg_ComponentFactory *cf, void* data ){
	TriSingleCellLayout	*self = (TriSingleCellLayout*)triSingleCellLayout;
	Dimension_Index		dim = 0;

	_CellLayout_AssignFromXML( self, cf, data );

	dim = Stg_ComponentFactory_GetRootDictUnsignedInt( cf, "dim", 0 );
	assert( dim );

	_TriSingleCellLayout_Init( (TriSingleCellLayout*)self, cf->rootDict, dim );
}
	
void _TriSingleCellLayout_Build( void* triSingleCellLayout, void* data ){
	
}
	
void _TriSingleCellLayout_Initialise( void* triSingleCellLayout, void* data ){
	
}
	
void _TriSingleCellLayout_Execute( void* triSingleCellLayout, void* data ){
	
}
	
void _TriSingleCellLayout_Destroy( void* triSingleCellLayout, void* data ){
	TriSingleCellLayout	*self = (TriSingleCellLayout*)triSingleCellLayout;

	_CellLayout_Destroy( self, data );	
}

Cell_Index _TriSingleCellLayout_CellLocalCount( void* triSingleCellLayout ) {
	/* There is only one cell... */
	return 1;
}


Cell_Index _TriSingleCellLayout_CellShadowCount( void* triSingleCellLayout ) {
	/* No shadow cells */
	return 0;
}


Cell_PointIndex _TriSingleCellLayout_PointCount( void* triSingleCellLayout, Cell_Index cellIndex ) {
	TriSingleCellLayout* self = (TriSingleCellLayout*)triSingleCellLayout;
	
	switch( self->dim ) {
		case 1:
			return 2;
		case 2:
			return 3;
		case 3:
			return 4;
		default:
			assert( 0 );
	}
	return 0;
}


void _TriSingleCellLayout_InitialisePoints( 
		void*			triSingleCellLayout, 
		Cell_Index		cellIndex, 
		Cell_PointIndex		pointCount, 
		Cell_Points		points ) 
{
	TriSingleCellLayout* self = (TriSingleCellLayout*)triSingleCellLayout;
	
	switch( self->dim ) {
		case 1:
			assert( 0 );
		case 2:
			*points[0] = Memory_Alloc_Array( double, self->dim, "points[0]" );
			*points[1] = Memory_Alloc_Array( double, self->dim, "points[1]" );
			*points[2] = Memory_Alloc_Array( double, self->dim, "points[2]" );
			
			(*points[0])[0] = 0.0f;
			(*points[0])[1] = 0.0f;
			(*points[1])[0] = 1.0f;
			(*points[1])[1] = 0.0f;
			(*points[2])[0] = 0.0f;
			(*points[2])[1] = 1.0f;
			break;
		case 3:
			assert( 0 );
		default:
			assert( 0 );
	}
}


Cell_Index _TriSingleCellLayout_MapElementIdToCellId( void* cellLayout, Element_DomainIndex element_dI ) {
	
	/* Always 0: see the header comment */
	return 0;
}


Bool _TriSingleCellLayout_IsInCell( void* triSingleCellLayout, Cell_Index cellIndex, void* particle ) {
	assert( 0 );
	return 0;
}


Cell_Index _TriSingleCellLayout_CellOf( void* triSingleCellLayout, void* particle ) {
	assert( 0 );
	return 0;
}


ShadowInfo* _TriSingleCellLayout_GetShadowInfo( void* triSingleCellLayout ) {
  /*TriSingleCellLayout*      self = (TriSingleCellLayout*)triSingleCellLayout; */

	/* TODO: this should return a shadow info with at least nbr info for my processors */
	Journal_Firewall( 0, Swarm_Warning, "Error: %s not implemented yet!\n", __func__ );
	return NULL;
}
