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
** $Id: SingleCellLayout.c 4081 2007-04-27 06:20:07Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>

#include <StgDomain/Geometry/Geometry.h>
#include <StgDomain/Shape/Shape.h>
#include <StgDomain/Mesh/Mesh.h>
#include <StgDomain/Utils/Utils.h>

#include "types.h"
#include "shortcuts.h"
#include "CellLayout.h"
#include "SingleCellLayout.h"

#include <stdio.h>
#include <string.h>
#include <assert.h>

const Type SingleCellLayout_Type = "SingleCellLayout";


SingleCellLayout* SingleCellLayout_New( Name name, AbstractContext* context, const Bool dimExists[3], const XYZ min, const XYZ max ) { 
	SingleCellLayout* self = (SingleCellLayout*) _SingleCellLayout_DefaultNew( name );

	self->isConstructed = True;
	_CellLayout_Init( (CellLayout*)self, context );
	_SingleCellLayout_Init( self, dimExists, min, max );

	return self;
}

void* _SingleCellLayout_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                                _sizeOfSelf = sizeof(SingleCellLayout);
	Type                                                        type = SingleCellLayout_Type;
	Stg_Class_DeleteFunction*                                _delete = _SingleCellLayout_Delete;
	Stg_Class_PrintFunction*                                  _print = _SingleCellLayout_Print;
	Stg_Class_CopyFunction*                                    _copy = _SingleCellLayout_Copy;
	Stg_Component_DefaultConstructorFunction*    _defaultConstructor = _SingleCellLayout_DefaultNew;
	Stg_Component_ConstructFunction*                      _construct = _SingleCellLayout_AssignFromXML;
	Stg_Component_BuildFunction*                              _build = _SingleCellLayout_Build;
	Stg_Component_InitialiseFunction*                    _initialise = _SingleCellLayout_Initialise;
	Stg_Component_ExecuteFunction*                          _execute = _SingleCellLayout_Execute;
	Stg_Component_DestroyFunction*                          _destroy = _SingleCellLayout_Destroy;
	AllocationType                                nameAllocationType = NON_GLOBAL;
	CellLayout_CellCountFunction*                    _cellLocalCount = _SingleCellLayout_CellLocalCount;
	CellLayout_CellCountFunction*                   _cellShadowCount = _SingleCellLayout_CellShadowCount;
	CellLayout_PointCountFunction*                       _pointCount = _SingleCellLayout_PointCount;
	CellLayout_InitialisePointsFunction*           _initialisePoints = _SingleCellLayout_InitialisePoints;
	CellLayout_MapElementIdToCellIdFunction*   _mapElementIdToCellId = _SingleCellLayout_MapElementIdToCellId;
	CellLayout_IsInCellFunction*                           _isInCell = _SingleCellLayout_IsInCell;
	CellLayout_CellOfFunction*                               _cellOf = _SingleCellLayout_CellOf;
	CellLayout_GetShadowInfoFunction*                 _getShadowInfo = _SingleCellLayout_GetShadowInfo;

	return (void*)_SingleCellLayout_New(  SINGLECELLLAYOUT_PASSARGS  );
}

SingleCellLayout* _SingleCellLayout_New(  SINGLECELLLAYOUT_DEFARGS  ) {
	SingleCellLayout* self;
	
	/* Allocate memory */
	self = (SingleCellLayout*)_CellLayout_New(  CELLLAYOUT_PASSARGS  );
	
	/* Virtual info */
	
	/* SingleCellLayout info */
	
	return self;
}

void _SingleCellLayout_Init( void* cellLayout, const Bool dimExists[3], const XYZ min, const XYZ max ) { 
	SingleCellLayout* self        = (SingleCellLayout*) cellLayout;
	Dimension_Index   dim_I;
	XYZ               minDefaults = { -1.0, -1.0, -1.0 };
	XYZ               maxDefaults = { 1.0, 1.0, 1.0 };
	
	/* General and Virtual info should already be set */

	self->isConstructed = True;
	for ( dim_I=0; dim_I < 3; dim_I++ ) {
		self->dimExists[dim_I] = dimExists[dim_I];
	}
	
	/* Get min and max values */
	if ( min )
		memcpy( self->min, min, sizeof( XYZ ));
	else
		memcpy( self->min, minDefaults, sizeof( XYZ ));
	if ( max )
		memcpy( self->max, max, sizeof( XYZ ));
	else
		memcpy( self->max, maxDefaults, sizeof( XYZ ));

	for ( dim_I=0; dim_I < 3; dim_I++ ) {
		if ( !dimExists[dim_I] ) {
			self->min[dim_I] = 0.0;
			self->max[dim_I] = 0.0;
		}
	}	

	/* Since pointcount and cell points are the same for all cells, calculate them now */
	_SingleCellLayout_CalculateGlobalPointCount( self );
	_SingleCellLayout_InitialiseGlobalCellPointPositions( self );
}


void _SingleCellLayout_Delete( void* singleCellLayout ) {
	SingleCellLayout* self = (SingleCellLayout*)singleCellLayout;
	
	/* Stg_Class_Delete parent class */
	_CellLayout_Delete( self );
}

void _SingleCellLayout_Print( void* singleCellLayout, Stream* stream ) {
	SingleCellLayout* self = (SingleCellLayout*)singleCellLayout;
	Index i;
	
	/* Set the Journal for printing informations */
	Stream* singleCellLayoutStream = stream;
	
	/* General info */
	Journal_Printf( singleCellLayoutStream, "SingleCellLayout (ptr): %p\n", self ); 
	
	/* Parent class info */
	_CellLayout_Print( self, stream );
	
	/* Virtual info */
	
	/* SingleCellLayout info */
	Journal_Printf( singleCellLayoutStream, "self->dimExists[" );
	for (i=I_AXIS; i < 3; i++ ) {
		Journal_Printf( singleCellLayoutStream, " %d,", self->dimExists[i] );
	}
	Journal_Printf( singleCellLayoutStream, "]\n" );
	
	Journal_Printf( singleCellLayoutStream, "self->min[" );
	for (i=I_AXIS; i < 3; i++ ) {
		Journal_Printf( singleCellLayoutStream, " %f,", self->min[i] );
	}	
	Journal_Printf( singleCellLayoutStream, "]\n" );
	
	Journal_Printf( singleCellLayoutStream, "self->max[" );
	for (i=I_AXIS; i < 3; i++ ) {
		Journal_Printf( singleCellLayoutStream, " %f,", self->max[i] );
	}	
	Journal_Printf( singleCellLayoutStream, "]\n" );
}


void* _SingleCellLayout_Copy( void* singleCellLayout, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	SingleCellLayout*	self = (SingleCellLayout*)singleCellLayout;
	SingleCellLayout*	newSingleCellLayout;
	PtrMap*			map = ptrMap;
	Bool			ownMap = False;
	
	if( !map ) {
		map = PtrMap_New( 10 );
		ownMap = True;
	}
	
	newSingleCellLayout = _CellLayout_Copy( self, dest, deep, nameExt, ptrMap );
	
	newSingleCellLayout->dimExists[0] = self->dimExists[0];
	newSingleCellLayout->dimExists[1] = self->dimExists[1];
	newSingleCellLayout->dimExists[2] = self->dimExists[2];
	newSingleCellLayout->min[0] = self->min[0];
	newSingleCellLayout->min[1] = self->min[1];
	newSingleCellLayout->min[2] = self->min[2];
	newSingleCellLayout->max[0] = self->max[0];
	newSingleCellLayout->max[1] = self->max[1];
	newSingleCellLayout->max[2] = self->max[2];
	newSingleCellLayout->pointCount = self->pointCount;
	
	if( deep ) {
		if( (newSingleCellLayout->cellPointCoords = PtrMap_Find( map, self->cellPointCoords )) == NULL && self->cellPointCoords ) {
			unsigned	p_i;

			newSingleCellLayout->cellPointCoords = Memory_Alloc_2DArray( double, newSingleCellLayout->pointCount, 3, "SingleCellLayout->cellPoints" );
			for( p_i = 0; p_i < newSingleCellLayout->pointCount; p_i++ )
				memcpy( newSingleCellLayout->cellPointCoords[p_i], self->cellPointCoords[p_i], 
					sizeof(double) * 3 * newSingleCellLayout->pointCount );
			PtrMap_Append( map, self->cellPointCoords, newSingleCellLayout->cellPointCoords );
		}
	}
	else {
		newSingleCellLayout->cellPointCoords = self->cellPointCoords;
	}
	
	if( ownMap ) {
		Stg_Class_Delete( map );
	}
	
	return (void*)newSingleCellLayout;
}

void _SingleCellLayout_AssignFromXML( void* singleCellLayout, Stg_ComponentFactory* cf, void* data ){
	SingleCellLayout* self              = (SingleCellLayout*)singleCellLayout;
	Bool              dimExists[]       = { False, False, False };
	Dimension_Index   dim;
	XYZ               min;
	XYZ               max;

	_CellLayout_AssignFromXML( self, cf, data );

	dim = Stg_ComponentFactory_GetRootDictUnsignedInt( cf, "dim", 0 );

	dimExists[ I_AXIS ] = Stg_ComponentFactory_GetBool( cf, self->name, "dimExistsI", True );
	dimExists[ J_AXIS ] = Stg_ComponentFactory_GetBool( cf, self->name, "dimExistsJ", True );
	dimExists[ K_AXIS ] = Stg_ComponentFactory_GetBool( cf, self->name, "dimExistsK", (dim == 3) );

	min[ I_AXIS ] = Stg_ComponentFactory_GetDouble( cf, self->name, "minX", -1.0 );
	min[ J_AXIS ] = Stg_ComponentFactory_GetDouble( cf, self->name, "minY", -1.0 );
	min[ K_AXIS ] = Stg_ComponentFactory_GetDouble( cf, self->name, "minZ", -1.0 );

	max[ I_AXIS ] = Stg_ComponentFactory_GetDouble( cf, self->name, "maxX", 1.0 );
	max[ J_AXIS ] = Stg_ComponentFactory_GetDouble( cf, self->name, "maxY", 1.0 );
	max[ K_AXIS ] = Stg_ComponentFactory_GetDouble( cf, self->name, "maxZ", 1.0 );

	_SingleCellLayout_Init( self, dimExists, min, max );
}
	
void _SingleCellLayout_Build( void* singleCellLayout, void* data ){
	
}
	
void _SingleCellLayout_Initialise( void* singleCellLayout, void* data ){
	
}
	
void _SingleCellLayout_Execute( void* singleCellLayout, void* data ){
	
}
	
void _SingleCellLayout_Destroy( void* singleCellLayout, void* data ){
	SingleCellLayout* self = (SingleCellLayout*)singleCellLayout;

	Memory_Free( self->cellPointCoords );

	_CellLayout_Destroy( self, data );
}

Cell_Index _SingleCellLayout_CellLocalCount( void* singleCellLayout ) {
	/* There is only one cell... */
	return 1;
}


Cell_Index _SingleCellLayout_CellShadowCount( void* singleCellLayout ) {
	/* No shadow cells */
	return 0;
}


void _SingleCellLayout_CalculateGlobalPointCount( SingleCellLayout* self ) {
	Index dim_I;

	self->pointCount = 1;
	for ( dim_I = I_AXIS; dim_I < 3; dim_I++ ) {
		if ( self->dimExists[dim_I] )
			self->pointCount *= 2;
	}
}


Cell_PointIndex _SingleCellLayout_PointCount( void* singleCellLayout, Cell_Index cellIndex ) {
	SingleCellLayout* self = (SingleCellLayout*)singleCellLayout;

	/* already calculated, just return that value */
	return self->pointCount;
}


void _SingleCellLayout_InitialiseGlobalCellPointPositions( SingleCellLayout* self ) {
	Cell_PointIndex	point_I = 0;
	Coord		tempCoord;
	double*		currPointCoord = NULL;
	Index		i, j, k;	/* loop iterators for each dimension */
	
	tempCoord[0] = self->min[I_AXIS];
	tempCoord[1] = self->min[J_AXIS];
	tempCoord[2] = self->min[K_AXIS];
	 
	self->cellPointCoords = Memory_Alloc_2DArray( double, self->pointCount, 3, "SingleCellLayout->cellPoints" );
	
	/* Now generate the coordinates */
	for ( k=0; k <= self->dimExists[K_AXIS]; k++ ) {
		for (j=0; j <= self->dimExists[J_AXIS]; j++ ) {
			for (i=0; i <= self->dimExists[I_AXIS]; i++ ) {
				currPointCoord = self->cellPointCoords[ RegularMeshUtils_ascendingIJK_ToHughesNodeNumberMap[point_I++] ];
				currPointCoord[I_AXIS] = tempCoord[I_AXIS];
				currPointCoord[J_AXIS] = tempCoord[J_AXIS];
				currPointCoord[K_AXIS] = tempCoord[K_AXIS];

				/* flip/flop the i for next time */
				tempCoord[I_AXIS] = ( self->min[I_AXIS] == tempCoord[I_AXIS] ) ? self->max[I_AXIS] : self->min[I_AXIS];
			}

			/* flip/flop the j for next time */
			tempCoord[J_AXIS] = ( self->min[J_AXIS] == tempCoord[J_AXIS] ) ? self->max[J_AXIS] : self->min[J_AXIS];
		}
		/* flip/flop the k for next time */
		tempCoord[K_AXIS] = ( self->min[K_AXIS] == tempCoord[K_AXIS] ) ? self->max[K_AXIS] : self->min[K_AXIS];
	}
}


void _SingleCellLayout_InitialisePoints( void* singleCellLayout, Cell_Index cellIndex, Cell_PointIndex pointCount, 
		Cell_Points points ) 
{
	SingleCellLayout* self = (SingleCellLayout*)singleCellLayout;
	Cell_PointIndex point_I = 0;

	/* since points have been pre-calculated, just return pointers to them */
	for ( point_I=0; point_I < pointCount; point_I++ ) {
		points[point_I] = &self->cellPointCoords[point_I];
	}
}


Cell_Index _SingleCellLayout_MapElementIdToCellId( void* cellLayout, Element_DomainIndex element_dI ) {
	
	/* Always 0: see the header comment */
	return 0;
}


Bool _SingleCellLayout_IsInCell( void* singleCellLayout, Cell_Index cellIndex, void* particle ) {
	SingleCellLayout* self = (SingleCellLayout*)singleCellLayout;
	double** coord = (double**)particle;
	Index dim_I = 0;
	
	for (dim_I=0; dim_I < 3; dim_I++ ) {
		if ( self->dimExists[dim_I] ) {
			if ( ((*coord)[dim_I] < self->min[dim_I]) || ((*coord)[dim_I] > self->max[dim_I]) ) {  
				return False;
			}	
		}
	}	
	
	return True;
}


Cell_Index _SingleCellLayout_CellOf( void* singleCellLayout, void* particle ) {
	
	/* in the single cell case, all particles belong to this cell */
	return 0;
}


ShadowInfo* _SingleCellLayout_GetShadowInfo( void* singleCellLayout ) {
	/* SingleCellLayout*      self = (SingleCellLayout*)singleCellLayout; */

	/* TODO: this should return a shadow info with at least nbr info for my processors */
	Journal_Firewall( 0, Swarm_Warning, "Error: %s not implemented yet!\n", __func__ );
	return NULL;
}


