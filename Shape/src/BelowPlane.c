/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street, Melbourne, 3053, Australia.
**
** Authors:
**	Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**	Patrick D. Sunter, Software Engineer, VPAC. (pds@vpac.org)
**	Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**	Siew-Ching Tan, Software Engineer, VPAC. (siew@vpac.org) )
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
** $Id: BelowPlane.c 3523 2006-04-11 06:42:09Z AlanLo $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/Geometry/Geometry.h>

#include "types.h"
#include "ShapeClass.h"
#include "BelowPlane.h"

#include <assert.h>
#include <string.h>
#include <math.h>


/* Textual name of this class */
const Type BelowPlane_Type = "BelowPlane";

/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/
BelowPlane* BelowPlane_New(
		Name                                  name,
		Dimension_Index                       dim,
		XYZ                                   centre, 
		double                                alpha,
		double                                beta,
		double                                gamma,
		double                                offset,
		XYZ                                   width,
		XYZ                                   minValue,
		XYZ                                   maxValue )
{
	BelowPlane* self = (BelowPlane*) _BelowPlane_DefaultNew( name );

   _Stg_Shape_Init( self, dim, centre, False, alpha, beta, gamma );
	_BelowPlane_Init( self, offset, width, minValue, maxValue );

	return self;
}

BelowPlane* _BelowPlane_New(  BELOWPLANE_DEFARGS  )
{
	BelowPlane* self;
	
	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(BelowPlane) );
	self = (BelowPlane*)_Stg_Shape_New(  STG_SHAPE_PASSARGS  );
	
	/* General info */
	return self;
}

void _BelowPlane_Init( void* belowPlane, double offset, XYZ width, XYZ minValue, XYZ maxValue ) {
	BelowPlane* self = (BelowPlane*)belowPlane;

	self->offset = offset;

	memcpy( self->width, width, sizeof(XYZ) );
	memcpy( self->minValue, maxValue, sizeof(XYZ) );
	memcpy( self->maxValue, maxValue, sizeof(XYZ) );
}

/*------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/
void _BelowPlane_Delete( void* belowPlane ) {
	BelowPlane* self = (BelowPlane*)belowPlane;
	
	/* Delete parent */
	_Stg_Shape_Delete( self );
}

void _BelowPlane_Print( void* belowPlane, Stream* stream ) {
	BelowPlane* self = (BelowPlane*)belowPlane;
	
	/* Print parent */
	_Stg_Shape_Print( self, stream );
}

void* _BelowPlane_Copy( void* belowPlane, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	BelowPlane*	self = (BelowPlane*)belowPlane;
	BelowPlane*	newBelowPlane;
	
	newBelowPlane = (BelowPlane*)_Stg_Shape_Copy( self, dest, deep, nameExt, ptrMap );

	newBelowPlane->offset = self->offset;
	
	return (void*)newBelowPlane;
}

void* _BelowPlane_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                                  _sizeOfSelf = sizeof(BelowPlane);
	Type                                                          type = BelowPlane_Type;
	Stg_Class_DeleteFunction*                                  _delete = _BelowPlane_Delete;
	Stg_Class_PrintFunction*                                    _print = _BelowPlane_Print;
	Stg_Class_CopyFunction*                                      _copy = _BelowPlane_Copy;
	Stg_Component_DefaultConstructorFunction*      _defaultConstructor = _BelowPlane_DefaultNew;
	Stg_Component_ConstructFunction*                        _construct = _BelowPlane_AssignFromXML;
	Stg_Component_BuildFunction*                                _build = _BelowPlane_Build;
	Stg_Component_InitialiseFunction*                      _initialise = _BelowPlane_Initialise;
	Stg_Component_ExecuteFunction*                            _execute = _BelowPlane_Execute;
	Stg_Component_DestroyFunction*                            _destroy = _BelowPlane_Destroy;
	Stg_Shape_IsCoordInsideFunction*                    _isCoordInside = _BelowPlane_IsCoordInside;
	Stg_Shape_CalculateVolumeFunction*                _calculateVolume = _BelowPlane_CalculateVolume;
	Stg_Shape_DistanceFromCenterAxisFunction*  _distanceFromCenterAxis = _BelowPlane_DistanceFromCenterAxis;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

	return (void*) _BelowPlane_New(  BELOWPLANE_PASSARGS  );
}


void _BelowPlane_AssignFromXML( void* belowPlane, Stg_ComponentFactory* cf, void* data ) {
	BelowPlane*	             self          = (BelowPlane*) belowPlane;
	double                       offset;
	XYZ                          minValue;
	XYZ                          maxValue;
	XYZ                          width;

	_Stg_Shape_AssignFromXML( self, cf, data );

	offset = Stg_ComponentFactory_GetDouble( cf, self->name, "offset", 0.5 );

	minValue[ I_AXIS ] = Stg_ComponentFactory_GetRootDictDouble( cf, "minX", 0.0 );
	minValue[ J_AXIS ] = Stg_ComponentFactory_GetRootDictDouble( cf, "minY", 0.0 );
	minValue[ K_AXIS ] = Stg_ComponentFactory_GetRootDictDouble( cf, "minZ", 0.0 );

	maxValue[ I_AXIS ] = Stg_ComponentFactory_GetRootDictDouble( cf, "maxX", 1.0 );
	maxValue[ J_AXIS ] = Stg_ComponentFactory_GetRootDictDouble( cf, "maxY", 1.0 );
	maxValue[ K_AXIS ] = Stg_ComponentFactory_GetRootDictDouble( cf, "maxZ", 1.0 );

	width[ I_AXIS ] = maxValue[ I_AXIS ] - minValue[ I_AXIS ] ;
	width[ J_AXIS ] = maxValue[ J_AXIS ] - minValue[ J_AXIS ] ;
	width[ K_AXIS ] = maxValue[ K_AXIS ] - minValue[ K_AXIS ] ;

	_BelowPlane_Init( self, offset, width, minValue, maxValue );
}

void _BelowPlane_Build( void* belowPlane, void* data ) {
	BelowPlane*	self = (BelowPlane*)belowPlane;

	_Stg_Shape_Build( self, data );
}
void _BelowPlane_Initialise( void* belowPlane, void* data ) {
	BelowPlane*	self = (BelowPlane*)belowPlane;
	
	_Stg_Shape_Initialise( self, data );
}
void _BelowPlane_Execute( void* belowPlane, void* data ) {
	BelowPlane*	self = (BelowPlane*)belowPlane;
	
	_Stg_Shape_Execute( self, data );
}
void _BelowPlane_Destroy( void* belowPlane, void* data ) {
	BelowPlane*	self = (BelowPlane*)belowPlane;
	
	_Stg_Shape_Destroy( self, data );
}

/*--------------------------------------------------------------------------------------------------------------------------
** Public Functions
*/
	
/*---------------------------------------------------------------------------------------------------------------------
** Private Member functions
*/

Bool _BelowPlane_IsCoordInside( void* belowPlane, Coord coord ) {
	BelowPlane*            self       = (BelowPlane*)belowPlane;
	Coord           newCoord;

	/* Transform coordinate into canonical reference frame */
	Stg_Shape_TransformCoord( self, coord, newCoord );

	if ( fabs( newCoord[ J_AXIS ] < self->offset ) ) {
		return True;
	}
	return False;
}

double _BelowPlane_CalculateVolume( void* belowPlane ) {
	BelowPlane* self = (BelowPlane*)belowPlane;
	double volume;

	if ( self->dim == 2 ) {
		volume = self->width[ I_AXIS ] * self->offset;
	}
	else {
		volume = self->width[ I_AXIS ] * self->width[ K_AXIS ] * self->offset;
	}

	return volume;
}

void _BelowPlane_DistanceFromCenterAxis( void* shape, Coord coord, double* disVec ){
	Stg_Shape* self = (Stg_Shape*)shape;
	Journal_Firewall( False, Journal_Register( Error_Type, self->type ),
	"Error in function %s: This functions hasn't been implemented.", 
	"Please inform underworld-dev@vpac.org you've received this error.\n", __func__ );
}



