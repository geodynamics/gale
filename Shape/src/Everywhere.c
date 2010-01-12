/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street, Melbourne, 3053, Australia.
**
** Authors:
**	Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**	Patrick D. Sunter, Software Engineer, VPAC. (pds@vpac.org)
**	Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**	Siew-Ching Tan, Software Engineer, VPAC. (siew@vpac.org) ) {
	IrregTopology* self = (IrregTopology*)ir
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
** $Id: Everywhere.c 3869 2006-10-16 13:42:59Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/Geometry/Geometry.h>

#include "types.h"
#include "ShapeClass.h"
#include "Everywhere.h"

#include <assert.h>
#include <string.h>
#include <math.h>

/* Textual name of this class */
const Type Everywhere_Type = "Everywhere";

/*-------------------------------------------------------------------------------------------------------------------------
** Constructors
*/
Everywhere* Everywhere_New(
		Name                                  name,
		Dimension_Index                       dim )
{
	Everywhere* self = (Everywhere*) _Everywhere_DefaultNew( name );
	XYZ         centre = { 0.0,0.0,0.0 };

	_Stg_Shape_Init( self, dim, centre, False, 0.0, 0.0, 0.0 );
	_Everywhere_Init( self );

	return self;
}

Everywhere* _Everywhere_New(  EVERYWHERE_DEFARGS  )
{
	Everywhere* self;
	
	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(Everywhere) );
	self = (Everywhere*)_Stg_Shape_New(  STG_SHAPE_PASSARGS  );
	
	/* General info */

	return self;
}

void _Everywhere_Init( void* everywhere ) {
}


/*------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _Everywhere_Delete( void* everywhere ) {
	Everywhere* self = (Everywhere*)everywhere;
	
	/* Delete parent */
	_Stg_Shape_Delete( self );
}

void _Everywhere_Print( void* everywhere, Stream* stream ) {
	Everywhere* self = (Everywhere*)everywhere;
	
	/* Print parent */
	_Stg_Shape_Print( self, stream );
}

void* _Everywhere_Copy( void* everywhere, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	Everywhere*	self = (Everywhere*)everywhere;
	Everywhere*	newEverywhere;
	
	newEverywhere = (Everywhere*)_Stg_Shape_Copy( self, dest, deep, nameExt, ptrMap );
	
	return (void*)newEverywhere;
}

void* _Everywhere_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                                  _sizeOfSelf = sizeof(Everywhere);
	Type                                                          type = Everywhere_Type;
	Stg_Class_DeleteFunction*                                  _delete = _Everywhere_Delete;
	Stg_Class_PrintFunction*                                    _print = _Everywhere_Print;
	Stg_Class_CopyFunction*                                      _copy = _Everywhere_Copy;
	Stg_Component_DefaultConstructorFunction*      _defaultConstructor = _Everywhere_DefaultNew;
	Stg_Component_ConstructFunction*                        _construct = _Everywhere_AssignFromXML;
	Stg_Component_BuildFunction*                                _build = _Everywhere_Build;
	Stg_Component_InitialiseFunction*                      _initialise = _Everywhere_Initialise;
	Stg_Component_ExecuteFunction*                            _execute = _Everywhere_Execute;
	Stg_Component_DestroyFunction*                            _destroy = _Everywhere_Destroy;
	Stg_Shape_IsCoordInsideFunction*                    _isCoordInside = _Everywhere_IsCoordInside;
	Stg_Shape_CalculateVolumeFunction*                _calculateVolume = _Everywhere_CalculateVolume;
	Stg_Shape_DistanceFromCenterAxisFunction*  _distanceFromCenterAxis = _Everywhere_DistanceFromCenterAxis;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

	return (void*) _Everywhere_New(  EVERYWHERE_PASSARGS  );
}


void _Everywhere_AssignFromXML( void* everywhere, Stg_ComponentFactory* cf, void* data ) {
	Everywhere*	self          = (Everywhere*) everywhere;

	_Stg_Shape_AssignFromXML( self, cf, data );
	_Everywhere_Init( self );
}

void _Everywhere_Build( void* everywhere, void* data ) {
	Everywhere*	self = (Everywhere*)everywhere;
	_Stg_Shape_Build( self, data );
}
void _Everywhere_Initialise( void* everywhere, void* data ) {
	Everywhere*	self = (Everywhere*)everywhere;
	_Stg_Shape_Initialise( self, data );
}
void _Everywhere_Execute( void* everywhere, void* data ) {
	Everywhere*	self = (Everywhere*)everywhere;
	_Stg_Shape_Execute( self, data );
}
void _Everywhere_Destroy( void* everywhere, void* data ) {
	Everywhere*	self = (Everywhere*)everywhere;
	_Stg_Shape_Destroy( self, data );
}

/*--------------------------------------------------------------------------------------------------------------------------
** Public Functions
*/
	
/*---------------------------------------------------------------------------------------------------------------------
** Private Member functions
*/

Bool _Everywhere_IsCoordInside( void* everywhere, Coord coord ) {
	return True;
}

double _Everywhere_CalculateVolume( void* everywhere ) {
	return 1.0;
}	
void _Everywhere_DistanceFromCenterAxis( void* shape, Coord coord, double* disVec ){
	Stg_Shape* self = (Stg_Shape*)shape;
	Journal_Firewall( False, Journal_Register( Error_Type, self->type ),
	"Error in function %s: This functions hasn't been implemented.", 
	"Please inform underworld-dev@vpac.org you've received this error.\n", __func__ );
}

/*----------------------------------------------------------------------------------------------------------------------------------
** Private Functions
*/



