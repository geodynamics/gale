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
** $Id: Cylinder.c 3869 2006-10-16 13:42:59Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/Geometry/Geometry.h>

#include "types.h"
#include "ShapeClass.h"
#include "Cylinder.h"

#include <assert.h>
#include <string.h>
#include <math.h>


/* Textual name of this class */
const Type Cylinder_Type = "Cylinder";

/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/
Cylinder* Cylinder_New(
		Name                                  name,
		Dimension_Index                       dim,
		XYZ                                   centre, 
		double                                alpha,
		double                                beta,
		double                                gamma,
		double                                radius, 
		XYZ                                   start, 
		XYZ                                   end, 
		Axis                                  alongAxis )
{
	Cylinder* self = (Cylinder*) _Cylinder_DefaultNew( name );

   _Stg_Shape_Init( self, dim, centre, False, alpha, beta, gamma );
	_Cylinder_Init( self, radius, start, end, alongAxis );

	return self;
}

Cylinder* _Cylinder_New(  CYLINDER_DEFARGS  )
{
	Cylinder* self;
	
	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(Cylinder) );
	self = (Cylinder*)_Stg_Shape_New(  STG_SHAPE_PASSARGS  );
	
	/* General info */

	/* Virtual Info */
	
	return self;
}

void _Cylinder_Init( Cylinder* self, double radius, XYZ start, XYZ end, Axis alongAxis ) {
	memcpy( self->start, start, sizeof(XYZ));
	memcpy( self->end, end, sizeof(XYZ));
	self->alongAxis = alongAxis;
	self->radius = radius;
}

/*------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _Cylinder_Delete( void* cylinder ) {
	Cylinder* self = (Cylinder*)cylinder;
	
	/* Delete parent */
	_Stg_Shape_Delete( self );
}

void _Cylinder_Print( void* cylinder, Stream* stream ) {
	Cylinder* self = (Cylinder*)cylinder;
	
	/* Print parent */
	_Stg_Shape_Print( self, stream );
}

void* _Cylinder_Copy( void* cylinder, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	Cylinder*	self = (Cylinder*)cylinder;
	Cylinder*	newCylinder;
	
	newCylinder = (Cylinder*)_Stg_Shape_Copy( self, dest, deep, nameExt, ptrMap );
	
	memcpy( newCylinder->start, self->start, sizeof(XYZ));
	memcpy( newCylinder->end, self->end, sizeof(XYZ));

	newCylinder->radius = self->radius;
	newCylinder->alongAxis = self->alongAxis;
	
	return (void*)newCylinder;
}

void* _Cylinder_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                                  _sizeOfSelf = sizeof(Cylinder);
	Type                                                          type = Cylinder_Type;
	Stg_Class_DeleteFunction*                                  _delete = _Cylinder_Delete;
	Stg_Class_PrintFunction*                                    _print = _Cylinder_Print;
	Stg_Class_CopyFunction*                                      _copy = _Cylinder_Copy;
	Stg_Component_DefaultConstructorFunction*      _defaultConstructor = _Cylinder_DefaultNew;
	Stg_Component_ConstructFunction*                        _construct = _Cylinder_AssignFromXML;
	Stg_Component_BuildFunction*                                _build = _Cylinder_Build;
	Stg_Component_InitialiseFunction*                      _initialise = _Cylinder_Initialise;
	Stg_Component_ExecuteFunction*                            _execute = _Cylinder_Execute;
	Stg_Component_DestroyFunction*                            _destroy = _Cylinder_Destroy;
	Stg_Shape_IsCoordInsideFunction*                    _isCoordInside = _Cylinder_IsCoordInside;
	Stg_Shape_CalculateVolumeFunction*                _calculateVolume = _Cylinder_CalculateVolume;
	Stg_Shape_DistanceFromCenterAxisFunction*  _distanceFromCenterAxis = _Cylinder_DistanceFromCenterAxis;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

	return (void*) _Cylinder_New(  CYLINDER_PASSARGS  );
}

#define BIG 1.0e99

void _Cylinder_AssignFromXML( void* cylinder, Stg_ComponentFactory* cf, void* data ) {
	Cylinder*            self                     = (Cylinder*) cylinder;
	XYZ                  start                    = { -BIG, -BIG, -BIG }; 
	XYZ                  end                      = {  BIG,  BIG,  BIG };
	double               radius                   = 0.0;
	Axis                 alongAxis        = I_AXIS;
	char*                perpendicularAxisName    = NULL;

	_Stg_Shape_AssignFromXML( self, cf, data );
	
	radius = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"radius", 0.0  );

	start[ I_AXIS ] = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"startX", -BIG  );
	start[ J_AXIS ] = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"startY", -BIG  );
	start[ K_AXIS ] = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"startZ", -BIG  );
	end[ I_AXIS ] = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"endX", BIG  );
	end[ J_AXIS ] = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"endY", BIG  );
	end[ K_AXIS ] = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"endZ", BIG  );

	perpendicularAxisName = Stg_ComponentFactory_GetString( cf, self->name, (Dictionary_Entry_Key)"perpendicularAxis", "x" );
	perpendicularAxisName = Stg_ComponentFactory_GetString( cf, self->name, (Dictionary_Entry_Key)"alongAxis", perpendicularAxisName );
	switch ( perpendicularAxisName[0] ) {
		case 'x': case 'X': case 'i': case 'I': case '0':
			alongAxis = I_AXIS; break;
		case 'y': case 'Y': case 'j': case 'J': case '1':
			alongAxis = J_AXIS; break;
		case 'z': case 'Z': case 'k': case 'K': case '2':
			alongAxis = K_AXIS; break;
		default:
			Journal_Firewall( False, Journal_Register( Error_Type, (Name)self->type  ),
					"Cannot understand alongAxis '%s'\n", perpendicularAxisName );
	}
	
	_Cylinder_Init( self, radius, start, end, alongAxis );
}

void _Cylinder_Build( void* cylinder, void* data ) {
	Cylinder*	self = (Cylinder*)cylinder;

	_Stg_Shape_Build( self, data );
}
void _Cylinder_Initialise( void* cylinder, void* data ) {
	Cylinder*	self = (Cylinder*)cylinder;
	
	_Stg_Shape_Initialise( self, data );
}
void _Cylinder_Execute( void* cylinder, void* data ) {
	Cylinder*	self = (Cylinder*)cylinder;
	
	_Stg_Shape_Execute( self, data );
}
void _Cylinder_Destroy( void* cylinder, void* data ) {
	Cylinder*	self = (Cylinder*)cylinder;

	_Stg_Shape_Destroy( self, data );
}

/*--------------------------------------------------------------------------------------------------------------------------
** Public Functions
*/
	
/*---------------------------------------------------------------------------------------------------------------------
** Private Member functions
*/

Bool _Cylinder_IsCoordInside( void* cylinder, Coord coord ) {
	Cylinder*       self       = (Cylinder*)cylinder;
	Coord           newCoord;
	double          insideOutsideValue;
	double          x, y, z;
	Axis            axis_I;

	/* Check whether coord is within min and max values */
	for ( axis_I = 0 ; axis_I < self->dim ; axis_I++ ) {
		if ( coord[ axis_I ] < self->start[ axis_I ] || coord[ axis_I ] > self->end[ axis_I ] )
			return False;
	}
	
	/* Transform coordinate into canonical reference frame */
	Stg_Shape_TransformCoord( self, coord, newCoord );
	
	newCoord[ self->alongAxis ] = 0.0;

	/* Check if coord is within radius */
	x = newCoord[ I_AXIS ];
	y = newCoord[ J_AXIS ];
	if(self->dim == 2)
		insideOutsideValue = x*x + y*y;
	else {
		z = newCoord[ K_AXIS ];
		insideOutsideValue = x*x + y*y + z*z;
	}
	if ( insideOutsideValue > (self->radius * self->radius) )
		return False;


	return True;
}

void _Cylinder_DistanceFromCenterAxis( void* cylinder, Coord coord, double* disVec ) {
	Cylinder*       self       = (Cylinder*)cylinder;
	Coord           newCoord;

	/* Transform coordinate into canonical reference frame */
	Stg_Shape_TransformCoord( self, coord, newCoord );
	
	newCoord[ self->alongAxis ] = 0.0;

	/* Check if coord is within radius */
	disVec[0] = newCoord[ I_AXIS ];
	disVec[1] = newCoord[ J_AXIS ];
	if(self->dim == 3)
		disVec[2] = newCoord[ K_AXIS ];

	return;
}


double _Cylinder_CalculateVolume( void* cylinder ) {
	assert( 0 /* unsure how this cylinder is setup...but shouldn't be hard to implement -- Alan */ );
	return 0.0;
}



