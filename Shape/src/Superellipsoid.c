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
** $Id: Superellipsoid.c 3869 2006-10-16 13:42:59Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/Geometry/Geometry.h>

#include "types.h"
#include "ShapeClass.h"
#include "Superellipsoid.h"

#include <assert.h>
#include <string.h>
#include <math.h>


/* Textual name of this class */
const Type Superellipsoid_Type = "Superellipsoid";

/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/

Superellipsoid* Superellipsoid_New(
		Name                                  name,
		Dimension_Index                       dim,
		XYZ                                   centre, 
		double                                alpha,
		double                                beta,
		double                                gamma,
		double                                epsilon1,
		double                                epsilon2,   
		XYZ                                   radius )
{
	Superellipsoid* self = (Superellipsoid*) _Superellipsoid_DefaultNew( name );

   _Stg_Shape_Init( self, dim, centre, False, alpha, beta, gamma );
   _Superellipsoid_Init( self, epsilon1, epsilon2, radius );

	return self;
}

Superellipsoid* _Superellipsoid_New(  SUPERELLIPSOID_DEFARGS  )
{
	Superellipsoid* self;
	
	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(Superellipsoid) );
	self = (Superellipsoid*)_Stg_Shape_New(  STG_SHAPE_PASSARGS  );
	
	/* General info */

	return self;
}

void _Superellipsoid_Init( void* superellipsoid, double epsilon1, double epsilon2, XYZ radius ) {
	Superellipsoid* self = (Superellipsoid*)superellipsoid;
	
	self->epsilon1 = epsilon1;
	self->epsilon2 = epsilon2;

	memcpy( self->radius, radius, sizeof(XYZ));
}



/*------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _Superellipsoid_Delete( void* superellipsoid ) {
	Superellipsoid* self = (Superellipsoid*)superellipsoid;
	
	/* Delete parent */
	_Stg_Shape_Delete( self );
}


void _Superellipsoid_Print( void* superellipsoid, Stream* stream ) {
	Superellipsoid* self = (Superellipsoid*)superellipsoid;
	
	/* Print parent */
	_Stg_Shape_Print( self, stream );
}



void* _Superellipsoid_Copy( void* superellipsoid, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	Superellipsoid*	self = (Superellipsoid*)superellipsoid;
	Superellipsoid*	newSuperellipsoid;
	
	newSuperellipsoid = (Superellipsoid*)_Stg_Shape_Copy( self, dest, deep, nameExt, ptrMap );
	
	newSuperellipsoid->epsilon1 = self->epsilon1;
	newSuperellipsoid->epsilon2 = self->epsilon2;
	
	memcpy( newSuperellipsoid->radius, self->radius, sizeof(XYZ));
	
	return (void*)newSuperellipsoid;
}

void* _Superellipsoid_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                                  _sizeOfSelf = sizeof(Superellipsoid);
	Type                                                          type = Superellipsoid_Type;
	Stg_Class_DeleteFunction*                                  _delete = _Superellipsoid_Delete;
	Stg_Class_PrintFunction*                                    _print = _Superellipsoid_Print;
	Stg_Class_CopyFunction*                                      _copy = _Superellipsoid_Copy;
	Stg_Component_DefaultConstructorFunction*      _defaultConstructor = _Superellipsoid_DefaultNew;
	Stg_Component_ConstructFunction*                        _construct = _Superellipsoid_AssignFromXML;
	Stg_Component_BuildFunction*                                _build = _Superellipsoid_Build;
	Stg_Component_InitialiseFunction*                      _initialise = _Superellipsoid_Initialise;
	Stg_Component_ExecuteFunction*                            _execute = _Superellipsoid_Execute;
	Stg_Component_DestroyFunction*                            _destroy = _Superellipsoid_Destroy;
	Stg_Shape_IsCoordInsideFunction*                    _isCoordInside = _Superellipsoid_IsCoordInside;
	Stg_Shape_CalculateVolumeFunction*                _calculateVolume = _Superellipsoid_CalculateVolume;
	Stg_Shape_DistanceFromCenterAxisFunction*  _distanceFromCenterAxis = _Superellipsoid_DistanceFromCenterAxis;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

	return (void*) _Superellipsoid_New(  SUPERELLIPSOID_PASSARGS  );
}


void _Superellipsoid_AssignFromXML( void* superellipsoid, Stg_ComponentFactory* cf, void* data ) {
	Superellipsoid*	self      = (Superellipsoid*) superellipsoid;
	XYZ             radius;
	double          epsilon1;
	double          epsilon2;

	_Stg_Shape_AssignFromXML( self, cf, data );

	radius[ I_AXIS ] = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"radiusX", 1.0  );
	radius[ J_AXIS ] = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"radiusY", 1.0  );
	radius[ K_AXIS ] = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"radiusZ", 1.0  );

	epsilon1 = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"epsilon1", 1.0  );
	epsilon2 = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"epsilon2", 1.0  );

	_Superellipsoid_Init( self, epsilon1, epsilon2, radius );
}

void _Superellipsoid_Build( void* superellipsoid, void* data ) {
	Superellipsoid*	self = (Superellipsoid*)superellipsoid;

	_Stg_Shape_Build( self, data );
}
void _Superellipsoid_Initialise( void* superellipsoid, void* data ) {
	Superellipsoid*	self = (Superellipsoid*)superellipsoid;
	
	_Stg_Shape_Initialise( self, data );
}
void _Superellipsoid_Execute( void* superellipsoid, void* data ) {
	Superellipsoid*	self = (Superellipsoid*)superellipsoid;
	
	_Stg_Shape_Execute( self, data );
}
void _Superellipsoid_Destroy( void* superellipsoid, void* data ) {
	Superellipsoid*	self = (Superellipsoid*)superellipsoid;
	
	_Stg_Shape_Destroy( self, data );
}

/*--------------------------------------------------------------------------------------------------------------------------
** Public Functions
*/
	
/*---------------------------------------------------------------------------------------------------------------------
** Private Member functions
*/

Bool _Superellipsoid_IsCoordInside( void* superellipsoid, Coord coord ) {
	Superellipsoid* self              = (Superellipsoid*)superellipsoid;
	Coord           newCoord;
	double          insideOutsideValue;
	double          epsilon1          = self->epsilon1;
	double          epsilon2          = self->epsilon2;
	double          x, y, z;

	/* Transform coordinate into canonical reference frame */
	Stg_Shape_TransformCoord( self, coord, newCoord );
	
	x = newCoord[ I_AXIS ]/self->radius[ I_AXIS ];
	y = newCoord[ J_AXIS ]/self->radius[ J_AXIS ];
	z = newCoord[ K_AXIS ]/self->radius[ K_AXIS ];

	/* Evaluate Inside-Outside Function */
	if (self->dim == 2)
		insideOutsideValue = pow( x*x, 1.0/epsilon1 ) + pow( y*y, 1.0/epsilon1 );
	else 
		insideOutsideValue = pow( pow( x*x, 1.0/epsilon2 ) + pow( y*y, 1.0/epsilon2 ) , epsilon2/epsilon1 )
			+ pow( z*z, 1.0/epsilon1 );

	/* Return True if coord is inside and False otherwise */
	return ( insideOutsideValue <= 1.0 );
}

double _Superellipsoid_CalculateVolume( void* superellipsoid ) {
	assert( 0  );
	return 0.0;
}

void _Superellipsoid_DistanceFromCenterAxis( void* shape, Coord coord, double* disVec ){
	Stg_Shape* self = (Stg_Shape*)shape;
	Journal_Firewall( False, Journal_Register( Error_Type, (Name)self->type  ),
	"Error in function %s: This functions hasn't been implemented.", 
	"Please inform underworld-dev@vpac.org you've received this error.\n", __func__ );
}	


