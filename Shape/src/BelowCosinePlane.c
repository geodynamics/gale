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
** $Id: BelowCosinePlane.c 3523 2006-04-11 06:42:09Z AlanLo $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/Geometry/Geometry.h>

#include "types.h"
#include "ShapeClass.h"
#include "BelowPlane.h"
#include "BelowCosinePlane.h"

#include <assert.h>
#include <string.h>
#include <math.h>


/* Textual name of this class */
const Type BelowCosinePlane_Type = "BelowCosinePlane";

/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/
BelowCosinePlane* BelowCosinePlane_New(
		Name                                  name,
		Dimension_Index                       dim,
		XYZ                                   centre, 
		double                                alpha,
		double                                beta,
		double                                gamma,
		double                                offset,
		XYZ                                   width,
		XYZ                                   minValue,
		XYZ                                   maxValue,
		double                                amplitude,
		double                                wavelength,
		double                                phase )
{
	BelowCosinePlane* self = (BelowCosinePlane*) _BelowCosinePlane_DefaultNew( name );

   _Stg_Shape_Init( self, dim, centre, False, alpha, beta, gamma );
   _BelowPlane_Init( self, offset, width, minValue, maxValue );
   _BelowCosinePlane_Init( self, width, amplitude, wavelength, phase );

	return self;
}

BelowCosinePlane* _BelowCosinePlane_New(  BELOWCOSINEPLANE_DEFARGS  )
{
	BelowCosinePlane* self;
	
	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(BelowCosinePlane) );
	self = (BelowCosinePlane*)_BelowPlane_New(  BELOWPLANE_PASSARGS  );
	
	/* General info */
	return self;
}

void _BelowCosinePlane_Init( void* belowPlane, XYZ width, double amplitude, double wavelength, double phase ) {
	BelowCosinePlane* self = (BelowCosinePlane*)belowPlane;

	self->amplitude = amplitude;
	self->wavelength = wavelength;
	self->phase = phase;
}


/*------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _BelowCosinePlane_Delete( void* belowPlane ) {
	BelowCosinePlane* self = (BelowCosinePlane*)belowPlane;
	
	/* Delete parent */
	_BelowPlane_Delete( self );
}


void _BelowCosinePlane_Print( void* belowPlane, Stream* stream ) {
	BelowCosinePlane* self = (BelowCosinePlane*)belowPlane;
	
	/* Print parent */
	_Stg_Shape_Print( self, stream );
}

void* _BelowCosinePlane_Copy( void* belowPlane, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	BelowCosinePlane*	self = (BelowCosinePlane*)belowPlane;
	BelowCosinePlane*	newBelowCosinePlane;
	
	newBelowCosinePlane = (BelowCosinePlane*)_BelowPlane_Copy( self, dest, deep, nameExt, ptrMap );

	newBelowCosinePlane->amplitude = self->amplitude;
	newBelowCosinePlane->wavelength = self->wavelength;
	newBelowCosinePlane->phase = self->phase;
	
	return (void*)newBelowCosinePlane;
}

void* _BelowCosinePlane_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                                  _sizeOfSelf = sizeof(BelowCosinePlane);
	Type                                                          type = BelowCosinePlane_Type;
	Stg_Class_DeleteFunction*                                  _delete = _BelowCosinePlane_Delete;
	Stg_Class_PrintFunction*                                    _print = _BelowCosinePlane_Print;
	Stg_Class_CopyFunction*                                      _copy = _BelowCosinePlane_Copy;
	Stg_Component_DefaultConstructorFunction*      _defaultConstructor = _BelowCosinePlane_DefaultNew;
	Stg_Component_ConstructFunction*                        _construct = _BelowCosinePlane_AssignFromXML;
	Stg_Component_BuildFunction*                                _build = _BelowCosinePlane_Build;
	Stg_Component_InitialiseFunction*                      _initialise = _BelowCosinePlane_Initialise;
	Stg_Component_ExecuteFunction*                            _execute = _BelowCosinePlane_Execute;
	Stg_Component_DestroyFunction*                            _destroy = _BelowCosinePlane_Destroy;
	Stg_Shape_IsCoordInsideFunction*                    _isCoordInside = _BelowCosinePlane_IsCoordInside;
	Stg_Shape_CalculateVolumeFunction*                _calculateVolume = _BelowCosinePlane_CalculateVolume;
	Stg_Shape_DistanceFromCenterAxisFunction*  _distanceFromCenterAxis = _BelowCosinePlane_DistanceFromCenterAxis;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

	return (void*) _BelowCosinePlane_New(  BELOWCOSINEPLANE_PASSARGS  );
}


void _BelowCosinePlane_AssignFromXML( void* belowPlane, Stg_ComponentFactory* cf, void* data ) {
	BelowCosinePlane*            self          = (BelowCosinePlane*) belowPlane;
	double                       amplitude;
	double                       wavelength;
	double                       phase;

	_BelowPlane_AssignFromXML( self, cf, data );

	amplitude = Stg_ComponentFactory_GetDouble( cf, self->name, "amplitude", 0.1 );
	wavelength = Stg_ComponentFactory_GetDouble( cf, self->name, "wavelength", 2*M_PI );
	phase = Stg_ComponentFactory_GetDouble( cf, self->name, "phase", 0.0 );

	_BelowCosinePlane_Init( self, self->width, amplitude, wavelength, phase );
}

void _BelowCosinePlane_Build( void* belowPlane, void* data ) {
	BelowCosinePlane*	self = (BelowCosinePlane*)belowPlane;

	_BelowPlane_Build( self, data );
}
void _BelowCosinePlane_Initialise( void* belowPlane, void* data ) {
	BelowCosinePlane*	self = (BelowCosinePlane*)belowPlane;
	
	_BelowPlane_Initialise( self, data );
}
void _BelowCosinePlane_Execute( void* belowPlane, void* data ) {
	BelowCosinePlane*	self = (BelowCosinePlane*)belowPlane;
	
	_BelowPlane_Execute( self, data );
}
void _BelowCosinePlane_Destroy( void* belowPlane, void* data ) {
	BelowCosinePlane*	self = (BelowCosinePlane*)belowPlane;
    
	_BelowPlane_Destroy( self, data );
}

/*--------------------------------------------------------------------------------------------------------------------------
** Public Functions
*/
	
/*---------------------------------------------------------------------------------------------------------------------
** Private Member functions
*/

Bool _BelowCosinePlane_IsCoordInside( void* belowPlane, Coord coord ) {
	BelowCosinePlane*            self       = (BelowCosinePlane*)belowPlane;
	Coord                        newCoord;

	double                       x, y;

	/* Transform coordinate into canonical reference frame */
	Stg_Shape_TransformCoord( self, coord, newCoord );

	x = newCoord[ I_AXIS ];

	y =  self->offset + self->amplitude * cos( (2*M_PI * x /self->wavelength)  + self->phase );

	if ( fabs( newCoord[ J_AXIS ] < y) ) {
		return True;
	}
	return False;
}

double _BelowCosinePlane_CalculateVolume( void* belowPlane ) {
	BelowCosinePlane* self = (BelowCosinePlane*)belowPlane;
	double volume;
	double wavelength = self->wavelength;
	double dx = self->width[ I_AXIS ];

	/* using the identity sin(u)-sin(v) = 2 * cos( (u+v)/2 ) * sin( (u-v)/2 ) */

	volume = self->offset*dx + 2 * self->amplitude * wavelength / (2*M_PI) *
				 (
				 cos( (M_PI/wavelength) * (self->maxValue[I_AXIS] - self->minValue[I_AXIS] ) + self->phase) *
				 sin( (M_PI/wavelength) * (self->maxValue[I_AXIS] - self->minValue[I_AXIS]) )
				 );

	if ( self->dim == 3 ) 
		volume = self->width[ K_AXIS ] * volume;
	
	return volume;
}
void _BelowCosinePlane_DistanceFromCenterAxis( void* shape, Coord coord, double* disVec ) {
	Stg_Shape* self = (Stg_Shape*)shape;
	Journal_Firewall( False, Journal_Register( Error_Type, self->type ),
	"Error in function %s: This functions hasn't been implemented.", 
	"Please inform underworld-dev@vpac.org you've received this error.\n", __func__ );
}



