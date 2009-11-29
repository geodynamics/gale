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
** $Id: Box.c 3869 2006-10-16 13:42:59Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/Geometry/Geometry.h>

#include "types.h"
#include "ShapeClass.h"
#include "Box.h"

#include <assert.h>
#include <string.h>
#include <math.h>


/* Textual name of this class */
const Type Box_Type = "Box";

/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/
Box* Box_New(
		Name                                  name,
		Dimension_Index                       dim,
		XYZ                                   centre, 
		double                                alpha,
		double                                beta,
		double                                gamma,
		XYZ                                   width )
{
	Box* self = (Box*) _Box_DefaultNew( name );

   _Stg_Shape_Init( self, dim, centre, False, alpha, beta, gamma );
	_Box_Init( self, width );
	return self;
}

Box* _Box_New(  BOX_DEFARGS  )
{
	Box* self;
	
	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(Box) );
	self = (Box*)_Stg_Shape_New(  STG_SHAPE_PASSARGS  );
	
	/* General info */

	return self;
}

void _Box_Init( void* shape, XYZ width ) {
	Box* self = (Box*)shape;
	
	memcpy( self->width, width, sizeof(XYZ));
}


/*------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _Box_Delete( void* shape ) {
	Box* self = (Box*)shape;
	
	/* Delete parent */
	_Stg_Shape_Delete( self );
}


void _Box_Print( void* shape, Stream* stream ) {
	Box* self = (Box*)shape;
	
	/* Print parent */
	_Stg_Shape_Print( self, stream );
}

void* _Box_Copy( void* shape, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	Box*	self = (Box*)shape;
	Box*	newBox;
	
	newBox = (Box*)_Stg_Shape_Copy( self, dest, deep, nameExt, ptrMap );
	
	memcpy( newBox->width, self->width, sizeof(XYZ));
	
	return (void*)newBox;
}

void* _Box_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                                  _sizeOfSelf = sizeof(Box);
	Type                                                          type = Box_Type;
	Stg_Class_DeleteFunction*                                  _delete = _Box_Delete;
	Stg_Class_PrintFunction*                                    _print = _Box_Print;
	Stg_Class_CopyFunction*                                      _copy = _Box_Copy;
	Stg_Component_DefaultConstructorFunction*      _defaultConstructor = _Box_DefaultNew;
	Stg_Component_ConstructFunction*                        _construct = _Box_AssignFromXML;
	Stg_Component_BuildFunction*                                _build = _Box_Build;
	Stg_Component_InitialiseFunction*                      _initialise = _Box_Initialise;
	Stg_Component_ExecuteFunction*                            _execute = _Box_Execute;
	Stg_Component_DestroyFunction*                            _destroy = _Box_Destroy;
	Stg_Shape_IsCoordInsideFunction*                    _isCoordInside = _Box_IsCoordInside;
	Stg_Shape_CalculateVolumeFunction*                _calculateVolume = _Box_CalculateVolume;
	Stg_Shape_DistanceFromCenterAxisFunction*  _distanceFromCenterAxis = _Box_DistanceFromCenterAxis;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = ZERO;

	return (void*) _Box_New(  BOX_PASSARGS  );
}


void _Box_AssignFromXML( void* shape, Stg_ComponentFactory* cf, void* data ) {
	Box*	             self          = (Box*) shape;
	Dictionary*          dictionary    = Dictionary_GetDictionary( cf->componentDict, self->name );
	XYZ                  width;
	double               start, end;
	Dictionary_Entry_Key startKey      = StG_Strdup("startX");
	Dictionary_Entry_Key endKey        = StG_Strdup("endX");
	Dictionary_Entry_Key widthKey      = StG_Strdup("widthX");
	char*                startCharPtr  = strchr( startKey, 'X' );
	char*                endCharPtr    = strchr( endKey, 'X' );
	char*                widthCharPtr  = strchr( widthKey, 'X' );
	char                 axisLetters[] = {'X','Y','Z'};
	Dimension_Index      dim_I;

	_Stg_Shape_AssignFromXML( self, cf, data );

	for ( dim_I = 0 ; dim_I < 3 ; dim_I++ ) {
		*startCharPtr = axisLetters[ dim_I ];
		*endCharPtr   = axisLetters[ dim_I ];
		*widthCharPtr = axisLetters[ dim_I ];

		/* Check to see whether the user wants to specify the start and end explicitly */
		if ( Dictionary_Get( dictionary, startKey ) && Dictionary_Get( dictionary, endKey ) ) {
			start = Stg_ComponentFactory_GetDouble( cf, self->name, startKey, 0.0 );
			end   = Stg_ComponentFactory_GetDouble( cf, self->name, endKey,   0.0 );

			width[ dim_I ] = end - start;
			self->centre[ dim_I ] = start + 0.5 * width[dim_I];
		}
		else 
			width[ dim_I ] = Stg_ComponentFactory_GetDouble( cf, self->name, widthKey, 0.0 );
	}

	Memory_Free( startKey );
	Memory_Free( endKey );
	Memory_Free( widthKey );

	_Box_Init( self, width );
}

void _Box_Build( void* shape, void* data ) {
	Box*	self = (Box*)shape;

	_Stg_Shape_Build( self, data );
}
void _Box_Initialise( void* shape, void* data ) {
	Box*	self = (Box*)shape;
	
	_Stg_Shape_Initialise( self, data );
}
void _Box_Execute( void* shape, void* data ) {
	Box*	self = (Box*)shape;
	
	_Stg_Shape_Execute( self, data );
}
void _Box_Destroy( void* shape, void* data ) {
	Box*	self = (Box*)shape;
    
	_Stg_Shape_Destroy( self, data );
}

/*--------------------------------------------------------------------------------------------------------------------------
** Public Functions
*/
	
/*---------------------------------------------------------------------------------------------------------------------
** Private Member functions
*/

Bool _Box_IsCoordInside( void* shape, Coord coord ) {
	Box*            self       = (Box*)shape;
	Coord           newCoord;
	Dimension_Index dim_I;

	/* Transform coordinate into canonical reference frame */
	Stg_Shape_TransformCoord( self, coord, newCoord );
	
	for ( dim_I = 0 ; dim_I < self->dim ; dim_I++ ) {
		if ( fabs( newCoord[ dim_I ] ) > 0.5 * self->width[ dim_I ] ) 
			return False;
	}
	return True;
}

double _Box_CalculateVolume( void* shape ) {
	Box* self = (Box*)shape;
	Dimension_Index dim_I;
	double result;
	result = 1.0;
	for ( dim_I = 0; dim_I < self->dim; dim_I++ ) {
		result *= self->width[dim_I];
	}
	return result;
}

void _Box_DistanceFromCenterAxis( void* shape, Coord coord, double* disVec ) {
	/* To be implemented */
	Stg_Shape* self = (Stg_Shape*)shape;
	Journal_Firewall( False, Journal_Register( Error_Type, self->type ),
	"Error in function %s: This functions hasn't been implemented.", 
	"Please inform underworld-dev@vpac.org you've received this error.\n", __func__ );
}



