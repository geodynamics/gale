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
** $Id: Union.c 4081 2007-04-27 06:20:07Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/Geometry/Geometry.h>

#include "types.h"
#include "ShapeClass.h"
#include "Union.h"

#include <assert.h>
#include <string.h>
#include <math.h>

/* Textual name of this class */
const Type Union_Type = "Union";

/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/

Union* Union_New(
		Name                                  name,
		Dimension_Index                       dim,
		XYZ                                   centre, 
		double                                alpha,
		double                                beta,
		double                                gamma,
		Stg_Shape**                           shapeList,
		Index                                 shapeCount,
		Bool*                                 isComplement
		)
{
	Union* self = (Union*)_Union_DefaultNew( name );

   _Stg_Shape_Init( self, dim, centre, False, alpha, beta, gamma);
	_Union_Init( self, shapeList, shapeCount, isComplement );
	return self;
}

Union* _Union_New(  UNION_DEFARGS  )
{
	Union* self;
	
	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(Union) );
	self = (Union*)_Stg_Shape_New(  STG_SHAPE_PASSARGS  );
	
	/* General info */
	
	return self;
}

void _Union_Init( void* combination,  Stg_Shape** shapeList, Index shapeCount, Bool* isComplement ) {
	Union* self = (Union*)combination;
	
	self->shapeList    = Memory_Alloc_Array( Stg_Shape* , shapeCount , "shapeList" );
	self->isComplement = Memory_Alloc_Array( Bool,        shapeCount , "isComplement" );

	memcpy( self->shapeList , shapeList, sizeof(Stg_Shape*) * shapeCount );
	memcpy( self->isComplement , isComplement, sizeof(Bool) * shapeCount );
	self->shapeCount = shapeCount;
}



/*------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _Union_Delete( void* combination ) {
	Union*       self = (Union*)combination;

	Memory_Free( self->shapeList );
	/* Delete parent */
	_Stg_Shape_Delete( self );
}


void _Union_Print( void* combination, Stream* stream ) {
	Union* self = (Union*)combination;
	
	/* Print parent */
	_Stg_Shape_Print( self, stream );
}

void* _Union_Copy( void* combination, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	Union*	self = (Union*)combination;
	Union*	newUnion;
	
	newUnion = (Union*)_Stg_Shape_Copy( self, dest, deep, nameExt, ptrMap );

	newUnion->shapeList = Memory_Alloc_Array( Stg_Shape*, self->shapeCount, "shapeList" );
	memcpy( newUnion->shapeList , self->shapeList, sizeof(Stg_Shape*) * self->shapeCount );
	
	newUnion->isComplement = Memory_Alloc_Array( Bool , self->shapeCount, "isComplement" );
	memcpy( newUnion->isComplement , self->isComplement, sizeof(Bool) * self->shapeCount );

	newUnion->isComplement = self->isComplement;
	newUnion->shapeList = self->shapeList;
	newUnion->shapeCount = self->shapeCount;
	
	return (void*)newUnion;
}

void* _Union_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                                  _sizeOfSelf = sizeof(Union);
	Type                                                          type = Union_Type;
	Stg_Class_DeleteFunction*                                  _delete = _Union_Delete;
	Stg_Class_PrintFunction*                                    _print = _Union_Print;
	Stg_Class_CopyFunction*                                      _copy = _Union_Copy;
	Stg_Component_DefaultConstructorFunction*      _defaultConstructor = _Union_DefaultNew;
	Stg_Component_ConstructFunction*                        _construct = _Union_AssignFromXML;
	Stg_Component_BuildFunction*                                _build = _Union_Build;
	Stg_Component_InitialiseFunction*                      _initialise = _Union_Initialise;
	Stg_Component_ExecuteFunction*                            _execute = _Union_Execute;
	Stg_Component_DestroyFunction*                            _destroy = _Union_Destroy;
	Stg_Shape_IsCoordInsideFunction*                    _isCoordInside = _Union_IsCoordInside;
	Stg_Shape_CalculateVolumeFunction*                _calculateVolume = _Union_CalculateVolume;
	Stg_Shape_DistanceFromCenterAxisFunction*  _distanceFromCenterAxis = _Union_DistanceFromCenterAxis;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

	return (void*) _Union_New(  UNION_PASSARGS  );
}


void _Union_AssignFromXML( void* combination, Stg_ComponentFactory* cf, void* data ) {
	Union*                  self       = (Union*)combination;
	Index                   shapeCount;
	Stg_Shape**             shapeList;
	Bool*                   isComplement;
	Index                   shape_I;
	Dictionary_Entry_Value* optionsList;
	Dictionary*             dictionary = Dictionary_GetDictionary( cf->componentDict, self->name );
	char*                   nameShape;
	Stream*                 stream     = Journal_Register( Info_Type, CURR_MODULE_NAME );
	
	_Stg_Shape_AssignFromXML( self, cf, data );

	optionsList = Dictionary_Get( dictionary, "shapes" );

	shapeCount = Dictionary_Entry_Value_GetCount(optionsList);

	/* Allocate space */
	shapeList     = Memory_Alloc_Array( Stg_Shape* , shapeCount, "Shape Array" );
	isComplement  = Memory_Alloc_Array( Bool, shapeCount, "Complement Array" );
	memset( shapeList,     0, shapeCount * sizeof(Stg_Shape*) );
	memset( isComplement,  0, shapeCount * sizeof(Bool) );
	
	Stream_Indent( stream );
	for ( shape_I = 0 ; shape_I < shapeCount ; shape_I++) { 
          /* gets the textual name corresponding to the shape elements	 */
		nameShape = Dictionary_Entry_Value_AsString( Dictionary_Entry_Value_GetElement( optionsList, shape_I));

		if ( nameShape[0] == '!' ) {
			shapeList[ shape_I ] =  Stg_ComponentFactory_ConstructByName( cf, &nameShape[1], Stg_Shape, True, data ) ;
			isComplement[ shape_I ] = True;
		}
		else {
			shapeList[ shape_I ] =  Stg_ComponentFactory_ConstructByName( cf, nameShape, Stg_Shape, True, data ) ;
			isComplement[ shape_I ] = False;
		}
		
	}
	Stream_UnIndent( stream );

	_Union_Init( self, shapeList, shapeCount, isComplement );

	Memory_Free( shapeList );
	Memory_Free( isComplement );
}

void _Union_Build( void* combination, void* data ) {
	Union*	self = (Union*)combination;
   unsigned shape_I = 0;
    
   for( shape_I = 0 ; shape_I < self->shapeCount ; shape_I++ ) {
      Stg_Component_Build( self->shapeList[shape_I], data, False );
   }

	_Stg_Shape_Build( self, data );
}
void _Union_Initialise( void* combination, void* data ) {
	Union*	self = (Union*)combination;
	unsigned shape_I = 0;
    
   for( shape_I = 0 ; shape_I < self->shapeCount ; shape_I++ ) {
      Stg_Component_Initialise( self->shapeList[shape_I], data, False );
   }
	_Stg_Shape_Initialise( self, data );
}
void _Union_Execute( void* combination, void* data ) {
	Union*	self = (Union*)combination;
	
	_Stg_Shape_Execute( self, data );
}
void _Union_Destroy( void* combination, void* data ) {
	Union*	self = (Union*)combination;
	unsigned shape_I = 0;
    
   for( shape_I = 0 ; shape_I < self->shapeCount ; shape_I++ ) {
      Stg_Component_Destroy( self->shapeList[shape_I], data, False );
   }
	Memory_Free( self->isComplement );

	_Stg_Shape_Destroy( self, data );
}

/*--------------------------------------------------------------------------------------------------------------------------
** Public Functions
*/
	
/*---------------------------------------------------------------------------------------------------------------------
** Private Member functions
*/
Bool _Union_IsCoordInside( void* combination, Coord coord ) {
	Union*    self           = (Union*)combination;
	Index           shapeCount     = self->shapeCount;
	Index           shape_I;
	Bool            value;

	for( shape_I = 0 ; shape_I < shapeCount ; shape_I++ ) {
		value = Stg_Shape_IsCoordInside( self->shapeList[ shape_I ], coord );

		
		if ( self->isComplement[ shape_I ] )
			value = !value;
			
		if ( value )
			return True;
	}
	return False;
}	


double _Union_CalculateVolume( void* combination ) {
	assert ( 0 );
	return 0.0;
}

void _Union_DistanceFromCenterAxis( void* shape, Coord coord, double* disVec ) {
	Stg_Shape* self = (Stg_Shape*)shape;
	Journal_Firewall( False, Journal_Register( Error_Type, self->type ),
	"Error in function %s: This functions hasn't been implemented.", 
	"Please inform underworld-dev@vpac.org you've received this error.\n", __func__ );
}
	


