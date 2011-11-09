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
** $Id: EquationShape.c 4056 2007-03-29 04:55:51Z JulianGiordani $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/Geometry/Geometry.h>
#include <StgDomain/Shape/Shape.h>
#include <StgDomain/Mesh/Mesh.h>
#include <StgDomain/Utils/Utils.h>

#include "types.h"
#include "ShapeClass.h"
#include "EquationShape.h"


#include <assert.h>
#include <string.h>
#include <math.h>

/* Textual name of this class */
const Type EquationShape_Type = "EquationShape";

/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/

EquationShape* EquationShape_New(Name name,
                                 Dimension_Index dim,
                                 XYZ centre, 
                                 double alpha,
                                 double beta,
                                 double gamma,
                                 char* equation)
{ 
  EquationShape* self = (EquationShape*)_EquationShape_DefaultNew( name );
  
  _Stg_Shape_Init( self, dim, centre, False, alpha, beta, gamma);
  _EquationShape_Init( self, equation );
  return self;
}

EquationShape* _EquationShape_New(  EQUATIONSHAPE_DEFARGS  )
{
  EquationShape* self;
	
  /* Allocate memory */
  assert( _sizeOfSelf >= sizeof(EquationShape) );
  self = (EquationShape*)_Stg_Shape_New(  STG_SHAPE_PASSARGS  );
	
  /* General info */
  return self;
}

void _EquationShape_Init( void* equation, char *eqn ) {
  EquationShape* self = (EquationShape*)equation;
	
  self->equation = StG_Strdup(eqn);
}
	
/*------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _EquationShape_Delete( void* equation ) {
  EquationShape*       self = (EquationShape*)equation;
	
  /* Delete parent */
  _Stg_Shape_Delete( self );
}


void _EquationShape_Print( void* equation, Stream* stream ) {
  EquationShape* self = (EquationShape*)equation;
	
  /* Print parent */
  _Stg_Shape_Print( self, stream );
}



void* _EquationShape_Copy( const void* equation, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
  EquationShape*	self = (EquationShape*)equation;
  EquationShape*	newEquationShape;
	
  newEquationShape = (EquationShape*)_Stg_Shape_Copy( self, dest, deep, nameExt, ptrMap );

  newEquationShape->equation = StG_Strdup(self->equation);
  return (void*)newEquationShape;
}

void* _EquationShape_DefaultNew( Name name ) {
  /* Variables set in this function */
  SizeT                                                  _sizeOfSelf = sizeof(EquationShape);
  Type                                                          type = EquationShape_Type;
  Stg_Class_DeleteFunction*                                  _delete = _EquationShape_Delete;
  Stg_Class_PrintFunction*                                    _print = _EquationShape_Print;
  Stg_Class_CopyFunction*                                      _copy = _EquationShape_Copy;
  Stg_Component_DefaultConstructorFunction*      _defaultConstructor = _EquationShape_DefaultNew;
  Stg_Component_ConstructFunction*                        _construct = _EquationShape_AssignFromXML;
  Stg_Component_BuildFunction*                                _build = _EquationShape_Build;
  Stg_Component_InitialiseFunction*                      _initialise = _EquationShape_Initialise;
  Stg_Component_ExecuteFunction*                            _execute = _EquationShape_Execute;
  Stg_Component_DestroyFunction*                            _destroy = _EquationShape_Destroy;
  Stg_Shape_IsCoordInsideFunction*                    _isCoordInside = _EquationShape_IsCoordInside;
  Stg_Shape_CalculateVolumeFunction*                _calculateVolume = _EquationShape_CalculateVolume;
  Stg_Shape_DistanceFromCenterAxisFunction*  _distanceFromCenterAxis = _EquationShape_DistanceFromCenterAxis;

  /* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
  AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

  return (void*) _EquationShape_New(  EQUATIONSHAPE_PASSARGS  );
}


void _EquationShape_AssignFromXML( void* equation, Stg_ComponentFactory* cf, void* data ) {
  EquationShape*           self       = (EquationShape*)equation;
	
  _Stg_Shape_AssignFromXML( self, cf, data );

  char *eqn = Stg_ComponentFactory_GetString( cf, self->name, (Dictionary_Entry_Key)"equation", "");
  _EquationShape_Init( self, eqn );
}

void _EquationShape_Build( void* equation, void* data ) {
  EquationShape*	self = (EquationShape*)equation;

  _Stg_Shape_Build( self, data );
}

void _EquationShape_Initialise( void* equation, void* data ) {
  EquationShape*	self = (EquationShape*)equation;
	
  _Stg_Shape_Initialise( self, data );
}

void _EquationShape_Execute( void* equation, void* data ) {
  EquationShape*	self = (EquationShape*)equation;
	
  _Stg_Shape_Execute( self, data );
}

void _EquationShape_Destroy( void* equation, void* data ) {
  EquationShape*	self = (EquationShape*)equation;

  Memory_Free(self->equation);
  _Stg_Shape_Destroy( self, data );
}

/*--------------------------------------------------------------------------------------------------------------------------
** Public Functions
*/
	
/*---------------------------------------------------------------------------------------------------------------------
** Private Member functions
*/

Bool _EquationShape_IsCoordInside( void* equation, const Coord coord ) {
  EquationShape* self=(EquationShape*) equation;
  Coord newCoord;

  /* Transform coordinate into canonical reference frame */
  Stg_Shape_TransformCoord( self, coord, newCoord );

  return Equation_eval(newCoord,(DomainContext*)(self->context),
                       self->equation)>=0 ? True : False;
}


double _EquationShape_CalculateVolume( void* equation ) {
  assert( 0 );
  return 0.0;
}

void _EquationShape_DistanceFromCenterAxis( void* shape, const Coord coord, double* disVec ){
  Stg_Shape* self = (Stg_Shape*)shape;
  Journal_Firewall( False, Journal_Register( Error_Type, (Name)self->type  ),
                    "Error in function %s: This functions hasn't been implemented.", 
                    "Please inform cig-long@geodynamics.org that you have received this error.\n", __func__ );
}



