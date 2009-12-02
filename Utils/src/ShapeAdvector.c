/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street, Melbourne, 3053, Australia.
**
** Authors:
**	Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**	Patrick D. Sunter, Software Engineer, VPAC. (pds@vpac.org)
**	Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**	Siew-Ching Tan, Software Engineer, VPAC. (siew@vpac.org) ) {
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
** $Id: ShapeAdvector.c 212 2005-11-08 23:50:02Z RobertTurnbull $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>

#include <StgDomain/Geometry/Geometry.h>
#include <StgDomain/Shape/Shape.h>
#include <StgDomain/Mesh/Mesh.h>

#include "types.h"
#include "DomainContext.h"
#include "ShapeAdvector.h"
#include "TimeIntegrator.h"
#include "TimeIntegratee.h"
#include "FieldVariable.h"

#include <assert.h>
#include <string.h>
#include <math.h>


/* Textual name of this class */
const Type ShapeAdvector_Type = "ShapeAdvector";

/*-------------------------------------------------------------------------------------------------------------------------
** Constructors
*/
ShapeAdvector* ShapeAdvector_New(
	Name					name,
	DomainContext*		context,
	TimeIntegrator*	timeIntegrator,
	FieldVariable*		velocityField,
	Stg_Shape*			shape,
	Bool					allowFallbackToFirstOrder )
{
	ShapeAdvector* self = (ShapeAdvector*) _ShapeAdvector_DefaultNew( name );

	self->isConstructed = True;
	_ShapeAdvector_Init( self, context, timeIntegrator, velocityField, shape, allowFallbackToFirstOrder );

	return self;
}

void* _ShapeAdvector_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                              _sizeOfSelf = sizeof(ShapeAdvector);
	Type                                                      type = ShapeAdvector_Type;
	Stg_Class_DeleteFunction*                              _delete = _ShapeAdvector_Delete;
	Stg_Class_PrintFunction*                                _print = _ShapeAdvector_Print;
	Stg_Class_CopyFunction*                                  _copy = _ShapeAdvector_Copy;
	Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _ShapeAdvector_DefaultNew;
	Stg_Component_ConstructFunction*                    _construct = _ShapeAdvector_AssignFromXML;
	Stg_Component_BuildFunction*                            _build = _ShapeAdvector_Build;
	Stg_Component_InitialiseFunction*                  _initialise = _ShapeAdvector_Initialise;
	Stg_Component_ExecuteFunction*                        _execute = _ShapeAdvector_Execute;
	Stg_Component_DestroyFunction*                        _destroy = _ShapeAdvector_Destroy;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = ZERO;

	return (void*) _ShapeAdvector_New(  SHAPEADVECTOR_PASSARGS  );
}

ShapeAdvector* _ShapeAdvector_New(  SHAPEADVECTOR_DEFARGS  )
{
	ShapeAdvector* self;
	
	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(ShapeAdvector) );
	/* The following terms are parameters that have been passed into this function but are being set before being passed onto the parent */
	/* This means that any values of these parameters that are passed into this function are not passed onto the parent function
	   and so should be set to ZERO in any children of this class. */
	nameAllocationType = NON_GLOBAL;

	self = (ShapeAdvector*)_Stg_Component_New(  STG_COMPONENT_PASSARGS  );
	
	/* General info */

	/* Virtual Info */
	
	return self;
}

void _ShapeAdvector_Init( 
	ShapeAdvector*		self,
	DomainContext*		context,
	TimeIntegrator*	timeIntegrator,
	FieldVariable*		velocityField,
	Stg_Shape*			shape,
	Bool					allowFallbackToFirstOrder )
{
	self->context = context;
	self->velocityField = velocityField;
	self->shape = shape;
	self->shapeCount = 1;
	self->shapeCentrePtr = shape->centre;

	self->shapeCentreVariable = Variable_NewVector( "shapeCentreVariable", (AbstractContext*)self->context,
		Variable_DataType_Double, shape->dim, &self->shapeCount, NULL,  &self->shapeCentrePtr, NULL );
	self->timeIntegratee = TimeIntegratee_New( "shapeTimeIntegratee", self->context, timeIntegrator, self->shapeCentreVariable, 1,
		(Stg_Component**) &velocityField, allowFallbackToFirstOrder );
}

/*------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _ShapeAdvector_Delete( void* shapeAdvector ) {
	ShapeAdvector* self = (ShapeAdvector*)shapeAdvector;

	/* Delete parent */
	_Stg_Component_Delete( self );
}


void _ShapeAdvector_Print( void* shapeAdvector, Stream* stream ) {
	ShapeAdvector* self = (ShapeAdvector*)shapeAdvector;
	
	/* Print parent */
	_Stg_Component_Print( self, stream );
}


void* _ShapeAdvector_Copy( void* shapeAdvector, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	ShapeAdvector*	self = (ShapeAdvector*)shapeAdvector;
	ShapeAdvector*	newShapeAdvector;
	
	newShapeAdvector = (ShapeAdvector*)_Stg_Component_Copy( self, dest, deep, nameExt, ptrMap );

	newShapeAdvector->velocityField = self->velocityField;
	newShapeAdvector->shape = self->shape;
	
	return (void*)newShapeAdvector;
}

void _ShapeAdvector_AssignFromXML( void* shapeAdvector, Stg_ComponentFactory* cf, void* data ) {
	ShapeAdvector*		self = (ShapeAdvector*) shapeAdvector;
	FieldVariable*		velocityField;
	Stg_Shape*			shape;
	TimeIntegrator*	timeIntegrator;
	Bool					allowFallbackToFirstOrder = False;
	DomainContext*		context;

	context = Stg_ComponentFactory_ConstructByKey( cf, self->name, "Context", DomainContext, False, data );
	if( !context )
		context = Stg_ComponentFactory_ConstructByName( cf, "context", DomainContext, True, data );

	timeIntegrator = Stg_ComponentFactory_ConstructByKey( cf, self->name, "TimeIntegrator", TimeIntegrator, True, data  ) ;
	velocityField = Stg_ComponentFactory_ConstructByKey( cf, self->name, "VelocityField", FieldVariable, True, data  ) ;
	shape = Stg_ComponentFactory_ConstructByKey( cf, self->name, "Shape", Stg_Shape, True, data ) ;
	allowFallbackToFirstOrder = Stg_ComponentFactory_GetBool( cf, self->name, "allowFallbackToFirstOrder", False );

	_ShapeAdvector_Init( self, context, timeIntegrator, velocityField, shape, allowFallbackToFirstOrder );
}

void _ShapeAdvector_Build( void* shapeAdvector, void* data ) {
}

void _ShapeAdvector_Initialise( void* shapeAdvector, void* data ) {
}

void _ShapeAdvector_Execute( void* shapeAdvector, void* data ) {
}

void _ShapeAdvector_Destroy( void* shapeAdvector, void* data ) {
	ShapeAdvector* self = (ShapeAdvector*)shapeAdvector;

	_Stg_Component_Delete( self->shapeCentreVariable );
	_Stg_Component_Delete( self->timeIntegratee );
}




