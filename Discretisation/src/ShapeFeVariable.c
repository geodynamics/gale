/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003-2006, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street,
**	Melbourne, 3053, Australia.
**
** Primary Contributing Organisations:
**	Victorian Partnership for Advanced Computing Ltd, Computational Software Development - http://csd.vpac.org
**	Australian Computational Earth Systems Simulator - http://www.access.edu.au
**	Monash Cluster Computing - http://www.mcc.monash.edu.au
**	Computational Infrastructure for Geodynamics - http://www.geodynamics.org
**
** Contributors:
**	Patrick D. Sunter, Software Engineer, VPAC. (pds@vpac.org)
**	Robert Turnbull, Research Assistant, Monash University. (robert.turnbull@sci.monash.edu.au)
**	Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**	David May, PhD Student, Monash University (david.may@sci.monash.edu.au)
**	Louis Moresi, Associate Professor, Monash University. (louis.moresi@sci.monash.edu.au)
**	Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**	Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
**	Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
**	Julian Giordani, Research Assistant, Monash University. (julian.giordani@sci.monash.edu.au)
**	Vincent Lemiale, Postdoctoral Fellow, Monash University. (vincent.lemiale@sci.monash.edu.au)
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
** $Id: ShapeFeVariable.c 532 2006-04-04 00:21:59Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <string.h>
#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>

#include "types.h"
#include "FeMesh.h"
#include "FeVariable.h"
#include "ShapeFeVariable.h"

#include <assert.h>

const Type ShapeFeVariable_Type = "ShapeFeVariable";

void* ShapeFeVariable_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                                       _sizeOfSelf = sizeof(ShapeFeVariable);
	Type                                                               type = ShapeFeVariable_Type;
	Stg_Class_DeleteFunction*                                       _delete = _ShapeFeVariable_Delete;
	Stg_Class_PrintFunction*                                         _print = _ShapeFeVariable_Print;
	Stg_Class_CopyFunction*                                           _copy = _ShapeFeVariable_Copy;
	Stg_Component_DefaultConstructorFunction*           _defaultConstructor = (Stg_Component_DefaultConstructorFunction*)ShapeFeVariable_DefaultNew;
	Stg_Component_ConstructFunction*                             _construct = _ShapeFeVariable_AssignFromXML;
	Stg_Component_BuildFunction*                                     _build = _ShapeFeVariable_Build;
	Stg_Component_InitialiseFunction*                           _initialise = _ShapeFeVariable_Initialise;
	Stg_Component_ExecuteFunction*                                 _execute = _ShapeFeVariable_Execute;
	Stg_Component_DestroyFunction*                                 _destroy = _ShapeFeVariable_Destroy;
	FieldVariable_InterpolateValueAtFunction*           _interpolateValueAt = _FeVariable_InterpolateValueAt;
	FieldVariable_GetValueFunction*             _getMinGlobalFieldMagnitude = _FeVariable_GetMinGlobalFieldMagnitude;
	FieldVariable_GetValueFunction*             _getMaxGlobalFieldMagnitude = _FeVariable_GetMaxGlobalFieldMagnitude;
	FieldVariable_GetCoordFunction*                _getMinAndMaxLocalCoords = _FeVariable_GetMinAndMaxLocalCoords;
	FieldVariable_GetCoordFunction*               _getMinAndMaxGlobalCoords = _FeVariable_GetMinAndMaxGlobalCoords;
	FeVariable_InterpolateWithinElementFunction*  _interpolateWithinElement = _FeVariable_InterpolateNodeValuesToElLocalCoord;
	FeVariable_GetValueAtNodeFunction*                      _getValueAtNode = _FeVariable_GetValueAtNode;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType                             nameAllocationType = ZERO;
	FeVariable_SyncShadowValuesFunc*            _syncShadowValues = ZERO;

	return (ShapeFeVariable*) _ShapeFeVariable_New(  SHAPEFEVARIABLE_PASSARGS  );
}			

ShapeFeVariable* _ShapeFeVariable_New(  SHAPEFEVARIABLE_DEFARGS  ) {
	ShapeFeVariable* self;
	
	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(ShapeFeVariable) );
	/* The following terms are parameters that have been passed into this function but are being set before being passed onto the parent */
	/* This means that any values of these parameters that are passed into this function are not passed onto the parent function
	   and so should be set to ZERO in any children of this class. */
	nameAllocationType = NON_GLOBAL;
	_syncShadowValues  = _FeVariable_SyncShadowValues;

	self = (ShapeFeVariable*) _FeVariable_New(  FEVARIABLE_PASSARGS  );

	return self;
}

void _ShapeFeVariable_Init( void* shapeFeVariable, Stg_Shape* shape ) {
	ShapeFeVariable* self = (ShapeFeVariable*) shapeFeVariable;

	self->shape = shape;
	
	/* EP_AppendClassHook( Context_GetEntryPoint( context, AbstractContext_EP_UpdateClass ),	ParticleFeVariable_Update, self ); */
}

void _ShapeFeVariable_Delete( void* _shapeFeVariable ) {
	ShapeFeVariable* self = (ShapeFeVariable*) _shapeFeVariable;

	_FeVariable_Delete( self );
}

void _ShapeFeVariable_Print( void* _shapeFeVariable, Stream* stream ) {
	ShapeFeVariable* self = (ShapeFeVariable*) _shapeFeVariable;

	/* Print parent */
	_FeVariable_Print( self, stream );

	Journal_PrintPointer( stream, self->shape );
}

void* _ShapeFeVariable_Copy( void* shapeFeVariable, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	/* ShapeFeVariable*	self = (ShapeFeVariable*)shapeFeVariable; */
	ShapeFeVariable* newShapeFeVariable;
	
	assert(0);
	return (void*)newShapeFeVariable;
}

void _ShapeFeVariable_AssignFromXML( void* shapeFeVariable, Stg_ComponentFactory* cf, void* data ) {
	ShapeFeVariable* self = (ShapeFeVariable*) shapeFeVariable;

	_FeVariable_AssignFromXML( self, cf, data );

	_ShapeFeVariable_Init( self, Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"Shape", Stg_Shape, True, data )   ) ;
}

void _ShapeFeVariable_Build( void* shapeFeVariable, void* data ) {
	ShapeFeVariable* self = (ShapeFeVariable*) shapeFeVariable;

	Stg_Component_Build( self->shape, data, False );
	_FeVariable_Build( self, data );
}

void _ShapeFeVariable_Initialise( void* shapeFeVariable, void* data ) {
	ShapeFeVariable* self = (ShapeFeVariable*) shapeFeVariable;
	Node_DomainIndex node_dI = 0;

	Stg_Component_Initialise( self->shape, data, False );
	_FeVariable_Initialise( self, data );

	/* Set up the basic "level set" describing if nodes are inside the shape or not */
	for ( node_dI = 0; node_dI < Mesh_GetDomainSize( self->feMesh, MT_VERTEX ); node_dI++ ) {
		if ( True == Stg_Shape_IsCoordInside( self->shape, Mesh_GetVertex( self->feMesh, node_dI ) ) ) {
			/* set value = 1 */
			FeVariable_SetComponentAtNode( self, node_dI, 0, 1 );
		}		
		else {
			/* set value = 0 */
			FeVariable_SetComponentAtNode( self, node_dI, 0, 0 );
		}
	}
}

void _ShapeFeVariable_Execute( void* shapeFeVariable, void* data ) {
	ShapeFeVariable* self = (ShapeFeVariable*) shapeFeVariable;
	
	_FeVariable_Execute( self, data );
}

void _ShapeFeVariable_Destroy( void* shapeFeVariable, void* data ) {
	ShapeFeVariable* self = (ShapeFeVariable*) shapeFeVariable;
	
	_FeVariable_Destroy( self, data );
}


