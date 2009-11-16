/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street, Melbourne, 3053, Australia.
**
** Authors:
**	Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**	Patrick D. Sunter, Software Engineer, VPAC. (pds@vpac.org)
**	Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**	Siew-Ching Tan, Software Engineer, VPAC. (siew@vpac.org)
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
** $Id: OperatorFieldVariable.c 4137 2007-06-07 05:46:46Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <string.h>
#include <mpi.h>
#include <StGermain/StGermain.h>

#include <StgDomain/Geometry/Geometry.h>
#include <StgDomain/Shape/Shape.h>
#include <StgDomain/Mesh/Mesh.h>

#include "types.h"
#include "FieldVariable.h"
#include "Operator.h"
#include "OperatorFieldVariable.h"
#include "DomainContext.h"
#include "FieldVariable_Register.h"

#include <assert.h>

const Type OperatorFieldVariable_Type = "OperatorFieldVariable";
const Name defaultOperatorFieldVariableName = "defaultOperatorFieldVariableName";

OperatorFieldVariable* OperatorFieldVariable_NewUnary( Name	name, DomainContext* context, void* _fieldVariable, Name operatorName ) {
	FieldVariable* fieldVariable = (FieldVariable*) _fieldVariable;
       	
	return OperatorFieldVariable_New( 
		name,
		context,
		OperatorFieldVariable_UnaryInterpolationFunc, 
		operatorName,
		1,
		&fieldVariable, 
		fieldVariable->dim,
		fieldVariable->isCheckpointedAndReloaded,
		fieldVariable->communicator,
		fieldVariable->fieldVariable_Register );
}

OperatorFieldVariable* OperatorFieldVariable_NewBinary( Name name, DomainContext* context, void* _fieldVariable1, void* _fieldVariable2, Name operatorName ) {
	FieldVariable* fieldVariableList[2];
       
	fieldVariableList[0] = (FieldVariable*) _fieldVariable1;
	fieldVariableList[1] = (FieldVariable*) _fieldVariable2;
	
	return OperatorFieldVariable_New( 
		name,
		context,
		OperatorFieldVariable_BinaryInterpolationFunc, 
		operatorName,
		2, 
		fieldVariableList, 
		fieldVariableList[0]->dim,
		fieldVariableList[0]->isCheckpointedAndReloaded,
		fieldVariableList[0]->communicator,
		fieldVariableList[0]->fieldVariable_Register );
}

OperatorFieldVariable* OperatorFieldVariable_New( 
	Name													name,
	DomainContext*										context,
	FieldVariable_InterpolateValueAtFunction*	interpolateValueAt,
	Name													operatorName,
	Index													fieldVariableCount,
	FieldVariable**									fieldVariableList,
	Dimension_Index									dim,
	Bool													isCheckpointedAndReloaded,
	MPI_Comm												communicator,
	FieldVariable_Register*							fieldVariable_Register ) 
{
	OperatorFieldVariable* self = _OperatorFieldVariable_DefaultNew( name );

	self->isConstructed = True;
	_FieldVariable_Init( self, context, fieldVariableCount, dim, isCheckpointedAndReloaded, communicator, fieldVariable_Register );
	_OperatorFieldVariable_Init( self, operatorName, fieldVariableCount, fieldVariableList  );

	return self;
}

OperatorFieldVariable* _OperatorFieldVariable_DefaultNew( Name name ) {
	return _OperatorFieldVariable_New( 
		sizeof(OperatorFieldVariable), 
		OperatorFieldVariable_Type, 
		_FieldVariable_Delete, 
		_OperatorFieldVariable_Print,
		_OperatorFieldVariable_Copy, 
		(Stg_Component_DefaultConstructorFunction*)_OperatorFieldVariable_DefaultNew,
		_OperatorFieldVariable_AssignFromXML,
		_OperatorFieldVariable_Build, 
		_OperatorFieldVariable_Initialise, 
		_OperatorFieldVariable_Execute,
		_OperatorFieldVariable_Destroy,
		name,
		NON_GLOBAL, 
		_OperatorFieldVariable_InterpolateValueAt,
		_OperatorFieldVariable_GetMinLocalFieldMagnitude,
		_OperatorFieldVariable_GetMaxLocalFieldMagnitude, 
		_OperatorFieldVariable_GetMinAndMaxLocalCoords,
		_OperatorFieldVariable_GetMinAndMaxGlobalCoords, 
		0,
		0,
		False, /* Setting default to false, as do not want to checkpoint OperatorFeVariables */
		MPI_COMM_WORLD,
		NULL,
		NULL,
		NULL);
}

OperatorFieldVariable* _OperatorFieldVariable_New( OPERATORFIELDVARIABLE_DEFARGS ) {
	OperatorFieldVariable* self;
	
	/* Allocate memory */
	assert( sizeOfSelf >= sizeof(OperatorFieldVariable) );
	self = (OperatorFieldVariable*) _FieldVariable_New( FIELDVARIABLE_PASSARGS );

	return self;
} 

void _OperatorFieldVariable_Delete( void* _fieldVariable ) {
	OperatorFieldVariable* self = (OperatorFieldVariable*) _fieldVariable;

	_FieldVariable_Delete( self );
}

void _OperatorFieldVariable_Print( void* _fieldVariable, Stream* stream ) {
	OperatorFieldVariable* self = (OperatorFieldVariable*) _fieldVariable;
	Index                  fieldVariable_I;

	_FieldVariable_Print( self, stream );

	Journal_PrintValue( stream, self->fieldVariableCount );
	for ( fieldVariable_I = 0 ; fieldVariable_I < self->fieldVariableCount ; fieldVariable_I++ ) 
		Journal_Printf( stream, "\tFieldVariable %u - '%s'\n", fieldVariable_I, self->fieldVariableList[ fieldVariable_I ]->name );

}

void _OperatorFieldVariable_Init( void* ofv, Name operatorName, Index fieldVariableCount, FieldVariable** fieldVariableList ) {
	OperatorFieldVariable*	self = (OperatorFieldVariable*)ofv;
	FieldVariable*              fieldVariable;
	Index                       fieldVariable_I;
	Stream*                     errorStream       = Journal_Register( Error_Type, self->type );

	/* Create operator */
	self->_operator = Operator_NewFromName( operatorName, fieldVariableList[0]->fieldComponentCount, self->dim );
	self->fieldComponentCount = self->_operator->resultDofs; /* Reset this value from the one generated from the operator */

	self->fieldVariableCount     = fieldVariableCount;

	/* Copy field variable list */
	self->fieldVariableList      = Memory_Alloc_Array( FieldVariable*, fieldVariableCount, "Array of Field Variables" );
	memcpy( self->fieldVariableList, fieldVariableList, fieldVariableCount * sizeof( FieldVariable* ) );

	for ( fieldVariable_I = 0 ; fieldVariable_I < fieldVariableCount ; fieldVariable_I++ ) {
		fieldVariable = fieldVariableList[ fieldVariable_I ];
		Journal_Firewall( fieldVariable != NULL, errorStream, 
				"In func %s: FieldVariable %u in list is NULL\n", __func__, fieldVariable_I );
		Journal_Firewall( fieldVariable->fieldComponentCount <= MAX_FIELD_COMPONENTS, errorStream, 
			"In func %s: Field Variable '%s' has too many components.\n", __func__, fieldVariable->name );
	}
}


void* _OperatorFieldVariable_Copy( void* fieldVariable, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	OperatorFieldVariable*	self = (OperatorFieldVariable*)fieldVariable;
	OperatorFieldVariable*	newOperatorFieldVariable;
	
	newOperatorFieldVariable = _FieldVariable_Copy( self, dest, deep, nameExt, ptrMap );
	
	newOperatorFieldVariable->_operator              = self->_operator;
	newOperatorFieldVariable->fieldVariableCount     = self->fieldVariableCount;
	
	if (deep) {
		newOperatorFieldVariable->fieldVariableList = Memory_Alloc_Array( FieldVariable*, self->fieldVariableCount, 
				"Array of Field Variables" );
		memcpy( newOperatorFieldVariable->fieldVariableList, self->fieldVariableList, 
				self->fieldVariableCount * sizeof( FieldVariable* ) );
	}
	else 
		newOperatorFieldVariable->fieldVariableList = self->fieldVariableList;
	
	return (void*)newOperatorFieldVariable;
}

void _OperatorFieldVariable_AssignFromXML( void* fieldVariable, Stg_ComponentFactory* cf, void* data ) {
	OperatorFieldVariable*	self = (OperatorFieldVariable*) fieldVariable;
	Index							fieldVariableCount = 0;
	Name							operatorName;
	FieldVariable**			fieldVariableList;
	
	/* Construct Parent */
	_FieldVariable_AssignFromXML( self, cf, data );

	operatorName = Stg_ComponentFactory_GetString( cf, self->name, "Operator", "" );

	fieldVariableList = Stg_ComponentFactory_ConstructByList( 
		cf, 
		self->name, 
		"FieldVariables", 
		Stg_ComponentFactory_Unlimited, 
		FieldVariable, 
		True, 
		&fieldVariableCount,
		data );

	_OperatorFieldVariable_Init( self, operatorName, fieldVariableCount, fieldVariableList );

	Memory_Free( fieldVariableList );
}

void _OperatorFieldVariable_Build( void* fieldVariable, void* data ) {
	OperatorFieldVariable* self = (OperatorFieldVariable*) fieldVariable;
	Index                  fieldVariable_I;

	for ( fieldVariable_I = 0 ; fieldVariable_I < self->fieldVariableCount ; fieldVariable_I++ ) 
		Stg_Component_Build( self->fieldVariableList[ fieldVariable_I ] , data, False );
}

void _OperatorFieldVariable_Execute( void* fieldVariable, void* data ) {}

void _OperatorFieldVariable_Destroy( void* fieldVariable, void* data ) {
	OperatorFieldVariable* self = (OperatorFieldVariable*) fieldVariable;

	Memory_Free( self->fieldVariableList );
	
	_FieldVariable_Destroy( self, data );
}

void _OperatorFieldVariable_Initialise( void* fieldVariable, void* data ) {
	OperatorFieldVariable* self = (OperatorFieldVariable*) fieldVariable;
	Index                  fieldVariable_I;

	for ( fieldVariable_I = 0 ; fieldVariable_I < self->fieldVariableCount ; fieldVariable_I++ ) 
		Stg_Component_Initialise( self->fieldVariableList[ fieldVariable_I ] , data, False );
}

/* TODO - Think of something clever for these */
double _OperatorFieldVariable_GetMinLocalFieldMagnitude( void* fieldVariable ) { return 0.0; }

double _OperatorFieldVariable_GetMaxLocalFieldMagnitude( void* fieldVariable ) { return 0.0; }

void  _OperatorFieldVariable_GetMinAndMaxLocalCoords( void* fieldVariable, Coord min, Coord max ) {
	OperatorFieldVariable* self            = (OperatorFieldVariable*) fieldVariable;

	FieldVariable_GetMinAndMaxLocalCoords( self->fieldVariableList[0], min, max );
}

void  _OperatorFieldVariable_GetMinAndMaxGlobalCoords( void* fieldVariable, Coord min, Coord max ) {
	OperatorFieldVariable* self            = (OperatorFieldVariable*) fieldVariable;

	FieldVariable_GetMinAndMaxGlobalCoords( self->fieldVariableList[0], min, max );
}

InterpolationResult _OperatorFieldVariable_InterpolateValueAt( void* fieldVariable, Coord coord, double* value ) {
	OperatorFieldVariable* self            = (OperatorFieldVariable*) fieldVariable;

	switch ( self->fieldVariableCount ) {
		case 1:
			self->_interpolateValueAt = OperatorFieldVariable_UnaryInterpolationFunc; break;
		case 2:
			self->_interpolateValueAt = OperatorFieldVariable_BinaryInterpolationFunc; break;
		default:
			Journal_Firewall( False, Journal_Register( Error_Type, self->type ),
					"Can't use func '%s' with fieldVariableCount = %d\n", __func__, self->fieldVariableCount );
	}

	return FieldVariable_InterpolateValueAt( self, coord, value );
}

InterpolationResult OperatorFieldVariable_UnaryInterpolationFunc( void* fieldVariable, Coord coord, double* value ) {
	OperatorFieldVariable* self            = (OperatorFieldVariable*) fieldVariable;
	FieldVariable*         field0          = self->fieldVariableList[0];
	InterpolationResult    result;
	double                 fieldValue[ MAX_FIELD_COMPONENTS ]; 

	result = field0->_interpolateValueAt( field0, coord, fieldValue );
	Operator_CarryOutUnaryOperation( self->_operator, fieldValue, value );
	return result;
}

InterpolationResult OperatorFieldVariable_BinaryInterpolationFunc( void* fieldVariable, Coord coord, double* value ) {
	OperatorFieldVariable* self            = (OperatorFieldVariable*) fieldVariable;
	FieldVariable*         field0          = self->fieldVariableList[0];
	FieldVariable*         field1          = self->fieldVariableList[1];
	double                 fieldValue0[ MAX_FIELD_COMPONENTS ]; 
	double                 fieldValue1[ MAX_FIELD_COMPONENTS ]; 
	InterpolationResult    result0;
	InterpolationResult    result1;

	result0 = field0->_interpolateValueAt( field0, coord, fieldValue0 );
	result1 = field1->_interpolateValueAt( field1, coord, fieldValue1 );

	Operator_CarryOutBinaryOperation( self->_operator, fieldValue0, fieldValue1, value ); 

	return result0;
}
