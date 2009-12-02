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
** $Id: DynamicVC.c 3881 2006-10-26 03:14:19Z KathleenHumble $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdarg.h>
#include "Base/Foundation/Foundation.h"
#include "Base/IO/IO.h"
#include "Base/Container/Container.h"
#include "Base/Automation/Automation.h"
#include "Base/Extensibility/Extensibility.h"

#include "types.h"
#include "shortcuts.h"
#include "Variable.h"
#include "Variable_Register.h"
#include "ConditionFunction.h"
#include "ConditionFunction_Register.h"
#include "VariableCondition.h"
#include "DynamicVC.h"

#include <string.h>
#include <assert.h>

const Type DynamicVC_Type = "DynamicVC";
const Name defaultDynamicVCName = "defaultDynamicVCName";

/*--------------------------------------------------------------------------------------------------------------------------
** Constructor
*/

VariableCondition* DynamicVC_Factory(
	AbstractContext*					context,
	Variable_Register*				varReg, 
	ConditionFunction_Register*	conFuncReg, 
	Dictionary*							dict, 
	void*									data )
{
	return (VariableCondition*)DynamicVC_New( defaultDynamicVCName, context, varReg, conFuncReg, dict );
}

DynamicVC* DynamicVC_New(
	Name									name,
	AbstractContext*					context,
	Variable_Register*				variable_Register, 
	ConditionFunction_Register*	conFunc_Register, 
	Dictionary*							dictionary )
{
	DynamicVC* self = _DynamicVC_DefaultNew( name );

	self->isConstructed = True;
	_VariableCondition_Init( self, context, variable_Register, conFunc_Register, dictionary );
	_DynamicVC_Init( self );	

	return self;
}

DynamicVC* _DynamicVC_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                              _sizeOfSelf = sizeof(DynamicVC);
	Type                                                      type = DynamicVC_Type;
	Stg_Class_DeleteFunction*                              _delete = _DynamicVC_Delete;
	Stg_Class_PrintFunction*                                _print = _DynamicVC_Print;
	Stg_Class_CopyFunction*                                  _copy = _DynamicVC_Copy;
	Stg_Component_DefaultConstructorFunction*  _defaultConstructor = (Stg_Component_DefaultConstructorFunction*)_DynamicVC_DefaultNew;
	Stg_Component_ConstructFunction*                    _construct = _VariableCondition_AssignFromXML;
	Stg_Component_BuildFunction*                            _build = _VariableCondition_Build;
	Stg_Component_InitialiseFunction*                  _initialise = _VariableCondition_Initialise;
	Stg_Component_ExecuteFunction*                        _execute = _VariableCondition_Execute;
	Stg_Component_DestroyFunction*                        _destroy = _VariableCondition_Destroy;
	AllocationType                              nameAllocationType = NON_GLOBAL;
	VariableCondition_BuildSelfFunc*                    _buildSelf = NULL;
	VariableCondition_PrintConciseFunc*              _printConcise = _DynamicVC_PrintConcise;
	VariableCondition_ReadDictionaryFunc*          _readDictionary = _DynamicVC_ReadDictionary;
	VariableCondition_GetSetFunc*                          _getSet = _DynamicVC_GetSet;
	VariableCondition_GetVariableCountFunc*      _getVariableCount = _DynamicVC_GetVariableCount;
	VariableCondition_GetVariableIndexFunc*      _getVariableIndex = _DynamicVC_GetVariableIndex;
	VariableCondition_GetValueIndexFunc*            _getValueIndex = _DynamicVC_GetValueIndex;
	VariableCondition_GetValueCountFunc*            _getValueCount = _DynamicVC_GetValueCount;
	VariableCondition_GetValueFunc*                      _getValue = _DynamicVC_GetValue;
	VariableCondition_ApplyFunc*                            _apply = DynamicVC_Apply;

	return (DynamicVC*)_DynamicVC_New(  DYNAMICVC_PASSARGS  );
}

DynamicVC* _DynamicVC_New(  DYNAMICVC_DEFARGS  ) {
	DynamicVC* self;

	/* Allocate memory/General info */
	assert( _sizeOfSelf >= sizeof(DynamicVC) );
	self = (DynamicVC*)_VariableCondition_New(  VARIABLECONDITION_PASSARGS  );
	
	/* Virtual info */
	
	/* Stg_Class info */
	
	return self;
}

void _DynamicVC_Init( void* vc ) {
	DynamicVC* self = (DynamicVC*)vc;

	assert( self );
	self->vcMap = IMap_New();
}	

/*--------------------------------------------------------------------------------------------------------------------------
** General virtual functions
*/

void _DynamicVC_AssignFromXML( void* vc, Stg_ComponentFactory* cf, void* data ) {
}

void _DynamicVC_Build( void* vc, void* data ) {
}

void _DynamicVC_Initialise( void* vc, void* data ) {
}

void _DynamicVC_Execute( void* vc, void* data ) {
}

void _DynamicVC_Destroy( void* vc, void* data ) {
	DynamicVC* self = (DynamicVC*)vc;

	NewClass_Delete( self->vcMap );

	__VariableCondition( self, data );
}

void _DynamicVC_ReadDictionary( void* vc, void* dict ) {
}

void _DynamicVC_Delete( void* vc ) {
	DynamicVC* self = (DynamicVC*)vc;

	/* Stg_Class_Delete parent */
	_VariableCondition_Delete( self );
}

void _DynamicVC_Print( void* vc, Stream* stream ) {
	_VariableCondition_Print( vc );
}

void* _DynamicVC_Copy( void* vc, void* dest, Bool deep, Name nameExt, struct PtrMap* ptrMap ) {
	return NULL;
}


/*--------------------------------------------------------------------------------------------------------------------------
** Macros
*/


/*--------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

IndexSet* _DynamicVC_GetSet( void* vc ) {
	return NULL;
}

VariableCondition_VariableIndex _DynamicVC_GetVariableCount( void* vc, Index globalIndex ) {
	return 0;
}

Variable_Index _DynamicVC_GetVariableIndex( void* vc, Index globalIndex, VariableCondition_VariableIndex varIndex ) {
	return 0;
}

VariableCondition_ValueIndex _DynamicVC_GetValueIndex( void* vc, 
						       Index globalIndex, 
						       VariableCondition_VariableIndex varIndex )
{
	return 0;
}

VariableCondition_ValueIndex _DynamicVC_GetValueCount( void* vc ) {
	return 0;
}

VariableCondition_Value _DynamicVC_GetValue( void* vc, VariableCondition_ValueIndex valIndex ) {
	VariableCondition_Value val;

	val.type = VC_ValueType_Double;
	val.as.typeDouble = 0.0;
	return val;
}

void _DynamicVC_PrintConcise( void* vc, Stream* stream ) {
}

void DynamicVC_Apply( void* vc, void* ctx ) {
	DynamicVC* self = (DynamicVC*)vc;
	IMapIter* iter;
	int ind, valInd, varInd;
	VariableCondition_Value* val;
	ConditionFunction* cf;
	Stream *errorStrm = Journal_Register( Error_Type, self->type );

	varInd = Variable_Register_GetIndex( self->variable_Register, self->var->name );
	iter = IMapIter_New();
	for( IMap_First( self->vcMap, iter );
	     Iter_IsValid( iter ); 
	     IMapIter_Next( iter ) )
	{
		ind = IMapIter_GetKey( iter );
		valInd = IMapIter_GetValue( iter );
		val = self->valueTbl + valInd;
		switch( val->type ) {
			case VC_ValueType_Double:
				Journal_Firewall( self->var->dataTypeCounts[0] == 1, errorStrm,
					"Error - in %s: while applying values for variable condition "
					"\"%s\", to index %d - asked to apply a scalar %s to Variable \"%s\" "
					"which has %d components. Specify a scalar Variable instead.\n",
					__func__, self->name, ind, "double",
					self->var->name, self->var->dataTypeCounts[0] );
				Variable_SetValueDouble( self->var, 
							 ind, 
							 val->as.typeDouble );
				break;

			case VC_ValueType_DoubleArray:
				Variable_SetValue( self->var, 
						   ind, 
						   val->as.typeArray.array );
				break;

			case VC_ValueType_CFIndex:
				Journal_Firewall( val->as.typeCFIndex != (unsigned)-1, errorStrm,
					"Error - in %s: trying to apply to index %d of variable \"%s\", which "
					"is supposed to be a condition function, but the cond. func. wasn't "
					"found in the c.f. register.\n", __func__, ind, self->var->name );
				cf = self->conFunc_Register->_cf[val->as.typeCFIndex];
				ConditionFunction_Apply( cf, 
							 ind, 
							 varInd, 
							 ctx, 
							 Variable_GetStructPtr( self->var, ind ) );
				break;

			default:
				assert( 0 );
				break;
		}
	}
	NewClass_Delete( iter );
}


/*--------------------------------------------------------------------------------------------------------------------------
** Build functions
*/


/*--------------------------------------------------------------------------------------------------------------------------
** Functions
*/

void DynamicVC_SetVariable( void* vc, Variable* var ) {
	DynamicVC* self = (DynamicVC*)vc;

	assert( self );
	self->var = var;
}

void DynamicVC_SetMaxEntries( void* vc, int maxEntries ) {
	DynamicVC* self = (DynamicVC*)vc;

	assert( self );
	/*IMap_Clear( self->vcMap );*/
	IMap_SetMaxSize( self->vcMap, maxEntries );
}

int DynamicVC_GetMaxEntries( void* vc ) {
   DynamicVC* self = (DynamicVC*)vc;

   return IMap_GetMaxSize( self->vcMap );
}

void DynamicVC_SetValues( void* vc, int nVals, VariableCondition_Value* vals ) {
	DynamicVC* self = (DynamicVC*)vc;

	assert( self );
	self->valueCount = nVals;
	self->valueTbl = MemRearray( self->valueTbl, VariableCondition_Value, nVals, DynamicVC_Type );
	memcpy( self->valueTbl, vals, nVals * sizeof(VariableCondition_Value) );
}

void DynamicVC_Insert( void* vc, int index, int valIndex ) {
	DynamicVC* self = (DynamicVC*)vc;

	assert( self );
	if( !IMap_Has( self->vcMap, index ) )
		IMap_Insert( self->vcMap, index, valIndex );
}

void DynamicVC_Remove( void* vc, int index ) {
	DynamicVC* self = (DynamicVC*)vc;

	assert( self );
	if( IMap_Has( self->vcMap, index ) )
		IMap_Remove( self->vcMap, index );
}

Bool DynamicVC_Has( void* _self, int index ) {
	DynamicVC* self = (DynamicVC*)_self;

	assert( self );
        return IMap_Has( self->vcMap, index );
}


