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

#include "types.h"
#include "shortcuts.h"
#include "Stg_Component.h"
#include "Stg_ComponentFactory.h"
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

VariableCondition* DynamicVC_Factory( Variable_Register* varReg, 
				      ConditionFunction_Register* conFuncReg, 
				      Dictionary* dict, 
				      void* data )
{
	return (VariableCondition*)DynamicVC_New( defaultDynamicVCName, varReg, conFuncReg, dict );
}

DynamicVC* DynamicVC_New( Name name,
			  Variable_Register* varReg, 
			  ConditionFunction_Register* conFuncReg, 
			  Dictionary* dict )
{
	return _DynamicVC_New( sizeof(DynamicVC), 
			       DynamicVC_Type, 
			       _DynamicVC_Delete, 
			       _DynamicVC_Print, 
			       _DynamicVC_Copy,
			       (Stg_Component_DefaultConstructorFunction*)DynamicVC_DefaultNew,
			       _DynamicVC_Construct,
			       _DynamicVC_Build,
			       _DynamicVC_Initialise,
			       _DynamicVC_Execute,
			       _DynamicVC_Destroy,
			       name,
			       True,	
			       NULL,
			       _DynamicVC_PrintConcise,
			       _DynamicVC_ReadDictionary,
			       _DynamicVC_GetSet, 
			       _DynamicVC_GetVariableCount, 
			       _DynamicVC_GetVariableIndex, 
			       _DynamicVC_GetValueIndex, 
			       _DynamicVC_GetValueCount, 
			       _DynamicVC_GetValue,
			       DynamicVC_Apply, 
			       varReg, 
			       conFuncReg, 
			       dict );
}

DynamicVC* DynamicVC_DefaultNew( Name name ) {
	return (DynamicVC*)_DynamicVC_New( sizeof(DynamicVC), 
					   DynamicVC_Type, 
					   _DynamicVC_Delete, 
					   _DynamicVC_Print, 
					   _DynamicVC_Copy,
					   (Stg_Component_DefaultConstructorFunction*)DynamicVC_DefaultNew,
					   _VariableCondition_Construct,
					   _VariableCondition_Build,
					   _VariableCondition_Initialise,
					   _VariableCondition_Execute,
					   _VariableCondition_Destroy,
					   name, 
					   False,
					   NULL,
					   _DynamicVC_PrintConcise,
					   _DynamicVC_ReadDictionary,
					   _DynamicVC_GetSet, 
					   _DynamicVC_GetVariableCount, 
					   _DynamicVC_GetVariableIndex, 
					   _DynamicVC_GetValueIndex, 
					   _DynamicVC_GetValueCount, 
					   _DynamicVC_GetValue,
					   DynamicVC_Apply, 
					   NULL, 
					   NULL, 
					   NULL );
}

void DynamicVC_Init( DynamicVC* self,
		     Name name,
		     Variable_Register* varReg, 
		     ConditionFunction_Register* conFuncReg, 
		     Dictionary* dict )
{
	/* General info */
	self->type = DynamicVC_Type;
	self->_sizeOfSelf = sizeof(DynamicVC);
	self->_deleteSelf = False;
	
	/* Virtual info */
	self->_delete = _DynamicVC_Delete;
	self->_print = _DynamicVC_Print;
	self->_copy = _DynamicVC_Copy;
	self->_build = _VariableCondition_Build;
	self->_initialise = _VariableCondition_Initialise;
	self->_execute = _VariableCondition_Execute;
	self->_buildSelf = NULL;
	self->_printConcise = _DynamicVC_PrintConcise;
	self->_readDictionary = _DynamicVC_ReadDictionary;
	self->_getSet = _DynamicVC_GetSet;
	self->_getVariableCount = _DynamicVC_GetVariableCount;
	self->_getVariableIndex = _DynamicVC_GetVariableIndex;
	self->_getValueIndex = _DynamicVC_GetValueIndex;
	self->_getValueCount = _DynamicVC_GetValueCount;
	self->_getValue = _DynamicVC_GetValue;
	self->_apply = DynamicVC_Apply;
	
	_Stg_Class_Init( (Stg_Class*)self );
	_Stg_Object_Init( (Stg_Object*)self, name, NON_GLOBAL );
	_Stg_Component_Init( (Stg_Component*)self );
	_VariableCondition_Init( (VariableCondition*)self, varReg, conFuncReg, dict );

	/* Stg_Class info */
	_DynamicVC_Init( self );
}

DynamicVC* _DynamicVC_New( SizeT _sizeOfSelf, 
			   Type type,
			   Stg_Class_DeleteFunction* _delete,
			   Stg_Class_PrintFunction* _print,
			   Stg_Class_CopyFunction* _copy, 
			   Stg_Component_DefaultConstructorFunction* _defaultConstructor,
			   Stg_Component_ConstructFunction* _construct,
			   Stg_Component_BuildFunction* _build,
			   Stg_Component_InitialiseFunction* _initialise,
			   Stg_Component_ExecuteFunction* _execute,
			   Stg_Component_DestroyFunction* _destroy,
			   Name name,
			   Bool initFlag,
			   VariableCondition_BuildSelfFunc* _buildSelf, 
			   VariableCondition_PrintConciseFunc* _printConcise,
			   VariableCondition_ReadDictionaryFunc* _readDictionary,
			   VariableCondition_GetSetFunc* _getSet,
			   VariableCondition_GetVariableCountFunc* _getVariableCount,
			   VariableCondition_GetVariableIndexFunc* _getVariableIndex,
			   VariableCondition_GetValueIndexFunc* _getValueIndex,
			   VariableCondition_GetValueCountFunc* _getValueCount,
			   VariableCondition_GetValueFunc* _getValue,
			   VariableCondition_ApplyFunc* _apply, 
			   Variable_Register* varReg, 
			   ConditionFunction_Register* conFuncReg, 
			   Dictionary* dict )
{
	DynamicVC* self;

	/* Allocate memory/General info */
	assert(_sizeOfSelf >= sizeof(DynamicVC));
	self = (DynamicVC*)_VariableCondition_New( _sizeOfSelf, 
						   type, 
						   _delete, 
						   _print, 
						   _copy, 
						   _defaultConstructor,
						   _construct,
						   _build,
						   _initialise,
						   _execute,
						   _destroy,
						   name,
						   initFlag,
						   _buildSelf, 
						   _printConcise,
						   _readDictionary,
						   _getSet, 
						   _getVariableCount, 
						   _getVariableIndex, 
						   _getValueIndex, 
						   _getValueCount, 
						   _getValue,
						   _apply, 
						   varReg, 
						   conFuncReg, 
						   dict );
	
	/* Virtual info */
	
	/* Stg_Class info */
	if( initFlag )
		_DynamicVC_Init( self );
	
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

void _DynamicVC_Construct( void* vc, Stg_ComponentFactory* cf, void* data ) {
}

void _DynamicVC_Build( void* vc, void* data ) {
}

void _DynamicVC_Initialise( void* vc, void* data ) {
}

void _DynamicVC_Execute( void* vc, void* data ) {
}

void _DynamicVC_Destroy( void* vc, void* data ) {
}

void _DynamicVC_ReadDictionary( void* vc, void* dict ) {
}

void _DynamicVC_Delete( void* vc ) {
	DynamicVC* self = (DynamicVC*)vc;

	NewClass_Delete( self->vcMap );

	/* Stg_Class_Delete parent */
	_VariableCondition_Delete(self);
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
