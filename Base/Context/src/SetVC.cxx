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
** $Id: SetVC.c 4153 2007-07-26 02:25:22Z LukeHodkinson $
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
#include "SetVC.h"

#include <string.h>
#include <assert.h>


const Type SetVC_Type = "SetVC";
const Name defaultSetVCName = "defaultSetVCName";

/*--------------------------------------------------------------------------------------------------------------------------
** Constructor
*/

VariableCondition* SetVC_Factory(
	AbstractContext*					context,
	Variable_Register*				variable_Register, 
	ConditionFunction_Register*	conFunc_Register, 
	Dictionary*							dictionary,
	void*									data )
{
  return (VariableCondition*)SetVC_New( (char*)defaultSetVCName, context, NULL, variable_Register, conFunc_Register, dictionary );
}

SetVC* SetVC_New(
	Name									name,
	AbstractContext*					context,
	const char*						_dictionaryEntryName, 
	Variable_Register*				variable_Register, 
	ConditionFunction_Register*	conFunc_Register,
	Dictionary*							dictionary )
{
	SetVC* self = _SetVC_DefaultNew( name );

	self->isConstructed = True;
	_VariableCondition_Init( self, context, variable_Register, conFunc_Register, dictionary );
	_SetVC_Init( self,  _dictionaryEntryName );

	return self;
}

SetVC* _SetVC_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                               _sizeOfSelf = sizeof(SetVC);
	Type                                                       type = SetVC_Type;
	Stg_Class_DeleteFunction*                               _delete = _SetVC_Delete;
	Stg_Class_PrintFunction*                                 _print = _SetVC_Print;
	Stg_Class_CopyFunction*                                   _copy = _SetVC_Copy;
	Stg_Component_DefaultConstructorFunction*   _defaultConstructor = (Stg_Component_DefaultConstructorFunction*)_SetVC_DefaultNew;
	Stg_Component_ConstructFunction*                     _construct = _VariableCondition_AssignFromXML;
	Stg_Component_BuildFunction*                             _build = _VariableCondition_Build;
	Stg_Component_InitialiseFunction*                   _initialise = _VariableCondition_Initialise;
	Stg_Component_ExecuteFunction*                         _execute = _VariableCondition_Execute;
	Stg_Component_DestroyFunction*                         _destroy = _SetVC_Destroy;
	AllocationType                               nameAllocationType = NON_GLOBAL;
	VariableCondition_BuildSelfFunc*                     _buildSelf = NULL;
	VariableCondition_PrintConciseFunc*               _printConcise = _SetVC_PrintConcise;
	VariableCondition_ReadDictionaryFunc*           _readDictionary = _SetVC_ReadDictionary;
	VariableCondition_GetSetFunc*                           _getSet = _SetVC_GetSet;
	VariableCondition_GetVariableCountFunc*       _getVariableCount = _SetVC_GetVariableCount;
	VariableCondition_GetVariableIndexFunc*       _getVariableIndex = _SetVC_GetVariableIndex;
	VariableCondition_GetValueIndexFunc*             _getValueIndex = _SetVC_GetValueIndex;
	VariableCondition_GetValueCountFunc*             _getValueCount = _SetVC_GetValueCount;
	VariableCondition_GetValueFunc*                       _getValue = _SetVC_GetValue;
	VariableCondition_ApplyFunc*                             _apply = _VariableCondition_Apply;

	return (SetVC*)_SetVC_New(  SETVC_PASSARGS  );
}

SetVC* _SetVC_New(  SETVC_DEFARGS  ) {
	SetVC* self;
	
	/* Allocate memory/General info */
	assert( _sizeOfSelf >= sizeof(SetVC) );
	self = (SetVC*)_VariableCondition_New(  VARIABLECONDITION_PASSARGS  );
	
	/* Virtual info */
	
	/* Stg_Class info */
	
	return self;
}

void _SetVC_Init( void* setVC, Name _dictionaryEntryName ) {
	SetVC* self = (SetVC*)setVC;
	
	self->_dictionaryEntryName = _dictionaryEntryName;
}	
	
/*--------------------------------------------------------------------------------------------------------------------------
** General virtual functions
*/

void _SetVC_ReadDictionary( void* setVC, void* dictionary ) {
	SetVC*			self = (SetVC*)setVC;
	Dictionary_Entry_Value*	vcDictVal;
	Dictionary_Entry_Value	_vcDictVal;
	Dictionary_Entry_Value*	varsVal;
	SetVC_Entry_Index	entry_I;
	
	
	/* Find dictionary entry */
	if (self->_dictionaryEntryName)
		vcDictVal = Dictionary_Get( dictionary, (Dictionary_Entry_Key)self->_dictionaryEntryName  );
	else {
		vcDictVal = &_vcDictVal;
		Dictionary_Entry_Value_InitFromStruct( vcDictVal, dictionary );
	}
	
	if (vcDictVal) {
		Dictionary_Entry_Value*		setVal = Dictionary_Entry_Value_GetMember( vcDictVal, (Dictionary_Entry_Key)"indices"  );
		Index				indexCnt = Dictionary_Entry_Value_AsUnsignedInt( 
							Dictionary_Entry_Value_GetMember( vcDictVal, (Dictionary_Entry_Key)"indexCount" )  );
		Index				i, cnt;
		
		self->_vcset = IndexSet_New( indexCnt );
		cnt = Dictionary_Entry_Value_GetCount( setVal );

		for( i = 0; i < cnt; i++ )
			IndexSet_Add( self->_vcset, Dictionary_Entry_Value_AsUnsignedInt( 
				Dictionary_Entry_Value_GetElement( setVal, i ) ) );
		
		/* Obtain the variable entries */
		varsVal = Dictionary_Entry_Value_GetMember( vcDictVal, (Dictionary_Entry_Key)"variables");
		self->_entryCount = Dictionary_Entry_Value_GetCount( varsVal  );
		self->_entryTbl = Memory_Alloc_Array( SetVC_Entry, self->_entryCount, "SetVC->_entryTbl");
		
		for (entry_I = 0; entry_I < self->_entryCount; entry_I++) {
			char* valType;
			Dictionary_Entry_Value*	valueEntry;
			Dictionary_Entry_Value*	varDictListVal;
			
			varDictListVal = Dictionary_Entry_Value_GetElement(varsVal, entry_I);
			valueEntry = Dictionary_Entry_Value_GetMember( varDictListVal, (Dictionary_Entry_Key)"value" );
			
			self->_entryTbl[entry_I].varName = Dictionary_Entry_Value_AsString(
				Dictionary_Entry_Value_GetMember( varDictListVal, (Dictionary_Entry_Key)"name") );
				
			valType = Dictionary_Entry_Value_AsString(Dictionary_Entry_Value_GetMember( varDictListVal, (Dictionary_Entry_Key)"type") );

			if (strlen(valType)==0 || 0==strcasecmp(valType, "equation"))
			{
                          self->_entryTbl[entry_I].value.type = VC_ValueType_Equation;
                          /* This leaks memory */
                          self->_entryTbl[entry_I].value.as.equation=
                            StG_Strdup(Dictionary_Entry_Value_AsString(valueEntry));
                        }
			else if (0 == strcasecmp(valType, "func"))
                          {
				char*	funcName = Dictionary_Entry_Value_AsString(valueEntry);
				
				self->_entryTbl[entry_I].value.type = VC_ValueType_CFIndex;
				self->_entryTbl[entry_I].value.as.typeCFIndex = ConditionFunction_Register_GetIndex(
					self->conFunc_Register, funcName);
                          }
			else if (!strcasecmp(valType, "array")) {
				Dictionary_Entry_Value*	valueElement;
				Index			i;

				self->_entryTbl[entry_I].value.type = VC_ValueType_DoubleArray;
				self->_entryTbl[entry_I].value.as.typeArray.size = Dictionary_Entry_Value_GetCount(valueEntry);
				self->_entryTbl[entry_I].value.as.typeArray.array = Memory_Alloc_Array( double,
					self->_entryTbl[entry_I].value.as.typeArray.size, "SetVC->_entryTbl[].value.as.typeArray.array" );
					
				for (i = 0; i < self->_entryTbl[entry_I].value.as.typeArray.size; i++) {
					valueElement = Dictionary_Entry_Value_GetElement(valueEntry, i);
					self->_entryTbl[entry_I].value.as.typeArray.array[i] = 
						Dictionary_Entry_Value_AsDouble(valueElement);
				}
			}
			else if( !strcasecmp( valType, "double" ) || !strcasecmp( valType, "d" ) || !strcasecmp( valType, "float" ) || !strcasecmp( valType, "f" ) ) {
				self->_entryTbl[entry_I].value.type = VC_ValueType_Double;
				self->_entryTbl[entry_I].value.as.typeDouble = Dictionary_Entry_Value_AsDouble( valueEntry );
			}
			else if( !strcasecmp( valType, "integer" ) || !strcasecmp( valType, "int" ) || !strcasecmp( valType, "i" ) ) {
				self->_entryTbl[entry_I].value.type = VC_ValueType_Int;
				self->_entryTbl[entry_I].value.as.typeInt = Dictionary_Entry_Value_AsUnsignedInt( valueEntry );
			}
			else if( !strcasecmp( valType, "short" ) || !strcasecmp( valType, "s" ) ) {
				self->_entryTbl[entry_I].value.type = VC_ValueType_Short;
				self->_entryTbl[entry_I].value.as.typeShort = Dictionary_Entry_Value_AsUnsignedInt( valueEntry );
			}
			else if( !strcasecmp( valType, "char" ) || !strcasecmp( valType, "c" ) ) {
				self->_entryTbl[entry_I].value.type = VC_ValueType_Char;
				self->_entryTbl[entry_I].value.as.typeChar = Dictionary_Entry_Value_AsUnsignedInt( valueEntry );
			}
			else if( !strcasecmp( valType, "pointer" ) || !strcasecmp( valType, "ptr" ) || !strcasecmp( valType, "p" ) ) {
				self->_entryTbl[entry_I].value.type = VC_ValueType_Ptr;
				self->_entryTbl[entry_I].value.as.typePtr = (void*) ((ArithPointer) Dictionary_Entry_Value_AsUnsignedInt( valueEntry ));
			}
			else {
                          Journal_Firewall(False,Journal_Register( Error_Type,SetVC_Type),
                                           "Unknown type for variable condition: %s\n",valType);
			}
		}
	}
	else {
		self->_entryCount = 0;
		self->_entryTbl = NULL;
	}
}

void _SetVC_Delete(void* setVC) {
	SetVC* self = (SetVC*)setVC;
	
	/* Stg_Class_Delete parent */
	_VariableCondition_Delete(self);
}

void _SetVC_Destroy(void* setVC, void* data) {
	SetVC* self = (SetVC*)setVC;
		
	if (self->_entryTbl) Memory_Free( self->_entryTbl );

	_VariableCondition_Destroy( self, data );
}

void _SetVC_Print(void* setVC, Stream* stream) {
	SetVC*				self = (SetVC*)setVC;
	SetVC_Entry_Index	entry_I;
	Index					i;
	
	/* Set the Journal for printing informations */
	Stream* info = stream;
	
	/* General info */
	Journal_Printf( info, "SetVC (ptr): %p\n", self);
	
	/* Virtual info */
	
	/* Stg_Class info */
	Journal_Printf( info, "\tdictionary (ptr): %p\n", self->dictionary);
	Journal_Printf( info, "\t_dictionaryEntryName (ptr): %p\n", self->_dictionaryEntryName);
	if (self->_dictionaryEntryName)
		Journal_Printf( info, "\t\t_dictionaryEntryName: %s\n", self->_dictionaryEntryName);
	Journal_Printf( info, "\t_entryCount: %u\n", self->_entryCount);
	Journal_Printf( info, "\t_entryTbl (ptr): %p\n", self->_entryTbl);
	if (self->_entryTbl)
		for (entry_I = 0; entry_I < self->_entryCount; entry_I++)
		{
			Journal_Printf( info, "\t\t_entryTbl[%u]:\n", entry_I);
			Journal_Printf( info, "\t\t\tvarName (ptr): %p\n", self->_entryTbl[entry_I].varName);
			if (self->_entryTbl[entry_I].varName)
				Journal_Printf( info, "\t\t\t\tvarName: %s\n", self->_entryTbl[entry_I].varName);
			Journal_Printf( info, "\t\t\tvalue:\n");
			switch (self->_entryTbl[entry_I].value.type)
			{
				case VC_ValueType_Double:
					Journal_Printf( info, "\t\t\t\ttype: VC_ValueType_Double\n" );
					Journal_Printf( info, "\t\t\t\tasDouble: %g\n", self->_entryTbl[entry_I].value.as.typeDouble );
					break;
					
				case VC_ValueType_Int:
					Journal_Printf( info, "\t\t\t\ttype: VC_ValueType_Int\n" );
					Journal_Printf( info, "\t\t\t\tasInt: %i\n", self->_entryTbl[entry_I].value.as.typeInt );
					break;
					
				case VC_ValueType_Short:
					Journal_Printf( info, "\t\t\t\ttype: VC_ValueType_Short\n" );
					Journal_Printf( info, "\t\t\t\tasShort: %i\n", self->_entryTbl[entry_I].value.as.typeShort );
					break;
					
				case VC_ValueType_Char:
					Journal_Printf( info, "\t\t\t\ttype: VC_ValueType_Char\n");
					Journal_Printf( info, "\t\t\t\tasChar: %c\n", self->_entryTbl[entry_I].value.as.typeChar );
					break;
					
				case VC_ValueType_Ptr:
					Journal_Printf( info, "\t\t\t\ttype: VC_ValueType_Ptr\n");
					Journal_Printf( info, "\t\t\t\tasPtr: %g\n", self->_entryTbl[entry_I].value.as.typePtr );
					break;
					
				case VC_ValueType_DoubleArray:
					Journal_Printf( info, "\t\t\t\ttype: VC_ValueType_DoubleArray\n");
					Journal_Printf( info, "\t\t\t\tarraySize: %u\n", self->_entryTbl[entry_I].value.as.typeArray.size);
					Journal_Printf( info, "\t\t\t\tasDoubleArray (ptr): %p\n", 
						self->_entryTbl[entry_I].value.as.typeArray.array);
					if (self->_entryTbl[entry_I].value.as.typeArray.array)
						for (i = 0; i < self->_entryTbl[entry_I].value.as.typeArray.size; i++)
							Journal_Printf( info, "\t\t\t\tasDoubleArray[%u]: %g\n", i,
								self->_entryTbl[entry_I].value.as.typeArray.array[i]);
					break;
					
				case VC_ValueType_CFIndex:
					Journal_Printf( info, "\t\t\t\ttype: VC_ValueType_CFIndex\n");
					Journal_Printf( info, "\t\t\t\tasCFIndex: %u\n", self->_entryTbl[entry_I].value.as.typeCFIndex);
					break;
				case VC_ValueType_Equation:
					Journal_Printf( info, "\t\t\t\ttype: VC_ValueType_Equation\n");
					Journal_Printf( info, "\t\t\t\tasEquation: %s\n", self->_entryTbl[entry_I].value.as.equation);
					break;
			}
		}
	
	/* Print parent */
	_VariableCondition_Print(self);
}


void* _SetVC_Copy( const void* setVC, void* dest, Bool deep, Name nameExt, struct PtrMap* ptrMap ) {
	SetVC*	self = (SetVC*)setVC;
	SetVC*	newSetVC;
	PtrMap*	map = ptrMap;
	Bool		ownMap = False;
	
	if( !map ) {
		map = PtrMap_New( 10 );
		ownMap = True;
	}
	
	newSetVC = (SetVC*)_VariableCondition_Copy( self, dest, deep, nameExt, map );
	
	newSetVC->_dictionaryEntryName = self->_dictionaryEntryName;
	newSetVC->_entryCount = self->_entryCount;
	
	if( deep ) {
		newSetVC->_vcset = (IndexSet*)Stg_Class_Copy( self->_vcset, NULL, deep, nameExt, map );
		
		if( (newSetVC->_entryTbl = (SetVC_Entry*)PtrMap_Find( map, self->_entryTbl )) == NULL && self->_entryTbl ) {
			newSetVC->_entryTbl = Memory_Alloc_Array( SetVC_Entry, newSetVC->_entryCount, "SetVC->_entryTbl");
			memcpy( newSetVC->_entryTbl, self->_entryTbl, sizeof(SetVC_Entry) * newSetVC->_entryCount );
			PtrMap_Append( map, newSetVC->_entryTbl, self->_entryTbl );
		}
	}
	else {
		newSetVC->_vcset = self->_vcset;
		newSetVC->_entryTbl = self->_entryTbl;
	}
	
	if( ownMap ) {
		Stg_Class_Delete( map );
	}
	
	return (void*)newSetVC;
}


/*--------------------------------------------------------------------------------------------------------------------------
** Macros
*/


/*--------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

IndexSet* _SetVC_GetSet( void* variableCondition ) {
	SetVC*		self = (SetVC*)variableCondition;
	
	return (IndexSet*) IndexSet_Duplicate( self->_vcset );
}


VariableCondition_VariableIndex _SetVC_GetVariableCount( void* variableCondition, Index globalIndex ) {
	SetVC*	self = (SetVC*)variableCondition;
	
	return self->_entryCount;
}


Variable_Index _SetVC_GetVariableIndex( void* variableCondition, Index globalIndex, VariableCondition_VariableIndex varIndex ) {
	SetVC*	self = (SetVC*)variableCondition;
	
	return Variable_Register_GetIndex(self->variable_Register, self->_entryTbl[varIndex].varName);
}


VariableCondition_ValueIndex _SetVC_GetValueIndex( void* variableCondition, Index globalIndex, VariableCondition_VariableIndex varIndex ) {
	return varIndex;
}


VariableCondition_ValueIndex _SetVC_GetValueCount( void* variableCondition ) {
	SetVC*	self = (SetVC*)variableCondition;
	
	return self->_entryCount;
}


VariableCondition_Value _SetVC_GetValue( void* variableCondition, VariableCondition_ValueIndex valIndex ) {
	SetVC*	self = (SetVC*)variableCondition;

	return self->_entryTbl[valIndex].value;
}


void _SetVC_PrintConcise( void* variableCondition, Stream* stream ) {
	SetVC*		self = (SetVC*)variableCondition;
	IndexSet_Index	set_I;
	
	Journal_Printf( stream, "\ttype: %s, set: {", self->type );
	for( set_I = 0; set_I < self->_set->size; set_I++ ) {
		if( IndexSet_IsMember( self->_set, set_I ) ) {
			Journal_Printf( stream, "%u ", set_I );
		}
	}
	Journal_Printf( stream, "}\n" );
}

/*--------------------------------------------------------------------------------------------------------------------------
** Build functions
*/


/*--------------------------------------------------------------------------------------------------------------------------
** Functions
*/


