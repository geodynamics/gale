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
** $Id: AllNodesVC.c 4153 2007-07-26 02:25:22Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>

#include <StgDomain/Geometry/Geometry.h>
#include <StgDomain/Shape/Shape.h>
#include <StgDomain/Mesh/Mesh.h>

#include "types.h"
#include "AllNodesVC.h"
#include "RegularMeshUtils.h"

#include <string.h>
#include <assert.h>


const Type AllNodesVC_Type = "AllNodesVC";
const Name defaultAllNodesVCName = "defaultAllNodesVCName";


/*--------------------------------------------------------------------------------------------------------------------------
** Constructor
*/

VariableCondition* AllNodesVC_Factory(
	AbstractContext*					context,
	Variable_Register*				variable_Register, 
	ConditionFunction_Register*	conFunc_Register, 
	Dictionary*							dictionary,
	void*									data )
{
	return (VariableCondition*)AllNodesVC_New( defaultAllNodesVCName, context, NULL, variable_Register, conFunc_Register, dictionary, data );
}


AllNodesVC*	AllNodesVC_New(
	Name									name,
	AbstractContext*					context,
	Name									_dictionaryEntryName, 
	Variable_Register*				variable_Register, 
	ConditionFunction_Register*	conFunc_Register,
	Dictionary*							dictionary,
	void*									mesh )
{
	AllNodesVC* self = _AllNodesVC_DefaultNew( name );

	self->isConstructed = True;
	_VariableCondition_Init( self, context, variable_Register, conFunc_Register, dictionary );	
	_AllNodesVC_Init( self, _dictionaryEntryName, mesh );

	return self;
}

AllNodesVC* _AllNodesVC_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                               _sizeOfSelf = sizeof(AllNodesVC);
	Type                                                       type = AllNodesVC_Type;
	Stg_Class_DeleteFunction*                               _delete = _AllNodesVC_Delete;
	Stg_Class_PrintFunction*                                 _print = _AllNodesVC_Print;
	Stg_Class_CopyFunction*                                   _copy = _AllNodesVC_Copy;
	Stg_Component_DefaultConstructorFunction*   _defaultConstructor = (Stg_Component_DefaultConstructorFunction*)_AllNodesVC_DefaultNew;
	Stg_Component_ConstructFunction*                     _construct = _AllNodesVC_AssignFromXML;
	Stg_Component_BuildFunction*                             _build = _AllNodesVC_Build;
	Stg_Component_InitialiseFunction*                   _initialise = _VariableCondition_Initialise;
	Stg_Component_ExecuteFunction*                         _execute = _VariableCondition_Execute;
	Stg_Component_DestroyFunction*                         _destroy = _AllNodesVC_Destroy;
	AllocationType                               nameAllocationType = NON_GLOBAL;
	VariableCondition_BuildSelfFunc*                     _buildSelf = _AllNodesVC_BuildSelf;
	VariableCondition_PrintConciseFunc*               _printConcise = _AllNodesVC_PrintConcise;
	VariableCondition_ReadDictionaryFunc*           _readDictionary = _AllNodesVC_ReadDictionary;
	VariableCondition_GetSetFunc*                           _getSet = _AllNodesVC_GetSet;
	VariableCondition_GetVariableCountFunc*       _getVariableCount = _AllNodesVC_GetVariableCount;
	VariableCondition_GetVariableIndexFunc*       _getVariableIndex = _AllNodesVC_GetVariableIndex;
	VariableCondition_GetValueIndexFunc*             _getValueIndex = _AllNodesVC_GetValueIndex;
	VariableCondition_GetValueCountFunc*             _getValueCount = _AllNodesVC_GetValueCount;
	VariableCondition_GetValueFunc*                       _getValue = _AllNodesVC_GetValue;
	VariableCondition_ApplyFunc*                             _apply = _VariableCondition_Apply;

	return _AllNodesVC_New(  ALLNODESVC_PASSARGS  );
}

AllNodesVC* _AllNodesVC_New(  ALLNODESVC_DEFARGS  ) {
	AllNodesVC*	self;
	
	/* Allocate memory/General info */
	assert( _sizeOfSelf >= sizeof(AllNodesVC) );
	self = (AllNodesVC*)_VariableCondition_New(  VARIABLECONDITION_PASSARGS  );
	
	/* Virtual info */
	
	/* Stg_Class info */
	
	return self;
}

void _AllNodesVC_Init(
	void*	allNodesVC, 
	Name	_dictionaryEntryName, 
	void*	mesh )
{
	AllNodesVC* self = (AllNodesVC*)allNodesVC;
	
	self->isConstructed = True;
	self->_dictionaryEntryName = _dictionaryEntryName;
	self->mesh = (Mesh*)mesh;
}

/*--------------------------------------------------------------------------------------------------------------------------
** General virtual functions
*/

void _AllNodesVC_ReadDictionary( void* variableCondition, void* dictionary ) {
	AllNodesVC*					self = (AllNodesVC*)variableCondition;
	Dictionary_Entry_Value*	vcDictVal;
	Dictionary_Entry_Value	_vcDictVal;
	Dictionary_Entry_Value*	varsVal;
	AllNodesVC_Entry_Index	entry_I;
	
	/* Find dictionary entry */
	if (self->_dictionaryEntryName)
		vcDictVal = Dictionary_Get( dictionary, (Dictionary_Entry_Key)self->_dictionaryEntryName  );
	else {
		vcDictVal = &_vcDictVal;
		Dictionary_Entry_Value_InitFromStruct( vcDictVal, dictionary );
	}
	
	if (vcDictVal) {
		/* Obtain the variable entries */
		self->_entryCount = Dictionary_Entry_Value_GetCount(Dictionary_Entry_Value_GetMember( vcDictVal, (Dictionary_Entry_Key)"variables") );
		self->_entryTbl = Memory_Alloc_Array( AllNodesVC_Entry, self->_entryCount, "AllNodesVC->_entryTbl" );
		varsVal = Dictionary_Entry_Value_GetMember( vcDictVal, (Dictionary_Entry_Key)"variables");
		
		for (entry_I = 0; entry_I < self->_entryCount; entry_I++ ) {
			char*							valType;
			Dictionary_Entry_Value*	valueEntry;
			Dictionary_Entry_Value*	varDictListVal;
			
			varDictListVal = Dictionary_Entry_Value_GetElement(varsVal, entry_I);
			valueEntry = Dictionary_Entry_Value_GetMember( varDictListVal, (Dictionary_Entry_Key)"value" );
			
			self->_entryTbl[entry_I].varName = Dictionary_Entry_Value_AsString( Dictionary_Entry_Value_GetMember( varDictListVal, (Dictionary_Entry_Key)"name") );
				
			valType = Dictionary_Entry_Value_AsString(Dictionary_Entry_Value_GetMember( varDictListVal, (Dictionary_Entry_Key)"type") );

			if (0 == strcasecmp(valType, "func")) {
				char*	funcName = Dictionary_Entry_Value_AsString(valueEntry);
				Index	cfIndex;
				
				self->_entryTbl[entry_I].value.type = VC_ValueType_CFIndex;
				cfIndex = ConditionFunction_Register_GetIndex( self->conFunc_Register, funcName);

				if ( cfIndex == (unsigned)-1 ) {	
					Stream*	errorStr = Journal_Register( Error_Type, (Name)self->type  );

					Journal_Printf( errorStr, "Error- in %s: While parsing "
						"definition of %s \"%s\", the cond. func. applied to "
						"variable \"%s\" - \"%s\" - wasn't found in the c.f. register.\n",
						__func__, self->type, self->_dictionaryEntryName, 
						self->_entryTbl[entry_I].varName, funcName );
					Journal_Printf( errorStr, "(Available functions in the C.F. register are: ");	
					ConditionFunction_Register_PrintNameOfEachFunc( self->conFunc_Register, errorStr );
					Journal_Printf( errorStr, ")\n");	
					assert(0);
				}	
				self->_entryTbl[entry_I].value.as.typeCFIndex = cfIndex;
			}	
			else if (!strcasecmp(valType, "array")) {
				Dictionary_Entry_Value*	valueElement;
				Index			i;

				self->_entryTbl[entry_I].value.type = VC_ValueType_DoubleArray;
				self->_entryTbl[entry_I].value.as.typeArray.size = Dictionary_Entry_Value_GetCount(valueEntry);
				self->_entryTbl[entry_I].value.as.typeArray.array = Memory_Alloc_Array( double, 
					self->_entryTbl[entry_I].value.as.typeArray.size, "AllNodesVC->_entryTbl[].value.as.typeArray.array" );
					
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
				self->_entryTbl[entry_I].value.as.typePtr = (void*) ((ArithPointer)Dictionary_Entry_Value_AsUnsignedInt( valueEntry ));
			}
			else {
				/* Assume double */
				Journal_Printf( 
					Journal_Register( InfoStream_Type, (Name)"myStream"  ), 
					"Type to variable on variable condition not given, assuming double\n" );
				self->_entryTbl[entry_I].value.type = VC_ValueType_Double;
				self->_entryTbl[entry_I].value.as.typeDouble = Dictionary_Entry_Value_AsDouble( valueEntry );
			}
		}
	}
	else
	{
		self->_entryCount = 0;
		self->_entryTbl = NULL;
	}
}

void _AllNodesVC_Delete( void* allNodesVC ) {
	AllNodesVC* self = (AllNodesVC*)allNodesVC;
	
	/* Stg_Class_Delete parent */
	_VariableCondition_Delete(self);
}

void _AllNodesVC_Destroy( void* allNodesVC, void* data ) {
	AllNodesVC* self = (AllNodesVC*)allNodesVC;

	if (self->_entryTbl) Memory_Free(self->_entryTbl);

	_VariableCondition_Destroy( self, data );
}

void _AllNodesVC_Print( void* allNodesVC, Stream* stream ) {
	AllNodesVC*					self = (AllNodesVC*)allNodesVC;
	AllNodesVC_Entry_Index	entry_I;
	Index							i;
	
	/* Set the Journal for printing informations */
	Stream* info = stream;
	
	/* General info */
	Journal_Printf( info, "AllNodesVC (ptr): %p\n", self);
	
	/* Virtual info */
	
	/* Stg_Class info */
	Journal_Printf( info, "\tdictionary (ptr): %p\n", self->dictionary);
	Journal_Printf( info, "\t_dictionaryEntryName (ptr): %p\n", self->_dictionaryEntryName);

	if (self->_dictionaryEntryName)
		Journal_Printf( info, "\t\t_dictionaryEntryName: %s\n", self->_dictionaryEntryName);

	Journal_Printf( info, "\t_entryCount: %u\n", self->_entryCount);
	Journal_Printf( info, "\t_entryTbl (ptr): %p\n", self->_entryTbl);

	if (self->_entryTbl)
		for (entry_I = 0; entry_I < self->_entryCount; entry_I++) {
			Journal_Printf( info, "\t\t_entryTbl[%u]:\n", entry_I);
			Journal_Printf( info, "\t\t\tvarName (ptr): %p\n", self->_entryTbl[entry_I].varName);

			if (self->_entryTbl[entry_I].varName)
				Journal_Printf( info, "\t\t\t\tvarName: %s\n", self->_entryTbl[entry_I].varName);
			Journal_Printf( info, "\t\t\tvalue:\n");

			switch (self->_entryTbl[entry_I].value.type) {
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
					Journal_Printf( info, "\t\t\t\tasDoubleArray (ptr): %p\n", self->_entryTbl[entry_I].value.as.typeArray.array);

					if (self->_entryTbl[entry_I].value.as.typeArray.array)
						for (i = 0; i < self->_entryTbl[entry_I].value.as.typeArray.size; i++)
							Journal_Printf( info, "\t\t\t\tasDoubleArray[%u]: %g\n", i,
								self->_entryTbl[entry_I].value.as.typeArray.array[i]);
					break;
					
				case VC_ValueType_CFIndex:
					Journal_Printf( info, "\t\t\t\ttype: VC_ValueType_CFIndex\n");
					Journal_Printf( info, "\t\t\t\tasCFIndex: %u\n", self->_entryTbl[entry_I].value.as.typeCFIndex);
					break;
			}
		}
	
	/* Print parent */
	_VariableCondition_Print(self);
}


void* _AllNodesVC_Copy( void* allNodesVC, void* dest, Bool deep, Name nameExt, struct PtrMap* ptrMap ) {
	AllNodesVC*		self = (AllNodesVC*)allNodesVC;
	AllNodesVC*		newAllNodesVC;
	PtrMap*			map = ptrMap;
	Bool			ownMap = False;
	
	if( !map ) {
		map = PtrMap_New( 10 );
		ownMap = True;
	}
	
	newAllNodesVC = (AllNodesVC*)_VariableCondition_Copy( self, dest, deep, nameExt, map );
	
	newAllNodesVC->_dictionaryEntryName = self->_dictionaryEntryName;
	newAllNodesVC->_entryCount = self->_entryCount;
	
	if( deep ) {
		newAllNodesVC->mesh = (Mesh*)Stg_Class_Copy( self->mesh, NULL, deep, nameExt, map );
		
		if( (newAllNodesVC->_entryTbl = PtrMap_Find( map, self->_entryTbl )) == NULL && self->_entryTbl ) {
			newAllNodesVC->_entryTbl = Memory_Alloc_Array( AllNodesVC_Entry, newAllNodesVC->_entryCount, "AllNodesVC->_entryTbl");
			memcpy( newAllNodesVC->_entryTbl, self->_entryTbl, sizeof(AllNodesVC_Entry) * newAllNodesVC->_entryCount );
			PtrMap_Append( map, newAllNodesVC->_entryTbl, self->_entryTbl );
		}
	}
	else {
		newAllNodesVC->mesh = self->mesh;
		newAllNodesVC->_entryTbl = self->_entryTbl;
	}
	
	if( ownMap ) {
		Stg_Class_Delete( map );
	}
	
	return (void*)newAllNodesVC;
}


void _AllNodesVC_Build( void* allNodesVC, void* data ) {
	AllNodesVC*		self = (AllNodesVC*)allNodesVC;
	
	_AllNodesVC_BuildSelf( self, data );
	
	_VariableCondition_Build( self, data );
}


/*--------------------------------------------------------------------------------------------------------------------------
** Macros
*/


/*--------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _AllNodesVC_AssignFromXML( void* allNodesVC, Stg_ComponentFactory* cf, void* data ) { 
}

void _AllNodesVC_BuildSelf( void* allNodesVC, void* data ) {
	AllNodesVC*		self = (AllNodesVC*)allNodesVC;
	
	if( self->mesh ) {
		Stg_Component_Build( self->mesh, data, False );
	}
}


IndexSet* _AllNodesVC_GetSet( void* variableCondition ) {
	AllNodesVC*		self = (AllNodesVC*)variableCondition;
	IndexSet*		set;

	set = IndexSet_New( Mesh_GetDomainSize( self->mesh, MT_VERTEX ) );
	IndexSet_AddAll( set );
	
	return set;
}


VariableCondition_VariableIndex _AllNodesVC_GetVariableCount( void* variableCondition, Index globalIndex ) {
	AllNodesVC*	self = (AllNodesVC*)variableCondition;
	
	return self->_entryCount;
}


Variable_Index _AllNodesVC_GetVariableIndex(
		void*				variableCondition, 
		Index				globalIndex, 
		VariableCondition_VariableIndex varIndex )
{
	AllNodesVC*	self = (AllNodesVC*)variableCondition;
	
	return Variable_Register_GetIndex(self->variable_Register, self->_entryTbl[varIndex].varName);
}


VariableCondition_ValueIndex _AllNodesVC_GetValueIndex(
		void*				variableCondition, 
		Index				globalIndex, 
		VariableCondition_VariableIndex	varIndex )
{
	return varIndex;
}


VariableCondition_ValueIndex _AllNodesVC_GetValueCount( void* variableCondition ) {
	AllNodesVC*	self = (AllNodesVC*)variableCondition;
	
	return self->_entryCount;
}


VariableCondition_Value _AllNodesVC_GetValue(void* variableCondition, VariableCondition_ValueIndex valIndex) {
	AllNodesVC*	self = (AllNodesVC*)variableCondition;

	return self->_entryTbl[valIndex].value;
}

void _AllNodesVC_PrintConcise( void* variableCondition, Stream* stream ) {
	AllNodesVC*		self = (AllNodesVC*)variableCondition;
	
	Journal_Printf( stream, "\ttype: %s, set: all\n", self->type );
}

/*--------------------------------------------------------------------------------------------------------------------------
** Build functions
*/


/*--------------------------------------------------------------------------------------------------------------------------
** Functions
*/


