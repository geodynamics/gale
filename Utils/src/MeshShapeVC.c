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
** $Id: MeshShapeVC.c 4160 2007-07-30 06:17:06Z DavidLee $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>

#include <StgDomain/Geometry/Geometry.h>
#include <StgDomain/Shape/Shape.h>
#include <StgDomain/Mesh/Mesh.h>

#include "types.h"
#include "MeshShapeVC.h"

#include <assert.h>
#include <string.h>

const Type MeshShapeVC_Type = "MeshShapeVC";
const Name defaultMeshShapeVCName = "defaultMeshShapeVCName";

/*-----------------------------------------------------------------------------------------------------------------
** Constructor
*/
VariableCondition* MeshShapeVC_Factory(
	Variable_Register*				variable_Register, 
	ConditionFunction_Register*	conFunc_Register, 
	Dictionary*							dictionary,
	void*									data )
{
	return (VariableCondition*) MeshShapeVC_New( defaultMeshShapeVCName, NULL, variable_Register, conFunc_Register, dictionary, (Mesh*)data );
}

MeshShapeVC* MeshShapeVC_New(
	Name									name,
	Name									_dictionaryEntryName, 
	Variable_Register*				variable_Register, 
	ConditionFunction_Register*	conFunc_Register, 
	Dictionary*							dictionary,
	void*									_mesh )
{
	MeshShapeVC* self = (MeshShapeVC*) _MeshShapeVC_DefaultNew( name );

	self->isConstructed = True;
	_VariableCondition_Init( self, variable_Register, conFunc_Register, dictionary );
	_MeshShapeVC_Init( self, _dictionaryEntryName, _mesh );

	return self;
}

void* _MeshShapeVC_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                               _sizeOfSelf = sizeof(MeshShapeVC);
	Type                                                       type = MeshShapeVC_Type;
	Stg_Class_DeleteFunction*                               _delete = _MeshShapeVC_Delete;
	Stg_Class_PrintFunction*                                 _print = _MeshShapeVC_Print;
	Stg_Class_CopyFunction*                                   _copy = _MeshShapeVC_Copy;
	Stg_Component_DefaultConstructorFunction*   _defaultConstructor = _MeshShapeVC_DefaultNew;
	Stg_Component_ConstructFunction*                     _construct = _MeshShapeVC_AssignFromXML;
	Stg_Component_BuildFunction*                             _build = _MeshShapeVC_Build;
	Stg_Component_InitialiseFunction*                   _initialise = _VariableCondition_Initialise;
	Stg_Component_ExecuteFunction*                         _execute = _VariableCondition_Execute;
	Stg_Component_DestroyFunction*                         _destroy = _MeshShapeVC_Destroy;
	AllocationType                               nameAllocationType = NON_GLOBAL;
	VariableCondition_BuildSelfFunc*                     _buildSelf = _MeshShapeVC_BuildSelf;
	VariableCondition_PrintConciseFunc*               _printConcise = _MeshShapeVC_PrintConcise;
	VariableCondition_ReadDictionaryFunc*           _readDictionary = _MeshShapeVC_ReadDictionary;
	VariableCondition_GetSetFunc*                           _getSet = _MeshShapeVC_GetSet;
	VariableCondition_GetVariableCountFunc*       _getVariableCount = _MeshShapeVC_GetVariableCount;
	VariableCondition_GetVariableIndexFunc*       _getVariableIndex = _MeshShapeVC_GetVariableIndex;
	VariableCondition_GetValueIndexFunc*             _getValueIndex = _MeshShapeVC_GetValueIndex;
	VariableCondition_GetValueCountFunc*             _getValueCount = _MeshShapeVC_GetValueCount;
	VariableCondition_GetValueFunc*                       _getValue = _MeshShapeVC_GetValue;
	VariableCondition_ApplyFunc*                             _apply = _VariableCondition_Apply;

	return (void*) _MeshShapeVC_New(  MESHSHAPEVC_PASSARGS  );
}

MeshShapeVC* _MeshShapeVC_New(  MESHSHAPEVC_DEFARGS  ) {
	MeshShapeVC* self;
	
	/* Allocate memory/General info */
	assert( _sizeOfSelf >= sizeof(MeshShapeVC) );
	self = (MeshShapeVC*)_VariableCondition_New(  VARIABLECONDITION_PASSARGS  );
	
	/* Virtual info */
	
	return self;
}

void _MeshShapeVC_Init( void* variableCondition, Name _dictionaryEntryName, void* _mesh ) {
	MeshShapeVC* self = (MeshShapeVC*) variableCondition;

	self->isConstructed        = True;
	self->_dictionaryEntryName = _dictionaryEntryName;
	self->_mesh                = (Mesh*)_mesh;
	self->_entryTbl            = 0;
	self->_entryCount          = 0;

	assert( _mesh && Stg_Class_IsInstance( _mesh, Mesh_Type ) );
}


/*--------------------------------------------------------------------------------------------------------------------------
** General virtual functions
*/


void _MeshShapeVC_Delete( void* variableCondition ) {
	MeshShapeVC* self = (MeshShapeVC*)variableCondition;
	
	/* Stg_Class_Delete parent */
	_VariableCondition_Delete(self);
}

void _MeshShapeVC_Destroy( void* variableCondition, void* data ) {
	MeshShapeVC* self = (MeshShapeVC*)variableCondition;
	
	if ( self->_entryTbl ) 
		Memory_Free(self->_entryTbl);

	if ( self->shapeName )
		Memory_Free( self->shapeName );
	
	/* Stg_Class_Delete parent */
	_VariableCondition_Destroy( self, data );
}

void _MeshShapeVC_Print(void* variableCondition, Stream* stream) {
	MeshShapeVC*                self = (MeshShapeVC*)variableCondition;
	MeshShapeVC_Entry_Index     entry_I;
	Index                   array_I;
	
	/* General info */
	Journal_Printf( stream, "MeshShapeVC (ptr): %p\n", self);
	
	/* Virtual info */
	
	/* Stg_Class info */
	Journal_Printf( stream, "\tdictionary (ptr): %p\n", self->dictionary);
	Journal_Printf( stream, "\t_dictionaryEntryName (ptr): %p\n", self->_dictionaryEntryName);
	if (self->_dictionaryEntryName)
		Journal_Printf( stream, "\t\t_dictionaryEntryName: %s\n", self->_dictionaryEntryName);
	if ( self->_shape )
		Journal_Printf( stream, "\t_shape: %s '%s'\n", self->_shape->type, self->_shape->name );
	
	Journal_Printf( stream, "\t_entryCount: %u\n", self->_entryCount);
	Journal_Printf( stream, "\t_entryTbl (ptr): %p\n", self->_entryTbl);
	if (self->_entryTbl) {
		for (entry_I = 0; entry_I < self->_entryCount; entry_I++) {
			Journal_Printf( stream, "\t\t_entryTbl[%u]:\n", entry_I);
			Journal_Printf( stream, "\t\t\tvarName (ptr): %p\n", self->_entryTbl[entry_I].varName);
			if (self->_entryTbl[entry_I].varName)
				Journal_Printf( stream, "\t\t\t\tvarName: %s\n", self->_entryTbl[entry_I].varName);
			Journal_Printf( stream, "\t\t\tvalue:\n");
			switch (self->_entryTbl[entry_I].value.type) {
				case VC_ValueType_Double:
					Journal_Printf( stream, "\t\t\t\ttype: VC_ValueType_Double\n" );
					Journal_Printf( stream, "\t\t\t\tasDouble: %g\n", self->_entryTbl[entry_I].value.as.typeDouble );
					break;
					
				case VC_ValueType_Int:
					Journal_Printf( stream, "\t\t\t\ttype: VC_ValueType_Int\n" );
					Journal_Printf( stream, "\t\t\t\tasInt: %i\n", self->_entryTbl[entry_I].value.as.typeInt );
					break;
					
				case VC_ValueType_Short:
					Journal_Printf( stream, "\t\t\t\ttype: VC_ValueType_Short\n" );
					Journal_Printf( stream, "\t\t\t\tasShort: %i\n", self->_entryTbl[entry_I].value.as.typeShort );
					break;
					
				case VC_ValueType_Char:
					Journal_Printf( stream, "\t\t\t\ttype: VC_ValueType_Char\n");
					Journal_Printf( stream, "\t\t\t\tasChar: %c\n", self->_entryTbl[entry_I].value.as.typeChar );
					break;
					
				case VC_ValueType_Ptr:
					Journal_Printf( stream, "\t\t\t\ttype: VC_ValueType_Ptr\n");
					Journal_Printf( stream, "\t\t\t\tasPtr: %g\n", self->_entryTbl[entry_I].value.as.typePtr );
					break;
					
				case VC_ValueType_DoubleArray:
					Journal_Printf( stream, "\t\t\t\ttype: VC_ValueType_DoubleArray\n");
					Journal_Printf( stream, "\t\t\t\tarraySize: %u\n", self->_entryTbl[entry_I].value.as.typeArray.size);
					Journal_Printf( stream, "\t\t\t\tasDoubleArray (ptr): %p\n", 
						self->_entryTbl[entry_I].value.as.typeArray.array);
					if (self->_entryTbl[entry_I].value.as.typeArray.array)
						for ( array_I = 0;  array_I < self->_entryTbl[entry_I].value.as.typeArray.size;  array_I++)
							Journal_Printf( stream, "\t\t\t\tasDoubleArray[%u]: %g\n",  array_I,
								self->_entryTbl[entry_I].value.as.typeArray.array[ array_I]);
					break;
					
				case VC_ValueType_CFIndex:
					Journal_Printf( stream, "\t\t\t\ttype: VC_ValueType_CFIndex\n");
					Journal_Printf( stream, "\t\t\t\tasCFIndex: %u\n", self->_entryTbl[entry_I].value.as.typeCFIndex);
					break;
			}
		}
	}
	Journal_Printf( stream, "\t_mesh (ptr): %p\n", self->_mesh);
	
	/* Print parent */
	_VariableCondition_Print( self );
}


void* _MeshShapeVC_Copy( void* variableCondition, void* dest, Bool deep, Name nameExt, struct PtrMap* ptrMap ) {
	MeshShapeVC*        self           = (MeshShapeVC*)variableCondition;
	MeshShapeVC*        newMeshShapeVC;
	PtrMap*         map            = ptrMap;
	Bool            ownMap         = False;
	
	if( !map ) {
		map = PtrMap_New( 10 );
		ownMap = True;
	}
	
	newMeshShapeVC = (MeshShapeVC*)_VariableCondition_Copy( self, dest, deep, nameExt, map );
	
	newMeshShapeVC->_dictionaryEntryName = self->_dictionaryEntryName;
	newMeshShapeVC->_shape = self->_shape;
	newMeshShapeVC->_entryCount = self->_entryCount;
	
	if( deep ) {
		newMeshShapeVC->_mesh = (Mesh*)Stg_Class_Copy( self->_mesh, NULL, deep, nameExt, map );
		
		if( (newMeshShapeVC->_entryTbl = PtrMap_Find( map, self->_entryTbl )) == NULL && self->_entryTbl ) {
			newMeshShapeVC->_entryTbl = Memory_Alloc_Array( MeshShapeVC_Entry, newMeshShapeVC->_entryCount, "MeshShapeVC->_entryTbl");
			memcpy( newMeshShapeVC->_entryTbl, self->_entryTbl, sizeof(MeshShapeVC_Entry) * newMeshShapeVC->_entryCount );
			PtrMap_Append( map, newMeshShapeVC->_entryTbl, self->_entryTbl );
		}
	}
	else {
		newMeshShapeVC->_mesh = self->_mesh;
		newMeshShapeVC->_entryTbl = self->_entryTbl;
	}
	
	if( ownMap ) {
		Stg_Class_Delete( map );
	}
	
	return (void*)newMeshShapeVC;
}
	
void _MeshShapeVC_AssignFromXML( void* variableCondition, Stg_ComponentFactory* cf, void* data ) {
}

void _MeshShapeVC_Build(  void* variableCondition, void* data ) {
	MeshShapeVC*			self = (MeshShapeVC*)variableCondition;

	_MeshShapeVC_BuildSelf( self, data );
	_VariableCondition_Build( self, data );
}
	
/****************** VariableCondition Virtual Functions ******************/
void _MeshShapeVC_BuildSelf(  void* variableCondition, void* data /* for build phase */ ) {
	MeshShapeVC*         self    = (MeshShapeVC*)variableCondition;
	AbstractContext* context = (AbstractContext*) data;

	assert( context && Stg_Class_IsInstance( context, AbstractContext_Type ) );
	assert( self->shapeName );
	Journal_Firewall( strlen( self->shapeName ) > 0, Journal_MyStream( Error_Type, self ),
			"You need to fill out the 'Shape' dictionary entry for this MeshShapeVC.\n" );
	assert( self->_mesh );

	self->_shape =  Stg_ComponentFactory_ConstructByName(  context->CF,  self->shapeName, Stg_Shape,  True, 0 /* dummy */  ) ;
	
	Stg_Component_Build( self->_mesh, data, False );
	Stg_Component_Build( self->_shape, data, False );
}

void _MeshShapeVC_PrintConcise( void* variableCondition, Stream* stream ) {
	MeshShapeVC* self = (MeshShapeVC*) variableCondition;
	
	Journal_Printf( stream, "\ttype: %s, Shape: %s '%s'", self->type, self->_shape->type, self->_shape->name );
}

void _MeshShapeVC_ReadDictionary( void* variableCondition, void* dictionary ) {
	MeshShapeVC*                  self = (MeshShapeVC*)variableCondition;
	Dictionary_Entry_Value*   vcDictVal;
	Dictionary_Entry_Value    _vcDictVal;
	Dictionary_Entry_Value*   varsVal;
	MeshShapeVC_Entry_Index	      entry_I;
	
	/* Find dictionary entry */
	if (self->_dictionaryEntryName)
		vcDictVal = Dictionary_Get(dictionary, self->_dictionaryEntryName);
	else {
		vcDictVal = &_vcDictVal;
		Dictionary_Entry_Value_InitFromStruct(vcDictVal, dictionary);
	}
	
	if (vcDictVal) {
		/* Get Name of Shape from dictionary - Grab pointer to shape later on */
		self->shapeName = StG_Strdup( 
				Dictionary_Entry_Value_AsString(Dictionary_Entry_Value_GetMember(vcDictVal, "Shape" )) );

		/* Obtain the variable entries */
		self->_entryCount = Dictionary_Entry_Value_GetCount(Dictionary_Entry_Value_GetMember(vcDictVal, "variables"));
		self->_entryTbl = Memory_Alloc_Array( MeshShapeVC_Entry, self->_entryCount, "MeshShapeVC->_entryTbl" );
		varsVal = Dictionary_Entry_Value_GetMember(vcDictVal, "variables");
		
		for (entry_I = 0; entry_I < self->_entryCount; entry_I++) {
			char*			valType;
			Dictionary_Entry_Value*	valueEntry;
			Dictionary_Entry_Value*	varDictListVal;
			
			varDictListVal = Dictionary_Entry_Value_GetElement(varsVal, entry_I);
			valueEntry = Dictionary_Entry_Value_GetMember(varDictListVal, "value");
			
			self->_entryTbl[entry_I].varName = Dictionary_Entry_Value_AsString(
				Dictionary_Entry_Value_GetMember(varDictListVal, "name"));
				
			valType = Dictionary_Entry_Value_AsString(Dictionary_Entry_Value_GetMember(varDictListVal, "type"));
			if (0 == strcasecmp(valType, "func")) {
				char*	funcName = Dictionary_Entry_Value_AsString(valueEntry);
				Index	cfIndex;
				
				self->_entryTbl[entry_I].value.type = VC_ValueType_CFIndex;
				cfIndex = ConditionFunction_Register_GetIndex( self->conFunc_Register, funcName);
				if ( cfIndex == (Index) -1 ) {	
					Stream*	errorStr = Journal_Register( Error_Type, self->type );

					Journal_Printf( errorStr, "Error- in %s: While parsing "
						"definition of MeshShapeVC \"%s\" (applies to shape \"%s\"), the cond. func. applied to "
						"variable \"%s\" - \"%s\" - wasn't found in the c.f. register.\n",
						__func__, self->_dictionaryEntryName, self->shapeName,
						self->_entryTbl[entry_I].varName, funcName );
					Journal_Printf( errorStr, "(Available functions in the C.F. register are: ");	
					ConditionFunction_Register_PrintNameOfEachFunc( self->conFunc_Register, errorStr );
					Journal_Printf( errorStr, ")\n");	
					assert(0);
				}	
				self->_entryTbl[entry_I].value.as.typeCFIndex = cfIndex;
			}
			else if (0 == strcasecmp(valType, "array"))
			{
				Dictionary_Entry_Value*	valueElement;
				Index			i;

				self->_entryTbl[entry_I].value.type = VC_ValueType_DoubleArray;
				self->_entryTbl[entry_I].value.as.typeArray.size = Dictionary_Entry_Value_GetCount(valueEntry);
				self->_entryTbl[entry_I].value.as.typeArray.array = Memory_Alloc_Array( double,
					self->_entryTbl[entry_I].value.as.typeArray.size, "MeshShapeVC->_entryTbl[].value.as.typeArray.array" );
					
				for (i = 0; i < self->_entryTbl[entry_I].value.as.typeArray.size; i++)
				{
					valueElement = Dictionary_Entry_Value_GetElement(valueEntry, i);
					self->_entryTbl[entry_I].value.as.typeArray.array[i] = 
						Dictionary_Entry_Value_AsDouble(valueElement);
				}
			}
			else if( 0 == strcasecmp( valType, "double" ) || 0 == strcasecmp( valType, "d" ) ||
				0 == strcasecmp( valType, "float" ) || 0 == strcasecmp( valType, "f" ) )
			{
				self->_entryTbl[entry_I].value.type = VC_ValueType_Double;
				self->_entryTbl[entry_I].value.as.typeDouble = Dictionary_Entry_Value_AsDouble( valueEntry );
			}
			else if( 0 == strcasecmp( valType, "integer" ) || 0 == strcasecmp( valType, "int" ) || 0 == strcasecmp( valType, "i" ) ) {
				self->_entryTbl[entry_I].value.type = VC_ValueType_Int;
				self->_entryTbl[entry_I].value.as.typeInt = Dictionary_Entry_Value_AsUnsignedInt( valueEntry );
			}
			else if( 0 == strcasecmp( valType, "short" ) || 0 == strcasecmp( valType, "s" ) ) {
				self->_entryTbl[entry_I].value.type = VC_ValueType_Short;
				self->_entryTbl[entry_I].value.as.typeShort = Dictionary_Entry_Value_AsUnsignedInt( valueEntry );
			}
			else if( 0 == strcasecmp( valType, "char" ) || 0 == strcasecmp( valType, "c" ) ) {
				self->_entryTbl[entry_I].value.type = VC_ValueType_Char;
				self->_entryTbl[entry_I].value.as.typeChar = Dictionary_Entry_Value_AsUnsignedInt( valueEntry );
			}
			else if( 0 == strcasecmp( valType, "pointer" ) || 0 == strcasecmp( valType, "ptr" ) || 0 == strcasecmp( valType, "p" ) ) {
				self->_entryTbl[entry_I].value.type = VC_ValueType_Ptr;
				self->_entryTbl[entry_I].value.as.typePtr = (void*) ( (ArithPointer)Dictionary_Entry_Value_AsUnsignedInt( valueEntry ));
			}
			else {
				/* Assume double */
				Journal_DPrintf( 
					Journal_Register( InfoStream_Type, "myStream" ), 
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

IndexSet* _MeshShapeVC_GetSet(void* variableCondition) {
	MeshShapeVC*	self = (MeshShapeVC*)variableCondition;
	Mesh*		mesh = self->_mesh;
	IndexSet*	set;
	unsigned	v_i;

	Stg_Component_Initialise( mesh, NULL, False );

	set = IndexSet_New( Mesh_GetDomainSize( mesh, MT_VERTEX ) );

	for( v_i = 0; v_i < Mesh_GetDomainSize( mesh, MT_VERTEX ); v_i++ ) {
		if( Stg_Shape_IsCoordInside( self->_shape, Mesh_GetVertex( mesh, v_i ) ) )
			IndexSet_Add( set, v_i );
	}

	return set;
}

VariableCondition_VariableIndex _MeshShapeVC_GetVariableCount(void* variableCondition, Index globalIndex) {
	MeshShapeVC*	self = (MeshShapeVC*)variableCondition;

	return self->_entryCount;
}

Variable_Index _MeshShapeVC_GetVariableIndex(void* variableCondition, Index globalIndex, VariableCondition_VariableIndex varIndex) {
	MeshShapeVC*        self          = (MeshShapeVC*)variableCondition;
	Variable_Index  searchedIndex = 0;
	Stream*         errorStr      = Journal_Register( Error_Type, self->type );
	Name            varName;
	
	varName = self->_entryTbl[varIndex].varName;
	searchedIndex = Variable_Register_GetIndex(self->variable_Register, varName );
	
	Journal_Firewall( 
			( searchedIndex < self->variable_Register->count ),
			errorStr,
			"Error- in %s: searching for index of varIndex %u (\"%s\") at global node number %u failed"
			" - register returned index %u, greater than count %u.\n",
			__func__, varIndex, varName, globalIndex, searchedIndex, self->variable_Register->count );

	return searchedIndex; 
}


VariableCondition_ValueIndex _MeshShapeVC_GetValueIndex(
		void*                                       variableCondition, 
		Index                                       globalIndex, 
		VariableCondition_VariableIndex             varIndex )
{
	return varIndex;
}


VariableCondition_ValueIndex _MeshShapeVC_GetValueCount(void* variableCondition) {
	MeshShapeVC*	self = (MeshShapeVC*)variableCondition;
	
	return self->_entryCount;
}


VariableCondition_Value _MeshShapeVC_GetValue(void* variableCondition, VariableCondition_ValueIndex valIndex) {
	MeshShapeVC*	self = (MeshShapeVC*)variableCondition;

	return self->_entryTbl[valIndex].value;
}



