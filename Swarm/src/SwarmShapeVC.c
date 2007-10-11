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
** $Id: SwarmShapeVC.c 4153 2007-07-26 02:25:22Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/Geometry/Geometry.h>
#include <StgDomain/Shape/Shape.h>
#include <StgDomain/Mesh/Mesh.h>
#include <StgDomain/Utils/Utils.h>

#include "types.h"
#include "SwarmShapeVC.h"

#include "Swarm.h"
#include <assert.h>
#include <string.h>

const Type SwarmShapeVC_Type = "SwarmShapeVC";
const Name defaultSwarmShapeVCName = "defaultSwarmShapeVCName";

/*-----------------------------------------------------------------------------------------------------------------
** Constructor
*/
VariableCondition* SwarmShapeVC_Factory(
		Variable_Register*                          variable_Register, 
		ConditionFunction_Register*                 conFunc_Register, 
		Dictionary*                                 dictionary,
		void*                                       data )
{
	return (VariableCondition*) 
		SwarmShapeVC_New( defaultSwarmShapeVCName, NULL, variable_Register, conFunc_Register, dictionary, (Swarm*)data );
}

SwarmShapeVC* SwarmShapeVC_New(
		Name                                        name,
		Name                                        _dictionaryEntryName, 
		Variable_Register*                          variable_Register, 
		ConditionFunction_Register*                 conFunc_Register, 
		Dictionary*	                                dictionary,
		void*                                       _swarm )
{
	SwarmShapeVC* self = (SwarmShapeVC*) _SwarmShapeVC_DefaultNew( name );

	_VariableCondition_Init( self, variable_Register, conFunc_Register, dictionary );
	_SwarmShapeVC_Init( self, _dictionaryEntryName, _swarm );

	return self;
}

SwarmShapeVC* _SwarmShapeVC_New( 
		SizeT                                       _sizeOfSelf, 
		Type                                        type,
		Stg_Class_DeleteFunction*                   _delete,
		Stg_Class_PrintFunction*                    _print,
		Stg_Class_CopyFunction*                     _copy,
		Stg_Component_DefaultConstructorFunction*   _defaultConstructor,
		Stg_Component_ConstructFunction*            _construct,
		Stg_Component_BuildFunction*                _build,
		Stg_Component_InitialiseFunction*           _initialise,
		Stg_Component_ExecuteFunction*              _execute,
		Stg_Component_DestroyFunction*              _destroy,
		VariableCondition_BuildSelfFunc*            _buildSelf, 
		VariableCondition_PrintConciseFunc*         _printConcise,
		VariableCondition_ReadDictionaryFunc*       _readDictionary,
		VariableCondition_GetSetFunc*               _getSet,
		VariableCondition_GetVariableCountFunc*     _getVariableCount,
		VariableCondition_GetVariableIndexFunc*     _getVariableIndex,
		VariableCondition_GetValueIndexFunc*        _getValueIndex,
		VariableCondition_GetValueCountFunc*        _getValueCount,
		VariableCondition_GetValueFunc*             _getValue,
		VariableCondition_ApplyFunc*	  	    _apply, 
		Name                                        name  )
{
	SwarmShapeVC*	self;
	
	/* Allocate memory/General info */
	assert(_sizeOfSelf >= sizeof(SwarmShapeVC));
	self = (SwarmShapeVC*)_VariableCondition_New(
		_sizeOfSelf, 
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
		False,
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
		NULL, 
		NULL,
		NULL );
	
	/* Virtual info */
	
	return self;
}


void _SwarmShapeVC_Init(
		void*                                       variableCondition, 
		Name                                        _dictionaryEntryName, 
		void*                                       _swarm )
{
	SwarmShapeVC*			self = (SwarmShapeVC*) variableCondition;

	self->isConstructed        = True;
	self->_dictionaryEntryName = _dictionaryEntryName;
	self->_swarm                = (Swarm*)_swarm;
	self->_entryTbl            = 0;
	self->_entryCount          = 0;

	assert( _swarm && Stg_Class_IsInstance( _swarm, Swarm_Type ) );
}


/*--------------------------------------------------------------------------------------------------------------------------
** General virtual functions
*/


void _SwarmShapeVC_Delete(void* variableCondition) {
	SwarmShapeVC*	self = (SwarmShapeVC*)variableCondition;
	
	if ( self->_entryTbl ) 
		Memory_Free(self->_entryTbl);

	if ( self->shapeName )
		Memory_Free( self->shapeName );
	
	/* Stg_Class_Delete parent */
	_VariableCondition_Delete(self);
}


void _SwarmShapeVC_Print(void* variableCondition, Stream* stream) {
	SwarmShapeVC*                self = (SwarmShapeVC*)variableCondition;
	SwarmShapeVC_Entry_Index     entry_I;
	Index                   array_I;
	
	/* General info */
	Journal_Printf( stream, "SwarmShapeVC (ptr): %p\n", self);
	
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
	Journal_Printf( stream, "\t_swarm (ptr): %p\n", self->_swarm);
	
	/* Print parent */
	_VariableCondition_Print( self );
}


void* _SwarmShapeVC_Copy( void* variableCondition, void* dest, Bool deep, Name nameExt, struct PtrMap* ptrMap ) {
	SwarmShapeVC*   self           = (SwarmShapeVC*)variableCondition;
	SwarmShapeVC*   newSwarmShapeVC;
	PtrMap*         map            = ptrMap;
	Bool            ownMap         = False;
	
	if( !map ) {
		map = PtrMap_New( 10 );
		ownMap = True;
	}
	
	newSwarmShapeVC = (SwarmShapeVC*)_VariableCondition_Copy( self, dest, deep, nameExt, map );
	
	newSwarmShapeVC->_dictionaryEntryName = self->_dictionaryEntryName;
	newSwarmShapeVC->_shape = self->_shape;
	newSwarmShapeVC->_entryCount = self->_entryCount;
	
	if( deep ) {
		newSwarmShapeVC->_swarm = (Swarm*)Stg_Class_Copy( self->_swarm, NULL, deep, nameExt, map );
		
		if( (newSwarmShapeVC->_entryTbl = PtrMap_Find( map, self->_entryTbl )) == NULL && self->_entryTbl ) {
			newSwarmShapeVC->_entryTbl = Memory_Alloc_Array( SwarmShapeVC_Entry, newSwarmShapeVC->_entryCount, "SwarmShapeVC->_entryTbl");
			memcpy( newSwarmShapeVC->_entryTbl, self->_entryTbl, sizeof(SwarmShapeVC_Entry) * newSwarmShapeVC->_entryCount );
			PtrMap_Append( map, newSwarmShapeVC->_entryTbl, self->_entryTbl );
		}
	}
	else {
		newSwarmShapeVC->_swarm = self->_swarm;
		newSwarmShapeVC->_entryTbl = self->_entryTbl;
	}
	
	if( ownMap ) {
		Stg_Class_Delete( map );
	}
	
	return (void*)newSwarmShapeVC;
}
	
/****************** Stg_Component Virtual Functions ******************/
void* _SwarmShapeVC_DefaultNew( Name name ) {
	return (void*) _SwarmShapeVC_New(
		sizeof(SwarmShapeVC), 
		SwarmShapeVC_Type, 
		_SwarmShapeVC_Delete, 
		_SwarmShapeVC_Print, 
		_SwarmShapeVC_Copy,
		_SwarmShapeVC_DefaultNew,
		_SwarmShapeVC_Construct,	
		_SwarmShapeVC_Build,
		/*_VariableCondition_Initialise,*/
		_SwarmShapeVC_Initialise,
		_VariableCondition_Execute,
		_VariableCondition_Destroy,
		_SwarmShapeVC_BuildSelf, 
		_SwarmShapeVC_PrintConcise,
		_SwarmShapeVC_ReadDictionary,
		_SwarmShapeVC_GetSet, 
		_SwarmShapeVC_GetVariableCount, 
		_SwarmShapeVC_GetVariableIndex, 
		_SwarmShapeVC_GetValueIndex, 
		_SwarmShapeVC_GetValueCount, 
		_SwarmShapeVC_GetValue,
		_VariableCondition_Apply, 
		name );
}

void _SwarmShapeVC_Construct( void* variableCondition, Stg_ComponentFactory* cf, void* data ) {
	SwarmShapeVC*	self    		= (SwarmShapeVC*)variableCondition;
	void*		conFunc_Register 	= NULL;
	void*		variable_Register 	= NULL;

	self = Stg_ComponentFactory_ConstructByName( cf, self->name, SwarmShapeVC, False, data );

	conFunc_Register = (void*)Stg_ObjectList_Get( cf->registerRegister, "ConditionFunction_Register" );
	assert( conFunc_Register );
	self->conFunc_Register = conFunc_Register; 
	
	variable_Register = (void*)Stg_ObjectList_Get( cf->registerRegister, "Variable_Register" );
	assert( variable_Register );
	self->variable_Register = variable_Register;

	_VariableCondition_Init( self, variable_Register, conFunc_Register, NULL );
	
	self->dictionary = Dictionary_GetDictionary( cf->componentDict, self->name );
}

void _SwarmShapeVC_Build(  void* variableCondition, void* data ) {
	SwarmShapeVC*		self 	= (SwarmShapeVC*)variableCondition;
	AbstractContext* 	context = (AbstractContext*) data;
	
	assert( context && Stg_Class_IsInstance( context, AbstractContext_Type ) );

	self->_swarm = Stg_ComponentFactory_ConstructByKey(  context->CF, self->name, "Swarm", Swarm,  True, 0  ) ; 
	assert( self->_swarm );
	
	self->_shape = Stg_ComponentFactory_ConstructByKey(  context->CF, self->name, "Shape", Stg_Shape,  True, 0 /* dummy */  ) ;
	assert( self->_shape );

	_VariableCondition_Build( self, data );

	Stg_Component_Build( self->_shape, data, False );
	/*_SwarmShapeVC_BuildSelf( self, data );*/
}
	
/****************** VariableCondition Virtual Functions ******************/
void _SwarmShapeVC_BuildSelf(  void* variableCondition, void* data /* for build phase */ ) {
	SwarmShapeVC*         	self    = (SwarmShapeVC*)variableCondition;
	AbstractContext* 	context = (AbstractContext*) data;

	assert( context && Stg_Class_IsInstance( context, AbstractContext_Type ) );

	/* dave - 06.08.07 */
	/*self->shapeName = Stg_ComponentFactory_GetString( context->CF, self->name, "Shape", "" );
*/
	/*Journal_Firewall( strlen( self->shapeName ) > 0, Journal_MyStream( Error_Type, self ),
			"You need to fill out the 'Shape' dictionary entry for this SwarmShapeVC.\n" );*/
	/*assert( self->shapeName );*/
	
	/*Stg_Component_Build( self->_swarm, data, False );*/ /* remove this call? */
	Stg_Component_Build( self->_shape, data, False );
}

/* added to call the porisity field standard condition plugin */
void _SwarmShapeVC_Initialise(  void* variableCondition, void* data ) {
	SwarmShapeVC*	self = (SwarmShapeVC*)variableCondition;

	_VariableCondition_Initialise( self, data );

	/* need to call the standard condition function here */
}
	
void _SwarmShapeVC_PrintConcise( void* variableCondition, Stream* stream ) {
	SwarmShapeVC* self = (SwarmShapeVC*) variableCondition;
	
	Journal_Printf( stream, "\ttype: %s, Shape: %s '%s'", self->type, self->_shape->type, self->_shape->name );
}

void _SwarmShapeVC_ReadDictionary( void* variableCondition, void* dict /**/ ) {
	SwarmShapeVC*                  	self = (SwarmShapeVC*)variableCondition;
	Dictionary_Entry_Value*   	vcDictVal;
	Dictionary_Entry_Value    	_vcDictVal;
	Dictionary_Entry_Value*   	varsVal;
	SwarmShapeVC_Entry_Index	entry_I;
	Dictionary*			dictionary = (Dictionary*)dict; /**/
	
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
		self->_entryTbl = Memory_Alloc_Array( SwarmShapeVC_Entry, self->_entryCount, "SwarmShapeVC->_entryTbl" );
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
						"definition of swarmShapeVC \"%s\" (applies to shape \"%s\"), the cond. func. applied to "
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
					self->_entryTbl[entry_I].value.as.typeArray.size, "SwarmShapeVC->_entryTbl[].value.as.typeArray.array" );
					
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

IndexSet* _SwarmShapeVC_GetSet(void* variableCondition) {
	SwarmShapeVC*		self = (SwarmShapeVC*)variableCondition;
	Swarm*			swarm = self->_swarm;
	IndexSet*		set;
	Index           	particleDomainCount;
	unsigned		particleIndex;
	GlobalParticle*		particle; /* check that this is ok for the given swarm type */

	/*Stg_Component_Initialise( swarm, NULL, False );*/
	/*Stg_Component_Build( swarm, NULL, False );*/

	particleDomainCount = swarm->particleLocalCount; /* sizeof(swarm particle array) - think this is right?? */
	set = IndexSet_New( particleDomainCount );

	for( particleIndex = 0; particleIndex < particleDomainCount; particleIndex++ ) {
		particle = (GlobalParticle*)Swarm_ParticleAt( swarm, particleIndex );
		if( Stg_Shape_IsCoordInside( self->_shape, particle->coord ) ) {
			IndexSet_Add( set, particleIndex );
		}	
	}

	return set;
}

VariableCondition_VariableIndex _SwarmShapeVC_GetVariableCount(void* variableCondition, Index globalIndex) {
	SwarmShapeVC*	self = (SwarmShapeVC*)variableCondition;

	return self->_entryCount;
}

Variable_Index _SwarmShapeVC_GetVariableIndex(void* variableCondition, Index globalIndex, VariableCondition_VariableIndex varIndex) {
	SwarmShapeVC*   self          = (SwarmShapeVC*)variableCondition;
	Variable_Index  searchedIndex = 0;
	Stream*         errorStr      = Journal_Register( Error_Type, self->type );
	Name            varName;

	Index		swarmVar_I;
	Swarm*		swarm		= self->_swarm;
	char*		swarmVarName;

	varName = self->_entryTbl[varIndex].varName;

	swarmVarName = (char*)calloc( strlen(swarm->name) + 1 + strlen(varName) + 1, sizeof(char) );
	strcat( swarmVarName, swarm->name );
	strcat( swarmVarName, "-" );
	strcat( swarmVarName, varName );
	/*searchedIndex = (Variable_Index)-1;
	for( swarmVar_I = 0; swarmVar_I < self->_swarm->nSwarmVars; swarmVar_I++ ) {
		if( swarm->swarmVars[swarmVar_I]->name && !strcmp(swarmVarName, swarm->swarmVars[swarmVar_I]->name) ) {
			searchedIndex = swarmVar_I;
		}
	}*/
	searchedIndex = Variable_Register_GetIndex(self->variable_Register, swarmVarName );
	
	Journal_Firewall( 
			( searchedIndex < self->variable_Register->count ),
			errorStr,
			"Error- in %s: searching for index of varIndex %u (\"%s\") at global node number %u failed"
			" - register returned index %u, greater than count %u.\n",
			__func__, varIndex, varName, globalIndex, searchedIndex, self->variable_Register->count );

	return searchedIndex; 
}


VariableCondition_ValueIndex _SwarmShapeVC_GetValueIndex(
		void*                                       variableCondition, 
		Index                                       globalIndex, 
		VariableCondition_VariableIndex             varIndex )
{
	return varIndex;
}


VariableCondition_ValueIndex _SwarmShapeVC_GetValueCount(void* variableCondition) {
	SwarmShapeVC*	self = (SwarmShapeVC*)variableCondition;
	
	return self->_entryCount;
}


VariableCondition_Value _SwarmShapeVC_GetValue(void* variableCondition, VariableCondition_ValueIndex valIndex) {
	SwarmShapeVC*	self = (SwarmShapeVC*)variableCondition;

	return self->_entryTbl[valIndex].value;
}

