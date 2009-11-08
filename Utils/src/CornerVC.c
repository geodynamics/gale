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
** $Id: CornerVC.c 3310 2005-10-26 07:10:18Z RobertTurnbull $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>

#include <StgDomain/Geometry/Geometry.h>
#include <StgDomain/Shape/Shape.h>
#include <StgDomain/Mesh/Mesh.h>

#include "types.h"
#include "CornerVC.h"
#include "RegularMeshUtils.h"

#include <string.h>
#include <assert.h>



const Type CornerVC_Type = "CornerVC";
const Name defaultCornerVCName = "defaultCornerVCName";
/* List of corner names */
const char* CornerVC_CornerEnumToStr[CornerVC_Corner_Size] = {
   "bottomLeftFront",
   "bottomRightFront",
   "topLeftFront",
   "topRightFront",
   "bottomLeftBack",
   "bottomRightBack",
   "topLeftBack",
   "topRightBack" };


/*--------------------------------------------------------------------------------------------------------------------------
** Constructor
*/

VariableCondition* CornerVC_Factory(
   Variable_Register*				variable_Register, 
   ConditionFunction_Register*	conFunc_Register, 
   Dictionary*							dictionary,
   void*									data )
{
   return (VariableCondition*)CornerVC_New( defaultCornerVCName, NULL, variable_Register, conFunc_Register, dictionary, (Mesh*)data );
}

CornerVC* CornerVC_New(
   Name									name,
   Name									_dictionaryEntryName, 
   Variable_Register*				variable_Register, 
   ConditionFunction_Register*	conFunc_Register, 
   Dictionary*							dictionary,
   void*									_mesh )
{
   CornerVC* self = CornerVC_DefaultNew( name );

	self->isConstructed = True;
	_VariableCondition_Init( self, variable_Register, conFunc_Register, dictionary );
	_CornerVC_Init( self, _dictionaryEntryName, _mesh );

	return self;
}

CornerVC* CornerVC_DefaultNew( Name name ) {
	return _CornerVC_New(
		sizeof(CornerVC), 
		CornerVC_Type, 
		_CornerVC_Delete, 
		_CornerVC_Print, 
		_CornerVC_Copy,
		(Stg_Component_DefaultConstructorFunction*)CornerVC_DefaultNew,
		_CornerVC_AssignFromXML,	
		_CornerVC_Build,
		_VariableCondition_Initialise,
		_VariableCondition_Execute,
		_VariableCondition_Destroy,
		name,
		NON_GLOBAL,
		_CornerVC_BuildSelf, 
		_CornerVC_PrintConcise,
		_CornerVC_ReadDictionary,
		_CornerVC_GetSet, 
		_CornerVC_GetVariableCount, 
		_CornerVC_GetVariableIndex, 
		_CornerVC_GetValueIndex, 
		_CornerVC_GetValueCount, 
		_CornerVC_GetValue,
		_VariableCondition_Apply, 
		NULL,
		NULL, 
		NULL, 
		NULL, 
		NULL);
}

CornerVC* _CornerVC_New( CORNERVC_DEFARGS ) {
   CornerVC* self;
	
   /* Allocate memory/General info */
   assert(sizeOfSelf >= sizeof(CornerVC));
   self = (CornerVC*)_VariableCondition_New( VARIABLECONDITION_PASSARGS );
	
   /* Virtual info */
	
   /* Stg_Class info */
	
   return self;
}


void _CornerVC_Init(
   void*	cornerVC, 
   Name	_dictionaryEntryName, 
   void*	_mesh )
{
   CornerVC* self = (CornerVC*)cornerVC;

   self->isConstructed =         True;
   self->_dictionaryEntryName =  _dictionaryEntryName;
   self->_mesh =                 (Mesh*)_mesh;
   self->_corner =               CornerVC_Corner_Size;
   self->_entryTbl =             0;
   self->_entryCount =           0;
}


/*--------------------------------------------------------------------------------------------------------------------------
** General virtual functions
*/

void _CornerVC_ReadDictionary( void* variableCondition, void* dictionary ) {
   CornerVC*                 self = (CornerVC*)variableCondition;
   Dictionary_Entry_Value*   vcDictVal;
   Dictionary_Entry_Value    _vcDictVal;
   Dictionary_Entry_Value*   varsVal;
   CornerVC_Entry_Index      entry_I;
	
   /* Find dictionary entry */
   if (self->_dictionaryEntryName) {
      vcDictVal = Dictionary_Get(dictionary, self->_dictionaryEntryName);
   }
   else
   {
      vcDictVal = &_vcDictVal;
      Dictionary_Entry_Value_InitFromStruct(vcDictVal, dictionary);
   }

   if (vcDictVal)
   {
      char*	cornerStr;
	
      /* Obtain which corner */
      cornerStr = Dictionary_Entry_Value_AsString(Dictionary_Entry_Value_GetMember(vcDictVal, "corner" ));

      if (!strcasecmp(cornerStr, "bottomLeftFront")){
	 self->_corner = CornerVC_Corner_BottomLeftFront;
      }
      else if (!strcasecmp(cornerStr, "bottomRightFront"))
	 self->_corner = CornerVC_Corner_BottomRightFront;
      else if (!strcasecmp(cornerStr, "topLeftFront"))
	 self->_corner = CornerVC_Corner_TopLeftFront;
      else if (!strcasecmp(cornerStr, "topRightFront"))
	 self->_corner = CornerVC_Corner_TopRightFront;
      else if (!strcasecmp(cornerStr, "bottomLeftBack"))
	 self->_corner = CornerVC_Corner_BottomLeftBack;
      else if (!strcasecmp(cornerStr, "bottomRightBack"))
	 self->_corner = CornerVC_Corner_BottomRightBack;
      else if (!strcasecmp(cornerStr, "topLeftBack"))
	 self->_corner = CornerVC_Corner_TopLeftBack;
      else if (!strcasecmp(cornerStr, "topRightBack"))
	 self->_corner = CornerVC_Corner_TopRightBack;
      else {
	 assert( 0 );
	 self->_corner = CornerVC_Corner_Size; /* invalid entry */
      }
		
      /* Obtain the variable entries */
      self->_entryCount = Dictionary_Entry_Value_GetCount(Dictionary_Entry_Value_GetMember(vcDictVal, "variables"));
      self->_entryTbl = Memory_Alloc_Array( CornerVC_Entry, self->_entryCount, "CornerVC->_entryTbl" );
      varsVal = Dictionary_Entry_Value_GetMember(vcDictVal, "variables");
		
      for (entry_I = 0; entry_I < self->_entryCount; entry_I++)
      {
	 char*			valType;
	 Dictionary_Entry_Value*	valueEntry;
	 Dictionary_Entry_Value*	varDictListVal;
			
	 varDictListVal = Dictionary_Entry_Value_GetElement(varsVal, entry_I);
	 valueEntry = Dictionary_Entry_Value_GetMember(varDictListVal, "value");
			
	 self->_entryTbl[entry_I].varName = Dictionary_Entry_Value_AsString(
	    Dictionary_Entry_Value_GetMember(varDictListVal, "name"));
				
	 valType = Dictionary_Entry_Value_AsString(Dictionary_Entry_Value_GetMember(varDictListVal, "type"));
	 if (0 == strcasecmp(valType, "func"))
	 {
	    char*	funcName = Dictionary_Entry_Value_AsString(valueEntry);
	    Index	cfIndex;
				
	    self->_entryTbl[entry_I].value.type = VC_ValueType_CFIndex;
	    cfIndex = ConditionFunction_Register_GetIndex( self->conFunc_Register, funcName);
	    if ( cfIndex == (unsigned)-1 ) {	
	       Stream*	errorStr = Journal_Register( Error_Type, self->type );

	       Journal_Printf( errorStr, "Error- in %s: While parsing "
			       "definition of cornerVC \"%s\" (applies to corner \"%s\"), the cond. func. applied to "
			       "variable \"%s\" - \"%s\" - wasn't found in the c.f. register.\n",
			       __func__, self->_dictionaryEntryName, CornerVC_CornerEnumToStr[self->_corner],
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
										    self->_entryTbl[entry_I].value.as.typeArray.size, "CornerVC->_entryTbl[].value.as.typeArray.array" );
					
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
      self->_corner =     CornerVC_Corner_Size;
      self->_entryCount = 0;
      self->_entryTbl =   NULL;
   }
}


void _CornerVC_Delete(void* cornerVC)
{
   CornerVC*	self = (CornerVC*)cornerVC;
	
   if (self->_entryTbl) Memory_Free(self->_entryTbl);
	
   /* Stg_Class_Delete parent */
   _VariableCondition_Delete(self);
}


void _CornerVC_Print(void* cornerVC, Stream* stream)
{
   CornerVC*                   self = (CornerVC*)cornerVC;
   CornerVC_Entry_Index        entry_I;
   Index                       i;
	
   /* Set the Journal for printing informations */
   Stream* info = stream;
	
   /* General info */
   Journal_Printf( info, "CornerVC (ptr): %p\n", self);
	
   /* Virtual info */
	
   /* Stg_Class info */
   Journal_Printf( info, "\tdictionary (ptr): %p\n", self->dictionary);
   Journal_Printf( info, "\t_dictionaryEntryName (ptr): %p\n", self->_dictionaryEntryName);
   if (self->_dictionaryEntryName)
      Journal_Printf( info, "\t\t_dictionaryEntryName: %s\n", self->_dictionaryEntryName);
	
   Journal_Printf( info, "\t_corner: %s\n", self->_corner == CornerVC_Corner_BottomLeftFront ? "BottomLeftFront" :
		   self->_corner == CornerVC_Corner_BottomRightFront ? "BottomRightFront" : self->_corner == CornerVC_Corner_TopLeftFront ? "TopLeftFront" :
		   self->_corner == CornerVC_Corner_TopRightFront ? "TopRightFront" : self->_corner == CornerVC_Corner_BottomLeftBack ? "BottomLeftBack" :
		   self->_corner == CornerVC_Corner_BottomRightBack ? "BottomRightBack" : self->_corner == CornerVC_Corner_TopLeftBack ? "TopLeftBack" :
		   self->_corner == CornerVC_Corner_TopRightBack ? "TopRightBack" : "None");
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
	 }
      }
   Journal_Printf( info, "\t_mesh (ptr): %p\n", self->_mesh);
	
   /* Print parent */
   _VariableCondition_Print(self);
}


void* _CornerVC_Copy( void* cornerVC, void* dest, Bool deep, Name nameExt, struct PtrMap* ptrMap ) {
   CornerVC*       self = (CornerVC*)cornerVC;
   CornerVC*       newCornerVC;
   PtrMap*         map = ptrMap;
   Bool            ownMap = False;
	
   if( !map ) {
      map = PtrMap_New( 10 );
      ownMap = True;
   }
	
   newCornerVC = (CornerVC*)_VariableCondition_Copy( self, dest, deep, nameExt, map );
	
   newCornerVC->_dictionaryEntryName = self->_dictionaryEntryName;
   newCornerVC->_corner = self->_corner;
   newCornerVC->_entryCount = self->_entryCount;
	
   if( deep ) {
      newCornerVC->_mesh = (Mesh*)Stg_Class_Copy( self->_mesh, NULL, deep, nameExt, map );
		
      if( (newCornerVC->_entryTbl = PtrMap_Find( map, self->_entryTbl )) == NULL && self->_entryTbl ) {
	 newCornerVC->_entryTbl = Memory_Alloc_Array( CornerVC_Entry, newCornerVC->_entryCount, "CornerVC->_entryTbl");
	 memcpy( newCornerVC->_entryTbl, self->_entryTbl, sizeof(CornerVC_Entry) * newCornerVC->_entryCount );
	 PtrMap_Append( map, newCornerVC->_entryTbl, self->_entryTbl );
      }
   }
   else {
      newCornerVC->_mesh = self->_mesh;
      newCornerVC->_entryTbl = self->_entryTbl;
   }
	
   if( ownMap ) {
      Stg_Class_Delete( map );
   }
	
   return (void*)newCornerVC;
}


void _CornerVC_Build(  void* cornerVC, void* data ) {
   CornerVC*         self = (CornerVC*)cornerVC;
	
   _CornerVC_BuildSelf( self, data );
	
   _VariableCondition_Build( self, data );
}
	

/*--------------------------------------------------------------------------------------------------------------------------
** Macros
*/


/*--------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _CornerVC_AssignFromXML( void* cornerVC, Stg_ComponentFactory* cf, void* data )
{
	
}

void _CornerVC_BuildSelf(  void* cornerVC, void* data ) {
   CornerVC*        self = (CornerVC*)cornerVC;
	
   if( self->_mesh ) {
      Stg_Component_Build( self->_mesh, data, False );
   }
}

IndexSet* _CornerVC_GetSet(void* variableCondition)
{
   CornerVC*   self = (CornerVC*)variableCondition;
   IndexSet    *set = NULL;
   Stream*     warningStr = Journal_Register( Error_Type, self->type );
   unsigned	nDims;
   Grid*	vertGrid;

   nDims = Mesh_GetDimSize( self->_mesh );
	vertGrid = *(Grid**)ExtensionManager_Get( self->_mesh->info, self->_mesh, 
						  ExtensionManager_GetHandle( self->_mesh->info, 
									      "vertexGrid" ) );

   switch (self->_corner) {
      case CornerVC_Corner_BottomLeftFront:
	 if ( nDims < 3 || !vertGrid->sizes[2] ) {
	    Journal_Printf( warningStr, "Warning - in %s: Can't build a %s corner VC "
			    "when mesh has no elements in the %s axis. Returning an empty set.\n", __func__,
			    CornerVC_CornerEnumToStr[self->_corner], "K" );
	    set = IndexSet_New( Mesh_GetDomainSize( self->_mesh, MT_VERTEX ) );
	 }
	 else {
	    set = RegularMeshUtils_CreateGlobalBottomLeftFrontSet(self->_mesh);
	 }
	 break;
			
      case CornerVC_Corner_BottomRightFront:
	 if ( !vertGrid->sizes[0] ) {
	    Journal_Printf( warningStr, "Warning - in %s: Can't build a %s corner VC "
			    "when mesh has no elements in the %s axis. Returning an empty set.\n", __func__,
			    CornerVC_CornerEnumToStr[self->_corner], "I" );
	    set = IndexSet_New( Mesh_GetDomainSize( self->_mesh, MT_VERTEX ) );	
	 }
	 else if ( nDims < 3 || !vertGrid->sizes[2] ) {
	    Journal_Printf( warningStr, "Warning - in %s: Can't build a %s corner VC "
			    "when mesh has no elements in the %s axis. Returning an empty set.\n", __func__,
			    CornerVC_CornerEnumToStr[self->_corner], "K" );
	    set = IndexSet_New( Mesh_GetDomainSize( self->_mesh, MT_VERTEX ) );	
	 }
	 else {
	    set = RegularMeshUtils_CreateGlobalBottomRightFrontSet(self->_mesh);
	 }	
	 break;
			
      case CornerVC_Corner_TopLeftFront:
	 if ( nDims < 2 || !vertGrid->sizes[1] ) {
	    Journal_Printf( warningStr, "Warning - in %s: Can't build a %s corner VC "
			    "when mesh has no elements in the %s axis. Returning an empty set.\n", __func__,
			    CornerVC_CornerEnumToStr[self->_corner], "J" );
	    set = IndexSet_New( Mesh_GetDomainSize( self->_mesh, MT_VERTEX ));	
	 }
	 else if ( nDims < 3 || !vertGrid->sizes[2] ) {
	    Journal_Printf( warningStr, "Warning - in %s: Can't build a %s corner VC "
			    "when mesh has no elements in the %s axis. Returning an empty set.\n", __func__,
			    CornerVC_CornerEnumToStr[self->_corner], "K" );
	    set = IndexSet_New( Mesh_GetDomainSize( self->_mesh, MT_VERTEX ) );	
	 }
	 else {
	    set = RegularMeshUtils_CreateGlobalTopLeftFrontSet(self->_mesh);
	 }	
	 break;
			
      case CornerVC_Corner_TopRightFront:
	 if ( nDims < 2 || !vertGrid->sizes[1] ) {
	    Journal_Printf( warningStr, "Warning - in %s: Can't build a %s corner VC "
			    "when mesh has no elements in the %s axis. Returning an empty set.\n", __func__,
			    CornerVC_CornerEnumToStr[self->_corner], "J" );
	    set = IndexSet_New( Mesh_GetDomainSize( self->_mesh, MT_VERTEX ) );	
	 }
	 else if ( nDims < 3 || !vertGrid->sizes[2] ) {
	    Journal_Printf( warningStr, "Warning - in %s: Can't build a %s corner VC "
			    "when mesh has no elements in the %s axis. Returning an empty set.\n", __func__,
			    CornerVC_CornerEnumToStr[self->_corner], "K" );
	    set = IndexSet_New( Mesh_GetDomainSize( self->_mesh, MT_VERTEX ) );	
	 }
	 else if ( !vertGrid->sizes[0] ) {
	    Journal_Printf( warningStr, "Warning - in %s: Can't build a %s corner VC "
			    "when mesh has no elements in the %s axis. Returning an empty set.\n", __func__,
			    CornerVC_CornerEnumToStr[self->_corner], "I" );
	    set = IndexSet_New( Mesh_GetDomainSize( self->_mesh, MT_VERTEX ) );	
	 }
	 else {
	    set = RegularMeshUtils_CreateGlobalTopRightFrontSet(self->_mesh);
	 }	
	 break;
			
      case CornerVC_Corner_BottomLeftBack:
	 set = RegularMeshUtils_CreateGlobalBottomLeftBackSet(self->_mesh);
	 break;
			
      case CornerVC_Corner_BottomRightBack:
	 if ( !vertGrid->sizes[0] ) {
	    Journal_Printf( warningStr, "Warning - in %s: Can't build a %s corner VC "
			    "when mesh has no elements in the %s axis. Returning an empty set.\n", __func__,
			    CornerVC_CornerEnumToStr[self->_corner], "I" );
	    set = IndexSet_New( Mesh_GetDomainSize( self->_mesh, MT_VERTEX ) );	
	 }
	 else {
	    set = RegularMeshUtils_CreateGlobalBottomRightBackSet(self->_mesh);
	 }
	 break;
			
      case CornerVC_Corner_TopLeftBack:
	 if ( nDims < 2 || !vertGrid->sizes[1] ) {
	    Journal_Printf( warningStr, "Warning - in %s: Can't build a %s corner VC "
			    "when mesh has no elements in the %s axis. Returning an empty set.\n", __func__,
			    CornerVC_CornerEnumToStr[self->_corner], "J" );
	    set = IndexSet_New( Mesh_GetDomainSize( self->_mesh, MT_VERTEX ) );	
	 }
	 else {
	    set = RegularMeshUtils_CreateGlobalTopLeftBackSet(self->_mesh);
	 }	
	 break;
			
      case CornerVC_Corner_TopRightBack:
	 if ( !vertGrid->sizes[0] ) {
	    Journal_Printf( warningStr, "Warning - in %s: Can't build a %s corner VC "
			    "when mesh has no elements in the %s axis. Returning an empty set.\n", __func__,
			    CornerVC_CornerEnumToStr[self->_corner], "I" );
	    set = IndexSet_New(Mesh_GetDomainSize( self->_mesh, MT_VERTEX ) );	
	 }
	 else if ( nDims < 2 || !vertGrid->sizes[1] ) {
	    Journal_Printf( warningStr, "Warning - in %s: Can't build a %s corner VC "
			    "when mesh has no elements in the %s axis. Returning an empty set.\n", __func__,
			    CornerVC_CornerEnumToStr[self->_corner], "J" );
	    set = IndexSet_New(Mesh_GetDomainSize( self->_mesh, MT_VERTEX ));
	 }
	 else {
	    set = RegularMeshUtils_CreateGlobalTopRightBackSet(self->_mesh);
	 }
	 break;
		
      case CornerVC_Corner_Size:
      default:
	 assert(0);
	 break;
   }

   return set;
}


VariableCondition_VariableIndex _CornerVC_GetVariableCount(void* variableCondition, Index globalIndex)
{
   CornerVC*   self = (CornerVC*)variableCondition;
	
   return self->_entryCount;
}


Variable_Index _CornerVC_GetVariableIndex(void* variableCondition, Index globalIndex, VariableCondition_VariableIndex varIndex)
{
   CornerVC*       self = (CornerVC*)variableCondition;
   Variable_Index  searchedIndex = 0;
   Stream*         errorStr = Journal_Register( Error_Type, self->type );
   Name            varName;
	
   varName = self->_entryTbl[varIndex].varName;
   searchedIndex = Variable_Register_GetIndex(self->variable_Register, varName );
	
   Journal_Firewall( ( searchedIndex < self->variable_Register->count ), errorStr, "Error- in %s: searching for index of "
		     "varIndex %u (\"%s\") at global node number %u failed - register returned index %u, greater than "
		     "count %u.\n", __func__, varIndex, varName, globalIndex, searchedIndex, self->variable_Register->count );

   return searchedIndex; 
}


VariableCondition_ValueIndex _CornerVC_GetValueIndex(
   void*				variableCondition, 
   Index				globalIndex, 
   VariableCondition_VariableIndex	varIndex)
{
   return varIndex;
}


VariableCondition_ValueIndex _CornerVC_GetValueCount(void* variableCondition)
{
   CornerVC*	self = (CornerVC*)variableCondition;
	
   return self->_entryCount;
}


VariableCondition_Value _CornerVC_GetValue(void* variableCondition, VariableCondition_ValueIndex valIndex)
{
   CornerVC*	self = (CornerVC*)variableCondition;

   return self->_entryTbl[valIndex].value;
}

void _CornerVC_PrintConcise( void* variableCondition, Stream* stream ) {
   CornerVC*	self = (CornerVC*)variableCondition;
	
   Journal_Printf( stream, "\ttype: %s, set: ", self->type );
   Journal_Printf( stream, "%s\n", 
		   self->_corner == CornerVC_Corner_BottomLeftFront ? "BottomLeftFront" :
		   self->_corner == CornerVC_Corner_BottomRightFront ? "BottomRightFront" : 
		   self->_corner == CornerVC_Corner_TopLeftFront ? "TopLeftFront" :
		   self->_corner == CornerVC_Corner_TopRightFront ? "TopRightFront" : 
		   self->_corner == CornerVC_Corner_BottomLeftBack ? "BottomLeftBack" :
		   self->_corner == CornerVC_Corner_BottomRightBack ? "BottomRightBack" : 
		   self->_corner == CornerVC_Corner_TopLeftBack ? "TopLeftBack" :
		   self->_corner == CornerVC_Corner_TopRightBack ? "TopRightBack" : "None" );
}

/*--------------------------------------------------------------------------------------------------------------------------
** Build functions
*/


/*--------------------------------------------------------------------------------------------------------------------------
** Functions
*/
