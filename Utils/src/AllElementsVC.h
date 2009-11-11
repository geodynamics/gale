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
*/
/** \file
**  Role:
**
** Assumptions:
**
** Comments:
**
** $Id: AllElementsVC.h 4153 2007-07-26 02:25:22Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __Domain_Utils_AllElementsVC_h__
#define __Domain_Utils_AllElementsVC_h__
	

	extern const Type AllElementsVC_Type;
	
	#define __AllElementsVC_Entry \
		Name							varName; \
		VariableCondition_Value	value; \
		
	struct _AllElementsVC_Entry { __AllElementsVC_Entry };
	
	#define __AllElementsVC \
		/* General info */ \
		__VariableCondition \
		\
		/* Virtual info */ \
		\
		/* Stg_Class info */ \
		Name								_dictionaryEntryName; \
		AllElementsVC_Entry_Index	_entryCount; \
		AllElementsVC_Entry*			_entryTbl; \
		Mesh*								mesh;

	struct _AllElementsVC { __AllElementsVC };
	
	#define ALLELEMENTSVC_DEFARGS \
      VARIABLECONDITION_DEFARGS, \
         Name  _dictionaryEntryName, \
         void* _mesh 

   #define ALLELEMENTSVC_PASSARGS \
      VARIABLECONDITION_PASSARGS, \
         _dictionaryEntryName, \
         _mesh 
	
	/*--------------------------------------------------------------------------------------------------------------------------
	** Constructor
	*/
	
	VariableCondition* AllElementsVC_Factory(
		Variable_Register*				variable_Register, 
		ConditionFunction_Register*	conFunc_Register, 
		Dictionary*							dictionary,
		void*									data );
	
	AllElementsVC*	AllElementsVC_DefaultNew( Name name );
	
	AllElementsVC* AllElementsVC_New(
		Name									name,
		Name									_dictionaryEntryName, 
		Variable_Register*				variable_Register, 
		ConditionFunction_Register*	conFunc_Register,
		Dictionary*							dictionary,
		void*									mesh );
	
	AllElementsVC* _AllElementsVC_New( ALLELEMENTSVC_DEFARGS );
	
	void _AllElementsVC_Init( void* allElementsVC, Name _dictionaryEntryName, void* mesh );
	
	/*--------------------------------------------------------------------------------------------------------------------------
	** General virtual functions
	*/
	
	void _AllElementsVC_AssignFromXML( void* allElementsVC, Stg_ComponentFactory *cf, void* data );
	
	void _AllElementsVC_Delete( void* allElementsVC );
	
	void _AllElementsVC_Print( void* allElementsVC, Stream* stream );
	
	void _AllElementsVC_Destroy( void* allElementsVC, void* data );

	/* Copy */
	#define AllElementsVC_Copy( self ) \
		(VariableCondition*)Stg_Class_Copy( self, NULL, False, NULL, NULL )
	#define AllElementsVC_Copy( self ) \
		(VariableCondition*)Stg_Class_Copy( self, NULL, False, NULL, NULL )
	
	void* _AllElementsVC_Copy( void* allElementsVC, void* dest, Bool deep, Name nameExt, struct PtrMap* ptrMap );
	
	void _AllElementsVC_Build( void* allElementsVC, void* data );
	
	
	/*--------------------------------------------------------------------------------------------------------------------------
	** Macros
	*/
	
	
	/*--------------------------------------------------------------------------------------------------------------------------
	** Virtual functions
	*/
	
	void _AllElementsVC_BuildSelf( void* allElementsVC, void* data );
	
	void _AllElementsVC_ReadDictionary( void* variableCondition, void* dictionary );
	
	IndexSet* _AllElementsVC_GetSet( void* variableCondition );
	
	VariableCondition_VariableIndex _AllElementsVC_GetVariableCount( void* variableCondition, Index globalIndex );
	
	Variable_Index _AllElementsVC_GetVariableIndex( void* variableCondition, Index globalIndex, VariableCondition_VariableIndex varIndex );
						
	VariableCondition_ValueIndex _AllElementsVC_GetValueIndex( void* variableCondition, Index globalIndex, VariableCondition_VariableIndex varIndex );
						
	VariableCondition_ValueIndex _AllElementsVC_GetValueCount( void* variableCondition );
	
	VariableCondition_Value _AllElementsVC_GetValue( void* variableCondition, VariableCondition_ValueIndex valIndex );
	
	void _AllElementsVC_PrintConcise( void* variableCondition, Stream* stream );
	
	/*--------------------------------------------------------------------------------------------------------------------------
	** Build functions
	*/
	
	
	/*--------------------------------------------------------------------------------------------------------------------------
	** Functions
	*/

	
#endif /* __Domain_Utils_AllElementsVC_h__ */
