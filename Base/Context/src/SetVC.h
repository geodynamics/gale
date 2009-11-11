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
** Role:
**
** Assumptions:
**
** Comments:
**
** $Id: SetVC.h 4153 2007-07-26 02:25:22Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __Base_Automation_SetVC_h__
#define __Base_Automation_SetVC_h__
	
	extern const Type SetVC_Type;
	
	#define __SetVC_Entry \
		Name							varName; \
		VariableCondition_Value	value; \
		
	struct _SetVC_Entry { __SetVC_Entry };
	
	#define __SetVC \
		/* General info */ \
		__VariableCondition \
		\
		/* Virtual info */ \
		\
		/* Stg_Class info */ \
		Name					_dictionaryEntryName; \
		SetVC_Entry_Index	_entryCount; \
		SetVC_Entry*		_entryTbl; \
		IndexSet*			_vcset;

	struct _SetVC { __SetVC };

	#define SETVC_DEFARGS \
    	VARIABLECONDITION_DEFARGS, \
			Name _dictionaryEntryName 

	#define SETVC_PASSARGS \
    	VARIABLECONDITION_PASSARGS, \
			_dictionaryEntryName
	
	/*--------------------------------------------------------------------------------------------------------------------------
	** Constructor
	*/
	
	VariableCondition* SetVC_Factory(
		Variable_Register*				variable_Register, 
		ConditionFunction_Register*	conFunc_Register, 
		Dictionary*							dictionary,
		void*									data );
	
	SetVC* SetVC_New(
		Name									name,
		Name									_dictionaryEntryName, 
		Variable_Register*				variable_Register, 
		ConditionFunction_Register*	conFunc_Register,
		Dictionary*							dictionary );
	
	SetVC* _SetVC_DefaultNew( Name name );
	
	SetVC* _SetVC_New( SETVC_DEFARGS );
	
	void _SetVC_Init( void* setVC, Name _dictionaryEntryName );
	
	/*--------------------------------------------------------------------------------------------------------------------------
	** General virtual functions
	*/
	
	void _SetVC_Delete( void* setVC );
	
	void _SetVC_Print( void* setVC, Stream* stream );

	void _SetVC_Destroy( void* setVC, void* data );
	
	/* Copy */
	#define SetVC_Copy( self ) \
		(VariableCondition*)Stg_Class_Copy( self, NULL, False, NULL, NULL )
	#define SetVC_Copy( self ) \
		(VariableCondition*)Stg_Class_Copy( self, NULL, False, NULL, NULL )
	
	void* _SetVC_Copy( void* setVC, void* dest, Bool deep, Name nameExt, struct PtrMap* ptrMap );
	
	/*--------------------------------------------------------------------------------------------------------------------------
	** Macros
	*/
	
	/*--------------------------------------------------------------------------------------------------------------------------
	** Virtual functions
	*/
	
	void _SetVC_ReadDictionary( void* setVC, void* dictionary );
	
	IndexSet* _SetVC_GetSet( void* variableCondition );
	
	VariableCondition_VariableIndex _SetVC_GetVariableCount( void* variableCondition, Index globalIndex );
	
	Variable_Index _SetVC_GetVariableIndex(
		void*										variableCondition,
		Index										globalIndex, 
		VariableCondition_VariableIndex	varIndex);
						
	VariableCondition_ValueIndex _SetVC_GetValueIndex(
		void*										variableCondition, 
		Index										globalIndex, 
		VariableCondition_VariableIndex	varIndex);
						
	VariableCondition_ValueIndex _SetVC_GetValueCount( void* variableCondition );
	
	VariableCondition_Value _SetVC_GetValue( void* variableCondition, VariableCondition_ValueIndex valIndex );
	
	void _SetVC_PrintConcise( void* variableCondition, Stream* stream );
	
	/*--------------------------------------------------------------------------------------------------------------------------
	** Build functions
	*/
	
	/*--------------------------------------------------------------------------------------------------------------------------
	** Functions
	*/

#endif /* __Base_Automation_SetVC_h__ */
