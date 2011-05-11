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
** $Id: CornerVC.h 3291 2005-10-18 00:05:33Z AlanLo $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __StgDomain_Utils_CornerVC_h__
#define __StgDomain_Utils_CornerVC_h__
	

	extern const Type CornerVC_Type;

	extern Name CornerVC_CornerEnumToStr[CornerVC_Corner_Size];
	
	#define __CornerVC_Entry \
		Name							varName; \
		VariableCondition_Value	value; \
		
	struct _CornerVC_Entry { __CornerVC_Entry };
	
	
	#define __CornerVC \
		/* General info */ \
		__VariableCondition \
		\
		/* Virtual info */ \
		\
		/* Stg_Class info */ \
		Name						_dictionaryEntryName; \
		CornerVC_Corner		_corner; \
		CornerVC_Entry_Index	_entryCount; \
		CornerVC_Entry*		_entryTbl; \
		Mesh*						_mesh;

	struct _CornerVC { __CornerVC };

	/*--------------------------------------------------------------------------------------------------------------------------
	** Constructor
	*/
	
	VariableCondition* CornerVC_Factory(
		AbstractContext*					context,
		Variable_Register*				variable_Register, 
		ConditionFunction_Register*	conFunc_Register, 
		Dictionary*							dictionary,
		void*									data );
	
	CornerVC* _CornerVC_DefaultNew( Name name );

	CornerVC* CornerVC_New(
		Name									name,
		AbstractContext*					context,
		Name									_dictionaryEntryName, 
		Variable_Register*				variable_Register, 
		ConditionFunction_Register*	conFunc_Register, 
		Dictionary*							dictionary,
		void*									_mesh );
	
	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define CORNERVC_DEFARGS \
                VARIABLECONDITION_DEFARGS

	#define CORNERVC_PASSARGS \
                VARIABLECONDITION_PASSARGS

	CornerVC* _CornerVC_New(  CORNERVC_DEFARGS  );
	
	void _CornerVC_Init( void* cornerVC, Name _dictionaryEntryName, void* _mesh );
	
	
	/*--------------------------------------------------------------------------------------------------------------------------
	** General virtual functions
	*/
	
	void _CornerVC_Delete( void* cornerVC );
	
	void _CornerVC_Print( void* cornerVC, Stream* stream );
	
	void _CornerVC_Destroy( void* cornerVC, void* data );
	
	/* Copy */
	#define CornerVC_Copy( self ) \
		(VariableCondition*)Stg_Class_Copy( self, NULL, False, NULL, NULL )
	#define CornerVC_Copy( self ) \
		(VariableCondition*)Stg_Class_Copy( self, NULL, False, NULL, NULL )
	
	void* _CornerVC_Copy( const void* cornerVC, void* dest, Bool deep, Name nameExt, struct PtrMap* ptrMap );
	
	void _CornerVC_Build(  void* cornerVC, void* data );
	
	
	/*--------------------------------------------------------------------------------------------------------------------------
	** Macros
	*/
	
	
	/*--------------------------------------------------------------------------------------------------------------------------
	** Virtual functions
	*/
	
	void _CornerVC_AssignFromXML( void* cornerVC, Stg_ComponentFactory* cf, void* data );
	
	void _CornerVC_BuildSelf( void* cornerVC, void* data );
	
	void _CornerVC_ReadDictionary( void* variableCondition, void* dictionary );
	
	IndexSet* _CornerVC_GetSet( void* variableCondition );
	
	VariableCondition_VariableIndex _CornerVC_GetVariableCount( void* variableCondition, Index globalIndex );
	
	Variable_Index _CornerVC_GetVariableIndex( void* variableCondition, Index globalIndex, VariableCondition_VariableIndex varIndex );
						
	VariableCondition_ValueIndex _CornerVC_GetValueIndex( void* variableCondition, Index globalIndex, VariableCondition_VariableIndex varIndex );
						
	VariableCondition_ValueIndex _CornerVC_GetValueCount( void* variableCondition );
	
	VariableCondition_Value _CornerVC_GetValue( void* variableCondition, VariableCondition_ValueIndex valIndex );
	
	void _CornerVC_PrintConcise( void* variableCondition, Stream* stream );
	
	/*--------------------------------------------------------------------------------------------------------------------------
	** Build functions
	*/
	
	
	/*--------------------------------------------------------------------------------------------------------------------------
	** Functions
	*/

	
#endif /* __StgDomain_Utils_CornerVC_h__ */

