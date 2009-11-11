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
** $Id: DynamicVC.h 3462 2006-02-19 06:53:24Z WalterLandry $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __Base_Automation_DynamicVC_h__
#define __Base_Automation_DynamicVC_h__

extern const Type DynamicVC_Type;

#define __DynamicVC				\
	/* General info */			\
	__VariableCondition			\
						\
	/* Virtual info */			\
						\
	/* Stg_Class info */			\
	Variable*	var;			\
	IMap*		vcMap;

	struct DynamicVC { __DynamicVC };

	#define DYNAMICVC_DEFARGS \
    	VARIABLECONDITION_DEFARGS

	#define DYNAMICVC_PASSARGS \
    	VARIABLECONDITION_PASSARGS

	/*--------------------------------------------------------------------------------------------------------------------------
	** Constructor
	*/

	VariableCondition* DynamicVC_Factory( Variable_Register* varReg, 
		ConditionFunction_Register*	conFuncReg, 
		Dictionary*							dict, 
		void*									data );

	DynamicVC* DynamicVC_New( Name name,
		Variable_Register*				variable_Register, 
		ConditionFunction_Register*	conFunc_Register,
		Dictionary*							dict );

	DynamicVC* _DynamicVC_DefaultNew( Name name );

	DynamicVC* _DynamicVC_New( DYNAMICVC_DEFARGS );

	void _DynamicVC_Init( void* vc );

	/*--------------------------------------------------------------------------------------------------------------------------
	** General virtual functions
	*/
	
	void _DynamicVC_Delete( void* vc );
	void _DynamicVC_Print( void* vc, Stream* stream );

	/* Copy */
	#define DynamicVC_Copy( self )						\
		(VariableCondition*)Stg_Class_Copy( self, NULL, False, NULL, NULL )
	#define DynamicVC_Copy( self )						\
		(VariableCondition*)Stg_Class_Copy( self, NULL, False, NULL, NULL )

	void* _DynamicVC_Copy( void* vc, void* dest, Bool deep, Name nameExt, struct PtrMap* ptrMap );

	void _DynamicVC_AssignFromXML( void* vc, Stg_ComponentFactory* cf, void* data );

	void _DynamicVC_Build( void* vc, void* data );

	void _DynamicVC_Initialise( void* vc, void* data );

	void _DynamicVC_Execute( void* vc, void* data );

	void _DynamicVC_Destroy( void* vc, void* data );

	void _DynamicVC_ReadDictionary( void* vc, void* dict );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Macros
	*/

	/*--------------------------------------------------------------------------------------------------------------------------
	** Virtual functions
	*/

	void _DynamicVC_ReadDictionary( void* vc, void* dict );

	IndexSet* _DynamicVC_GetSet( void* vc );

	VariableCondition_VariableIndex _DynamicVC_GetVariableCount( void* vc, Index globalIndex );

	Variable_Index _DynamicVC_GetVariableIndex( void* vc,
		Index globalIndex, 
		VariableCondition_VariableIndex varIndex );

	VariableCondition_ValueIndex _DynamicVC_GetValueIndex( void* vc, 
		Index globalIndex, 
		VariableCondition_VariableIndex varIndex );

	VariableCondition_ValueIndex _DynamicVC_GetValueCount( void* vc );

	VariableCondition_Value _DynamicVC_GetValue( void* vc, VariableCondition_ValueIndex valIndex );

	void _DynamicVC_PrintConcise( void* vc, Stream* stream );

	void DynamicVC_Apply( void* vc, void* ctx );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Build functions
	*/

	/*--------------------------------------------------------------------------------------------------------------------------
	** Functions
	*/

	void DynamicVC_SetVariable( void* vc, Variable* var );
	void DynamicVC_SetMaxEntries( void* vc, int maxEntries );
	int DynamicVC_GetMaxEntries( void* vc );
	void DynamicVC_SetValues( void* vc, int nVals, VariableCondition_Value* vals );
	void DynamicVC_Insert( void* vc, int index, int valIndex );
	void DynamicVC_Remove( void* vc, int index );
	Bool DynamicVC_Has( void* _self, int index );

#endif /* __Base_Automation_DynamicVC_h__ */
