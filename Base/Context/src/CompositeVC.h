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
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/** \file
** <b>Role:</b>
**
** <b>Assumptions:</b>
**
** <b>Comments:</b>
**
** $Id: CompositeVC.h 4297 2008-08-19 16:54:05Z LukeHodkinson $
**
**/

#ifndef __StGermain_Base_Context_CompositeVC_h__
#define __StGermain_Base_Context_CompositeVC_h__
	
	extern const Type CompositeVC_Type;
	
	#define __CompositeVC \
		/* General info */ \
		__VariableCondition \
		\
		/* Virtual info */ \
		\
		/* Stg_Class info */ \
		CompositeVC_ItemIndex	itemCount; \
		VariableCondition**		itemTbl; \
		int							nIndepItems;	\
   	VariableCondition**		indepItems; \
		Bool*							iOwnTbl; \
		SizeT							_size; \
		SizeT							_delta; \
		void*							data; \
		IndexSet**					attachedSets; \
		Bool							hasReadDictionary;
	
	struct _CompositeVC { __CompositeVC };
	
	/*--------------------------------------------------------------------------------------------------------------------------
	** Constructor
	*/
	
	VariableCondition* CompositeVC_Factory(
		AbstractContext*					context,
		Variable_Register*				variable_Register, 
		ConditionFunction_Register*	conFunc_Register,
		Dictionary*							dictionary,
		void*									data );
	
	CompositeVC* CompositeVC_New(
		Name									name,
		AbstractContext*					context,
		Variable_Register*				variable_Register, 
		ConditionFunction_Register*	conFunc_Register,
		Dictionary*							dictionary,
		void*									data );
	
	CompositeVC* _CompositeVC_DefaultNew( Name name );
	
	void CompositeVC_Init(
		CompositeVC*						self, 
		Name									name,
		Variable_Register*				variable_Register, 
		ConditionFunction_Register*	conFunc_Register,
		Dictionary*							dictionary,
		void*									data );
	
	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define COMPOSITEVC_DEFARGS \
                VARIABLECONDITION_DEFARGS

	#define COMPOSITEVC_PASSARGS \
                VARIABLECONDITION_PASSARGS

	CompositeVC* _CompositeVC_New(  COMPOSITEVC_DEFARGS  );
	
	void _CompositeVC_Init(
		void* compositeVC,
		void* data );
	
	/*--------------------------------------------------------------------------------------------------------------------------
	** General virtual functions
	*/
	
	void _CompositeVC_Delete( void* compositeVC );
	
	void _CompositeVC_Print( void* compositeVC, Stream* stream );
	
	void _CompositeVC_Destroy( void* compositeVC, void* data );

	/* Copy */
	#define CompositeVC_Copy( self ) \
		(VariableCondition*)Stg_Class_Copy( self, NULL, False, NULL, NULL )
	#define CompositeVC_Copy( self ) \
		(VariableCondition*)Stg_Class_Copy( self, NULL, False, NULL, NULL )
	
	void* _CompositeVC_Copy( void* compositeVC, void* dest, Bool deep, Name nameExt, struct PtrMap* ptrMap );
	
	/*--------------------------------------------------------------------------------------------------------------------------
	** Macros
	*/
	
	/*--------------------------------------------------------------------------------------------------------------------------
	** Virtual functions
	*/
	
	void _CompositeVC_Build( void* compositeVC, void* data );
	
	void _CompositeVC_AssignFromXML( void* compositeVC, Stg_ComponentFactory* cf, void* data );
	
	void _CompositeVC_ReadDictionary( void* compositeVC, void* dictionary );
	
	IndexSet* _CompositeVC_GetSet( void* compositeVC );
	
	VariableCondition_VariableIndex	_CompositeVC_GetVariableCount( void* compositeVC, Index globalIndex );
	
	Variable_Index _CompositeVC_GetVariableIndex(
		void*										compositeVC, 
		Index										globalIndex, 
		VariableCondition_VariableIndex	varIndex );
	
	VariableCondition_ValueIndex _CompositeVC_GetValueIndex(
		void*										compositeVC, 
		Index										globalIndex, 
		VariableCondition_VariableIndex	varIndex );
	
	VariableCondition_ValueIndex	_CompositeVC_GetValueCount( void* compositeVC );
	
	VariableCondition_Value _CompositeVC_GetValue( void* compositeVC, VariableCondition_ValueIndex valIndex );
	
	void _CompositeVC_PrintConcise( void* variableCondition, Stream* stream );
	
	/*--------------------------------------------------------------------------------------------------------------------------
	** Build functions
	*/
	
	CompositeVC_ItemIndex CompositeVC_Add( void* compositeVC, void* variableCondition, Bool iOwn );
	
	/*--------------------------------------------------------------------------------------------------------------------------
	** Functions
	*/

void _CompositeVC_Apply( void* _self, void* _ctx );

	
#endif /* __StGermain_Base_Context_CompositeVC_h__ */

