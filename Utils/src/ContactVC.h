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
**  Role: Allows variable conditions to be defined on the walls of a regular mesh.
**
** Assumptions:
**
** Comments:
**
** $Id: ContactVC.h 4153 2007-07-26 02:25:22Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __Domain_Utils_ContactVC_h__
#define __Domain_Utils_ContactVC_h__
	

	extern const Type ContactVC_Type;
	
	
	#define __ContactVC \
		/* General info */ \
		__WallVC \
		\
		/* Virtual info */ \
		\
		/* Stg_Class info */ \
		Bool deep;

	struct _ContactVC { __ContactVC };
	
	
	/*--------------------------------------------------------------------------------------------------------------------------
	** Constructor
	*/
	
	VariableCondition*	ContactVC_Factory(
						Variable_Register*				variable_Register, 
						ConditionFunction_Register*			conFunc_Register, 
						Dictionary*					dictionary,
						void*						data );
	
	ContactVC*				ContactVC_DefaultNew( Name name );

	ContactVC*				ContactVC_New(
						Name						name,
						Name						_dictionaryEntryName, 
						Variable_Register*				variable_Register, 
						ConditionFunction_Register*			conFunc_Register, 
						Dictionary*					dictionary,
						void*						_mesh );
	
	void				ContactVC_Init(
						ContactVC*						self,
						Name						name,
						Name						_dictionaryEntryName, 
						Variable_Register*				variable_Register, 
						ConditionFunction_Register*			conFunc_Register, 
						Dictionary*					dictionary,
						void*						_mesh );
	
	ContactVC*				_ContactVC_New( 
						SizeT						_sizeOfSelf, 
						Type						type,
						Stg_Class_DeleteFunction*				_delete,
						Stg_Class_PrintFunction*				_print,
						Stg_Class_CopyFunction*				_copy,
						Stg_Component_DefaultConstructorFunction*	_defaultConstructor,
						Stg_Component_ConstructFunction*			_construct,
						Stg_Component_BuildFunction*			_build,
						Stg_Component_InitialiseFunction*			_initialise,
						Stg_Component_ExecuteFunction*			_execute,
						Stg_Component_DestroyFunction*			_destroy,
						Name						name, 
						Bool						initFlag,
						VariableCondition_BuildSelfFunc*		_buildSelf, 
						VariableCondition_PrintConciseFunc*		_printConcise,
						VariableCondition_ReadDictionaryFunc*		_readDictionary,
						VariableCondition_GetSetFunc*			_getSet,
						VariableCondition_GetVariableCountFunc*		_getVariableCount,
						VariableCondition_GetVariableIndexFunc*		_getVariableIndex,
						VariableCondition_GetValueIndexFunc*		_getValueIndex,
						VariableCondition_GetValueCountFunc*		_getValueCount,
						VariableCondition_GetValueFunc*			_getValue,
						VariableCondition_ApplyFunc*			_apply, 
						Name						_dictionaryEntryName, 
						Variable_Register*				variable_Register, 
						ConditionFunction_Register*			conFunc_Register, 
						Dictionary*					dictionary,
						void*						_mesh );
	
	void				_ContactVC_Init(
						void*						wallVC, 
						Name						_dictionaryEntryName, 
						void*						_mesh );
	
	
	/*--------------------------------------------------------------------------------------------------------------------------
	** General virtual functions
	*/
	
	void				_ContactVC_Delete( void* wallVC );
	
	void				_ContactVC_Print( void* wallVC, Stream* stream );
	
	/* Copy */
	#define ContactVC_Copy( self ) \
		(VariableCondition*)Stg_Class_Copy( self, NULL, False, NULL, NULL )
	#define ContactVC_Copy( self ) \
		(VariableCondition*)Stg_Class_Copy( self, NULL, False, NULL, NULL )
	
	void* _ContactVC_Copy( void* wallVC, void* dest, Bool deep, Name nameExt, struct PtrMap* ptrMap );
	
	void				_ContactVC_Build(  void* wallVC, void* data );
	
	
	/*--------------------------------------------------------------------------------------------------------------------------
	** Macros
	*/
	
	
	/*--------------------------------------------------------------------------------------------------------------------------
	** Virtual functions
	*/
	
	void				_ContactVC_Construct( void* wallVC, Stg_ComponentFactory* cf, void* data );
	
	void				_ContactVC_ReadDictionary( void* variableCondition, void* dictionary );
	
	IndexSet*			_ContactVC_GetSet( void* variableCondition );
	
	/*--------------------------------------------------------------------------------------------------------------------------
	** Build functions
	*/
	
	
	/*--------------------------------------------------------------------------------------------------------------------------
	** Functions
	*/

	
#endif /* __Domain_Utils_ContactVC_h__ */
