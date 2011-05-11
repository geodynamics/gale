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

#ifndef __StgDomain_Utils_ContactVC_h__
#define __StgDomain_Utils_ContactVC_h__
	

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
	
	VariableCondition* ContactVC_Factory(
		AbstractContext*					context,
		Variable_Register*				variable_Register, 
		ConditionFunction_Register*	conFunc_Register, 
		Dictionary*							dictionary,
		void*									data );
	
	ContactVC* _ContactVC_DefaultNew( Name name );

	ContactVC* ContactVC_New(
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

	#define CONTACTVC_DEFARGS \
                WALLVC_DEFARGS

	#define CONTACTVC_PASSARGS \
                WALLVC_PASSARGS

	ContactVC* _ContactVC_New(  CONTACTVC_DEFARGS  );
	
	void _ContactVC_Init( void* wallVC, Name _dictionaryEntryName, void* _mesh );
	
	
	/*--------------------------------------------------------------------------------------------------------------------------
	** General virtual functions
	*/
	
	void _ContactVC_Delete( void* wallVC );
	
	void _ContactVC_Destroy( void* wallVC, void* data );

	void _ContactVC_Print( void* wallVC, Stream* stream );
	
	/* Copy */
	#define ContactVC_Copy( self ) \
		(VariableCondition*)Stg_Class_Copy( self, NULL, False, NULL, NULL )
	#define ContactVC_Copy( self ) \
		(VariableCondition*)Stg_Class_Copy( self, NULL, False, NULL, NULL )
	
	void* _ContactVC_Copy( const void* wallVC, void* dest, Bool deep, Name nameExt, struct PtrMap* ptrMap );
	
	void _ContactVC_Build(  void* wallVC, void* data );
	
	
	/*--------------------------------------------------------------------------------------------------------------------------
	** Macros
	*/
	
	
	/*--------------------------------------------------------------------------------------------------------------------------
	** Virtual functions
	*/
	
	void _ContactVC_AssignFromXML( void* wallVC, Stg_ComponentFactory* cf, void* data );
	
	void _ContactVC_ReadDictionary( void* variableCondition, void* dictionary );
	
	IndexSet* _ContactVC_GetSet( void* variableCondition );
	
	/*--------------------------------------------------------------------------------------------------------------------------
	** Build functions
	*/
	
	
	/*--------------------------------------------------------------------------------------------------------------------------
	** Functions
	*/

	
#endif /* __StgDomain_Utils_ContactVC_h__ */

