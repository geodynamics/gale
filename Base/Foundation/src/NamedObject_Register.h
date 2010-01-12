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
**  <b>Role:</b>
**	A generic register class for any Object.
**
** Assumptions:
**
** Comments:
**	Just uses a Stg_ObjectList to do the grunt work.
**
** $Id: NamedObject_Register.h 2428 2004-12-16 03:33:16Z AlanLo $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __StGermain_Base_Foundation_NamedObject_Register_h__
#define __StGermain_Base_Foundation_NamedObject_Register_h__ 

	extern const Type NamedObject_Register_Type;
	
	#define __NamedObject_Register \
		/* General info */ \
		__Stg_Class \
		\
		/* Virtual info */ \
		\
		/* Stg_Class info */ \
		Stg_ObjectList*	objects;
		

	struct NamedObject_Register { __NamedObject_Register };
	
	
	/* Stg_Class Administration members ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	
	
	NamedObject_Register*	NamedObject_Register_New( void );
	
	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define NAMEDOBJECT_REGISTER_DEFARGS \
                STG_CLASS_DEFARGS

	#define NAMEDOBJECT_REGISTER_PASSARGS \
                STG_CLASS_PASSARGS

	NamedObject_Register*	_NamedObject_Register_New(  NAMEDOBJECT_REGISTER_DEFARGS  );
	
	void _NamedObject_Register_Init( NamedObject_Register* self );
	
	void _NamedObject_Register_Delete( void* nameObjectRegister );
	
	void _NamedObject_Register_Print( void* nameObjectRegister, struct Stream* stream );
	
	void* _NamedObject_Register_Copy( void* namedObjectRegister, void* dest, Bool deep, Name nameExt, struct PtrMap* ptrMap );
	
	
	/* Public member functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	
   /* Deletes all elements from register and then
      deletes the register */
   void NamedObject_Register_DeleteAll( void* reg );
	
	#define NamedObject_Register_Add( self, nameObject ) \
		( Stg_ObjectList_Append( (self)->objects, nameObject ) )

	#define NamedObject_Register_GetIndex( self, name ) \
		( Stg_ObjectList_GetIndex( (self)->objects, name ) )

	#define NamedObject_Register_GetByName( self, name ) \
		( Stg_ObjectList_Get( (self)->objects, name ) )

	#define NamedObject_Register_GetByIndex( self, index ) \
		( (self)->objects->data[(index)] )
	
	#define NamedObject_Register_PrintAllEntryNames( self, stream ) \
		( Stg_ObjectList_PrintAllEntryNames( (self)->objects, stream ) )
	
	/* Private member functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
#endif /* __StGermain_Base_Foundation_NamedObject_Register_h__ */

