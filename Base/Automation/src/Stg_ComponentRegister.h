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
**	Abstract class for objects that are named.
**
** Assumptions:
**
** Comments:
**
** $Id: Stg_ComponentRegister.h 2745 2005-05-11 08:12:18Z RaquibulHassan $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __Base_Automation_Stg_ComponentRegister_h__
#define __Base_Automation_Stg_ComponentRegister_h__
	
	/* Textual name of this class */
	extern const Type Stg_ComponentRegister_Type;
	extern Stg_ComponentRegister *stgComponentRegister;

	/*struct Stg_Component_DefaultConstructorFunction;*/
	#define __Stg_ComponentRegisterElement \
		__Stg_Object						\
		Type								componentType; \
		Stg_Component_DefaultConstructorFunction*			defaultConstructor; \
		Stg_Component_MetaAsDictionaryFunction*                         metadata; \
		Name								version;

	struct Stg_ComponentRegisterElement{ __Stg_ComponentRegisterElement };

	extern const Type Stg_ComponentRegisterElement_Type;

	/** ComponentRegisterElement Constructor interface. */
	Stg_ComponentRegisterElement* Stg_ComponentRegisterElement_New(
		Type			type,
		Stg_Class_DeleteFunction*	_delete,
		Stg_Class_PrintFunction*	_print,
      Type        componentType,
		Stg_Component_DefaultConstructorFunction*		defaultConstructor,
		Stg_Component_MetaAsDictionaryFunction*      metadata,
		Name								version
      );
	
   	/** Stg_Class_Delete interface. */
   	void _Stg_ComponentRegisterElement_Delete( void* element );

   	/** Print interaface. */
	   void _Stg_ComponentRegisterElement_Print( void* element, Stream* paramStream );	

	/* Stg_ComponentRegister information */
	#define __Stg_ComponentRegister \
		/* General info */ \
		__Stg_Class \
		\
		/* Virtual info */ \
		\
		/* Stg_ComponentRegister info */ \
		Stg_ObjectList*									constructors;  \
      Stream*                                   debugStream;   \
	
	struct Stg_ComponentRegister { __Stg_ComponentRegister };
	
	/* Class Administration members ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	
	/** Constructor Implementation */
	Stg_ComponentRegister *_Stg_ComponentRegister_New(
		SizeT					_sizeOfSelf, 
		Type					type,
		Stg_Class_DeleteFunction*		_delete,
		Stg_Class_PrintFunction*		_print, 
		Stg_Class_CopyFunction*			_copy );
	
	Stg_ComponentRegister *Stg_ComponentRegister_New(  );
	
	/* Initialisation implementation */
	void _Stg_ComponentRegister_Init( Stg_ComponentRegister* self );
	
	void Stg_ComponentRegister_Init( Stg_ComponentRegister* self );
	
	/* Delete boundary condition layout implementation */
	void _Stg_ComponentRegister_Delete( void* componentRegister );
	
	/* Print boundary condition layout implementation */
	void _Stg_ComponentRegister_Print( void* componentRegister, Stream* stream );
	
	
	/* Private member functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	int constructorElementCompareFunction( void *data1, void *data2 );
	
	int constructorElementCompareFunction1( void *data1, void *data2 );
	
	void constructorElementPrintFunction( void *nodeData, Stream *printStream );

	void constructorElementDeleteFunction( void *nodeData );
	
	/* Public member functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	int Stg_ComponentRegister_AddFunc( 
			Stg_ComponentRegister *self,
			Name componentType,
			Name version,
			Stg_Component_DefaultConstructorFunction *func,
			Stg_Component_MetaAsDictionaryFunction* metadata );

	/* Adds a Component to the database/register. It is a macro because of the auto argument generation based
		ont the "componentType" argument, to load in the associated meta data function */
	#define Stg_ComponentRegister_Add( self, componentType, version, func ) \
	{ \
		Dictionary* componentType ##_MetaAsDictionary(); \
		Dictionary* componentType ##_Type_MetaAsDictionary(); \
		Stg_ComponentRegister_AddFunc( self, componentType, version, func, componentType ##_MetaAsDictionary ); \
	}

   /* Remove and free a component in the register */
   Bool Stg_ComponentRegister_RemoveEntry(
      Stg_ComponentRegister* self,
      Name                   componentType,
      Name                   version );

	Stg_Component_DefaultConstructorFunction* Stg_ComponentRegister_Get( 
			Stg_ComponentRegister* self,
			Name                   componentType,
			Name                   version );
	
	/* Same function as above, but it asserts and gives a nice message if it cannot find the default constructor */
	Stg_Component_DefaultConstructorFunction* Stg_ComponentRegister_AssertGet( 
			Stg_ComponentRegister* self,
			Name                   componentType,
			Name                   version ); 

	Dictionary* Stg_ComponentRegister_GetMetadata(
			Stg_ComponentRegister* self,
			Name                   componentType,
			Name                   version );

	/** This function returns the pointer to the singleton "stgComponentRegister" */
	Stg_ComponentRegister *Stg_ComponentRegister_Get_ComponentRegister( );

	/** This function prints the types registered that are most similar to the 'name' passed in */
	void Stg_ComponentRegister_PrintSimilar( void* componentRegister, Name name, void* _stream, unsigned int number ) ;

	Stg_Component_DefaultConstructorFunction* Stg_ComponentRegister_AssertGet( 
		Stg_ComponentRegister* self,
		Name                   componentType,
		Name                   version );
	void Stg_ComponentRegister_PrintAllTypes( void* componentRegister, void* stream );

	/* Functions for iterating through the component element list ---------------------------------------------------*/
   int Stg_ComponentRegister_GetCount( void* componentRegister );
   Stg_ComponentRegisterElement* Stg_ComponentRegister_GetByIndex( void* componentRegister, int index );

	/** Obtain the component type from the component list element */
	Type Stg_ComponentRegisterElement_GetType( Stg_ComponentRegisterElement* element );

	/** Obtain the component version from the component list element */
	Name Stg_ComponentRegisterElement_GetVersion( Stg_ComponentRegisterElement* element );

	/** Obtain the component constructor function from the component list element */
	Stg_Component_DefaultConstructorFunction* Stg_ComponentRegisterElement_GetConstructor( Stg_ComponentRegisterElement* element );

	/** Obtain the component metadata from the component list element */
	Dictionary* Stg_ComponentRegisterElement_GetMetadata( Stg_ComponentRegisterElement* element );
	
#endif /* __Base_Automation_Stg_ComponentRegister_h__ */
