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
** $Id: Stg_ComponentRegister.c 2745 2005-05-1 08:12:18Z RaquibulHassan $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdarg.h>
#include "Base/Foundation/Foundation.h"
#include "Base/IO/IO.h"
#include "Base/Container/Container.h"

#include "types.h"
#include "shortcuts.h"
#include "Stg_Component.h"
#include "Stg_ComponentRegister.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>


/* Textual name of this class */
const Type Stg_ComponentRegister_Type = "Stg_ComponentRegister";
const Type Stg_ComponentRegisterElement_Type = "Stg_ComponentRegisterElement";

const Name Version = "0";
Stg_ComponentRegister *stgComponentRegister = NULL;

Stg_ComponentRegister *_Stg_ComponentRegister_New(  STG_COMPONENTREGISTER_DEFARGS  )
{
	Stg_ComponentRegister *self = NULL;

	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(Stg_ComponentRegister) );
	self = (Stg_ComponentRegister*)_Stg_Class_New(  STG_CLASS_PASSARGS  );

	return self;
}
	
Stg_ComponentRegister *Stg_ComponentRegister_New(  )
{
	/* Variables set in this function */
	SizeT                      _sizeOfSelf = sizeof( Stg_ComponentRegister );
	Type                              type = Stg_ComponentRegister_Type;
	Stg_Class_DeleteFunction*      _delete = _Stg_ComponentRegister_Delete;
	Stg_Class_PrintFunction*        _print = _Stg_ComponentRegister_Print;

	Stg_ComponentRegister *self = NULL;

	if( stgComponentRegister == NULL ){
		
	/* The following terms are parameters that have been passed into or defined in this function but are being set before being passed onto the parent */
	Stg_Class_CopyFunction*  _copy = NULL;

		self = _Stg_ComponentRegister_New(  STG_COMPONENTREGISTER_PASSARGS  );
		Stg_ComponentRegister_Init( self );
	}
	else{
		self = stgComponentRegister;
	}

	return self;
}
	
/* Initialisation implementation */
void _Stg_ComponentRegister_Init( Stg_ComponentRegister* self )
{
	assert( self );
	
	self->constructors = Stg_ObjectList_New();
   self->debugStream = Journal_Register( Debug_Type, "ComponentRegisterDebug" );
}
	
void Stg_ComponentRegister_Init( Stg_ComponentRegister* self )
{
	assert( self );
	_Stg_ComponentRegister_Init( self );
}
	
/* Delete boundary condition layout implementation */
void _Stg_ComponentRegister_Delete( void* componentRegister )
{
	Stg_ComponentRegister *self = NULL;

	self = (Stg_ComponentRegister*) componentRegister;
	assert( self );

   Stg_ObjectList_DeleteAllObjects( self->constructors );
	Stg_Class_Delete( self->constructors );
	_Stg_Class_Delete( self );
}
	
void _Stg_ComponentRegister_Print( void* componentRegister, Stream* stream )
{
	Stg_ComponentRegister *self = NULL;

	self = ( Stg_ComponentRegister* ) componentRegister;

	assert( self );
	
		/* General info */
	Journal_Printf( (void*) stream, "Stg_ComponentRegister (ptr): %p\n", self );
	_Stg_Class_Print( self, stream );
	Journal_Printf( stream, "Constructors:\n" );
	Stg_Class_Print( self->constructors, stream );
}
	
/* Public member functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
int Stg_ComponentRegister_AddFunc( 
		Stg_ComponentRegister *self,
		const Name componentType,
		const Name version,
		Stg_Component_DefaultConstructorFunction *func,
		Stg_Component_MetaAsDictionaryFunction* metadata )
{
	Stg_ComponentRegisterElement *element = NULL;

	assert( self );

   element = Stg_ComponentRegisterElement_New(Stg_ComponentRegisterElement_Type,
                                              _Stg_ComponentRegisterElement_Delete,
                                              _Stg_ComponentRegisterElement_Print,
                                              componentType,
                                              func,
                                              metadata,
                                              version
      );

   Journal_Printf( self->debugStream, "Adding [%s] to ComponentRegister\n", componentType);

   /* Copy component type into name field - used as key in objectList */
   element->name = StG_Strdup(componentType);   

   /* search object list to avoid duplicates... */
	Journal_Firewall( Stg_ObjectList_Get(self->constructors, element->name ) == NULL,
			Journal_Register( Error_Type, Stg_ComponentRegister_Type ), 
			"Error in func %s: Attempting to enter duplicate constructors for type '%s' (version '%s').\n"
			"This should only be done once per type.\n",
			__func__, componentType, version );

   /* Append element to list */
 	Stg_ObjectList_Append(self->constructors, element );
	return 1;
}

Bool Stg_ComponentRegister_RemoveEntry(
		Stg_ComponentRegister* self,
		Name                   componentType,
		Name                   version ) 
{
	assert( self );
   Stg_ObjectList_Remove( self->constructors, componentType, DELETE);
	return True;
}

Stg_Component_DefaultConstructorFunction* Stg_ComponentRegister_Get( 
		Stg_ComponentRegister* self,
		Name                   componentType,
		Name                   version ) 
{
	assert( self );
   /* Get the element object */
	Stg_ComponentRegisterElement *element = Stg_ObjectList_Get(self->constructors, componentType); 
	if ( element )
      /* Return the constructor function pointer */
	   return element->defaultConstructor;
   else
      return NULL;
}

Stg_Component_DefaultConstructorFunction* Stg_ComponentRegister_AssertGet( 
		Stg_ComponentRegister* self,
		Name                   componentType,
		Name                   version ) 
{
	Stg_Component_DefaultConstructorFunction* componentConstructorFunction;
	
	componentConstructorFunction = Stg_ComponentRegister_Get( self, componentType, version );

	/* If we cannot find the default construct for this componentType - then abort() with a nice message */
	if ( !componentConstructorFunction ) {
		Stream* errorStream = Journal_Register( Error_Type, self->type );

		Journal_Printf( errorStream, "Cannot find default constructor function for type '%s'\n", componentType );
		Journal_Printf( errorStream, "Could you have meant one of these?\n" );

		Stream_Indent( errorStream );
		Stg_ComponentRegister_PrintSimilar( self, componentType, errorStream, 5 );
		abort();
	}	

	return componentConstructorFunction;
}

Dictionary* Stg_ComponentRegister_GetMetadata(
		Stg_ComponentRegister* self,
		Name                   componentType,
		Name                   version ) 
{
	Stg_ComponentRegisterElement *element = NULL;
	assert( self );
	element = Stg_ObjectList_Get(self->constructors, componentType); 
	if( element ){
		return (void*) element->metadata();
	}

	return NULL;
}

Stg_ComponentRegister *Stg_ComponentRegister_Get_ComponentRegister() {
	return stgComponentRegister;
}

int Stg_ComponentRegister_GetCount( void* componentRegister ) {
	Stg_ComponentRegister* self = (Stg_ComponentRegister*)componentRegister;
	return self->constructors->count;
}

/** Stg_ComponentRegisterElement methods: 
 * Constructor interface... */
Stg_ComponentRegisterElement* Stg_ComponentRegisterElement_New(
		Type			type,
		Stg_Class_DeleteFunction*	_delete,
		Stg_Class_PrintFunction*	_print,
      Type        componentType,
		Stg_Component_DefaultConstructorFunction*		defaultConstructor,
		Stg_Component_MetaAsDictionaryFunction*      metadata,
		Name								version
      )
{
   Stg_ComponentRegisterElement* self = ( Stg_ComponentRegisterElement*)_Stg_Class_New( 
                                          sizeof(Stg_ComponentRegisterElement),
                                          type, _delete, _print, _Stg_Class_Copy);

   self->componentType = StG_Strdup(componentType);
   self->defaultConstructor = defaultConstructor;
   self->metadata = metadata;
   self->version = StG_Strdup(version);
   return self;   /* How in the hell this worked previously without returning anything I have no idea */ 
}

void _Stg_ComponentRegisterElement_Delete( void* self )
{
   Stg_ComponentRegisterElement* element = (Stg_ComponentRegisterElement*)self; 
	if( element ){
		if (element->componentType) Memory_Free( element->componentType );
		if (element->version) Memory_Free( element->version );
		if (element->name) Memory_Free( element->name);
      _Stg_Class_Delete( element ); /* element's parent is a class so delete it */
	}
}
Stg_ComponentRegisterElement* Stg_ComponentRegister_GetByIndex( void* componentRegister, int index ) {
	Stg_ComponentRegister* self = (Stg_ComponentRegister*)componentRegister;
   assert(index < self->constructors->count);
	return (Stg_ComponentRegisterElement*)self->constructors->data[index]; 
}

Type Stg_ComponentRegisterElement_GetType( Stg_ComponentRegisterElement* element ) {
	return element->componentType;
}

Name Stg_ComponentRegisterElement_GetVersion( Stg_ComponentRegisterElement* element ) {
	return element->version;
}

Stg_Component_DefaultConstructorFunction* Stg_ComponentRegisterElement_GetConstructor( Stg_ComponentRegisterElement* element ) {
	return element->defaultConstructor;
}

Dictionary* Stg_ComponentRegisterElement_GetMetadata( Stg_ComponentRegisterElement* element ) {
	return element->metadata();
}

void _Stg_ComponentRegisterElement_Print( void* self, Stream* printStream )	
{
   Stg_ComponentRegisterElement* element = (Stg_ComponentRegisterElement*)self; 
	Journal_Printf( printStream, "Constructor Information\n");
	Journal_Printf( printStream, "\tStg_ComponentType                : %s\n", element->componentType );
	Journal_Printf( printStream, "\tStg_Component Default Constructor: %p\n", element->defaultConstructor );
	Journal_Printf( printStream, "\tVersion                      : %s\n", element->version );
}



void Stg_ComponentRegister_PrintSimilar( void* componentRegister, Name name, void* _stream, unsigned int number ) {
	Stg_ComponentRegister*                  self               = (Stg_ComponentRegister*) componentRegister;
   Stg_ObjectList_PrintSimilar( self->constructors, name, _stream, number );
}

void Stg_ComponentRegister_PrintAllTypes( void* componentRegister, void* stream ) {
	Stg_ComponentRegister*                  self               = (Stg_ComponentRegister*) componentRegister;

	/* Parse the list, printing all the names */
   int i;
	Stg_ComponentRegisterElement* element;
	for(i = 0; i < Stg_ComponentRegister_GetCount(self); i++) {
      element = Stg_ComponentRegister_GetByIndex(self, i);
	   Journal_Printf( stream, "%s\n", element->componentType );
	}
}


