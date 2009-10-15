/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003-2006, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street,
**	Melbourne, 3053, Australia.
**
** Primary Contributing Organisations:
**	Victorian Partnership for Advanced Computing Ltd, Computational Software Development - http://csd.vpac.org
**	Australian Computational Earth Systems Simulator - http://www.access.edu.au
**	Monash Cluster Computing - http://www.mcc.monash.edu.au
**	Computational Infrastructure for Geodynamics - http://www.geodynamics.org
**
** Contributors:
**	Patrick D. Sunter, Software Engineer, VPAC. (pds@vpac.org)
**	Robert Turnbull, Research Assistant, Monash University. (robert.turnbull@sci.monash.edu.au)
**	Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**	David May, PhD Student, Monash University (david.may@sci.monash.edu.au)
**	Louis Moresi, Associate Professor, Monash University. (louis.moresi@sci.monash.edu.au)
**	Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**	Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
**	Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
**	Julian Giordani, Research Assistant, Monash University. (julian.giordani@sci.monash.edu.au)
**	Vincent Lemiale, Postdoctoral Fellow, Monash University. (vincent.lemiale@sci.monash.edu.au)
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
** $Id: ElementType_Register.c 964 2007-10-11 08:03:06Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include "units.h"
#include "types.h"
#include "shortcuts.h"
#include "ElementType.h"
#include "ElementType_Register.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

/* Textual name of this class */
const Type ElementType_Register_Type = "ElementType_Register";

ElementType_Register* elementType_Register = 0;

void* ElementType_Register_DefaultNew( Name name ) {
	return (void*) _ElementType_Register_New( 
			sizeof(ElementType_Register), 
			ElementType_Register_Type,
			_ElementType_Register_Delete, 
			_ElementType_Register_Print, 
			NULL,
			ElementType_Register_DefaultNew,
			_ElementType_Register_Construct,
			_ElementType_Register_Build,
			_ElementType_Register_Initialise,
			_ElementType_Register_Execute,
			_ElementType_Register_Destroy,
			name,
			False);
}

ElementType_Register* ElementType_Register_New( Name name ) {
	return _ElementType_Register_New( sizeof(ElementType_Register), ElementType_Register_Type,
			_ElementType_Register_Delete, _ElementType_Register_Print, NULL,
			ElementType_Register_DefaultNew,
			_ElementType_Register_Construct,
			_ElementType_Register_Build,
			_ElementType_Register_Initialise,
			_ElementType_Register_Execute,
			_ElementType_Register_Destroy,
			name,
			True
			);
}

void ElementType_Register_Init( void* elementType_Register, Name name ) {
	ElementType_Register* self = (ElementType_Register*)elementType_Register;
	
	/* General info */
	self->type = ElementType_Register_Type;
	self->_sizeOfSelf = sizeof(ElementType_Register);
	self->_deleteSelf = False;
	
	/* Virtual info */
	self->_delete = _ElementType_Register_Delete;
	self->_print = _ElementType_Register_Print;
	self->_copy = NULL;
	self->_defaultConstructor = ElementType_Register_DefaultNew;
	self->_construct = _ElementType_Register_Construct;
	self->_build = _ElementType_Register_Build;
	self->_initialise = _ElementType_Register_Initialise;
	self->_execute = _ElementType_Register_Execute;
	self->_destroy = _ElementType_Register_Destroy;

	_Stg_Class_Init( (Stg_Class*)self );
	_Stg_Object_Init( (Stg_Object*)self, name, NON_GLOBAL );
	_Stg_Component_Init( (Stg_Component*)self );
	
	/* ElementType_Register info */
	_ElementType_Register_Init( self );
}

ElementType_Register* _ElementType_Register_New(
		SizeT				_sizeOfSelf,
		Type				type,
		Stg_Class_DeleteFunction*		_delete,
		Stg_Class_PrintFunction*		_print,
		Stg_Class_CopyFunction*		_copy,
		Stg_Component_DefaultConstructorFunction*	_defaultConstructor,
		Stg_Component_ConstructFunction*			_construct,
		Stg_Component_BuildFunction*		_build,
		Stg_Component_InitialiseFunction*		_initialise,
		Stg_Component_ExecuteFunction*		_execute,
		Stg_Component_DestroyFunction*		_destroy,
		Name							name,
		Bool							initFlag )
{
	ElementType_Register* self;
	
	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(ElementType_Register) );
	self = (ElementType_Register*)_Stg_Component_New( _sizeOfSelf, type, _delete, _print, _copy, _defaultConstructor, _construct, _build, 
			_initialise, _execute, _destroy, name, NON_GLOBAL );
	
	/* General info */
	
	/* Virtual info */
	
	if( initFlag ){
		_ElementType_Register_Init( self );
	}
	
	return self;
}

void _ElementType_Register_Init( void* elementType_Register ) {
	ElementType_Register* self = (ElementType_Register*)elementType_Register;
	
	/* General and Virtual info should already be set */
	
	/* ElementType_Register info */
	self->isConstructed = True;
	self->count = 0;
	self->_size = 8;
	self->_delta = 8;
	self->elementType = Memory_Alloc_Array( ElementType*, self->_size, "ElementType_Register->elementType" );
	memset( self->elementType, 0, sizeof(ElementType*) * self->_size );
	self->debug = Stream_RegisterChild( StgFEM_Discretisation_Debug, ElementType_Register_Type );
}

void _ElementType_Register_Delete( void* elementType_Register ) {
	ElementType_Register* self = (ElementType_Register*)elementType_Register;

	Journal_DPrintf( self->debug, "In %s\n", __func__ );
	Stream_IndentBranch( StgFEM_Debug );

	/* Assumes ownerships of the element types */
	if( self->elementType ) {
		ElementType_Index elementType_I;
		
		for( elementType_I = 0; elementType_I < self->count; elementType_I++ ) {
			_Stg_Object_Delete( self->elementType[elementType_I] );
		}
		
		Memory_Free( self->elementType );
	}
	
	/* Stg_Class_Delete parent */
	_Stg_Class_Delete( self );
	Stream_UnIndentBranch( StgFEM_Debug );
}

void _ElementType_Register_Print( void* elementType_Register, Stream* stream ) {
	ElementType_Register* self = (ElementType_Register*)elementType_Register;
	ElementType_Index elementType_I;
	
	/* Set the Journal for printing informations */
	Stream* elementType_RegisterStream = stream;
	
	/* General info */
	Journal_Printf( stream, "ElementType_Register (ptr): %p\n", self );
	
	/* Print parent */
	_Stg_Class_Print( self, elementType_RegisterStream );
	
	/* Virtual info */
	
	/* ElementType_Register info */
	Journal_Printf( stream, "\tcount: %u\n", self->count );
	Journal_Printf( stream, "\t_size: %lu\n", self->_size );
	Journal_Printf( stream, "\t_delta: %lu\n", self->_delta );
	
	Journal_Printf( stream, "\telementType (ptr): %p\n", self->elementType );
	Journal_Printf( stream, "\telementType[0-%u]:\n", self->count );
	for( elementType_I = 0; elementType_I < self->count; elementType_I++ ) {
		Journal_Printf( stream, "elementType[%u]: ", elementType_I );
		Stg_Class_Print( self->elementType[elementType_I], elementType_RegisterStream );
	}
	Journal_Printf( stream, "\t]\n" );
}

void _ElementType_Register_Construct( void* elementType_Register, Stg_ComponentFactory *cf, void* data ){
	ElementType_Register*	self = (ElementType_Register*)elementType_Register;
	
	self->context = Stg_ComponentFactory_ConstructByKey( cf, self->name, "Context", DomainContext, False, data );
	if( !self->context ) 
		self->context = Stg_ComponentFactory_ConstructByName( cf, "context", DomainContext, True, data );
}
	
void _ElementType_Register_Build( void* elementType_Register, void *data ){
	
}
	
void _ElementType_Register_Initialise( void* elementType_Register, void *data ){
	
}
	
void _ElementType_Register_Execute( void* elementType_Register, void *data ){
	
}

void _ElementType_Register_Destroy( void* elementType_Register, void *data ){
	
}

ElementType_Index ElementType_Register_Add( void* elementType_Register, void* elementType ) {
	ElementType_Register*	self = (ElementType_Register*)elementType_Register;
	ElementType_Index	handle;
	
	if( self->count >= self->_size ) {
		ElementType**	newElementType;
		
		self->_size += self->_delta;
		newElementType = Memory_Alloc_Array( ElementType*, self->_size, "ElementType_Register->elementType" );
		memcpy( newElementType, self->elementType, sizeof(ElementType*) * self->count );
		Memory_Free( self->elementType );
		self->elementType = newElementType;
	}
	
	handle = self->count;
	self->elementType[handle] = (ElementType*)elementType;
	self->count++;
	
	/* Build the elementType... i.e assume it hasn't been built already */
	ElementType_Build( self->elementType[handle], NULL );
	
	return handle;
}

ElementType_Index ElementType_Register_GetIndex( void* elementType_Register, Type type ) {
	ElementType_Register*	self = (ElementType_Register*)elementType_Register;
	ElementType_Index elementType_I;
	
	for( elementType_I = 0; elementType_I < self->count; elementType_I++ ) {
		if( self->elementType[elementType_I]->type == type ) {
			return elementType_I;
		}
	}
	return (unsigned)-1;
}

ElementType* _ElementType_Register_At( void* elementType_Register, ElementType_Index handle ) {
	ElementType_Register*	self = (ElementType_Register*)elementType_Register;
	
	return ElementType_Register_At( self, handle );
}
