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
** $Id: SimpleExtensionInfo.c 2915 2005-05-08 12:51:46Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdarg.h>
#include "Base/Foundation/Foundation.h"
#include "Base/IO/IO.h"
#include "Base/Container/Container.h"
#include "Base/Automation/Automation.h"

#include "types.h"
#include "shortcuts.h"
#include "ExtensionManager.h"
#include "ExtensionInfo.h"
#include "SimpleExtensionInfo.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>


/* Class Administration members ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


/* Textual name of this class */
const Type SimpleExtensionInfo_Type = "SimpleExtensionInfo";


/** allocate and initialise a new SimpleExtensionInfo. */
SimpleExtensionInfo* SimpleExtensionInfo_New( Name name, SizeT size, Index count ) {
	/* Variables set in this function */
	SizeT                            _sizeOfSelf = sizeof(SimpleExtensionInfo);
	Type                                    type = SimpleExtensionInfo_Type;
	Stg_Class_DeleteFunction*            _delete = _SimpleExtensionInfo_Delete;
	Stg_Class_PrintFunction*              _print = _SimpleExtensionInfo_Print;
	Stg_Class_CopyFunction*                _copy = _SimpleExtensionInfo_Copy;
	ExtensionInfo_DataCopyFunction*    _dataCopy = _SimpleExtensionInfo_DataCopy;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

	return _SimpleExtensionInfo_New(  SIMPLEEXTENSIONINFO_PASSARGS  );
}


void SimpleExtensionInfo_Init( void* simpleExtensionInfo, Name name, SizeT size, Index count ) {
	SimpleExtensionInfo* self = (SimpleExtensionInfo*)simpleExtensionInfo;

	/* General info */
	self->type = SimpleExtensionInfo_Type;
	self->_sizeOfSelf = sizeof(SimpleExtensionInfo);
	self->_deleteSelf = False;
	
	/* Virtual info */
	self->_delete = _SimpleExtensionInfo_Delete;
	self->_print = _SimpleExtensionInfo_Print;
	self->_copy = _SimpleExtensionInfo_Copy;
	self->_dataCopy = _SimpleExtensionInfo_DataCopy;
	
	/* SimpleExtensionInfo info */
	_SimpleExtensionInfo_Init( self, (Name)name, size, count );
}


SimpleExtensionInfo* _SimpleExtensionInfo_New(  SIMPLEEXTENSIONINFO_DEFARGS  )
{
	SimpleExtensionInfo* self;
	
	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(SimpleExtensionInfo) );
	self = (SimpleExtensionInfo*)_ExtensionInfo_New(  EXTENSIONINFO_PASSARGS  );
	
	/* General info */
	
	/* Virtual info */
	
	/* SimpleExtensionInfo info */
	_SimpleExtensionInfo_Init( self, name, size, count );
	
	return self;
}


void _SimpleExtensionInfo_Init( SimpleExtensionInfo* self, Name name, SizeT size, Index count ) {
	/* General and Virtual info should already be set */
	/* SimpleExtensionInfo info */
	_ExtensionInfo_Init( (ExtensionInfo*)self, name, size, count );
}


void _SimpleExtensionInfo_Delete( void* simpleExtensionInfo ) {
	SimpleExtensionInfo* self = (SimpleExtensionInfo*)simpleExtensionInfo;

	Journal_DFirewall( 
		self != NULL,  
		Journal_Register( ErrorStream_Type, ExtensionInfo_Type ), 
		"Attempting to \"%s\" delete a pointer of value NULL\n", SimpleExtensionInfo_Type );
	
	/* Purposely do nothing... the extension is assumed to be static (non-pointer) data */
	
	/* Delete parent */
	_ExtensionInfo_Delete( self );
}


void _SimpleExtensionInfo_Print( void* simpleExtensionInfo, Stream* stream ) {
	SimpleExtensionInfo* self = (SimpleExtensionInfo*)simpleExtensionInfo;
		
	Journal_Firewall( 
		self != NULL,  
		Journal_Register( ErrorStream_Type, ExtensionInfo_Type ), 
		"Attempting to \"%s\" print a pointer of value NULL\n", SimpleExtensionInfo_Type );
	
	/* General info */
	Journal_Printf( (void*)stream, "SimpleExtensionInfo (%p):\n", self );
	_ExtensionInfo_Print( self, stream );
	Journal_Printf( (void*)stream, "\tNo print function implemented/assigned for this extension\n" );
}


void* _SimpleExtensionInfo_Copy( const void* simpleExtensionInfo, void* dest, Bool deep, Name nameExt, struct PtrMap* ptrMap  ) {
	SimpleExtensionInfo*	self = (SimpleExtensionInfo*)simpleExtensionInfo;
	SimpleExtensionInfo*	newSimpleExtensionInfo;
		
	Journal_Firewall( 
		self != NULL,  
		Journal_Register( ErrorStream_Type, SimpleExtensionInfo_Type ), 
		"Attempting to \"%s\" copy a pointer of value NULL\n", SimpleExtensionInfo_Type );
	
	/* Copy parent */
	newSimpleExtensionInfo = (SimpleExtensionInfo*)_ExtensionInfo_Copy( self, dest, deep, nameExt, ptrMap );
	
	return newSimpleExtensionInfo;
}

void* _SimpleExtensionInfo_DataCopy( 
	void* simpleExtensionInfo, 
	void* source,
	void* dest, 
	Bool deep, 
	Name nameExt, 
	PtrMap* ptrMap )
{
	SimpleExtensionInfo* self = (SimpleExtensionInfo*)simpleExtensionInfo;
	
	memcpy( dest, source, self->size );

	return dest;
}

/* Public member functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/* Private member functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


