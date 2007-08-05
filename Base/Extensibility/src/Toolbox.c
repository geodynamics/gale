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
** $Id: Toolbox.c 3192 2005-08-25 01:45:42Z AlanLo $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdarg.h>
#include "Base/Foundation/Foundation.h"
#include "Base/IO/IO.h"
#include "Base/Container/Container.h"
#include "Base/Automation/Automation.h"

#include "types.h"
#include "shortcuts.h"
#include "Module.h"
#include "Toolbox.h"

#include <string.h>

const Type Toolbox_Type = "Toolbox";

static const char* TOOLBOX_REGISTER_SUFFIX = "_Register";
static const char* TOOLBOX_INITIALISE_SUFFIX = "_Initialise";
static const char* TOOLBOX_FINALISE_SUFFIX = "_Finalise";
static const char* TOOLBOX_MODULE_SUFFIX = "_Toolbox";
#ifdef MEMORY_STATS
	static const char* TOOLBOX_MANGLEDNAME = "mangledName";
#endif

Toolbox* Toolbox_New( Name name, Stg_ObjectList* directories ) {
	return _Toolbox_New( 
		sizeof(Toolbox), 
		Toolbox_Type, 
		_Toolbox_Delete, 
		_Toolbox_Print, 
		NULL,
		name,
		_Toolbox_MangleName,
		directories );
}

Module* Toolbox_Factory( Name name, Stg_ObjectList* directories ) {
	return (Module*)Toolbox_New( name, directories );
}
	
Toolbox* _Toolbox_New( 
		SizeT                        _sizeOfSelf,
		Type                         type,
		Stg_Class_DeleteFunction*    _delete,
		Stg_Class_PrintFunction*     _print,
		Stg_Class_CopyFunction*      _copy, 
		Name                         name,
		Module_MangleNameFunction    MangleName,
		Stg_ObjectList*              directories )
{
	Toolbox* self;

	assert( _sizeOfSelf >= sizeof(Toolbox) );

	self = (Toolbox*)_Module_New( _sizeOfSelf, type, _delete, _print, _copy, name, MangleName, directories );
	
	_Toolbox_Init( self );

	return self;
}
	
void _Toolbox_Init( Toolbox* self ) {
	Stream* stream = Journal_Register( Info_Type, self->type );

	if( self->dllPtr != NULL ) {
		self->Register = (Toolbox_RegisterFunction*)Module_LoadSymbol( self, TOOLBOX_REGISTER_SUFFIX );
		self->Initialise = (Toolbox_InitialiseFunction*)Module_LoadSymbol( self, TOOLBOX_INITIALISE_SUFFIX );
		self->Finalise = (Toolbox_FinaliseFunction*)Module_LoadSymbol( self, TOOLBOX_FINALISE_SUFFIX );
	}
	/* If the register function is not found, then unload the module... it's not a toolbox. */
	if( self->Register == NULL || self->Initialise == NULL || self->Finalise == NULL ) {
		Journal_Printf( stream, "Toolbox %s is not a toolbox, unloading.\n", self->name );
		Module_UnLoad( self );
	}
}
	
void _Toolbox_Delete( void* toolbox ) {
	Toolbox* self = (Toolbox*)toolbox;

	/* Delete parent */
	_Module_Delete( self );
}
	
void _Toolbox_Print( void* toolbox, Stream* stream ) {
	Toolbox* self = (Toolbox*)toolbox;

	Journal_Printf( stream, "Toolbox: %s\n", self->name );
	Stream_Indent( stream );
	
	/* Print parent */
	_Module_Print( self, stream );
	
	Stream_UnIndent( stream );
}

char* _Toolbox_MangleName( char* name ) {
	char* mangledName = Memory_Alloc_Array( char, strlen( name ) + strlen( TOOLBOX_MODULE_SUFFIX ) + 1, TOOLBOX_MANGLEDNAME );
	sprintf( mangledName, "%s%s", name, TOOLBOX_MODULE_SUFFIX );
	return mangledName;
}


Toolbox_RegisterFunction* Toolbox_GetRegisterFunc( void* toolbox ) {
	Toolbox* self = (Toolbox*)toolbox;

	return self->Register;
}

Toolbox_InitialiseFunction* Toolbox_GetInitialiseFunc( void* toolbox ) {
	Toolbox* self = (Toolbox*)toolbox;

	return self->Initialise;
}

Toolbox_FinaliseFunction* Toolbox_GetFinaliseFunc( void* toolbox ) {
	Toolbox* self = (Toolbox*)toolbox;

	return self->Finalise;
}
