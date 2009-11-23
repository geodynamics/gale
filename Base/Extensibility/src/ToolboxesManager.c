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
** $Id: ToolboxesManager.c 4081 2007-04-27 06:20:07Z LukeHodkinson $
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
#include "ModulesManager.h"
#include "ToolboxesManager.h"

#include <stdlib.h>
#include <string.h>

/* Textual name of this class */
const Type ToolboxesManager_Type = "ToolboxesManager";


ToolboxesManager* ToolboxesManager_New( int* argc, char*** argv ) {
	return _ToolboxesManager_New( 
		sizeof(ToolboxesManager), 
		ToolboxesManager_Type, 
		_ToolboxesManager_Delete, 
		_ToolboxesManager_Print, 
		NULL, 
		_ToolboxesManager_GetToolboxesList,
		_ToolboxesManager_LoadToolbox,
		_ToolboxesManager_UnloadToolbox,
		Toolbox_Factory,
		argc,
		argv );
}

ToolboxesManager* _ToolboxesManager_New(
		SizeT                                   _sizeOfSelf,
		Type                                    type,
		Stg_Class_DeleteFunction*               _delete,
		Stg_Class_PrintFunction*                _print,
		Stg_Class_CopyFunction*                 _copy, 
		ModulesManager_GetModulesListFunction*  _getModulesList,
		ModulesManager_LoadModuleFunction*	_loadModule,
		ModulesManager_UnloadModuleFunction*	_unloadModule,
		ModulesManager_ModuleFactoryFunction*   _moduleFactory,
		int*					argc,
		char***					argv )
{
	ToolboxesManager* self;
	
	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(ToolboxesManager) );
	self = (ToolboxesManager*)_ModulesManager_New( 
		_sizeOfSelf, 
		type, 
		_delete, 
		_print, 
		_copy, 
		_getModulesList, 
		_loadModule, 
		_unloadModule,
		_moduleFactory );
	
	/* General info */
	
	/* Virtual info */
	
	_ToolboxesManager_Init( self, argc, argv );
	
	return self;
}

void _ToolboxesManager_Init( void* toolboxesManager, int* argc, char*** argv ) {
	ToolboxesManager*         self = (ToolboxesManager*)toolboxesManager;
	
	self->argc = argc;
	self->argv = argv;
	self->_initialisedSize = 8;
	self->initialised = Memory_Alloc_Array( char*, self->_initialisedSize, ToolboxesManager_Type );
	self->_initialisedCount = 0;
}


void _ToolboxesManager_Delete( void* toolboxesManager ) {
	ToolboxesManager*         self = (ToolboxesManager*)toolboxesManager;

	Memory_Free( self->initialised );
	self->_initialisedSize = 0;
	self->_initialisedCount = 0;
	self->initialised = 0;
	
	/* Delete parent */
	_ModulesManager_Delete( self );
}

void _ToolboxesManager_Print( void* toolboxesManager, Stream* stream ) {
	ToolboxesManager* self = (ToolboxesManager*)toolboxesManager;
	
	/* General info */
	Journal_Printf( (void*) stream, "Toolboxes (ptr): %p\n", self );
	
	if( self->_initialisedCount > 0 ) {
		Index i;
		
		Journal_Printf( stream, "Initialised Modules:\n" );
		Stream_Indent( stream );
		for( i = 0; i < self->_initialisedCount; ++i ) {
			Journal_Printf( stream, "%s\n", self->initialised[i] );
		}
		Stream_UnIndent( stream );
	}
	
	/* Print parent */
	_ModulesManager_Print( self, stream );
}


Dictionary_Entry_Value* _ToolboxesManager_GetToolboxesList( void* toolboxesManager, void* _dictionary ) {
	/*ToolboxesManager*		self = (ToolboxesManager*)toolboxesManager;*/
	Dictionary*			dictionary = (Dictionary*)_dictionary;
	Dictionary_Entry_Value*		pluginsList = NULL;
	
	pluginsList = Dictionary_Get( dictionary, "import" );

	return pluginsList;
}

Bool _ToolboxesManager_LoadToolbox( void* toolboxesManager, Module* toolbox ) {
	ToolboxesManager* self = (ToolboxesManager*)toolboxesManager;
	
	((Toolbox*)toolbox)->Initialise( self, self->argc, self->argv );
	((Toolbox*)toolbox)->Register( self );
    
	return True;
}

Bool _ToolboxesManager_UnloadToolbox( void* toolboxesManager, Module* toolbox ) {
	ToolboxesManager* self = (ToolboxesManager*)toolboxesManager;
	
	((Toolbox*)toolbox)->Finalise( self );
    
	return True;
}

Index ToolboxesManager_SetInitialised( void* initRegister, char* label ) {
	ToolboxesManager* self = (ToolboxesManager*)initRegister;
	
	if( self->_initialisedCount == self->_initialisedSize ) {
		Index oldSize = self->_initialisedSize;
		char** tmp;
		
		self->_initialisedSize += 8;
		tmp = Memory_Alloc_Array( char*, self->_initialisedSize, ToolboxesManager_Type );
		memcpy( tmp, self->initialised, sizeof( char* ) * oldSize );
		Memory_Free( self->initialised );
		self->initialised = tmp;
	}
	
	self->initialised[self->_initialisedCount] = StG_Strdup( label );
	self->_initialisedCount += 1;
	return self->_initialisedCount - 1;
}

Bool ToolboxesManager_IsInitialised( void* initRegister, char* label ) {
	ToolboxesManager* self = (ToolboxesManager*)initRegister;
	Index i;
	
	for( i = 0; i < self->_initialisedCount; i++ ) {
		if( strcmp( label, self->initialised[i] ) == 0 ) {
			return True;
		}
	}
	return False;
}
