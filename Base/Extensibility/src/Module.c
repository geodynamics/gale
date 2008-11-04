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
** $Id: Module.c 3192 2005-08-25 01:45:42Z AlanLo $
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


#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <errno.h>

#ifdef NOSHARED
	#define dlopen( x, y ) (NULL)
	#define dlsym( x, y ) (NULL)
	#define dlclose( x )
	#define dlerror() "blah"
	#define RTLD_LAZY 0
	#define RTLD_GLOBAL 0
#else
        #include <dlfcn.h>
#endif

#include <limits.h>


const Type Module_Type = "Module";

static const char* MODULE_SUFFIX = "module.";
static const char* MODULE_GETMETADATA_SUFFIX = "_MetaAsDictionary";
static const char* MODULE_GETNAME_SUFFIX = "_GetName";
static const char* MODULE_GETVERSION_SUFFIX = "_GetVersion";

const char* PLUGIN_DEPENDENCY_KEY = "associations";
const char* PLUGIN_DEPENDENCY_NAME_KEY = "type";
const char* PLUGIN_DEPENDENCY_VERSION_KEY = "version";
const char* PLUGIN_DEPENDENCY_URL_KEY = "url";

#ifdef MEMORY_STATS
	static const char* MODULE_FILENAME = "fileName";
	static const char* MODULE_SYMBOLNAME = "symbolName";
#endif


Module* _Module_New( 
		SizeT                        _sizeOfSelf,
		Type                         type,
		Stg_Class_DeleteFunction*    _delete,
		Stg_Class_PrintFunction*     _print,
		Stg_Class_CopyFunction*      _copy, 
		Name                         name,
		Module_MangleNameFunction    MangleName,
		Stg_ObjectList*              directories )
{
	Module* self;

	assert( _sizeOfSelf >= sizeof(Module) );

	self = (Module*)_Stg_Object_New( _sizeOfSelf, type, _delete, _print, _copy, name, NON_GLOBAL );
	
	_Module_Init( self, MangleName, directories );

	return self;
}
	
void _Module_Init(
		Module*                      self,
		Module_MangleNameFunction    MangleName,
		Stg_ObjectList*              directories )
{
	char*                           mangledName;
	char*                           fileName = NULL;
	char*                           fullPathName = NULL;
	int                             fullPathLength = 0;
	int                             length;
	Index                           dir_i;

	Stream*                         stream;
	Stream*                         debug;
	Stream*                         error;

	XML_IO_Handler*	                ioHandler;
	const char*	                metadata;
	
	stream =  Journal_Register( Info_Type, self->type );
	debug =  Journal_Register( Debug_Type, self->type );
	error =  Journal_Register( Error_Type, self->type );
	
	self->MangleName = MangleName;
	
	
	Journal_Printf( debug, "Finding module: \"%s\"... ", self->name );

	/* Try the plugin name by itself (allows LD_LIBRARY_PATH) to take precendence */
	mangledName = self->MangleName( self->name );
	fileName = Memory_Alloc_Array( char, strlen(mangledName) + strlen(MODULE_SUFFIX) + strlen(MODULE_EXT) +1, MODULE_FILENAME );
	sprintf( fileName, "%s%s%s", mangledName, MODULE_SUFFIX, MODULE_EXT );
		
	self->dllPtr = dlopen( fileName, RTLD_LAZY | RTLD_GLOBAL );
	if( self->dllPtr ) {
		Journal_RPrintf( stream, "%s \"%s\" found using %s\n", self->type, self->name, fileName );
	}
	else {
		Journal_RPrintf( stream, "%s \"%s\" failed: %s\n", self->type, self->name, dlerror() );
		
		for ( dir_i = 0; dir_i < directories->count; ++dir_i ) {
			length = strlen(Stg_ObjectList_ObjectAt( directories, dir_i )) + 1 + strlen(fileName) + 1;
			if ( fullPathLength < length ) {
				fullPathLength = length;
				fullPathName = Memory_Realloc_Array( fullPathName, char, fullPathLength );
			}
			PathJoin( fullPathName, 2, Stg_ObjectList_ObjectAt( directories, dir_i ), fileName );
			self->dllPtr = dlopen( fullPathName, RTLD_LAZY | RTLD_GLOBAL );
			if( self->dllPtr ) {
				Journal_RPrintf( stream, "%s \"%s\" found using %s\n", self->type, self->name, fileName );
				break;
			}
			else {
				Journal_RPrintf( stream, "%s \"%s\" failed: %s\n", self->type, self->name, dlerror() );
			}
		}
		/* If it failed alltogether, print a error message. */
		if ( dir_i == directories->count) {
			Journal_Printf( debug, "failed to find in any of directories above, or had error.\n" );
		}
	}

	/* Load the symbols */
	if( self->dllPtr ) {
		self->GetMetadata = (Module_GetMetadataFunction*)Module_LoadSymbol( self, MODULE_GETMETADATA_SUFFIX );
		self->GetName = (Module_GetNameFunction*)Module_LoadSymbol( self, MODULE_GETNAME_SUFFIX );
		self->GetVersion = (Module_GetVersionFunction*)Module_LoadSymbol( self, MODULE_GETVERSION_SUFFIX );
	}
	
	/* Load the meta data */
	if( self->GetMetadata ) {
		self->_meta = self->GetMetadata();
		/* TODO: Needs to be replaced with the meta xsd equivalent SQ-20081103
		if( metadata ) {
			ioHandler = XML_IO_Handler_New();
			IO_Handler_ReadAllFromBuffer( ioHandler, metadata, self->_meta );
			Stg_Class_Delete( ioHandler );
		}*/
	}
	else {
		self->_meta = Dictionary_New();
	}

	if ( fullPathName ) {
		Memory_Free( fullPathName );
	}
	Memory_Free( fileName );
	Memory_Free( mangledName );
}
	
void _Module_Delete( void* module ) {
	Module* self = (Module*)module;

	Module_UnLoad( self );
	Stg_Class_Delete( self->_meta );
	
	/* Delete parent */
	_Stg_Object_Delete( self );
}
	
void _Module_Print( void* module, Stream* stream ) {
	Module* self = (Module*)module;

	Dictionary* depList;
	Index count = 0;
	Index i;
	const char* version;

	Journal_RPrintf( stream, "Module: %s\n", self->name );
	Stream_Indent( stream );
	
	version = Module_GetVersion( self );
	if ( version == NULL ) {
		version = "Unknown";
	}
	Journal_Printf( stream, "Version: (version) %s\n", version );

	depList = Module_GetDependencies( self );
	if ( depList ) {
		count = Dictionary_GetCount( depList );
	}
	for ( i = 0; i < count; ++i ) {
		Dictionary* dep;
		char* depName;
		char* depVersion;
		char* depUrl;
		
		dep = Dictionary_Entry_Value_AsDictionary( Dictionary_GetByIndex( depList, i ) );
		depName = Dictionary_GetString_WithDefault( dep, (char*)PLUGIN_DEPENDENCY_NAME_KEY, "NULL" );
		depVersion = Dictionary_GetString_WithDefault( dep, (char*)PLUGIN_DEPENDENCY_VERSION_KEY, "Unknown" );
		depUrl = Dictionary_GetString_WithDefault( dep, (char*)PLUGIN_DEPENDENCY_URL_KEY, "None" );
	
		Journal_RPrintf( stream, "Depends on: %s, version=%s, url=%s\n", depName, depVersion, depUrl );
	}
	
	/* Print parent */
	_Stg_Object_Print( self, stream );
	
	Stream_UnIndent( stream );
}


Dictionary* Module_GetMetadata( void* module ) {
	Module* self = (Module*)module;
	
	return self->_meta;
}

const char* Module_GetName( void* module ) {
	Module* self = (Module*)module;

	if ( self->GetName ) {
		return self->GetName();
	}
	return self->name;
}

const char* Module_GetVersion( void* module ) {
	Module* self = (Module*)module;

	if ( self->GetVersion ) {
		return self->GetVersion();
	}
	return NULL;
}


Dictionary* Module_GetDependencies( void* module ) {
	Module* self = (Module*)module;
	Dictionary_Entry_Value* dev;

	/* Kinda important... we want to only return a real dictionary if the meta data has a dependency entry, even if it is
	   empty. Otherwise, indicate the effective error state by returning null */
	dev = Dictionary_Get( self->_meta, (char*)PLUGIN_DEPENDENCY_KEY );
	if( dev ) {
		return Dictionary_Entry_Value_AsDictionary( dev );
	}
	else {
		return NULL;
	}
}

Dictionary_Entry_Value* Module_GetValue( void* module, char* key ) {
	Module* self = (Module*)module;

	return Dictionary_Get( self->_meta, key );
}
	

void* Module_LoadSymbol( void* module, const char* suffix ) {
	Module* self = (Module*)module;
	char*   mangledName;
	char*   symbolText;
	void*   result;
	
	mangledName = self->MangleName( self->name );
	symbolText = Memory_Alloc_Array( char, strlen( mangledName ) + strlen( suffix ) + 2, MODULE_SYMBOLNAME );

	sprintf( symbolText, "%s%s",  mangledName, suffix );

	result = dlsym( self->dllPtr, symbolText );
	if( result == NULL ) {
		/* Try with a leading "_"... this is because on macx the dlcompat library can work either placing
		   this "_" for you and without and there is no easy way to know */
		sprintf( symbolText, "_%s%s", mangledName, suffix );
		result = dlsym( self->dllPtr, symbolText );
	}

	Memory_Free( symbolText );
	Memory_Free( mangledName );

	return result;
}


void Module_UnLoad( void* module ) {
	Module* self = (Module*)module;
	
	if( self->dllPtr ) {
		dlclose( self->dllPtr );
	}
	self->dllPtr = 0;
	self->GetMetadata = 0;
	self->GetName = 0;
	self->GetVersion = 0;
}

char* Module_MangledName( void* module ) {
	Module* self = (Module*)module;
	
	return self->MangleName( self->name );
}
