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

#ifndef NOSHARED
#include <dlfcn.h>
#else
extern int stg_num_modules;
extern char **stg_module_names;
extern int *stg_num_module_syms;
extern char ***stg_module_syms;
extern void ***stg_module_funcs;

Bool Module_CheckStaticName(Module *self);
void* Module_MapStaticSymbol(Module *self, const char *sym);
#endif

#include <limits.h>


const Type Module_Type = "Module";

static Name MODULE_GETMETADATA_SUFFIX = "_MetaAsDictionary";
static Name MODULE_GETNAME_SUFFIX = "_GetName";
static Name MODULE_GETVERSION_SUFFIX = "_GetVersion";

#ifdef MEMORY_STATS
	static Name MODULE_FILENAME = "fileName";
	static Name MODULE_SYMBOLNAME = "symbolName";
#endif

#ifdef NOSHARED
int stg_num_modules = 0;
char **stg_module_names = NULL;
int *stg_num_module_syms = NULL;
char ***stg_module_syms = NULL;
void ***stg_module_funcs = NULL;
#endif

Module* _Module_New(  MODULE_DEFARGS  )
{
	Module* self;

	assert( _sizeOfSelf >= sizeof(Module) );

	/* The following terms are parameters that have been passed into this function but are being set before being passed onto the parent */
	/* This means that any values of these parameters that are passed into this function are not passed onto the parent function
	   and so should be set to ZERO in any children of this class. */
	nameAllocationType = NON_GLOBAL;

	self = (Module*)_Stg_Object_New(  STG_OBJECT_PASSARGS  );
	
	_Module_Init( self, MangleName, directories );

	return self;
}
	
void _Module_Init(
		Module*                      self,
		Module_MangleNameFunction    MangleName,
		Stg_ObjectList*              directories )
{
	char*                           fileName = NULL;
	char*                           fullPathName = NULL;

#ifndef NOSHARED
	int                             fullPathLength = 0;
	int                             length;
	Index                           dir_i;
	
	Stream *stream =  Journal_Register( Info_Type, self->type );
	Stream *error =  Journal_Register( Error_Type, self->type );
#endif
	
	Stream *debug =  Journal_Register( Debug_Type, self->type );
	self->MangleName = MangleName;
	
	
	Journal_Printf( debug, "Finding module: \"%s\"... ", self->name );

	/* Try the plugin name by itself (allows LD_LIBRARY_PATH) to take precendence */
	self->mangledName = self->MangleName( self->name );
#ifndef NOSHARED
	fileName = Memory_Alloc_Array( char, strlen(self->mangledName) + strlen(MODULE_SUFFIX) + strlen(MODULE_EXT) +1, MODULE_FILENAME );
	sprintf( fileName, "%s%s%s", self->mangledName, MODULE_SUFFIX, MODULE_EXT );
		
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
#else
	/* This is a static build, so we just need to make sure we can find the name
	   of the module in the module list. */
	self->dllPtr = (DLL_Handle)Module_CheckStaticName(self);
#endif

	/* Load the symbols */
	if( self->dllPtr ) {
		self->GetMetadata = (Module_GetMetadataFunction*)Module_LoadSymbol( self, MODULE_GETMETADATA_SUFFIX );
		self->GetName = (Module_GetNameFunction*)Module_LoadSymbol( self, MODULE_GETNAME_SUFFIX );
		self->GetVersion = (Module_GetVersionFunction*)Module_LoadSymbol( self, MODULE_GETVERSION_SUFFIX );
	}
	
	/* Load the meta data */
	if( self->GetMetadata ) {
		self->_meta = self->GetMetadata();
	}
	else {
		self->_meta = Dictionary_New();
	}

	if ( fullPathName ) {
		Memory_Free( fullPathName );
	}
	Memory_Free( fileName );
}
	
void _Module_Delete( void* module ) {
	Module* self = (Module*)module;

	Module_UnLoad( self );
	Stg_Class_Delete( self->_meta );
	Memory_Free( self->mangledName );
	
	/* Delete parent */
	_Stg_Object_Delete( self );
}
	
void _Module_Print( void* module, Stream* stream ) {
	Module* self = (Module*)module;

	Index count = 0;
	Index i;
	Name version;

	Journal_RPrintf( stream, "Module: %s\n", self->name );
	Stream_Indent( stream );
	
	version = Module_GetVersion( self );
	if ( version == NULL ) {
		version = "Unknown";
	}
	Journal_Printf( stream, "Version: (version) %s\n", version );

	count = Stg_Meta_GetAssociationCount( self->_meta );
	for( i = 0; i < count; ++i ) {
		char* name;
		
		name = Stg_Meta_GetAssociationType( self->_meta, i );
		Journal_RPrintf( stream, "Depends on: %s\n", name ? name : "(type missing - erronous meta data)" );
	}
	
	/* Print parent */
	_Stg_Object_Print( self, stream );
	
	Stream_UnIndent( stream );
}


Dictionary* Module_GetMetadata( void* module ) {
	Module* self = (Module*)module;
	
	return self->_meta;
}

Name Module_GetName( void* module ) {
	Module* self = (Module*)module;

	if ( self->GetName ) {
		return self->GetName();
	}
	return self->name;
}

Name Module_GetVersion( void* module ) {
	Module* self = (Module*)module;

	if ( self->GetVersion ) {
		return self->GetVersion();
	}
	return NULL;
}


void* Module_LoadSymbol( void* module, Name suffix ) {
	Module* self = (Module*)module;
	char*   mangledName;
	char*   symbolText;
	void*   result;
	
	mangledName = self->MangleName( self->name );
	symbolText = Memory_Alloc_Array( char, strlen( mangledName ) + strlen( suffix ) + 3, MODULE_SYMBOLNAME );

	sprintf( symbolText, "%s%s",  mangledName, suffix );

#ifndef NOSHARED
	result = dlsym( self->dllPtr, symbolText );
#else
	/* Try to map the symbol. */
	result = Module_MapStaticSymbol(self, symbolText);
#endif

	if( result == NULL ) {
		/* Try with a leading "_"... this is because on macx the dlcompat library can work either placing
		   this "_" for you and without and there is no easy way to know */
		sprintf( symbolText, "_%s%s", mangledName, suffix );
#ifndef NOSHARED
		result = dlsym( self->dllPtr, symbolText );
#else
		result = Module_MapStaticSymbol(self, symbolText);
#endif
	}

	Memory_Free( symbolText );
	Memory_Free( mangledName );

	return result;
}


void Module_UnLoad( void* module ) {
	Module* self = (Module*)module;
	
#ifndef NOSHARED
	if( self->dllPtr ) {
		dlclose( self->dllPtr );
	}
#endif
	self->dllPtr = 0;
	self->GetMetadata = 0;
	self->GetName = 0;
	self->GetVersion = 0;
}

char* Module_MangledName( void* module ) {
	Module* self = (Module*)module;
	
	return self->MangleName( self->name );
}

#ifdef NOSHARED
Bool Module_CheckStaticName(Module *self) {
    int ii;

    for(ii = 0; ii < stg_num_modules; ii++) {
	if(!strcmp(self->mangledName, stg_module_names[ii]))
	    break;
    }
    if(ii == stg_num_modules)
	return False;
    return True;
}

void* Module_MapStaticSymbol(Module *self, const char *sym) {
    int ii, jj;

    for(ii = 0; ii < stg_num_modules; ii++) {
	if(!strcmp(self->mangledName, stg_module_names[ii]))
	    break;
    }
    if(ii == stg_num_modules)
	return NULL;

    for(jj = 0; jj < stg_num_module_syms[ii]; jj++) {
	if(!strcmp(sym, stg_module_syms[ii][jj]))
	    break;
    }
    if(jj == stg_num_module_syms[ii])
	return NULL;

    return stg_module_funcs[ii][jj];
}
#endif




