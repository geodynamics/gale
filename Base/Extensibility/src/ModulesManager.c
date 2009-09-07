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
** $Id: ModulesManager.c 4081 2007-04-27 06:20:07Z LukeHodkinson $
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
#include "ModulesManager.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#if NOSHARED
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

/* Textual name of this class */
const Type ModulesManager_Type = "ModulesManager";

/* For when compiling a single a single executable (no loading of modules at runtime */
/* This is mainly for platforms which does not support dynamic libs */
void SingleRegister();


const char* Codelet_GetMetadata();

/* Define memory for global pointer to moduleDirectories */
Stg_ObjectList*  moduleDirectories = NULL;	

ModulesManager* _ModulesManager_New(
		SizeT                                   _sizeOfSelf,
		Type                                    type,
		Stg_Class_DeleteFunction*               _delete,
		Stg_Class_PrintFunction*                _print,
		Stg_Class_CopyFunction*                 _copy, 
		ModulesManager_GetModulesListFunction*  _getModulesList,
		ModulesManager_LoadModuleFunction*	_loadModule,
		ModulesManager_UnloadModuleFunction*	_unloadModule,
		ModulesManager_ModuleFactoryFunction*   _moduleFactory,
		ModulesManager_CheckContextFunction*	_checkContext )
{
	ModulesManager* self;
	
	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(ModulesManager) );
	self = (ModulesManager*)_Stg_Class_New( _sizeOfSelf, type, _delete, _print, _copy );
	
	/* General info */
	
	/* Virtual info */
	self->_getModulesList = _getModulesList;
	self->_loadModule = _loadModule;
	self->_unloadModule = _unloadModule;
	self->_moduleFactory = _moduleFactory;
	self->_checkContext = _checkContext;
	
	_ModulesManager_Init( self );
	
	return self;
}

void _ModulesManager_Init( void* modulesManager ) {
	ModulesManager* self = (ModulesManager*)modulesManager;
	
	/* General and Virtual info should already be set */
	/* Modules info */

	self->modules = Stg_ObjectList_New();
	self->codelets = Stg_ObjectList_New();
}

void _ModulesManager_Delete( void* modulesManager ) {
	ModulesManager*         self = (ModulesManager*)modulesManager;
	LiveComponentRegister*  lcRegister = NULL;

	/* Note : We have to delete the codelets here rather than let the 
		component factory do it as they refer to static data inside the
		loaded module dlls, which are no longer available once we delete
		all in self->modules. -- Main.PatrickSunter 18 May 2006 */
	if ( (lcRegister = LiveComponentRegister_GetLiveComponentRegister()) ) {
		/* sweep through codelets list and delete each from the lcRegister */
		Index codelet_I;
		
		for ( codelet_I = 0; codelet_I < self->codelets->count; ++codelet_I ) {
			LiveComponentRegister_RemoveOneComponentsEntry( lcRegister,
				((Stg_Object*)self->codelets->data[codelet_I])->name );
		}
	}
	Stg_ObjectList_DeleteAllObjects( self->codelets );
	Stg_Class_Delete( self->codelets );
	ModulesManager_Unload( self ); 
	Stg_Class_Delete( self->modules );
	
	/* Delete parent */
	_Stg_Class_Delete( self );
}

void _ModulesManager_Print( void* modulesManager, Stream* stream ) {
	ModulesManager* self = (ModulesManager*)modulesManager;
	Index i;
	
	/* General info */
	Journal_Printf( (void*) stream, "%s (ptr): %p\n", ModulesManager_Type, self );
	
	/* Print parent */
	_Stg_Class_Print( self, stream );

	if ( moduleDirectories->count > 0 ) {
		Journal_Printf( stream, "Search Path:\n" );
		Stream_Indent( stream );
		for ( i = 0; i < moduleDirectories->count; ++i ) {
			Journal_Printf( stream, "(path) %s\n", Stg_ObjectList_ObjectAt( moduleDirectories, i ) );
		}
		Stream_UnIndent( stream );
	}

	if ( self->modules->count > 0 ) {
		Journal_Printf( stream, "Loaded modules:\n" );
		Stream_Indent( stream );
		for ( i = 0; i < self->modules->count; ++i ) {
			Module* module = (Module*)Stg_ObjectList_At( self->modules, i );
			Journal_Printf( stream, "%s\n", module->name );
			if ( Module_GetVersion( module ) ) {
				Journal_Printf( stream, "\t(version) %s\n", Module_GetVersion( module ) );
			}
		}
		Stream_UnIndent( stream );
	}
}

Dictionary_Entry_Value* ModulesManager_GetModulesList( void* modulesManager, void* _dictionary ) {
	ModulesManager*			self = (ModulesManager*)modulesManager;
	Dictionary*			dictionary = (Dictionary*)_dictionary;
	
	return self->_getModulesList( self, dictionary );
}

Bool ModulesManager_CheckContext( void* modulesManager, void* _dictionary, Name moduleName, Name contextName ) {
	ModulesManager*			self 		= (ModulesManager*)modulesManager;
	Dictionary*			dictionary 	= (Dictionary*)_dictionary;

	return self->_checkContext( self, dictionary, moduleName, contextName );
}

void ModulesManager_Load( void* modulesManager, void* _dictionary, Name contextName ) {
	ModulesManager*			self 		= (ModulesManager*)modulesManager;
	Dictionary*			dictionary 	= (Dictionary*)_dictionary;
	unsigned int			entryCount;
	unsigned int			entry_I;
	Dictionary_Entry_Value*		modulesVal;

	/* First add the directory list onto LD_LIBRARY_PATH so that it can potentially
	 * resolve the unknown symbols */
	char* curEnvPath;
	char* newEnvPath;
	Index newEnvPathLength = 0;
	Index dir_I;
	
#ifdef SINGLE_EXE
	SingleRegister();
#else
	if( dictionary ) {
		Dictionary_Entry_Value* localLibDirList = Dictionary_Get( dictionary, "LD_LIBRARY_PATH" );
		if( localLibDirList ) {
			Index i, count;
			char* dir;
			count = Dictionary_Entry_Value_GetCount( localLibDirList );
			for ( i = 0; i < count; ++i ) {
				dir = Dictionary_Entry_Value_AsString( Dictionary_Entry_Value_GetElement( localLibDirList, i ) );
				ModulesManager_AddDirectory( "FromDictionary", dir );
			}
		}
	}
	
	for ( dir_I = 0; dir_I < moduleDirectories->count; ++dir_I ) {
		newEnvPathLength += strlen( Stg_ObjectList_ObjectAt( moduleDirectories, dir_I ) );
		/* Add one make space for the ':' inbetween the directories */
		newEnvPathLength += 1; 
	}
	curEnvPath = getenv("LD_LIBRARY_PATH");
	if ( curEnvPath ) {
		newEnvPathLength += strlen( curEnvPath );
	}

	if ( newEnvPathLength > 0 ) {
		/* Add one to make space for the Null Terminator '\0' */
		newEnvPathLength += 1;
		
		newEnvPath = Memory_Alloc_Array( char, newEnvPathLength, "LD_LIBRARY_PATH" );
		newEnvPath[0] = '\0';
		for ( dir_I = 0; dir_I < moduleDirectories->count; ++dir_I ) {
			strcat( newEnvPath, Stg_ObjectList_ObjectAt( moduleDirectories, dir_I ) );
			strcat( newEnvPath, ":" );
		}
		if ( curEnvPath ) {
			strcat( newEnvPath, curEnvPath );
		}
		setenv( "LD_LIBRARY_PATH", newEnvPath, 1 );
		Journal_Printf(
			Journal_Register( Debug_Type, self->type ),
			"Using LD_LIBRARY_PATH=%s\n",
			newEnvPath );
		Memory_Free( newEnvPath );
	}

	modulesVal = ModulesManager_GetModulesList( self, dictionary );
	
	if( !modulesVal ) {
		entryCount = 0;
	}
	else {
		entryCount = Dictionary_Entry_Value_GetCount( modulesVal );
	}
	
	for( entry_I = 0; entry_I < entryCount; entry_I++ ) {
		Name		moduleName;
		moduleName = Dictionary_Entry_Value_AsString( Dictionary_Entry_Value_GetElement( modulesVal, entry_I ) );

		if( !ModulesManager_CheckContext( self, dictionary, moduleName, contextName ) )
			continue;

		if ( ! ModulesManager_LoadModule( self, moduleName ) ) {
			Journal_Firewall(
				0,
				Journal_Register( Info_Type, self->type ),
				"Error. Unable to load module %s\n"
				"Please check that:-\n"
				"   1. Module %s(%smodule.so) is built\n"
				"   2. Module has function 'Index %s_Register( ... )'\n"
				"   3. Module does not have undefined symbols.\n"
				"   4. $LD_LIBRARY_PATH is set\n"
				"   5. $LD_LIBRARY_PATH is set correctly\n"
				"You can get more information about what went wrong by adding the following two lines to your input file:\n"
				"	<param name=\"journal.debug\">true</param>\n"
				"	<param name=\"journal.debug.%s\">true</param>\n"
				"For more help visit http://csd.vpac.org/twiki/bin/view/Stgermain/FrequentlyAskedQuestions#Plugins_Won_t_Load_in_Parallel\n",
				moduleName, moduleName, moduleName, moduleName, self->type );
		}
	}
#endif /* ifdef SINGLE_EXE */
}

Bool ModulesManager_LoadModule( void* modulesManager, Name moduleName ) {
	ModulesManager*	self = (ModulesManager*)modulesManager;
	Module* theModule;

	Dictionary* metadata;
	Index depCount = 0;
	Index dep_I;
	
	Stream* stream;
	Stream* debug;
	Stream* error;

    
	stream =  Journal_Register( Info_Type, self->type );
	debug =  Journal_Register( Debug_Type, self->type );
	error =  Journal_Register( Error_Type, self->type );

	if ( Stg_ObjectList_Get( self->modules, moduleName ) != NULL ) {
		Journal_Printf( debug, "Module %s already loaded\n", moduleName );
		return True;
	}

	/* Try local search first */
	theModule = ModulesManager_CreateModule( self, moduleName, moduleDirectories );

	/* remote search? */

	if ( theModule == NULL ) {
		return False;
	}

	/* version check? */
	
	metadata = Module_GetMetadata( theModule );
	depCount = metadata ? Stg_Meta_GetAssociationCount( metadata ) : 0;
	for( dep_I = 0; dep_I < depCount; ++dep_I ) {
		char* name;

		name = Stg_Meta_GetAssociationType( metadata, dep_I );
		Journal_Firewall( (Bool)name, stream, "Module name/type not found in dependency/associations meta entry\n" );

		if( !ModulesManager_LoadModule( self, name ) ) {
			Journal_Printf( stream, "Dependency %s failed to load\n", name );
			Stg_Class_Delete( theModule );
			return False;
		}
	}
	
	/* Do the actual loading of the module */
	if( self->_loadModule( self, theModule ) ) {
		Stg_ObjectList_Append( self->modules, theModule );
	}
	else {
		return False;
	}

	return True;
}


Bool ModulesManager_UnloadModule( void* modulesManager, Name moduleName ) {
	ModulesManager* self = (ModulesManager*)modulesManager;
	
	if( self->_unloadModule( self, Stg_ObjectList_Get( self->modules, moduleName ) ) ) {
		if( Stg_ObjectList_Remove( self->modules, moduleName, DELETE ) ) {
			return True;
		}
	}
	return False;
}

void ModulesManager_Unload( void* modulesManager ) {
	ModulesManager* self = (ModulesManager*)modulesManager;
	
	while( self->modules->count ) {
		Module* module = (Module*)Stg_ObjectList_At( self->modules, self->modules->count - 1 ); /* reverse order deletion */
		
		ModulesManager_UnloadModule( self, module->name );
	}
}

void ModulesManager_AddDirectory( Name name, char* directory ) {
	Bool				found;
	Index				i;
	
	/* Check if global list is initialised */
	if (moduleDirectories == NULL) {
		moduleDirectories = Stg_ObjectList_New();
	}
	
	/* Check that the directory isn't already added */
	found = False;
	for( i =  0; i < moduleDirectories->count; i++ ){
		if( strcmp( directory, Stg_ObjectList_ObjectAt( moduleDirectories, i ) ) == 0 ) {
			found = True;
		}
	}
	
	if( !found ) {
		Stg_ObjectList_PointerAppend( moduleDirectories, StG_Strdup( directory ), name, 0, 0, 0 ); 
	}
}

/* This prototype is provided here as its generated by the build system from the meta files, and the call to 
   ModulesManager_Add in the function below is the raw one, and hence more info is needed */
Dictionary* Codelet_MetaAsDictionary();

Index ModulesManager_Submit( 
		void* modulesManager, 
		Name codeletName, 
		Name version,
		Stg_Component_DefaultConstructorFunction* defaultNew ) 
{
	ModulesManager* self = (ModulesManager*)modulesManager;

	void* codeletInstance = NULL;
	Index result;
	Stg_ComponentRegister* componentRegister;
	Stg_Component_DefaultConstructorFunction*  functionCheck=NULL;
	
	/* Check whether the module has already been added to the component register */
	
	componentRegister = Stg_ComponentRegister_Get_ComponentRegister();
	functionCheck = Stg_ComponentRegister_Get(componentRegister, codeletName, version);
	
	/* Add module if it isn't already there. */
	if( functionCheck == NULL ) {
		Stg_ComponentRegister_AddFunc( 
			Stg_ComponentRegister_Get_ComponentRegister(),
			codeletName,
			version,
			defaultNew,
			Codelet_MetaAsDictionary );
	
		codeletInstance = defaultNew( codeletName );
		result = Stg_ObjectList_Append( self->codelets, codeletInstance );
		if( LiveComponentRegister_GetLiveComponentRegister() ) {
			LiveComponentRegister_Add( LiveComponentRegister_GetLiveComponentRegister(), codeletInstance );
		}
		
		return result;
	}
	else {
		// TODO check if an existing self->codelets has the same Type - if so error & exit
		
		/* Otherwise since this Codelet has already been added to the global ComponentRegister,
		 * but not to this ModulesManager instance, then just instantiate another Codelet instance
		 * and add it to this ModulesManager's codelet array. */
		codeletInstance = defaultNew( codeletName );
		result = Stg_ObjectList_Append( self->codelets, codeletInstance );
		
		return result;
	}
}

Module* ModulesManager_CreateModule( void* modulesManager, Name name, Stg_ObjectList* directories ) {
	ModulesManager* self = (ModulesManager*)modulesManager;
	Module* module = self->_moduleFactory( name, directories );

	if( module->dllPtr ) {
		return module;
	}
	else {
		return 0;
	}
}


void ModulesManager_ConstructModules( void* modulesManager, Stg_ComponentFactory* cf, void* data ) {
	ModulesManager* self = (ModulesManager*)modulesManager;
	int i;

	for( i = 0; i < self->codelets->count; ++i ) {
		Stg_Component_Construct( self->codelets->data[i], cf, data, False );
	}
}


Bool ModulesManager_ConstructModule( void* modulesManager, Name moduleName, Stg_ComponentFactory* cf, void* data ) {
	ModulesManager* self = (ModulesManager*)modulesManager;

	/* Get the codelet for the module we're wanting to construct. */
	void* codelet = Stg_ObjectList_Get( self->codelets, moduleName );
				
	if( codelet == NULL )
		return False;
	else
		Stg_Component_Construct( codelet, cf, data, False );

	return True;
}
