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
** $Id: PluginsManager.c 4163 2007-08-02 08:32:40Z SteveQuenette $
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
#include "Plugin.h"
#include "ModulesManager.h"
#include "PluginsManager.h"


/* Textual name of this class */
const Type PluginsManager_Type = "PluginsManager";


PluginsManager* PluginsManager_New( void ) {
	return _PluginsManager_New( 
		sizeof(PluginsManager), 
		PluginsManager_Type, 
		_PluginsManager_Delete, 
		_PluginsManager_Print, 
		NULL, 
		_PluginsManager_GetPluginsList,
		_PluginsManager_LoadPlugin,
		_PluginsManager_UnloadPlugin,
		Plugin_Factory,
		_PluginsManager_CheckContext,
		_PluginsManager_GetModuleName );
}

PluginsManager* _PluginsManager_New(
		SizeT                                   _sizeOfSelf,
		Type                                    type,
		Stg_Class_DeleteFunction*               _delete,
		Stg_Class_PrintFunction*                _print,
		Stg_Class_CopyFunction*                 _copy, 
		ModulesManager_GetModulesListFunction*  _getModulesList,
		ModulesManager_LoadModuleFunction*	_loadModule,
		ModulesManager_UnloadModuleFunction*	_unloadModule,
		ModulesManager_ModuleFactoryFunction*   _moduleFactory,
		ModulesManager_CheckContextFunction*	_checkContext,
		ModulesManager_GetModuleNameFunction*	_getModuleName )
{
	PluginsManager* self;
	
	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(PluginsManager) );
	self = (PluginsManager*)_ModulesManager_New( 
		_sizeOfSelf, 
		type, 
		_delete, 
		_print, 
		_copy, 
		_getModulesList, 
		_loadModule, 
		_unloadModule,
		_moduleFactory,
		_checkContext,
		_getModuleName );
	
	/* General info */
	
	/* Virtual info */
	
	_PluginsManager_Init( self );
	
	return self;
}

void _PluginsManager_Init( void* pluginsManager ) {
}


void _PluginsManager_Delete( void* pluginsManager ) {
	PluginsManager*         self = (PluginsManager*)pluginsManager;

	/* Delete parent */
	_ModulesManager_Delete( self );
}

void _PluginsManager_Print( void* pluginsManager, Stream* stream ) {
	PluginsManager* self = (PluginsManager*)pluginsManager;
	
	/* General info */
	Journal_Printf( (void*) stream, "Plugins (ptr): %p\n", self );
	
	/* Print parent */
	_ModulesManager_Print( self, stream );
}


Dictionary_Entry_Value* _PluginsManager_GetPluginsList( void* pluginsManager, void* _dictionary ) {
	/*PluginsManager*			self = (PluginsManager*)pluginsManager;*/
	Dictionary*			dictionary = (Dictionary*)_dictionary;
	Dictionary_Entry_Value*		pluginsList = NULL;
	
	pluginsList = Dictionary_Get( dictionary, "plugins" );

	return pluginsList;
}


Bool _PluginsManager_LoadPlugin( void* pluginsManager, Module* plugin ) {
	PluginsManager* self = (PluginsManager*)pluginsManager;
	
	((Plugin*)plugin)->Register( self );
    
	return True;
}

Bool _PluginsManager_UnloadPlugin( void* pluginsManager, Module* plugin ) {
	/*PluginsManager* self = (PluginsManager*)pluginsManager;*/
	
	return True;
}

Bool _PluginsManager_CheckContext( void* pluginsManager, Dictionary_Entry_Value* modulesVal, unsigned int entry_I, Name contextName ) {
	PluginsManager* 	self 		= (PluginsManager*)pluginsManager;
	Dictionary_Entry_Value*	pluginDEV	= Dictionary_Entry_Value_GetElement( modulesVal, entry_I );
	Dictionary*		pluginDict;
	Name			componentName;

	pluginDict = Dictionary_Entry_Value_AsDictionary( pluginDEV );
	if( !pluginDict )
		return False;

	componentName = Dictionary_GetString_WithDefault( pluginDict, "Context", "context" );

	if( !strcmp( componentName, contextName ) )
		return True;

	return False;
}

Name _PluginsManager_GetModuleName( void* pluginsManager, Dictionary_Entry_Value* moduleVal, unsigned int entry_I ) {
	PluginsManager* 	self 		= (PluginsManager*)pluginsManager;
	Dictionary_Entry_Value*	pluginDEV	= Dictionary_Entry_Value_GetElement( moduleVal, entry_I );
	Dictionary*		pluginDict 	= Dictionary_Entry_Value_AsDictionary( pluginDEV );
	Name			pluginName 	= Dictionary_GetString( pluginDict, "Type" );

	return pluginName;	
}


