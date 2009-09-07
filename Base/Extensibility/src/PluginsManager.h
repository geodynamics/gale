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
**	Handles the loading of "Plugin" modules, which can extend the functionality or data structures of 
**	a main StGermain program.
**
** Assumptions:
**
** Comments:
**
** $Id: PluginsManager.h 4163 2007-08-02 08:32:40Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __Base_Extensibility_PluginsManager_h__
#define __Base_Extensibility_PluginsManager_h__
	

	/* Textual name of this class */
	extern const Type PluginsManager_Type;

	
	/* Plugins info */
	#define __PluginsManager \
		/* General info */ \
		__ModulesManager \
		\
		/* Virtual info */ \
		\
		/* Plugins info */
		
	struct PluginsManager { __PluginsManager };
	
    /** Define a global list of plugin directories*/
     extern Stg_ObjectList*  pluginDirectories;	

	/* Create a new Plugins */
	PluginsManager* PluginsManager_New( void );
	
	/* Creation implementation / Virtual constructor */
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
		ModulesManager_CheckContextFunction*	_checkContext );
	
	/* Initialisation implementation */
	void _PluginsManager_Init( void* pluginsManager );
	
	/* Stg_Class_Delete implementation */
	void _PluginsManager_Delete( void* pluginsManager );
	
	/* Print implementation */
	void _PluginsManager_Print( void* pluginsManager, Stream* stream );
	
	/** Get the plugins list from the dictionary */
	Dictionary_Entry_Value* _PluginsManager_GetPluginsList( void* pluginsManager, void* dictionary );

	/** Exactly what to do to load the plugin */
	Bool _PluginsManager_LoadPlugin( void* pluginsManager, Module* plugin );
	
	/** Exactly what to do to unload the plugin */
	Bool _PluginsManager_UnloadPlugin( void* pluginsManager, Module* plugin );

	Bool _PluginsManager_CheckContext( void* pluginsManager, void* dictionary, Name pluginName, Name contextName );

	#define PluginsManager_Submit ModulesManager_Submit
	
#endif /* __Base_Extensibility_PluginsManager_h__ */
