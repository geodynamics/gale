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

#ifndef __StGermain_Base_Extensibility_PluginsManager_h__
#define __StGermain_Base_Extensibility_PluginsManager_h__
	

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
	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define PLUGINSMANAGER_DEFARGS \
                MODULESMANAGER_DEFARGS

	#define PLUGINSMANAGER_PASSARGS \
                MODULESMANAGER_PASSARGS

	PluginsManager* _PluginsManager_New(  PLUGINSMANAGER_DEFARGS  );
	
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

   /** unload all plugins, this includes dlclosing when dynamic libs are used */
   Bool PluginsManager_UnloadAll( void* pluginsManager );

	/** Remove all plugins from register */
   void PluginsManager_RemoveAllFromComponentRegister( void* pluginsManager );
	
	Bool _PluginsManager_CheckContext( void* pluginsManager, Dictionary_Entry_Value* modulesVal, unsigned int entry_I, Name contextName );

	Name _PluginsManager_GetModuleName( void* pluginsManager, Dictionary_Entry_Value* moduleVal, unsigned int entry_I );

	#define PluginsManager_Submit ModulesManager_Submit
	
#endif /* __StGermain_Base_Extensibility_PluginsManager_h__ */

