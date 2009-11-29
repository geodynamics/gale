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
** $Id: ModulesManager.h 4014 2007-02-23 02:15:16Z KathleenHumble $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __Base_Extensibility_ModulesManager_h__
#define __Base_Extensibility_ModulesManager_h__
	

	/* Textual name of this class */
	extern const Type ModulesManager_Type;
	
	/* Templates of virtual functions */
	typedef Dictionary_Entry_Value* (ModulesManager_GetModulesListFunction) ( void* modulesManager, void* dictionary );
	typedef Bool                    (ModulesManager_LoadModuleFunction)     ( void* modulesManager, Module* module );
	typedef Bool                    (ModulesManager_UnloadModuleFunction)   ( void* modulesManager, Module* module );
	typedef Module*                 (ModulesManager_ModuleFactoryFunction)  ( Name name, Stg_ObjectList* directories );
	typedef Bool			(ModulesManager_CheckContextFunction)	( void* modulesManager, 
										  Dictionary_Entry_Value* modulesVal, 
										  unsigned int entry_I, 
										  Name contextName );
	typedef Name			(ModulesManager_GetModuleNameFunction)	( void* modulesManager, Dictionary_Entry_Value* modulesVal, unsigned int entry_I );

	/* Modules info */
	#define __ModulesManager \
		/* General info */ \
		__Stg_Class \
		\
		/* Virtual info */ \
		ModulesManager_GetModulesListFunction* _getModulesList; \
		ModulesManager_LoadModuleFunction*     _loadModule; \
		ModulesManager_UnloadModuleFunction*   _unloadModule; \
		ModulesManager_ModuleFactoryFunction*  _moduleFactory; \
		ModulesManager_CheckContextFunction*   _checkContext; \
		ModulesManager_GetModuleNameFunction*  _getModuleName; \
		\
		/* Modules info */ \
		Stg_ObjectList* modules; \
		Stg_ObjectList* codelets; 
		
	struct ModulesManager { __ModulesManager };
	
	/** Define a global list of module directories*/
	extern Stg_ObjectList* moduleDirectories;	

	/* Creation implementation / Virtual constructor */
	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define MODULESMANAGER_DEFARGS \
                STG_CLASS_DEFARGS, \
                ModulesManager_GetModulesListFunction*  _getModulesList, \
                ModulesManager_LoadModuleFunction*          _loadModule, \
                ModulesManager_UnloadModuleFunction*      _unloadModule, \
                ModulesManager_ModuleFactoryFunction*    _moduleFactory, \
                ModulesManager_CheckContextFunction*      _checkContext, \
                ModulesManager_GetModuleNameFunction*    _getModuleName

	#define MODULESMANAGER_PASSARGS \
                STG_CLASS_PASSARGS, \
	        _getModulesList, \
	        _loadModule,     \
	        _unloadModule,   \
	        _moduleFactory,  \
	        _checkContext,   \
	        _getModuleName 

	ModulesManager* _ModulesManager_New(  MODULESMANAGER_DEFARGS  );
	
	/* Initialisation implementation */
	void _ModulesManager_Init( void* modulesManager );
	
	/* Stg_Class_Delete implementation */
	void _ModulesManager_Delete( void* modulesManager );
	
	/* Print implementation */
	void _ModulesManager_Print( void* modulesManager, Stream* stream );
	
	/** Get the modules list from the dictionary */
	Dictionary_Entry_Value* ModulesManager_GetModulesList( void* modulesManager, void* dictionary );

	/** check that a given module is associated with a given context (from the XML) - always true for toolboxes, 
	    whereas plugins are associated with a specific context */
	Bool ModulesManager_CheckContext( void* modulesManager, Dictionary_Entry_Value* modulesVal, unsigned int entry_I, Name contextName );

	Name ModulesManaget_GetModuleName( void* modulesManager, Dictionary_Entry_Value* moduleVal, unsigned int entry_I );
	
	/** Load the specified modules from the dictionary attaching the codelets (components) */ 
	void ModulesManager_Load( void* modulesManager, void* dictionary, Name contextName );

	/** Load a module attaching the codelet (component) */
	Bool ModulesManager_LoadModule( void* modulesManager, Name moduleName );

	/** Unload a module */
	Bool ModulesManager_UnloadModule( void* modulesManager, Name moduleName );
	
	/** Unload all the modules */
	void ModulesManager_Unload( void* modulesManager );
	
	/** Add a directory to the list of prefixes it will try if the "dlopen" fails on a module. It won't add the same dir
	    twice. */
	void ModulesManager_AddDirectory( Name name, char* directory );
	
	Index ModulesManager_Submit( 
			void* modulesManager, 
			Name codeletName, 
			Name version,
			Stg_Component_DefaultConstructorFunction* defaultNew );

	Module* ModulesManager_CreateModule( void* modulesManager, Name name, Stg_ObjectList* directories );

	void ModulesManager_ConstructModules( void* modulesManager, Stg_ComponentFactory* cf, void* data );

	Bool ModulesManager_ConstructModule( void* modulesManager, Name moduleName, Stg_ComponentFactory* cf, void* data );

#endif /* __Base_Extensibility_ModulesManager_h__ */

