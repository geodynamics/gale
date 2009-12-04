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
**	Handles the loading of "Toolbox" modules, which can extend the functionality or data structures of 
**	a main StGermain program.
**
** Assumptions:
**
** Comments:
**
** $Id: ToolboxesManager.h 4014 2007-02-23 02:15:16Z KathleenHumble $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __StGermain_Base_Extensibility_ToolboxesManager_h__
#define __StGermain_Base_Extensibility_ToolboxesManager_h__
	

	/* Textual name of this class */
	extern const Type ToolboxesManager_Type;

	
	/* Toolboxes info */
	#define __ToolboxesManager \
		/* General info */ \
		__ModulesManager \
		\
		/* Virtual info */ \
		\
		/* Toolboxes info */ \
		int*       argc; \
		char***    argv; \
		char**     initialised; \
		Index      _initialisedSize; \
		Index      _initialisedCount;
		
	struct ToolboxesManager { __ToolboxesManager };
	
    /** Define a global list of plugin directories*/
     extern Stg_ObjectList*  pluginDirectories;	

	/* Create a new Toolboxes */
	ToolboxesManager* ToolboxesManager_New( int* argc, char*** argv );
	
	/* Creation implementation / Virtual constructor */
	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define TOOLBOXESMANAGER_DEFARGS \
                MODULESMANAGER_DEFARGS, \
                int*     argc, \
                char***  argv

	#define TOOLBOXESMANAGER_PASSARGS \
                MODULESMANAGER_PASSARGS, \
	        argc, \
	        argv

	ToolboxesManager* _ToolboxesManager_New(  TOOLBOXESMANAGER_DEFARGS  );

	/* Initialisation implementation */
	void _ToolboxesManager_Init( void* toolboxesManager, int* argc, char*** argv );
	
	/* Stg_Class_Delete implementation */
	void _ToolboxesManager_Delete( void* toolboxesManager );
	
	/* Print implementation */
	void _ToolboxesManager_Print( void* toolboxesManager, Stream* stream );
	
	/** Get the toolbox list from the dictionary */
	Dictionary_Entry_Value* _ToolboxesManager_GetToolboxesList(  void* toolboxesManager, void* dictionary );
	
	/** Exactly what to do to load the toolbox */
	Bool _ToolboxesManager_LoadToolbox( void* toolboxesManager, Module* toolbox );

	/** Exactly what to do to unload the toolbox */
	Bool _ToolboxesManager_UnloadToolbox( void* toolboxesManager, Module* toolbox );

	Bool _ToolboxesManager_CheckContext( void* toolboxesManager, Dictionary_Entry_Value* modulesVal, unsigned int entry_I, Name contextName );

	Name _ToolboxesManager_GetModuleName( void* toolboxesManager, Dictionary_Entry_Value* moduleVal, unsigned int entry_I );

	#define ToolboxesManager_Submit ModulesManager_Submit

	/** Let StGermain know that the "Init" function of a module has been called. This exists to handle the case where a module
	   is linked into a binary and the user attempts to module load the module too. */
	Index ToolboxesManager_SetInitialised( void* toolboxesManager, char* label );
	
	/** This exists to handle the case where a module is linked into a binary and the user attempts to module load the module too. 
	   Its expected at modules will check to see if they have been inited already before doing initialisation work in the init 
	   function. */
	Bool ToolboxesManager_IsInitialised( void* toolboxesManager, char* label );

#endif /* __StGermain_Base_Extensibility_ToolboxesManager_h__ */

