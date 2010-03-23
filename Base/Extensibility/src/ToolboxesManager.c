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
	/* Variables set in this function */
	SizeT                                       _sizeOfSelf = sizeof(ToolboxesManager);
	Type                                               type = ToolboxesManager_Type;
	Stg_Class_DeleteFunction*                       _delete = _ToolboxesManager_Delete;
	Stg_Class_PrintFunction*                         _print = _ToolboxesManager_Print;
	Stg_Class_CopyFunction*                           _copy = NULL;
	ModulesManager_GetModulesListFunction*  _getModulesList = _ToolboxesManager_GetToolboxesList;
	ModulesManager_LoadModuleFunction*          _loadModule = _ToolboxesManager_LoadToolbox;
	ModulesManager_UnloadModuleFunction*      _unloadModule = _ToolboxesManager_UnloadToolbox;
	ModulesManager_ModuleFactoryFunction*    _moduleFactory = Toolbox_Factory;
	ModulesManager_CheckContextFunction*      _checkContext = _ToolboxesManager_CheckContext;
	ModulesManager_GetModuleNameFunction*    _getModuleName = _ToolboxesManager_GetModuleName;

	return _ToolboxesManager_New(  TOOLBOXESMANAGER_PASSARGS  );
}

ToolboxesManager* _ToolboxesManager_New(  TOOLBOXESMANAGER_DEFARGS  )
{
	ToolboxesManager* self;
	
	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(ToolboxesManager) );
	self = (ToolboxesManager*)_ModulesManager_New(  MODULESMANAGER_PASSARGS  );
	
	/* General info */
	
	/* Virtual info */
	
	_ToolboxesManager_Init( self, argc, argv );
	
	return self;
}

void _ToolboxesManager_Init( void* toolboxesManager, int* argc, char*** argv ) {
	ToolboxesManager*         self = (ToolboxesManager*)toolboxesManager;
	
   self->initTB = Stg_ObjectList_New();
	self->argc = argc;
	self->argv = argv;
}


void _ToolboxesManager_Delete( void* toolboxesManager ) {
   ToolboxesManager*         self = (ToolboxesManager*)toolboxesManager;
   int ii, originalListSize;

	Stg_ObjectList_DeleteAllObjects( self->codelets );
	Stg_Class_Delete( self->codelets );
	ModulesManager_Unload( self );  /* this will unload all toolboxes implicitly */
	Stg_Class_Delete( self->modules );
	
	/* Delete parent */
	_Stg_Class_Delete( self );
}

void _ToolboxesManager_Print( void* toolboxesManager, Stream* stream ) {
	ToolboxesManager* self = (ToolboxesManager*)toolboxesManager;
	
	/* General info */
	Journal_Printf( (void*) stream, "Toolboxes (ptr): %p\n", self );
	
	if( Stg_ObjectList_Count(self->initTB) > 0 ) {
		Index i;
		
		Journal_Printf( stream, "Initialised Modules:\n" );
		Stream_Indent( stream );
		for( i = 0; i < Stg_ObjectList_Count(self->initTB) ; ++i ) {
			Journal_Printf( stream, "%s\n", self->initTB->data[i]->name );
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
	
   /* if not Loaded call the Initialise() and Register() */
   if( !Stg_ObjectList_Get( self->initTB, toolbox->name ) ) {

      ((Toolbox*)toolbox)->Initialise( self, self->argc, self->argv );
      ((Toolbox*)toolbox)->Register( self );

      Stg_ObjectList_Append( self->initTB, toolbox );
   }
	return True;
}

Bool _ToolboxesManager_UnloadToolbox( void* toolboxesManager, Module* toolbox ) {
	ToolboxesManager* self = (ToolboxesManager*)toolboxesManager;
	
   if( Stg_ObjectList_Get( self->initTB, toolbox->name ) ) {
      ((Toolbox*)toolbox)->Finalise( self );

      /* remove the toolbox from the initTB list, but don't actually Delete it's memory */
      Stg_ObjectList_Remove( self->initTB, toolbox->name, KEEP );
   }
    
	return True;
}

/* toolboxes do not need to be associated with contexts, so just return true */
Bool _ToolboxesManager_CheckContext( void* toolboxesManager, Dictionary_Entry_Value* modulesVal, unsigned int entry_I, Name contextName ) {
	return True;
}

Name _ToolboxesManager_GetModuleName( void* toolboxesManager, Dictionary_Entry_Value* moduleVal, unsigned int entry_I ) {
	return Dictionary_Entry_Value_AsString( Dictionary_Entry_Value_GetElement( moduleVal, entry_I ) );
}

Bool ToolboxesManager_IsInitialised( void* initRegister, char* label ) {
	ToolboxesManager* self = (ToolboxesManager*)initRegister;

   if( Stg_ObjectList_Get( self->initTB, label ) )
      return True;
   else
      return False;
}



