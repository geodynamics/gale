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
** $Id: Stg_Component.c 4288 2008-07-24 05:32:42Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdarg.h>
#include "Base/Foundation/Foundation.h"
#include "Base/IO/IO.h"
#include "Base/Container/Container.h"

#include "types.h"
#include "shortcuts.h"
#include "Meta.h"
#include "Stg_Component.h"
#include "Stg_ComponentFactory.h"
#include "LiveComponentRegister.h"
#include "CallGraph.h"

#include "Stg_ComponentRegister.h"

#include <assert.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>


/* Textual name of this class */
const Type Stg_Component_Type = "Stg_Component";


Stg_Component* _Stg_Component_New( 
		SizeT						_sizeOfSelf, 
		Type						type,
		Stg_Class_DeleteFunction*			_delete,
		Stg_Class_PrintFunction*			_print, 
		Stg_Class_CopyFunction*				_copy, 
		Stg_Component_DefaultConstructorFunction*	_defaultConstructor,
		Stg_Component_ConstructFunction*		_construct,
		Stg_Component_BuildFunction*			_build,
		Stg_Component_InitialiseFunction*		_initialise,
		Stg_Component_ExecuteFunction*			_execute,
		Stg_Component_DestroyFunction*			_destroy,
		Name						name,
		AllocationType					nameAllocationType )
{
	Stg_Component* self;
	
	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(Stg_Component) );
	self = (Stg_Component*)_Stg_Object_New( _sizeOfSelf, type, _delete, _print, _copy, name, nameAllocationType );
	
	/* General info */
	
	/* Virtual functions */
	self->_defaultConstructor = _defaultConstructor;
	self->_construct = _construct;
	self->_build = _build;
	self->_initialise =  _initialise;
	self->_execute = _execute;
	self->_destroy = _destroy;
	
	/* Stg_Component info */
	_Stg_Component_Init( self );
	
	return self;
}


void _Stg_Component_Init( Stg_Component* self ) {
	/* General and Virtual info should already be set */
	
	/* Stg_Component info */
	self->isConstructed = False;
	self->isBuilt = False;
	self->isInitialised = False;
	self->hasExecuted = False;
	self->isDestroyed = False;
	self->constructType = Stg_Object_AppendSuffix( self, "Construct" );
	self->buildType = Stg_Object_AppendSuffix( self, "Build" );
	self->initialiseType = Stg_Object_AppendSuffix( self, "Initialise" );
	self->executeType = Stg_Object_AppendSuffix( self, "Execute" );
	self->destroyType = Stg_Object_AppendSuffix( self, "Destroy" );
}


void _Stg_Component_Delete( void* component ) {
	Stg_Component* self = (Stg_Component*)component;
	
   /* just make sure the component is destroyed */
 	Stg_Component_Destroy( self, NULL, False );
	
	/* Stg_Class_Delete parent class */
	_Stg_Object_Delete( self );
}


void _Stg_Component_Print( void* component, Stream* stream ) {
	Stg_Component* self = (Stg_Component*)component;
	
	/* General info */
	Journal_Printf( (void*) stream, "Stg_Component (ptr): %p\n", self );
	
	/* Print parent class */
	_Stg_Object_Print( self, stream );
	
	/* Virtual info */
	Journal_Printf( (void*) stream, "\t_defaultConstructor (func ptr): %p\n", self->_defaultConstructor );
	Journal_Printf( (void*) stream, "\t_construct (func ptr): %p\n", self->_construct );
	
	Journal_Printf( (void*) stream, "\t_initialise (func ptr): %p\n", self->_initialise );
	Journal_Printf( (void*) stream, "\t_execute (func ptr): %p\n", self->_execute );
	Journal_Printf( (void*) stream, "\t_destroy (func ptr): %p\n", self->_destroy );
	
	/* Stg_Component */
	Journal_Printf( (void*) stream, "\tisConstructed: %s\n", self->isConstructed ? "True" : "False" );
	Journal_Printf( (void*) stream, "\tisBuilt: %s\n", self->isBuilt ? "True" : "False" );
	Journal_Printf( (void*) stream, "\tisInitialised: %s\n", self->isInitialised ? "True" : "False" );
	Journal_Printf( (void*) stream, "\thasExecuted: %s\n", self->hasExecuted ? "True" : "False" );
	Journal_Printf( (void*) stream, "\tisDestroyed: %s\n", self->isDestroyed ? "True" : "False" );

	Journal_Printf( (void*) stream, "\tconstructor function name: %s\n", self->constructType );
	Journal_Printf( (void*) stream, "\tbuild function name: %s\n", self->buildType );
	Journal_Printf( (void*) stream, "\tinitialise function name: %s\n", self->initialiseType );
	Journal_Printf( (void*) stream, "\texecute function name: %s\n", self->executeType );
	Journal_Printf( (void*) stream, "\tdestroy function name: %s\n", self->destroyType );
}

void* _Stg_Component_Copy( void* component, void* dest, Bool deep, Name nameExt, struct PtrMap* ptrMap ) {
	Stg_Component*	self = (Stg_Component*)component;
	Stg_Component*	newComponent;
	
	newComponent = _Stg_Object_Copy( component, dest, deep, nameExt, ptrMap );
	
	/* Virtual methods */
	newComponent->_defaultConstructor = self->_defaultConstructor;
	newComponent->_construct = self->_construct;
	newComponent->_build = self->_build;
	newComponent->_initialise = self->_initialise;
	newComponent->_execute = self->_execute;
	newComponent->_destroy = self->_destroy;
	
	newComponent->isConstructed = self->isConstructed;
	newComponent->isBuilt = self->isBuilt;
	newComponent->isInitialised = self->isInitialised;
	newComponent->hasExecuted = self->hasExecuted;
	newComponent->isDestroyed = self->isDestroyed;
	
	newComponent->constructType = StG_Strdup( self->constructType );
	newComponent->buildType = StG_Strdup( self->buildType );
	newComponent->initialiseType = StG_Strdup( self->initialiseType );
	newComponent->executeType = StG_Strdup( self->executeType );
	newComponent->destroyType = StG_Strdup( self->destroyType );

	if ( LiveComponentRegister_GetLiveComponentRegister() ) {
		LiveComponentRegister_Add( LiveComponentRegister_GetLiveComponentRegister(), newComponent );
		Memory_CountInc( newComponent );	/* Live register now has a refernce to it */
	}

	return newComponent;
}

void Stg_Component_AssignFromXML( void* component, Stg_ComponentFactory* cf, void* data, Bool force ) {
	Stg_Component* self   = (Stg_Component*)component;

	assert( self );
	
	if( force || !self->isConstructed ) {
		/*Stg_TimeMonitor* tm;
		Stg_MemMonitor*  mm;
		char*            buf;*/
		
		Stream* stream = Journal_Register( Info_Type, "Stg_ComponentFactory" );

		Journal_RPrintf( stream, "Constructing %s..\n", self->name );
		Stream_Indent( stream );
		
		/*buf = Stg_Object_AppendSuffix( self, "Construct" );
		tm = Stg_TimeMonitor_New( buf, True, True, MPI_COMM_WORLD );
		mm = Stg_MemMonitor_New( buf, True, True, MPI_COMM_WORLD );
	
		Stg_MemMonitor_Begin( mm );
		Stg_TimeMonitor_Begin( tm );*/
		
		#ifdef USE_PROFILE
			Stg_CallGraph_Push( stgCallGraph, self->_construct, self->constructType );
		#endif
		
		if ( !self->isConstructed ) {
			Stg_Component_SetupStreamFromDictionary( 
				self, 
				Dictionary_GetDictionary( cf->componentDict, self->name ) );
		}
		
		self->isConstructed = True;
		self->_construct( self, cf, data );
		
		#ifdef USE_PROFILE
			Stg_CallGraph_Pop( stgCallGraph );
		#endif
		
		/*Stg_TimeMonitor_End( tm );
		Stg_MemMonitor_End( mm );
	
		Stg_TimeMonitor_Delete( tm );
		Stg_MemMonitor_Delete( mm );

		Memory_Free( buf );*/
		
		Stream_UnIndent( stream );
	}
}

void Stg_Component_Build( void* component, void* data, Bool force ) {
	Stg_Component* self = (Stg_Component*)component;
	
	if( force || !self->isBuilt ) {
		/*Stg_TimeMonitor* tm;
		Stg_MemMonitor*  mm;
		char*            buf;
		
		buf = Stg_Object_AppendSuffix( self, "Build" );
		
		tm = Stg_TimeMonitor_New( buf, True, True, MPI_COMM_WORLD );
		mm = Stg_MemMonitor_New( buf, True, True, MPI_COMM_WORLD );
	
		Stg_MemMonitor_Begin( mm );
		Stg_TimeMonitor_Begin( tm );*/
		
		#ifdef USE_PROFILE
			Stg_CallGraph_Push( stgCallGraph, self->_build, self->buildType );
		#endif
		
		self->_build( self, data );
		self->isBuilt = True;
		
		#ifdef USE_PROFILE
			Stg_CallGraph_Pop( stgCallGraph );
		#endif
		
		/*Stg_TimeMonitor_End( tm );
		Stg_MemMonitor_End( mm );
	
		Stg_TimeMonitor_Delete( tm );
		Stg_MemMonitor_Delete( mm );

		Memory_Free( buf );*/
	}
}

void Stg_Component_Initialise( void* component, void* data, Bool force ) {
	Stg_Component* self = (Stg_Component*)component;

	if( force || !self->isInitialised ) {
		/*Stg_TimeMonitor* tm;
		Stg_MemMonitor*  mm;
		char*            buf;
		
		buf = Stg_Object_AppendSuffix( self, "Initialise" );
		
		tm = Stg_TimeMonitor_New( buf, True, True, MPI_COMM_WORLD );
		mm = Stg_MemMonitor_New( buf, True, True, MPI_COMM_WORLD );
	
		Stg_MemMonitor_Begin( mm );
		Stg_TimeMonitor_Begin( tm );*/
		
		#ifdef USE_PROFILE
			Stg_CallGraph_Push( stgCallGraph, self->_initialise, self->initialiseType );
		#endif
		
		self->isInitialised = True;
		self->_initialise( self, data );
		
		
		#ifdef USE_PROFILE
			Stg_CallGraph_Pop( stgCallGraph );
		#endif
		
		/*Stg_TimeMonitor_End( tm );
		Stg_MemMonitor_End( mm );
	
		Stg_TimeMonitor_Delete( tm );
		Stg_MemMonitor_Delete( mm );
		
		Memory_Free( buf );*/
	}
}

void Stg_Component_Execute( void* component, void* data, Bool force ) {
	Stg_Component* self = (Stg_Component*)component;
	if( force || !self->hasExecuted ) {
		/*Stg_TimeMonitor* tm;
		Stg_MemMonitor*  mm;
		char*            buf;
		
		buf = Stg_Object_AppendSuffix( self, "Execute" );
		
		tm = Stg_TimeMonitor_New( buf, True, True, MPI_COMM_WORLD );
		mm = Stg_MemMonitor_New( buf, True, True, MPI_COMM_WORLD );
	
		Stg_MemMonitor_Begin( mm );
		Stg_TimeMonitor_Begin( tm );*/
		
		#ifdef USE_PROFILE
			Stg_CallGraph_Push( stgCallGraph, self->_execute, self->executeType );
		#endif
		
		self->_execute( self, data );
		self->hasExecuted = True;
		
		#ifdef USE_PROFILE
			Stg_CallGraph_Pop( stgCallGraph );
		#endif
		
		/*Stg_TimeMonitor_End( tm );
		Stg_MemMonitor_End( mm );
	
		Stg_TimeMonitor_Delete( tm );
		Stg_MemMonitor_Delete( mm );	
		
		Memory_Free( buf );*/
	}
}

void Stg_Component_Destroy( void* component, void* data, Bool force ) {
	Stg_Component* self = (Stg_Component*)component;

	if( /* force TODO: remove this comment after pcu_rejig period|| */ !self->isDestroyed ) {
		/*Stg_TimeMonitor*	tm;
		Stg_MemMonitor*		mm;
		char*						buf;
		
		buf = Stg_Object_AppendSuffix( self, "Destroy" );
		
		tm = Stg_TimeMonitor_New( buf, True, True, MPI_COMM_WORLD );
		mm = Stg_MemMonitor_New( buf, True, True, MPI_COMM_WORLD );
	
		Stg_MemMonitor_Begin( mm );
		Stg_TimeMonitor_Begin( tm );*/
		
		#ifdef USE_PROFILE
			Stg_CallGraph_Push( stgCallGraph, self->_destroy, self->destroyType );
		#endif
		
		self->isDestroyed = True;

      if( self->_destroy )
         self->_destroy( self, data );

		Memory_Free( self->destroyType );
		Memory_Free( self->executeType );
		Memory_Free( self->initialiseType );
		Memory_Free( self->buildType );
		Memory_Free( self->constructType );

		#ifdef USE_PROFILE
			Stg_CallGraph_Pop( stgCallGraph );
		#endif
		
		/*Stg_TimeMonitor_End( tm );
		Stg_MemMonitor_End( mm );
	
		Stg_TimeMonitor_Delete( tm );
		Stg_MemMonitor_Delete( mm );
		Memory_Free( buf );*/
	}
}

/* TODO: UPT TO HERE */

Bool Stg_Component_IsConstructed( void* component ) {
	Stg_Component* self = (Stg_Component*)component;
	
	return self->isConstructed;
}

Bool Stg_Component_IsBuilt( void* component ) {
	Stg_Component* self = (Stg_Component*)component;
	
	return self->isBuilt;
}

Bool Stg_Component_IsInitialised( void* component ) {
	Stg_Component* self = (Stg_Component*)component;
	
	return self->isInitialised;
}

Bool Stg_Component_HasExecuted( void* component ) {
	Stg_Component* self = (Stg_Component*)component;
	
	return self->hasExecuted;
}

Bool Stg_Component_IsDestroyed( void* component ) {
	Stg_Component* self = (Stg_Component*)component;
	
	return self->isDestroyed;
}


void Stg_Component_SetupStreamFromDictionary( void* component, Dictionary* dictionary ) {
	Stg_Component* self = (Stg_Component*)component;

	Index index;
	Bool valid;

	char* keyCopy;
	char* operation;
	char* streamType;

	Stream* stream;

	if ( dictionary == NULL ) {
		return;
	}
	
	/* Iterate through the whole dictionary, checking for journal related commands */
	for ( index = 0; index < dictionary->count; ++index ) {
		/* Check to see if it is a journal operation by searching for the JOURNAL_KEY. */
		if ( (strncasecmp( dictionary->entryPtr[index]->key, JOURNAL_KEY, strlen(JOURNAL_KEY))) == 0 ) {
			
			keyCopy = StG_Strdup( dictionary->entryPtr[index]->key );

			operation = strtok( keyCopy, JOURNAL_DELIMITER );
			streamType = strtok( NULL, JOURNAL_DELIMITER );

			if ( streamType == NULL ) {
				streamType = Info_Type;
			}

			stream = Journal_MyStream( streamType, self );

			valid = True;
			if ( strcmp( operation, JOURNAL_ENABLE_KEY ) == 0 ) {
				Bool enable = Dictionary_Entry_Value_AsBool( dictionary->entryPtr[index]->value );
				Stream_Enable( stream, enable );
			}
			else if ( strcmp( operation, JOURNAL_ENABLE_BRANCH_KEY ) == 0 ) {
				Bool enable = Dictionary_Entry_Value_AsBool( dictionary->entryPtr[index]->value );
				Stream_EnableBranch( stream, enable );
			}
			else if ( strcmp( operation, JOURNAL_LEVEL_KEY ) == 0 ) {
				JournalLevel level = Dictionary_Entry_Value_AsUnsignedInt( dictionary->entryPtr[index]->value );
				Stream_SetLevelBranch( stream, level );
			}
			else if ( strcmp( operation, JOURNAL_RANK_KEY ) == 0 ) {
				int rank = Dictionary_Entry_Value_AsUnsignedInt( dictionary->entryPtr[index]->value );
				Stream_SetPrintingRank( stream, rank );
			}
			else if ( strcmp( operation, JOURNAL_FLUSH_KEY ) == 0 ) {
				Bool flush = Dictionary_Entry_Value_AsBool( dictionary->entryPtr[index]->value );
				Stream_SetAutoFlush( stream, flush );
			}
			else {
				valid = False;
			}

			if ( !valid ) {
				Journal_Printf(
					Journal_Register( Info_Type, "Journal" ),
					"Warning - unknown journal operation %s for component %s\n",
					dictionary->entryPtr[index]->key,
					self->name );
			}
						
			Memory_Free( keyCopy );
		}
	}
											
}

