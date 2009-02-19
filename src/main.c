/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003-2006, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street,
**	Melbourne, 3053, Australia.
**
** Primary Contributing Organisations:
**	Victorian Partnership for Advanced Computing Ltd, Computational Software Development - http://csd.vpac.org
**	Australian Computational Earth Systems Simulator - http://www.access.edu.au
**	Monash Cluster Computing - http://www.mcc.monash.edu.au
**	Computational Infrastructure for Geodynamics - http://www.geodynamics.org
**
** Contributors:
**	Patrick D. Sunter, Software Engineer, VPAC. (pds@vpac.org)
**	Robert Turnbull, Research Assistant, Monash University. (robert.turnbull@sci.monash.edu.au)
**	Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**	David May, PhD Student, Monash University (david.may@sci.monash.edu.au)
**	Louis Moresi, Associate Professor, Monash University. (louis.moresi@sci.monash.edu.au)
**	Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**	Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
**	Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
**	Julian Giordani, Research Assistant, Monash University. (julian.giordani@sci.monash.edu.au)
**	Vincent Lemiale, Postdoctoral Fellow, Monash University. (vincent.lemiale@sci.monash.edu.au)
**	Kent Humphries, Software Engineer, VPAC. (kenth@vpac.org)
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
** $Id: main.c 532 2006-04-04 00:21:59Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifdef HAVE_PYTHON
	#include <Python.h>
#endif
#include <mpi.h>
#include <StGermain/StGermain.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

const Type StGermain_Type = "StGermain";

static void deleteListArgItem( void* ptr ) {
	/* Correct way to delete result of stgParseList[All]CmdLineArg */
	Memory_Free( ptr );
}

int main( int argc, char* argv[] ) {
	/* StGermain standard bits & pieces */
	MPI_Comm			CommWorld;
	int				rank;
	int				numProcessors;
	Dictionary*			dictionary;
	XML_IO_Handler*			ioHandler;
	Stream*                         stream;
	char*                           helpTopic;
	char*                           listAllTopic;
	Stg_ObjectList*                 listAllTopics;
	char*                           listTopic;
	Stg_ObjectList*                 listTopics;

	/* Initialise PETSc, get world info */
	MPI_Init( &argc, &argv );
	MPI_Comm_dup( MPI_COMM_WORLD, &CommWorld );
	MPI_Comm_size( CommWorld, &numProcessors );
	MPI_Comm_rank( CommWorld, &rank );
	StGermain_Init( &argc, &argv );
	stream = Journal_Register( Info_Type, StGermain_Type );
	#ifdef HAVE_PYTHON
		Py_Initialize();
	#endif	
	MPI_Barrier( CommWorld ); /* Ensures copyright info always come first in output */
	
	
	/* Create the application's dictionary & read input */
	dictionary = Dictionary_New();
	ioHandler = XML_IO_Handler_New();
	IO_Handler_ReadAllFromCommandLine( ioHandler, argc, argv, dictionary );
	Journal_ReadFromDictionary( dictionary );
	
	/* Parse the additional (stgermain exe specific) commandline arguments */
	helpTopic = stgParseHelpCmdLineArg( &argc, &argv );
	listAllTopics = Stg_ObjectList_New();
	while( (listAllTopic = stgParseListAllCmdLineArg( &argc, &argv )) != 0 ) {
		Stg_ObjectList_Append( listAllTopics, Stg_ObjectAdaptor_NewOfPointer( listAllTopic, 0, True, False, deleteListArgItem, 0, 0 ) );
	}
	listTopics = Stg_ObjectList_New();
	while( (listTopic = stgParseListCmdLineArg( &argc, &argv )) != 0 ) {
		Stg_ObjectList_PointerAppend( listTopics, listTopic, 0, deleteListArgItem, 0, 0 );
	}

	/* if the command line arguments ask for "help" or "list-all", then load the toolboxes and plugins (to include all 
		symbols/components, and print the help for the selected component, but don't run the wholse app. There can be
		multiple "list-all"s. */
	if( helpTopic || Stg_ObjectList_Count( listAllTopics ) ) {
		PluginsManager* plugins = PluginsManager_New();
		Dictionary* metadata;
		Stg_ComponentRegister* cr = Stg_ComponentRegister_Get_ComponentRegister();
		Index i;

		ModulesManager_Load( stgToolboxesManager, dictionary );
		ModulesManager_Load( plugins, dictionary );

		if( helpTopic ) {
			metadata = Stg_ComponentRegister_GetMetadata( cr, helpTopic, "0" );
			if( metadata ) {
				Stg_Meta_Print( metadata, stream );
			}
			else {
				Journal_Printf( stream, "Help topic '%s' not found.\n", helpTopic );
			}
			Memory_Free( helpTopic );
			helpTopic = 0;
		}
		for( i = 0; i < Stg_ObjectList_Count( listAllTopics ); i++ ) {
			listAllTopic = (char*)Stg_ObjectAdaptor_Object( (Stg_ObjectAdaptor*)Stg_ObjectList_At( listAllTopics, i ) );

			if( strcmp( listAllTopic, "components" ) == 0 ) {
				BTreeIterator* i;
				Stg_ComponentRegisterElement* cre;

				Journal_Printf( stream, "Registered components are:\n" );
				i = Stg_ComponentRegister_GetIterator( cr );
				for( 
					cre = Stg_ComponentRegisterIterator_First( i ); 
					cre != NULL; 
					cre = Stg_ComponentRegisterIterator_Next( i ) )
				{
					Journal_Printf( stream, "\t'%s'\n", Stg_ComponentRegisterElement_GetType( cre ) );
				}

				Stg_Class_Delete( i );
			}
			else if( strcmp( listAllTopic, "references" ) == 0 ) {
				BTreeIterator* i;
				Stg_ComponentRegisterElement* cre;

				Journal_Printf( stream, "Component references are:\n" );
				i = Stg_ComponentRegister_GetIterator( cr );
				for( 
					cre = Stg_ComponentRegisterIterator_First( i ); 
					cre != NULL; 
					cre = Stg_ComponentRegisterIterator_Next( i ) )
				{
					char* reference = Stg_Meta_GetReference( Stg_ComponentRegisterElement_GetMetadata( cre ) );
					Journal_Printf( 
						stream, 
						"\t'%s': %s\n", 
						Stg_ComponentRegisterElement_GetType( cre ),
						( reference && reference[0] ) ? reference : "(None provided)" );
						/* i.e. if not null and not an empty string print the value else default */
				}

				Stg_Class_Delete( i );
			}
			else if( strcmp( listAllTopic, "equations" ) == 0 ) {
				BTreeIterator* i;
				Stg_ComponentRegisterElement* cre;

				Journal_Printf( stream, "Component equations are:\n" );
				i = Stg_ComponentRegister_GetIterator( cr );
				for( 
					cre = Stg_ComponentRegisterIterator_First( i ); 
					cre != NULL; 
					cre = Stg_ComponentRegisterIterator_Next( i ) )
				{
					char* equation = Stg_Meta_GetEquation( Stg_ComponentRegisterElement_GetMetadata( cre ) );
					Journal_Printf( 
						stream, 
						"\t'%s': %s\n", 
						Stg_ComponentRegisterElement_GetType( cre ),
						( equation && equation[0] ) ? equation : "(None provided)" );
						/* i.e. if not null and not an empty string print the value else default */
				}

				Stg_Class_Delete( i );
			}
			else if( strcmp( listAllTopic, "rights" ) == 0 ) {
				BTreeIterator* i;
				Stg_ComponentRegisterElement* cre;

				Journal_Printf( stream, "Component references are:\n" );
				i = Stg_ComponentRegister_GetIterator( cr );
				for( 
					cre = Stg_ComponentRegisterIterator_First( i ); 
					cre != NULL; 
					cre = Stg_ComponentRegisterIterator_Next( i ) )
				{
					char* rights = Stg_Meta_GetRights( Stg_ComponentRegisterElement_GetMetadata( cre ) );
					Journal_Printf( 
						stream, 
						"\t'%s': %s\n", 
						Stg_ComponentRegisterElement_GetType( cre ),
						( rights && rights[0] ) ? rights : "(None provided)" );
						/* i.e. if not null and not an empty string print the value else default */
				}

				Stg_Class_Delete( i );
			}
			else {
				Journal_Printf( stream, "List topic '%s' not found.\n", listAllTopic );
				Journal_Printf( stream, "Available topics are:\n\t'all-components'\n\t'all-references'\n\t'all-equations'\n\t'all-rights'\n", listAllTopic );
			}
		}

		/* metadata is provided as a reference - don't delete */
		Stg_Class_Delete( plugins );
	}
	else {  /* ... run the app */
		Index i;
		AbstractContext*                context;

		/* Magic happens here! */
		context = stgMainInit( dictionary, CommWorld );
		stgMainLoop( context );

		for( i = 0; i < Stg_ObjectList_Count( listTopics ); i++ ) {
			listTopic = (char*)Stg_ObjectAdaptor_Object( (Stg_ObjectAdaptor*)Stg_ObjectList_At( listTopics, i ) );

			if( strcmp( listTopic, "components" ) == 0 ) {
				Stg_ObjectList* uniqueComponentTypes;
				Index i;

				/* Add each instantiated component type to a list, ensuring the list is of unique entries */
				uniqueComponentTypes = Stg_ObjectList_New();
				for( i = 0; i < LiveComponentRegister_GetCount( stgLiveComponentRegister ); i ++ ) {
					Type componentType = Stg_Class_GetType( LiveComponentRegister_At( stgLiveComponentRegister, i ) );
					Index j;
					Bool found;

					found = False;
					for( j = 0; componentType && j < Stg_ObjectList_Count( uniqueComponentTypes ); j++ ) {
						char* added = (char*)Stg_ObjectAdaptor_Object( (Stg_ObjectAdaptor*)Stg_ObjectList_At( uniqueComponentTypes, j ) );
						if( strcmp( componentType, added ) == 0 ) {
							found = True;
						}
					}
					if( componentType && !found ) {
						Stg_ObjectList_PointerAppend( uniqueComponentTypes, componentType, 0, 0, 0, 0 );
					}
				}

				Journal_Printf( stream, "Instantiated/used components are:\n" );
				for( i = 0; i < Stg_ObjectList_Count( uniqueComponentTypes ); i ++ ) {
					Type componentType = (char*)Stg_ObjectAdaptor_Object( (Stg_ObjectAdaptor*)Stg_ObjectList_At( uniqueComponentTypes, i ) );
					Journal_Printf( stream, "\t'%s'\n", componentType );
				}
				
				Stg_Class_Delete( uniqueComponentTypes );
			}
		}

		stgMainDestroy( context );
	}
	
	Stg_Class_Delete( listAllTopics );
	Stg_Class_Delete( listTopics );
	Stg_Class_Delete( dictionary );
	
	/* Close off everything */
	#ifdef HAVE_PYTHON
		Py_Finalize();
	#endif
	StGermain_Finalise();
//	Journal_Printf( stream, "Finalised: StGermain Framework.\n");
	MPI_Finalize();
	
	return 0; /* success */
}
