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

#include <mpi.h>
#include "StGermain/Base/Base.h"

#include "main.h"

#include <stdio.h>


void stgMain( Dictionary* dictionary, MPI_Comm communicator ) {
	Stg_ComponentFactory*	cf;

	cf = stgMainConstruct( dictionary, communicator, NULL );
	stgMainLoop( cf );
	stgMainDestroy( cf );

	Stg_Class_Delete( cf );
}

/* TODO: need to find a way to add different communicators for different contexts */
Stg_ComponentFactory* stgMainConstruct( Dictionary* dictionary, MPI_Comm communicator, void* _context ) {
	Stg_ComponentFactory*	cf;
	Dictionary*					componentDict;
	Stg_Component*				component;
	AbstractContext*			context;
	unsigned						component_I;

	if( ( componentDict = Dictionary_Entry_Value_AsDictionary( Dictionary_Get( dictionary, "components" ) ) ) == NULL )
		componentDict = Dictionary_New();
	
	CheckDictionaryKeys( componentDict, "Component dictionary must have unique names\n" );
	cf = Stg_ComponentFactory_New( dictionary, componentDict );

	if( _context ) {
		context = (AbstractContext*)_context;
		context->CF = cf;
		context->dictionary = dictionary;
		context->communicator = communicator;
		LiveComponentRegister_Add( cf->LCRegister, (void*)context );
	}

	/* Instantion phase -------------------------------------------------------------------------------------------------*/
	Stg_ComponentFactory_CreateComponents( cf );

	/* Assign the dictionary, componentFactory & the communicator for the contexts */
	/* TODO: if different contexts require different communicators, then StG. components will be required for these, and 
	 * they should be passed in from the XML */
	/* Also, this is a little hacky, as nothing is known about the other layers of StG or their associated contexts here */
	for( component_I = 0; component_I < LiveComponentRegister_GetCount( cf->LCRegister ); component_I++ ) {
		component = LiveComponentRegister_At( cf->LCRegister, component_I );
		if( Stg_CheckContext( component, AbstractContext ) ) { 
			Journal_Firewall( dictionary->count, Journal_Register( Error_Type, "Error Stream" ), 
				 "Error in %s: The dictionary is empty, meaning no input parameters have been feed into your program. Perhaps you've forgot to pass any input files ( or command-line arguments ) in.\n", __func__); 	
			context = (AbstractContext*)component;
			context->dictionary = dictionary;
			context->CF = cf;
			context->communicator = communicator;
		}
	}
	/* Construction phase -----------------------------------------------------------------------------------------------*/
	Stg_ComponentFactory_ConstructComponents( cf, NULL );
	
	return cf;
}

void stgMainBuildAndInitialise( Stg_ComponentFactory* cf ) {
	/* Building phase ---------------------------------------------------------------------------------------------------*/
	Stg_ComponentFactory_BuildComponents( cf, NULL );
	
	/* Initialisaton phase ----------------------------------------------------------------------------------------------*/
	Stg_ComponentFactory_InitialiseComponents( cf, NULL );
}

Stg_ComponentFactory* stgMainInitFromXML( char* xmlInputFilename, MPI_Comm communicator, void* _context ) {
   Dictionary*       		dictionary = NULL;
   Bool              		result;
   XML_IO_Handler*   		ioHandler;
   Stg_ComponentFactory*	cf;

   dictionary = Dictionary_New();
   ioHandler = XML_IO_Handler_New();
   result = IO_Handler_ReadAllFromFile( ioHandler, xmlInputFilename, dictionary );
   /* In case the user has put any journal configuration in the XML, read here */
   Journal_ReadFromDictionary( dictionary );

   cf = stgMainConstruct( dictionary, communicator, _context );
   
   /* We don't need the XML IO handler again (however don't delete the dictionary as it's 
    * 'owned' by the context from hereon */
   Stg_Class_Delete( ioHandler );

   return cf;
}

	
void stgMainLoop( Stg_ComponentFactory* cf ) {
	/* Run (Solve) phase ------------------------------------------------------------------------------------------------*/
	/* do this by running the contexts, which manage the entry points which call the _Execute() funcs for the other components */
	unsigned component_i;
	Stg_Component* component;
	AbstractContext* context;
	
	for( component_i = 0; component_i < LiveComponentRegister_GetCount( cf->LCRegister ); component_i++ ) {
		component = LiveComponentRegister_At( cf->LCRegister, component_i );
		if( Stg_CheckContext( component, AbstractContext ) ) { 
			context = (AbstractContext*)component;
			AbstractContext_Dump( context );
			Stg_Component_Execute( context, 0, True );
		}
	}
}

void stgMainDestroy( Stg_ComponentFactory* cf ) {
	/* Destruct phase ---------------------------------------------------------------------------------------------------*/
	Stg_ComponentFactory_DestroyComponents( cf, NULL );
	Stg_Class_Delete( cf->LCRegister );
	_Stg_Class_Delete( cf );
}

void stgImportToolbox( Dictionary* dictionary, char* toolboxName ) {
	Dictionary_Entry_Value* dev = Dictionary_Get( dictionary, "import" );
	
	if( !dev ) {
		dev = Dictionary_Entry_Value_NewList();
		Dictionary_Add( dictionary, "import", dev );
	}
	
	Dictionary_Entry_Value_AddElement( dev, Dictionary_Entry_Value_FromString( toolboxName ) );
}
