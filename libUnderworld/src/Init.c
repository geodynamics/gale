/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
** Copyright (c) 2005, Monash Cluster Computing 
** All rights reserved.
** Redistribution and use in source and binary forms, with or without modification,
** are permitted provided that the following conditions are met:
**
** 		* Redistributions of source code must retain the above copyright notice, 
** 			this list of conditions and the following disclaimer.
** 		* Redistributions in binary form must reproduce the above copyright 
**			notice, this list of conditions and the following disclaimer in the 
**			documentation and/or other materials provided with the distribution.
** 		* Neither the name of the Monash University nor the names of its contributors 
**			may be used to endorse or promote products derived from this software 
**			without specific prior written permission.
**
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
** THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR 
** PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS 
** BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
** CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
** SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) 
** HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT 
** LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT 
** OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
**
**
** Contact:
*%		Louis Moresi - Louis.Moresi@sci.monash.edu.au
*%
** Contributors:
*+		Robert Turnbull
*+		Vincent Lemiale
*+		Louis Moresi
*+		David May
*+		David Stegman
*+		Mirko Velic
*+		Patrick Sunter
*+		Julian Giordani
*+
** $Id: Init.c 610 2007-10-11 08:09:29Z SteveQuenette $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>

#include "Underworld/Rheology/Rheology.h"
#include "Underworld/Utils/Utils.h"

#include "Init.h"

#include <stdio.h>
#include <signal.h>
#include <string.h>

/** Initialises this package, then any init for this package
such as streams etc */

Stream* Underworld_Debug = NULL;
Stream* Underworld_Info = NULL;
Stream* Underworld_Error = NULL;

Bool Underworld_Init( int* argc, char** argv[] ) {
	/* This init function tells StGermain of all the component types, etc this module contributes. Because it can be linked at compile
	   time or linked in by a toolbox at runtime, we need to make sure it isn't run twice (compiled in and loaded through a toolbox.*/
	if( !ToolboxesManager_IsInitialised( stgToolboxesManager, "Underworld" ) ) {
		char* argString;
		int arg_I;
		int tmp;
		Bool useSignalHandler = True;
		char* directory;

		for ( arg_I = 0; argc && arg_I < *argc; arg_I++ ) {
			argString = (*argv)[arg_I];
			/* Leverage from PETSC's no signal flag */
			if ( strcmp( argString, "-no_signal_handler" ) == 0 ) {
				useSignalHandler = False;
			}
		}

		if ( useSignalHandler ) {
			signal( SIGSEGV, Underworld_SignalHandler );
			signal( SIGTERM, Underworld_SignalHandler );
		}

		Underworld_Rheology_Init( argc, argv );
		Underworld_Utils_Init( argc, argv );
	
		Journal_Printf( Journal_Register( DebugStream_Type, "Context" ), "In: %s\n", __func__ ); /* DO NOT CHANGE OR REMOVE */
		tmp = Stream_GetPrintingRank( Journal_Register( InfoStream_Type, "Context" ) );
		Stream_SetPrintingRank( Journal_Register( InfoStream_Type, "Context" ), 0 );
		Journal_Printf( /* DO NOT CHANGE OR REMOVE */
			Journal_Register( InfoStream_Type, "Context" ), 
			"Underworld (Geodynamics framework) revision %s. Copyright (C) 2005 Monash Cluster Computing.\n", VERSION );
		Stream_Flush( Journal_Register( InfoStream_Type, "Context" ) );
		Stream_SetPrintingRank( Journal_Register( InfoStream_Type, "Context" ), tmp );

		/* Create Streams */
		Underworld_Debug  = Journal_Register( Debug_Type, "Context" );
		Underworld_Info   = Journal_Register( Info_Type,  "Context" );
		Underworld_Error  = Journal_Register( Error_Type, "Context" );
	
		/* Add the Underworld path to the global xml path dictionary */
		directory = Memory_Alloc_Array( char, 400, "xmlDirectory" ) ;
		sprintf(directory, "%s%s", LIB_DIR, "/StGermain" );
		XML_IO_Handler_AddDirectory( "Underworld", directory );
		Memory_Free(directory);

		/* Add the plugin path to the global plugin list */
		#ifdef GLUCIFER_LIBDIR
			ModulesManager_AddDirectory( "gLucifer", GLUCIFER_LIBDIR );
		#endif
	
		ModulesManager_AddDirectory( "Underworld", LIB_DIR );

		ToolboxesManager_SetInitialised( stgToolboxesManager, "Underworld" );
		return True;
	}
	return False;
}


void Underworld_SignalHandler( int signal ) {
	fprintf( stderr, 
			"\n\n=====================================================================================\n"
			"Error running Underworld (revision %s) - Signal %d ",
			VERSION, signal );

	switch ( signal ) {
		case SIGSEGV:
			fprintf( stderr, 
					"'SIGSEGV' (Segmentation Fault).\n" 
					"This is probably caused by an illegal access of memory.\n"
					"We recommend running the code in a debugger to work out where the problem is (e.g. 'gdb')\n"
					"and also to contact the developers.\n" );
			break;
		case SIGTERM:
			fprintf( stderr, 
					"'SIGTERM' (Termination Request).\n" 
					"This is caused by an external call to terminate the code.\n"
					"This could have happened by a queueing system (e.g. if the code has run longer than allowed),\n"
					"the code might have been killed on another processor or it may have been killed by the user.\n" );
			break;
	}
	exit(EXIT_FAILURE);
}
