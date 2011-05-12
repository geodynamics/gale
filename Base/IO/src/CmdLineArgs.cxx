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
** $Id: Dictionary.c 4096 2007-05-16 00:54:10Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <string.h>

#include "Base/Foundation/Foundation.h"


void stgRemoveCmdLineArg( int* argc, char** argv[], int index ) {
	if( index > 0 && index < *argc ) {
		int i;
		char* tmpPtr = (*argv)[index];

		for( i = index; i < *argc - 1; i++ ) {
			(*argv)[i] = (*argv)[i+1];
		}

		(*argv)[*argc-1] = tmpPtr;
		*argc -= 1;
	}
}


char* stgParseHelpCmdLineArg( int* argc, char** argv[] ) {
	int                   arg_I;

	/* Loop over all the arguments from command line and reads all arguments of form "--help topic" */
	for( arg_I = 1 ; arg_I < *argc; arg_I++ ) {
		char*                   valueString = 0;
		char*                   argumentString = (*argv)[arg_I];
		Name             preceedingString = "--help";
		unsigned int            preceedingStringLength = strlen( preceedingString );

		/* Check is string has preceeding string is "--help" if not then continue in loop */
		if( strncmp( preceedingString, argumentString , preceedingStringLength ) != 0 ) {
			continue;
		}
		if( strlen( argumentString ) != strlen( preceedingString ) ) {
			continue;
		}
		if( arg_I >= (*argc - 1) ) {
			continue;
		}

		valueString = StG_Strdup( (*argv)[arg_I+1] );
		stgRemoveCmdLineArg( argc, argv, arg_I ); /* first argument: --help */
		stgRemoveCmdLineArg( argc, argv, arg_I ); /* second argument: topic */
		return valueString;		
	}

	return 0;
}

char* stgParseListCmdLineArg( int* argc, char** argv[] ) {
	int                   arg_I;

	/* Loop over all the arguments from command line and reads all arguments of form "--help topic" */
	for( arg_I = 1 ; arg_I < *argc; arg_I++ ) {
		char*                   valueString = 0;
		char*                   argumentString = (*argv)[arg_I];
		Name             preceedingString = "--list";
		unsigned int            preceedingStringLength = strlen( preceedingString );

		/* Check is string has preceeding string is "--list" if not then continue in loop */
		if( strncmp( preceedingString, argumentString , preceedingStringLength ) != 0 ) {
			continue;
		}
		if( strlen( argumentString ) != strlen( preceedingString ) ) {
			continue;
		}
		if( arg_I >= (*argc - 1) ) { /* "--list" is the last commandline argument */
			valueString = StG_Strdup( "" );
			stgRemoveCmdLineArg( argc, argv, arg_I ); /* first argument: --list */
			return valueString;
		}

		valueString = StG_Strdup( (*argv)[arg_I+1] );
		stgRemoveCmdLineArg( argc, argv, arg_I ); /* first argument: --list */
		stgRemoveCmdLineArg( argc, argv, arg_I ); /* second argument: topic */
		return valueString;		
	}

	return 0;
}

char* stgParseListAllCmdLineArg( int* argc, char** argv[] ) {
	int                   arg_I;

	/* Loop over all the arguments from command line and reads all arguments of form "--help topic" */
	for( arg_I = 1 ; arg_I < *argc; arg_I++ ) {
		char*                   valueString = 0;
		char*                   argumentString = (*argv)[arg_I];
		Name             preceedingString = "--list-all";
		unsigned int            preceedingStringLength = strlen( preceedingString );

		/* Check is string has preceeding string is "--list" if not then continue in loop */
		if( strncmp( preceedingString, argumentString , preceedingStringLength ) != 0 ) {
			continue;
		}
		if( strlen( argumentString ) != strlen( preceedingString ) ) {
			continue;
		}
		if( arg_I >= (*argc - 1) ) { /* "--list" is the last commandline argument */
			valueString = StG_Strdup( "" );
			stgRemoveCmdLineArg( argc, argv, arg_I ); /* first argument: --list */
			return valueString;
		}

		valueString = StG_Strdup( (*argv)[arg_I+1] );
		stgRemoveCmdLineArg( argc, argv, arg_I ); /* first argument: --list */
		stgRemoveCmdLineArg( argc, argv, arg_I ); /* second argument: topic */
		return valueString;		
	}

	return 0;
}



