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
** $Id: main.c 645 2008-02-12 01:18:41Z DavidMay $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#ifdef HAVE_PYTHON
	#include <Python.h>
#endif
#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

const Type Underworld_Type = "Underworld";

int main( int argc, char* argv[] ) {
	/* StGermain standard bits & pieces */
	MPI_Comm			CommWorld;
	int				rank;
	int				numProcessors;
	Dictionary*			dictionary;
	XML_IO_Handler*			ioHandler;
	Stream*                         stream;

	/* Initialise PETSc, get world info */
	MPI_Init( &argc, &argv );
	MPI_Comm_dup( MPI_COMM_WORLD, &CommWorld );
	MPI_Comm_size( CommWorld, &numProcessors );
	MPI_Comm_rank( CommWorld, &rank );
	StGermain_Init( &argc, &argv );
	stream = Journal_Register( Info_Type, Underworld_Type );
	#ifdef HAVE_PYTHON
		Py_Initialize();
	#endif	
	MPI_Barrier( CommWorld ); /* Ensures copyright info always come first in output */
	
	/* Create the application's dictionary & read input */
	dictionary = Dictionary_New();
	ioHandler = XML_IO_Handler_New();
	IO_Handler_ReadAllFromCommandLine( ioHandler, argc, argv, dictionary );
	Journal_ReadFromDictionary( dictionary );
	stgImportToolbox( dictionary, "Underworld" );
	
	stgMainLoop( dictionary, CommWorld );
	
	Stg_Class_Delete( dictionary );
	
	/* Close off everything */
	#ifdef HAVE_PYTHON
		Py_Finalize();
	#endif
	StGermain_Finalise();
	Journal_RPrintf( stream, "Finalised: StGermain Framework.\n");
	MPI_Finalize();
	
	return 0; /* success */
}
