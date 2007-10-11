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
*%		Cecile Duboz - Cecile.Duboz@sci.monash.edu.au
*%
** Contributors:
*+		Cecile Duboz
*+		Robert Turnbull
*+		Alan Lo
*+		Louis Moresi
*+		David Stegman
*+		David May
*+		Stevan Quenette
*+		Patrick Sunter
*+		Greg Watson
*+
** $Id: main.c 740 2007-10-11 08:05:31Z SteveQuenette $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/



#ifdef HAVE_PYTHON
#include <Python.h>
#endif

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <glucifer/glucifer.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

int main( int argc, char* argv[] ) {
	MPI_Comm               CommWorld;
	int                    rank;
	int                    numProcessors;
	Dictionary*            dictionary;
	XML_IO_Handler*        ioHandler;
	DomainContext* context         = NULL;
	
	/* Initialise PETSc, get world info */
	MPI_Init( &argc, &argv );
	MPI_Comm_dup( MPI_COMM_WORLD, &CommWorld );
	MPI_Comm_size( CommWorld, &numProcessors );
	MPI_Comm_rank( CommWorld, &rank );
	
	StGermain_Init( &argc, &argv );
	StgDomain_Init( &argc, &argv );
	glucifer_Init();
	#ifdef HAVE_PYTHON
	Py_Initialize();
	#endif

	MPI_Barrier( CommWorld ); /* Ensures copyright info always come first in output */
	
	/* Create the application's dictionary */
	dictionary = Dictionary_New();

	/* Read input */
	ioHandler = XML_IO_Handler_New();
	IO_Handler_ReadAllFromCommandLine( ioHandler, argc, argv, dictionary );
	Journal_ReadFromDictionary( dictionary );

	/* Construction phase -----------------------------------------------------------------------------------------------*/
	context = DomainContext_New( "context", 0, 0, CommWorld, dictionary );
	Stg_Component_Construct( context, 0 /* dummy */, &context, True	);

	if( rank == 0 ) 
		Context_PrintConcise( context, context->verbose );

	if ( True == Dictionary_GetBool_WithDefault( dictionary, "showJournalStatus", False ) ) {
		Journal_PrintConcise();	
	}			

	/* Building phase ---------------------------------------------------------------------------------------------------*/
	Stg_Component_Build( context, 0 /* dummy */, False );
	
	/* Initialisaton phase ----------------------------------------------------------------------------------------------*/
	Stg_Component_Initialise( context, 0 /* dummy */, False );
	
	/* Run (Solve) phase ------------------------------------------------------------------------------------------------*/
	AbstractContext_Dump( context );
	
	/* Destruct phase ---------------------------------------------------------------------------------------------------*/
	Stg_Component_Destroy( context, 0 /* dummy */, False );
	Stg_Class_Delete( context );
	Stg_Class_Delete( dictionary );

	#ifdef HAVE_PYTHON
	Py_Finalize();
	#endif

	lucBase_Finalise();
	StgDomain_Finalise();
	StGermain_Finalise();
		
	/* Close off MPI */
	MPI_Finalize();
	
	return 0; /* success */
}
