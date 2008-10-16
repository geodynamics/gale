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
** $Id: testContext.c 964 2007-10-11 08:03:06Z SteveQuenette $
**
** Purpose: primarily in this test we're going to check the timestep
**	limiting works.
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifdef HAVE_PYTHON
	#include <Python.h>
#endif
#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include "StgFEM/Discretisation/Discretisation.h"
//#include "StgFEM/SLE/LinearAlgebra/LinearAlgebra.h"
#include "StgFEM/SLE/SystemSetup/SystemSetup.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct _Element {
	__FiniteElement_Element
};

struct _Particle {
	__IntegrationPoint
};


double testCalcDtFunc( FiniteElementContext* context, FiniteElementContext* contextArg2 ) {
	if ( context->timeStep == 1 ) {
		return 1.0;
	}
	else {
		return 1.5 * context->dt;
	}
}

int main( int argc, char* argv[] ) 
{
	/* StGermain standard bits & pieces */
	MPI_Comm			CommWorld;
	int				rank;
	int				numProcessors;
	Dictionary*			dictionary;
	Dictionary*			componentDict;
	XML_IO_Handler*			ioHandler;
	Stream*                         stream;
	
	/* context */
	FiniteElementContext*		context = NULL;
	
	/* Initialise PETSc, get world info */
	MPI_Init( &argc, &argv );
	MPI_Comm_dup( MPI_COMM_WORLD, &CommWorld );
	MPI_Comm_size( CommWorld, &numProcessors );
	MPI_Comm_rank( CommWorld, &rank );
	
	StGermain_Init( &argc, &argv );
	StgFEM_Discretisation_Init( &argc, &argv );
	//StgFEM_SLE_LinearAlgebra_Init( &argc, &argv );
	StgFEM_SLE_SystemSetup_Init( &argc, &argv );
	MPI_Barrier( CommWorld ); /* Ensures copyright info always come first in output */
	
	stream = Journal_Register( Info_Type, "test" );
	Stream_SetPrintingRank( stream, 0 );
	
	/* Create the application's dictionary */
	dictionary = Dictionary_New();

	Dictionary_Add( dictionary, "outputPath", Dictionary_Entry_Value_FromString( "output" ) );

	/* Read input */
	ioHandler = XML_IO_Handler_New();
	IO_Handler_ReadAllFromCommandLine( ioHandler, argc, argv, dictionary );

	Journal_ReadFromDictionary( dictionary );

	/* Construction phase ----------------------------------------------------------------------------------------------*/
	context = FiniteElementContext_New( "context", 0, 0, CommWorld, dictionary );
	
	Stream_SetPrintingRank( context->info, 0 );

	componentDict = Dictionary_GetDictionary( dictionary, "components" );

	if ( componentDict == NULL ) {
		componentDict = Dictionary_New();
	}
	context->CF = Stg_ComponentFactory_New( dictionary, componentDict, context->register_Register );
	
	/* This is where we'd normally construct components if it was real main.
		instead, we'll just set the dt function so we can test it */
	EP_AppendClassHook( context->calcDtEP, testCalcDtFunc, context );
	
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

	context->maxTimeSteps = 10;

	Journal_Printf( stream, "Running with no timestep braking, using  "
		"dt that increases 50%% each step:\n" );

	Stg_Component_Execute( context, 0 /* dummy */, False );
	context->currentTime=0;
	context->dt = 0;

	Journal_Printf( stream, "\nTurning on timestep braking, at default "
		"level, running again:\n" );
	context->limitTimeStepIncreaseRate = True;
	Stg_Component_Execute( context, 0 /* dummy */, True );
	context->currentTime=0;
	context->dt = 0;

	Journal_Printf( stream, "\nTurning on timestep braking, at 80%% "
		"level, running again - expect same as original:\n" );
	context->maxTimeStepIncreasePercentage = 80;
	Stg_Component_Execute( context, 0 /* dummy */, True );

	/* Destruct phase ---------------------------------------------------------------------------------------------------*/
	Stg_Component_Destroy( context, 0 /* dummy */, False );
	Stg_Class_Delete( context );
	Stg_Class_Delete( dictionary );

	StgFEM_SLE_SystemSetup_Finalise();
	//StgFEM_SLE_LinearAlgebra_Finalise();
	StgFEM_Discretisation_Finalise();
	StGermain_Finalise();
		
	/* Close off MPI */
	MPI_Finalize();
	
	return 0; /* success */
}
