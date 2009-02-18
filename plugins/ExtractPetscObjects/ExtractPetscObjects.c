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
** $Id: ExtractPetscObjects.c 358 2007-06-02 06:17:30Z Mayhem $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
#include <string.h>
#include <mpi.h>
#include <petsc.h>
#include <petscmat.h>
#include <petscvec.h>

#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include <Underworld/Underworld.h>
#include <StgFEM/FrequentOutput/FrequentOutput.h>



#include "ExtractPetscObjects.h"

const Type Underworld_ExtractPetscObjects_Type = "Underworld_ExtractPetscObjects";

void _Underworld_ExtractPetscObjects_Construct( void* component, Stg_ComponentFactory* cf, void* data ) 
{
	UnderworldContext*  context;
	
	context = Stg_ComponentFactory_ConstructByName( cf, "context", UnderworldContext, True, data ); 
	
	Underworld_ExtractPetscObjects_PrintHeaderToFile( context );
	ContextEP_Append( context, AbstractContext_EP_FrequentOutput     , Underworld_ExtractPetscObjects_Dump );
}

void* _Underworld_ExtractPetscObjects_DefaultNew( Name name ) 
{
	return Codelet_New(
		Underworld_ExtractPetscObjects_Type,
		_Underworld_ExtractPetscObjects_DefaultNew,
		_Underworld_ExtractPetscObjects_Construct,
		_Codelet_Build,
		_Codelet_Initialise,
		_Codelet_Execute,
		_Codelet_Destroy,
		name );
}

Index Underworld_ExtractPetscObjects_Register( PluginsManager* pluginsManager ) 
{
	return PluginsManager_Submit( pluginsManager, Underworld_ExtractPetscObjects_Type, "0", _Underworld_ExtractPetscObjects_DefaultNew );
}


/* 
Dump the matrices and vectors to binary files
The file suffix is extracted from the xml variable <ProbDescription>

The file,
    ProbDescription-xxxx.info,
(where xxxx is the time step counter) contains the list of file names associated with 
the current problem. The list defines the names for the matrices K,G,D,C and vectors f,h.
If the name is NULL, there is no saved petsc object.

*/

void Underworld_ExtractPetscObjects_Dump( void* _context ) 
{
	UnderworldContext *context = (UnderworldContext*) _context;
	
	int		length;
	char	machine_name[256];
	int		rank;
	char	*probName, fileName[256];
	FILE	*info;
	char	kName[256];
	char	GradName[256];
	char	DivName[256];
	char	CName[256];
	char	FName[256];
	char	HName[256];
	char	SpcName[256];
	
	char			mat_name[256];
	char			vec_name[256];
	PetscViewer		mat_view_file;
	PetscViewer		vec_view_file;
	Stokes_SLE		*stokesSLE = context->stokesSLE;
	Mat A;
	Vec b;
	MPI_Comm comm;
	int step;
	Stokes_SLE_UzawaSolver *uzawaSolver;
	Bool uzawa_present;
	int ex_step, default_step;
	
	/* Init names */
	strcpy( kName,    "NULL" );
	strcpy( GradName, "NULL" );
	strcpy( DivName,  "NULL" );
	strcpy( CName,    "NULL" );
	strcpy( FName,    "NULL" );
	strcpy( HName,    "NULL" );
	
	
	step = context->timeStep;

	/* Get time step to perform extraction */
	default_step =  -6699;
	ex_step = Dictionary_Entry_Value_AsInt(
			Dictionary_GetDefault( context->dictionary, "ExtractMatricesAtStep", Dictionary_Entry_Value_FromInt(default_step) )
						);
	if( ex_step != default_step ) { /* Then we want to dump a specific step */
		if( ex_step != step ) { /* If current step does not match required step */
			return;
		}
	}	

	
	/* Write the matrices and vectors to disk */
	comm = context->communicator;
	
	// get filename from problem description
	probName = strdup( Dictionary_Entry_Value_AsString( Dictionary_Get( context->dictionary, "ProbDescription" ) ));
	printf("\n\n");
	printf("**********************************************\n" );
	printf("******     Extracting PetscObjects      ******\n" );
	printf("**********************************************\n");
	printf("  ProbDescription: %s \n", probName );
	
	MPI_Get_processor_name( machine_name, &length);
	MPI_Comm_rank( context->communicator, &rank );

	
	
	/* Write K to file */
	if( stokesSLE->kStiffMat != NULL ) {
		sprintf( kName, "k--%s-step%d-%s", probName,step, machine_name );
		sprintf( mat_name, "%s/%s", context->outputPath, kName );
		printf("  Writing kMatrix:                    %s \n",mat_name );
		
/* 		if( !stokesSLE->kStiffMat->useShellMatrix ) */
			A = stokesSLE->kStiffMat->matrix;
/* 		else */
/* 			A = stokesSLE->kStiffMat->shellMatrix->matrix; */
		PetscViewerBinaryOpen( comm, mat_name, FILE_MODE_WRITE, &mat_view_file );
		MatView( A, mat_view_file );
		PetscViewerDestroy( mat_view_file );
	}
	/* Write G to file */
	if( stokesSLE->gStiffMat != NULL ) {
		sprintf( GradName, "g--%s-step%d-%s", probName,step, machine_name );
		sprintf( mat_name, "%s/%s", context->outputPath, GradName );
		printf("  Writing Grad:                       %s \n",mat_name );
		
/* 		if( !stokesSLE->gStiffMat->useShellMatrix ) */
			A = stokesSLE->gStiffMat->matrix;
/* 		else */
/* 			A = stokesSLE->gStiffMat->shellMatrix->matrix; */
		PetscViewerBinaryOpen( comm, mat_name, FILE_MODE_WRITE, &mat_view_file );
		MatView( A, mat_view_file );
		PetscViewerDestroy( mat_view_file );
	}
	
	/* Write Div to file */	
	if( stokesSLE->dStiffMat != NULL ) {
		sprintf( DivName, "div--%s-step%d-%s", probName, step, machine_name );
		sprintf( mat_name, "%s/%s", context->outputPath, DivName );
		printf("  Writing Div:                        %s \n",mat_name );
		
/* 		if( !stokesSLE->dStiffMat->useShellMatrix ) */
			A = stokesSLE->dStiffMat->matrix;
/* 		else */
/* 			A = stokesSLE->dStiffMat->shellMatrix->matrix; */
		PetscViewerBinaryOpen( comm, mat_name, FILE_MODE_WRITE, &mat_view_file );
		MatView( A, mat_view_file );
		PetscViewerDestroy( mat_view_file );
	}
	
	/* Write C to file */
	if( stokesSLE->cStiffMat != NULL ) {
		sprintf( CName, "c--%s-step%d-%s", probName, step, machine_name );
		sprintf( mat_name, "%s/%s", context->outputPath, CName );
		printf("  Writing C:                          %s \n",mat_name );
		
/* 		if( !stokesSLE->cStiffMat->useShellMatrix ) */
			A = stokesSLE->cStiffMat->matrix;
/* 		else */
/* 			A = stokesSLE->cStiffMat->shellMatrix->matrix; */
		PetscViewerBinaryOpen( comm, mat_name, FILE_MODE_WRITE, &mat_view_file );
		MatView( A, mat_view_file );
		PetscViewerDestroy( mat_view_file );
	}
	
	
	/* Momentum rhs */
	if( stokesSLE->fForceVec != NULL ) {
		sprintf( FName, "F--%s-step%d-%s", probName, step, machine_name );
		sprintf( vec_name, "%s/%s", context->outputPath, FName );
		printf("  Writing F:                          %s \n", vec_name );
		
		b = stokesSLE->fForceVec->vector;
		PetscViewerBinaryOpen( comm, vec_name, FILE_MODE_WRITE, &vec_view_file );
		VecView( b, vec_view_file );
		PetscViewerDestroy( vec_view_file );
	}
	
	/* Continuity rhs */
	if( stokesSLE->hForceVec != NULL ) {
		sprintf( HName, "H--%s-step%d-%s", probName, step, machine_name );
		sprintf( vec_name, "%s/%s", context->outputPath, HName );
		printf("  Writing H:                          %s \n", vec_name );
		
		b = stokesSLE->hForceVec->vector;
		PetscViewerBinaryOpen( comm, vec_name, FILE_MODE_WRITE, &vec_view_file );
		VecView( b, vec_view_file );
		PetscViewerDestroy( vec_view_file );
	}
	
	/* Check if solver is Uzawa and spit out the preconditioner used */
	uzawa_present = False;
	uzawa_present = Stg_Class_IsInstance( stokesSLE->solver, Stokes_SLE_UzawaSolver_Type );
	if( uzawa_present == True ) {
		sprintf( SpcName, "Spc--%s-step%d-%s", probName, step, machine_name );
		sprintf( mat_name, "%s/%s", context->outputPath, SpcName );
		printf("  Writing Schur pc operator (UW_Q22): %s \n",mat_name );
	
		uzawaSolver = (Stokes_SLE_UzawaSolver*)stokesSLE->solver;
		A = uzawaSolver->preconditioner->matrix;
		PetscViewerBinaryOpen( comm, mat_name, FILE_MODE_WRITE, &mat_view_file );
		MatView( A, mat_view_file );
		PetscViewerDestroy( mat_view_file );
	}
	
	
	
	if( rank == 0 ) {
		
		sprintf( fileName, "%s/%s-step%d.info", context->outputPath, probName, context->timeStep );
		printf("  Info file:                          %s \n", fileName );			
		
		if (( info = fopen (fileName, "w")) == NULL)  {
			(void) fprintf(stderr, "Error: couldn't open file %s for writing. Exiting.\n", fileName);
			PetscFinalize();
			exit(1);
		}
		
		/* Write the object info to file */
		fprintf( info, "%s %s %s %s %s %s",
				kName, GradName, DivName, CName,
				FName, HName );
		
		fclose( info );
	}
	
//	b = (Vec)stokesSLE->pSolnVec->vector;
//	VecView( b, PETSC_VIEWER_STDOUT_WORLD );
	
	
}

void Underworld_ExtractPetscObjects_PrintHeaderToFile( void* context ) 
{
	StgFEM_FrequentOutput_PrintString( context, "ExtractPetscObjects" );
}

