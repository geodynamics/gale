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
** $Id: MatrixSolver.c 656 2006-10-18 06:45:50Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgFEM/Discretisation/Discretisation.h>
#include "SLE/LinearAlgebra/LinearAlgebra.h"

#include <petsc.h>
#include <petscmat.h>
#include "petscksp.h"
#include "petscpc.h"
#include "petscmg.h"
#include <StGermain/compatibility/petsccompat.h>

#include "ErrorChecking.h"

MatrixSolver* MatrixSolver_Build( MPI_Comm comm, void* matrix ) {
	SLES           solverCtx;
	PetscErrorCode errorFlag;
	
	errorFlag = SLESCreate( comm, &solverCtx );
	CheckPETScError( errorFlag );

	errorFlag = SLESSetFromOptions( (SLES)solverCtx );
	CheckPETScError( errorFlag );

	return (MatrixSolver*)solverCtx;
}

/*
This is the simplist way to set up a PETSc solver. Many options are available
to specialise and tune the solver setup.
*/
void MatrixSolver_Setup( void* solverCtx, void* matrix ) {
	PetscErrorCode errorFlag;
	
	errorFlag = SLESSetOperators( (SLES)solverCtx, matrix, matrix, DIFFERENT_NONZERO_PATTERN );
	CheckPETScError( errorFlag );
}


void MatrixSolver_SetupKSP( void* solver ) {
	PetscErrorCode	ec;
	KSP		ksp;

	SLESGetKSP( (SLES)solver, &ksp );
	ec = KSPSetUp( ksp );	CheckPETScError( ec );
}


void MatrixSolver_Setup_SameNonZeroPattern( void* solverCtx, void* matrix ) {
	PetscErrorCode errorFlag;
	
	errorFlag = SLESSetOperators( (SLES)solverCtx, matrix, matrix, SAME_NONZERO_PATTERN );
	CheckPETScError( errorFlag );

	errorFlag = SLESSetFromOptions( (SLES)solverCtx );
	CheckPETScError( errorFlag );
}

void MatrixSolver_SetPC_Type( void* matSolver, char* pcType ) {
	PetscErrorCode	ef;
	PC		pc;
	
	ef = SLESGetPC( (SLES)matSolver, &pc );
	ef = PCSetType( pc, pcType );
	
	CheckPETScError( ef );
}


void MatrixSolver_SetKSP_Type( void* matSolver, char* kspType ) {
	PetscErrorCode	ef;
	KSP		ksp;
	
	SLESGetKSP( (SLES)matSolver, &ksp );
	ef = KSPSetType( ksp, kspType );
	
	CheckPETScError( ef );
}


/* 	The following routines should not clobber each other (the PETSc code does nothing at all where PETSC_DEFAULT
 	is given as an argument - it does not restore a default value each time)
*/
 
void MatrixSolver_SetMaxIts( void* solverCtx, Iteration_Index maxIts ) {
	KSP            ksp;
	PetscErrorCode errorFlag;
 
	SLESGetKSP( (SLES)solverCtx, &ksp );

	errorFlag = KSPSetTolerances( ksp, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, maxIts );
	CheckPETScError( errorFlag );
}


void MatrixSolver_SetRelativeTolerance( void* solverCtx, PetscReal relativeTolerance ) {
	KSP            ksp;
	PetscErrorCode errorFlag;
 
	SLESGetKSP( (SLES)solverCtx, &ksp );

	errorFlag = KSPSetTolerances( ksp, relativeTolerance, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT);
	CheckPETScError( errorFlag );
}

void MatrixSolver_SetAbsoluteTolerance( void* solverCtx, PetscReal absoluteTolerance ) {
	KSP            ksp;
	PetscErrorCode errorFlag;
 
	SLESGetKSP( (SLES)solverCtx, &ksp );

	errorFlag = KSPSetTolerances( ksp, PETSC_DEFAULT, absoluteTolerance, PETSC_DEFAULT, PETSC_DEFAULT);
	CheckPETScError( errorFlag );
}



Iteration_Index MatrixSolver_Solve( MatrixSolver* matSolver, Vector* solnVec, Vector* rhsVec ) {
	PetscErrorCode errorFlag;
	int iterationCount;

	#ifdef DUMP_FOR_PETSCGROUP
		MPI_Comm comm;
		PetscViewer viewer;
		Mat mat;
		char filename[256];
		static int counter = 0;
	#endif
	
	#ifdef DUMP_FOR_PETSCGROUP
		errorFlag = PetscObjectGetComm((PetscObject) matSolver, &comm); CHKERRQ(errorFlag);
		sprintf( filename, "petsc-%u.bin", counter );
		errorFlag = PetscViewerBinaryOpen(comm, filename, PETSC_FILE_CREATE, &viewer); CHKERRQ(errorFlag);
		errorFlag = KSPGetOperators((KSP)matSolver, &mat, PETSC_NULL, PETSC_NULL); CHKERRQ(errorFlag);
		errorFlag = MatView( mat, viewer ); CHKERRQ(errorFlag);
		errorFlag = VecView( (Vec)rhsVec, viewer ); CHKERRQ(errorFlag);
		errorFlag = PetscViewerDestroy(viewer); CHKERRQ(errorFlag);
		counter += 1;
	#endif

	errorFlag = SLESSolve( (SLES)matSolver, (Vec)rhsVec, (Vec)solnVec, &iterationCount );
	CheckPETScError( errorFlag );

	return (Iteration_Index) iterationCount;
}


Iteration_Index MatrixSolver_StatSolve( MatrixSolver* solver, 
								Vector* sol, Vector* rhs, 
								unsigned nReps )
{
	PetscInt		nIts, maxIts;
	PetscTruth	guess;
	Vec			tmpSol, tmpRes;
	double*		itTimes;
	double*		netTimes;
	double*		rNorms;
	unsigned		it_i;
	Mat			aMat;
	KSP			ksp;
	PetscErrorCode	ec;
	
	/* Backwards compatibility. */
	SLESGetKSP( (SLES)solver, &ksp );
	
	/* Backup PETSc state. */
	ec = KSPGetInitialGuessNonzero( ksp, &guess );						CheckPETScError( ec );
	ec = KSPGetTolerances( ksp, PETSC_NULL, PETSC_NULL, PETSC_NULL, &maxIts );	CheckPETScError( ec );
	
	/* Make a copy of the solution vector and a residual. */
	VecDuplicate( (Vec)sol, &tmpRes );
	VecDuplicate( (Vec)sol, &tmpSol );
	VecCopy( (Vec)sol, tmpSol );
	
	/* Extract the KSP's operators. */
	KSPGetOperators( ksp, &aMat, PETSC_NULL, PETSC_NULL );	CheckPETScError( ec );
	
	/* We need to find out how many iterations this solve will use. */
	ec = SLESSolve( (SLES)solver, (Vec)rhs, tmpSol, &nIts );	CheckPETScError( ec );
	if( nReps < nIts ) {
		nIts = nReps;
	}
	
	/* Allocate space for the results. */
	itTimes = Memory_Alloc_Array( double, nIts, "MatrixSolver" );
	netTimes = Memory_Alloc_Array( double, nIts, "MatrixSolver" );
	rNorms = Memory_Alloc_Array( double, nIts, "MatrixSolver" );
	
	/* Calculate values for each level of iterations. */
	for( it_i = 0; it_i < nIts; it_i++ ) {
		/* Cap the number of iterations. */
		ec = KSPSetTolerances( ksp, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, it_i + 1 );	CheckPETScError( ec );
		
		/* Clear the accumulative values. */
		netTimes[it_i] = 0.0;
		rNorms[it_i] = 0.0;
		nReps = 0;
		
		/* Perform the solve the number of times specified. */
		do {
			double	sTime = MPI_Wtime();
			PetscReal	resNorm;
			unsigned	tmpIts;
			
			/* Revert the temporary solution back to the original solution. */
			VecCopy( (Vec)sol, tmpSol );
			
			/* Solve and store results. */
			ec = SLESSolve( (SLES)solver, (Vec)rhs, tmpSol, (PetscInt*)(&tmpIts) );	CheckPETScError( ec );
			netTimes[it_i] += MPI_Wtime() - sTime;
			Vector_ScaleContents( (Vector*) tmpSol, -1.0 );						CheckPETScError( ec );
			ec = MatMultAdd( aMat, tmpSol, (Vec)rhs, tmpRes );		CheckPETScError( ec );
			ec = VecNorm( tmpRes, NORM_2, &resNorm );				CheckPETScError( ec );
			rNorms[it_i] += resNorm;
			nReps++;
		}
		while( netTimes[it_i] < 30.0 );
		
		/* Average the results. */
		itTimes[it_i] = netTimes[it_i] / (double)nReps;
		rNorms[it_i] /= (double)nReps;
	}
	
	/* Kill the temporary vectors. */
	VecDestroy( tmpSol );
	VecDestroy( tmpRes );
	
	/* Restore original settings. */
	ec = KSPSetInitialGuessNonzero( ksp, guess );								CheckPETScError( ec );
	ec = KSPSetTolerances( ksp, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, maxIts );	CheckPETScError( ec );
	
	/* Output the results. */
	printf( "*** Stat: Solved up to %d iterations using %d solve repetitions:\n", nIts, nReps );
	for( it_i = 0; it_i < nIts; it_i++ ) {
		printf( "*** Stat: \tResults for %d iterations:\n", it_i + 1 );
		printf( "*** Stat: \t\tAverage residual norm: %g\n", rNorms[it_i] );
		printf( "*** Stat: \t\tAverage iteration time: %g\n", itTimes[it_i] );
		printf( "*** Stat: \t\tTotal time: %g\n", netTimes[it_i] );
	}
	
	/* Free storage arrays. */
	FreeArray( rNorms );
	FreeArray( itTimes );
	FreeArray( netTimes );
	
	return nIts;
}


void MatrixSolver_Destroy( MatrixSolver* matSolver ) {
	PetscErrorCode errorFlag;
	
	errorFlag = SLESDestroy( (SLES)matSolver );
	CheckPETScError( errorFlag );
}


double MatrixSolver_GetResidualNorm( MatrixSolver* solver ) {
	PetscErrorCode	ec;
	PetscReal		rNorm;
	KSP			ksp;
	
	SLESGetKSP( (SLES)solver, &ksp );
	ec = KSPGetResidualNorm( ksp, &rNorm );	CheckPETScError( ec );
	
	return rNorm;
}


void MatrixSolver_HaveInitialGuess( MatrixSolver* solver ) {
   KSP	ksp;

   SLESGetKSP( (SLES)solver, &ksp );
   KSPSetInitialGuessNonzero( ksp, PETSC_TRUE );
}


void MatrixSolver_PC_GetBJacobiSubBlocks( MatrixSolver* solver, MatrixSolver*** blocks, unsigned* nBlocks ) {
	PetscErrorCode	ec;
	KSP		ksp;
	PC		pc;

	SLESGetKSP( (SLES)solver, &ksp );
	ec = KSPGetPC( ksp, &pc );						CheckPETScError( ec );
	ec = PCBJacobiGetSubKSP( pc, (PetscInt*)nBlocks, PETSC_NULL, (KSP**)blocks );	CheckPETScError( ec );
}


void MatrixSolver_PC_SetSORIts( MatrixSolver* solver, unsigned nIts ) {
	PetscErrorCode	ec;
	KSP		ksp;
	PC		pc;

	SLESGetKSP( (SLES)solver, &ksp );
	ec = KSPGetPC( ksp, &pc );		CheckPETScError( ec );
	ec = PCSORSetIterations( pc, 1, nIts );	CheckPETScError( ec );
}


/*
** Multi-grid functions.
*/

void MatrixSolver_MG_SetLevels( MatrixSolver* solver, unsigned nLevels ) {
	PetscErrorCode	ec;
	PC		pc;
	KSP		ksp;
	
	SLESGetKSP( (SLES)solver, &ksp );
	ec = KSPGetPC( ksp, &pc );			CheckPETScError( ec );
	ec = PCMGSetLevels( pc, nLevels, NULL );	CheckPETScError( ec );
	ec = PCMGSetType( pc, PC_MG_MULTIPLICATIVE );	CheckPETScError( ec );
}


void MatrixSolver_MG_SetCycles( MatrixSolver* solver, unsigned level, unsigned nCycles ) {
	PetscErrorCode	ec;
	PC		pc;
	KSP		ksp;
	
	SLESGetKSP( (SLES)solver, &ksp );
	ec = KSPGetPC( ksp, &pc );				CheckPETScError( ec );
	ec = PCMGSetCyclesOnLevel( pc, level, nCycles );	CheckPETScError( ec );
}


void MatrixSolver_MG_SetLevelMatrix( MatrixSolver* solver, unsigned level, Matrix* mat ) {
	PetscErrorCode	ec;
	PC		pc;
	KSP		ksp;

	SLESGetKSP( (SLES)solver, &ksp );
	ec = KSPGetPC( ksp, &pc );
	ec = PCMGSetResidual( pc, level, PCMGDefaultResidual, (Mat)mat );
}


void MatrixSolver_MG_SetSmoothUpIts( MatrixSolver* solver, unsigned nUpIts ) {
	PetscErrorCode	ec;
	PC		pc;
	KSP		ksp;

	SLESGetKSP( (SLES)solver, &ksp );
	ec = KSPGetPC( ksp, &pc );			CheckPETScError( ec );
	ec = PCMGSetNumberSmoothUp( pc, (int)nUpIts );	CheckPETScError( ec );
}


void MatrixSolver_MG_SetSmoothDownIts( MatrixSolver* solver, unsigned nDownIts ) {
	PetscErrorCode	ec;
	PC		pc;
	KSP		ksp;

	SLESGetKSP( (SLES)solver, &ksp );
	ec = KSPGetPC( ksp, &pc );				CheckPETScError( ec );
	ec = PCMGSetNumberSmoothDown( pc, (int)nDownIts );	CheckPETScError( ec );
}


void MatrixSolver_MG_GetDownSmoother( MatrixSolver* solver, unsigned level, MatrixSolver** smoother ) {
	PetscErrorCode	ec;
	PC		pc;
	KSP		ksp;
	
	SLESGetKSP( (SLES)solver, &ksp );
	ec = KSPGetPC( ksp, &pc );					CheckPETScError( ec );
	ec = PCMGGetSmootherDown( pc, (int)level, (KSP*)smoother );	CheckPETScError( ec );
}


void MatrixSolver_MG_GetUpSmoother( MatrixSolver* solver, unsigned level, MatrixSolver** smoother ) {
	PetscErrorCode	ec;
	PC		pc;
	KSP		ksp;
	
	SLESGetKSP( (SLES)solver, &ksp );
	ec = KSPGetPC( ksp, &pc );					CheckPETScError( ec );
	ec = PCMGGetSmootherUp( pc, (int)level, (KSP*)smoother );	CheckPETScError( ec );
}


void MatrixSolver_MG_SetRestrictionOp( MatrixSolver* solver, unsigned level, Matrix* op ) {
	PetscErrorCode	ec;
	PC		pc;
	KSP		ksp;
	
	SLESGetKSP( (SLES)solver, &ksp );
	ec = KSPGetPC( ksp, &pc );				CheckPETScError( ec );
	ec = PCMGSetRestriction( pc, (int)level, (Mat)op );	CheckPETScError( ec );
}


void MatrixSolver_MG_SetInterpolationOp( MatrixSolver* solver, unsigned level, Matrix* op ) {
	PetscErrorCode	ec;
	PC		pc;
	KSP		ksp;
	
	SLESGetKSP( (SLES)solver, &ksp );
	ec = KSPGetPC( ksp, &pc );				CheckPETScError( ec );
	ec = PCMGSetInterpolate( pc, (int)level, (Mat)op );	CheckPETScError( ec );
}


void MatrixSolver_MG_SetWorkSpace( MatrixSolver* solver, unsigned level, Vector* rhsVec, Vector* xVec, Vector* rVec ) {
	PetscErrorCode	ec;
	PC		pc;
	KSP		ksp;
	
	SLESGetKSP( (SLES)solver, &ksp );
	ec = KSPGetPC( (KSP)solver, &pc );				CheckPETScError( ec );
	if( rhsVec ) {
		ec = PCMGSetRhs( pc, (PetscInt)level, (Vec)rhsVec );	CheckPETScError( ec );
	}
	if( xVec ) {
		ec = PCMGSetX( pc, (PetscInt)level, (Vec)xVec );	CheckPETScError( ec );
	}
	if( rVec ) {
		ec = PCMGSetR( pc, (PetscInt)level, (Vec)rVec );	CheckPETScError( ec );
	}
}


void MatrixSolver_MG_ForceMonitor( MatrixSolver* solver ) {
	PetscErrorCode	ec;
	KSP		ksp;

	SLESGetKSP( (SLES)solver, &ksp );
	ec = KSPSetMonitor( ksp, KSPDefaultMonitor, PETSC_NULL, PETSC_NULL );
	CheckPETScError( ec );
}
