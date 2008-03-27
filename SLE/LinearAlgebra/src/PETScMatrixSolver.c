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
** $Id: PETScMatrixSolver.c 3584 2006-05-16 11:11:07Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include "Discretisation/Discretisation.h"
#include "LinearAlgebra.h"


/* Textual name of this class */
const Type PETScMatrixSolver_Type = "PETScMatrixSolver";


/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/

PETScMatrixSolver* PETScMatrixSolver_New( Name name ) {
	return _PETScMatrixSolver_New( sizeof(PETScMatrixSolver), 
				       PETScMatrixSolver_Type, 
				       _PETScMatrixSolver_Delete, 
				       _PETScMatrixSolver_Print, 
				       NULL, 
				       (void* (*)(Name))PETScMatrixSolver_New, 
				       _PETScMatrixSolver_Construct, 
				       _PETScMatrixSolver_Build, 
				       _PETScMatrixSolver_Initialise, 
				       _PETScMatrixSolver_Execute, 
				       _PETScMatrixSolver_Destroy, 
				       name, 
				       NON_GLOBAL, 
				       PETScMatrixSolver_SetComm, 
				       PETScMatrixSolver_SetMatrix, 
				       PETScMatrixSolver_SetMaxIterations, 
				       PETScMatrixSolver_SetRelativeTolerance, 
				       PETScMatrixSolver_SetAbsoluteTolerance, 
				       PETScMatrixSolver_SetUseInitialSolution, 
				       PETScMatrixSolver_Solve, 
				       PETScMatrixSolver_Setup, 
				       PETScMatrixSolver_GetSolveStatus, 
				       PETScMatrixSolver_GetIterations, 
				       PETScMatrixSolver_GetMaxIterations, 
				       PETScMatrixSolver_GetResidualNorm );
}

PETScMatrixSolver* _PETScMatrixSolver_New( PETSCMATRIXSOLVER_DEFARGS ) {
	PETScMatrixSolver*	self;

	/* Allocate memory */
	assert( sizeOfSelf >= sizeof(PETScMatrixSolver) );
	self = (PETScMatrixSolver*)_MatrixSolver_New( MATRIXSOLVER_PASSARGS );

	/* Virtual info */

	/* PETScMatrixSolver info */
	_PETScMatrixSolver_Init( self );

	return self;
}

void _PETScMatrixSolver_Init( PETScMatrixSolver* self ) {
	assert( self && Stg_CheckType( self, PETScMatrixSolver ) );

	self->ksp = PETSC_NULL;
	PETScMatrixSolver_SetComm( self, MPI_COMM_WORLD );
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _PETScMatrixSolver_Delete( void* matrixSolver ) {
	PETScMatrixSolver*	self = (PETScMatrixSolver*)matrixSolver;

	assert( self && Stg_CheckType( self, PETScMatrixSolver ) );

	KSPDestroy( self->ksp );

	/* Delete the parent. */
	_MatrixSolver_Delete( self );
}

void _PETScMatrixSolver_Print( void* matrixSolver, Stream* stream ) {
	PETScMatrixSolver*	self = (PETScMatrixSolver*)matrixSolver;
	
	/* Set the Journal for printing informations */
	Stream* matrixSolverStream;
	matrixSolverStream = Journal_Register( InfoStream_Type, "PETScMatrixSolverStream" );

	assert( self && Stg_CheckType( self, PETScMatrixSolver ) );

	/* Print parent */
	Journal_Printf( stream, "PETScMatrixSolver (ptr): (%p)\n", self );
	_MatrixSolver_Print( self, stream );
}

void _PETScMatrixSolver_Construct( void* matrixSolver, Stg_ComponentFactory* cf, void* data ) {
	PETScMatrixSolver*		self = (PETScMatrixSolver*)matrixSolver;

	assert( self && Stg_CheckType( self, PETScMatrixSolver ) );
	assert( cf );

	_MatrixSolver_Construct( self, cf, data );
}

void _PETScMatrixSolver_Build( void* matrixSolver, void* data ) {
}

void _PETScMatrixSolver_Initialise( void* matrixSolver, void* data ) {
}

void _PETScMatrixSolver_Execute( void* matrixSolver, void* data ) {
}

void _PETScMatrixSolver_Destroy( void* matrixSolver, void* data ) {
}

void PETScMatrixSolver_SetComm( void* matrixSolver, MPI_Comm comm ) {
	PETScMatrixSolver*	self = (PETScMatrixSolver*)matrixSolver;
	PetscErrorCode		ec;

	assert( self && Stg_CheckType( self, PETScMatrixSolver ) );

	_MatrixSolver_SetComm( self, comm );

	if( self->ksp != PETSC_NULL )
		KSPDestroy( self->ksp );
	ec = KSPCreate( self->comm, &self->ksp );
	CheckPETScError( ec );
}

void PETScMatrixSolver_SetMatrix( void* matrixSolver, void* _matrix ) {
	PETScMatrixSolver*	self = (PETScMatrixSolver*)matrixSolver;
	PETScMatrix*		matrix = (PETScMatrix*)_matrix;
	PetscErrorCode		ec;

	assert( self && Stg_CheckType( self, PETScMatrixSolver ) );
	assert( matrix && Stg_CheckType( matrix, PETScMatrix ) );

	_MatrixSolver_SetMatrix( self, matrix );

	ec = KSPSetOperators( self->ksp, matrix->petscMat, matrix->petscMat, DIFFERENT_NONZERO_PATTERN );
	CheckPETScError( ec );
	if( !self->optionsReady ) {
		ec = KSPSetFromOptions( self->ksp );
		CheckPETScError( ec );
		self->optionsReady = True;
	}
}

void PETScMatrixSolver_SetMaxIterations( void* matrixSolver, unsigned nIterations ) {
	PETScMatrixSolver*	self = (PETScMatrixSolver*)matrixSolver;
	PetscErrorCode		ec;

	assert( self && Stg_CheckType( self, PETScMatrixSolver ) );

	ec = KSPSetTolerances( self->ksp, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, (PetscInt)nIterations );
	CheckPETScError( ec );
}

void PETScMatrixSolver_SetRelativeTolerance( void* matrixSolver, double tolerance ) {
	PETScMatrixSolver*	self = (PETScMatrixSolver*)matrixSolver;
	PetscErrorCode		ec;

	assert( self && Stg_CheckType( self, PETScMatrixSolver ) );

	ec = KSPSetTolerances( self->ksp, tolerance, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT );
	CheckPETScError( ec );
}

void PETScMatrixSolver_SetAbsoluteTolerance( void* matrixSolver, double tolerance ) {
	PETScMatrixSolver*	self = (PETScMatrixSolver*)matrixSolver;
	PetscErrorCode		ec;

	assert( self && Stg_CheckType( self, PETScMatrixSolver ) );

	ec = KSPSetTolerances( self->ksp, PETSC_DEFAULT, tolerance, PETSC_DEFAULT, PETSC_DEFAULT );
	CheckPETScError( ec );
}

void PETScMatrixSolver_SetUseInitialSolution( void* matrixSolver, Bool state ) {
	PETScMatrixSolver*	self = (PETScMatrixSolver*)matrixSolver;
	PetscErrorCode		ec;

	assert( self && Stg_CheckType( self, PETScMatrixSolver ) );

	ec = KSPSetInitialGuessNonzero( self->ksp, (PetscTruth)state );
	CheckPETScError( ec );
}

void PETScMatrixSolver_Solve( void* matrixSolver, void* _rhs, void* _solution ) {
	PETScMatrixSolver*	self = (PETScMatrixSolver*)matrixSolver;
	PETScVector*		rhs = (PETScVector*)_rhs;
	PETScVector*		solution = (PETScVector*)_solution;
	PetscErrorCode		ec;

	assert( self && Stg_CheckType( self, PETScMatrixSolver ) );
	assert( solution && Stg_CheckType( solution, PETScVector ) );
	assert( rhs && Stg_CheckType( rhs, PETScVector ) );

	MatrixSolver_Setup( self, rhs, solution );
	ec = KSPSolve( self->ksp, rhs->petscVec, solution->petscVec );
	CheckPETScError( ec );
}

void PETScMatrixSolver_Setup( void* matrixSolver, void* rhs, void* solution ) {
	PETScMatrixSolver*	self = (PETScMatrixSolver*)matrixSolver;
	PetscErrorCode		ec;

	assert( self && Stg_CheckType( self, PETScMatrixSolver ) );

	_MatrixSolver_Setup( self, rhs, solution );
}

MatrixSolver_Status PETScMatrixSolver_GetSolveStatus( void* matrixSolver ) {
	PETScMatrixSolver*	self = (PETScMatrixSolver*)matrixSolver;
	PC			pc;
	KSPType			kspType;
	PCType			pcType;
	KSPConvergedReason	reason;
	PetscErrorCode		ec;

	assert( self && Stg_CheckType( self, PETScMatrixSolver ) );

	ec = KSPGetType( self->ksp, &kspType );
	CheckPETScError( ec );
	ec = KSPGetPC( self->ksp, &pc );
	CheckPETScError( ec );
	ec = PCGetType( pc, &pcType );
	CheckPETScError( ec );

	if( !strcmp( kspType, KSPRICHARDSON ) && !strcmp( pcType, PCSOR ) ) {
		double		rnorm;
		unsigned	curIt;

		rnorm = PETScMatrixSolver_GetResidualNorm( self );
		curIt = PETScMatrixSolver_GetIterations( self );
		PETScMatrixSolver_SetNormType( self, PETScMatrixSolver_NormType_Preconditioned );
		ec = KSPDefaultConverged( self->ksp, (PetscInt)curIt, (PetscScalar)rnorm, &reason, PETSC_NULL );
		CheckPETScError( ec );
	}
	else {
		ec = KSPGetConvergedReason( self->ksp, &reason );
		CheckPETScError( ec );
	}

	return reason;
}

unsigned PETScMatrixSolver_GetIterations( void* matrixSolver ) {
	PETScMatrixSolver*	self = (PETScMatrixSolver*)matrixSolver;
	PetscInt		nIts;
	PetscErrorCode		ec;

	assert( self && Stg_CheckType( self, PETScMatrixSolver ) );

	ec = KSPGetIterationNumber( self->ksp, &nIts );
	CheckPETScError( ec );

	return (unsigned)nIts;
}

unsigned PETScMatrixSolver_GetMaxIterations( void* matrixSolver ) {
	PETScMatrixSolver*	self = (PETScMatrixSolver*)matrixSolver;
	PetscInt		nIts;
	PetscErrorCode		ec;

	assert( self && Stg_CheckType( self, PETScMatrixSolver ) );

	ec = KSPGetTolerances( self->ksp, PETSC_NULL, PETSC_NULL, PETSC_NULL, &nIts );
	CheckPETScError( ec );

	return (unsigned)nIts;
}

double PETScMatrixSolver_GetResidualNorm( void* matrixSolver ) {
	PETScMatrixSolver*	self = (PETScMatrixSolver*)matrixSolver;
	PC			pc;
	KSPType			kspType;
	PCType			pcType;
	PetscScalar		rnorm;
	PetscErrorCode		ec;

	assert( self && Stg_CheckType( self, PETScMatrixSolver ) );

	ec = KSPGetType( self->ksp, &kspType );
	CheckPETScError( ec );
	ec = KSPGetPC( self->ksp, &pc );
	CheckPETScError( ec );
	ec = PCGetType( pc, &pcType );
	CheckPETScError( ec );

	if( !strcmp( kspType, KSPRICHARDSON ) && !strcmp( pcType, PCSOR ) ) {
		Vector*	residual;

		residual = MatrixSolver_GetResidual( self );
		rnorm = (PetscScalar)Vector_L2Norm( residual );
	}
	else {
		ec = KSPGetResidualNorm( self->ksp, &rnorm );
		CheckPETScError( ec );
	}

	return (double)rnorm;
}


/*--------------------------------------------------------------------------------------------------------------------------
** Public Functions
*/

void PETScMatrixSolver_SetKSPType( void* matrixSolver, PETScMatrixSolver_KSPType type ) {
	PETScMatrixSolver*	self = (PETScMatrixSolver*)matrixSolver;
	PetscErrorCode		ec;

	assert( self && Stg_CheckType( self, PETScMatrixSolver ) );

	switch( type ) {
	case PETScMatrixSolver_KSPType_Richardson:
		ec = KSPSetType( self->ksp, KSPRICHARDSON );
		CheckPETScError( ec );
		break;
	case PETScMatrixSolver_KSPType_GMRes:
		ec = KSPSetType( self->ksp, KSPGMRES );
		CheckPETScError( ec );
		break;
	case PETScMatrixSolver_KSPType_FGMRes:
		ec = KSPSetType( self->ksp, KSPFGMRES );
		CheckPETScError( ec );
		break;
	case PETScMatrixSolver_KSPType_CG:
		ec = KSPSetType( self->ksp, KSPCG );
		CheckPETScError( ec );
		break;
	case PETScMatrixSolver_KSPType_PreOnly:
		ec = KSPSetType( self->ksp, KSPPREONLY );
		CheckPETScError( ec );
		break;
	default:
		assert( 0 );
	}
}

void PETScMatrixSolver_SetPCType( void* matrixSolver, PETScMatrixSolver_PCType type ) {
	PETScMatrixSolver*	self = (PETScMatrixSolver*)matrixSolver;
	PC			pc;
	PetscErrorCode		ec;

	assert( self && Stg_CheckType( self, PETScMatrixSolver ) );

	ec = KSPGetPC( self->ksp, &pc );
	CheckPETScError( ec );
	if( type == PETScMatrixSolver_PCType_Jacobi ) {
		ec = PCSetType( pc, PCJACOBI );
		CheckPETScError( ec );
	}
	else if( type == PETScMatrixSolver_PCType_BlockJacobi ) {
		ec = PCSetType( pc, PCBJACOBI );
		CheckPETScError( ec );
	}
	else if( type == PETScMatrixSolver_PCType_SOR ) {
		ec = PCSetType( pc, PCSOR );
		CheckPETScError( ec );
	}
	else if( type == PETScMatrixSolver_PCType_LU ) {
		ec = PCSetType( pc, PCLU );
		CheckPETScError( ec );
	}
	else if( type == PETScMatrixSolver_PCType_RedundantLU ) {
		ec = PCSetType( pc, PCREDUNDANT );
		CheckPETScError( ec );
		ec = PCRedundantGetPC( pc, &pc );
		CheckPETScError( ec );
		ec = PCSetType( pc, PCLU );
		CheckPETScError( ec );
	}
	else if( type == PETScMatrixSolver_PCType_ILU ) {
		ec = PCSetType( pc, PCILU );
		CheckPETScError( ec );
	}
	else if( type == PETScMatrixSolver_PCType_Multigrid ) {
		ec = PCSetType( pc, PCMG );
		CheckPETScError( ec );
	}
	else if( type == PETScMatrixSolver_PCType_None ) {
		ec = PCSetType( pc, PCNONE );
		CheckPETScError( ec );
	}
	else {
		assert( 0 ); /* TODO */
	}
}

void PETScMatrixSolver_GetSubBlocks( void* matrixSolver, unsigned* nBlocks, KSP** ksps ) {
	PETScMatrixSolver*	self = (PETScMatrixSolver*)matrixSolver;
	PC			pc;
	PetscErrorCode		ec;

	assert( self && Stg_CheckType( self, PETScMatrixSolver ) );

	ec = KSPGetPC( self->ksp, &pc );
	CheckPETScError( ec );
	ec = PCBJacobiGetSubKSP( pc, (PetscInt*)nBlocks, PETSC_NULL, ksps );
	CheckPETScError( ec );
}

void PETScMatrixSolver_EnableShifting( void* matrixSolver, Bool state ) {
	PETScMatrixSolver*	self = (PETScMatrixSolver*)matrixSolver;
	PC			pc;
	PetscErrorCode		ec;

	assert( self && Stg_CheckType( self, PETScMatrixSolver ) );

	ec = KSPGetPC( self->ksp, &pc );
	CheckPETScError( ec );
	ec = PCFactorSetShiftPd( pc, (PetscTruth)state );
	CheckPETScError( ec );
}

void PETScMatrixSolver_SetNormType( void* matrixSolver, PETScMatrixSolver_NormType normType ) {
	PETScMatrixSolver*	self = (PETScMatrixSolver*)matrixSolver;
	PetscErrorCode		ec;

	assert( self && Stg_CheckType( self, PETScMatrixSolver ) );

	ec = KSPSetNormType( self->ksp, (KSPNormType)normType );
	CheckPETScError( ec );
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Private Functions
*/
