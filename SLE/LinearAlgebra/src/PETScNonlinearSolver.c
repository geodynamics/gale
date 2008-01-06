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
** $Id: PETScNonlinearSolver.c 3584 2006-05-16 11:11:07Z PatrickSunter $
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
const Type PETScNonlinearSolver_Type = "PETScNonlinearSolver";


/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/

PETScNonlinearSolver* PETScNonlinearSolver_New( Name name ) {
	return _PETScNonlinearSolver_New( sizeof(PETScNonlinearSolver), 
				       PETScNonlinearSolver_Type, 
				       _PETScNonlinearSolver_Delete, 
				       _PETScNonlinearSolver_Print, 
				       NULL, 
				       (void* (*)(Name))PETScNonlinearSolver_New, 
				       _PETScNonlinearSolver_Construct, 
				       _PETScNonlinearSolver_Build, 
				       _PETScNonlinearSolver_Initialise, 
				       _PETScNonlinearSolver_Execute, 
				       _PETScNonlinearSolver_Destroy, 
				       name, 
				       NON_GLOBAL, 
				       PETScNonlinearSolver_SetComm, 
				       //PETScNonlinearSolver_Create,
			               PETScNonlinearSolver_Destroy,
			               PETScNonlinearSolver_SetFunction,
			               PETScNonlinearSolver_GetJacobian,
			               PETScNonlinearSolver_Solve,
			               PETScNonlinearSolver_SetSolution,
				       PETScNonlinearSolver_SetRhs,       
				       NULL,//PETScNonlinearSolver_GetSolveStatus, 
				       PETScNonlinearSolver_GetIterations, 
				       PETScNonlinearSolver_GetMaxIterations, 
				       NULL/*PETScNonlinearSolver_GetResidualNorm*/ );
}

PETScNonlinearSolver* _PETScNonlinearSolver_New( PETSCNONLINEARSOLVER_DEFARGS ) {
	PETScNonlinearSolver*	self;

	/* Allocate memory */
	assert( sizeOfSelf >= sizeof(PETScNonlinearSolver) );
	self = (PETScNonlinearSolver*)_NonlinearSolver_New( NONLINEARSOLVER_PASSARGS );

	/* Virtual info */

	/* PETScNonlinearSolver info */
	_PETScNonlinearSolver_Init( self );

	return self;
}

void _PETScNonlinearSolver_Init( PETScNonlinearSolver* self ) {
	PETScMatrix*		J	= (PETScMatrix*)self->J;
	PetscErrorCode		ec;
	
	assert( self && Stg_CheckType( self, PETScNonlinearSolver ) );
	assert( J && Stg_CheckType( J, PETScMatrix ) );

	self->snes = PETSC_NULL;
	self->ksp = PETSC_NULL;
	PETScNonlinearSolver_SetComm( self, MPI_COMM_WORLD );

	ec = SNESCreate( self->comm, &self->snes );

	/* i assume this is ok?? */
	ec = KSPSetOperators( self->ksp, J->petscMat, J->petscMat, DIFFERENT_NONZERO_PATTERN );
	CheckPETScError( ec );
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _PETScNonlinearSolver_Delete( void* nls ) {
	PETScNonlinearSolver*	self = (PETScNonlinearSolver*)nls;

	assert( self && Stg_CheckType( self, PETScNonlinearSolver ) );

	KSPDestroy( self->ksp );

	/* Delete the parent. */
	_NonlinearSolver_Delete( self );
}

void _PETScNonlinearSolver_Print( void* nls, Stream* stream ) {
	PETScNonlinearSolver*	self = (PETScNonlinearSolver*)nls;
	
	/* Set the Journal for printing informations */
	Stream* nlsStream;
	nlsStream = Journal_Register( InfoStream_Type, "PETScNonlinearSolverStream" );

	assert( self && Stg_CheckType( self, PETScNonlinearSolver ) );

	/* Print parent */
	Journal_Printf( stream, "PETScNonlinearSolver (ptr): (%p)\n", self );
	_NonlinearSolver_Print( self, stream );
}

void _PETScNonlinearSolver_Construct( void* nls, Stg_ComponentFactory* cf, void* data ) {
	PETScNonlinearSolver*		self = (PETScNonlinearSolver*)nls;

	assert( self && Stg_CheckType( self, PETScNonlinearSolver ) );
	assert( cf );

	_NonlinearSolver_Construct( self, cf, data );
}

void _PETScNonlinearSolver_Build( void* nls, void* data ) {
}

void _PETScNonlinearSolver_Initialise( void* nls, void* data ) {
}

void _PETScNonlinearSolver_Execute( void* nls, void* data ) {
}

void _PETScNonlinearSolver_Destroy( void* nls, void* data ) {
	PETScNonlinearSolver*		self = (PETScNonlinearSolver*)nls;
	
	PETScNonlinearSolver_Destroy( self );
	_NonlinearSolver_Destroy( self, data );
}

void PETScNonlinearSolver_SetComm( void* nls, MPI_Comm comm ) {
	PETScNonlinearSolver*	self = (PETScNonlinearSolver*)nls;
	PetscErrorCode		ec;

	assert( self && Stg_CheckType( self, PETScNonlinearSolver ) );

	_NonlinearSolver_SetComm( self, comm );

	if( self->snes != PETSC_NULL )
		SNESDestroy( self->snes );
	ec = SNESCreate( self->comm, &self->snes );
	CheckPETScError( ec );

	if( self->ksp != PETSC_NULL )
		KSPDestroy( self->ksp );
	ec = KSPCreate( self->comm, &self->ksp );
	CheckPETScError( ec );
}

void PETScNonlinearSolver_Destroy( void* nls ) {
	PETScNonlinearSolver*	self = (PETScNonlinearSolver*)nls;
	PetscErrorCode		ec;

	assert( self && Stg_CheckType( self, PETScNonlinearSolver ) );

	ec = SNESDestroy( self->snes );
	CheckPETScError( ec );
	ec = KSPDestroy( self->ksp );
	CheckPETScError( ec );

}

void PETScNonlinearSolver_SetFunction( void* nls, void* _f, void* _func, void* context ) {
	PETScNonlinearSolver*		self 	= (PETScNonlinearSolver*)nls;
	PETScVector*			f	= (PETScVector*)_f;
	PETScNonlinearSolver_Func*	func	= (PETScNonlinearSolver_Func*)_func;	
	PetscErrorCode			ec;

	assert( self && Stg_CheckType( self, PETScNonlinearSolver ) );

	ec = SNESSetFunction( self->snes, f->petscVec, func, context );
	CheckPETScError( ec );
}

void PETScNonlinearSolver_GetJacobian( void* nls, void* _J, void* _pc, void** context ) {
	PETScNonlinearSolver*	self 	= (PETScNonlinearSolver*)nls;
	PETScMatrix*		J	= (PETScMatrix*)_J;
	PETScMatrix*		pc	= (PETScMatrix*)_pc;
	PetscErrorCode		ec;

	assert( self && Stg_CheckType( self, PETScNonlinearSolver ) );
	assert( J && Stg_CheckType( J, PETScMatrix ) );
	assert( pc && Stg_CheckType( pc, PETScMatrix ) );

	/* just using the matrix free PETSc routines for generating J for now, so not passing in any function here.
	 * will need to specify this argument when manually computing J however... 28.12.07 */
	ec = SNESGetJacobian( self->snes, &J->petscMat, &pc->petscMat, PETSC_NULL, context );
	CheckPETScError( ec );
}

void PETScNonlinearSolver_Solve( void* nls, void* _b, void* _x ) {
	PETScNonlinearSolver*	self 	= (PETScNonlinearSolver*)nls;
	PETScVector*		b 	= (PETScVector*)_b;
	PETScVector*		x 	= (PETScVector*)_x;
	PetscErrorCode		ec;

	assert( self && Stg_CheckType( self, PETScNonlinearSolver ) );
	assert( b && Stg_CheckType( b, PETScVector ) );
	assert( x && Stg_CheckType( x, PETScVector ) );

	//NonlinearSolver_Setup( self, rhs, solution );
	//need to setup for solve??
	ec = SNESSolve( self->snes, b->petscVec, x->petscVec );
	CheckPETScError( ec );
}

void PETScNonlinearSolver_SetSolution( void* nls, void* _x ) {
	PETScNonlinearSolver*	self 	= (PETScNonlinearSolver*)nls;
	PETScVector*		x 	= (PETScVector*)_x;
	PetscErrorCode		ec;

	assert( self && Stg_CheckType( self, PETScNonlinearSolver ) );
	assert( x && Stg_CheckType( x, PETScVector ) );

	ec = SNESSetSolution( self->snes, x->petscVec );
	CheckPETScError( ec );
}

void PETScNonlinearSolver_SetRhs( void* nls, void* _rhs ) {
	PETScNonlinearSolver*	self 	= (PETScNonlinearSolver*)nls;
	PETScVector*		rhs 	= (PETScVector*)_rhs;
	PetscErrorCode		ec;

	assert( self && Stg_CheckType( self, PETScNonlinearSolver ) );
	assert( rhs && Stg_CheckType( rhs, PETScVector ) );

	ec = SNESSetRhs( self->snes, rhs->petscVec );
	CheckPETScError( ec );
}

/*
void PETScNonlinearSolver_Setup( void* nls, void* rhs, void* solution ) {
	PETScNonlinearSolver*	self = (PETScNonlinearSolver*)nls;
	PetscErrorCode		ec;

	assert( self && Stg_CheckType( self, PETScNonlinearSolver ) );

	_NonlinearSolver_Setup( self, rhs, solution );

	ec = KSPSetFromOptions( self->ksp );
	CheckPETScError( ec );
}*/
/*
MatrixSolver_Status PETScNonlinearSolver_GetSolveStatus( void* nls ) {
	PETScNonlinearSolver*	self = (PETScNonlinearSolver*)nls;
	PC			pc;
	KSPType			kspType;
	PCType			pcType;
	KSPConvergedReason	reason;
	PetscErrorCode		ec;

	assert( self && Stg_CheckType( self, PETScNonlinearSolver ) );

	ec = KSPGetType( self->ksp, &kspType );
	CheckPETScError( ec );
	ec = KSPGetPC( self->ksp, &pc );
	CheckPETScError( ec );
	ec = PCGetType( pc, &pcType );
	CheckPETScError( ec );

	if( !strcmp( kspType, KSPRICHARDSON ) && !strcmp( pcType, PCSOR ) ) {
		double		rnorm;
		unsigned	curIt;

		rnorm = PETScNonlinearSolver_GetResidualNorm( self );
		curIt = PETScNonlinearSolver_GetIterations( self );
		PETScNonlinearSolver_SetNormType( self, PETScNonlinearSolver_NormType_Preconditioned );
		ec = KSPDefaultConverged( self->ksp, (PetscInt)curIt, (PetscScalar)rnorm, &reason, PETSC_NULL );
		CheckPETScError( ec );
	}
	else {
		ec = KSPGetConvergedReason( self->ksp, &reason );
		CheckPETScError( ec );
	}

	return reason;
}*/

unsigned PETScNonlinearSolver_GetIterations( void* nls ) {
	PETScNonlinearSolver*	self = (PETScNonlinearSolver*)nls;
	PetscInt		nIts;
	PetscErrorCode		ec;

	assert( self && Stg_CheckType( self, PETScNonlinearSolver ) );

	ec = KSPGetIterationNumber( self->ksp, &nIts );
	CheckPETScError( ec );

	return (unsigned)nIts;
}

unsigned PETScNonlinearSolver_GetMaxIterations( void* nls ) {
	PETScNonlinearSolver*	self = (PETScNonlinearSolver*)nls;
	PetscInt		nIts;
	PetscErrorCode		ec;

	assert( self && Stg_CheckType( self, PETScNonlinearSolver ) );

	ec = KSPGetTolerances( self->ksp, PETSC_NULL, PETSC_NULL, PETSC_NULL, &nIts );
	CheckPETScError( ec );

	return (unsigned)nIts;
}
/*
double PETScNonlinearSolver_GetResidualNorm( void* nls ) {
	PETScNonlinearSolver*	self = (PETScNonlinearSolver*)nls;
	PC			pc;
	KSPType			kspType;
	PCType			pcType;
	PetscScalar		rnorm;
	PetscErrorCode		ec;

	assert( self && Stg_CheckType( self, PETScNonlinearSolver ) );

	ec = KSPGetType( self->ksp, &kspType );
	CheckPETScError( ec );
	ec = KSPGetPC( self->ksp, &pc );
	CheckPETScError( ec );
	ec = PCGetType( pc, &pcType );
	CheckPETScError( ec );

	if( !strcmp( kspType, KSPRICHARDSON ) && !strcmp( pcType, PCSOR ) ) {
		Vector*	residual;

		residual = NonlinearSolver_GetResidual( self );
		rnorm = (PetscScalar)Vector_L2Norm( residual );
	}
	else {
		ec = KSPGetResidualNorm( self->ksp, &rnorm );
		CheckPETScError( ec );
	}

	return (double)rnorm;
}*/


/*--------------------------------------------------------------------------------------------------------------------------
** Public Functions
*/
/*
void PETScNonlinearSolver_SetKSPType( void* nls, PETScNonlinearSolver_KSPType type ) {
	PETScNonlinearSolver*	self = (PETScNonlinearSolver*)nls;
	PetscErrorCode		ec;

	assert( self && Stg_CheckType( self, PETScNonlinearSolver ) );

	switch( type ) {
	case PETScNonlinearSolver_KSPType_Richardson:
		ec = KSPSetType( self->ksp, KSPRICHARDSON );
		CheckPETScError( ec );
		break;
	case PETScNonlinearSolver_KSPType_GMRes:
		ec = KSPSetType( self->ksp, KSPGMRES );
		CheckPETScError( ec );
		break;
	case PETScNonlinearSolver_KSPType_CG:
		ec = KSPSetType( self->ksp, KSPCG );
		CheckPETScError( ec );
		break;
	case PETScNonlinearSolver_KSPType_PreOnly:
		ec = KSPSetType( self->ksp, KSPPREONLY );
		CheckPETScError( ec );
		break;
	default:
		assert( 0 );
	}
}*/
/*
void PETScNonlinearSolver_SetPCType( void* nls, PETScNonlinearSolver_PCType type ) {
	PETScNonlinearSolver*	self = (PETScNonlinearSolver*)nls;
	PC			pc;
	PetscErrorCode		ec;

	assert( self && Stg_CheckType( self, PETScNonlinearSolver ) );

	ec = KSPGetPC( self->ksp, &pc );
	CheckPETScError( ec );
	if( type == PETScNonlinearSolver_PCType_Jacobi ) {
		ec = PCSetType( pc, PCJACOBI );
		CheckPETScError( ec );
	}
	else if( type == PETScNonlinearSolver_PCType_BlockJacobi ) {
		ec = PCSetType( pc, PCBJACOBI );
		CheckPETScError( ec );
	}
	else if( type == PETScNonlinearSolver_PCType_SOR ) {
		ec = PCSetType( pc, PCSOR );
		CheckPETScError( ec );
	}
	else if( type == PETScNonlinearSolver_PCType_LU ) {
		ec = PCSetType( pc, PCLU );
		CheckPETScError( ec );
	}
	else if( type == PETScNonlinearSolver_PCType_RedundantLU ) {
		ec = PCSetType( pc, PCREDUNDANT );
		CheckPETScError( ec );
		ec = PCRedundantGetPC( pc, &pc );
		CheckPETScError( ec );
		ec = PCSetType( pc, PCLU );
		CheckPETScError( ec );
	}
	else if( type == PETScNonlinearSolver_PCType_ILU ) {
		ec = PCSetType( pc, PCILU );
		CheckPETScError( ec );
	}
	else if( type == PETScNonlinearSolver_PCType_Multigrid ) {
		ec = PCSetType( pc, PCMG );
		CheckPETScError( ec );
	}
	else if( type == PETScNonlinearSolver_PCType_None ) {
		ec = PCSetType( pc, PCNONE );
		CheckPETScError( ec );
	}
	else {
		assert( 0 ); 
	}
}*/
/*
void PETScNonlinearSolver_GetSubBlocks( void* nls, unsigned* nBlocks, KSP** ksps ) {
	PETScNonlinearSolver*	self = (PETScNonlinearSolver*)nls;
	PC			pc;
	PetscErrorCode		ec;

	assert( self && Stg_CheckType( self, PETScNonlinearSolver ) );

	ec = KSPGetPC( self->ksp, &pc );
	CheckPETScError( ec );
	ec = PCBJacobiGetSubKSP( pc, (PetscInt*)nBlocks, PETSC_NULL, ksps );
	CheckPETScError( ec );
}

void PETScNonlinearSolver_EnableShifting( void* nls, Bool state ) {
	PETScNonlinearSolver*	self = (PETScNonlinearSolver*)nls;
	PC			pc;
	PetscErrorCode		ec;

	assert( self && Stg_CheckType( self, PETScNonlinearSolver ) );

	ec = KSPGetPC( self->ksp, &pc );
	CheckPETScError( ec );
	ec = PCFactorSetShiftPd( pc, (PetscTruth)state );
	CheckPETScError( ec );
}

void PETScNonlinearSolver_SetNormType( void* nls, PETScNonlinearSolver_NormType normType ) {
	PETScNonlinearSolver*	self = (PETScNonlinearSolver*)nls;
	PetscErrorCode		ec;

	assert( self && Stg_CheckType( self, PETScNonlinearSolver ) );

	ec = KSPSetNormType( self->ksp, (KSPNormType)normType );
	CheckPETScError( ec );
}*/


/*----------------------------------------------------------------------------------------------------------------------------------
** Private Functions
*/
