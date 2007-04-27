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
** $Id: PETScMGSolver.c 3584 2006-05-16 11:11:07Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <mpi.h>
#include "StGermain/StGermain.h"
#include "Discretisation/Discretisation.h"
#include "LinearAlgebra.h"


/* Textual name of this class */
const Type PETScMGSolver_Type = "PETScMGSolver";


/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/

PETScMGSolver* PETScMGSolver_New( Name name ) {
	return _PETScMGSolver_New( sizeof(PETScMGSolver), 
				   PETScMGSolver_Type, 
				   _PETScMGSolver_Delete, 
				   _PETScMGSolver_Print, 
				   NULL, 
				   (void* (*)(Name))_PETScMGSolver_New, 
				   _PETScMGSolver_Construct, 
				   _PETScMGSolver_Build, 
				   _PETScMGSolver_Initialise, 
				   _PETScMGSolver_Execute, 
				   _PETScMGSolver_Destroy, 
				   name, 
				   NON_GLOBAL, 
				   PETScMGSolver_SetComm, 
				   PETScMatrixSolver_SetMatrix, 
				   PETScMatrixSolver_SetMaxIterations, 
				   PETScMatrixSolver_SetRelativeTolerance, 
				   PETScMatrixSolver_SetAbsoluteTolerance, 
				   PETScMatrixSolver_SetUseInitialSolution, 
				   PETScMatrixSolver_Solve, 
				   PETScMGSolver_Setup, 
				   PETScMatrixSolver_GetSolveStatus, 
				   PETScMatrixSolver_GetIterations, 
				   PETScMatrixSolver_GetMaxIterations, 
				   PETScMatrixSolver_GetResidualNorm );
}

PETScMGSolver* _PETScMGSolver_New( MULTIGRIDSOLVER_DEFARGS ) {
	PETScMGSolver*	self;

	/* Allocate memory */
	assert( sizeOfSelf >= sizeof(PETScMGSolver) );
	self = (PETScMGSolver*)_PETScMatrixSolver_New( PETSCMATRIXSOLVER_PASSARGS );

	/* Virtual info */

	/* PETScMGSolver info */
	_PETScMGSolver_Init( self );

	return self;
}

void _PETScMGSolver_Init( PETScMGSolver* self ) {
	assert( self && Stg_CheckType( self, PETScMGSolver ) );

	self->nLevels = 0;
	self->levels = NULL;
	self->opGen = NULL;
	self->solversChanged = True;
	self->opsChanged = True;
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _PETScMGSolver_Delete( void* matrixSolver ) {
	PETScMGSolver*	self = (PETScMGSolver*)matrixSolver;

	assert( self && Stg_CheckType( self, PETScMGSolver ) );

	PETScMGSolver_Destruct( self );

	/* Delete the parent. */
	_PETScMatrixSolver_Delete( self );
}

void _PETScMGSolver_Print( void* matrixSolver, Stream* stream ) {
	PETScMGSolver*	self = (PETScMGSolver*)matrixSolver;
	
	/* Set the Journal for printing informations */
	Stream* matrixSolverStream;
	matrixSolverStream = Journal_Register( InfoStream_Type, "PETScMGSolverStream" );

	assert( self && Stg_CheckType( self, PETScMGSolver ) );

	/* Print parent */
	Journal_Printf( stream, "PETScMGSolver (ptr): (%p)\n", self );
	_PETScMatrixSolver_Print( self, stream );
}

void _PETScMGSolver_Construct( void* matrixSolver, Stg_ComponentFactory* cf, void* data ) {
	PETScMGSolver*	self = (PETScMGSolver*)matrixSolver;
	Bool		pure;
	unsigned	nLevels;
	unsigned	nCycles;
	unsigned	nDownIts, nUpIts;

	assert( self && Stg_CheckType( self, PETScMGSolver ) );
	assert( cf );

	_PETScMatrixSolver_Construct( self, cf, data );

	pure = Stg_ComponentFactory_GetBool( cf, self->name, "pure", True );
	nLevels = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, "levels", 1 );
	nCycles = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, "cycles", 1 );
	nDownIts = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, "downIterations", 1 );
	nUpIts = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, "upIterations", 1 );

	self->pure = pure;
	PETScMGSolver_SetLevels( self, nLevels );
	PETScMGSolver_SetLevelCycles( self, nLevels - 1, nCycles );
	PETScMGSolver_SetAllDownIterations( self, nDownIts );
	PETScMGSolver_SetAllUpIterations( self, nUpIts );

	self->opGen = Stg_ComponentFactory_ConstructByKey( cf, self->name, "opGenerator", MGOpGenerator, 
							   True, data );
	MGOpGenerator_SetMatrixSolver( self->opGen, self );
	MGOpGenerator_SetNumLevels( self->opGen, nLevels );
}

void _PETScMGSolver_Build( void* matrixSolver, void* data ) {
	PETScMGSolver*	self = (PETScMGSolver*)matrixSolver;

	assert( self && Stg_CheckType( self, PETScMGSolver ) );

	if( self->opGen )
		Build( self->opGen, data, False );
}

void _PETScMGSolver_Initialise( void* matrixSolver, void* data ) {
	PETScMGSolver*	self = (PETScMGSolver*)matrixSolver;

	assert( self && Stg_CheckType( self, PETScMGSolver ) );

	if( self->opGen )
		Initialise( self->opGen, data, False );
}

void _PETScMGSolver_Execute( void* matrixSolver, void* data ) {
}

void _PETScMGSolver_Destroy( void* matrixSolver, void* data ) {
}

void PETScMGSolver_SetComm( void* matrixSolver, MPI_Comm comm ) {
	PETScMGSolver*	self = (PETScMGSolver*)matrixSolver;

	assert( self && Stg_CheckType( self, PETScMGSolver ) );

	PETScMatrixSolver_SetComm( self, comm );

	if( self->pure )
		PETScMatrixSolver_SetKSPType( self, PETScMatrixSolver_KSPType_Richardson );
	else
		PETScMatrixSolver_SetKSPType( self, PETScMatrixSolver_KSPType_GMRes );
	PETScMatrixSolver_SetPCType( self, PETScMatrixSolver_PCType_Multigrid );
}

void PETScMGSolver_Setup( void* matrixSolver, void* rhs, void* solution ) {
	PETScMGSolver*	self = (PETScMGSolver*)matrixSolver;
	Bool		rebuildOps;

	assert( self && Stg_CheckType( self, PETScMGSolver ) );

	_MatrixSolver_Setup( self, rhs, solution );

	if( self->opGen )
		rebuildOps = MGOpGenerator_HasExpired( self->opGen );
	else
		rebuildOps = False;
	if( !rebuildOps ) {
		unsigned		l_i;

		for( l_i = 1; l_i < self->nLevels; l_i++ ) {
			if( !self->levels[l_i].R || !self->levels[l_i].P ) {
				rebuildOps = True;
				break;
			}
		}
	}

	if( self->solversChanged || rebuildOps || self->opsChanged )
		PETScMGSolver_UpdateSolvers( self );
	if( rebuildOps )
		PETScMGSolver_UpdateOps( self );
	if( self->matrixChanged || self->solversChanged || rebuildOps || self->opsChanged ) {
		PETScMGSolver_UpdateMatrices( self );
		PETScMGSolver_UpdateWorkVectors( self );
	}

	self->solversChanged = False;
	self->matrixChanged = False;
	self->opsChanged = False;
}


/*--------------------------------------------------------------------------------------------------------------------------
** Public Functions
*/

void PETScMGSolver_SetLevels( void* matrixSolver, unsigned nLevels ) {
	PETScMGSolver*	self = (PETScMGSolver*)matrixSolver;
	unsigned	l_i;

	assert( self && Stg_CheckType( self, PETScMGSolver ) );

	PETScMGSolver_DestructLevels( self );

	self->nLevels = nLevels;
	self->levels = AllocArray( PETScMGSolver_Level, nLevels );
	for( l_i = 0; l_i < nLevels; l_i++ ) {
		PETScMGSolver_Level*	level = self->levels + l_i;

		level->nDownIts = 1;
		level->nUpIts = (l_i == 0) ? 0 : 1;
		level->nCycles = 1;
		level->R = NULL;
		level->P = NULL;
		level->A = NULL;
		level->workRes = NULL;
		level->workSol = NULL;
		level->workRHS = NULL;
	}
}

void PETScMGSolver_SetRestriction( void* matrixSolver, unsigned levelInd, void* R ) {
	PETScMGSolver*		self = (PETScMGSolver*)matrixSolver;
	PETScMGSolver_Level*	level;

	assert( self && Stg_CheckType( self, PETScMGSolver ) );
	assert( levelInd < self->nLevels && levelInd > 0 );
	assert( !R || Stg_CheckType( R, PETScMatrix ) );

	level = self->levels + levelInd;
	if( level->R != R )
		self->opsChanged = True;
	if( R )
		Stg_Class_AddRef( R );
	if( level->R )
		Stg_Class_RemoveRef( level->R );
	level->R = R;
}

void PETScMGSolver_SetProlongation( void* matrixSolver, unsigned levelInd, void* P ) {
	PETScMGSolver*		self = (PETScMGSolver*)matrixSolver;
	PETScMGSolver_Level*	level;

	assert( self && Stg_CheckType( self, PETScMGSolver ) );
	assert( levelInd < self->nLevels && levelInd > 0 );
	assert( !P || Stg_CheckType( P, PETScMatrix ) );

	level = self->levels + levelInd;
	if( level->P != P )
		self->opsChanged = True;
	if( P )
		Stg_Class_AddRef( P );
	if( level->P )
		Stg_Class_RemoveRef( level->P );
	level->P = P;
}

void PETScMGSolver_SetLevelDownIterations( void* matrixSolver, unsigned level, unsigned nIts ) {
	PETScMGSolver*	self = (PETScMGSolver*)matrixSolver;

	assert( self && Stg_CheckType( self, PETScMGSolver ) );
	assert( level < self->nLevels );

	self->levels[level].nDownIts = nIts;
	self->solversChanged = True;
}

void PETScMGSolver_SetLevelUpIterations( void* matrixSolver, unsigned level, unsigned nIts ) {
	PETScMGSolver*	self = (PETScMGSolver*)matrixSolver;

	assert( self && Stg_CheckType( self, PETScMGSolver ) );
	assert( level < self->nLevels );
	assert( level > 0 );

	self->levels[level].nUpIts = nIts;
	self->solversChanged = True;
}

void PETScMGSolver_SetLevelCycles( void* matrixSolver, unsigned level, unsigned nCycles ) {
	PETScMGSolver*	self = (PETScMGSolver*)matrixSolver;

	assert( self && Stg_CheckType( self, PETScMGSolver ) );
	assert( level < self->nLevels );
	assert( level > 0 );

	self->levels[level].nCycles = nCycles;
	self->solversChanged = True;
}

void PETScMGSolver_SetAllDownIterations( void* matrixSolver, unsigned nIts ) {
	PETScMGSolver*	self = (PETScMGSolver*)matrixSolver;
	unsigned	l_i;

	assert( self && Stg_CheckType( self, PETScMGSolver ) );

	for( l_i = 1; l_i < self->nLevels; l_i++ )
		PETScMGSolver_SetLevelDownIterations( self, l_i, nIts );
}

void PETScMGSolver_SetAllUpIterations( void* matrixSolver, unsigned nIts ) {
	PETScMGSolver*	self = (PETScMGSolver*)matrixSolver;
	unsigned	l_i;

	assert( self && Stg_CheckType( self, PETScMGSolver ) );

	for( l_i = 1; l_i < self->nLevels; l_i++ )
		PETScMGSolver_SetLevelUpIterations( self, l_i, nIts );
}

unsigned PETScMGSolver_GetNumLevels( void* matrixSolver ) {
	PETScMGSolver*	self = (PETScMGSolver*)matrixSolver;

	assert( self && Stg_CheckType( self, PETScMGSolver ) );

	return self->nLevels;
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Private Functions
*/

void PETScMGSolver_UpdateOps( PETScMGSolver* self ) {
	PC		pc;
	PETScMatrix	**pOps, **rOps;
	PetscErrorCode	ec;
	unsigned	l_i;

	assert( self && Stg_CheckType( self, PETScMGSolver ) );

	ec = KSPGetPC( self->ksp, &pc );
	CheckPETScError( ec );

	MGOpGenerator_Generate( self->opGen, (Matrix***)&pOps, (Matrix***)&rOps );

	for( l_i = 1; l_i < self->nLevels; l_i++ ) {
		assert( Stg_CheckType( pOps[l_i], PETScMatrix ) );
		assert( Stg_CheckType( rOps[l_i], PETScMatrix ) );

		PETScMGSolver_SetProlongation( self, l_i, pOps[l_i] );
		ec = PCMGSetInterpolate( pc, l_i, pOps[l_i]->petscMat );
		CheckPETScError( ec );

		PETScMGSolver_SetRestriction( self, l_i, rOps[l_i] );
		ec = PCMGSetRestriction( pc, l_i, rOps[l_i]->petscMat );
		CheckPETScError( ec );
	}

	FreeArray( pOps );
	FreeArray( rOps );
}

void PETScMGSolver_UpdateMatrices( PETScMGSolver* self ) {
	Stream*			stream;
	PETScMGSolver_Level*	level;
	PC			pc;
	KSP			levelKSP;
	PetscErrorCode		ec;
	unsigned		l_i;

	assert( self && Stg_CheckType( self, PETScMGSolver ) );
	assert( self->matrix && Stg_CheckType( self->matrix, PETScMatrix ) );

	stream = Journal_Register( InfoStream_Type, "general" ); assert( stream );
	Journal_Printf( stream, "Updating MG matrices ...\n" );

	ec = KSPGetPC( self->ksp, &pc );
	CheckPETScError( ec );

	for( l_i = self->nLevels - 1; l_i < self->nLevels; l_i-- ) {
		level = self->levels + l_i;

		if( l_i == self->nLevels - 1 )
			level->A = (PETScMatrix*)self->matrix;
		else {
			KillObject( level->A );
			Matrix_PtAP( self->levels[l_i + 1].A, self->levels[l_i + 1].P, (void**)&level->A );
		}

		ec = PCMGGetSmootherDown( pc, l_i, &levelKSP );
		CheckPETScError( ec );
		ec = KSPSetOperators( levelKSP, level->A->petscMat, level->A->petscMat, DIFFERENT_NONZERO_PATTERN );
		CheckPETScError( ec );
		ec = PCMGSetResidual( pc, l_i, PCMGDefaultResidual, level->A->petscMat );
		CheckPETScError( ec );

		if( l_i > 0 ) {
			PCMGGetSmootherUp( pc, l_i, &levelKSP );
			ec = KSPSetOperators( levelKSP, level->A->petscMat, level->A->petscMat, DIFFERENT_NONZERO_PATTERN );
			CheckPETScError( ec );
		}
	}

	Journal_Printf( stream, "done\n" );
}

void PETScMGSolver_UpdateWorkVectors( PETScMGSolver* self ) {
	PETScMGSolver_Level*	level;
	PC			pc;
	unsigned		size;
	PetscErrorCode		ec;
	unsigned		l_i;

	assert( self && Stg_CheckType( self, PETScMGSolver ) );

	if( self->nLevels == 1 )
		return;

	ec = KSPGetPC( self->ksp, &pc );
	CheckPETScError( ec );

	for( l_i = 0; l_i < self->nLevels; l_i++ ) {
		level = self->levels + l_i;

		Matrix_GetLocalSize( level->A, &size, NULL );

		if( l_i > 0 && (!level->workRes || Vector_GetLocalSize( level->workRes ) != size) ) {
			if( level->workRes )
				FreeObject( level->workRes );
			Vector_Duplicate( self->curSolution, (void**)&level->workRes );
			Vector_SetLocalSize( level->workRes, size );
			ec = PCMGSetR( pc, l_i, level->workRes->petscVec );
			CheckPETScError( ec );
		}

		if( l_i < self->nLevels - 1 ) {
			if( !level->workSol || Vector_GetLocalSize( level->workSol ) != size ) {
				if( level->workSol )
					FreeObject( level->workSol );
				Vector_Duplicate( self->curSolution, (void**)&level->workSol );
				Vector_SetLocalSize( level->workSol, size );
				ec = PCMGSetX( pc, l_i, level->workSol->petscVec );
				CheckPETScError( ec );
			}

			if( !level->workRHS || Vector_GetLocalSize( level->workRHS ) != size ) {
				if( level->workRHS )
					FreeObject( level->workRHS );
				Vector_Duplicate( self->curSolution, (void**)&level->workRHS );
				Vector_SetLocalSize( level->workRHS, size );
				ec = PCMGSetRhs( pc, l_i, level->workRHS->petscVec );
				CheckPETScError( ec );
			}
		}
	}
}

void PETScMGSolver_UpdateSolvers( PETScMGSolver* self ) {
	PETScMGSolver_Level*	level;
	PC			pc;
	KSP			levelKSP;
	PC			levelPC;
	PetscErrorCode		ec;
	unsigned		l_i;

	assert( self && Stg_CheckType( self, PETScMGSolver ) );

	ec = KSPGetPC( self->ksp, &pc );
	CheckPETScError( ec );

	ec = PCMGSetLevels( pc, self->nLevels, PETSC_NULL );
	CheckPETScError( ec );
	ec = PCMGSetType( pc, PC_MG_MULTIPLICATIVE );
	CheckPETScError( ec );

	for( l_i = 1; l_i < self->nLevels; l_i++ ) {
		level = self->levels + l_i;

		ec = PCMGGetSmootherDown( pc, l_i, &levelKSP );
		CheckPETScError( ec );
		ec = KSPSetType( levelKSP, KSPRICHARDSON );
		CheckPETScError( ec );
		ec = KSPGetPC( levelKSP, &levelPC );
		CheckPETScError( ec );
		ec = PCSetType( levelPC, PCSOR );
		CheckPETScError( ec );
		ec = KSPSetTolerances( levelKSP, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, level->nDownIts );
		CheckPETScError( ec );
		if( l_i == self->nLevels - 1 )
			ec = KSPSetInitialGuessNonzero( levelKSP, PETSC_TRUE );
		else
			ec = KSPSetInitialGuessNonzero( levelKSP, PETSC_FALSE );
		CheckPETScError( ec );

		ec = PCMGGetSmootherUp( pc, l_i, &levelKSP );
		CheckPETScError( ec );
		ec = KSPSetType( levelKSP, KSPRICHARDSON );
		CheckPETScError( ec );
		ec = KSPGetPC( levelKSP, &levelPC );
		CheckPETScError( ec );
		ec = PCSetType( levelPC, PCSOR );
		CheckPETScError( ec );
		ec = KSPSetTolerances( levelKSP, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, level->nUpIts );
		CheckPETScError( ec );
		ec = KSPSetInitialGuessNonzero( levelKSP, PETSC_TRUE );
		CheckPETScError( ec );

		ec = PCMGSetCyclesOnLevel( pc, l_i, level->nCycles );
		CheckPETScError( ec );
	}
}

void PETScMGSolver_Destruct( PETScMGSolver* self ) {
	assert( self && Stg_CheckType( self, PETScMGSolver ) );

	PETScMGSolver_DestructLevels( self );
	KillObject( self->opGen );
}

void PETScMGSolver_DestructLevels( PETScMGSolver* self ) {
	unsigned	l_i;

	assert( self && Stg_CheckType( self, PETScMGSolver ) );

	for( l_i = 0; l_i < self->nLevels; l_i++ ) {
		PETScMGSolver_Level*	level = self->levels + l_i;

		if( level->R )
			Stg_Class_RemoveRef( level->R );
		if( level->P )
			Stg_Class_RemoveRef( level->P );
		if( level->A )
			Stg_Class_RemoveRef( level->A );

		FreeObject( level->workRes );
		FreeObject( level->workSol );
		FreeObject( level->workRHS );
	}

	KillArray( self->levels );
	self->nLevels = 0;
	self->solversChanged = True;
	self->opsChanged = True;
}
