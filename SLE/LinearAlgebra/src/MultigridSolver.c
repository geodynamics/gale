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
** $Id: MultigridSolver.c 3584 2006-05-16 11:11:07Z PatrickSunter $
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
const Type MultigridSolver_Type = "MultigridSolver";


/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/

MultigridSolver* MultigridSolver_New( Name name ) {
	return _MultigridSolver_New( sizeof(MultigridSolver), 
				     MultigridSolver_Type, 
				     _MultigridSolver_Delete, 
				     _MultigridSolver_Print, 
				     NULL, 
				     (void* (*)(Name))_MultigridSolver_New, 
				     _MultigridSolver_Construct, 
				     _MultigridSolver_Build, 
				     _MultigridSolver_Initialise, 
				     _MultigridSolver_Execute, 
				     _MultigridSolver_Destroy, 
				     name, 
				     NON_GLOBAL, 
				     _MatrixSolver_SetComm, 
				     _MatrixSolver_SetMatrix, 
				     MultigridSolver_SetMaxIterations, 
				     MultigridSolver_SetRelativeTolerance, 
				     MultigridSolver_SetAbsoluteTolerance, 
				     MultigridSolver_SetUseInitialSolution, 
				     MultigridSolver_Solve, 
				     MultigridSolver_Setup, 
				     MultigridSolver_GetSolveStatus, 
				     MultigridSolver_GetIterations, 
				     MultigridSolver_GetMaxIterations, 
				     MultigridSolver_GetResidualNorm );
}

MultigridSolver* _MultigridSolver_New( MULTIGRIDSOLVER_DEFARGS ) {
	MultigridSolver*	self;

	/* Allocate memory */
	assert( sizeOfSelf >= sizeof(MultigridSolver) );
	self = (MultigridSolver*)_MatrixSolver_New( MATRIXSOLVER_PASSARGS );

	/* Virtual info */

	/* MultigridSolver info */
	_MultigridSolver_Init( self );

	return self;
}

void _MultigridSolver_Init( MultigridSolver* self ) {
	assert( self && Stg_CheckType( self, MultigridSolver ) );

	self->nLevels = 0;
	self->levels = NULL;
	self->opGen = NULL;
	self->solversChanged = True;
	self->opsChanged = True;

	self->curIt = 0;
	self->maxIts = 500;
	self->relTol = 1e-5;
	self->rnorm = 0.0;
	self->useInitial = False;
	self->outerSolver = NULL;
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _MultigridSolver_Delete( void* matrixSolver ) {
	MultigridSolver*	self = (MultigridSolver*)matrixSolver;

	assert( self && Stg_CheckType( self, MultigridSolver ) );

	MultigridSolver_Destruct( self );

	/* Delete the parent. */
	_MatrixSolver_Delete( self );
}

void _MultigridSolver_Print( void* matrixSolver, Stream* stream ) {
	MultigridSolver*	self = (MultigridSolver*)matrixSolver;
	
	/* Set the Journal for printing informations */
	Stream* matrixSolverStream;
	matrixSolverStream = Journal_Register( InfoStream_Type, "MultigridSolverStream" );

	assert( self && Stg_CheckType( self, MultigridSolver ) );

	/* Print parent */
	Journal_Printf( stream, "MultigridSolver (ptr): (%p)\n", self );
	_MatrixSolver_Print( self, stream );
}

void _MultigridSolver_Construct( void* matrixSolver, Stg_ComponentFactory* cf, void* data ) {
	MultigridSolver*	self = (MultigridSolver*)matrixSolver;
	unsigned		nLevels;

	assert( self && Stg_CheckType( self, MultigridSolver ) );
	assert( cf );

	_MatrixSolver_Construct( self, cf, data );

	nLevels = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, "levels", 1 );
	MultigridSolver_SetLevels( self, nLevels );
	self->opGen = Stg_ComponentFactory_ConstructByKey( cf, self->name, "opGenerator", MGOpGenerator, 
							   True, data );
	MGOpGenerator_SetMatrixSolver( self->opGen, self );
	MGOpGenerator_SetNumLevels( self->opGen, nLevels );
}

void _MultigridSolver_Build( void* matrixSolver, void* data ) {
	MultigridSolver*	self = (MultigridSolver*)matrixSolver;

	assert( self && Stg_CheckType( self, MultigridSolver ) );

	self->stream = Journal_Register( InfoStream_Type, "general" );
	if( self->opGen )
		Stg_Component_Build( self->opGen, data, False );
}

void _MultigridSolver_Initialise( void* matrixSolver, void* data ) {
	MultigridSolver*	self = (MultigridSolver*)matrixSolver;

	assert( self && Stg_CheckType( self, MultigridSolver ) );

	if( self->opGen )
		Stg_Component_Initialise( self->opGen, data, False );
}

void _MultigridSolver_Execute( void* matrixSolver, void* data ) {
}

void _MultigridSolver_Destroy( void* matrixSolver, void* data ) {
}

void MultigridSolver_SetMaxIterations( void* matrixSolver, unsigned nIterations ) {
	MultigridSolver*	self = (MultigridSolver*)matrixSolver;

	assert( self && Stg_CheckType( self, MultigridSolver ) );

	self->maxIts = nIterations;
}

void MultigridSolver_SetRelativeTolerance( void* matrixSolver, double tolerance ) {
	MultigridSolver*	self = (MultigridSolver*)matrixSolver;

	assert( self && Stg_CheckType( self, MultigridSolver ) );

	self->relTol = tolerance;
	self->solversChanged = True;
}

void MultigridSolver_SetAbsoluteTolerance( void* matrixSolver, double tolerance ) {
	MultigridSolver*	self = (MultigridSolver*)matrixSolver;

	assert( self && Stg_CheckType( self, MultigridSolver ) );

	MatrixSolver_SetAbsoluteTolerance( self->outerSolver, tolerance );
	self->solversChanged = True;
}

void MultigridSolver_SetUseInitialSolution( void* matrixSolver, Bool state ) {
	MultigridSolver*	self = (MultigridSolver*)matrixSolver;

	assert( self && Stg_CheckType( self, MultigridSolver ) );

	self->useInitial = state;
}

void MultigridSolver_Solve( void* matrixSolver, void* _rhs, void* _solution ) {
	MultigridSolver*	self = (MultigridSolver*)matrixSolver;
	Vector*			rhs = (Vector*)_rhs;
	Vector*			solution = (Vector*)_solution;

	assert( self && Stg_CheckType( self, MultigridSolver ) );
	assert( rhs && Stg_CheckType( rhs, Vector ) );
	assert( solution && Stg_CheckType( solution, Vector ) );

	Journal_Printf( self->stream, "MultigridSolver: Starting solve ...\n" );
	Stream_Indent( self->stream );

	if( !self->useInitial ) {
		Journal_Printf( self->stream, "Zeroing initial solution\n" );
		Vector_Zero( solution );
	}
	else
		Journal_Printf( self->stream, "Keeping initial solution\n" );

	MultigridSolver_Setup( self, rhs, solution );
	self->curIt = 0;
	while( MatrixSolver_GetSolveStatus( self ) == MatrixSolver_Status_Iterating && self->curIt < self->maxIts ) {
		Journal_Printf( self->stream, "Iteration %d: residual %.10lf\n", 
				self->curIt, MatrixSolver_GetResidualNorm( self->outerSolver ) );
		MultigridSolver_LevelCycle( self, self->nLevels - 1, rhs, solution );
		self->curIt++;
	}
	Journal_Printf( self->stream, "Iteration %d: residual %.10lf\n", 
			self->curIt, MatrixSolver_GetResidualNorm( self->outerSolver ) );

	Stream_UnIndent( self->stream );
	Journal_Printf( self->stream, "done.\n" );
}

void MultigridSolver_Setup( void* matrixSolver, void* rhs, void* solution ) {
	MultigridSolver*	self = (MultigridSolver*)matrixSolver;
	Bool			rebuildOps;

	assert( self && Stg_CheckType( self, MultigridSolver ) );

	_MatrixSolver_Setup( self, rhs, solution );

	/* Need to rebuild the operators? */
	if( self->opGen )
		rebuildOps = MGOpGenerator_HasExpired( self->opGen );
	else
		rebuildOps = False;
	if( !rebuildOps ) {
		unsigned		l_i;

		for( l_i = 1; l_i < self->nLevels; l_i++ ) {
			if( !MultigridSolver_GetRestriction( self, l_i ) || !MultigridSolver_GetProlongation( self, l_i ) ) {
				rebuildOps = True;
				break;
			}
		}
	}

	if( rebuildOps )
		MultigridSolver_UpdateOps( self );
	if( self->matrixChanged || rebuildOps || self->opsChanged || self->solversChanged ) {
		MultigridSolver_UpdateMatrices( self );
		MultigridSolver_UpdateWorkVectors( self );
	}
	self->solversChanged = False;
	self->matrixChanged = False;
	self->opsChanged = False;
}

MatrixSolver_Status MultigridSolver_GetSolveStatus( void* matrixSolver ) {
	MultigridSolver*	self = (MultigridSolver*)matrixSolver;
	MatrixSolver_Status	outerStatus;

	assert( self && Stg_CheckType( self, MultigridSolver ) );

	MatrixSolver_Solve( self->outerSolver, self->curRHS, self->curSolution );
	outerStatus = MatrixSolver_GetSolveStatus( self->outerSolver );

	if( outerStatus == MatrixSolver_Status_ConvergedIterations || 
	    outerStatus == MatrixSolver_Status_DivergedIterations || 
	    outerStatus == MatrixSolver_Status_Iterating )
	{
		return MatrixSolver_Status_Iterating;
	}
	else if( outerStatus == MatrixSolver_Status_ConvergedRelative )
		return MatrixSolver_Status_ConvergedRelative;
	else
		return outerStatus;
}

unsigned MultigridSolver_GetIterations( void* matrixSolver ) {
	MultigridSolver*	self = (MultigridSolver*)matrixSolver;

	assert( self && Stg_CheckType( self, MultigridSolver ) );

	return self->curIt;
}

unsigned MultigridSolver_GetMaxIterations( void* matrixSolver ) {
	MultigridSolver*	self = (MultigridSolver*)matrixSolver;

	assert( self && Stg_CheckType( self, MultigridSolver ) );

	return self->maxIts;
}

double MultigridSolver_GetResidualNorm( void* matrixSolver ) {
	MultigridSolver*	self = (MultigridSolver*)matrixSolver;

	assert( self && Stg_CheckType( self, MultigridSolver ) );

	return MatrixSolver_GetResidualNorm( self->outerSolver );
}


/*--------------------------------------------------------------------------------------------------------------------------
** Public Functions
*/

void MultigridSolver_SetLevels( void* matrixSolver, unsigned nLevels ) {
	MultigridSolver*	self = (MultigridSolver*)matrixSolver;
	unsigned		nProcs;
	unsigned		l_i;

	assert( self && Stg_CheckType( self, MultigridSolver ) );

	MultigridSolver_DestructLevels( self );

	self->nLevels = nLevels;
	self->levels = AllocArray( MultigridSolver_Level, nLevels );
	MPI_Comm_size( self->comm, (int*)&nProcs );

	for( l_i = 0; l_i < nLevels; l_i++ ) {
		MultigridSolver_Level*	level = self->levels + l_i;

		level->downSolver = NULL;
		level->nDownIts = 1;
		level->upSolver = NULL;
		level->nUpIts = l_i ? 1 : 0;
		level->nCycles = 1;

		level->A = NULL;
		level->R = NULL;
		level->P = NULL;

		level->workRHS = NULL;
		level->workSol = NULL;
	}
}

void MultigridSolver_SetRestriction( void* matrixSolver, unsigned levelInd, void* _R ) {
	MultigridSolver*	self = (MultigridSolver*)matrixSolver;
	MultigridSolver_Level*	level;
	Matrix*			R = (Matrix*)_R;

	assert( self && Stg_CheckType( self, MultigridSolver ) );
	assert( levelInd < self->nLevels && levelInd > 0 );
	assert( !R || Stg_CheckType( R, Matrix ) );

	level = self->levels + levelInd;
	if( level->R != R )
		self->opsChanged = True;
	if( level->R )
		Stg_Class_RemoveRef( level->R );
	level->R = R;
	if( R )
		Stg_Class_AddRef( R );
}

void MultigridSolver_SetProlongation( void* matrixSolver, unsigned levelInd, void* P ) {
	MultigridSolver*	self = (MultigridSolver*)matrixSolver;
	MultigridSolver_Level*	level;

	assert( self && Stg_CheckType( self, MultigridSolver ) );
	assert( levelInd < self->nLevels && levelInd > 0 );
	assert( !P || Stg_CheckType( P, Matrix ) );

	level = self->levels + levelInd;
	if( level->P != P )
		self->opsChanged = True;
	if( level->P )
		Stg_Class_RemoveRef( level->P );
	level->P = P;
	if( P )
		Stg_Class_AddRef( P );
}

void MultigridSolver_SetLevelDownSolver( void* matrixSolver, unsigned levelInd, void* solver ) {
	MultigridSolver*	self = (MultigridSolver*)matrixSolver;
	MultigridSolver_Level*	level;

	assert( self && Stg_CheckType( self, MultigridSolver ) );
	assert( levelInd < self->nLevels );
	assert( !solver || Stg_CheckType( solver, MatrixSolver ) );

	level = self->levels + levelInd;
	if( level->downSolver != solver )
		self->solversChanged = True;
	if( level->downSolver )
		Stg_Class_RemoveRef( level->downSolver );
	level->downSolver = solver;
	if( solver )
		Stg_Class_AddRef( solver );
}

void MultigridSolver_SetLevelDownIterations( void* matrixSolver, unsigned level, unsigned nIts ) {
	MultigridSolver*	self = (MultigridSolver*)matrixSolver;

	assert( self && Stg_CheckType( self, MultigridSolver ) );
	assert( level < self->nLevels );

	self->levels[level].nDownIts = nIts;
	self->solversChanged = True;
}

void MultigridSolver_SetLevelUpSolver( void* matrixSolver, unsigned levelInd, void* solver ) {
	MultigridSolver*	self = (MultigridSolver*)matrixSolver;
	MultigridSolver_Level*	level;

	assert( self && Stg_CheckType( self, MultigridSolver ) );
	assert( levelInd < self->nLevels && levelInd > 0 );

	level = self->levels + levelInd;
	if( level->downSolver != solver )
		self->solversChanged = True;
	if( level->downSolver )
		Stg_Class_RemoveRef( level->downSolver );
	level->downSolver = solver;
	if( solver )
		Stg_Class_AddRef( solver );
}

void MultigridSolver_SetLevelUpIterations( void* matrixSolver, unsigned level, unsigned nIts ) {
	MultigridSolver*	self = (MultigridSolver*)matrixSolver;

	assert( self && Stg_CheckType( self, MultigridSolver ) );
	assert( level < self->nLevels );
	assert( level > 0 );

	self->levels[level].nUpIts = nIts;
	self->solversChanged = True;
}

void MultigridSolver_SetLevelCycles( void* matrixSolver, unsigned level, unsigned nCycles ) {
	MultigridSolver*	self = (MultigridSolver*)matrixSolver;

	assert( self && Stg_CheckType( self, MultigridSolver ) );
	assert( level < self->nLevels );
	assert( level > 0 );

	self->levels[level].nCycles = nCycles;
	self->solversChanged = True;
}

void MultigridSolver_SetAllDownSolver( void* matrixSolver, void* solver ) {
	MultigridSolver*	self = (MultigridSolver*)matrixSolver;
	unsigned		l_i;

	assert( self && Stg_CheckType( self, MultigridSolver ) );

	for( l_i = 1; l_i < self->nLevels; l_i++ )
		MultigridSolver_SetLevelDownSolver( self, l_i, solver );
}

void MultigridSolver_SetAllDownIterations( void* matrixSolver, unsigned level, unsigned nIts ) {
	MultigridSolver*	self = (MultigridSolver*)matrixSolver;
	unsigned		l_i;

	assert( self && Stg_CheckType( self, MultigridSolver ) );

	for( l_i = 1; l_i < self->nLevels; l_i++ )
		MultigridSolver_SetLevelDownIterations( self, l_i, nIts );
}

void MultigridSolver_SetAllUpSolver( void* matrixSolver, void* solver ) {
	MultigridSolver*	self = (MultigridSolver*)matrixSolver;
	unsigned		l_i;

	assert( self && Stg_CheckType( self, MultigridSolver ) );

	for( l_i = 1; l_i < self->nLevels; l_i++ )
		MultigridSolver_SetLevelUpSolver( self, l_i, solver );
}

void MultigridSolver_SetAllUpIterations( void* matrixSolver, unsigned level, unsigned nIts ) {
	MultigridSolver*	self = (MultigridSolver*)matrixSolver;
	unsigned		l_i;

	assert( self && Stg_CheckType( self, MultigridSolver ) );

	for( l_i = 1; l_i < self->nLevels; l_i++ )
		MultigridSolver_SetLevelUpIterations( self, l_i, nIts );
}

void MultigridSolver_SetAllSolver( void* matrixSolver, void* solver ) {
	MultigridSolver*	self = (MultigridSolver*)matrixSolver;

	assert( self );

	MultigridSolver_SetAllDownSolver( self, solver );
	MultigridSolver_SetAllUpSolver( self, solver );
}

void MultigridSolver_SetCoarseSolver( void* matrixSolver, void* solver ) {
	MultigridSolver*	self = (MultigridSolver*)matrixSolver;

	assert( self );

	MultigridSolver_SetLevelDownSolver( self, 0, solver );
}

unsigned MultigridSolver_GetNumLevels( void* matrixSolver ) {
	MultigridSolver*	self = (MultigridSolver*)matrixSolver;

	assert( self && Stg_CheckType( self, MultigridSolver ) );

	return self->nLevels;
}

Matrix* MultigridSolver_GetRestriction( void* matrixSolver, unsigned level ) {
	MultigridSolver*	self = (MultigridSolver*)matrixSolver;

	assert( self && Stg_CheckType( self, MultigridSolver ) );
	assert( level < self->nLevels && level > 0 );

	return self->levels[level].R;
}

Matrix* MultigridSolver_GetProlongation( void* matrixSolver, unsigned level ) {
	MultigridSolver*	self = (MultigridSolver*)matrixSolver;

	assert( self && Stg_CheckType( self, MultigridSolver ) );
	assert( level < self->nLevels && level > 0 );

	return self->levels[level].P;
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Private Functions
*/

void MultigridSolver_RestrictMatrix( MultigridSolver* self, MultigridSolver_Level* level, Matrix** dstMatrix ) {
	Matrix		*srcMat;
	unsigned	nSrcRows, nSrcCols;
	unsigned	nRRows, nRCols;

	assert( self );
	assert( level );
	assert( dstMatrix );

	srcMat = MatrixSolver_GetMatrix( level->downSolver );
	Matrix_GetGlobalSize( srcMat, &nSrcRows, &nSrcCols );
	Matrix_GetGlobalSize( level->R, &nRRows, &nRCols );
	if( nRCols == nSrcRows )
		Matrix_PAPt( srcMat, level->R, (void**)dstMatrix );
	else
		Matrix_PtAP( srcMat, level->R, (void**)dstMatrix );
}

void MultigridSolver_LevelCycle( MultigridSolver* self, unsigned levelInd, Vector* rhs, Vector* solution ) {
	MultigridSolver_Level*	level;

	assert( self );
	assert( levelInd < self->nLevels );

	Stream_Indent( self->stream );
	level = self->levels + levelInd;

	if( level->nDownIts ) {
		Journal_Printf( self->stream, "Down-solve on level %d... ", levelInd );
		MatrixSolver_Solve( level->downSolver, rhs, solution );
		Journal_Printf( self->stream, "residual: %.10lf\n", MatrixSolver_GetResidualNorm( level->downSolver ) );
	}

	if( levelInd > 0 ) {
		MultigridSolver_Level*	nextLevel;
		unsigned		c_i;

		nextLevel = self->levels + (levelInd - 1);
		if( level->R == level->P )
			Matrix_TransposeMultiply( level->R, MatrixSolver_GetResidual( level->downSolver ), nextLevel->workRHS );
		else
			Matrix_Multiply( level->R, MatrixSolver_GetResidual( level->downSolver ), nextLevel->workRHS );
		for( c_i = 0; c_i < level->nCycles; c_i++ )
			MultigridSolver_LevelCycle( self, levelInd - 1, nextLevel->workRHS, nextLevel->workSol );
		Matrix_MultiplyAdd( level->P, nextLevel->workSol, solution, solution );
	}

	if( level->nUpIts ) {
		Journal_Printf( self->stream, "Up-solve on level %d... ", levelInd );
		MatrixSolver_Solve( level->upSolver, rhs, solution );
		Journal_Printf( self->stream, "residual: %.10lf\n", MatrixSolver_GetResidualNorm( level->upSolver ) );
	}

	Stream_UnIndent( self->stream );
}

void MultigridSolver_UpdateSolvers( MultigridSolver* self ) {
	MultigridSolver_Level*	level;
	unsigned		l_i;

	assert( self && Stg_CheckType( self, MultigridSolver ) );

	for( l_i = self->nLevels - 1; l_i < self->nLevels; l_i-- ) {
		level = self->levels + l_i;

		if( l_i > 0 ) {
			if( MatrixSolver_GetMaxIterations( level->downSolver ) != level->nDownIts )
				MatrixSolver_SetMaxIterations( level->downSolver, level->nDownIts );
			if( MatrixSolver_GetMaxIterations( level->upSolver ) != level->nUpIts )
				MatrixSolver_SetMaxIterations( level->upSolver, level->nUpIts );

			if( l_i == self->nLevels - 1 )
				MatrixSolver_SetUseInitialSolution( level->downSolver, True );
			else
				MatrixSolver_SetUseInitialSolution( level->downSolver, False );
			MatrixSolver_SetUseInitialSolution( level->upSolver, True );
		}
	}

	FreeObject( self->outerSolver );
	self->outerSolver = MultigridSolver_CreateOuterSolver( self, self->matrix );
	MatrixSolver_SetRelativeTolerance( self->outerSolver, self->relTol );
}

void MultigridSolver_UpdateMatrices( MultigridSolver* self ) {
	MultigridSolver_Level*	level;
	unsigned		l_i;

	assert( self && Stg_CheckType( self, MultigridSolver ) );

	for( l_i = self->nLevels - 1; l_i < self->nLevels; l_i-- ) {
		level = self->levels + l_i;

		Journal_Printf( self->stream, "Updating matrix on level %d...\n", l_i );

		if( l_i < self->nLevels - 1 ) {
			Matrix*	mat;

			if( level->downSolver )
				mat = MatrixSolver_GetMatrix( level->downSolver );
			else
				mat = NULL;
			if( mat ) {
				KillObject( mat );
			}
			MultigridSolver_RestrictMatrix( self, self->levels + l_i + 1, &mat );
			level->A = mat;
			if( level->downSolver )
				MatrixSolver_SetMatrix( level->downSolver, mat );
			else {
				if( l_i > 0 )
					level->downSolver = MultigridSolver_CreateSmoother( self, mat );
				else
					level->downSolver = MultigridSolver_CreateCoarseSolver( self, mat );
			}
			if( l_i > 0 ) {
				if( level->upSolver )
					MatrixSolver_SetMatrix( level->upSolver, mat );
				else
					level->upSolver = MultigridSolver_CreateSmoother( self, mat );
			}
		}
		else {
			level->A = self->matrix;
			if( level->downSolver )
				MatrixSolver_SetMatrix( level->downSolver, self->matrix );
			else {
				if( l_i > 0 )
					level->downSolver = MultigridSolver_CreateSmoother( self, self->matrix );
				else
					level->downSolver = MultigridSolver_CreateCoarseSolver( self, self->matrix );
			}
			if( l_i > 0 ) {
				if( level->upSolver )
					MatrixSolver_SetMatrix( level->upSolver, self->matrix );
				else
					level->upSolver = MultigridSolver_CreateSmoother( self, self->matrix );
			}
		}

		Journal_Printf( self->stream, "done\n" );
	}

	MultigridSolver_UpdateSolvers( self );
	MatrixSolver_SetMatrix( self->outerSolver, self->matrix );
}

void MultigridSolver_UpdateOps( MultigridSolver* self ) {
	MultigridSolver_Level*	level;
	Matrix			**pOps, **rOps;
	unsigned		l_i;

	assert( self && Stg_CheckType( self, MultigridSolver ) );

	MGOpGenerator_Generate( self->opGen, &pOps, &rOps );
	for( l_i = 1; l_i < self->nLevels; l_i++ ) {
		level = self->levels + l_i;

		if( !level->P ) {
			level->P = pOps[l_i];
			Stg_Class_AddRef( level->P );
		}
		else
			Stg_Class_RemoveRef( pOps[l_i] );

		if( !level->R ) {
			level->R = rOps[l_i];
			Stg_Class_AddRef( level->R );
		}
		else
			Stg_Class_RemoveRef( rOps[l_i] );
	}

	FreeArray( pOps );
	FreeArray( rOps );
}

void MultigridSolver_UpdateWorkVectors( MultigridSolver* self ) {
	MultigridSolver_Level*	level;
	unsigned		rowSize, colSize;
	unsigned		l_i;

	assert( self && Stg_CheckType( self, MultigridSolver ) );

	for( l_i = 0; l_i < self->nLevels - 1; l_i++ ) {
		level = self->levels + l_i;

		Matrix_GetLocalSize( MatrixSolver_GetMatrix( level->downSolver ), &rowSize, &colSize );

		if( !level->workSol || Vector_GetLocalSize( level->workSol ) != rowSize ) {
			if( level->workSol )
				Stg_Class_RemoveRef( level->workSol );
			Vector_Duplicate( self->curSolution, (void**)&level->workSol );
			Vector_SetLocalSize( level->workSol, rowSize );
		}

		if( !level->workRHS || Vector_GetLocalSize( level->workRHS ) != rowSize ) {
			if( level->workRHS )
				Stg_Class_RemoveRef( level->workRHS );
			Vector_Duplicate( self->curSolution, (void**)&level->workRHS );
			Vector_SetLocalSize( level->workRHS, rowSize );
		}
	}
}

MatrixSolver* MultigridSolver_CreateOuterSolver( MultigridSolver* self, Matrix* matrix ) {
	MatrixSolver*	outerSolver;

#ifdef HAVE_PETSC
	outerSolver = (MatrixSolver*)PETScMatrixSolver_New( "" );
	PETScMatrixSolver_SetKSPType( outerSolver, PETScMatrixSolver_KSPType_Richardson );
	PETScMatrixSolver_SetPCType( outerSolver, PETScMatrixSolver_PCType_SOR );
	PETScMatrixSolver_SetMatrix( outerSolver, matrix );
	PETScMatrixSolver_SetMaxIterations( outerSolver, 3 );
	PETScMatrixSolver_SetUseInitialSolution( outerSolver, True );
	PETScMatrixSolver_SetNormType( outerSolver, PETScMatrixSolver_NormType_Preconditioned );
#else
#error Oi! Where is PETSc?
#endif

	return outerSolver;
}

MatrixSolver* MultigridSolver_CreateSmoother( MultigridSolver* self, Matrix* matrix ) {
	MatrixSolver*	smoother;
	unsigned	nBlocks;
	KSP*		ksps;
	PC		pc;
	PetscErrorCode	ec;

#ifdef HAVE_PETSC
	smoother = (MatrixSolver*)PETScMatrixSolver_New( "" );
	PETScMatrixSolver_SetKSPType( smoother, 
				      PETScMatrixSolver_KSPType_Richardson );
	PETScMatrixSolver_SetPCType( smoother, 
				     PETScMatrixSolver_PCType_SOR );
	MatrixSolver_SetMatrix( smoother, matrix );
#else
#error *** Oi! Where is PETSc?
#endif

	return smoother;
}

MatrixSolver* MultigridSolver_CreateCoarseSolver( MultigridSolver* self, Matrix* matrix ) {
	MatrixSolver*	coarseSolver;
	unsigned	nProcs;

	MPI_Comm_size( self->comm, (int*)&nProcs );

#ifdef HAVE_PETSC
	coarseSolver = (MatrixSolver*)PETScMatrixSolver_New( "" );
	PETScMatrixSolver_SetKSPType( coarseSolver, 
				      PETScMatrixSolver_KSPType_PreOnly );
	if( nProcs == 1 ) {
		PETScMatrixSolver_SetPCType( coarseSolver, 
					     PETScMatrixSolver_PCType_LU );
	}
	else {
		PETScMatrixSolver_SetPCType( coarseSolver, 
					     PETScMatrixSolver_PCType_RedundantLU );
	}
	MatrixSolver_SetMatrix( coarseSolver, matrix );
#else
#error *** Oi! Where is PETSc?
#endif

	return coarseSolver;
}

void MultigridSolver_Destruct( MultigridSolver* self ) {
	assert( self && Stg_CheckType( self, MultigridSolver ) );

	MultigridSolver_DestructLevels( self );
	KillObject( self->opGen );
}

void MultigridSolver_DestructLevels( MultigridSolver* self ) {
	unsigned	l_i;

	assert( self && Stg_CheckType( self, MultigridSolver ) );

	for( l_i = 0; l_i < self->nLevels; l_i++ ) {
		MultigridSolver_Level*	level = self->levels + l_i;

		if( level->downSolver ) {
			Stg_Class_RemoveRef( MatrixSolver_GetMatrix( level->downSolver ) );
			Stg_Class_RemoveRef( level->downSolver );
		}
		if( level->upSolver ) {
			Stg_Class_RemoveRef( MatrixSolver_GetMatrix( level->upSolver ) );
			Stg_Class_RemoveRef( level->upSolver );
		}
		if( level->R )
			Stg_Class_RemoveRef( level->R );
		if( level->P )
			Stg_Class_RemoveRef( level->P );
		if( level->workRHS )
			Stg_Class_RemoveRef( level->workRHS );
		if( level->workSol )
			Stg_Class_RemoveRef( level->workSol );
	}

	KillArray( self->levels );
	self->nLevels = 0;
	self->solversChanged = True;
	self->opsChanged = True;

	/* Temporary. */
	KillObject( self->outerSolver );
}
