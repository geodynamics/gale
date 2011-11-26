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
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/Discretisation/Discretisation.h>

#include "SystemSetup.h"


/* Textual name of this class */
const Type MultigridSolver_Type = "MultigridSolver";



/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/

MultigridSolver* MultigridSolver_New( Name name ) {
	/* Variables set in this function */
	SizeT                                                    _sizeOfSelf = sizeof(MultigridSolver);
	Type                                                            type = MultigridSolver_Type;
	Stg_Class_DeleteFunction*                                    _delete = _MultigridSolver_Delete;
	Stg_Class_PrintFunction*                                      _print = _MultigridSolver_Print;
	Stg_Class_CopyFunction*                                        _copy = NULL;
	Stg_Component_DefaultConstructorFunction*        _defaultConstructor = (void* (*)(Name))_MultigridSolver_New;
	Stg_Component_ConstructFunction*                          _construct = _MultigridSolver_AssignFromXML;
	Stg_Component_BuildFunction*                                  _build = _MultigridSolver_Build;
	Stg_Component_InitialiseFunction*                        _initialise = _MultigridSolver_Initialise;
	Stg_Component_ExecuteFunction*                              _execute = _MultigridSolver_Execute;
	Stg_Component_DestroyFunction*                              _destroy = _MultigridSolver_Destroy;
	AllocationType                                    nameAllocationType = NON_GLOBAL;
	MGSolver_SetCommFunc*                                    setCommFunc = MultigridSolver_SetComm;
	MGSolver_SetMatrixFunc*                                setMatrixFunc = MultigridSolver_SetMatrix;
	MGSolver_SetMaxIterationsFunc*                  setMaxIterationsFunc = MultigridSolver_SetMaxIterations;
	MGSolver_SetRelativeToleranceFunc*          setRelativeToleranceFunc = MultigridSolver_SetRelativeTolerance;
	MGSolver_SetAbsoluteToleranceFunc*          setAbsoluteToleranceFunc = MultigridSolver_SetAbsoluteTolerance;
	MGSolver_SetUseInitialSolutionFunc*        setUseInitialSolutionFunc = MultigridSolver_SetUseInitialSolution;
	MGSolver_SolveFunc*                                        solveFunc = MultigridSolver_Solve;
	MGSolver_SetupFunc*                                        setupFunc = MultigridSolver_Setup;
	MGSolver_GetSolveStatusFunc*                      getSolveStatusFunc = MultigridSolver_GetSolveStatus;
	MGSolver_GetIterationsFunc*                        getIterationsFunc = MultigridSolver_GetIterations;
	MGSolver_GetMaxIterationsFunc*                  getMaxIterationsFunc = MultigridSolver_GetMaxIterations;
	MGSolver_GetResidualNormFunc*                    getResidualNormFunc = MultigridSolver_GetResidualNorm;

	return _MultigridSolver_New(  MULTIGRIDSOLVER_PASSARGS  );
}

MultigridSolver* _MultigridSolver_New(  MULTIGRIDSOLVER_DEFARGS  ) {
	MultigridSolver*	self;

	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(MultigridSolver) );

	self = (MultigridSolver*)_Stg_Component_New(  STG_COMPONENT_PASSARGS  );

	/* function assignments previously in the MatrixSolver_New func */
	self->setCommFunc = setCommFunc;
	self->setMatrixFunc = setMatrixFunc;
	self->setMaxIterationsFunc = setMaxIterationsFunc;
	self->setRelativeToleranceFunc = setRelativeToleranceFunc;
	self->setAbsoluteToleranceFunc = setAbsoluteToleranceFunc;
	self->setUseInitialSolutionFunc = setUseInitialSolutionFunc;

	self->solveFunc = solveFunc;
	self->setupFunc = setupFunc;

	self->getSolveStatusFunc = getSolveStatusFunc;
	self->getIterationsFunc = getIterationsFunc;
	self->getMaxIterationsFunc = getMaxIterationsFunc;
	self->getResidualNormFunc = getResidualNormFunc;

	/* Virtual info */

	/* MultigridSolver info */
	_MultigridSolver_Init( self );

	return self;
}

void _MultigridSolver_Init( MultigridSolver* self ) {
	assert( self && Stg_CheckType( self, MultigridSolver ) );

	/* these initialisations previously done in the MatrixSolver_Init func */
	self->mgData = (MGSolver_PETScData*)malloc( sizeof( MGSolver_PETScData ) );

	self->mgData->comm = MPI_COMM_WORLD;
	KSPCreate( MPI_COMM_WORLD, &self->mgData->ksp );
	self->mgData->matrix = PETSC_NULL;
	self->mgData->inversion = PETSC_NULL;
	self->mgData->residual = PETSC_NULL;
	self->mgData->expiredResidual = True;
	self->mgData->matrixChanged = True;

	self->mgData->curRHS = PETSC_NULL;
	self->mgData->curSolution = PETSC_NULL;
	/* end of old MatrixSolver_Init initialisations */

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
	//self->outerSolver = NULL;
	self->outerSolver = (MGSolver_PETScData*)malloc( sizeof( MGSolver_PETScData ) );
	self->outerSolver->ksp = PETSC_NULL;
	self->outerSolver->matrix = PETSC_NULL;
	self->outerSolver->inversion = PETSC_NULL;
	self->outerSolver->residual = PETSC_NULL;
	self->outerSolver->curRHS = PETSC_NULL;
	self->outerSolver->curSolution = PETSC_NULL;
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _MultigridSolver_Delete( void* matrixSolver ) {
	MultigridSolver*	self = (MultigridSolver*)matrixSolver;

	assert( self && Stg_CheckType( self, MultigridSolver ) );

	MultigridSolver_Destruct( self );

	/* Delete the parent. */
	_Stg_Component_Delete( self );
}

void _MultigridSolver_Print( void* matrixSolver, Stream* stream ) {
	MultigridSolver*	self = (MultigridSolver*)matrixSolver;
	
	assert( self && Stg_CheckType( self, MultigridSolver ) );

	/* Print parent */
	Journal_Printf( stream, "MultigridSolver (ptr): (%p)\n", self );
	_Stg_Component_Print( self, stream );
}

void _MultigridSolver_AssignFromXML( void* matrixSolver, Stg_ComponentFactory* cf, void* data ) {
	MultigridSolver*	self = (MultigridSolver*)matrixSolver;
	unsigned		nLevels;

	assert( self && Stg_CheckType( self, MultigridSolver ) );
	assert( cf );

	nLevels = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, (Dictionary_Entry_Key)"levels", 1  );
	MultigridSolver_SetLevels( self, nLevels );
	self->opGen = Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"opGenerator", MGOpGenerator, True, data  );
	MGOpGenerator_SetMatrixSolver( self->opGen, self );
	MGOpGenerator_SetNumLevels( self->opGen, nLevels );
}

void _MultigridSolver_Build( void* matrixSolver, void* data ) {
	MultigridSolver*	self = (MultigridSolver*)matrixSolver;

	assert( self && Stg_CheckType( self, MultigridSolver ) );

	self->stream = Journal_Register( InfoStream_Type, (Name)"general" );
	if( self->opGen  )
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

/* copied from the PETScMatrixSolver class, has been depreciated */
void MultigridSolver_SetComm( void* matrixSolver, MPI_Comm comm ) {
	MultigridSolver*	self = (MultigridSolver*)matrixSolver;

	self->mgData->comm = comm;

	if( self->mgData->ksp != PETSC_NULL )
		KSPDestroy( self->mgData->ksp );
	KSPCreate( comm, &self->mgData->ksp );
}

/* copied from the PETScMatrixSolver class, has been depreciated */
void MultigridSolver_SetMatrix( void* matrixSolver, void* _matrix ) {
	MultigridSolver*	self 	= (MultigridSolver*)matrixSolver;
	Mat			matrix	= (Mat)_matrix;

	self->mgData->matrix = matrix;
	KSPSetOperators( self->mgData->ksp, matrix, matrix, DIFFERENT_NONZERO_PATTERN );
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

	//MatrixSolver_SetAbsoluteTolerance( self->outerSolver, tolerance );
	KSPSetTolerances( self->outerSolver->ksp, PETSC_DEFAULT, tolerance, PETSC_DEFAULT, PETSC_DEFAULT );

	self->solversChanged = True;
}

void MultigridSolver_SetUseInitialSolution( void* matrixSolver, Bool state ) {
	MultigridSolver*	self = (MultigridSolver*)matrixSolver;

	assert( self && Stg_CheckType( self, MultigridSolver ) );

	self->useInitial = state;
}

/* these functions were associated with the PETScMatrixSolver class, before it was
 * depreciated */
Vec _GetResidual( MGSolver_PETScData* mgData ) {
	if( mgData->expiredResidual ) {
		VecDuplicate( mgData->curSolution, &mgData->residual );	
		VecSetFromOptions( mgData->residual );
#if( PETSC_VERSION_MAJOR <= 2 && PETSC_VERSION_MINOR >= 3 && PETSC_VERSION_SUBMINOR >= 3 )
		VecSetOption( mgData->residual, VEC_IGNORE_NEGATIVE_INDICES );
#elif( PETSC_VERSION_MAJOR >= 3 )
		VecSetOption( mgData->residual, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE );
#endif
		MatMult( mgData->matrix, mgData->curSolution, mgData->residual );
		VecAYPX( mgData->residual, -1.0, mgData->curRHS );

		mgData->expiredResidual = False;
	}

	return mgData->residual;
}

double _GetResidualNorm( MGSolver_PETScData* mgData ) {
	PC			pc;
	const KSPType		kspType;
	const PCType		pcType;
	PetscScalar		rnorm;
	PetscErrorCode		ec;

	ec = KSPGetType( mgData->ksp, &kspType );
	CheckPETScError( ec );
	ec = KSPGetPC( mgData->ksp, &pc );
	CheckPETScError( ec );
	ec = PCGetType( pc, &pcType );
	CheckPETScError( ec );

	if( !strcmp( kspType, KSPRICHARDSON ) && !strcmp( pcType, PCSOR ) ) {
		Vec	residual;

		//residual = MatrixSolver_GetResidual( mgData );
		//rnorm = (PetscScalar)Vector_L2Norm( residual );
		residual = _GetResidual( mgData );
		VecNorm( residual, NORM_2, &rnorm );
	}
	else {
		ec = KSPGetResidualNorm( mgData->ksp, &rnorm );
		CheckPETScError( ec );
	}

	return (double)rnorm;
}

#define MultigridSolver_NormType_Preconditioned 1

MGSolver_Status _GetSolveStatus( MGSolver_PETScData* mgData ) {
	PC			pc;
	const KSPType		kspType;
	const PCType		pcType;
	KSPConvergedReason	reason;
	PetscErrorCode		ec;

	ec = KSPGetType( mgData->ksp, &kspType );
	CheckPETScError( ec );
	ec = KSPGetPC( mgData->ksp, &pc );
	CheckPETScError( ec );
	ec = PCGetType( pc, &pcType );
	CheckPETScError( ec );

	if( !strcmp( kspType, KSPRICHARDSON ) && !strcmp( pcType, PCSOR ) ) {
		double		rnorm;
		PetscInt	curIt;

		//rnorm = PETScMatrixSolver_GetResidualNorm( self );
		//curIt = PETScMatrixSolver_GetIterations( self );
		rnorm = _GetResidualNorm( mgData );
		KSPGetIterationNumber( mgData->ksp, &curIt );
		//PETScMatrixSolver_SetNormType( self, MultigridSolver_NormType_Preconditioned );
		KSPSetNormType( mgData->ksp, (KSPNormType)MultigridSolver_NormType_Preconditioned );
		ec = KSPDefaultConverged( mgData->ksp, curIt, (PetscScalar)rnorm, &reason, PETSC_NULL );
		CheckPETScError( ec );
	}
	else {
		ec = KSPGetConvergedReason( mgData->ksp, &reason );
		CheckPETScError( ec );
	}

	return (MGSolver_Status)reason;
}

/* end of functions formerly defined in the PETScMatrixSolver class */

void MultigridSolver_Solve( void* matrixSolver, void* _rhs, void* _solution ) {
	MultigridSolver*	self = (MultigridSolver*)matrixSolver;
	Vec			rhs = (Vec)_rhs;
	Vec			solution = (Vec)_solution;
	double wallTime;

	assert( self && Stg_CheckType( self, MultigridSolver ) );
	//assert( rhs && Stg_CheckType( rhs, Vector ) );
	//assert( solution && Stg_CheckType( solution, Vector ) );

	Journal_Printf( self->stream, "MultigridSolver: Starting solve ...\n" );
	Stream_Indent( self->stream );

	if( !self->useInitial ) {
		Journal_Printf( self->stream, "Zeroing initial solution\n" );
		//Vector_Zero( solution );
		VecSet( solution, 0.0 );
	}
	else
		Journal_Printf( self->stream, "Keeping initial solution\n" );

	wallTime = MPI_Wtime();
	MultigridSolver_Setup( self, rhs, solution );
	stg_profile_Func( "MultigridSolver_Setup", MPI_Wtime() - wallTime);
	self->curIt = 0;
	//while( MatrixSolver_GetSolveStatus( self ) == MultigridSolver_Status_Iterating && self->curIt < self->maxIts ) {
	while( _GetSolveStatus( self->mgData ) == MGSolver_Status_Iterating && self->curIt < self->maxIts ) {
		Journal_Printf( self->stream, "Iteration %d: residual %.10lf\n", 
				//self->curIt, MatrixSolver_GetResidualNorm( self->outerSolver ) );
				self->curIt, _GetResidualNorm( self->outerSolver ) );
		MultigridSolver_LevelCycle( self, self->nLevels - 1, rhs, solution );
		self->curIt++;
	}
	Journal_Printf( self->stream, "Iteration %d: residual %.10lf\n", 
			//self->curIt, MatrixSolver_GetResidualNorm( self->outerSolver ) );
			self->curIt, _GetResidualNorm( self->outerSolver ) );

	Stream_UnIndent( self->stream );
	Journal_Printf( self->stream, "done.\n" );
}

void MultigridSolver_Setup( void* matrixSolver, void* rhs, void* solution ) {
	MultigridSolver*	self = (MultigridSolver*)matrixSolver;
	Bool			rebuildOps;
	double wallTime;

	assert( self && Stg_CheckType( self, MultigridSolver ) );

	//_MatrixSolver_Setup( self, rhs, solution );
	self->mgData->curRHS = (struct _p_Vec*)rhs;
	self->mgData->curSolution = (struct _p_Vec*)solution;
	self->mgData->expiredResidual = True;

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
		wallTime = MPI_Wtime();
		MultigridSolver_UpdateOps( self );
		stg_profile_Func( "MultigridSolver_UpdateOps", MPI_Wtime() - wallTime);
	//if( self->matrixChanged || rebuildOps || self->opsChanged || self->solversChanged ) {
	if( self->mgData->matrixChanged || rebuildOps || self->opsChanged || self->solversChanged ) {
		wallTime = MPI_Wtime();
		MultigridSolver_UpdateMatrices( self );
		stg_profile_Func( "MultigridSolver_UpdateMatrices", MPI_Wtime() - wallTime);
		
		wallTime = MPI_Wtime();
		MultigridSolver_UpdateWorkVectors( self );
		stg_profile_Func( "MultigridSolver_UpdateWorkVectors", MPI_Wtime() - wallTime);


	}
	self->solversChanged = False;
	//self->matrixChanged = False;
	self->mgData->matrixChanged = False;
	self->opsChanged = False;
}

MGSolver_Status MultigridSolver_GetSolveStatus( void* matrixSolver ) {
	MultigridSolver*	self = (MultigridSolver*)matrixSolver;
	MGSolver_Status		outerStatus;

	assert( self && Stg_CheckType( self, MultigridSolver ) );

	//MatrixSolver_Solve( self->outerSolver, self->curRHS, self->curSolution );
	//outerStatus = MatrixSolver_GetSolveStatus( self->outerSolver );
	KSPSolve( self->outerSolver->ksp, self->mgData->curRHS, self->mgData->curSolution );
	outerStatus = _GetSolveStatus( self->outerSolver );

	if( outerStatus == MGSolver_Status_ConvergedIterations || 
	    outerStatus == MGSolver_Status_DivergedIterations || 
	    outerStatus == MGSolver_Status_Iterating )
	{
		return MGSolver_Status_Iterating;
	}
	else if( outerStatus == MGSolver_Status_ConvergedRelative )
		return MGSolver_Status_ConvergedRelative;
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

	return _GetResidualNorm( self->outerSolver );
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
	//MPI_Comm_size( self->comm, (int*)&nProcs );
	MPI_Comm_size( self->mgData->comm, (int*)&nProcs );

	for( l_i = 0; l_i < nLevels; l_i++ ) {
		MultigridSolver_Level*	level = self->levels + l_i;

		//level->downSolver = NULL;
		level->downSolver = (MGSolver_PETScData*)malloc( sizeof( MGSolver_PETScData ) );
		level->nDownIts = 1;
		//level->upSolver = NULL;
		level->upSolver = (MGSolver_PETScData*)malloc( sizeof( MGSolver_PETScData ) );
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
	//Matrix*			R = (Matrix*)_R;
	Mat			R = (Mat)_R;

	assert( self && Stg_CheckType( self, MultigridSolver ) );
	assert( levelInd < self->nLevels && levelInd > 0 );
	//assert( !R || Stg_CheckType( R, Matrix ) );

	level = self->levels + levelInd;
	if( level->R != R )
		self->opsChanged = True;
	//if( level->R )
	//	Stg_Class_RemoveRef( level->R );
	level->R = R;
	//if( R )
	//	Stg_Class_AddRef( R );
}

void MultigridSolver_SetProlongation( void* matrixSolver, unsigned levelInd, void* PP ) {
	MultigridSolver*	self 	= (MultigridSolver*)matrixSolver;
	MultigridSolver_Level*	level;
	Mat			P	= (Mat)PP;

	assert( self && Stg_CheckType( self, MultigridSolver ) );
	assert( levelInd < self->nLevels && levelInd > 0 );
	//assert( !P || Stg_CheckType( P, Matrix ) );

	level = self->levels + levelInd;
	if( level->P != P )
		self->opsChanged = True;
	//if( level->P )
	//	Stg_Class_RemoveRef( level->P );
	if( level->P != PETSC_NULL )
		MatDestroy( level->P );
	level->P = P;
	//if( P )
	//	Stg_Class_AddRef( P );
}

void MultigridSolver_SetLevelDownSolver( void* matrixSolver, unsigned levelInd, void* solver ) {
	MultigridSolver*	self = (MultigridSolver*)matrixSolver;
	MultigridSolver_Level*	level;

	assert( self && Stg_CheckType( self, MultigridSolver ) );
	assert( levelInd < self->nLevels );
	//assert( !solver || Stg_CheckType( solver, MatrixSolver ) );

	level = self->levels + levelInd;
	if( level->downSolver != solver )
		self->solversChanged = True;
	if( level->downSolver ) {
	//	Stg_Class_RemoveRef( level->downSolver );
		if( level->downSolver->ksp != PETSC_NULL )         KSPDestroy( level->downSolver->ksp );
		if( level->downSolver->matrix != PETSC_NULL )      MatDestroy( level->downSolver->matrix );
		if( level->downSolver->inversion != PETSC_NULL )   MatDestroy( level->downSolver->inversion );
		if( level->downSolver->residual != PETSC_NULL )    VecDestroy( level->downSolver->residual );
		if( level->downSolver->curRHS != PETSC_NULL )      VecDestroy( level->downSolver->curRHS );
		if( level->downSolver->curSolution != PETSC_NULL ) VecDestroy( level->downSolver->curSolution );
		free( level->downSolver );
	}
	level->downSolver = (MGSolver_PETScData*)solver;
	//if( solver )
	//	Stg_Class_AddRef( solver );
}

void MultigridSolver_SetLevelDownIterations( void* matrixSolver, unsigned level, unsigned nIts ) {
	MultigridSolver*	self = (MultigridSolver*)matrixSolver;

	assert( self && Stg_CheckType( self, MultigridSolver ) );
	assert( level < self->nLevels );

	self->levels[level].nDownIts = nIts;
	self->solversChanged = True;
}

/* this function was applying to the down solver before - pretty sure this is a mistake */
void MultigridSolver_SetLevelUpSolver( void* matrixSolver, unsigned levelInd, void* solver ) {
	MultigridSolver*	self = (MultigridSolver*)matrixSolver;
	MultigridSolver_Level*	level;

	assert( self && Stg_CheckType( self, MultigridSolver ) );
	assert( levelInd < self->nLevels && levelInd > 0 );

	level = self->levels + levelInd;
	if( level->upSolver != solver )
		self->solversChanged = True;
	if( level->upSolver ) {
	//	Stg_Class_RemoveRef( level->downSolver );
		if( level->upSolver->ksp != PETSC_NULL )         KSPDestroy( level->upSolver->ksp );
		if( level->upSolver->matrix != PETSC_NULL )      MatDestroy( level->upSolver->matrix );
		if( level->upSolver->inversion != PETSC_NULL )   MatDestroy( level->upSolver->inversion );
		if( level->upSolver->residual != PETSC_NULL )    VecDestroy( level->upSolver->residual );
		if( level->upSolver->curRHS != PETSC_NULL )      VecDestroy( level->upSolver->curRHS );
		if( level->upSolver->curSolution != PETSC_NULL ) VecDestroy( level->upSolver->curSolution );
		free( level->upSolver );
	}
	level->upSolver = (MGSolver_PETScData*)solver;
	//if( solver )
	//	Stg_Class_AddRef( solver );
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

Mat MultigridSolver_GetRestriction( void* matrixSolver, unsigned level ) {
	MultigridSolver*	self = (MultigridSolver*)matrixSolver;

	assert( self && Stg_CheckType( self, MultigridSolver ) );
	assert( level < self->nLevels && level > 0 );

	return self->levels[level].R;
}

Mat MultigridSolver_GetProlongation( void* matrixSolver, unsigned level ) {
	MultigridSolver*	self = (MultigridSolver*)matrixSolver;

	assert( self && Stg_CheckType( self, MultigridSolver ) );
	assert( level < self->nLevels && level > 0 );

	return self->levels[level].P;
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Private Functions
*/

void MultigridSolver_RestrictMatrix( MultigridSolver* self, MultigridSolver_Level* level, Mat* dstMatrix ) {
	Mat		srcMat;
	PetscInt	nSrcRows, nSrcCols;
	PetscInt	nRRows, nRCols;

#if( PETSC_VERSION_MAJOR > 2 )
	PetscInt dummy;
#endif

	MatInfo		mInfo;
	PetscInt	nRowsA, nRowsC;

	assert( self );
	assert( level );
	assert( dstMatrix );

	srcMat = level->downSolver->matrix;
	MatGetSize( srcMat, &nSrcRows, &nSrcCols );
	MatGetSize( level->R, &nRRows, &nRCols );
	if( nRCols == nSrcRows ) {
		//Matrix_PAPt( srcMat, level->R, (void**)dstMatrix );
		/* this function implementation was commented out in the PETScMatrix function - not sure if it should 
		 * be implemented here?? */
		MatGetInfo( srcMat, MAT_GLOBAL_SUM, &mInfo );
#if( PETSC_VERSION_MAJOR <= 2 )
		nRowsA = mInfo.rows_global;
#else
		MatGetSize( srcMat, &nRowsA, &dummy );
#endif
                // PetscScalar nzA = mInfo.nz_used;
		MatGetInfo( level->R, MAT_GLOBAL_SUM, &mInfo );
#if( PETSC_VERSION_MAJOR <= 2 )
		nRowsC = mInfo.columns_global;
#else
		MatGetSize( level->R, &dummy, &nRowsC );
#endif
                // PetscScalar nzP = mInfo.nz_used;
                // PetscScalar nzC = ( nRowsC / nRowsA ) * nzA;
		// PetscScalar fillRatio = nzC / ( nzA + nzP );

		/* this function doesn't exist! */
		//MatPAPt( srcMat, level->R, MAT_REUSE_MATRIX, fillRatio, dstMatrix );
	}
	else {
		//Matrix_PtAP( srcMat, level->R, (void**)dstMatrix );
		if( dstMatrix )
			MatPtAP( srcMat, level->R, MAT_REUSE_MATRIX, 1.0, dstMatrix );
		else 
			MatPtAP( srcMat, level->R, MAT_INITIAL_MATRIX, 1.0, dstMatrix );
	}
}

//void MultigridSolver_LevelCycle( MultigridSolver* self, unsigned levelInd, Vector* rhs, Vector* solution ) {
void MultigridSolver_LevelCycle( MultigridSolver* self, unsigned levelInd, Vec rhs, Vec solution ) {
	MultigridSolver_Level*	level;

	assert( self );
	assert( levelInd < self->nLevels );

	Stream_Indent( self->stream );
	level = self->levels + levelInd;

	if( level->nDownIts ) {
		Journal_Printf( self->stream, "Down-solve on level %d... ", levelInd );
		//MatrixSolver_Solve( level->downSolver, rhs, solution );
		//Journal_Printf( self->stream, "residual: %.10lf\n", MatrixSolver_GetResidualNorm( level->downSolver ) );
		KSPSolve( level->downSolver->ksp, rhs, solution );
		Journal_Printf( self->stream, "residual: %.10lf\n", _GetResidualNorm( level->downSolver ) );
	}

	if( levelInd > 0 ) {
		MultigridSolver_Level*	nextLevel;
		unsigned		c_i;

		nextLevel = self->levels + (levelInd - 1);
		if( level->R == level->P )
			//Matrix_TransposeMultiply( level->R, MatrixSolver_GetResidual( level->downSolver ), nextLevel->workRHS );
			MatMultTranspose( level->R, _GetResidual( level->downSolver ), nextLevel->workRHS );
		else
			//Matrix_Multiply( level->R, MatrixSolver_GetResidual( level->downSolver ), nextLevel->workRHS );
			MatMult( level->R, _GetResidual( level->downSolver ), nextLevel->workRHS );
		for( c_i = 0; c_i < level->nCycles; c_i++ )
			MultigridSolver_LevelCycle( self, levelInd - 1, nextLevel->workRHS, nextLevel->workSol );
		//Matrix_MultiplyAdd( level->P, nextLevel->workSol, solution, solution );
		MatMultAdd( level->P, nextLevel->workSol, solution, solution );
	}

	if( level->nUpIts ) {
		Journal_Printf( self->stream, "Up-solve on level %d... ", levelInd );
		//MatrixSolver_Solve( level->upSolver, rhs, solution );
		//Journal_Printf( self->stream, "residual: %.10lf\n", MatrixSolver_GetResidualNorm( level->upSolver ) );
		KSPSolve( level->upSolver->ksp, rhs, solution );
		Journal_Printf( self->stream, "residual: %.10lf\n", _GetResidualNorm( level->upSolver ) );
	}

	Stream_UnIndent( self->stream );
}

void MultigridSolver_UpdateSolvers( MultigridSolver* self ) {
	MultigridSolver_Level*	level;
	unsigned		l_i;
	PetscInt		nDownIts, nUpIts;

	assert( self && Stg_CheckType( self, MultigridSolver ) );

	for( l_i = self->nLevels - 1; l_i < self->nLevels; l_i-- ) {
		level = self->levels + l_i;

		if( l_i > 0 ) {
			KSPGetTolerances( level->downSolver->ksp, PETSC_NULL, PETSC_NULL, PETSC_NULL, &nDownIts );
			KSPGetTolerances( level->upSolver->ksp, PETSC_NULL, PETSC_NULL, PETSC_NULL, &nUpIts );
			if( /*MatrixSolver_GetMaxIterations( level->downSolver )*/nDownIts != level->nDownIts )
				//MatrixSolver_SetMaxIterations( level->downSolver, level->nDownIts );
				KSPSetTolerances( level->downSolver->ksp, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, level->nDownIts );
			if( /*MatrixSolver_GetMaxIterations( level->upSolver )*/nUpIts != level->nUpIts )
				//MatrixSolver_SetMaxIterations( level->upSolver, level->nUpIts );
				KSPSetTolerances( level->upSolver->ksp, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, level->nUpIts );

			if( l_i == self->nLevels - 1 )
				//MatrixSolver_SetUseInitialSolution( level->downSolver, True );
				KSPSetInitialGuessNonzero( level->downSolver->ksp, (PetscTruth)True );
			else
				//MatrixSolver_SetUseInitialSolution( level->downSolver, False );
				KSPSetInitialGuessNonzero( level->downSolver->ksp, (PetscTruth)False );
			//MatrixSolver_SetUseInitialSolution( level->upSolver, True );
			KSPSetInitialGuessNonzero( level->upSolver->ksp, (PetscTruth)True );
		}
	}

	//FreeObject( self->outerSolver );
	self->outerSolver = MultigridSolver_CreateOuterSolver( self, self->mgData->matrix );
	//MatrixSolver_SetRelativeTolerance( self->outerSolver, self->relTol );
	KSPSetTolerances( self->outerSolver->ksp, self->relTol, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT );
}

void MultigridSolver_UpdateMatrices( MultigridSolver* self ) {
	MultigridSolver_Level*	level;
	unsigned		l_i;

	assert( self && Stg_CheckType( self, MultigridSolver ) );

	for( l_i = self->nLevels - 1; l_i < self->nLevels; l_i-- ) {
		level = self->levels + l_i;

		Journal_Printf( self->stream, "Updating matrix on level %d...\n", l_i );

		if( l_i < self->nLevels - 1 ) {
			//Matrix*	mat;
			Mat	mat;

			if( level->downSolver )
				//mat = MatrixSolver_GetMatrix( level->downSolver );
				mat = level->downSolver->matrix; 
			else
				mat = NULL;
			//if( mat ) {
			if( mat != PETSC_NULL ) {
				//KillObject( mat );
				MatDestroy( mat );
			}
			MultigridSolver_RestrictMatrix( self, self->levels + l_i + 1, &mat );
			level->A = mat;
			if( level->downSolver ) {
				//MatrixSolver_SetMatrix( level->downSolver, mat );
				level->downSolver->matrix = mat;
				KSPSetOperators( level->downSolver->ksp, mat, mat, DIFFERENT_NONZERO_PATTERN );
			}
			else {
				if( l_i > 0 )
					level->downSolver = MultigridSolver_CreateSmoother( self, mat );
				else
					level->downSolver = MultigridSolver_CreateCoarseSolver( self, mat );
			}
			if( l_i > 0 ) {
				if( level->upSolver ) {
					//MatrixSolver_SetMatrix( level->upSolver, mat );
					level->upSolver->matrix = mat;
					KSPSetOperators( level->upSolver->ksp, mat, mat, DIFFERENT_NONZERO_PATTERN );
				}
				else
					level->upSolver = MultigridSolver_CreateSmoother( self, mat );
			}
		}
		else {
			level->A = self->mgData->matrix;
			if( level->downSolver ) {
				//MatrixSolver_SetMatrix( level->downSolver, self->matrix );
				level->downSolver->matrix = self->mgData->matrix;
				KSPSetOperators( level->downSolver->ksp, self->mgData->matrix, self->mgData->matrix, DIFFERENT_NONZERO_PATTERN );
			}
			else {
				if( l_i > 0 )
					level->downSolver = MultigridSolver_CreateSmoother( self, self->mgData->matrix );
				else
					level->downSolver = MultigridSolver_CreateCoarseSolver( self, self->mgData->matrix );
			}
			if( l_i > 0 ) {
				if( level->upSolver ) {
					//MatrixSolver_SetMatrix( level->upSolver, self->matrix );
					level->upSolver->matrix = self->mgData->matrix;
					KSPSetOperators( level->upSolver->ksp, self->mgData->matrix, self->mgData->matrix, DIFFERENT_NONZERO_PATTERN );
				}
				else
					level->upSolver = MultigridSolver_CreateSmoother( self, self->mgData->matrix );
			}
		}

		Journal_Printf( self->stream, "done\n" );
	}

	MultigridSolver_UpdateSolvers( self );
	//MatrixSolver_SetMatrix( self->outerSolver, self->matrix );
	self->outerSolver->matrix = self->mgData->matrix;
	KSPSetOperators( self->outerSolver->ksp, self->mgData->matrix, self->mgData->matrix, DIFFERENT_NONZERO_PATTERN );
}

void MultigridSolver_UpdateOps( MultigridSolver* self ) {
	MultigridSolver_Level*	level;
	//Matrix			**pOps, **rOps;
	Mat			*pOps, *rOps;
	unsigned		l_i;

	assert( self && Stg_CheckType( self, MultigridSolver ) );

	MGOpGenerator_Generate( self->opGen, &pOps, &rOps );
	for( l_i = 1; l_i < self->nLevels; l_i++ ) {
		level = self->levels + l_i;

		//if( !level->P ) {
		if( level->P == PETSC_NULL ) {
			level->P = pOps[l_i];
			//Stg_Class_AddRef( level->P );
		}
		else
			MatDestroy( pOps[l_i] );
			//Stg_Class_RemoveRef( pOps[l_i] );


		//if( !level->R ) {
		if( level->R == PETSC_NULL ) {
			level->R = rOps[l_i];
			//Stg_Class_AddRef( level->R );
		}
		else
			MatDestroy( rOps[l_i] );
			//Stg_Class_RemoveRef( rOps[l_i] );
	}

	FreeArray( pOps );
	FreeArray( rOps );
}

void MultigridSolver_UpdateWorkVectors( MultigridSolver* self ) {
	MultigridSolver_Level*	level;
	//unsigned		rowSize, colSize;
	PetscInt		rowSize, colSize, vecSize;
	unsigned		l_i;

	assert( self && Stg_CheckType( self, MultigridSolver ) );

	for( l_i = 0; l_i < self->nLevels - 1; l_i++ ) {
		level = self->levels + l_i;

		//Matrix_GetLocalSize( MatrixSolver_GetMatrix( level->downSolver ), &rowSize, &colSize );
		MatGetLocalSize( level->downSolver->matrix, &rowSize, &colSize );

		VecGetLocalSize( level->workSol, &vecSize );

		//if( !level->workSol || Vector_GetLocalSize( level->workSol ) != rowSize ) {
		if( !level->workSol || vecSize != rowSize ) {
			//if( level->workSol )
			//	Stg_Class_RemoveRef( level->workSol );
			//Vector_Duplicate( self->curSolution, (void**)&level->workSol );
			//Vector_SetLocalSize( level->workSol, rowSize );
			if( level->workSol != PETSC_NULL )
				VecDestroy( level->workSol );
			VecCreate( self->mgData->comm, &level->workSol );
			VecSetSizes( level->workSol, rowSize, PETSC_DECIDE );
			VecSetFromOptions( level->workSol );
#if( PETSC_VERSION_MAJOR <= 2 && PETSC_VERSION_MINOR >= 3 && PETSC_VERSION_SUBMINOR >= 3 )
			VecSetOption( level->workSol, VEC_IGNORE_NEGATIVE_INDICES );
#elif( PETSC_VERSION_MAJOR >= 3 )
			VecSetOption( level->workSol, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE );
#endif
		}

		VecGetLocalSize( level->workRHS, &vecSize );

		if( !level->workRHS || /*Vector_GetLocalSize( level->workRHS )*/vecSize != rowSize ) {
			//if( level->workRHS )
			//	Stg_Class_RemoveRef( level->workRHS );
			//Vector_Duplicate( self->curSolution, (void**)&level->workRHS );
			//Vector_SetLocalSize( level->workRHS, rowSize );
			if( level->workRHS != PETSC_NULL )
				VecDestroy( level->workRHS );
			VecCreate( self->mgData->comm, &level->workRHS );
			VecSetSizes( level->workRHS, rowSize, PETSC_DECIDE );
			VecSetFromOptions( level->workRHS );
#if( PETSC_VERSION_MAJOR <= 2 && PETSC_VERSION_MINOR >= 3 && PETSC_VERSION_SUBMINOR >= 3 )
			VecSetOption( level->workRHS, VEC_IGNORE_NEGATIVE_INDICES );
#elif( PETSC_VERSION_MAJOR >= 3 )
			VecSetOption( level->workRHS, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE );
#endif
		}
	}
}

MGSolver_PETScData* MultigridSolver_CreateOuterSolver( MultigridSolver* self, Mat matrix ) {
	//MatrixSolver*	outerSolver;
  MGSolver_PETScData* 	outerSolver = (MGSolver_PETScData*)malloc( sizeof(MGSolver_PETScData) );
	PC			pc;

	/*
	outerSolver = (MatrixSolver*)PETScMatrixSolver_New( "" );
	PETScMatrixSolver_SetKSPType( outerSolver, PETScMatrixSolver_KSPType_Richardson );
	PETScMatrixSolver_SetPCType( outerSolver, PETScMatrixSolver_PCType_SOR );
	PETScMatrixSolver_SetMatrix( outerSolver, matrix );
	PETScMatrixSolver_SetMaxIterations( outerSolver, 3 );
	PETScMatrixSolver_SetUseInitialSolution( outerSolver, True );
	PETScMatrixSolver_SetNormType( outerSolver, PETScMatrixSolver_NormType_Preconditioned );
	*/
	KSPCreate( MPI_COMM_WORLD, &outerSolver->ksp );
	KSPSetType( outerSolver->ksp, KSPRICHARDSON );
	KSPGetPC( outerSolver->ksp, &pc );
	PCSetType( pc, PCSOR );
	if( outerSolver->matrix != PETSC_NULL )
		MatDestroy( outerSolver->matrix );
	outerSolver->matrix = matrix;
	KSPSetOperators( outerSolver->ksp, matrix, matrix, DIFFERENT_NONZERO_PATTERN );
	KSPSetTolerances( outerSolver->ksp, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, (PetscInt)3 );
	KSPSetInitialGuessNonzero( outerSolver->ksp, (PetscTruth)True );
	KSPSetNormType( outerSolver->ksp, (KSPNormType)MultigridSolver_NormType_Preconditioned );

	return outerSolver;
}

MGSolver_PETScData* MultigridSolver_CreateSmoother( MultigridSolver* self, Mat matrix ) {
	//MatrixSolver*	smoother;
  MGSolver_PETScData* smoother = (MGSolver_PETScData*)malloc( sizeof( MGSolver_PETScData ) );
	//unsigned	nBlocks;
	//KSP*		ksps;
	PC		pc;
	//PetscErrorCode	ec;

	/*
	smoother = (MatrixSolver*)PETScMatrixSolver_New( "" );
	PETScMatrixSolver_SetKSPType( smoother, 
				      PETScMatrixSolver_KSPType_Richardson );
	PETScMatrixSolver_SetPCType( smoother, 
				     PETScMatrixSolver_PCType_SOR );
	MatrixSolver_SetMatrix( smoother, matrix );
	*/
	KSPCreate( MPI_COMM_WORLD, &smoother->ksp );
	KSPSetType( smoother->ksp, KSPRICHARDSON );
	KSPGetPC( smoother->ksp, &pc );
	PCSetType( pc, PCSOR );
	if( smoother->matrix != PETSC_NULL )
		MatDestroy( smoother->matrix );
	smoother->matrix = matrix;
	KSPSetOperators( smoother->ksp, matrix, matrix, DIFFERENT_NONZERO_PATTERN );

	return smoother;
}

MGSolver_PETScData* MultigridSolver_CreateCoarseSolver( MultigridSolver* self, Mat matrix ) {

  abort();
  return (MGSolver_PETScData*)NULL;
	// //MatrixSolver*	coarseSolver;
	// MGSolver_PETScData* courseSolver;
	// unsigned	nProcs;
	// PC		pc;

	// MPI_Comm_size( self->mgData->comm, (int*)&nProcs );

	// /*
	// coarseSolver = (MatrixSolver*)PETScMatrixSolver_New( "" );
	// PETScMatrixSolver_SetKSPType( coarseSolver, 
	// 			      PETScMatrixSolver_KSPType_PreOnly );
	// if( nProcs == 1 ) {
	// 	PETScMatrixSolver_SetPCType( coarseSolver, 
	// 				     PETScMatrixSolver_PCType_LU );
	// }
	// else {
	// 	PETScMatrixSolver_SetPCType( coarseSolver, 
	// 				     PETScMatrixSolver_PCType_RedundantLU );
	// }
	// MatrixSolver_SetMatrix( coarseSolver, matrix );
	// */
	// KSPCreate( MPI_COMM_WORLD, &courseSolver->ksp );
	// KSPSetType( courseSolver->ksp, KSPPREONLY );
	// KSPGetPC( courseSolver->ksp, &pc );
	// if( nProcs == 1 )
	// 	PCSetType( pc, PCLU );
	// else {
	// 	PCSetType( pc, PCREDUNDANT );
	// 	PCRedundantGetPC( pc, &pc );
	// 	PCSetType( pc, PCLU );
	// }
	// if( courseSolver->matrix != PETSC_NULL )
	// 	MatDestroy( courseSolver->matrix );
	// courseSolver->matrix = matrix;
	// KSPSetOperators( courseSolver->ksp, matrix, matrix, DIFFERENT_NONZERO_PATTERN );

	// return courseSolver;
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
			//Stg_Class_RemoveRef( MatrixSolver_GetMatrix( level->downSolver ) );
			//Stg_Class_RemoveRef( level->downSolver );
			if( level->downSolver->ksp != PETSC_NULL )         KSPDestroy( level->downSolver->ksp );
			if( level->downSolver->matrix != PETSC_NULL )      MatDestroy( level->downSolver->matrix );
			if( level->downSolver->inversion != PETSC_NULL )   MatDestroy( level->downSolver->inversion );
			if( level->downSolver->residual != PETSC_NULL )    VecDestroy( level->downSolver->residual );
			if( level->downSolver->curRHS != PETSC_NULL )      VecDestroy( level->downSolver->curRHS );
			if( level->downSolver->curSolution != PETSC_NULL ) VecDestroy( level->downSolver->curSolution );
			free( level->downSolver );
		}
		if( level->upSolver ) {
			//Stg_Class_RemoveRef( MatrixSolver_GetMatrix( level->upSolver ) );
			//Stg_Class_RemoveRef( level->upSolver );
			if( level->upSolver->ksp != PETSC_NULL )         KSPDestroy( level->upSolver->ksp );
			if( level->upSolver->matrix != PETSC_NULL )      MatDestroy( level->upSolver->matrix );
			if( level->upSolver->inversion != PETSC_NULL )   MatDestroy( level->upSolver->inversion );
			if( level->upSolver->residual != PETSC_NULL )    VecDestroy( level->upSolver->residual );
			if( level->upSolver->curRHS != PETSC_NULL )      VecDestroy( level->upSolver->curRHS );
			if( level->upSolver->curSolution != PETSC_NULL ) VecDestroy( level->upSolver->curSolution );
			free( level->upSolver );
		}
		if( level->R != PETSC_NULL )
			//Stg_Class_RemoveRef( level->R );
			MatDestroy( level->R );
		if( level->P != PETSC_NULL )
			//Stg_Class_RemoveRef( level->P );
			MatDestroy( level->P );
		if( level->workRHS != PETSC_NULL )
			//Stg_Class_RemoveRef( level->workRHS );
			VecDestroy( level->workRHS );
		if( level->workSol != PETSC_NULL )
			//Stg_Class_RemoveRef( level->workSol );
			VecDestroy( level->workSol );
	}

	KillArray( self->levels );
	self->nLevels = 0;
	self->solversChanged = True;
	self->opsChanged = True;

	/* Temporary. */
	//KillObject( self->outerSolver );
	if( self->outerSolver->ksp )         KSPDestroy( self->outerSolver->ksp );
	if( self->outerSolver->matrix )      MatDestroy( self->outerSolver->matrix );
	if( self->outerSolver->inversion )   MatDestroy( self->outerSolver->inversion );
	if( self->outerSolver->residual )    VecDestroy( self->outerSolver->residual );
	if( self->outerSolver->curRHS )      VecDestroy( self->outerSolver->curRHS );
	if( self->outerSolver->curSolution ) VecDestroy( self->outerSolver->curSolution );
	free( self->outerSolver );
}



