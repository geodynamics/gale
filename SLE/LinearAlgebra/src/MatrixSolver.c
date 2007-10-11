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
** $Id: MatrixSolver.c 3584 2006-05-16 11:11:07Z PatrickSunter $
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
const Type MatrixSolver_Type = "MatrixSolver";


/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/

MatrixSolver* _MatrixSolver_New( MATRIXSOLVER_DEFARGS ) {
	MatrixSolver*	self;

	/* Allocate memory */
	assert( sizeOfSelf >= sizeof(MatrixSolver) );
	self = (MatrixSolver*)_Stg_Component_New( STG_COMPONENT_PASSARGS );

	/* Virtual info */
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

	/* MatrixSolver info */
	_MatrixSolver_Init( self );

	return self;
}

void _MatrixSolver_Init( MatrixSolver* self ) {
	assert( self && Stg_CheckType( self, MatrixSolver ) );

	self->comm = MPI_COMM_WORLD;
	self->matrix = NULL;
	self->inversion = NULL;
	self->residual = NULL;
	self->expiredResidual = True;
	self->matrixChanged = True;

	self->curRHS = NULL;
	self->curSolution = NULL;
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _MatrixSolver_Delete( void* matrixSolver ) {
	MatrixSolver*	self = (MatrixSolver*)matrixSolver;

	assert( self && Stg_CheckType( self, MatrixSolver ) );

	/* Delete the parent. */
	_Stg_Component_Delete( self );
}

void _MatrixSolver_Print( void* matrixSolver, Stream* stream ) {
	MatrixSolver*	self = (MatrixSolver*)matrixSolver;
	
	/* Set the Journal for printing informations */
	Stream* matrixSolverStream;
	matrixSolverStream = Journal_Register( InfoStream_Type, "MatrixSolverStream" );

	assert( self && Stg_CheckType( self, MatrixSolver ) );

	/* Print parent */
	Journal_Printf( stream, "MatrixSolver (ptr): (%p)\n", self );
	_Stg_Component_Print( self, stream );
}

void _MatrixSolver_Construct( void* matrixSolver, Stg_ComponentFactory* cf, void* data ) {
	MatrixSolver*		self = (MatrixSolver*)matrixSolver;

	assert( self && Stg_CheckType( self, MatrixSolver ) );
	assert( cf );
}

void _MatrixSolver_Build( void* matrixSolver, void* data ) {
}

void _MatrixSolver_Initialise( void* matrixSolver, void* data ) {
}

void _MatrixSolver_Execute( void* matrixSolver, void* data ) {
}

void _MatrixSolver_Destroy( void* matrixSolver, void* data ) {
}

void _MatrixSolver_SetComm( void* matrixSolver, MPI_Comm comm ) {
	MatrixSolver*	self = (MatrixSolver*)matrixSolver;

	assert( self && Stg_CheckType( self, MatrixSolver ) );

	self->comm = comm;
}

void _MatrixSolver_SetMatrix( void* matrixSolver, void* _matrix ) {
	MatrixSolver*	self = (MatrixSolver*)matrixSolver;
	Matrix*		matrix = (Matrix*)_matrix;

	assert( self && Stg_CheckType( self, MatrixSolver ) );
	assert( !matrix || Stg_CheckType( matrix, Matrix ) );

	if( matrix == self->matrix )
		return;

	if( self->matrix ) {
		List_Remove( self->matrix->solvers, &self );
		Stg_Class_RemoveRef( self->matrix );
	}
	self->matrix = matrix;
	if( matrix ) {
		Stg_Class_AddRef( matrix );
		List_Append( matrix->solvers, &self );
	}
	self->matrixChanged = True;
}

void _MatrixSolver_Setup( void* matrixSolver, void* rhs, void* solution ) {
	MatrixSolver*	self = (MatrixSolver*)matrixSolver;

	assert( self && Stg_CheckType( self, MatrixSolver ) );
	assert( rhs && Stg_CheckType( rhs, Vector ) );
	assert( solution && Stg_CheckType( solution, Vector ) );

	self->curRHS = rhs;
	self->curSolution = solution;
	self->expiredResidual = True;
}


/*--------------------------------------------------------------------------------------------------------------------------
** Public Functions
*/

MPI_Comm MatrixSolver_GetComm( void* matrixSolver ) {
	MatrixSolver*	self = (MatrixSolver*)matrixSolver;

	assert( self && Stg_CheckType( self, MatrixSolver ) );

	return self->comm;
}

Matrix* MatrixSolver_GetMatrix( void* matrixSolver ) {
	MatrixSolver*	self = (MatrixSolver*)matrixSolver;

	assert( self && Stg_CheckType( self, MatrixSolver ) );

	return self->matrix;
}

Vector* MatrixSolver_GetResidual( void* matrixSolver ) {
	MatrixSolver*	self = (MatrixSolver*)matrixSolver;

	assert( self && Stg_CheckType( self, MatrixSolver ) );

	if( !self->residual || self->expiredResidual ) {
		MatrixSolver_CalcResidual( self );
		self->expiredResidual = False;
	}

	return self->residual;
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Private Functions
*/

void MatrixSolver_CalcResidual( MatrixSolver* self ) {
	assert( self && Stg_CheckType( self, MatrixSolver ) );
	assert( self->curSolution && self->curRHS && self->matrix );

	if( !self->residual ) {
		Vector_Duplicate( self->curSolution, (void**)&self->residual );
		Vector_SetLocalSize( self->residual, Vector_GetLocalSize( self->curSolution ) );
	}
	Matrix_Multiply( self->matrix, self->curSolution, self->residual );
	Vector_ScaleAdd( self->residual, -1.0, self->curRHS );
}
