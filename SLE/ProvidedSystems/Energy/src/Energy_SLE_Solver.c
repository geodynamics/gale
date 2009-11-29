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
** $Id: Energy_SLE_Solver.c 964 2007-10-11 08:03:06Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include "StgFEM/Discretisation/Discretisation.h"
#include "StgFEM/SLE/SystemSetup/SystemSetup.h"
#include "types.h"
#include "Energy_SLE_Solver.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

/* Textual name of this class */
const Type Energy_SLE_Solver_Type = "Energy_SLE_Solver";

PetscTruth Energy_SLE_HasNullSpace( Mat A );

void* Energy_SLE_Solver_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                              _sizeOfSelf = sizeof(Energy_SLE_Solver);
	Type                                                      type = Energy_SLE_Solver_Type;
	Stg_Class_DeleteFunction*                              _delete = _Energy_SLE_Solver_Delete;
	Stg_Class_PrintFunction*                                _print = _Energy_SLE_Solver_Print;
	Stg_Class_CopyFunction*                                  _copy = _Energy_SLE_Solver_Copy;
	Stg_Component_DefaultConstructorFunction*  _defaultConstructor = Energy_SLE_Solver_DefaultNew;
	Stg_Component_ConstructFunction*                    _construct = _Energy_SLE_Solver_AssignFromXML;
	Stg_Component_BuildFunction*                            _build = _Energy_SLE_Solver_Build;
	Stg_Component_InitialiseFunction*                  _initialise = _Energy_SLE_Solver_Initialise;
	Stg_Component_ExecuteFunction*                        _execute = _SLE_Solver_Execute;
	Stg_Component_DestroyFunction*                        _destroy = _SLE_Solver_Destroy;
	SLE_Solver_SolverSetupFunction*                   _solverSetup = _Energy_SLE_Solver_SolverSetup;
	SLE_Solver_SolveFunction*                               _solve = _Energy_SLE_Solver_Solve;
	SLE_Solver_GetResidualFunc*                       _getResidual = _Energy_SLE_GetResidual;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = ZERO;

	return _Energy_SLE_Solver_New(  ENERGY_SLE_SOLVER_PASSARGS  );
}

Energy_SLE_Solver* Energy_SLE_Solver_New( Name name, Bool useStatSolve, int statReps ) {
	Energy_SLE_Solver* self = (Energy_SLE_Solver*) Energy_SLE_Solver_DefaultNew( name );

	Energy_SLE_Solver_InitAll( self, useStatSolve, statReps ) ;

	return self;
}

/* Creation implementation / Virtual constructor */
Energy_SLE_Solver* _Energy_SLE_Solver_New(  ENERGY_SLE_SOLVER_DEFARGS  )
{
	Energy_SLE_Solver* self;

	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(Energy_SLE_Solver) );
	self = (Energy_SLE_Solver*) _SLE_Solver_New(  SLE_SOLVER_PASSARGS  );
	
	/* Virtual info */
	return self;
}
	
/* Initialisation implementation */
void _Energy_SLE_Solver_Init( Energy_SLE_Solver* self ) {
	self->isConstructed = True;
	self->residual = NULL;
	//self->matrixSolver = NULL;
	self->matrixSolver = PETSC_NULL;
}
void Energy_SLE_Solver_InitAll( Energy_SLE_Solver* solver, Bool useStatSolve, int statReps ) {
	Energy_SLE_Solver* self = (Energy_SLE_Solver*)solver;

	SLE_Solver_InitAll( self, useStatSolve, statReps );
}

void _Energy_SLE_Solver_Delete( void* sle ) {
	Energy_SLE_Solver* self = (Energy_SLE_Solver*)sle;

	//FreeObject( self->matrixSolver );
	KSPDestroy( self->matrixSolver );

	if( self->residual != PETSC_NULL ) {
		//FreeObject( self->residual );
		VecDestroy( self->residual );
	}
}

void _Energy_SLE_Solver_Print( void* solver, Stream* stream ) {
	Energy_SLE_Solver* self = (Energy_SLE_Solver*)solver;
	
	/* Set the Journal for printing informations */
	Stream* standard_SLE_SolverStream = stream;
	
	/* General info */
	Journal_Printf( standard_SLE_SolverStream, "Energy_SLE_Solver (ptr): %p\n", self );

	/* other info */
	Journal_Printf( standard_SLE_SolverStream, "\tmatrixSolver (ptr): %p", self->matrixSolver );
}


void* _Energy_SLE_Solver_Copy( void* standardSleSolver, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	Energy_SLE_Solver*	self = (Energy_SLE_Solver*)standardSleSolver;
	Energy_SLE_Solver*	newEnergySleSolver;

	newEnergySleSolver = _SLE_Solver_Copy( self, dest, deep, nameExt, ptrMap );
	
	newEnergySleSolver->matrixSolver = self->matrixSolver;
	
	return (void*)newEnergySleSolver;
}

void _Energy_SLE_Solver_AssignFromXML( void* sleSolver, Stg_ComponentFactory* cf, void* data ) {
	Energy_SLE_Solver	*self = (Energy_SLE_Solver*)sleSolver;

	assert( self && Stg_CheckType( self, Energy_SLE_Solver ) );
	assert( cf && Stg_CheckType( cf, Stg_ComponentFactory ) );

	_SLE_Solver_AssignFromXML( self, cf, data );
}

/* Build */
void _Energy_SLE_Solver_Build( void* sleSolver, void* standardSLE ) {
	Energy_SLE_Solver*	self = (Energy_SLE_Solver*)sleSolver;
	SystemLinearEquations*	sle = (SystemLinearEquations*) standardSLE;
	StiffnessMatrix*	stiffMat = (StiffnessMatrix*)sle->stiffnessMatrices->data[0];

	Journal_DPrintf( self->debug, "In %s()\n", __func__ );
	Stream_IndentBranch( StgFEM_SLE_ProvidedSystems_Energy_Debug );
	Journal_DPrintf( self->debug, "building a standard L.A. solver for the \"%s\" matrix.\n", stiffMat->name );

	Stg_Component_Build( stiffMat, standardSLE, False );

	if( self->matrixSolver == PETSC_NULL )
		KSPCreate( sle->comm, &self->matrixSolver );

	Stream_UnIndentBranch( StgFEM_SLE_ProvidedSystems_Energy_Debug );
}


void _Energy_SLE_Solver_Initialise( void* sleSolver, void* standardSLE ) {
	Energy_SLE_Solver*	self = (Energy_SLE_Solver*)sleSolver;
	SystemLinearEquations*	sle = (SystemLinearEquations*) standardSLE;
	
	/* Initialise parent. */
	_SLE_Solver_Initialise( self, sle );
}

void _Energy_SLE_Solver_Execute( void* sleSolver, void* data ) {
}
	
void _Energy_SLE_Solver_Destroy( void* sleSolver, void* data ) {
}

void _Energy_SLE_Solver_SolverSetup( void* sleSolver, void* standardSLE ) {
	Energy_SLE_Solver* self = (Energy_SLE_Solver*)sleSolver;
	SystemLinearEquations* sle = (SystemLinearEquations*) standardSLE;
	StiffnessMatrix*	stiffMat = (StiffnessMatrix*)sle->stiffnessMatrices->data[0];
	
	Journal_DPrintf( self->debug, "In %s()\n", __func__ );
	Stream_IndentBranch( StgFEM_SLE_ProvidedSystems_Energy_Debug );
	
	Journal_DPrintf( self->debug, "Initialising the L.A. solver for the \"%s\" matrix.\n", stiffMat->name );
	KSPSetOperators( self->matrixSolver, stiffMat->matrix, stiffMat->matrix, DIFFERENT_NONZERO_PATTERN );
	Stream_UnIndentBranch( StgFEM_SLE_ProvidedSystems_Energy_Debug );
	
	if( self->maxIterations > 0 ) {
		KSPSetTolerances( self->matrixSolver, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT, self->maxIterations );
	}

	VecDuplicate( ((ForceVector**)sle->forceVectors->data)[0]->vector, &self->residual );
}


void _Energy_SLE_Solver_Solve( void* sleSolver, void* standardSLE ) {
	Energy_SLE_Solver*     self       = (Energy_SLE_Solver*)sleSolver;
	SystemLinearEquations* sle        = (SystemLinearEquations*) standardSLE;
	Iteration_Index        iterations;
        PetscTruth isNull;
        MatNullSpace nsp;

	
	Journal_DPrintf( self->debug, "In %s - for standard SLE solver\n", __func__ );
	Stream_IndentBranch( StgFEM_SLE_ProvidedSystems_Energy_Debug );
	
	if ( (sle->stiffnessMatrices->count > 1 ) ||
		(sle->forceVectors->count > 1 ) )
	{
		Stream* warning = Stream_RegisterChild( StgFEM_Warning, self->type );
		Journal_Printf( warning, "Warning: Energy solver unable to solve more that one matrix/vector.\n" );
	}

	if ( sle->solutionVectors->count > 1 ) {
		Stream* warning = Stream_RegisterChild( StgFEM_Warning, self->type );
		Journal_Printf( warning, "Warning: More than 1 solution vector provided to standard sle. Ignoring second and subsequent"
			" solution vectors.\n" );
	}

        isNull = Energy_SLE_HasNullSpace(((StiffnessMatrix**)sle->stiffnessMatrices->data)[0]->matrix);
        if(isNull) {
            MatNullSpaceCreate(PETSC_COMM_WORLD, PETSC_TRUE, 0, PETSC_NULL, &nsp);
            KSPSetNullSpace(self->matrixSolver, nsp);
        }

	KSPSolve( self->matrixSolver,
		    ((ForceVector*) sle->forceVectors->data[0])->vector, 
		    ((SolutionVector*) sle->solutionVectors->data[0])->vector );
	KSPGetIterationNumber( self->matrixSolver, &iterations );

	Journal_DPrintf( self->debug, "Solved after %u iterations.\n", iterations );
	Stream_UnIndentBranch( StgFEM_SLE_ProvidedSystems_Energy_Debug );
	
	/* calculate the residual */
	/* TODO: need to make this optional */
	MatMult( ((StiffnessMatrix**)sle->stiffnessMatrices->data)[0]->matrix, 
		 ((SolutionVector**)sle->solutionVectors->data)[0]->vector, 
		 self->residual );
	VecAYPX( self->residual, -1.0, ((ForceVector**)sle->forceVectors->data)[0]->vector );

        if(isNull)
            MatNullSpaceDestroy(nsp);
}


Vec _Energy_SLE_GetResidual( void* sleSolver, Index fv_I ) {
	Energy_SLE_Solver*	self = (Energy_SLE_Solver*)sleSolver;

	return self->residual;
}

PetscTruth Energy_SLE_HasNullSpace( Mat A ) {
    PetscInt N;
    PetscScalar sum;
    PetscReal nrm;
    PetscTruth isNull;
    Vec r, l;

    MatGetVecs(A, &r, &l); /* l = A r */

    VecGetSize(r, &N);
    sum = 1.0/(PetscScalar)N;
    VecSet(r, sum);

    MatMult(A, r, l); /* {l} = [A] {r} */

    VecNorm(l, NORM_2, &nrm);
    if(nrm < 1.e-7)
	isNull = PETSC_TRUE;
    else
	isNull = PETSC_FALSE;

    VecDestroy(l);
    VecDestroy(r);

    return isNull;
}


