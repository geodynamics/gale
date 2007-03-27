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


MatrixSolver* MatrixSolver_Build( MPI_Comm comm, void* matrix ) {
	return (MatrixSolver*)0;
}


void MatrixSolver_Setup( void* matSolver, void* matrix )
{
}

void MatrixSolver_SetupKSP( void* solver ) {
}

void MatrixSolver_Setup_SameNonZeroPattern( void* matSolver, void* matrix )
{
}

void MatrixSolver_SetPC_Type( void* matSolver, char* pcType ) {
}

void MatrixSolver_SetKSP_Type( void* matSolver, char* kspType ) {
}


void MatrixSolver_SetMaxIts( void* solverCtx, Iteration_Index maxIts )
{
}

void MatrixSolver_SetRelativeTolerance( void* matSolver, double relativeTolerance )
{
}

void MatrixSolver_SetAbsoluteTolerance( void* matSolver, double absoluteTolerance )
{
}




Iteration_Index MatrixSolver_Solve( MatrixSolver* matSolver, Vector* solnVec, Vector* rhsVec ) {
	return 0;
}


Iteration_Index MatrixSolver_StatSolve( MatrixSolver* solver, 
								Vector* sol, Vector* rhs, 
								unsigned nReps )
{
	return 0;
}


double MatrixSolver_GetResidualNorm( MatrixSolver* solver ) {
	return 0.0;
}


void MatrixSolver_HaveInitialGuess( MatrixSolver* solver ) {
}


void MatrixSolver_PC_GetBJacobiSubBlocks( MatrixSolver* solver, MatrixSolver*** blocks, unsigned* nBlocks ) {
}


void MatrixSolver_PC_SetSORIts( MatrixSolver* solver, unsigned nIts ) {
}


void MatrixSolver_MG_SetLevels( MatrixSolver* solver, unsigned nLevels ) {
}


void MatrixSolver_MG_SetCycles( MatrixSolver* solver, unsigned level, unsigned nLevels ) {
}


void MatrixSolver_MG_SetLevelMatrix( MatrixSolver* solver, unsigned level, Matrix* mat ) {
}


void MatrixSolver_MG_SetSmoothUpIts( MatrixSolver* solver, unsigned nUpIts ) {
}


void MatrixSolver_MG_SetSmoothDownIts( MatrixSolver* solver, unsigned nDownIts ) {
}


void MatrixSolver_MG_GetDownSmoother( MatrixSolver* solver, unsigned level, MatrixSolver** smoother ) {
}


void MatrixSolver_MG_GetUpSmoother( MatrixSolver* solver, unsigned level, MatrixSolver** smoother ) {
}


void MatrixSolver_MG_SetRestrictionOp( MatrixSolver* solver, unsigned level, Matrix* op ) {
}


void MatrixSolver_MG_SetInterpolationOp( MatrixSolver* solver, unsigned level, Matrix* op ) {
}


void MatrixSolver_MG_SetWorkSpace( MatrixSolver* solver, unsigned level, Vector* rhsVec, Vector* xVec, Vector* rVec ) {
}


void MatrixSolver_MG_ForceMonitor( MatrixSolver* solver ) {
}
