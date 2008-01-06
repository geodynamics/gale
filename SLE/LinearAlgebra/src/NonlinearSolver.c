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
** $Id: NonlinearSolver.c 3584 2006-05-16 11:11:07Z PatrickSunter $
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
const Type NonlinearSolver_Type = "NonlinearSolver";


/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/

NonlinearSolver* _NonlinearSolver_New( NONLINEARSOLVER_DEFARGS ) {
	NonlinearSolver*	self;

	/* Allocate memory */
	assert( sizeOfSelf >= sizeof(NonlinearSolver) );
	self = (NonlinearSolver*)_Stg_Component_New( STG_COMPONENT_PASSARGS );

	/* Virtual info */
	self->setCommFunc = setCommFunc;
	self->destroyFunc = destroyFunc;
	self->setFunctionFunc = setFunctionFunc;
	self->getJacobianFunc = getJacobianFunc;
	self->solveFunc = solveFunc;
	self->setSolutionFunc = setSolutionFunc;
	self->setRhsFunc = setRhsFunc;

	self->getSolveStatusFunc = getSolveStatusFunc;
	self->getIterationsFunc = getIterationsFunc;
	self->getMaxIterationsFunc = getMaxIterationsFunc;
	self->getResidualNormFunc = getResidualNormFunc;

	/* NonlinearSolver info */
	_NonlinearSolver_Init( self );

	return self;
}

void _NonlinearSolver_Init( NonlinearSolver* self ) {
	assert( self && Stg_CheckType( self, NonlinearSolver ) );

	self->comm = MPI_COMM_WORLD;
	self->J = NULL;
	self->Jinv = NULL;
	self->residual = NULL;
	self->expiredResidual = True;

	self->curRHS = NULL;
	self->curSolution = NULL;
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _NonlinearSolver_Delete( void* nls ) {
	NonlinearSolver*	self = (NonlinearSolver*)nls;

	assert( self && Stg_CheckType( self, NonlinearSolver ) );

	/* Delete the parent. */
	_Stg_Component_Delete( self );
}

void _NonlinearSolver_Print( void* nls, Stream* stream ) {
	NonlinearSolver*	self = (NonlinearSolver*)nls;
	
	/* Set the Journal for printing informations */
	Stream* nlsStream;
	nlsStream = Journal_Register( InfoStream_Type, "NonlinearSolverStream" );

	assert( self && Stg_CheckType( self, NonlinearSolver ) );

	/* Print parent */
	Journal_Printf( stream, "NonlinearSolver (ptr): (%p)\n", self );
	_Stg_Component_Print( self, stream );
}

void _NonlinearSolver_Construct( void* nls, Stg_ComponentFactory* cf, void* data ) {
	NonlinearSolver*		self = (NonlinearSolver*)nls;

	assert( self && Stg_CheckType( self, NonlinearSolver ) );
	assert( cf );
}

void _NonlinearSolver_Build( void* nls, void* data ) {
}

void _NonlinearSolver_Initialise( void* nls, void* data ) {
}

void _NonlinearSolver_Execute( void* nls, void* data ) {
}

void _NonlinearSolver_Destroy( void* nls, void* data ) {
}

void _NonlinearSolver_SetComm( void* nls, MPI_Comm comm ) {
	NonlinearSolver*	self = (NonlinearSolver*)nls;

	assert( self && Stg_CheckType( self, NonlinearSolver ) );

	self->comm = comm;
}

/*--------------------------------------------------------------------------------------------------------------------------
** Public Functions
*/


/*----------------------------------------------------------------------------------------------------------------------------------
** Private Functions
*/

