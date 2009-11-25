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
** $Id: MGOpGenerator.c 3584 2006-05-16 11:11:07Z PatrickSunter $
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
const Type MGOpGenerator_Type = "MGOpGenerator";


/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/

MGOpGenerator* _MGOpGenerator_New( MGOPGENERATOR_DEFARGS ) {
	MGOpGenerator*	self;

	/* Allocate memory */
	assert( sizeOfSelf >= sizeof(MGOpGenerator) );
	self = (MGOpGenerator*)_Stg_Component_New( STG_COMPONENT_PASSARGS );

	/* stuff previously housed in the MatrixSolver class */
	self->solver = malloc( sizeof( MGSolver_PETScData ) );
	self->solver->ksp = PETSC_NULL;
	self->solver->matrix = PETSC_NULL;
	self->solver->inversion = PETSC_NULL;
	self->solver->residual = PETSC_NULL;
	self->solver->curRHS = PETSC_NULL;
	self->solver->curSolution = PETSC_NULL;
	self->solver->expiredResidual = True;
	self->solver->matrixChanged = True;

	/* Virtual info */
	self->setNumLevelsFunc = setNumLevelsFunc;
	self->hasExpiredFunc = hasExpiredFunc;
	self->generateFunc = generateFunc;

	/* MGOpGenerator info */
	_MGOpGenerator_Init( self );

	return self;
}

void _MGOpGenerator_Init( MGOpGenerator* self ) {
	assert( self && Stg_CheckType( self, MGOpGenerator ) );

	self->nLevels = 0;
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _MGOpGenerator_Delete( void* mgOpGenerator ) {
	MGOpGenerator*	self = (MGOpGenerator*)mgOpGenerator;

	assert( self && Stg_CheckType( self, MGOpGenerator ) );

	/* this stuff was previously taken care of in the MatrixSolver class */
	if( self->solver ) {
            /*if( self->solver->ksp != PETSC_NULL )         KSPDestroy( self->solver->ksp );*/
		if( self->solver->matrix != PETSC_NULL )      MatDestroy( self->solver->matrix );
		if( self->solver->inversion != PETSC_NULL )   MatDestroy( self->solver->inversion );
		if( self->solver->residual != PETSC_NULL )    VecDestroy( self->solver->residual );
		if( self->solver->curRHS != PETSC_NULL )      VecDestroy( self->solver->curRHS );
		if( self->solver->curSolution != PETSC_NULL ) VecDestroy( self->solver->curSolution );
		free( self->solver );
	}

	/* Delete the parent. */
	_Stg_Component_Delete( self );
}

void _MGOpGenerator_Print( void* mgOpGenerator, Stream* stream ) {
	MGOpGenerator*	self = (MGOpGenerator*)mgOpGenerator;
	
	/* Set the Journal for printing informations */
	Stream* mgOpGeneratorStream;
	mgOpGeneratorStream = Journal_Register( InfoStream_Type, "MGOpGeneratorStream" );

	assert( self && Stg_CheckType( self, MGOpGenerator ) );

	/* Print parent */
	Journal_Printf( stream, "MGOpGenerator (ptr): (%p)\n", self );
	_Stg_Component_Print( self, stream );
}

void _MGOpGenerator_AssignFromXML( void* mgOpGenerator, Stg_ComponentFactory* cf, void* data ) {
	MGOpGenerator*		self = (MGOpGenerator*)mgOpGenerator;

	assert( self && Stg_CheckType( self, MGOpGenerator ) );
	assert( cf );
}

void _MGOpGenerator_Build( void* mgOpGenerator, void* data ) {
}

void _MGOpGenerator_Initialise( void* mgOpGenerator, void* data ) {
}

void _MGOpGenerator_Execute( void* mgOpGenerator, void* data ) {
}

void _MGOpGenerator_Destroy( void* mgOpGenerator, void* data ) {
}

void _MGOpGenerator_SetNumLevels( void* mgOpGenerator, unsigned nLevels ) {
	MGOpGenerator*		self = (MGOpGenerator*)mgOpGenerator;

	assert( self && Stg_CheckType( self, MGOpGenerator ) );

	self->nLevels = nLevels;
}


/*--------------------------------------------------------------------------------------------------------------------------
** Public Functions
*/

void MGOpGenerator_SetMatrixSolver( void* mgOpGenerator, void* _solver ) {
	MGOpGenerator*	self = (MGOpGenerator*)mgOpGenerator;
	//MatrixSolver*	solver = (MatrixSolver*)_solver;
	KSP		solver = (KSP)_solver;

	assert( self && Stg_CheckType( self, MGOpGenerator ) );
	//assert( solver && Stg_CheckType( solver, MatrixSolver ) );

	//self->solver = solver;
	self->solver->ksp = solver;
}

unsigned MGOpGenerator_GetNumLevels( void* mgOpGenerator ) {
	MGOpGenerator*		self = (MGOpGenerator*)mgOpGenerator;

	assert( self && Stg_CheckType( self, MGOpGenerator ) );

	return self->nLevels;
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Private Functions
*/
