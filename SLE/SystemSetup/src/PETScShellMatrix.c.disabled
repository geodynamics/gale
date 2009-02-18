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
** $Id: PETScShellMatrix.c 3584 2006-05-16 11:11:07Z PatrickSunter $
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
#include "SystemSetup.h"


/* Textual name of this class */
const Type PETScShellMatrix_Type = "PETScShellMatrix";


/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/

PETScShellMatrix* PETScShellMatrix_New( Name name ) {
	return _PETScShellMatrix_New( sizeof(PETScShellMatrix), 
				      PETScShellMatrix_Type, 
				      _PETScShellMatrix_Delete, 
				      _PETScShellMatrix_Print, 
				      NULL, 
				      (void* (*)(Name))PETScShellMatrix_New, 
				      _PETScShellMatrix_Construct, 
				      _PETScShellMatrix_Build, 
				      _PETScShellMatrix_Initialise, 
				      _PETScShellMatrix_Execute, 
				      _PETScShellMatrix_Destroy, 
				      name, 
				      NON_GLOBAL, 
				      PETScShellMatrix_SetComm, 
				      PETScShellMatrix_SetGlobalSize, 
				      PETScShellMatrix_SetLocalSize, 
				      PETScShellMatrix_SetNonZeroStructure );
}

PETScShellMatrix* _PETScShellMatrix_New( PETSCSHELLMATRIX_DEFARGS ) {
	PETScShellMatrix*	self;

	/* Allocate memory */
	assert( sizeOfSelf >= sizeof(PETScShellMatrix) );
	//self = (PETScShellMatrix*)_PETScMatrix_New( PETSCMATRIX_PASSARGS );
	self = (PETScShellMatrix*)_Stg_Component_New( STG_COMPONENT_PASSARGS );

	/* Virtual info */
	self->setCommFunc = setCommFunc;
	self->setGlobalSizeFunc = setGlobalSizeFunc;
	self->setLocalSizeFunc = setLocalSizeFunc;
	self->setNonZeroStructure = setNonZeroStructure;

	/* PETScShellMatrix info */
	_PETScShellMatrix_Init( self );

	return self;
}

void _PETScShellMatrix_Init( PETScShellMatrix* self ) {
	assert( self && Stg_CheckType( self, PETScShellMatrix ) );

	self->nRowDofs = 0;
	self->nColDofs = 0;
	self->elStiffMat = NULL;

	self->matrix = PETSC_NULL;
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _PETScShellMatrix_Delete( void* matrix ) {
	PETScShellMatrix*	self = (PETScShellMatrix*)matrix;

	assert( self && Stg_CheckType( self, PETScShellMatrix ) );

	/* Delete the parent. */
	//_PETScMatrix_Delete( self );
	if( self->matrix != PETSC_NULL )
		MatDestroy( self->matrix );
}

void _PETScShellMatrix_Print( void* matrix, Stream* stream ) {
	PETScShellMatrix*	self = (PETScShellMatrix*)matrix;
	
	/* Set the Journal for printing informations */
	Stream* matrixStream;
	matrixStream = Journal_Register( InfoStream_Type, "PETScShellMatrixStream" );

	assert( self && Stg_CheckType( self, PETScShellMatrix ) );

	/* Print parent */
	Journal_Printf( stream, "PETScShellMatrix (ptr): (%p)\n", self );
	//_PETScMatrix_Print( self, stream );
	_Stg_Component_Print( self, stream );
}

void _PETScShellMatrix_Construct( void* matrix, Stg_ComponentFactory* cf, void* data ) {
	PETScShellMatrix*	self = (PETScShellMatrix*)matrix;

	assert( self && Stg_CheckType( self, PETScShellMatrix ) );
	assert( cf );

	//_PETScMatrix_Construct( self, cf, data );

	self->stiffMat = Stg_ComponentFactory_ConstructByKey( cf, self->name, "stiffnessMatrix", StiffnessMatrix, True, data );
	self->sle = Stg_ComponentFactory_ConstructByKey( cf, self->name, "sle", SystemLinearEquations, True, data );
}

void _PETScShellMatrix_Build( void* matrix, void* data ) {
	PETScShellMatrix*	self = (PETScShellMatrix*)matrix;

	assert( self && Stg_CheckType( self, PETScShellMatrix ) );

	Stg_Component_Build( self->sle, data, False );
	Stg_Component_Build( self->stiffMat, data, False );
}

void _PETScShellMatrix_Initialise( void* matrix, void* data ) {
	PETScShellMatrix*	self = (PETScShellMatrix*)matrix;

	assert( self && Stg_CheckType( self, PETScShellMatrix ) );

	Stg_Component_Initialise( self->sle, data, False );
	Stg_Component_Initialise( self->stiffMat, data, False );

	PETScShellMatrix_UpdateAssembly( self );
}

void _PETScShellMatrix_Execute( void* matrix, void* data ) {
}

void _PETScShellMatrix_Destroy( void* matrix, void* data ) {
}

void PETScShellMatrix_SetComm( void* matrix, MPI_Comm comm ) {
	PETScShellMatrix*	self = (PETScShellMatrix*)matrix;
	PetscErrorCode	ec;

	assert( self && Stg_CheckType( self, PETScShellMatrix ) );

	//_Matrix_SetComm( self, comm );
	self->comm = comm;

	if( self->matrix != PETSC_NULL ) {
		ec = MatDestroy( self->matrix );
		CheckPETScError( ec );
	}

	self->hasChanged = True;
}

void PETScShellMatrix_SetGlobalSize( void* matrix, unsigned nRows, unsigned nColumns ) {
	PETScShellMatrix*	self = (PETScShellMatrix*)matrix;
	PetscErrorCode	ec;

	assert( self && Stg_CheckType( self, PETScShellMatrix ) );

	ec = MatCreateShell( self->comm, PETSC_DETERMINE, PETSC_DETERMINE, (PetscInt)nRows, (PetscInt)nColumns, 
			     self, &self->matrix );
	CheckPETScError( ec );
	ec = MatSetFromOptions( self->matrix );
	CheckPETScError( ec );
	ec = MatShellSetOperation( self->matrix, MATOP_MULT, (void (*)(void))PETScShellMatrix_MatMult );
	CheckPETScError( ec );
	ec = MatShellSetOperation( self->matrix, MATOP_MULT_ADD, (void (*)(void))PETScShellMatrix_MatMultAdd );
	CheckPETScError( ec );
	ec = MatShellSetOperation( self->matrix, MATOP_MULT_TRANSPOSE, (void (*)(void))PETScShellMatrix_MatMultTranspose );
	CheckPETScError( ec );
	ec = MatShellSetOperation( self->matrix, MATOP_GET_DIAGONAL, (void (*)(void))PETScShellMatrix_MatGetDiagonal );
	CheckPETScError( ec );

	self->hasChanged = True;
}

void PETScShellMatrix_SetLocalSize( void* matrix, unsigned nRows, unsigned nColumns ) {
	PETScShellMatrix*	self = (PETScShellMatrix*)matrix;
	PetscErrorCode	ec;

	assert( self && Stg_CheckType( self, PETScShellMatrix ) );

	ec = MatCreateShell( self->comm, (PetscInt)nRows, (PetscInt)nColumns, PETSC_DETERMINE, PETSC_DETERMINE, 
			     self, &self->matrix );
	CheckPETScError( ec );
	ec = MatSetFromOptions( self->matrix );
	CheckPETScError( ec );
	ec = MatShellSetOperation( self->matrix, MATOP_MULT, (void (*)(void))PETScShellMatrix_MatMult );
	CheckPETScError( ec );
	ec = MatShellSetOperation( self->matrix, MATOP_MULT_ADD, (void (*)(void))PETScShellMatrix_MatMultAdd );
	CheckPETScError( ec );
	ec = MatShellSetOperation( self->matrix, MATOP_MULT_TRANSPOSE, (void (*)(void))PETScShellMatrix_MatMultTranspose );
	CheckPETScError( ec );
	ec = MatShellSetOperation( self->matrix, MATOP_GET_DIAGONAL, (void (*)(void))PETScShellMatrix_MatGetDiagonal );
	CheckPETScError( ec );

	self->hasChanged = True;
}

void PETScShellMatrix_SetNonZeroStructure( void* matrix, unsigned nNonZeros, 
					   unsigned* diagonalNonZeroIndices, unsigned* offDiagonalNonZeroIndices)
{
	/* Don't do anything. */
}


/*--------------------------------------------------------------------------------------------------------------------------
** Public Functions
*/

void PETScShellMatrix_UpdateAssembly( void* matrix ) {
	PETScShellMatrix*	self = (PETScShellMatrix*)matrix;
	StiffnessMatrix*	stiffMat;
	SystemLinearEquations*	sle;
	FeVariable		*rowVar, *colVar;
	FeMesh			*rowMesh, *colMesh;
	FeEquationNumber	*rowEqNum, *colEqNum;
	DofLayout		*rowDofs, *colDofs;
	unsigned		nRowNodes, *rowNodes;
	unsigned		nColNodes, *colNodes;
	unsigned		nDofs, nRowDofs, nColDofs;
	double			**elStiffMat, **elStiffMatTrans;
	unsigned		n_i, dof_i, dof_j;

	assert( self && Stg_CheckType( self, PETScShellMatrix ) );

	FreeArray( self->elStiffMat );

	stiffMat = self->stiffMat;							assert( stiffMat );
	sle = self->sle;
	rowVar = stiffMat->rowVariable;							assert( rowVar );
	colVar = stiffMat->columnVariable ? stiffMat->columnVariable : rowVar;		assert( colVar );
	rowEqNum = rowVar->eqNum;
	colEqNum = colVar->eqNum;
	rowMesh = rowVar->feMesh;
	colMesh = colVar->feMesh;
	rowDofs = rowVar->dofLayout;							assert( rowDofs );
	colDofs = colVar->dofLayout;							assert( colDofs );

	assert( FeMesh_GetElementLocalSize( rowMesh ) );
	/*
	** Taken out while adding in AMR, need to fix this.
	**
	FeMesh_GetElementNodes( rowMesh, 0, &nRowNodes, &rowNodes );
	FeMesh_GetElementNodes( colMesh, 0, &nColNodes, &colNodes );
	*/
	abort();

	nRowDofs = 0;
	for( n_i = 0; n_i < nRowNodes; n_i++ )
		nRowDofs += rowDofs->dofCounts[rowNodes[n_i]];
	nColDofs = 0;
	for( n_i = 0; n_i < nColNodes; n_i++ )
		nColDofs += colDofs->dofCounts[colNodes[n_i]];
	nDofs = nRowDofs * nColDofs;
	elStiffMat = AllocArray2D( double, nRowDofs, nColDofs );

	memset( elStiffMat[0], 0, nDofs * sizeof(double) );
	StiffnessMatrix_AssembleElement( stiffMat, 0, sle, NULL, elStiffMat );

	elStiffMatTrans = AllocArray2D( double, nColDofs, nRowDofs );
	for( dof_i = 0; dof_i < nRowDofs; dof_i++ ) {
		for( dof_j = 0; dof_j < nColDofs; dof_j++ )
			elStiffMatTrans[dof_j][dof_i] = elStiffMat[dof_i][dof_j];
	}

	self->nRowDofs = nRowDofs;
	self->nColDofs = nColDofs;
	self->elStiffMat = elStiffMat;
	self->elStiffMatTrans = elStiffMatTrans;
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Private Functions
*/

PetscErrorCode PETScShellMatrix_MatMult( Mat A, Vec x, Vec y ) {
	PetscErrorCode	ec;

	ec = VecSet( y, 0.0 );
	CheckPETScError( ec );

	return PETScShellMatrix_MatMultAdd( A, x, y );
}

PetscErrorCode PETScShellMatrix_MatMultAdd( Mat A, Vec x, Vec y ) {
	PETScShellMatrix*	self;
	StiffnessMatrix*	stiffMat;
	SystemLinearEquations*	sle;
	FeVariable		*rowVar, *colVar;
	FeMesh			*rowMesh;
	FeEquationNumber	*rowEqNum, *colEqNum;
	DofLayout		*rowDofs, *colDofs;
	unsigned		nRowEls;
	unsigned		nRowDofs, nColDofs;
	int			*rowEqs, *colEqs;
	int			rowEq, colEq;
	PetscScalar		*xVals, *yVals;
	double**		elStiffMat;
	PetscErrorCode		ec;
	unsigned		e_i, dof_i, dof_j;

	assert( A );
	assert( x );
	assert( y );

	ec = MatShellGetContext( A, (void*)&self );
	CheckPETScError( ec );
	assert( self );

	ec = VecGetArray( x, &xVals );
	CheckPETScError( ec );

	stiffMat = self->stiffMat;
	sle = self->sle;
	rowVar = stiffMat->rowVariable;
	colVar = stiffMat->columnVariable ? stiffMat->columnVariable : rowVar;
	rowEqNum = rowVar->eqNum;
	colEqNum = colVar->eqNum;
	rowMesh = rowVar->feMesh;
	rowDofs = rowVar->dofLayout;
	colDofs = colVar->dofLayout;
	nRowEls = FeMesh_GetElementLocalSize( rowMesh );			assert( nRowEls );
	nRowDofs = self->nRowDofs;						assert( nRowDofs );
	nColDofs = self->nColDofs;						assert( nColDofs );
	elStiffMat = self->elStiffMat;						assert( elStiffMat );

	yVals = AllocArray( double, nRowDofs );

	for( e_i = 0; e_i < nRowEls; e_i++ ) {
		rowEqs = rowEqNum->locationMatrix[e_i][0];
		colEqs = colEqNum->locationMatrix[e_i][0];

		memset( yVals, 0, nRowDofs * sizeof(double) );

		for( dof_i = 0; dof_i < nRowDofs; dof_i++ ) {
			yVals[dof_i] = 0.0;

			rowEq = rowEqs[dof_i];
			if( rowEq == -1 )
				continue;

			for( dof_j = 0; dof_j < nColDofs; dof_j++ ) {
				colEq = colEqs[dof_j];
				if( colEq == -1 )
					continue;

				yVals[dof_i] += elStiffMat[dof_i][dof_j] * xVals[colEq];
			}
		}

		ec = VecSetValues( y, (PetscInt)nRowDofs, (PetscInt*)rowEqs, (PetscScalar*)yVals, ADD_VALUES );
		CheckPETScError( ec );
	}

	FreeArray( yVals );

	ec = VecRestoreArray( x, &xVals );
	CheckPETScError( ec );

	ec = VecAssemblyEnd( y );
	CheckPETScError( ec );

	return ec;
}

PetscErrorCode PETScShellMatrix_MatMultTranspose( Mat A, Vec x, Vec y ) {
	PETScShellMatrix*	self;
	StiffnessMatrix*	stiffMat;
	SystemLinearEquations*	sle;
	FeVariable		*rowVar, *colVar;
	FeMesh			*rowMesh;
	FeEquationNumber	*rowEqNum, *colEqNum;
	DofLayout		*rowDofs, *colDofs;
	unsigned		nRowEls;
	unsigned		nRowDofs, nColDofs;
	int			*rowEqs, *colEqs;
	int			rowEq, colEq;
	PetscScalar		*xVals, *yVals;
	double**		elStiffMat;
	PetscErrorCode		ec;
	unsigned		e_i, dof_i, dof_j;

	assert( A );
	assert( x );
	assert( y );

	ec = MatShellGetContext( A, (void*)&self );
	CheckPETScError( ec );
	assert( self );

	ec = VecSet( y, 0.0 );
	CheckPETScError( ec );

	ec = VecGetArray( x, &xVals );
	CheckPETScError( ec );

	stiffMat = self->stiffMat;
	sle = self->sle;
	colVar = stiffMat->rowVariable;
	rowVar = stiffMat->columnVariable ? stiffMat->columnVariable : colVar;
	rowEqNum = rowVar->eqNum;
	colEqNum = colVar->eqNum;
	rowMesh = rowVar->feMesh;
	rowDofs = rowVar->dofLayout;
	colDofs = colVar->dofLayout;
	nRowEls = FeMesh_GetElementLocalSize( rowMesh );			assert( nRowEls );
	nRowDofs = self->nColDofs;						assert( nRowDofs );
	nColDofs = self->nRowDofs;						assert( nColDofs );
	elStiffMat = self->elStiffMatTrans;					assert( elStiffMat );

	yVals = AllocArray( double, nRowDofs );

	for( e_i = 0; e_i < nRowEls; e_i++ ) {
		rowEqs = rowEqNum->locationMatrix[e_i][0];
		colEqs = colEqNum->locationMatrix[e_i][0];

		memset( yVals, 0, nRowDofs * sizeof(double) );

		for( dof_i = 0; dof_i < nRowDofs; dof_i++ ) {
			yVals[dof_i] = 0.0;

			rowEq = rowEqs[dof_i];
			if( rowEq == -1 )
				continue;

			for( dof_j = 0; dof_j < nColDofs; dof_j++ ) {
				colEq = colEqs[dof_j];
				if( colEq == -1 )
					continue;

				yVals[dof_i] += elStiffMat[dof_i][dof_j] * xVals[colEq];
			}
		}

		ec = VecSetValues( y, (PetscInt)nRowDofs, (PetscInt*)rowEqs, (PetscScalar*)yVals, ADD_VALUES );
		CheckPETScError( ec );
	}

	FreeArray( yVals );

	ec = VecRestoreArray( x, &xVals );
	CheckPETScError( ec );

	ec = VecAssemblyEnd( y );
	CheckPETScError( ec );

	return ec;
}

PetscErrorCode PETScShellMatrix_MatGetDiagonal( Mat A, Vec x ) {
	PETScShellMatrix*	self;
	StiffnessMatrix*	stiffMat;
	SystemLinearEquations*	sle;
	FeVariable		*rowVar, *colVar;
	FeMesh			*rowMesh;
	FeEquationNumber	*rowEqNum, *colEqNum;
	DofLayout		*rowDofs, *colDofs;
	unsigned		nRowEls;
	unsigned		nRowDofs, nColDofs;
	int			*rowEqs, *colEqs;
	int			rowEq, colEq;
	PetscScalar		*xVals;
	double**		elStiffMat;
	PetscErrorCode		ec;
	unsigned		e_i, dof_i, dof_j;

	assert( A );
	assert( x );

	ec = MatShellGetContext( A, (void*)&self );
	CheckPETScError( ec );
	assert( self );

	ec = VecSet( x, 0.0 );
	CheckPETScError( ec );

	stiffMat = self->stiffMat;
	sle = self->sle;
	rowVar = stiffMat->rowVariable;
	colVar = stiffMat->columnVariable ? stiffMat->columnVariable : rowVar;
	rowEqNum = rowVar->eqNum;
	colEqNum = colVar->eqNum;
	rowMesh = rowVar->feMesh;
	rowDofs = rowVar->dofLayout;
	colDofs = colVar->dofLayout;
	nRowEls = FeMesh_GetElementLocalSize( rowMesh );			assert( nRowEls );
	nRowDofs = self->nRowDofs;						assert( nRowDofs );
	nColDofs = self->nColDofs;						assert( nColDofs );
	elStiffMat = self->elStiffMat;						assert( elStiffMat );

	xVals = AllocArray( double, nRowDofs );

	for( e_i = 0; e_i < nRowEls; e_i++ ) {
		rowEqs = rowEqNum->locationMatrix[e_i][0];
		colEqs = colEqNum->locationMatrix[e_i][0];

		for( dof_i = 0; dof_i < nRowDofs; dof_i++ ) {
			rowEq = rowEqs[dof_i];
			if( rowEq == -1 )
				continue;

			for( dof_j = 0; dof_j < nColDofs; dof_j++ ) {
				colEq = colEqs[dof_j];
				if( colEq == -1 || colEq != rowEq )
					continue;

				xVals[dof_i] = elStiffMat[dof_i][dof_j];
				break;
			}
		}

		ec = VecSetValues( x, (PetscInt)nRowDofs, (PetscInt*)rowEqs, (PetscScalar*)xVals, ADD_VALUES );
		CheckPETScError( ec );
	}

	ec = VecAssemblyBegin( x );
	CheckPETScError( ec );

	return ec;
}
