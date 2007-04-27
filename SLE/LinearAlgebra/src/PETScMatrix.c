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
** $Id: PETScMatrix.c 3584 2006-05-16 11:11:07Z PatrickSunter $
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
const Type PETScMatrix_Type = "PETScMatrix";


/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/

PETScMatrix* PETScMatrix_New( Name name ) {
	return _PETScMatrix_New( sizeof(PETScMatrix), 
				 PETScMatrix_Type, 
				 _PETScMatrix_Delete, 
				 _PETScMatrix_Print, 
				 NULL, 
				 (void* (*)(Name))PETScMatrix_New, 
				 _PETScMatrix_Construct, 
				 _PETScMatrix_Build, 
				 _PETScMatrix_Initialise, 
				 _PETScMatrix_Execute, 
				 _PETScMatrix_Destroy, 
				 name, 
				 NON_GLOBAL, 
				 PETScMatrix_SetComm, 
				 PETScMatrix_SetGlobalSize, 
				 PETScMatrix_SetLocalSize, 
				 PETScMatrix_AddEntries, 
				 PETScMatrix_InsertEntries, 
				 PETScMatrix_DiagonalAddEntries, 
				 PETScMatrix_DiagonalInsertEntries, 
				 PETScMatrix_Zero, 
				 PETScMatrix_Load, 
				 PETScMatrix_AssemblyBegin, 
				 PETScMatrix_AssemblyEnd, 
				 PETScMatrix_Scale, 
				 PETScMatrix_AddScaled, 
				 PETScMatrix_DiagonalScale, 
				 PETScMatrix_Multiply, 
				 PETScMatrix_TransposeMultiply, 
				 PETScMatrix_MultiplyAdd, 
				 PETScMatrix_PAPt, 
				 PETScMatrix_PtAP, 
				 PETScMatrix_MatrixMultiply, 
				 PETScMatrix_L2Norm, 
				 PETScMatrix_Transpose, 
				 PETScMatrix_GetGlobalSize, 
				 PETScMatrix_GetLocalSize, 
				 PETScMatrix_GetRow, 
				 PETScMatrix_RestoreRow, 
				 PETScMatrix_GetDiagonal, 
				 PETScMatrix_Duplicate, 
				 PETScMatrix_CopyEntries, 
				 _PETScMatrix_SetNonZeroStructure );
}

PETScMatrix* _PETScMatrix_New( PETSCMATRIX_DEFARGS ) {
	PETScMatrix*	self;

	/* Allocate memory */
	assert( sizeOfSelf >= sizeof(PETScMatrix) );
	self = (PETScMatrix*)_Matrix_New( MATRIX_PASSARGS );

	/* Virtual info */
	self->setNonZeroStructure = setNonZeroStructure;

	/* PETScMatrix info */
	_PETScMatrix_Init( self );

	return self;
}

void _PETScMatrix_Init( PETScMatrix* self ) {
	assert( self && Stg_CheckType( self, PETScMatrix ) );

	self->petscMat = PETSC_NULL;
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _PETScMatrix_Delete( void* matrix ) {
	PETScMatrix*	self = (PETScMatrix*)matrix;

	assert( self && Stg_CheckType( self, PETScMatrix ) );

	if( self->petscMat != PETSC_NULL ) {
		PetscErrorCode	ec;

		ec = MatDestroy( self->petscMat );
		CheckPETScError( ec );
	}

	/* Delete the parent. */
	_Matrix_Delete( self );
}

void _PETScMatrix_Print( void* matrix, Stream* stream ) {
	PETScMatrix*	self = (PETScMatrix*)matrix;
	
	/* Set the Journal for printing informations */
	Stream* matrixStream;
	matrixStream = Journal_Register( InfoStream_Type, "PETScMatrixStream" );

	assert( self && Stg_CheckType( self, PETScMatrix ) );

	/* Print parent */
	Journal_Printf( stream, "PETScMatrix (ptr): (%p)\n", self );
	_Matrix_Print( self, stream );
}

void _PETScMatrix_Construct( void* matrix, Stg_ComponentFactory* cf, void* data ) {
	PETScMatrix*	self = (PETScMatrix*)matrix;

	assert( self && Stg_CheckType( self, PETScMatrix ) );
	assert( cf );
}

void _PETScMatrix_Build( void* matrix, void* data ) {
}

void _PETScMatrix_Initialise( void* matrix, void* data ) {
}

void _PETScMatrix_Execute( void* matrix, void* data ) {
}

void _PETScMatrix_Destroy( void* matrix, void* data ) {
}

void PETScMatrix_SetComm( void* matrix, MPI_Comm comm ) {
	PETScMatrix*	self = (PETScMatrix*)matrix;
	PetscErrorCode	ec;

	assert( self && Stg_CheckType( self, PETScMatrix ) );

	_Matrix_SetComm( self, comm );

	if( self->petscMat != PETSC_NULL )
		MatDestroy( self->petscMat );
	ec = MatCreate( self->comm, &self->petscMat );
	CheckPETScError( ec );

	self->hasChanged = True;
}

void PETScMatrix_SetGlobalSize( void* matrix, unsigned nRows, unsigned nColumns ) {
	PETScMatrix*	self = (PETScMatrix*)matrix;
	PetscErrorCode	ec;

	assert( self && Stg_CheckType( self, PETScMatrix ) );

	ec = MatSetSizes( self->petscMat, PETSC_DETERMINE, PETSC_DETERMINE, (PetscInt)nRows, (PetscInt)nColumns );
	CheckPETScError( ec );
	ec = MatSetFromOptions( self->petscMat );
	CheckPETScError( ec );

	self->hasChanged = True;
}

void PETScMatrix_SetLocalSize( void* matrix, unsigned nRows, unsigned nColumns ) {
	PETScMatrix*	self = (PETScMatrix*)matrix;
	PetscErrorCode	ec;

	assert( self && Stg_CheckType( self, PETScMatrix ) );

	ec = MatSetSizes( self->petscMat, (PetscInt)nRows, (PetscInt)nColumns, PETSC_DETERMINE, PETSC_DETERMINE );
	CheckPETScError( ec );
	ec = MatSetFromOptions( self->petscMat );
	CheckPETScError( ec );

	self->hasChanged = True;
}

void PETScMatrix_AddEntries( void* matrix, unsigned nRows, unsigned* rowIndices, 
			     unsigned nColumns, unsigned* columnIndices, 
			     double* values )
{
	PETScMatrix*	self = (PETScMatrix*)matrix;
	PetscErrorCode	ec;

	assert( self && Stg_CheckType( self, PETScMatrix ) );
	assert( !nRows || rowIndices );
	assert( !nColumns || columnIndices );
	assert( (!nRows && !nColumns) || values );
	assert( sizeof(PetscInt) == sizeof(unsigned) );

	ec = MatSetValues( self->petscMat, (PetscInt)nRows, (PetscInt*)rowIndices, 
			   (PetscInt)nColumns, (PetscInt*)columnIndices, 
			   (PetscScalar*)values, ADD_VALUES );
	CheckPETScError( ec );

	self->hasChanged = True;
}

void PETScMatrix_InsertEntries( void* matrix, unsigned nRows, unsigned* rowIndices, 
				unsigned nColumns, unsigned* columnIndices, 
				double* values )
{
	PETScMatrix*	self = (PETScMatrix*)matrix;
	PetscErrorCode	ec;

	assert( self && Stg_CheckType( self, PETScMatrix ) );
	assert( !nRows || rowIndices );
	assert( !nColumns || columnIndices );
	assert( (!nRows && !nColumns) || values );

	ec = MatSetValues( self->petscMat, (PetscInt)nRows, (PetscInt*)rowIndices, 
			   (PetscInt)nColumns, (PetscInt*)columnIndices, 
			   (PetscScalar*)values, INSERT_VALUES );
	CheckPETScError( ec );

	self->hasChanged = True;
}

void PETScMatrix_DiagonalAddEntries( void* matrix, void* _vector ) {
	PETScMatrix*	self = (PETScMatrix*)matrix;
	PETScVector*	vector = (PETScVector*)_vector;
	PetscErrorCode	ec;

	assert( self && Stg_CheckType( self, PETScMatrix ) );
	assert( vector && Stg_CheckType( vector, PETScVector ) );

	ec = MatDiagonalSet( self->petscMat, vector->petscVec, ADD_VALUES );
	CheckPETScError( ec );

	self->hasChanged = True;
}

void PETScMatrix_DiagonalInsertEntries( void* matrix, void* _vector ) {
	PETScMatrix*	self = (PETScMatrix*)matrix;
	PETScVector*	vector = (PETScVector*)_vector;
	PetscErrorCode	ec;

	assert( self && Stg_CheckType( self, PETScMatrix ) );
	assert( vector && Stg_CheckType( vector, PETScVector ) );

	ec = MatDiagonalSet( self->petscMat, vector->petscVec, INSERT_VALUES );
	CheckPETScError( ec );

	self->hasChanged = True;
}

void PETScMatrix_Zero( void* matrix ) {
	PETScMatrix*	self = (PETScMatrix*)matrix;
	PetscErrorCode	ec;

	assert( self && Stg_CheckType( self, PETScMatrix ) );

	ec = MatZeroEntries( self->petscMat );
	CheckPETScError( ec );

	self->hasChanged = True;
}

void PETScMatrix_Load( void* matrix, char* filename ) {
	PETScMatrix*	self = (PETScMatrix*)matrix;
	PetscViewer	viewer;
	PetscErrorCode	ec;

	assert( self && Stg_CheckType( self, PETScMatrix ) );

	ec = PetscViewerBinaryOpen( PETSC_COMM_SELF, filename, FILE_MODE_READ, &viewer );
	CheckPETScError( ec );
	if( self->petscMat != PETSC_NULL )
		MatDestroy( self->petscMat );
	ec = MatLoad( viewer, MATSEQAIJ, &self->petscMat );
	CheckPETScError( ec );
	ec = PetscViewerDestroy( viewer );
	CheckPETScError( ec );

	self->hasChanged = True;
}

void PETScMatrix_AssemblyBegin( void* matrix ) {
	PETScMatrix*	self = (PETScMatrix*)matrix;
	PetscErrorCode	ec;

	assert( self && Stg_CheckType( self, PETScMatrix ) );

	_Matrix_AssemblyBegin( self );

	ec = MatAssemblyBegin( self->petscMat, MAT_FINAL_ASSEMBLY );
	CheckPETScError( ec );
}

void PETScMatrix_AssemblyEnd( void* matrix ) {
	PETScMatrix*	self = (PETScMatrix*)matrix;
	PetscErrorCode	ec;

	assert( self && Stg_CheckType( self, PETScMatrix ) );

	ec = MatAssemblyEnd( self->petscMat, MAT_FINAL_ASSEMBLY );
	CheckPETScError( ec );
}

void PETScMatrix_Scale( void* matrix, double factor ) {
	PETScMatrix*	self = (PETScMatrix*)matrix;
	PetscErrorCode	ec;

	assert( self && Stg_CheckType( self, PETScMatrix ) );

	ec = MatScale( self->petscMat, (PetscScalar)factor );
	CheckPETScError( ec );

	self->hasChanged = True;
}

void PETScMatrix_AddScaled( void* matrix, double factor, void* _matrix0 ) {
	PETScMatrix*	self = (PETScMatrix*)matrix;
	PETScMatrix*	matrix0 = (PETScMatrix*)_matrix0;
	PetscErrorCode	ec;

	assert( self && Stg_CheckType( self, PETScMatrix ) );
	assert( matrix0 && Stg_CheckType( matrix0, PETScMatrix ) );

	ec = MatAXPY( self->petscMat, (PetscScalar)factor, matrix0->petscMat, SAME_NONZERO_PATTERN );
	CheckPETScError( ec );

	self->hasChanged = True;
}

void PETScMatrix_DiagonalScale( void* matrix, void* _leftVector, void* _rightVector ) {
	PETScMatrix*	self = (PETScMatrix*)matrix;
	PETScVector*	leftVector = (PETScVector*)_leftVector;
	PETScVector*	rightVector = (PETScVector*)_rightVector;
	PetscErrorCode	ec;

	assert( self && Stg_CheckType( self, PETScMatrix ) );
	assert( leftVector ? (leftVector && Stg_CheckType( leftVector, PETScVector )) : (rightVector != NULL) );
	assert( rightVector ? (rightVector && Stg_CheckType( rightVector, PETScVector )) : (leftVector != NULL) );

	ec = MatDiagonalScale( self->petscMat, leftVector ? leftVector->petscVec : PETSC_NULL, 
			       rightVector ? rightVector->petscVec : PETSC_NULL );
	CheckPETScError( ec );

	self->hasChanged = True;
}

void PETScMatrix_Multiply( void* matrix, void* _vector, void* _dstVector ) {
	PETScMatrix*	self = (PETScMatrix*)matrix;
	PETScVector*	vector = (PETScVector*)_vector;
	PETScVector*	dstVector = (PETScVector*)_dstVector;
	PetscErrorCode	ec;

	assert( self && Stg_CheckType( self, PETScMatrix ) );
	assert( vector && Stg_CheckType( vector, PETScVector ) );
	assert( dstVector && Stg_CheckType( dstVector, PETScVector ) );

	ec = MatMult( self->petscMat, vector->petscVec, dstVector->petscVec );
	CheckPETScError( ec );
}

void PETScMatrix_TransposeMultiply( void* matrix, void* _vector, void* _dstVector ) {
	PETScMatrix*	self = (PETScMatrix*)matrix;
	PETScVector*	vector = (PETScVector*)_vector;
	PETScVector*	dstVector = (PETScVector*)_dstVector;
	PetscErrorCode	ec;

	assert( self && Stg_CheckType( self, PETScMatrix ) );
	assert( vector && Stg_CheckType( vector, PETScVector ) );
	assert( dstVector && Stg_CheckType( dstVector, PETScVector ) );

	ec = MatMultTranspose( self->petscMat, vector->petscVec, dstVector->petscVec );
	CheckPETScError( ec );
}

void PETScMatrix_MultiplyAdd( void* matrix, void* _vector0, void* _vector1, void* _dstVector ) {
	PETScMatrix*	self = (PETScMatrix*)matrix;
	PETScVector*	vector0 = (PETScVector*)_vector0;
	PETScVector*	vector1 = (PETScVector*)_vector1;
	PETScVector*	dstVector = (PETScVector*)_dstVector;
	PetscErrorCode	ec;

	assert( self && Stg_CheckType( self, PETScMatrix ) );
	assert( vector0 && Stg_CheckType( vector0, PETScVector ) );
	assert( vector1 && Stg_CheckType( vector1, PETScVector ) );
	assert( dstVector && Stg_CheckType( dstVector, PETScVector ) );

	ec = MatMultAdd( self->petscMat, vector0->petscVec, vector1->petscVec, dstVector->petscVec );
	CheckPETScError( ec );
}

void PETScMatrix_PAPt( void* matrix, void* _P, void** _dstMatrix ) {
#if 0
	PETScMatrix*	self = (PETScMatrix*)matrix;
	PETScMatrix*	P = (PETScMatrix*)_P;
	PETScMatrix*	dstMatrix = (PETScMatrix*)_dstMatrix;
	PetscErrorCode	ec;
	double		fillRatio;
	double		nRowsA, nRowsC;
	double		nzA, nzP, nzC;
	MatInfo		mInfo;

	assert( self && Stg_CheckType( self, PETScMatrix ) );
	assert( P && Stg_CheckType( P, PETScMatrix ) );
	assert( dstMatrix && Stg_CheckType( dstMatrix, PETScMatrix ) );

	MatGetInfo( self->petscMat, MAT_GLOBAL_SUM, &mInfo );
	nRowsA = mInfo.rows_global;
	nzA = mInfo.nz_used;

	MatGetInfo( P->petscMat, MAT_GLOBAL_SUM, &mInfo );
	nRowsC = mInfo.columns_global;
	nzP = mInfo.nz_used;

	nzC = (nRowsC / nRowsA) * nzA;
	fillRatio = nzC / (nzA + nzP);

	assert( 0 );
	/*ec = MatPAPt( self->petscMat, P->petscMat, MAT_REUSE_MATRIX, (PetscScalar)fillRatio, &dstMatrix->petscMat );*/
	CheckPETScError( ec );

	dstMatrix->hasChanged = True;
#endif
}

void PETScMatrix_PtAP( void* matrix, void* _P, void** _dstMatrix ) {
	PETScMatrix*	self = (PETScMatrix*)matrix;
	PETScMatrix*	P = (PETScMatrix*)_P;
	PETScMatrix**	dstMatrix = (PETScMatrix**)_dstMatrix;
	PetscErrorCode	ec;
#if 0
	double		fillRatio;
	double		nRowsA, nRowsC;
	double		nzA, nzP, nzC;
	MatInfo		mInfo;
#endif

	assert( self && Stg_CheckType( self, PETScMatrix ) );
	assert( P && Stg_CheckType( P, PETScMatrix ) );
	assert( dstMatrix && (!(*dstMatrix) || Stg_CheckType( *dstMatrix, PETScMatrix )) );

#if 0
	MatGetInfo( self->petscMat, MAT_GLOBAL_SUM, &mInfo );
	nRowsA = mInfo.rows_global;
	nzA = mInfo.nz_used;

	MatGetInfo( P->petscMat, MAT_GLOBAL_SUM, &mInfo );
	nRowsC = mInfo.columns_global;
	nzP = mInfo.nz_used;

	nzC = (nRowsC / nRowsA) * nzA;
	fillRatio = nzC / (nzA + nzP);
#endif

	if( *dstMatrix )
		ec = MatPtAP( self->petscMat, P->petscMat, MAT_REUSE_MATRIX, 1.0, &(*dstMatrix)->petscMat );
	else {
		Matrix_Duplicate( self, (void**)dstMatrix );
		if( (*dstMatrix)->petscMat )
			MatDestroy( (*dstMatrix)->petscMat );
		ec = MatPtAP( self->petscMat, P->petscMat, MAT_INITIAL_MATRIX, 1.0, &(*dstMatrix)->petscMat );
	}
	CheckPETScError( ec );

	(*dstMatrix)->hasChanged = True;
}

void PETScMatrix_MatrixMultiply( void* matrix, void* _matrix0, void* _dstMatrix ) {
	PETScMatrix*	self = (PETScMatrix*)matrix;
	PETScMatrix*	matrix0 = (PETScMatrix*)_matrix0;
	PETScMatrix*	dstMatrix = (PETScMatrix*)_dstMatrix;
	PetscErrorCode	ec;
#if 0
	double		fillRatio;
	double		nRowsA, nRowsC;
	double		nzA, nzP, nzC;
	MatInfo		mInfo;
#endif

	assert( self && Stg_CheckType( self, PETScMatrix ) );
	assert( matrix0 && Stg_CheckType( matrix0, PETScMatrix ) );
	assert( dstMatrix && Stg_CheckType( dstMatrix, PETScMatrix ) );

#if 0
	MatGetInfo( self->petscMat, MAT_GLOBAL_SUM, &mInfo );
	nRowsA = mInfo.rows_global;
	nzA = mInfo.nz_used;

	MatGetInfo( matrix0->petscMat, MAT_GLOBAL_SUM, &mInfo );
	nRowsC = mInfo.columns_global;
	nzP = mInfo.nz_used;

	nzC = (nRowsC / nRowsA) * nzA;
	fillRatio = nzC / (nzA + nzP);
#endif

	if( dstMatrix->petscMat != PETSC_NULL )
		ec = MatDestroy( dstMatrix->petscMat );
	CheckPETScError( ec );
	ec = MatPtAP( self->petscMat, matrix0->petscMat, MAT_INITIAL_MATRIX, 1.0, &dstMatrix->petscMat );
	CheckPETScError( ec );

	dstMatrix->hasChanged = True;
}

double PETScMatrix_L2Norm( void* matrix ) {
	PETScMatrix*	self = (PETScMatrix*)matrix;
	PetscErrorCode	ec;
	PetscScalar	norm;

	assert( self && Stg_CheckType( self, PETScMatrix ) );

	ec = MatNorm( self->petscMat, NORM_FROBENIUS, &norm );
	CheckPETScError( ec );

	return (double)norm;
}

void PETScMatrix_Transpose( void* matrix, void* _dstMatrix ) {
	PETScMatrix*	self = (PETScMatrix*)matrix;
	PetscErrorCode	ec;
	PETScMatrix*	dstMatrix = (PETScMatrix*)_dstMatrix;

	assert( self && Stg_CheckType( self, PETScMatrix ) );
	assert( dstMatrix && Stg_CheckType( dstMatrix, PETScMatrix ) );

	if( self != dstMatrix ) {
		if( dstMatrix->petscMat != PETSC_NULL )
			ec = MatDestroy( dstMatrix->petscMat );
		CheckPETScError( ec );
		ec = MatTranspose( self->petscMat, &dstMatrix->petscMat );
	}
	else
		ec = MatTranspose( self->petscMat, PETSC_NULL );
	CheckPETScError( ec );

	dstMatrix->hasChanged = True;
}

void PETScMatrix_GetGlobalSize( void* matrix, unsigned* nRows, unsigned* nColumns ) {
	PETScMatrix*	self = (PETScMatrix*)matrix;
	PetscErrorCode	ec;
	PetscInt	nr, nc;

	assert( self && Stg_CheckType( self, PETScMatrix ) );
	assert( nRows );
	assert( nColumns );

	ec = MatGetSize( self->petscMat, &nr, &nc );
	CheckPETScError( ec );

	*nRows = (unsigned)nr;
	*nColumns = (unsigned)nc;
}

void PETScMatrix_GetLocalSize( void* matrix, unsigned* nRows, unsigned* nColumns ) {
	PETScMatrix*	self = (PETScMatrix*)matrix;
	PetscErrorCode	ec;
	PetscInt	nr, nc;

	assert( self && Stg_CheckType( self, PETScMatrix ) );
	assert( nRows || nColumns );

	ec = MatGetLocalSize( self->petscMat, &nr, &nc );
	CheckPETScError( ec );

	if( nRows )
		*nRows = (unsigned)nr;
	if( nColumns )
		*nColumns = (unsigned)nc;
}

void PETScMatrix_GetRow( void* matrix, unsigned rowIndex, 
			 unsigned* nEntries, unsigned** columnIndices, double** values )
{
	PETScMatrix*	self = (PETScMatrix*)matrix;
	PetscErrorCode	ec;

	assert( self && Stg_CheckType( self, PETScMatrix ) );
	assert( nEntries );
	assert( columnIndices );
	assert( values );
	assert( sizeof(PetscScalar) == sizeof(double) );
	assert( sizeof(PetscInt) == sizeof(unsigned) );

	ec = MatGetRow( self->petscMat, (PetscInt)rowIndex, 
			(PetscInt*)nEntries, (const PetscInt**)columnIndices, (const PetscScalar**)values );
	CheckPETScError( ec );
}

void PETScMatrix_RestoreRow( void* matrix, unsigned rowIndex, 
			     unsigned* nEntries, unsigned** columnIndices, double** values )
{
	PETScMatrix*	self = (PETScMatrix*)matrix;
	PetscErrorCode	ec;

	assert( self && Stg_CheckType( self, PETScMatrix ) );
	assert( nEntries );
	assert( columnIndices );
	assert( values );

	ec = MatRestoreRow( self->petscMat, (PetscInt)rowIndex, 
			    (PetscInt*)nEntries, (const PetscInt**)columnIndices, (const PetscScalar**)values );
	CheckPETScError( ec );
}

void PETScMatrix_GetDiagonal( void* matrix, void* _vector ) {
	PETScMatrix*	self = (PETScMatrix*)matrix;
	PETScVector*	vector = (PETScVector*)_vector;
	PetscErrorCode	ec;

	assert( self && Stg_CheckType( self, PETScMatrix ) );
	assert( vector && Stg_CheckType( vector, PETScVector ) );

	ec = MatGetDiagonal( self->petscMat, vector->petscVec );
	CheckPETScError( ec );
}

void PETScMatrix_Duplicate( void* matrix, void** _dstMatrix ) {
	PETScMatrix*	self = (PETScMatrix*)matrix;
	PETScMatrix**	dstMatrix = (PETScMatrix**)_dstMatrix;
	PetscErrorCode	ec;

	assert( self && Stg_CheckType( self, PETScMatrix ) );
	assert( dstMatrix );

	*dstMatrix = self->_defaultConstructor( "" );
	ec = MatCreate( self->comm, &(*dstMatrix)->petscMat );
	CheckPETScError( ec );

	(*dstMatrix)->hasChanged = True;
}

void PETScMatrix_CopyEntries( void* matrix, void* _dstMatrix ) {
	PETScMatrix*	self = (PETScMatrix*)matrix;
	PETScMatrix*	dstMatrix = (PETScMatrix*)dstMatrix;
	PetscErrorCode	ec;

	assert( self && Stg_CheckType( self, PETScMatrix ) );
	assert( dstMatrix && Stg_CheckType( dstMatrix, PETScMatrix ) );

	ec = MatCopy( self->petscMat, dstMatrix->petscMat, DIFFERENT_NONZERO_PATTERN );
	CheckPETScError( ec );

	dstMatrix->hasChanged = True;
}

void _PETScMatrix_SetNonZeroStructure( void* matrix, unsigned nNonZeros, 
				       unsigned* diagonalNonZeroIndices, unsigned* offDiagonalNonZeroIndices )
{
	PETScMatrix*	self = (PETScMatrix*)matrix;
	PetscErrorCode	ec;
	unsigned	nProcs;

	assert( self && Stg_CheckType( self, PETScMatrix ) );

	MPI_Comm_size( self->comm, (int*)&nProcs );
	if( diagonalNonZeroIndices || offDiagonalNonZeroIndices ) {
		if( nProcs > 1 ) {
			ec = MatMPIAIJSetPreallocation( self->petscMat, 
							PETSC_NULL, (int*)diagonalNonZeroIndices,
							PETSC_NULL, (int*)offDiagonalNonZeroIndices );
		}
		else
			ec = MatSeqAIJSetPreallocation( self->petscMat, PETSC_NULL, (int*)diagonalNonZeroIndices );
	}
	else {
		if( nProcs > 1 )
			ec = MatMPIAIJSetPreallocation( self->petscMat, nNonZeros, PETSC_NULL, nNonZeros, PETSC_NULL );
		else
			ec = MatSeqAIJSetPreallocation( self->petscMat, nNonZeros, PETSC_NULL );
	}
	CheckPETScError( ec );

	self->hasChanged = True;
}


/*--------------------------------------------------------------------------------------------------------------------------
** Public Functions
*/

void PETScMatrix_Draw( void* matrix ) {
	PETScMatrix*	self = (PETScMatrix*)matrix;
	PetscErrorCode	ec;

	assert( self && Stg_CheckType( self, PETScMatrix ) );

	ec = MatView( self->petscMat, PETSC_VIEWER_DRAW_WORLD );
	CheckPETScError( ec );
}

void PETScMatrix_View( void* matrix ) {
	PETScMatrix*	self = (PETScMatrix*)matrix;
	PetscErrorCode	ec;

	assert( self && Stg_CheckType( self, PETScMatrix ) );

	ec = MatView( self->petscMat, PETSC_VIEWER_STDOUT_SELF );
	CheckPETScError( ec );
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Private Functions
*/
