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
** $Id: Matrix.c 744 2007-02-15 00:34:22Z DavidMay $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <assert.h>
#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgFEM/Discretisation/Discretisation.h>
#include "SLE/LinearAlgebra/LinearAlgebra.h"

#include <petsc.h>
#include <petscmat.h>
#include <StGermain/compatibility/petsccompat.h>

#include "ErrorChecking.h"

Matrix* Matrix_New( 
		MPI_Comm			comm,
		Index				rowLocalSize,
		Index				colLocalSize,
		Index				nonZeroCount )
{
	Mat                matrix;
	PetscErrorCode     errorFlag;
	Processor_Index    numProcs;

	MPI_Comm_size( comm, (int*)(&numProcs) );
	
	errorFlag = MatCreate( comm, rowLocalSize, colLocalSize, PETSC_DETERMINE, PETSC_DETERMINE, &matrix );
	CheckPETScError( errorFlag );

	errorFlag = MatSetFromOptions( matrix ); 
	CheckPETScError( errorFlag );
	
	if ( numProcs > 1 ) {
		errorFlag = MatMPIAIJSetPreallocation( matrix, nonZeroCount, PETSC_NULL, nonZeroCount, PETSC_NULL );
	}
	else {
		errorFlag = MatSeqAIJSetPreallocation( matrix, nonZeroCount, PETSC_NULL );
	}
	CheckPETScError( errorFlag );
	
	return (void*)matrix;
}


Matrix* Matrix_NewEmpty( MPI_Comm comm, unsigned rowLocalSize, unsigned colLocalSize ) {
	Mat                mat;
	PetscErrorCode     ef;
	
	ef = MatCreate( comm, rowLocalSize, colLocalSize, PETSC_DETERMINE, PETSC_DETERMINE, &mat );	CheckPETScError( ef );
	ef = MatSetFromOptions( mat );													CheckPETScError( ef );
	
	return (Matrix*)mat;
}


Matrix* Matrix_New2( 
		MPI_Comm			comm,
		Index				rowLocalSize,
		Index				colLocalSize, 
		Index*				diagonalNonZeroIndices,
		Index*				offDiagonalNonZeroIndices )
{
	Mat                matrix;
	PetscErrorCode     errorFlag;
	Processor_Index    numProcs;
	
	MPI_Comm_size( comm, (int*)(&numProcs) );

	errorFlag = MatCreate( comm, rowLocalSize, colLocalSize, PETSC_DETERMINE, PETSC_DETERMINE, &matrix );
	CheckPETScError( errorFlag );

	errorFlag = MatSetFromOptions( matrix ); 
	CheckPETScError( errorFlag );

	if ( numProcs > 1 ) {
		errorFlag = MatMPIAIJSetPreallocation( matrix,
			PETSC_NULL, (int*)diagonalNonZeroIndices,
			PETSC_NULL, (int*)offDiagonalNonZeroIndices );
	}
	else {
		errorFlag = MatSeqAIJSetPreallocation( matrix, PETSC_NULL, (int*)diagonalNonZeroIndices );
	}
	CheckPETScError( errorFlag );
	
	return (void*)matrix;
}

Matrix* Matrix_NewFromFile( Name filename ) {
	PetscErrorCode     errorFlag;
	PetscViewer        viewer;
	Matrix*            newMatrix;

	errorFlag = PetscViewerBinaryOpen( PETSC_COMM_SELF, filename, PETSC_FILE_RDONLY, &viewer );
	CheckPETScError( errorFlag );

	errorFlag = MatLoad( viewer, MATSEQAIJ, (Mat*)&newMatrix );
	CheckPETScError( errorFlag );
	
	errorFlag = PetscViewerDestroy( viewer );
	CheckPETScError( errorFlag );

	return newMatrix;
}

void Matrix_Dump( Matrix* matrix, Name filename ) {
	PetscErrorCode     errorFlag;
	PetscViewer        viewer;

	errorFlag = PetscViewerBinaryOpen( PETSC_COMM_SELF, filename, PETSC_FILE_CREATE, &viewer );
	CheckPETScError( errorFlag );

	errorFlag = MatView( (Mat) matrix, viewer ); 
	CheckPETScError( errorFlag );

	errorFlag = PetscViewerDestroy( viewer );
	CheckPETScError( errorFlag );
}
	
void Matrix_View( void* matrix, Stream* stream ) {
	PetscErrorCode     errorFlag;

	if ( ! Stream_IsEnable( stream ) )
		return;
	
	errorFlag = PetscViewerSetFormat( PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_DENSE );
	CheckPETScError( errorFlag );

	errorFlag = MatView( matrix, PETSC_VIEWER_STDOUT_WORLD ); 
	CheckPETScError( errorFlag );
}

void Matrix_Zero( void* matrix ) {
	PetscErrorCode     errorFlag;
	
	errorFlag = MatZeroEntries( matrix ); 
	CheckPETScError( errorFlag );
}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
#ifdef PETSC_USE_COMPLEX

void Matrix_AddTo( void *matrix, Index rowCount, Index rowIndices[], Index colCount, Index colIndices[], double *values )
{
	PetscInt i;
	PetscScalar *_values;
	PetscInt size;
	
	size = rowCount*colCount;
	
	PetscMalloc( size*sizeof(PetscScalar), &_values );
	PetscMemzero( _values, size*sizeof(PetscScalar) );
	
	for( i=0; i<rowCount*colCount; i++ ) {
		PetscRealPart(_values[i]) = (PetscReal)values[i];
	}
/*	PetscScalarView( rowCount*colCount, _values, PETSC_VIEWER_STDOUT_WORLD ); */
	
	MatSetValues( matrix, rowCount, (int*)rowIndices, colCount, (int*)colIndices, _values, ADD_VALUES ); 
	PetscFree( _values );
}

#else

void Matrix_AddTo( 
	void*				matrix,
	Index				rowCount,
	Index				rowIndices[],
	Index				colCount,
	Index				colIndices[],
	double*				values )
{
	PetscErrorCode     errorFlag;
	
	errorFlag = MatSetValues( matrix, rowCount, (int*)rowIndices, colCount, (int*)colIndices, (PetscScalar*)values, ADD_VALUES ); 
	#ifdef DEBUG
		CheckPETScError( errorFlag );
	#endif
}

#endif


/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
#ifdef PETSC_USE_COMPLEX

void Matrix_Insert( void *matrix, Index rowCount, Index rowIndices[], Index colCount, Index colIndices[], double *values )
{
	PetscInt i;
	PetscScalar *_values;
	PetscInt size;
	
	size = rowCount*colCount;
	PetscMalloc( size*sizeof(PetscScalar), &_values );
	PetscMemzero( _values, size*sizeof(PetscScalar) );
	
	for( i=0; i<size; i++ ) {
		PetscRealPart(_values[i]) = (PetscReal)values[i];
	}
/*	PetscScalarView( rowCount*colCount, _values, PETSC_VIEWER_STDOUT_WORLD ); */
	MatSetValues( matrix, rowCount, (int*)rowIndices, colCount, (int*)colIndices, (PetscScalar*)values, INSERT_VALUES ); 
	
	PetscFree( _values );
}

#else

void Matrix_Insert( 
	void*				matrix,
	Index				rowCount,
	Index				rowIndices[],
	Index				colCount,
	Index				colIndices[],
	double*				values )
{
	PetscErrorCode     errorFlag;
	
	errorFlag = MatSetValues( matrix, rowCount, (int*)rowIndices, colCount, (int*)colIndices, (PetscScalar*)values, INSERT_VALUES ); 
	#ifdef DEBUG
		CheckPETScError( errorFlag );
	#endif
}

#endif



/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
#ifdef PETSC_USE_COMPLEX

void Matrix_SetValue( void* matrix, unsigned rowInd, unsigned colInd, double val ) 
{
	PetscScalar _val;
	
	PetscRealPart(_val) = (PetscReal)val;
	PetscImaginaryPart(_val) = 0.0;
	
	MatSetValue( (Mat)matrix, rowInd, colInd, _val, INSERT_VALUES );
}

#else

void Matrix_SetValue( void* matrix, unsigned rowInd, unsigned colInd, double val ) {
	PetscErrorCode	ec;

	ec = MatSetValue( (Mat)matrix, rowInd, colInd, val, INSERT_VALUES );
	CheckPETScError( ec );
}

#endif



/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
#ifdef PETSC_USE_COMPLEX

void Matrix_AddValue( void* matrix, unsigned rowInd, unsigned colInd, double val ) 
{
	PetscScalar _val;
	
	PetscRealPart(_val) = (PetscReal)val;
	PetscImaginaryPart(_val) = 0.0;
	MatSetValue( (Mat)matrix, rowInd, colInd, _val, ADD_VALUES );
}

#else

void Matrix_AddValue( void* matrix, unsigned rowInd, unsigned colInd, double val ) {
	PetscErrorCode	ec;

	ec = MatSetValue( (Mat)matrix, rowInd, colInd, val, ADD_VALUES );
	CheckPETScError( ec );
}

#endif

void Matrix_ZeroRow( void* matrix, unsigned rowInd ) {
	PetscErrorCode	ec;
	unsigned	nCols;
	unsigned*	colInds;
	double*		vals;
	unsigned	col_i;

	ec = MatGetSize( (Mat)matrix, PETSC_NULL, (PetscInt*)(&nCols) );	CheckPETScError( ec );

	colInds = Memory_Alloc_Array( unsigned, nCols, "Matrix_ZeroRow" );
	vals = Memory_Alloc_Array( double, nCols, "Matrix_ZeroRow" );
	for( col_i = 0; col_i < nCols; col_i++ ) {
		colInds[col_i] = col_i;
		vals[col_i] = 0.0;
	}

	ec = MatSetValues( (Mat)matrix, 1, (PetscInt*)(&rowInd), nCols, (PetscInt*)(colInds), (PetscScalar*)vals, INSERT_VALUES );	CheckPETScError( ec );

	Memory_Free( colInds );
	Memory_Free( vals );
}


void Matrix_ZeroCol( void* matrix, unsigned colInd ) {
	PetscErrorCode	ec;
	unsigned	nRows;
	unsigned*	rowInds;
	double*		vals;
	unsigned	row_i;

	ec = MatGetSize( (Mat)matrix, (PetscInt*)(&nRows), PETSC_NULL );	CheckPETScError( ec );

	rowInds = Memory_Alloc_Array( unsigned, nRows, "Matrix_ZeroCol" );
	vals = Memory_Alloc_Array( double, nRows, "Matrix_ZeroCol" );
	for( row_i = 0; row_i < nRows; row_i++ ) {
		rowInds[row_i] = row_i;
		vals[row_i] = 0.0;
	}

	ec = MatSetValues( (Mat)matrix, nRows, (PetscInt*)rowInds, 1, (PetscInt*)(&colInd), (PetscScalar*)vals, INSERT_VALUES );	CheckPETScError( ec );

	Memory_Free( rowInds );
	Memory_Free( vals );
}


void Matrix_AssemblyBegin( void* matrix ) {
	PetscErrorCode     errorFlag;
	
	errorFlag = MatAssemblyBegin( matrix, MAT_FINAL_ASSEMBLY );
	CheckPETScError( errorFlag );
}

void Matrix_AssemblyEnd( void* matrix ) {
	PetscErrorCode     errorFlag;

	errorFlag = MatAssemblyEnd( matrix, MAT_FINAL_ASSEMBLY );
	CheckPETScError( errorFlag );
}

void Matrix_DiagonalSet_Add( void* mat, void* vec) {
	PetscErrorCode     errorFlag;

	errorFlag = MatDiagonalSet( mat, vec, ADD_VALUES );
	CheckPETScError( errorFlag );
}

void Matrix_DiagonalSet_Insert( void* mat, void* vec) {
	PetscErrorCode     errorFlag;

	errorFlag = MatDiagonalSet( mat, vec, INSERT_VALUES );
	CheckPETScError( errorFlag );
}

void Matrix_Scale( Matrix* mat, double scaleFactor ) {
	PetscErrorCode     errorFlag;

	errorFlag = MatScale( &scaleFactor, (Mat)mat );
	CheckPETScError( errorFlag );
}


void Matrix_Destroy( Matrix* mat ) {
	PetscErrorCode     errorFlag;

	errorFlag = MatDestroy( (Mat)mat );
	CheckPETScError( errorFlag );
}

void MatrixMultiply( Matrix* A, Vector* x, Vector* b ) {
	PetscErrorCode     errorFlag;

	errorFlag = MatMult( (Mat)A, (Vec)x, (Vec)b );
	CheckPETScError( errorFlag );
}

void MatrixTransposeMultiply( Matrix* A, Vector* x, Vector* y ) {
	PetscErrorCode     errorFlag;
	
	errorFlag = MatMultTranspose( (Mat)A, (Vec)x, (Vec)y );
	CheckPETScError( errorFlag );
}
        
void MatrixMultiplyAdd( Matrix* A, Vector* x, Vector* y, Vector* z ) {
	PetscErrorCode     errorFlag;
	
	errorFlag = MatMultAdd( (Mat)A, (Vec)x, (Vec)y, (Vec)z );
	CheckPETScError( errorFlag );
}


void MatrixDuplicate_DataStructureOnly( Matrix* A, Matrix** B ) {
	PetscErrorCode     errorFlag;
	
	errorFlag = MatDuplicate( (Mat)A, MAT_DO_NOT_COPY_VALUES, (Mat*)B );
	CheckPETScError( errorFlag );
}

void MatrixDuplicate_IncludingValues( Matrix* A, Matrix** B ) {
	PetscErrorCode     errorFlag;
	
	errorFlag = MatDuplicate( (Mat)A, MAT_COPY_VALUES, (Mat*)B );
	CheckPETScError( errorFlag );
}

void Matrix_GetDiagonal( Matrix* matrix, Vector* vec ) {
	PetscErrorCode     errorFlag;
	
	errorFlag = MatGetDiagonal( (Mat)matrix, (Vec)vec );
	CheckPETScError( errorFlag );
}	
	
void Matrix_DiagonalScale( Matrix* matrix, Vector* leftScale, Vector* rightScale ) {
	PetscErrorCode     errorFlag;
	
	errorFlag = MatDiagonalScale( (Mat)matrix, (Vec)leftScale, (Vec)rightScale );
	CheckPETScError( errorFlag );
}	

void Matrix_AddScaledMatrix( Matrix* matrix, double scaleFactor, Matrix* matrixToAdd ) {
	PetscErrorCode errorFlag;

	errorFlag = MatAXPY( scaleFactor, (Mat) matrixToAdd, (Mat) matrix );
	CheckPETScError( errorFlag );
}

void Matrix_ScaleAndAddMatrix( Matrix* matrix, double scaleFactor, Matrix* matrixToAdd ) {
	PetscErrorCode errorFlag;

	errorFlag = MatAYPX( scaleFactor, (Mat) matrixToAdd, (Mat) matrix );
	CheckPETScError( errorFlag );
}

double Matrix_L2_Norm( Matrix* matrix ) {
	PetscErrorCode errorFlag;
	PetscReal    norm;
	
	errorFlag = MatNorm( (Mat) matrix, NORM_FROBENIUS, &norm );
	CheckPETScError(errorFlag);
	
	return (double)norm;
}


void Matrix_Transpose( Matrix* srcMat, Matrix** dstMat ) {
	PetscErrorCode	ef;
	
	ef = MatTranspose( (Mat)srcMat, (Mat*)dstMat );
	CheckPETScError( ef );
}


void Matrix_MatrixMult( Matrix* matA, Matrix* matB, Matrix** dst, double fillRatio ) {
	PetscErrorCode	ef;
	
	/* A guess at how this works... */
/*
	MatGetInfo( matA, MAT_GLOBAL_SUM, &mInfo );
	fill = (double)mInfo.nz_used;
	MatGetInfo( matB, MAT_GLOBAL_SUM, &mInfo );
	fill += (double)mInfo.nz_used;
	fill = expNonZeros / fill;
*/
	
	/* TODO: Figure out how this call works. */
	ef = MatMatMult( (Mat)matA, (Mat)matB, MAT_INITIAL_MATRIX, fillRatio, (Mat*)dst );
	CheckPETScError( ef );
}

void Matrix_MatrixMult_General( Matrix** cMat, Matrix* aMat, Matrix* bMat )
{
	int 		i,j;
	Vec 		tmpColVec, newColVec;
	PetscScalar*	col2Add;  
	int 		rows_A, cols_A;
	int 		rows_B, cols_B;
	int*		rowId;
	int 		start, end, nRows;
	Mat*		C = (Mat*)cMat;
	Mat		A = (Mat)aMat;
	Mat		B = (Mat)bMat;
	
	/*
	PetscViewerSetFormat (PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_DENSE);
	MatView (B, PETSC_VIEWER_STDOUT_WORLD);
	*/
	
	//Journal_DPrintf( snarkDebug, "In %s():\n", __func__ );
	MatGetSize( A, &rows_A, &cols_A );
	MatGetSize( B, &rows_B, &cols_B );
	
	if( cols_A != rows_B)  {
		Stream*	errorStr = Journal_Register( Error_Type, "PETSc Matrix" );
		Journal_Printf( errorStr, " Size of matrix A is not compatible with matrix B "
			"for multiplication.\n");
		exit( EXIT_FAILURE );
	}
	
	MatCreateMPIAIJ( MPI_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, rows_A, cols_B, PETSC_NULL,PETSC_NULL,
		PETSC_NULL, PETSC_NULL, C ); 
	MatSetFromOptions( *C );
	
	VecCreateMPI( MPI_COMM_WORLD, PETSC_DECIDE, rows_B, &tmpColVec );
	VecCreateMPI( MPI_COMM_WORLD, PETSC_DECIDE, rows_A, &newColVec );
	VecGetOwnershipRange( newColVec, &start, &end );
	//Journal_DPrintf( snarkDebug, "FirstRow = %d : LastRow = %d \n",start, end );
	nRows = end - start;
	
	PetscMalloc( nRows * sizeof(int), &rowId );
	
	
	for( i=0; i<cols_B; i++ )  { 
		MatGetColumnVector( B, tmpColVec, i );
		/* VecView(TempColVec, PETSC_VIEWER_STDOUT_WORLD); */
		MatMult( A, tmpColVec, newColVec );
		/* VecView(NewColVec, PETSC_VIEWER_STDOUT_WORLD);  */
		VecGetArray( newColVec, &col2Add );
		for(j=0; j<nRows; j++)  {
			rowId[j] = j + start;
			if( fabs(col2Add[j]) < 1.0e-15 ) {
				rowId[j] = -1;
			}
		}
		MatSetValues( *C, nRows, rowId, 1, &i, (PetscScalar*)col2Add, INSERT_VALUES );
		VecRestoreArray( newColVec, &col2Add );
		
		//Journal_DPrintf( snarkDebugThree, "done row %d \r",i);
	}  
	MatAssemblyBegin( *C, MAT_FINAL_ASSEMBLY );
	MatAssemblyEnd( *C, MAT_FINAL_ASSEMBLY );
	
	/*
	 PetscViewerSetFormat (PETSC_VIEWER_STDOUT_WORLD, PETSC_VIEWER_ASCII_DENSE);  
	 MatView (C, PETSC_VIEWER_STDOUT_WORLD);  
	*/
	PetscFree( rowId );
	VecDestroy( tmpColVec );
	VecDestroy( newColVec );
}


void Matrix_PtAP( Matrix* matA, Matrix* matP, Matrix** dst, double fillRatio ) {
	PetscErrorCode	ec;
	
	/* Assumes resulting matrix has a similar structure to A. */
	if( fillRatio == 0.0 ) {
		double	nRowsA, nRowsC;
		double	nzA, nzP, nzC;
		MatInfo	mInfo;
		
		MatGetInfo( (Mat)matA, MAT_GLOBAL_SUM, &mInfo );
		nRowsA = mInfo.rows_global;
		nzA = mInfo.nz_used;
		
		MatGetInfo( (Mat)matP, MAT_GLOBAL_SUM, &mInfo );
		nRowsC = mInfo.columns_global;
		nzP = mInfo.nz_used;
		
		nzC = (nRowsC / nRowsA) * nzA;
		fillRatio = nzC / (nzA + nzP);
	}
	
	if( *dst ) {
		ec = MatPtAP( (Mat)matA, (Mat)matP, MAT_REUSE_MATRIX, fillRatio, (Mat*)dst );	CheckPETScError( ec );
	}
	else {
		ec = MatPtAP( (Mat)matA, (Mat)matP, MAT_INITIAL_MATRIX, fillRatio, (Mat*)dst );	CheckPETScError( ec );
	}
}


void Matrix_GetLocalSize( Matrix* mat, unsigned* nRowsDst, unsigned* nColsDst ) {
	PetscErrorCode	ec;
	unsigned	nRows, nCols;

	ec = MatGetLocalSize( (Mat)mat, (PetscInt*)(&nRows), (PetscInt*)(&nCols) );	CheckPETScError( ec );
	if( nRowsDst ) {
		*nRowsDst = nRows;
	}
	if( nColsDst ) {
		*nColsDst = nCols;
	}
}


void Matrix_GetRow( Matrix* mat, unsigned rowInd, unsigned* nEntries, unsigned** cols, double** entries ) {
	PetscErrorCode	ec;

	assert( nEntries );
	assert( cols );
	assert( entries );

	ec = MatGetRow( (Mat)mat, rowInd, (PetscInt*)nEntries, (const PetscInt**)cols, (const PetscScalar**)entries );
	CheckPETScError( ec );
}


void Matrix_RestoreRow( Matrix* mat, unsigned rowInd, unsigned nEntries, unsigned** cols, double** entries ) {
	PetscErrorCode	ec;

	assert( cols );
	assert( entries );

	ec = MatRestoreRow( (Mat)mat, rowInd, (PetscInt*)(&nEntries), (const PetscInt**)cols, (const PetscScalar**)entries );
	CheckPETScError( ec );
}
