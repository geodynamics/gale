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
** $Id: Matrix.c 3584 2006-05-16 11:11:07Z PatrickSunter $
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
const Type Matrix_Type = "Matrix";


/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/

Matrix* _Matrix_New( MATRIX_DEFARGS ) {
	Matrix*	self;

	/* Allocate memory */
	assert( sizeOfSelf >= sizeof(Matrix) );
	self = (Matrix*)_Stg_Component_New( STG_COMPONENT_PASSARGS );

	/* Virtual info */
	self->setCommFunc = setCommFunc;
	self->setGlobalSizeFunc = setGlobalSizeFunc;
	self->setLocalSizeFunc = setLocalSizeFunc;
	self->addEntriesFunc = addEntriesFunc;
	self->insertEntriesFunc = insertEntriesFunc;
	self->diagonalAddEntriesFunc = diagonalAddEntriesFunc;
	self->diagonalInsertEntriesFunc = diagonalInsertEntriesFunc;
	self->zeroFunc = zeroFunc;
	self->loadFunc = loadFunc;
	self->assemblyBeginFunc = assemblyBeginFunc;
	self->assemblyEndFunc = assemblyEndFunc;

	self->scaleFunc = scaleFunc;
	self->addScaledFunc = addScaledFunc;
	self->diagonalScaleFunc = diagonalScaleFunc;
	self->multiplyFunc = multiplyFunc;
	self->transposeMultiplyFunc = transposeMultiplyFunc;
	self->multiplyAddFunc = multiplyAddFunc;
	self->paptFunc = paptFunc;
	self->ptapFunc = ptapFunc;
	self->matrixMultiplyFunc = matrixMultiplyFunc;
	self->l2NormFunc = l2NormFunc;
	self->transposeFunc = transposeFunc;

	self->getGlobalSizeFunc = getGlobalSizeFunc;
	self->getLocalSizeFunc = getLocalSizeFunc;
	self->getRowFunc = getRowFunc;
	self->restoreRowFunc = restoreRowFunc;
	self->getDiagonalFunc = getDiagonalFunc;
	self->duplicateFunc = duplicateFunc;
	self->copyEntriesFunc = copyEntriesFunc;

	/* Matrix info */
	_Matrix_Init( self );

	return self;
}

void _Matrix_Init( Matrix* self ) {
	assert( self && Stg_CheckType( self, Matrix ) );

	self->comm = MPI_COMM_WORLD;
	self->hasChanged = True;
	self->solvers = List_New();
	List_SetItemSize( self->solvers, sizeof(MatrixSolver*) );
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _Matrix_Delete( void* matrix ) {
	Matrix*		self = (Matrix*)matrix;
	MatrixSolver*	solver;
	unsigned	s_i;

	assert( self && Stg_CheckType( self, Matrix ) );

	/* Remove from solvers. */
	for( s_i = 0; s_i < List_GetSize( self->solvers ); s_i++ ) {
		solver = *List_Get( self->solvers, s_i, MatrixSolver* );
		solver->matrix = NULL;
	}
	FreeObject( self->solvers );

	/* Delete the parent. */
	_Stg_Component_Delete( self );
}

void _Matrix_Print( void* matrix, Stream* stream ) {
	Matrix*	self = (Matrix*)matrix;
	
	/* Set the Journal for printing informations */
	Stream* matrixStream;
	matrixStream = Journal_Register( InfoStream_Type, "MatrixStream" );

	assert( self && Stg_CheckType( self, Matrix ) );

	/* Print parent */
	Journal_Printf( stream, "Matrix (ptr): (%p)\n", self );
	_Stg_Component_Print( self, stream );
}

void _Matrix_Construct( void* matrix, Stg_ComponentFactory* cf, void* data ) {
	Matrix*	self = (Matrix*)matrix;

	assert( self && Stg_CheckType( self, Matrix ) );
	assert( cf );
}

void _Matrix_Build( void* matrix, void* data ) {
}

void _Matrix_Initialise( void* matrix, void* data ) {
}

void _Matrix_Execute( void* matrix, void* data ) {
}

void _Matrix_Destroy( void* matrix, void* data ) {
}

void _Matrix_SetComm( void* matrix, MPI_Comm comm ) {
	Matrix*	self = (Matrix*)matrix;

	assert( self && Stg_CheckType( self, Matrix ) );

	self->comm = comm;
}

void _Matrix_AssemblyBegin( void* matrix ) {
	Matrix*	self = (Matrix*)matrix;

	assert( self && Stg_CheckType( self, Matrix ) );

	if( self->hasChanged ) {
		Matrix_InvalidateSolvers( self );
		self->hasChanged = False;
	}
}


/*--------------------------------------------------------------------------------------------------------------------------
** Public Functions
*/

void Matrix_Dump( void* matrix, const char* filename ) {
	Matrix*		self = (Matrix*)matrix;
	FILE*		fp;
	unsigned	nRows, nEntries, *entries;
	double*		values;
	unsigned	row_i, entry_i, b_i;

	assert( self && Stg_CheckType( self, Matrix ) );

	insist( fp = fopen( filename, "w" ), != NULL );

	Matrix_GetLocalSize( self, &nRows, NULL );
	for( row_i = 0; row_i < nRows; row_i++ ) {
		Matrix_GetRow( self, row_i, &nEntries, &entries, &values );
		for( entry_i = 0; entry_i < nEntries; entry_i++ ) {
			fprintf( fp, "(%d, %d): ", row_i, entries[entry_i] );
			for( b_i = 0; b_i < sizeof(double) / sizeof(unsigned); b_i++ )
				fprintf( fp, "%x", ((unsigned*)(values + entry_i))[b_i] );
			fprintf( fp, "\n" );
		}
		Matrix_RestoreRow( self, row_i, &nEntries, &entries, &values );
	}

	fclose( fp );
}

void Matrix_InvalidateSolvers( void* matrix ) {
	Matrix*		self = (Matrix*)matrix;
	MatrixSolver*	solver;
	unsigned	s_i;

	assert( self && Stg_CheckType( self, Matrix ) );

	for( s_i = 0; s_i < List_GetSize( self->solvers ); s_i++ ) {
		solver = *List_Get( self->solvers, s_i, MatrixSolver* );
		solver->matrixChanged = True;
	}
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Private Functions
*/
