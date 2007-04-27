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
*/
/** \file
**  Role:
**
** Assumptions:
**
** Invariants:
**
** Comments:
**
** $Id: PETScMatrix.h 3584 2006-05-16 11:11:07Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __StgFEM_SLE_LinearAlgebra_PETScMatrix_h__
#define __StgFEM_SLE_LinearAlgebra_PETScMatrix_h__

	/** Textual name of this class */
	extern const Type PETScMatrix_Type;

	/** Virtual function types */
	typedef void (PETScMatrix_SetNonZeroStructure)( void* matrix, unsigned nNonZeros, 
							unsigned* diagonalNonZeroIndices, unsigned* offDiagonalNonZeroIndices);

	/** PETScMatrix class contents */
	#define __PETScMatrix							\
		/* General info */						\
		__Matrix							\
										\
		/* Virtual info */						\
		PETScMatrix_SetNonZeroStructure*	setNonZeroStructure;	\
										\
		/* PETScMatrix info */						\
		Mat		petscMat;

	struct PETScMatrix { __PETScMatrix };

	/*--------------------------------------------------------------------------------------------------------------------------
	** Constructors
	*/

	#define PETSCMATRIX_DEFARGS						\
		MATRIX_DEFARGS, 						\
		PETScMatrix_SetNonZeroStructure*	setNonZeroStructure

	#define PETSCMATRIX_PASSARGS	\
		MATRIX_PASSARGS, 	\
		setNonZeroStructure

	PETScMatrix* PETScMatrix_New( Name name );
	PETScMatrix* _PETScMatrix_New( PETSCMATRIX_DEFARGS );
	void _PETScMatrix_Init( PETScMatrix* self );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Virtual functions
	*/

	void _PETScMatrix_Delete( void* matrix );
	void _PETScMatrix_Print( void* matrix, Stream* stream );
	void _PETScMatrix_Construct( void* matrix, Stg_ComponentFactory* cf, void* data );
	void _PETScMatrix_Build( void* matrix, void* data );
	void _PETScMatrix_Initialise( void* matrix, void* data );
	void _PETScMatrix_Execute( void* matrix, void* data );
	void _PETScMatrix_Destroy( void* matrix, void* data );

	void PETScMatrix_SetComm( void* matrix, MPI_Comm comm );
	void PETScMatrix_SetGlobalSize( void* matrix, unsigned nRows, unsigned nColumns );
	void PETScMatrix_SetLocalSize( void* matrix, unsigned nRows, unsigned nColumns );
	void PETScMatrix_AddEntries( void* matrix, unsigned nRows, unsigned* rowIndices, 
				     unsigned nColumns, unsigned* columnIndices, 
				     double* values );
	void PETScMatrix_InsertEntries( void* matrix, unsigned nRows, unsigned* rowIndices, 
					unsigned nColumns, unsigned* columnIndices, 
					double* values );
	void PETScMatrix_DiagonalAddEntries( void* matrix, void* vector );
	void PETScMatrix_DiagonalInsertEntries( void* matrix, void* vector );
	void PETScMatrix_Zero( void* matrix );
	void PETScMatrix_Load( void* matrix, char* filename );
	void PETScMatrix_AssemblyBegin( void* matrix );
	void PETScMatrix_AssemblyEnd( void* matrix );

	void PETScMatrix_Scale( void* matrix, double factor );
	void PETScMatrix_AddScaled( void* matrix, double factor, void* matrix0 );
	void PETScMatrix_DiagonalScale( void* matrix, void* leftVector, void* rightVector );
	void PETScMatrix_Multiply( void* matrix, void* vector, void* dstVector );
	void PETScMatrix_TransposeMultiply( void* matrix, void* vector, void* dstVector );
	void PETScMatrix_MultiplyAdd( void* matrix, void* vector0, void* vector1, void* dstVector );
	void PETScMatrix_PAPt( void* matrix, void* P, void** dstMatrix );
	void PETScMatrix_PtAP( void* matrix, void* P, void** dstMatrix );
	void PETScMatrix_MatrixMultiply( void* matrix, void* matrix0, void* dstMatrix );
	double PETScMatrix_L2Norm( void* matrix );
	void PETScMatrix_Transpose( void* matrix, void* dstMatrix );

	void PETScMatrix_GetGlobalSize( void* matrix, unsigned* nRows, unsigned* nColumns );
	void PETScMatrix_GetLocalSize( void* matrix, unsigned* nRows, unsigned* nColumns );
	void PETScMatrix_GetRow( void* matrix, unsigned rowIndex, 
				 unsigned* nEntries, unsigned** columnIndices, double** values );
	void PETScMatrix_RestoreRow( void* matrix, unsigned rowIndex, 
				     unsigned* nEntries, unsigned** columnIndices, double** values );
	void PETScMatrix_GetDiagonal( void* matrix, void* vector );
	void PETScMatrix_Duplicate( void* matrix, void** dstMatrix );
	void PETScMatrix_CopyEntries( void* matrix, void* dstMatrix );

	void _PETScMatrix_SetNonZeroStructure( void* matrix, unsigned nNonZeros, 
					       unsigned* diagonalNonZeroIndices, unsigned* offDiagonalNonZeroIndices);

	#define PETScMatrix_SetNonZeroStructure( self, nNonZeros, diagonalNonZeroIndices, offDiagonalNonZeroIndices)		\
		VirtualCall( self, setNonZeroStructure, self, nNonZeros, diagonalNonZeroIndices, offDiagonalNonZeroIndices)

	/*--------------------------------------------------------------------------------------------------------------------------
	** Public functions
	*/

	void PETScMatrix_Draw( void* matrix );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Private Member functions
	*/

#endif /* __StgFEM_SLE_LinearAlgebra_PETScMatrix_h__ */
