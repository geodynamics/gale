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
** $Id: Matrix.h 983 2007-11-14 04:04:09Z JulianGiordani $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __StgFEM_SLE_LinearAlgebra_Matrix_h__
#define __StgFEM_SLE_LinearAlgebra_Matrix_h__

	/** Textual name of this class */
	extern const Type Matrix_Type;

	/** Virtual function types */
	typedef void (Matrix_SetCommFunc)( void* matrix, MPI_Comm comm );
	typedef void (Matrix_SetGlobalSizeFunc)( void* matrix, unsigned nRows, unsigned nColumns );
	typedef void (Matrix_SetLocalSizeFunc)( void* matrix, unsigned nRows, unsigned nColumns );
	typedef void (Matrix_AddEntriesFunc)( void* matrix, unsigned nRows, unsigned* rowIndices, 
					      unsigned nColumns, unsigned* columnIndices, 
					      double* values );
	typedef void (Matrix_InsertEntriesFunc)( void* matrix, unsigned nRows, unsigned* rowIndices, 
						 unsigned nColumns, unsigned* columnIndices, 
						 double* values );
	typedef void (Matrix_DiagonalAddEntriesFunc)( void* matrix, void* vector );
	typedef void (Matrix_DiagonalInsertEntriesFunc)( void* matrix, void* vector );
	typedef void (Matrix_ZeroFunc)( void* matrix );
	typedef void (Matrix_DumpFunc)( void* matrix, const char* filename );
	typedef void (Matrix_LoadFunc)( void* matrix, char* filename );
	typedef void (Matrix_AssemblyBeginFunc)( void* matrix );
	typedef void (Matrix_AssemblyEndFunc)( void* matrix );

	typedef void (Matrix_ScaleFunc)( void* matrix, double factor );
	typedef void (Matrix_AddScaledFunc)( void* matrix, double factor, void* matrix0 );
	typedef void (Matrix_DiagonalScaleFunc)( void* matrix, void* leftVector, void* rightVector );
	typedef void (Matrix_MultiplyFunc)( void* matrix, void* vector, void* dstVector );
	typedef void (Matrix_TransposeMultiplyFunc)( void* matrix, void* vector, void* dstVector );
	typedef void (Matrix_MultiplyAddFunc)( void* matrix, void* vector0, void* vector1, void* dstVector );
	typedef void (Matrix_PAPtFunc)( void* matrix, void* P, void** dstMatrix );
	typedef void (Matrix_PtAPFunc)( void* matrix, void* P, void** dstMatrix );
	typedef void (Matrix_MatrixMultiplyFunc)( void* matrix, void* matrix0, void* dstMatrix );
	typedef double (Matrix_L2NormFunc)( void* matrix );
	typedef void (Matrix_TransposeFunc)( void* matrix, void* dstMatrix );

	typedef void (Matrix_GetGlobalSizeFunc)( void* matrix, unsigned* nRows, unsigned* nColumns );
	typedef void (Matrix_GetLocalSizeFunc)( void* matrix, unsigned* nRows, unsigned* nColumns );
	typedef void (Matrix_GetRowFunc)( void* matrix, unsigned rowIndex, 
					  unsigned* nEntries, unsigned** columnIndices, double** values );
	typedef void (Matrix_RestoreRowFunc)( void* matrix, unsigned rowIndex, 
					      unsigned* nEntries, unsigned** columnIndices, double** values );
	typedef void (Matrix_GetDiagonalFunc)( void* matrix, void* vector );
	typedef void (Matrix_DuplicateFunc)( void* matrix, void** dstMatrix );
	typedef void (Matrix_CopyEntriesFunc)( void* matrix, void* dstMatrix );

	/** Matrix class contents */
	#define __Matrix								\
		/* General info */							\
		__Stg_Component								\
											\
		/* Virtual info */							\
		Matrix_SetCommFunc*			setCommFunc;			\
		Matrix_SetGlobalSizeFunc*		setGlobalSizeFunc;		\
		Matrix_SetLocalSizeFunc*		setLocalSizeFunc;		\
		Matrix_AddEntriesFunc*			addEntriesFunc;			\
		Matrix_InsertEntriesFunc*		insertEntriesFunc;		\
		Matrix_DiagonalAddEntriesFunc*		diagonalAddEntriesFunc;		\
		Matrix_DiagonalInsertEntriesFunc*	diagonalInsertEntriesFunc;	\
		Matrix_ZeroFunc*			zeroFunc;			\
		Matrix_DumpFunc*			dumpFunc;			\
		Matrix_LoadFunc*			loadFunc;			\
		Matrix_AssemblyBeginFunc*		assemblyBeginFunc;		\
		Matrix_AssemblyEndFunc*			assemblyEndFunc;		\
											\
		Matrix_ScaleFunc*			scaleFunc;			\
		Matrix_AddScaledFunc*			addScaledFunc;			\
		Matrix_DiagonalScaleFunc*		diagonalScaleFunc;		\
		Matrix_MultiplyFunc*			multiplyFunc;			\
		Matrix_TransposeMultiplyFunc*		transposeMultiplyFunc;		\
		Matrix_MultiplyAddFunc*			multiplyAddFunc;		\
		Matrix_PAPtFunc*			paptFunc;			\
		Matrix_PtAPFunc*			ptapFunc;			\
		Matrix_MatrixMultiplyFunc*		matrixMultiplyFunc;		\
		Matrix_L2NormFunc*			l2NormFunc;			\
		Matrix_TransposeFunc*			transposeFunc;			\
											\
		Matrix_GetGlobalSizeFunc*		getGlobalSizeFunc;		\
		Matrix_GetLocalSizeFunc*		getLocalSizeFunc;		\
		Matrix_GetRowFunc*			getRowFunc;			\
		Matrix_RestoreRowFunc*			restoreRowFunc;			\
		Matrix_GetDiagonalFunc*			getDiagonalFunc;		\
		Matrix_DuplicateFunc*			duplicateFunc;			\
		Matrix_CopyEntriesFunc*			copyEntriesFunc;		\
											\
		/* Matrix info */					       		\
		MPI_Comm		comm;						\
		Bool			hasChanged;					\
		List*			solvers;

	struct Matrix { __Matrix };

	/*--------------------------------------------------------------------------------------------------------------------------
	** Constructors
	*/

	#define MATRIX_DEFARGS								\
		STG_COMPONENT_DEFARGS,							\
		Matrix_SetCommFunc*			setCommFunc, 			\
		Matrix_SetGlobalSizeFunc*		setGlobalSizeFunc,		\
		Matrix_SetLocalSizeFunc*		setLocalSizeFunc,		\
		Matrix_AddEntriesFunc*			addEntriesFunc,			\
		Matrix_InsertEntriesFunc*		insertEntriesFunc,		\
		Matrix_DiagonalAddEntriesFunc*		diagonalAddEntriesFunc,		\
		Matrix_DiagonalInsertEntriesFunc*	diagonalInsertEntriesFunc, 	\
		Matrix_ZeroFunc*			zeroFunc,			\
		Matrix_DumpFunc*			dumpFunc,			\
		Matrix_LoadFunc*			loadFunc,			\
		Matrix_AssemblyBeginFunc*		assemblyBeginFunc,		\
		Matrix_AssemblyEndFunc*			assemblyEndFunc,		\
											\
		Matrix_ScaleFunc*			scaleFunc,			\
		Matrix_AddScaledFunc*			addScaledFunc,			\
		Matrix_DiagonalScaleFunc*		diagonalScaleFunc,		\
		Matrix_MultiplyFunc*			multiplyFunc,			\
		Matrix_TransposeMultiplyFunc*		transposeMultiplyFunc,		\
		Matrix_MultiplyAddFunc*			multiplyAddFunc,		\
		Matrix_PAPtFunc*			paptFunc,			\
		Matrix_PtAPFunc*			ptapFunc,			\
		Matrix_MatrixMultiplyFunc*		matrixMultiplyFunc,		\
		Matrix_L2NormFunc*			l2NormFunc,			\
		Matrix_TransposeFunc*			transposeFunc,			\
											\
		Matrix_GetGlobalSizeFunc*		getGlobalSizeFunc,		\
		Matrix_GetLocalSizeFunc*		getLocalSizeFunc,		\
		Matrix_GetRowFunc*			getRowFunc,			\
		Matrix_RestoreRowFunc*			restoreRowFunc,			\
		Matrix_GetDiagonalFunc*			getDiagonalFunc,		\
		Matrix_DuplicateFunc*			duplicateFunc, 			\
		Matrix_CopyEntriesFunc*			copyEntriesFunc

	#define MATRIX_PASSARGS			\
		STG_COMPONENT_PASSARGS, 	\
		setCommFunc, 			\
		setGlobalSizeFunc,		\
		setLocalSizeFunc,		\
		addEntriesFunc,			\
		insertEntriesFunc,		\
		diagonalAddEntriesFunc,		\
		diagonalInsertEntriesFunc,	\
		zeroFunc,			\
		dumpFunc,			\
		loadFunc,			\
		assemblyBeginFunc,		\
		assemblyEndFunc,		\
						\
		scaleFunc,			\
		addScaledFunc,			\
		diagonalScaleFunc, 		\
		multiplyFunc,			\
		transposeMultiplyFunc,		\
		multiplyAddFunc,		\
		paptFunc,			\
		ptapFunc,			\
		matrixMultiplyFunc,		\
		l2NormFunc,			\
		transposeFunc,			\
						\
		getGlobalSizeFunc,		\
		getLocalSizeFunc,		\
		getRowFunc,			\
		restoreRowFunc,			\
		getDiagonalFunc,		\
		duplicateFunc, 			\
		copyEntriesFunc

	Matrix* _Matrix_New( MATRIX_DEFARGS );
	void _Matrix_Init( Matrix* self );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Virtual functions
	*/

	void _Matrix_Delete( void* matrix );
	void _Matrix_Print( void* matrix, Stream* stream );
	void _Matrix_Construct( void* matrix, Stg_ComponentFactory* cf, void* data );
	void _Matrix_Build( void* matrix, void* data );
	void _Matrix_Initialise( void* matrix, void* data );
	void _Matrix_Execute( void* matrix, void* data );
	void _Matrix_Destroy( void* matrix, void* data );

	void _Matrix_SetComm( void* matrix, MPI_Comm comm );
	void _Matrix_AssemblyBegin( void* matrix );

	#define Matrix_SetComm( self, comm )									\
		VirtualCall( self, setCommFunc, self, comm )

	#define Matrix_SetGlobalSize( self, nRows, nColumns )							\
		VirtualCall( self, setGlobalSizeFunc, self, nRows, nColumns )

	#define Matrix_SetLocalSize( self, nRows, nColumns )							\
		VirtualCall( self, setLocalSizeFunc, self, nRows, nColumns )

	#define Matrix_AddEntries( self, nRows, rowIndices, nColumns, columnIndices, values )			\
		VirtualCall( self, addEntriesFunc, self, nRows, rowIndices, nColumns, columnIndices, values )

	#define Matrix_InsertEntries( self, nRows, rowIndices, nColumns, columnIndices, values )		\
		VirtualCall( self, insertEntriesFunc, self, nRows, rowIndices, nColumns, columnIndices, values )

	#define Matrix_DiagonalAddEntries( self, vector )							\
		VirtualCall( self, diagonalAddEntriesFunc, self, vector )

	#define Matrix_DiagonalInsertEntries( self, vector )							\
		VirtualCall( self, diagonalInsertEntriesFunc, self, vector )

	#define Matrix_Zero( self )										\
		VirtualCall( self, zeroFunc, self )

	#define Matrix_Dump( self, filename )									\
		VirtualCall( self, dumpFunc, self, filename )

	#define Matrix_Load( self, filename )									\
		VirtualCall( self, loadFunc, self, filename )

	#define Matrix_AssemblyBegin( self )								\
		VirtualCall( self, assemblyBeginFunc, self )

	#define Matrix_AssemblyEnd( self )								\
		VirtualCall( self, assemblyEndFunc, self )

	#define Matrix_Scale( self, factor )								\
		VirtualCall( self, scaleFunc, self, factor )

	#define Matrix_AddScaled( self, factor, matrix0 )						\
		VirtualCall( self, addScaledFunc, self, factor, matrix0 )

	#define Matrix_DiagonalScale( self, leftVector, rightVector )					\
		VirtualCall( self, diagonalScaleFunc, self, leftVector, rightVector )

	#define Matrix_Multiply( self, vector, dstVector )						\
		VirtualCall( self, multiplyFunc, self, vector, dstVector )

	#define Matrix_TransposeMultiply( self, vector, dstVector )					\
		VirtualCall( self, transposeMultiplyFunc, self, vector, dstVector )

	#define Matrix_MultiplyAdd( self, vector0, vector1, dstVector )					\
		VirtualCall( self, multiplyAddFunc, self, vector0, vector1, dstVector )

	#define Matrix_PAPt( self, P, dstMatrix )							\
		VirtualCall( self, paptFunc, self, P, dstMatrix )

	#define Matrix_PtAP( self, P, dstMatrix )							\
		VirtualCall( self, ptapFunc, self, P, dstMatrix )

	#define Matrix_MatrixMultiply( self, matrix0, dstMatrix )					\
		VirtualCall( self, matrixMultiplyFunc, self, matrix0, dstMatrix )

	#define Matrix_L2Norm( self )									\
		VirtualCall( self, l2NormFunc, self )

	#define Matrix_Transpose( self, dstMatrix )							\
		VirtualCall( self, transposeFunc, self, dstMatrix )

	#define Matrix_GetGlobalSize( self, nRows, nColumns )						\
		VirtualCall( self, getGlobalSizeFunc, self, nRows, nColumns )

	#define Matrix_GetLocalSize( self, nRows, nColumns )						\
		VirtualCall( self, getLocalSizeFunc, self, nRows, nColumns )

	#define Matrix_GetRow( self, rowIndex, nEntries, columnIndices, values )			\
		VirtualCall( self, getRowFunc, self, rowIndex, nEntries, columnIndices, values )

	#define Matrix_RestoreRow( self, rowIndex, nEntries, columnIndices, values )			\
		VirtualCall( self, restoreRowFunc, self, rowIndex, nEntries, columnIndices, values )

	#define Matrix_GetDiagonal( self, vector )							\
		VirtualCall( self, getDiagonalFunc, self, vector )

	#define Matrix_Duplicate( self, dstMatrix )							\
		VirtualCall( self, duplicateFunc, self, dstMatrix )

	#define Matrix_CopyEntries( self, dstMatrix )							\
		VirtualCall( self, copyEntriesFunc, self, dstMatrix )

	/*--------------------------------------------------------------------------------------------------------------------------
	** Public functions
	*/

	void _Matrix_Dump( void* matrix, const char* filename );
	void Matrix_InvalidateSolvers( void* matrix );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Private Member functions
	*/

#endif /* __StgFEM_SLE_LinearAlgebra_Matrix_h__ */
