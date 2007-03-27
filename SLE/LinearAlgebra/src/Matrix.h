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
*/
/** \file
**  Role:
**
** Assumptions:
**
** Comments:
**
** $Id: Matrix.h 656 2006-10-18 06:45:50Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __StgFEM_SLE_LinearAlgebra_Matrix_h__
#define __StgFEM_SLE_LinearAlgebra_Matrix_h__
	
	Matrix* Matrix_New( 
		MPI_Comm			comm,
		Index				rowSize,
		Index				colSize,
		Index				nonZeroCount );
	
	
	Matrix* Matrix_NewEmpty( MPI_Comm comm, unsigned rowLocalSize, unsigned colLocalSize );
	
	
	Matrix* Matrix_New2( 
		MPI_Comm			comm,
		Index				rowSize,
		Index				colSize, 
		Index*				diagnolNonZeroIndices,
		Index*				offDiagnolNonZeroIndices );

	Matrix* Matrix_NewFromFile( Name filename );
	
	void Matrix_Dump( Matrix* matrix, Name filename ) ;
	void Matrix_View( void* matrix, Stream* stream );
	
	void Matrix_Zero( void* matrix );
	void Matrix_AddTo( 
		void*				matrix,
		Index				rowCount,
		Index				rowIndices[],
		Index				colCount,
		Index				colIndices[],
		double*				values );
	void Matrix_Insert( 
		void*				matrix,
		Index				rowCount,
		Index				rowIndices[],
		Index				colCount,
		Index				colIndices[],
		double*				values );
	void Matrix_SetValue( void* matrix, unsigned rowInd, unsigned colInd, double val );
	void Matrix_AddValue( void* matrix, unsigned rowInd, unsigned colInd, double val );
	void Matrix_ZeroRow( void* matrix, unsigned rowInd );
	void Matrix_ZeroCol( void* matrix, unsigned colInd );
	void Matrix_AssemblyBegin( void* matrix );
	void Matrix_AssemblyEnd( void* matrix );

	void Matrix_DiagonalSet_Add( void* mat, void* vec);
	void Matrix_DiagonalSet_Insert( void* mat, void* vec);

	void Matrix_Scale( Matrix* mat, double scaleFactor );
	
	void MatrixMultiply( Matrix* A, Vector* x, Vector* b );
	void MatrixTransposeMultiply( Matrix* A, Vector* x, Vector* y );
	void MatrixMultiplyAdd( Matrix* A, Vector* x, Vector* y, Vector* z );
	void MatrixDuplicate_DataStructureOnly( Matrix* mat, Matrix** mat_copy );
	void MatrixDuplicate_IncludingValues( Matrix* A, Matrix** B );

	void Matrix_Destroy( Matrix* mat ); 

	void Matrix_AddScaledMatrix( Matrix* matrix, double scaleFactor, Matrix* matrixToAdd ) ;
	void Matrix_ScaleAndAddMatrix( Matrix* matrix, double scaleFactor, Matrix* matrixToAdd ) ;

	void Matrix_GetDiagonal( Matrix* matrix, Vector* vec );
	void Matrix_DiagonalScale( Matrix* matrix, Vector* leftScale, Vector* rightScale );
	
	double Matrix_L2_Norm( Matrix* matrix ) ;
	
	void Matrix_Transpose( Matrix* srcMat, Matrix** dstMat );
	
	void Matrix_MatrixMult( Matrix* matA, Matrix* matB, Matrix** dst, double fillRatio );
	void Matrix_MatrixMult_General( Matrix** cMat, Matrix* aMat, Matrix* bMat );

	void Matrix_PtAP( Matrix* matA, Matrix* matP, Matrix** dst, double fillRatio );
	void Matrix_GetLocalSize( Matrix* mat, unsigned* nRowsDst, unsigned* nColsDst );
	void Matrix_GetRow( Matrix* mat, unsigned rowInd, unsigned* nEntries, unsigned** cols, double** entries );
	void Matrix_RestoreRow( Matrix* mat, unsigned rowInd, unsigned nEntries, unsigned** cols, double** entries );

#endif /* __StgFEM_SLE_LinearAlgebra_Matrix_h__ */
