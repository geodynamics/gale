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
** $Id: PETScMatrixSolver.h 672 2006-12-14 00:58:36Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __StgFEM_SLE_LinearAlgebra_PETScMatrixSolver_h__
#define __StgFEM_SLE_LinearAlgebra_PETScMatrixSolver_h__

	/** Textual name of this class */
	extern const Type PETScMatrixSolver_Type;

	/** Virtual function types */

	/** PETScMatrixSolver class contents */
	#define __PETScMatrixSolver		\
		/* General info */		\
		__MatrixSolver			\
						\
		/* Virtual info */		\
						\
		/* PETScMatrixSolver info */	\
		KSP		ksp;		\
		Bool            optionsReady;

	struct PETScMatrixSolver { __PETScMatrixSolver };

	/*--------------------------------------------------------------------------------------------------------------------------
	** Constructors
	*/

	#define PETSCMATRIXSOLVER_DEFARGS \
		MATRIXSOLVER_DEFARGS

	#define PETSCMATRIXSOLVER_PASSARGS \
		MATRIXSOLVER_PASSARGS

	PETScMatrixSolver* PETScMatrixSolver_New( Name name );
	PETScMatrixSolver* _PETScMatrixSolver_New( PETSCMATRIXSOLVER_DEFARGS );
	void _PETScMatrixSolver_Init( PETScMatrixSolver* self );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Virtual functions
	*/

	void _PETScMatrixSolver_Delete( void* matrixSolver );
	void _PETScMatrixSolver_Print( void* matrixSolver, Stream* stream );
	void _PETScMatrixSolver_Construct( void* matrixSolver, Stg_ComponentFactory* cf, void* data );
	void _PETScMatrixSolver_Build( void* matrixSolver, void* data );
	void _PETScMatrixSolver_Initialise( void* matrixSolver, void* data );
	void _PETScMatrixSolver_Execute( void* matrixSolver, void* data );
	void _PETScMatrixSolver_Destroy( void* matrixSolver, void* data );

	void PETScMatrixSolver_SetComm( void* matrixSolver, MPI_Comm comm );
	void PETScMatrixSolver_SetMatrix( void* matrixSolver, void* matrix );
	void PETScMatrixSolver_SetMaxIterations( void* matrixSolver, unsigned nIterations );
	void PETScMatrixSolver_SetRelativeTolerance( void* matrixSolver, double tolerance );
	void PETScMatrixSolver_SetAbsoluteTolerance( void* matrixSolver, double tolerance );
	void PETScMatrixSolver_SetUseInitialSolution( void* matrixSolver, Bool state );

	void PETScMatrixSolver_Solve( void* matrixSolver, void* rhs, void* solution );
	void PETScMatrixSolver_Setup( void* matrixSolver, void* rhs, void* solution );

	MatrixSolver_Status PETScMatrixSolver_GetSolveStatus( void* matrixSolver );
	unsigned PETScMatrixSolver_GetIterations( void* matrixSolver );
	unsigned PETScMatrixSolver_GetMaxIterations( void* matrixSolver );
	double PETScMatrixSolver_GetResidualNorm( void* matrixSolver );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Public functions
	*/

	void PETScMatrixSolver_SetKSPType( void* matrixSolver, PETScMatrixSolver_KSPType type );
	void PETScMatrixSolver_SetPCType( void* matrixSolver, PETScMatrixSolver_PCType type );
	void PETScMatrixSolver_GetSubBlocks( void* matrixSolver, unsigned* nBlocks, KSP** ksps );
	void PETScMatrixSolver_EnableShifting( void* matrixSolver, Bool state );
	void PETScMatrixSolver_SetNormType( void* matrixSolver, PETScMatrixSolver_NormType normType );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Private Member functions
	*/

#endif /* __StgFEM_SLE_LinearAlgebra_PETScMatrixSolver_h__ */
