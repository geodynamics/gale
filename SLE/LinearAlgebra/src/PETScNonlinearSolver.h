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
** $Id: PETScNonlinearSolver.h 672 2006-12-14 00:58:36Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __StgFEM_SLE_LinearAlgebra_PETScNonlinearSolver_h__
#define __StgFEM_SLE_LinearAlgebra_PETScNonlinearSolver_h__

	/** Textual name of this class */
	extern const Type PETScNonlinearSolver_Type;

	/** Virtual function types */
	typedef PetscErrorCode (PETScNonlinearSolver_Func)( SNES snes, Vec x, Vec f, void* context );

	/** PETScNonlinearSolver class contents */
	#define __PETScNonlinearSolver		\
		/* General info */		\
		__NonlinearSolver		\
						\
		/* Virtual info */		\
						\
		/* PETScNonlinearSolver info */	\
		SNES		snes;		\
		KSP		ksp;		\

	struct PETScNonlinearSolver { __PETScNonlinearSolver };

	/*--------------------------------------------------------------------------------------------------------------------------
	** Constructors
	*/

	#define PETSCNONLINEARSOLVER_DEFARGS \
		NONLINEARSOLVER_DEFARGS

	#define PETSCNONLINEARSOLVER_PASSARGS \
		NONLINEARSOLVER_PASSARGS

	PETScNonlinearSolver* PETScNonlinearSolver_New( Name name ); /* wraps up the SNES create func. */
	PETScNonlinearSolver* _PETScNonlinearSolver_New( PETSCNONLINEARSOLVER_DEFARGS );
	void _PETScNonlinearSolver_Init( PETScNonlinearSolver* self );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Virtual functions
	*/

	void _PETScNonlinearSolver_Delete( void* nls );
	void _PETScNonlinearSolver_Print( void* nls, Stream* stream );
	void _PETScNonlinearSolver_Construct( void* nls, Stg_ComponentFactory* cf, void* data );
	void _PETScNonlinearSolver_Build( void* nls, void* data );
	void _PETScNonlinearSolver_Initialise( void* nls, void* data );
	void _PETScNonlinearSolver_Execute( void* nls, void* data );
	void _PETScNonlinearSolver_Destroy( void* nls, void* data );

	void PETScNonlinearSolver_SetComm( void* nls, MPI_Comm comm );

	MatrixSolver_Status PETScMatrixSolver_GetSolveStatus( void* nls );
	unsigned PETScNonlinearSolver_GetIterations( void* nls );
	unsigned PETScNonlinearSolver_GetMaxIterations( void* nls );
	double PETScNonlinearSolver_GetResidualNorm( void* nls );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Public functions
	*/

	void PETScNonlinearSolver_Create( void* nls );
	void PETScNonlinearSolver_Destroy( void* nls );
	void PETScNonlinearSolver_SetFunction( void* nls, void* f, void* func, void* context );
	void PETScNonlinearSolver_GetJacobian( void* nls, void* J, void* pc, /*no function to build J for now*/ void** context );
	void PETScNonlinearSolver_Solve( void* nls, void* b, void* x );
	void PETScNonlinearSolver_SetSolution( void* nls, void* x );
	void PETScNonlinearSolver_SetRhs( void* nls, void* rhs );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Private Member functions
	*/

#endif /* __StgFEM_SLE_LinearAlgebra_PETScNonlinearSolver_h__ */
