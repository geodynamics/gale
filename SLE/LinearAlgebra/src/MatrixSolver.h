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
** $Id: MatrixSolver.h 822 2007-04-27 06:20:35Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __StgFEM_SLE_LinearAlgebra_MatrixSolver_h__
#define __StgFEM_SLE_LinearAlgebra_MatrixSolver_h__

	/** Textual name of this class */
	extern const Type MatrixSolver_Type;

	/** Virtual function types */
	typedef void (MatrixSolver_SetCommFunc)( void* matrixSolver, MPI_Comm comm );
	typedef void (MatrixSolver_SetMatrixFunc)( void* matrixSolver, void* matrix );
	typedef void (MatrixSolver_SetMaxIterationsFunc)( void* matrixSolver, unsigned nIterations );
	typedef void (MatrixSolver_SetRelativeToleranceFunc)( void* matrixSolver, double tolerance );
	typedef void (MatrixSolver_SetAbsoluteToleranceFunc)( void* matrixSolver, double tolerance );
	typedef void (MatrixSolver_SetUseInitialSolutionFunc)( void* matrixSolver, Bool state );

	typedef void (MatrixSolver_SolveFunc)( void* matrixSolver, void* rhs, void* solution );
	typedef void (MatrixSolver_SetupFunc)( void* matrixSolver, void* rhs, void* solution );

	typedef MatrixSolver_Status (MatrixSolver_GetSolveStatusFunc)( void* matrixSolver );
	typedef unsigned (MatrixSolver_GetIterationsFunc)( void* matrixSolver );
	typedef unsigned (MatrixSolver_GetMaxIterationsFunc)( void* matrixSolver );
	typedef double (MatrixSolver_GetResidualNormFunc)( void* matrixSolver );

	/** MatrixSolver class contents */
	#define __MatrixSolver									\
		/* General info */								\
		__Stg_Component									\
												\
		/* Virtual info */								\
		MatrixSolver_SetCommFunc*			setCommFunc;			\
		MatrixSolver_SetMatrixFunc*			setMatrixFunc;			\
		MatrixSolver_SetMaxIterationsFunc*		setMaxIterationsFunc;		\
		MatrixSolver_SetRelativeToleranceFunc*		setRelativeToleranceFunc;	\
		MatrixSolver_SetAbsoluteToleranceFunc*		setAbsoluteToleranceFunc;	\
		MatrixSolver_SetUseInitialSolutionFunc*		setUseInitialSolutionFunc;	\
												\
		MatrixSolver_SolveFunc*				solveFunc;			\
		MatrixSolver_SetupFunc*				setupFunc;			\
												\
		MatrixSolver_GetSolveStatusFunc*		getSolveStatusFunc;		\
		MatrixSolver_GetIterationsFunc*			getIterationsFunc;		\
		MatrixSolver_GetMaxIterationsFunc*		getMaxIterationsFunc;		\
		MatrixSolver_GetResidualNormFunc*		getResidualNormFunc;		\
												\
		/* MatrixSolver info */								\
		MPI_Comm		comm;							\
		Matrix*			matrix;							\
		Matrix*			inversion;						\
		Vector*			residual;						\
		Bool			expiredResidual;					\
		Bool			matrixChanged;						\
												\
		Vector*			curRHS;							\
		Vector*			curSolution;

	struct MatrixSolver { __MatrixSolver };

	/*--------------------------------------------------------------------------------------------------------------------------
	** Constructors
	*/

	#define MATRIXSOLVER_DEFARGS								\
		STG_COMPONENT_DEFARGS,								\
		MatrixSolver_SetCommFunc*			setCommFunc,			\
		MatrixSolver_SetMatrixFunc*			setMatrixFunc,			\
		MatrixSolver_SetMaxIterationsFunc*		setMaxIterationsFunc,		\
		MatrixSolver_SetRelativeToleranceFunc*		setRelativeToleranceFunc,	\
		MatrixSolver_SetAbsoluteToleranceFunc*		setAbsoluteToleranceFunc,	\
		MatrixSolver_SetUseInitialSolutionFunc*		setUseInitialSolutionFunc,	\
												\
		MatrixSolver_SolveFunc*				solveFunc,			\
		MatrixSolver_SetupFunc*				setupFunc,			\
												\
		MatrixSolver_GetSolveStatusFunc*		getSolveStatusFunc,		\
		MatrixSolver_GetIterationsFunc*			getIterationsFunc,		\
		MatrixSolver_GetMaxIterationsFunc*		getMaxIterationsFunc,		\
		MatrixSolver_GetResidualNormFunc*		getResidualNormFunc

	#define MATRIXSOLVER_PASSARGS		\
		STG_COMPONENT_PASSARGS,		\
		setCommFunc, 			\
		setMatrixFunc,			\
		setMaxIterationsFunc,		\
		setRelativeToleranceFunc,	\
		setAbsoluteToleranceFunc,	\
		setUseInitialSolutionFunc,	\
						\
		solveFunc,			\
		setupFunc,			\
						\
		getSolveStatusFunc,		\
		getIterationsFunc,		\
		getMaxIterationsFunc,		\
		getResidualNormFunc

	MatrixSolver* _MatrixSolver_New( MATRIXSOLVER_DEFARGS );
	void _MatrixSolver_Init( MatrixSolver* self );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Virtual functions
	*/

	void _MatrixSolver_Delete( void* matrixSolver );
	void _MatrixSolver_Print( void* matrixSolver, Stream* stream );
	void _MatrixSolver_Construct( void* matrixSolver, Stg_ComponentFactory* cf, void* data );
	void _MatrixSolver_Build( void* matrixSolver, void* data );
	void _MatrixSolver_Initialise( void* matrixSolver, void* data );
	void _MatrixSolver_Execute( void* matrixSolver, void* data );
	void _MatrixSolver_Destroy( void* matrixSolver, void* data );

	void _MatrixSolver_SetComm( void* matrixSolver, MPI_Comm comm );
	void _MatrixSolver_SetMatrix( void* matrixSolver, void* matrix );
	void _MatrixSolver_Setup( void* matrixSolver, void* rhs, void* solution );

	#define MatrixSolver_SetComm( self, comm )			\
		VirtualCall( self, setCommFunc, self, comm )

	#define MatrixSolver_SetMatrix( self, matrix )			\
		VirtualCall( self, setMatrixFunc, self, matrix )

	#define MatrixSolver_SetMaxIterations( self, nIterations )	\
		VirtualCall( self, setMaxIterationsFunc, self, nIterations )

	#define MatrixSolver_SetRelativeTolerance( self, tolerance )	\
		VirtualCall( self, setRelativeToleranceFunc, self, tolerance )

	#define MatrixSolver_SetAbsoluteTolerance( self, tolerance )	\
		VirtualCall( self, setAbsoluteToleranceFunc, self, tolerance )

	#define MatrixSolver_SetUseInitialSolution( self, state )	\
		VirtualCall( self, setUseInitialSolutionFunc, self, state )

	#define MatrixSolver_Solve( self, rhs, solution )			\
		VirtualCall( self, solveFunc, self, rhs, solution )

	#define MatrixSolver_Setup( self, rhs, solution )		\
		VirtualCall( self, setupFunc, self, rhs, solution )

	#define MatrixSolver_GetSolveStatus( self )			\
		VirtualCall( self, getSolveStatusFunc, self )

	#define MatrixSolver_GetIterations( self )			\
		VirtualCall( self, getIterationsFunc, self )

	#define MatrixSolver_GetMaxIterations( self )			\
		VirtualCall( self, getMaxIterationsFunc, self )

	#define MatrixSolver_GetResidualNorm( self )			\
		VirtualCall( self, getResidualNormFunc, self )

	/*--------------------------------------------------------------------------------------------------------------------------
	** Public functions
	*/

	MPI_Comm MatrixSolver_GetComm( void* matrixSolver );
	Matrix* MatrixSolver_GetMatrix( void* matrixSolver );
	Vector* MatrixSolver_GetResidual( void* matrixSolver );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Private Member functions
	*/

	void MatrixSolver_CalcResidual( MatrixSolver* self );

#endif /* __StgFEM_SLE_LinearAlgebra_MatrixSolver_h__ */
