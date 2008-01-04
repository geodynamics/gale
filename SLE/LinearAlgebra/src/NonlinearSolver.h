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
** $Id: NonlinearSolver.h 822 2007-04-27 06:20:35Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __StgFEM_SLE_LinearAlgebra_NonlinearSolver_h__
#define __StgFEM_SLE_LinearAlgebra_NonlinearSolver_h__

	/** Textual name of this class */
	extern const Type NonlinearSolver_Type;

	typedef void* (NonlinearSolver_Func)( void* nls, Vector x, Vector f, void* context );

	/** Virtual function types */
	typedef void (NonlinearSolver_SetCommFunc)( void* nls, MPI_Comm comm );
	typedef void (NonlinearSolver_CreateFunc)( void* nls );
	typedef void (NonlinearSolver_DestroyFunc)( void* nls );	
	typedef void (NonlinearSolver_SetFunctionFunc)( void* nls, void* f, NonlinearSolver_Func* func, void* context );
	typedef void (NonlinearSolver_GetJacobianFunc)( void* nls, void* J, void* pc, /*no Jacobian func for now*/  void** context );
	typedef void (NonlinearSolver_SolveFunc)( void* nls, void* b, void* x );
	typedef void (NonlinearSolver_SetSolutionFunc)( void* nls, void* x );
	typedef void (NonlinearSolver_SetRhsFunc)( void* nls, void* rhs );

	typedef MatrixSolver_Status (NonlinearSolver_GetSolveStatusFunc)( void* nls );
	typedef unsigned (NonlinearSolver_GetIterationsFunc)( void* nls );
	typedef unsigned (NonlinearSolver_GetMaxIterationsFunc)( void* nls );
	typedef double (NonlinearSolver_GetResidualNormFunc)( void* nls );

	/** NonlinearSolver class contents */
	#define __NonlinearSolver								\
		/* General info */								\
		__Stg_Component									\
												\
		/* Virtual info */								\
		NonlinearSolver_SetCommFunc*			setCommFunc;			\
		NonlinearSolver_CreateFunc*			createFunc;			\
		NonlinearSolver_DestroyFunc*			destroyFunc;			\
		NonlinearSolver_SetFunctionFunc*		setFunctionFunc;		\
		NonlinearSolver_GetJacobianFunc*		getJacobianFunc;		\
		NonlinearSolver_SolveFunc*			solveFunc;			\
		NonlinearSolver_SetSolutionFunc*		setSolutionFunc;		\
		NonlinearSolver_SetRhsFunc*			setRhsFunc;			\
												\
		NonlinearSolver_GetSolveStatusFunc*		getSolveStatusFunc;		\
		NonlinearSolver_GetIterationsFunc*		getIterationsFunc;		\
		NonlinearSolver_GetMaxIterationsFunc*		getMaxIterationsFunc;		\
		NonlinearSolver_GetResidualNormFunc*		getResidualNormFunc;		\
												\
		/* NonlinearSolver info */							\
		MPI_Comm		comm;							\
		Matrix*			J;							\
		Matrix*			Jinv;							\
		Matrix*			pc;							\
		Vector*			residual;						\
		Bool			expiredResidual;					\
		double			toleranceNorm;						\
		double			residualNorm;						\
		unsigned		maxIterations;						\
												\
		Vector*			curRHS;							\
		Vector*			curSolution;

	struct NonlinearSolver { __NonlinearSolver };

	/*--------------------------------------------------------------------------------------------------------------------------
	** Constructors
	*/

	#define NONLINEARSOLVER_DEFARGS								\
		STG_COMPONENT_DEFARGS,								\
		NonlinearSolver_SetCommFunc*			setCommFunc,			\
		NonlinearSolver_CreateFunc*			createFunc,			\
		NonlinearSolver_DestroyFunc*			destroyFunc,			\
		NonlinearSolver_SetFunctionFunc*		setFunctionFunc,		\
		NonlinearSolver_GetJacobianFunc*		getJacobianFunc,		\
		NonlinearSolver_SolveFunc*			solveFunc,			\
		NonlinearSolver_SetSolutionFunc*		setSolutionFunc,		\
		NonlinearSolver_SetRhsFunc*			setRhsFunc,			\
												\
		NonlinearSolver_GetSolveStatusFunc*		getSolveStatusFunc,		\
		NonlinearSolver_GetIterationsFunc*		getIterationsFunc,		\
		NonlinearSolver_GetMaxIterationsFunc*		getMaxIterationsFunc,		\
		NonlinearSolver_GetResidualNormFunc*		getResidualNormFunc		\

	#define NONLINEARSOLVER_PASSARGS	\
		STG_COMPONENT_PASSARGS,		\
		setCommFunc, 			\
		createFunc,			\
		destroyFunc,			\
		setFunctionFunc,		\
		getJacobianFunc,		\
		solveFunc,			\
		setSolutionFunc,		\
		setRhsFunc,			\
						\
		getSolveStatusFunc,		\
		getIterationsFunc,		\
		getMaxIterationsFunc,		\
		getResidualNormFunc

	NonlinearSolver* _NonlinearSolver_New( NONLINEARSOLVER_DEFARGS );
	void _NonlinearSolver_Init( NonlinearSolver* self );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Virtual functions
	*/

	void _NonlinearSolver_Delete( void* matrixSolver );
	void _NonlinearSolver_Print( void* matrixSolver, Stream* stream );
	void _NonlinearSolver_Construct( void* matrixSolver, Stg_ComponentFactory* cf, void* data );
	void _NonlinearSolver_Build( void* matrixSolver, void* data );
	void _NonlinearSolver_Initialise( void* matrixSolver, void* data );
	void _NonlinearSolver_Execute( void* matrixSolver, void* data );
	void _NonlinearSolver_Destroy( void* matrixSolver, void* data );

	void _NonlinearSolver_SetComm( void* matrixSolver, MPI_Comm comm );

	#define NonlinearSolver_SetComm( self, comm )			\
		VirtualCall( self, setCommFunc, self, comm )

	#define NonlinearSolver_SetFunction( self, f, func, context )	\
		VirtualCall( self, setFunctionFunc, self, f, func, context )

	#define NonlinearSolver_GetJacobian( self, J, pc, context )	\
		VirtualCall( self, getJacobianFunc, self, j, pc, context )

	#define NonlinearSolver_Solve( self, b, x )			\
		VirtualCall( self, solveFunc, self, b, x )

	#define NonlinearSolver_SetSolution( self, x )			\
		VirtualCall( self, setSolutionFunc, self, x )

	#define NonlinearSolver_SetRhs( self, rhs )			\
		VirtualCall( self, setRhsFunc, self, rhs )
	/*
	#define NonlinearSolver_SetMaxIterations( self, nIterations )	\
		VirtualCall( self, setMaxIterationsFunc, self, nIterations )

	#define NonlinearSolver_SetRelativeTolerance( self, tolerance )	\
		VirtualCall( self, setRelativeToleranceFunc, self, tolerance )

	#define NonlinearSolver_SetAbsoluteTolerance( self, tolerance )	\
		VirtualCall( self, setAbsoluteToleranceFunc, self, tolerance )

	#define NonlinearSolver_SetUseInitialSolution( self, state )	\
		VirtualCall( self, setUseInitialSolutionFunc, self, state )*/

	#define NonlinearSolver_GetSolveStatus( self )			\
		VirtualCall( self, getSolveStatusFunc, self )

	#define NonlinearSolver_GetIterations( self )			\
		VirtualCall( self, getIterationsFunc, self )

	#define NonlinearSolver_GetMaxIterations( self )		\
		VirtualCall( self, getMaxIterationsFunc, self )

	#define NonlinearSolver_GetResidualNorm( self )			\
		VirtualCall( self, getResidualNormFunc, self )

	/*--------------------------------------------------------------------------------------------------------------------------
	** Public functions
	*/

	MPI_Comm NonlinearSolver_GetComm( void* nls );
	//Matrix* NonlinearSolver_GetJacobian( void* nls );
	//Vector* NonlinearSolver_GetResidual( void* nls );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Private Member functions
	*/

	void NonlinearSolver_CalcResidual( MatrixSolver* self );

#endif /* __StgFEM_SLE_LinearAlgebra_NonlinearSolver_h__ */
