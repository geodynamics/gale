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
** $Id: MultigridSolver.h 3584 2006-05-16 11:11:07Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __StgFEM_SLE_SystemSetup_MultigridSolver_h__
#define __StgFEM_SLE_SystemSetup_MultigridSolver_h__

	/** Textual name of this class */
	extern const Type MultigridSolver_Type;

	/** MultigridSolver class contents */
	typedef struct {
		MGSolver_PETScData*	downSolver;
		int		nDownIts;
		MGSolver_PETScData*	upSolver;
		int		nUpIts;
		unsigned		nCycles;

		Mat			A;
		Mat			R;
		Mat			P;

		Vec			workRHS;
		Vec			workSol;
	} MultigridSolver_Level;

	#define __MultigridSolver				\
		/* General info */				\
		__Stg_Component					\
								\
		/* Virtual info */				\
		MGSolver_SetCommFunc*			setCommFunc;			\
		MGSolver_SetMatrixFunc*			setMatrixFunc;			\
		MGSolver_SetMaxIterationsFunc*		setMaxIterationsFunc;		\
		MGSolver_SetRelativeToleranceFunc*	setRelativeToleranceFunc;	\
		MGSolver_SetAbsoluteToleranceFunc*	setAbsoluteToleranceFunc;	\
		MGSolver_SetUseInitialSolutionFunc*	setUseInitialSolutionFunc;	\
											\
		MGSolver_SolveFunc*			solveFunc;			\
		MGSolver_SetupFunc*			setupFunc;			\
											\
		MGSolver_GetSolveStatusFunc*		getSolveStatusFunc;		\
		MGSolver_GetIterationsFunc*		getIterationsFunc;		\
		MGSolver_GetMaxIterationsFunc*		getMaxIterationsFunc;		\
		MGSolver_GetResidualNormFunc*		getResidualNormFunc;		\
								\
		/* MultigridSolver info */			\
		MGSolver_PETScData*		mgData;		\
		Stream*				stream;		\
								\
		unsigned			nLevels;	\
		MultigridSolver_Level*		levels;		\
		MGOpGenerator*			opGen;		\
		Bool				solversChanged;	\
		Bool				opsChanged;	\
								\
		unsigned			curIt;		\
		unsigned			maxIts;		\
		double				relTol;		\
		double				rnorm;		\
		Bool				useInitial;	\
		MGSolver_PETScData*		outerSolver;

	struct MultigridSolver { __MultigridSolver };

	/*--------------------------------------------------------------------------------------------------------------------------
	** Constructors
	*/



	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define MULTIGRIDSOLVER_DEFARGS \
                STG_COMPONENT_DEFARGS, \
                MGSolver_SetCommFunc*                              setCommFunc, \
                MGSolver_SetMatrixFunc*                          setMatrixFunc, \
                MGSolver_SetMaxIterationsFunc*            setMaxIterationsFunc, \
                MGSolver_SetRelativeToleranceFunc*    setRelativeToleranceFunc, \
                MGSolver_SetAbsoluteToleranceFunc*    setAbsoluteToleranceFunc, \
                MGSolver_SetUseInitialSolutionFunc*  setUseInitialSolutionFunc, \
                MGSolver_SolveFunc*                                  solveFunc, \
                MGSolver_SetupFunc*                                  setupFunc, \
                MGSolver_GetSolveStatusFunc*                getSolveStatusFunc, \
                MGSolver_GetIterationsFunc*                  getIterationsFunc, \
                MGSolver_GetMaxIterationsFunc*            getMaxIterationsFunc, \
                MGSolver_GetResidualNormFunc*              getResidualNormFunc

	#define MULTIGRIDSOLVER_PASSARGS \
                STG_COMPONENT_PASSARGS, \
	        setCommFunc,               \
	        setMatrixFunc,             \
	        setMaxIterationsFunc,      \
	        setRelativeToleranceFunc,  \
	        setAbsoluteToleranceFunc,  \
	        setUseInitialSolutionFunc, \
	        solveFunc,                 \
	        setupFunc,                 \
	        getSolveStatusFunc,        \
	        getIterationsFunc,         \
	        getMaxIterationsFunc,      \
	        getResidualNormFunc      

	MultigridSolver* MultigridSolver_New( Name name );
	MultigridSolver* _MultigridSolver_New(  MULTIGRIDSOLVER_DEFARGS  );
	void _MultigridSolver_Init( MultigridSolver* self );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Virtual functions
	*/

	void _MultigridSolver_Delete( void* matrixSolver );
	void _MultigridSolver_Print( void* matrixSolver, Stream* stream );
	void _MultigridSolver_AssignFromXML( void* matrixSolver, Stg_ComponentFactory* cf, void* data );
	void _MultigridSolver_Build( void* matrixSolver, void* data );
	void _MultigridSolver_Initialise( void* matrixSolver, void* data );
	void _MultigridSolver_Execute( void* matrixSolver, void* data );
	void _MultigridSolver_Destroy( void* matrixSolver, void* data );

	void MultigridSolver_SetComm( void* matrixSolver, MPI_Comm comm );
	void MultigridSolver_SetMatrix( void* matrixSolver, void* _matrix );
	void MultigridSolver_SetMaxIterations( void* matrixSolver, unsigned nIterations );
	void MultigridSolver_SetRelativeTolerance( void* matrixSolver, double tolerance );
	void MultigridSolver_SetAbsoluteTolerance( void* matrixSolver, double tolerance );
	void MultigridSolver_SetUseInitialSolution( void* matrixSolver, Bool state );

	void MultigridSolver_Solve( void* matrixSolver, void* rhs, void* solution ); //need this
	void MultigridSolver_Setup( void* matrixSolver, void* rhs, void* solution ); //need this

	MGSolver_Status MultigridSolver_GetSolveStatus( void* matrixSolver ); //need this
	unsigned MultigridSolver_GetIterations( void* matrixSolver ); //need this
	unsigned MultigridSolver_GetMaxIterations( void* matrixSolver ); //need this
	double MultigridSolver_GetResidualNorm( void* matrixSolver ); //need this - i think??

	/*--------------------------------------------------------------------------------------------------------------------------
	** Public functions
	*/

	void MultigridSolver_SetLevels( void* matrixSolver, unsigned nLevels );
	void MultigridSolver_SetRestriction( void* matrixSolver, unsigned level, void* R );
	void MultigridSolver_SetProlongation( void* matrixSolver, unsigned level, void* P );
	void MultigridSolver_SetLevelDownSolver( void* matrixSolver, unsigned level, void* solver );
	void MultigridSolver_SetLevelDownIterations( void* matrixSolver, unsigned level, unsigned nIts );
	void MultigridSolver_SetLevelUpSolver( void* matrixSolver, unsigned level, void* solver );
	void MultigridSolver_SetLevelUpIterations( void* matrixSolver, unsigned level, unsigned nIts );
	void MultigridSolver_SetLevelCycles( void* matrixSolver, unsigned level, unsigned nCycles );
	void MultigridSolver_SetAllDownSolver( void* matrixSolver, void* solver );
	void MultigridSolver_SetAllDownIterations( void* matrixSolver, unsigned level, unsigned nIts );
	void MultigridSolver_SetAllUpSolver( void* matrixSolver, void* solver );
	void MultigridSolver_SetAllUpIterations( void* matrixSolver, unsigned level, unsigned nIts );
	void MultigridSolver_SetAllSolver( void* matrixSolver, void* solver );
	void MultigridSolver_SetCoarseSolver( void* matrixSolver, void* solver );

	unsigned MultigridSolver_GetNumLevels( void* matrixSolver );
	Mat MultigridSolver_GetRestriction( void* matrixSolver, unsigned level );
	Mat MultigridSolver_GetProlongation( void* matrixSolver, unsigned level );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Private Member functions
	*/

	//void MultigridSolver_RestrictMatrix( MultigridSolver* self, MultigridSolver_Level* level, Matrix** dstMatrix );
	//void MultigridSolver_LevelCycle( MultigridSolver* self, unsigned levelInd, Vector* rhs, Vector* solution );
	void MultigridSolver_RestrictMatrix( MultigridSolver* self, MultigridSolver_Level* level, Mat* dstMatrix );
	void MultigridSolver_LevelCycle( MultigridSolver* self, unsigned levelInd, Vec rhs, Vec solution );
	void MultigridSolver_UpdateWorkVectors( MultigridSolver* self );
	void MultigridSolver_UpdateSolvers( MultigridSolver* self );
	void MultigridSolver_UpdateMatrices( MultigridSolver* self );
	void MultigridSolver_UpdateOps( MultigridSolver* self );

	//MultigridSolver* MultigridSolver_CreateOuterSolver( MultigridSolver* self, Matrix* matrix );
	//MultigridSolver* MultigridSolver_CreateSmoother( MultigridSolver* self, Matrix* matrix );
	//MultigridSolver* MultigridSolver_CreateCoarseSolver( MultigridSolver* self, Matrix* matrix );
	MGSolver_PETScData* MultigridSolver_CreateOuterSolver( MultigridSolver* self, Mat matrix );
	MGSolver_PETScData* MultigridSolver_CreateSmoother( MultigridSolver* self, Mat matrix );
	MGSolver_PETScData* MultigridSolver_CreateCoarseSolver( MultigridSolver* self, Mat matrix );

	void MultigridSolver_Destruct( MultigridSolver* self );
	void MultigridSolver_DestructLevels( MultigridSolver* self );

#endif /* __StgFEM_SLE_SystemSetup_MultigridSolver_h__ */

