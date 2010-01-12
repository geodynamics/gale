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
** $Id: PETScMGSolver.h 3584 2006-05-16 11:11:07Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __StgFEM_SLE_SystemSetup_PETScMGSolver_h__
#define __StgFEM_SLE_SystemSetup_PETScMGSolver_h__

	/** Textual name of this class */
	extern const Type PETScMGSolver_Type;

	/** Virtual function types */

	/** PETScMGSolver class contents */
	typedef struct {
		unsigned	nDownIts;
		unsigned	nUpIts;
		unsigned	nCycles;

		Mat		R;
		Mat		P;
		Mat		A;

		Vec		workRes;
		Vec		workSol;
		Vec		workRHS;
	} PETScMGSolver_Level;

	#define __PETScMGSolver				\
		/* General info */			\
		__Stg_Component				\
							\
		/* Virtual info */			\
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
		/* PETScMGSolver info */		\
		Bool			pure;		\
		unsigned		nLevels;	\
		PETScMGSolver_Level*	levels;		\
		MGOpGenerator*		opGen;		\
		Bool			solversChanged;	\
		Bool			opsChanged;	\
		/* this stuff was previously stored in the */ \
		/* multigridSolver class, from which this inherited */ \
		MGSolver_PETScData*	mgData;

	struct PETScMGSolver { __PETScMGSolver };

	/*--------------------------------------------------------------------------------------------------------------------------
	** Constructors
	*/



	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define PETSCMGSOLVER_DEFARGS \
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

	#define PETSCMGSOLVER_PASSARGS \
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

	PETScMGSolver* PETScMGSolver_New( Name name );
	PETScMGSolver* _PETScMGSolver_New(  PETSCMGSOLVER_DEFARGS  );
	void _PETScMGSolver_Init( PETScMGSolver* self );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Virtual functions
	*/

	void _PETScMGSolver_Delete( void* matrixSolver );
	void _PETScMGSolver_Print( void* matrixSolver, Stream* stream );
	void _PETScMGSolver_AssignFromXML( void* matrixSolver, Stg_ComponentFactory* cf, void* data );
	void _PETScMGSolver_Build( void* matrixSolver, void* data );
	void _PETScMGSolver_Initialise( void* matrixSolver, void* data );
	void _PETScMGSolver_Execute( void* matrixSolver, void* data );
	void _PETScMGSolver_Destroy( void* matrixSolver, void* data );

	void PETScMGSolver_SetComm( void* matrixSolver, MPI_Comm comm );
	void PETScMGSolver_Setup( void* matrixSolver, void* rhs, void* solution );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Public functions
	*/

	void PETScMGSolver_SetLevels( void* matrixSolver, unsigned nLevels );
	void PETScMGSolver_SetRestriction( void* matrixSolver, unsigned levelInd, void* R );
	void PETScMGSolver_SetProlongation( void* matrixSolver, unsigned levelInd, void* P );
	void PETScMGSolver_SetLevelDownIterations( void* matrixSolver, unsigned level, unsigned nIts );
	void PETScMGSolver_SetLevelUpIterations( void* matrixSolver, unsigned level, unsigned nIts );
	void PETScMGSolver_SetLevelCycles( void* matrixSolver, unsigned level, unsigned nCycles );
	void PETScMGSolver_SetAllDownIterations( void* matrixSolver, unsigned nIts );
	void PETScMGSolver_SetAllUpIterations( void* matrixSolver, unsigned nIts );

	unsigned PETScMGSolver_GetNumLevels( void* matrixSolver );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Private Member functions
	*/

	void PETScMGSolver_UpdateOps( PETScMGSolver* self );
	void PETScMGSolver_UpdateMatrices( PETScMGSolver* self );
	void PETScMGSolver_UpdateWorkVectors( PETScMGSolver* self );
	void PETScMGSolver_UpdateSolvers( PETScMGSolver* self );

	void PETScMGSolver_Destruct( PETScMGSolver* self );
	void PETScMGSolver_DestructLevels( PETScMGSolver* self );

#endif /* __StgFEM_SLE_SystemSetup_PETScMGSolver_h__ */

