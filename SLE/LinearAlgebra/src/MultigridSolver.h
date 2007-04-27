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

#ifndef __StgFEM_SLE_LinearAlgebra_MultigridSolver_h__
#define __StgFEM_SLE_LinearAlgebra_MultigridSolver_h__

	/** Textual name of this class */
	extern const Type MultigridSolver_Type;

	/** Virtual function types */

	/** MultigridSolver class contents */
	typedef struct {
		MatrixSolver*	downSolver;
		unsigned	nDownIts;
		MatrixSolver*	upSolver;
		unsigned	nUpIts;
		unsigned	nCycles;

		Matrix*		A;
		Matrix*		R;
		Matrix*		P;

		Vector*		workRHS;
		Vector*		workSol;
	} MultigridSolver_Level;

	#define __MultigridSolver			\
		/* General info */			\
		__MatrixSolver				\
							\
		/* Virtual info */			\
							\
		/* MultigridSolver info */		\
		Stream*			stream;		\
							\
		unsigned		nLevels;	\
		MultigridSolver_Level*	levels;		\
		MGOpGenerator*		opGen;		\
		Bool			solversChanged;	\
		Bool			opsChanged;	\
							\
		unsigned		curIt;		\
		unsigned		maxIts;		\
		double			relTol;		\
		double			rnorm;		\
		Bool			useInitial;	\
		MatrixSolver*		outerSolver;

	struct MultigridSolver { __MultigridSolver };

	/*--------------------------------------------------------------------------------------------------------------------------
	** Constructors
	*/

	#define MULTIGRIDSOLVER_DEFARGS	\
		MATRIXSOLVER_DEFARGS

	#define MULTIGRIDSOLVER_PASSARGS \
		MATRIXSOLVER_PASSARGS

	MultigridSolver* MultigridSolver_New( Name name );
	MultigridSolver* _MultigridSolver_New( MATRIXSOLVER_DEFARGS );
	void _MultigridSolver_Init( MultigridSolver* self );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Virtual functions
	*/

	void _MultigridSolver_Delete( void* matrixSolver );
	void _MultigridSolver_Print( void* matrixSolver, Stream* stream );
	void _MultigridSolver_Construct( void* matrixSolver, Stg_ComponentFactory* cf, void* data );
	void _MultigridSolver_Build( void* matrixSolver, void* data );
	void _MultigridSolver_Initialise( void* matrixSolver, void* data );
	void _MultigridSolver_Execute( void* matrixSolver, void* data );
	void _MultigridSolver_Destroy( void* matrixSolver, void* data );

	void MultigridSolver_SetMaxIterations( void* matrixSolver, unsigned nIterations );
	void MultigridSolver_SetRelativeTolerance( void* matrixSolver, double tolerance );
	void MultigridSolver_SetAbsoluteTolerance( void* matrixSolver, double tolerance );
	void MultigridSolver_SetUseInitialSolution( void* matrixSolver, Bool state );

	void MultigridSolver_Solve( void* matrixSolver, void* rhs, void* solution );
	void MultigridSolver_Setup( void* matrixSolver, void* rhs, void* solution );

	MatrixSolver_Status MultigridSolver_GetSolveStatus( void* matrixSolver );
	unsigned MultigridSolver_GetIterations( void* matrixSolver );
	unsigned MultigridSolver_GetMaxIterations( void* matrixSolver );
	double MultigridSolver_GetResidualNorm( void* matrixSolver );

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
	Matrix* MultigridSolver_GetRestriction( void* matrixSolver, unsigned level );
	Matrix* MultigridSolver_GetProlongation( void* matrixSolver, unsigned level );

	/*--------------------------------------------------------------------------------------------------------------------------
	** Private Member functions
	*/

	void MultigridSolver_RestrictMatrix( MultigridSolver* self, MultigridSolver_Level* level, Matrix** dstMatrix );
	void MultigridSolver_LevelCycle( MultigridSolver* self, unsigned levelInd, Vector* rhs, Vector* solution );
	void MultigridSolver_UpdateWorkVectors( MultigridSolver* self );
	void MultigridSolver_UpdateSolvers( MultigridSolver* self );
	void MultigridSolver_UpdateMatrices( MultigridSolver* self );
	void MultigridSolver_UpdateOps( MultigridSolver* self );

	MatrixSolver* MultigridSolver_CreateOuterSolver( MultigridSolver* self, Matrix* matrix );
	MatrixSolver* MultigridSolver_CreateSmoother( MultigridSolver* self, Matrix* matrix );
	MatrixSolver* MultigridSolver_CreateCoarseSolver( MultigridSolver* self, Matrix* matrix );

	void MultigridSolver_Destruct( MultigridSolver* self );
	void MultigridSolver_DestructLevels( MultigridSolver* self );

#endif /* __StgFEM_SLE_LinearAlgebra_MultigridSolver_h__ */
