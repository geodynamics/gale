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
**		Public types for the FiniteElement module.
**
** Assumptions:
**
** Comments:
**
** $Id: types.h 1210 2008-08-25 01:17:12Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __StgFEM_SLE_SystemSetup_types_h__
#define __StgFEM_SLE_SystemSetup_types_h__
	
	/* FE types/classes */
	typedef struct FeEntryPoint                 FeEntryPoint;
	typedef struct StiffnessMatrix              StiffnessMatrix;
	typedef struct StiffnessMatrixTerm          StiffnessMatrixTerm;
	typedef struct ForceVector                  ForceVector;
	typedef struct ForceTerm                    ForceTerm;
	typedef struct SolutionVector               SolutionVector;
	typedef struct SystemLinearEquations        SystemLinearEquations;
	typedef struct SLE_Solver                   SLE_Solver;
	typedef struct FiniteElementContext         FiniteElementContext;
	typedef struct Assembler		Assembler;
        typedef struct MultigridSolver			MultigridSolver;
        typedef struct MGOpGenerator			MGOpGenerator;
        typedef struct SROpGenerator			SROpGenerator;
        typedef struct PETScMGSolver			PETScMGSolver;
	
	/* types for lists etc ... for readability */
	typedef Index                       StiffnessMatrix_Index;
	typedef StiffnessMatrix*            StiffnessMatrixPtr;
	typedef Stg_ObjectList              StiffnessMatrixList;
	typedef Index                       ForceVector_Index;
	typedef ForceVector*                ForceVectorPtr;
	typedef Stg_ObjectList              ForceVectorList;
	typedef Index                       SolutionVector_Index;
	typedef SolutionVector*             SolutionVectorPtr;
	typedef Stg_ObjectList              SolutionVectorList;
	typedef Index                       SLE_Solver_Index;
	typedef SolutionVector*             SLE_SolverPtr;
	typedef Stg_ObjectList              SLE_SolverList;
	typedef Index                       SystemLinearEquations_Index;
	typedef SystemLinearEquations*      SystemLinearEquationsPtr;
	typedef Stg_ObjectList              SystemLinearEquationList;

	/* output streams: initialised in StgFEM_SLE_SystemSetup_Init() */
	extern Stream* StgFEM_SLE_Debug;
	extern Stream* StgFEM_SLE_SystemSetup_Debug;

typedef struct {
      void* callback;
      void* object;
} Callback;

typedef enum {
	MGSolver_Status_ConvergedRelative = 2, 
	MGSolver_Status_ConvergedAbsolute = 3, 
	MGSolver_Status_ConvergedIterations = 4, 
	MGSolver_Status_DivergedNull = -2, 
	MGSolver_Status_DivergedIterations = -3, 
	MGSolver_Status_DivergedTolerance = -4, 
	MGSolver_Status_Iterating = 0
} MGSolver_Status;

/** Virtual function types */
typedef void (MGSolver_SetCommFunc)( void* matrixSolver, MPI_Comm comm );
typedef void (MGSolver_SetMatrixFunc)( void* matrixSolver, void* matrix );
typedef void (MGSolver_SetMaxIterationsFunc)( void* matrixSolver, unsigned nIterations );
typedef void (MGSolver_SetRelativeToleranceFunc)( void* matrixSolver, double tolerance );
typedef void (MGSolver_SetAbsoluteToleranceFunc)( void* matrixSolver, double tolerance );
typedef void (MGSolver_SetUseInitialSolutionFunc)( void* matrixSolver, Bool state );

typedef void (MGSolver_SolveFunc)( void* matrixSolver, void* rhs, void* solution );
typedef void (MGSolver_SetupFunc)( void* matrixSolver, void* rhs, void* solution );

typedef MGSolver_Status (MGSolver_GetSolveStatusFunc)( void* matrixSolver );
typedef unsigned (MGSolver_GetIterationsFunc)( void* matrixSolver );
typedef unsigned (MGSolver_GetMaxIterationsFunc)( void* matrixSolver );
typedef double (MGSolver_GetResidualNormFunc)( void* matrixSolver );


/* MatrixSolver class has been depreciated, so this class can no
 * longer inherit from it. all the data previously encapsulated in 
 * the MatrixSolver class is now wrapped up here */
typedef struct {
	MPI_Comm	comm;
	KSP		ksp;
	Mat		matrix;
	Mat		inversion;
	Vec		residual;
	Bool		expiredResidual;
	Bool		matrixChanged;
	Vec		curRHS;
	Vec		curSolution;
	Bool		optionsReady;
} MGSolver_PETScData;

#endif /* __StgFEM_SLE_SystemSetup_types_h__ */
