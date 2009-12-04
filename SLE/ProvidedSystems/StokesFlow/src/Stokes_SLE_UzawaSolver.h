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
**	Solves a Stokes SLE using the Uzawa-Preconditioned Conjugate Gradient method.
**
** Assumptions:
**
** Comments:
**	If C isn't active in the SLE we are solving, we assume its an incompressible system, and adjust accordingly.
**	Temporary Vector and MatrixSolver objects are re-used between iterations, and are only deleted once
**	Class_Delete() is called on this class. If something else in the code is using lots of memory temporarily,
**	these could be allocated and deleted just before and after the solve function each timestep.
**
**	TODO: Dave needs to find out if its ok to reuse PETSc solvers between timesteps if the matrix sizes change
**	between timesteps due to changes in BCs etc.
**
**	TODO: Need to handle Non-symmetric problems where Divergence isn't just G transpose. This will require a 
**	fn pointer or entry point to MatMult certain things based on D.
**
** $Id: Stokes_SLE_UzawaSolver.h 1125 2008-05-12 14:22:02Z DavidMay $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __StgFEM_SLE_ProvidedSystems_StokesFlow_Stokes_SLE_UzawaSolver_h__
#define __StgFEM_SLE_ProvidedSystems_StokesFlow_Stokes_SLE_UzawaSolver_h__

	/** Textual name of this class */
	extern const Type Stokes_SLE_UzawaSolver_Type;

	#define __Stokes_SLE_UzawaSolver \
		/* General info */ \
		__SLE_Solver \
		/* Virtual info */ \
		\
		/* Stokes_SLE_UzawaSolver info */ \
		\
		/* Matrix solvers */ \
		KSP                 velSolver;  /** Solver for the velocity system */ \
		KSP                 pcSolver;   /** Solver for the preconditioner system */\
		/* Auxiliary matrixes */ \
		StiffnessMatrix*    preconditioner; 	\
		/* auxillary vectors */ \
		Vec                 pTempVec; \
		Vec                 rVec; \
		Vec                 sVec; \
		Vec                 fTempVec; \
		Vec                 vStarVec; \
		/* parameters */ \
		Iteration_Index     maxUzawaIterations;	/** Max number of iterations the solver can run for */ \
		Iteration_Index     minUzawaIterations; \
		double              tolerance;	/** Tolerance that residual must be within to finish */ \
		Bool                useAbsoluteTolerance; \
		Bool		    monitor;

	/** Solves a Stokes SLE using the Uzawa-Preconditioned Conjugate Gradient method - see Stokes_SLE_UzawaSolver.h */
	struct Stokes_SLE_UzawaSolver { __Stokes_SLE_UzawaSolver };
	
	/* --- Constructors / Destructor --- */

	/** Constructor */
	void* _Stokes_SLE_UzawaSolver_DefaultNew( Name name );
	
	Stokes_SLE_UzawaSolver* Stokes_SLE_UzawaSolver_New( 
		Name                                        name,
		Bool                                        useStatSolve, 
		int                                         statReps,
		StiffnessMatrix*                            preconditioner,
		Iteration_Index                             maxUzawaIterations,
		Iteration_Index                             minUzawaIterations,
		double                                      tolerance,
		Bool                                        useAbsoluteTolerance,
		Bool                                        monitor );
		
	/** Creation implementation / Virtual constructor */
	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define STOKES_SLE_UZAWASOLVER_DEFARGS \
                SLE_SOLVER_DEFARGS

	#define STOKES_SLE_UZAWASOLVER_PASSARGS \
                SLE_SOLVER_PASSARGS

	Stokes_SLE_UzawaSolver* _Stokes_SLE_UzawaSolver_New(  STOKES_SLE_UZAWASOLVER_DEFARGS  );

	/** Class member variable initialisation */
	void _Stokes_SLE_UzawaSolver_Init(
		Stokes_SLE_UzawaSolver*      self,
		StiffnessMatrix*             preconditioner, 
		Iteration_Index              maxUzawaIterations,
		Iteration_Index              minUzawaIterations,
		double                       tolerance,
		Bool                         useAbsoluteTolerance, 
                Bool                         monitor );
		
	void Stokes_SLE_UzawaSolver_InitAll( 
		void*                        solver,
		Bool                         useStatSolve,
		int                          statReps, 
		StiffnessMatrix*             preconditioner, 
		Iteration_Index              maxUzawaIterations,
		Iteration_Index              minUzawaIterations,
		double                       tolerance,
		Bool                         useAbsoluteTolerance,
		Bool                         monitor );
		
	/** Class_Delete() implementation */
	void _Stokes_SLE_UzawaSolver_Delete( void* solver );

	/* --- Virtual Function Implementations --- */
	void _Stokes_SLE_UzawaSolver_Print( void* solver, Stream* stream );
	
	/* Copy */
	#define Stokes_SLE_UzawaSolver_Copy( self ) \
		(SLE_Solver*)Stg_Class_Copy( self, NULL, False, NULL, NULL )
	#define Stokes_SLE_UzawaSolver_DeepCopy( self ) \
		(SLE_Solver*)Class_DeepCopy( self, NULL, True, NULL, NULL )
	
	void* _Stokes_SLE_UzawaSolver_Copy( void* stokesSleUzawaSolver, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );

	/** Stg_Component_Build() implementations: allocates the 2 MatrixSolvers and additional Vectors */
	void _Stokes_SLE_UzawaSolver_Build( void* solver, void* stokesSLE );
	
	void _Stokes_SLE_UzawaSolver_AssignFromXML( void* solver, Stg_ComponentFactory* cf, void* data );
	
	void _Stokes_SLE_UzawaSolver_Initialise( void* solver, void* stokesSLE ) ;
	
	void _Stokes_SLE_UzawaSolver_Execute( void* solver, void* data );
	
	void _Stokes_SLE_UzawaSolver_Destroy( void* solver, void* data );

	/** SolverSetup: sets up the 2 MatrixSolvers */
	void _Stokes_SLE_UzawaSolver_SolverSetup( void* stokesSle, void* stokesSLE );

	/** Solve: - Applies Preconditioned Conjuagte Gradient - Pressure Corrected Algorithm. 
	REFERENCES:
 	* Ramage A. & Wathen J., "Iterative solution techniques for the stokes and navier-stokes equations", 
	Int. J. Num. Meth. Fluids, 19, (1994), pp 67-83	(See figure 1)

	* Atanga J. & Silvester D., "Iterative methods for stabilized mixed velocity-pressure finite elements",
	Int. J. Num. Meth. Fluids, 14, (1992), pp 71-81.	(See algorithm 2)
	
	This routine should NOT modify F or H!
	H should not be modified because G is only build once and H is derived from GTrans * bc values
	F cannot be modified because this func. is called multiple times by Uzawa manager and K is not
	rebuilt in the between each call to this function. 
	*/
	void _Stokes_SLE_UzawaSolver_Solve( void* solver, void* stokesSLE );

	/* Get residual implementation */
	//Vector* _Stokes_SLE_UzawaSolver_GetResidual( void* solver, Index fv_I );
	Vec _Stokes_SLE_UzawaSolver_GetResidual( void* solver, Index fv_I );

	//void _Stokes_SLE_UzawaSolver_FormResidual( void *stokesSLE, void *solver, Vector *stg_r );
	//void _Stokes_SLE_UzawaSolver_GetRhs( void *stokesSLE, void *solver, Vector *stg_rhs );
	//void _Stokes_SLE_UzawaSolver_GetSolution( void *stokesSLE, void *solver, Vector **x );
	void _Stokes_SLE_UzawaSolver_FormResidual( void *stokesSLE, void *solver, Vec r );
	void _Stokes_SLE_UzawaSolver_GetRhs( void *stokesSLE, void *solver, Vec rhs );
	void _Stokes_SLE_UzawaSolver_GetSolution( void *stokesSLE, void *solver, Vec *x );

#endif	

