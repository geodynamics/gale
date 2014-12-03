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
**	Abstract class defining the interface for a System of Linear Equations solver.
**
** Assumptions:
**
** Comments:
**	Note - as 1 September 2004 (rev 1994), the functioality for building Matrices etc
**	that was in this class is back in the SLE_Solver class.
**
** $Id: SLE_Solver.h 1125 2008-05-12 14:22:02Z DavidMay $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
	
#ifndef __StgFEM_SLE_SystemSetup_SLE_Solver_h__
#define __StgFEM_SLE_SystemSetup_SLE_Solver_h__
	

	/** Textual name of this class */
	extern const Type SLE_Solver_Type;
	
	/* virtual function interface */
	typedef void (SLE_Solver_SolverSetupFunction) ( void* sleSolver, void* sle );
	typedef void (SLE_Solver_SolveFunction) ( void* sleSolver, void* sle );
	typedef Vec (SLE_Solver_GetResidualFunc) ( void* sleSolver, Index fvIndex );

	typedef void (SLE_Solver_FormResidualFunc) ( void *someSLE, void *someSolver, Vec );
	typedef void (SLE_Solver_GetRhsFunc) ( void *someSLE, void *someSovler, Vec );
	typedef void (SLE_Solver_GetSolutionFunc) ( void *someSLE, void *someSolver, Vec* );

	/** SLE_Solver class contents */
	#define __SLE_Solver \
		__Stg_Component \
		FiniteElementContext*				context; \
		ExtensionManager*						extensionManager; \
		\
		/* Virtual info */ \
		SLE_Solver_SolverSetupFunction*	_solverSetup; \
		SLE_Solver_SolveFunction*			_solve; \
		SLE_Solver_GetResidualFunc*		_getResidual; \
  		SLE_Solver_FormResidualFunc*		_formResidual; \
		SLE_Solver_GetRhsFunc*				_getRhs; \
		SLE_Solver_GetSolutionFunc*		_getSolution; \
		\
		/* SLE_Solver info */ \
		Stream*                            debug; \
		Stream*                            info; \
		Iteration_Index                    maxIterations; \
                                                                  \
		/* Timing variables for solvers */ \
		double				inneritsinitialtime; \
		double				outeritsinitialtime; \
		double				nonlinearitsinitialtime; \
		double				inneritsendtime; \
		double				outeritsendtime; \
		double				nonlinearitsendtime; \
		double				totalinneritstime; \
		double				totalouteritstime; \
		double				totalnonlinearitstime; \
		int					totalnuminnerits; \
		int				totalnumouterits; \
		int				totalnumnonlinearits; \
		int				avgnuminnerits; \
		int				avgnumouterits; \
		double			avgtimeinnerits; \
		double			avgtimeouterits; \
		double			avgtimenonlinearits; \
		int				currenttimestep; \
		int				previoustimestep; \
		\
		Bool                               useStatSolve; \
		unsigned                           nStatReps;
		
	/** Abstract class defining the interface for a SLE_Solver solver - see SLE_Solver.h */
	struct SLE_Solver { __SLE_Solver };

	/* No SLE_Solver_New() or SLE_Solver_Init() as this is an abstract class. */

	/* --- Constructor functions --- */

	/** Creation implementation */
	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define SLE_SOLVER_DEFARGS \
                STG_COMPONENT_DEFARGS, \
                SLE_Solver_SolverSetupFunction*  _solverSetup, \
                SLE_Solver_SolveFunction*              _solve, \
                SLE_Solver_GetResidualFunc*      _getResidual

	#define SLE_SOLVER_PASSARGS \
                STG_COMPONENT_PASSARGS, \
	        _solverSetup, \
	        _solve,       \
	        _getResidual

	SLE_Solver* _SLE_Solver_New(  SLE_SOLVER_DEFARGS  );

	/** class member initialisation */
	void _SLE_Solver_Init( SLE_Solver* self, Bool useStatSolve, int statReps ) ;

	void SLE_Solver_InitAll( void* sleSolver, Bool useStatSolve, int statReps ) ;

	/* --- Virtual function implementations --- */
	
	/** Class Virtual Functions Implementations */
	void* _SLE_Solver_Copy( const void* sleSolver, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );

	void _SLE_Solver_Delete( void* sleSolver );

	void _SLE_Solver_Print( void* sleSolver, Stream* stream ) ;

	/** Stg_Component_Build() implementation: does nothing by default as some solvers may not need it. */
	void _SLE_Solver_Build( void* sleSolver, void* data );
	
	void _SLE_Solver_AssignFromXML( void* sleSolver, Stg_ComponentFactory* cf, void* data );

	/** Stg_Component_Initialise() implementation: does nothing by default as some solvers may not neet it. */
	void _SLE_Solver_Initialise( void* sleSolver, void* data );

	/** Stg_Component_Execute() implementation: Calls SolverSetup() for any per-iteration setup, then 
	calls Solve() to calculate the new solutions. */
	void _SLE_Solver_Execute( void* sleSolver, void* data );
	
	void _SLE_Solver_Destroy( void* sleSolver, void* data );

	/* --- Public functions --- */

	/** Does any required solver setup beyond assembly of the matrices to be solved: e.g. priming the Matrix solvers
	etc. */
	void SLE_Solver_SolverSetup( void* sleSolver, void* sle );
	
	/** Solve:- calculate the new values for all solution vectors in the system. */
	void SLE_Solver_Solve( void* sleSolver, void* sle );

	/** Set the maximum number of iterations */
	#define SLE_Solver_SetMaxIterations( self, its ) \
		(self)->maxIterations = its

	/** Get the systems most recent residual */
	#define SLE_Solver_GetResidual( self, fvIndex ) \
		(self)->_getResidual( self, fvIndex )
	

#endif

