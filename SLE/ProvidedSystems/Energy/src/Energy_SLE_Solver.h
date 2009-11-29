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
**	Solves a basic SLE consisting of only one matrix, one force vector and one soln vector.**
** Assumptions:
**
** Comments:
**	Uses a straightforward MatrixSolver to solve the system. Can be used for
**	problems such as temperature diffusion etc. Uses the general interface
**	to StiffnessMatrices, ForceVectors and SolutionVectors provided by the
**	SystemLinearEquations class.
**
** $Id: Energy_SLE_Solver.h 822 2007-04-27 06:20:35Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __StgFEM_SLE_ProvidedSystems_Energy_Energy_SLE_Solver_h__
#define __StgFEM_SLE_ProvidedSystems_Energy_Energy_SLE_Solver_h__

	/** Textual name of this class */
	extern const Type Energy_SLE_Solver_Type;
	
	/** Energy_SLE_Solver class contents */
	#define __Energy_SLE_Solver \
		/* General info */ \
		__SLE_Solver \
		\
		/* Virtual info */ \
		\
		/* Energy_SLE_Solver info */ \
		KSP			matrixSolver; \
		Vec    			residual; \

	/** Solves a basic SLE consisting of only one matrix, one force vector and one soln vector - see
	Energy_SLE_Solver.h */
	struct Energy_SLE_Solver { __Energy_SLE_Solver };	

	/** Constructor */
	void* Energy_SLE_Solver_DefaultNew( Name name );
	
	Energy_SLE_Solver* Energy_SLE_Solver_New( Name name, Bool useStatSolve, int statReps ) ;

	/** Creation implementation / Virtual constructor */
	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define ENERGY_SLE_SOLVER_DEFARGS \
                SLE_SOLVER_DEFARGS

	#define ENERGY_SLE_SOLVER_PASSARGS \
                SLE_SOLVER_PASSARGS

	Energy_SLE_Solver* _Energy_SLE_Solver_New(  ENERGY_SLE_SOLVER_DEFARGS  );
	
	/** Member variable initialisation */
	void _Energy_SLE_Solver_Init( Energy_SLE_Solver* self );
	void Energy_SLE_Solver_InitAll( Energy_SLE_Solver* solver, Bool useStatSolve, int statReps ) ;
	
	/** Class_Delete() implementation */
	void _Energy_SLE_Solver_Delete( void* sleSolver );
	
	/* --- Virtual Function Implementations --- */

	/** Stg_Class_Print() implementation */
	void _Energy_SLE_Solver_Print( void* sleSolver, Stream* stream );
	
	/* Copy */
	#define Energy_SLE_Solver_Copy( self ) \
		(Energy_SLE_Solver*)Stg_Class_Copy( self, NULL, False, NULL, NULL )
	#define Energy_SLE_Solver_DeepCopy( self ) \
		(Energy_SLE_Solver*)Stg_Class_Copy( self, NULL, True, NULL, NULL )
	
	void* _Energy_SLE_Solver_Copy( void* standardSleSolver, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );

	/* Stg_Component_Build() implementation */
	void _Energy_SLE_Solver_Build( void* sleSolver, void* standardSLE );
	
	void _Energy_SLE_Solver_Initialise( void* sleSolver, void* standardSLE );
	
	void _Energy_SLE_Solver_AssignFromXML( void* sleSolver, Stg_ComponentFactory* cf, void* data );
	
	void _Energy_SLE_Solver_Execute( void* sleSolver, void* data );
	
	void _Energy_SLE_Solver_Destroy( void* sleSolver, void* data );

	/* SLE_Solver_SolverSetup() implementation. */
	void _Energy_SLE_Solver_SolverSetup( void* sleSolver, void* standardSLE );

	/* SLE_Solver_Solve() implementation */
	void _Energy_SLE_Solver_Solve( void* sleSolver, void* standardSLE );

	/* Get residual implementation */
	//Vector* _Energy_SLE_GetResidual( void* sleSolver, Index fv_I );
	Vec _Energy_SLE_GetResidual( void* sleSolver, Index fv_I );

#endif

