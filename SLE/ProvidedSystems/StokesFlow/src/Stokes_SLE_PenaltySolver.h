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
**	Solves a Stokes SLE using the Penalty method.
**
** Assumptions:
**
** Comments:
**
** $Id: Stokes_SLE_PenaltySolver.h 822 2007-04-27 06:20:35Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __FE_Stokes_SLE_PenaltySolver_h__
#define __FE_Stokes_SLE_PenaltySolver_h__

	/** Textual name of this class */
	extern const Type Stokes_SLE_PenaltySolver_Type;

	#define __Stokes_SLE_PenaltySolver \
		/* General info */ \
		__SLE_Solver \
		\
		/* Virtual info */ \
		\
		/* Stokes_SLE_PenaltySolver info */ \
		\
		/* Matrix solvers */ \

	/** Solves a Stokes SLE using the Penalty Method */
	struct Stokes_SLE_PenaltySolver { __Stokes_SLE_PenaltySolver };
	
	/* --- Constructors / Destructor --- */

	/** Constructor */
	void* Stokes_SLE_PenaltySolver_DefaultNew( Name name );
	
	Stokes_SLE_PenaltySolver* Stokes_SLE_PenaltySolver_New(
		Name                                        name,
		Bool                                        useStatSolve, 
		int                                         statReps );

	/** Creation implementation / Virtual constructor */
	Stokes_SLE_PenaltySolver* _Stokes_SLE_PenaltySolver_New( 
		SizeT                                       sizeOfSelf,
		Type                                        type,
		Stg_Class_DeleteFunction*                   _delete,
		Stg_Class_PrintFunction*                    _print,
		Stg_Class_CopyFunction*                     _copy, 
		Stg_Component_DefaultConstructorFunction*   _defaultConstructor,
		Stg_Component_ConstructFunction*            _construct,
		Stg_Component_BuildFunction*                _build,
		Stg_Component_InitialiseFunction*           _initialise,
		Stg_Component_ExecuteFunction*              _execute,
		Stg_Component_DestroyFunction*              _destroy,
		SLE_Solver_SolverSetupFunction*             _solverSetup,
		SLE_Solver_SolveFunction*                   _solve,
		SLE_Solver_GetResidualFunc*                 _getResidual, 
		Name                                        name );

	/** Class member variable initialisation */
	void _Stokes_SLE_PenaltySolver_Init( void* solver ) ;
	void Stokes_SLE_PenaltySolver_InitAll( 
		void*                        solver,
		Bool                         useStatSolve,
		int                          statReps );

	/** Class_Delete() implementation */
	void _Stokes_SLE_PenaltySolver_Delete( void* solver );

	/* --- Virtual Function Implementations --- */
	void _Stokes_SLE_PenaltySolver_Print( void* solver, Stream* stream );
	
	/* Copy */
	#define Stokes_SLE_PenaltySolver_Copy( self ) \
		(SLE_Solver*)Stg_Class_Copy( self, NULL, False, NULL, NULL )
	#define Stokes_SLE_PenaltySolver_DeepCopy( self ) \
		(SLE_Solver*)Class_DeepCopy( self, NULL, True, NULL, NULL )
	
	void* _Stokes_SLE_PenaltySolver_Copy( void* stokesSlePenaltySolver, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );

	/** Stg_Component_Build() implementations: allocates the 2 MatrixSolvers and additional Vectors */
	void _Stokes_SLE_PenaltySolver_Build( void* solver, void* stokesSLE );
	
	void _Stokes_SLE_PenaltySolver_Construct( void* solver, Stg_ComponentFactory* cf, void* data );
	
	void _Stokes_SLE_PenaltySolver_Initialise( void* solver, void* stokesSLE ) ;
	
	void _Stokes_SLE_PenaltySolver_Execute( void* solver, void* data );
	
	void _Stokes_SLE_PenaltySolver_Destroy( void* solver, void* data );

	/** SolverSetup: sets up the 2 MatrixSolvers */
	void _Stokes_SLE_PenaltySolver_SolverSetup( void* stokesSle, void* stokesSLE );

	/**
	Solves
	Ku + Grad p = f
	Div u + C p = h
	by eliminating pressure via the penalty method.

	Hence we obtain the velocity by solving
	(K - Grad CInv Div )u = kHat u = f - Grad CInv h
	and recover pressure from
	p = CInv( h - Div u )
	*/
	void _Stokes_SLE_PenaltySolver_Solve( void* solver, void* stokesSLE );

	/* Get residual implementation */
	Vector* _Stokes_SLE_PenaltySolver_GetResidual( void* solver, Index fv_I );

#endif	

