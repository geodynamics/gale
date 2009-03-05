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


#ifndef __StgFEM_AdvectionDiffusion_Multicorrector_h__
#define __StgFEM_AdvectionDiffusion_Multicorrector_h__
	
	/** Textual name of this class */
	extern const Type AdvDiffMulticorrector_Type;

	/** AdvDiffMulticorrector class contents */
	#define __AdvDiffMulticorrector \
		/* General info */ \
		__SLE_Solver \
		\
		/* Virtual info */ \
		\
		/* AdvDiffMulticorrector info */ \
		double                                              gamma;                          \
		Iteration_Index                                     multiCorrectorIterations;       \
		KSP                                                 matrixSolver;


	struct AdvDiffMulticorrector { __AdvDiffMulticorrector };	

	AdvDiffMulticorrector* AdvDiffMulticorrector_New( 
		Name                                                name,
		double                                              gamma,
		Iteration_Index                                     multiCorrectorIterations );

	AdvDiffMulticorrector* _AdvDiffMulticorrector_New( 
		SizeT                                               sizeOfSelf,  
		Type                                                type,
		Stg_Class_DeleteFunction*                           _delete,
		Stg_Class_PrintFunction*                            _print,
		Stg_Class_CopyFunction*                             _copy, 
		Stg_Component_DefaultConstructorFunction*           _defaultConstructor,
		Stg_Component_ConstructFunction*                    _construct,
		Stg_Component_BuildFunction*                        _build,
		Stg_Component_InitialiseFunction*                   _initialise,
		Stg_Component_ExecuteFunction*                      _execute,
		Stg_Component_DestroyFunction*                      _destroy,
		SLE_Solver_SolverSetupFunction*                     _solverSetup,
		SLE_Solver_SolveFunction*                           _solve,
		SLE_Solver_GetResidualFunc*                         _getResidual, 
		Name                                                name );

	void _AdvDiffMulticorrector_Init( 
		AdvDiffMulticorrector*                              self, 
		double                                              gamma,
		Iteration_Index                                     multiCorrectorIterations );
	
	void AdvDiffMulticorrector_InitAll( 
		void*                                               solver,
		double                                              gamma,
		Iteration_Index                                     multiCorrectorIterations );

	void _AdvDiffMulticorrector_Delete( void* solver );
	void _AdvDiffMulticorrector_Print( void* solver, Stream* stream );

	void* _AdvDiffMulticorrector_DefaultNew( Name name ) ;
	void _AdvDiffMulticorrector_Construct( void* solver, Stg_ComponentFactory* cf, void* data ) ;
	void _AdvDiffMulticorrector_Build( void* solver, void* data ) ;
	void _AdvDiffMulticorrector_Initialise( void* solver, void* data ) ;
	void _AdvDiffMulticorrector_Execute( void* solver, void* data ) ;
	void _AdvDiffMulticorrector_Destroy( void* solver, void* data ) ;

	void _AdvDiffMulticorrector_SolverSetup( void* solver, void* data ) ;
	void _AdvDiffMulticorrector_Solve( void* solver, void* _sle ) ;

	void AdvDiffMulticorrector_Predictors( AdvDiffMulticorrector* self, AdvectionDiffusionSLE* sle, double dt ) ;

	void AdvDiffMulticorrector_Solution( AdvDiffMulticorrector* self, AdvectionDiffusionSLE* sle, Vec deltaPhiDot ) ;
	void AdvDiffMulticorrector_Correctors( AdvDiffMulticorrector* self, AdvectionDiffusionSLE* sle, Vec deltaPhiDot, double dt ) ;

	void AdvDiffMulticorrector_CalculatePhiDot( AdvDiffMulticorrector* self, AdvectionDiffusionSLE* sle, Vec deltaPhiDot ) ;
	void _AdvDiffMulticorrector_CalculatePhiDot_Explicit( AdvDiffMulticorrector* self, AdvectionDiffusionSLE* sle, Vec deltaPhiDot ) ;
	void _AdvDiffMulticorrector_CalculatePhiDot_Implicit( AdvDiffMulticorrector* self, AdvectionDiffusionSLE* sle, Vec deltaPhiDot ) ;

#endif
