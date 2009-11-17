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
**	Sets up and manages a Energy system of linear equations.
**
** Assumptions:
**
** Comments:
**
** $Id: Energy_SLE.h 654 2006-10-12 08:58:49Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __SLE_ProvidedSystems_Energy_SLE_h__
#define __SLE_ProvidedSystems_Energy_SLE_h__

	/** Textual name of this class */
	extern const Type Energy_SLE_Type;

	/** Energy_SLE class contents */
	#define __Energy_SLE \
		/* General info */ \
		__SystemLinearEquations \
		\
		/* Virtual info */ \
		\
		/* Energy_SLE info */ \
		\
		StiffnessMatrix*	stiffMat; \
		SolutionVector*	solutionVec; \
		ForceVector*		fVector; 

	struct Energy_SLE { __Energy_SLE };	

	Energy_SLE* Energy_SLE_New( 		
		Name							name,
		SLE_Solver*					solver,
		FiniteElementContext*	context,
		Bool							isNonLinear,
		double						nonLinearTolerance,
		Iteration_Index			nonLinearMaxIterations,
		Bool							killNonConvergent,
		EntryPoint_Register*		entryPoint_Register,
		MPI_Comm						comm,
		StiffnessMatrix*			stiffMat,
		SolutionVector*			solutionVec,
		ForceVector*				fVector );

	/* Creation implementation / Virtual constructor */
	Energy_SLE* _Energy_SLE_New(
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
		SystemLinearEquations_LM_SetupFunction*             _LM_Setup,
		SystemLinearEquations_MatrixSetupFunction*          _matrixSetup,
		SystemLinearEquations_VectorSetupFunction*          _vectorSetup,
		SystemLinearEquations_UpdateSolutionOntoNodesFunc*	_updateSolutionOntoNodes, 
		SystemLinearEquations_MG_SelectStiffMatsFunc*		_mgSelectStiffMats, 
		Name                                                name );

	void _Energy_SLE_Init( 
		void*					sle,
		StiffnessMatrix*	stiffMat,
		SolutionVector*	solutionVec,
		ForceVector*		fVector );

	/** 'Class' Virtual Functions */
	void _Energy_SLE_Delete( void* stokesSleSolver );

	void _Energy_SLE_Print( void* stokesSleSolver, Stream* stream );

	/* 'Stg_Component' Virtual Functions */
	void* _Energy_SLE_DefaultNew( Name name );

	void _Energy_SLE_AssignFromXML( void* sle, Stg_ComponentFactory* cf, void* data );

#endif
