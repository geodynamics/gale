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
**	Sets up and manages a Stokes system of linear equations.
**
** Assumptions:
**
** Comments:
**	The particular solver algorithm must be provided through the Solver Stg_Component that the user chooses
**	More on the system we are solving:
**	Ku + Grad p = f
**	Div u + C p = h
**	If C isn't provided, we assume its an incompressible system, and adjust accordingly.
**
** $Id: Stokes_SLE.h 654 2006-10-12 08:58:49Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __FE_Stokes_SLE_h__
#define __FE_Stokes_SLE_h__

	/** Textual name of this class */
	extern const Type Stokes_SLE_Type;

	/** Stokes_SLE class contents */
	#define __Stokes_SLE \
		/* General info */ \
		__SystemLinearEquations \
		\
		/* Virtual info */ \
		\
		/* Stokes_SLE info */ \
		\
		StiffnessMatrix*    kStiffMat;	/** Stress tensor matrix */ \
		StiffnessMatrix*    gStiffMat;	/** Gradient matrix */ \
		StiffnessMatrix*    dStiffMat;	/** Divergence matrix */ \
		StiffnessMatrix*    cStiffMat;	/** Compressibility matrix */\
		SolutionVector*     uSolnVec;	/** velocity vector */\
		SolutionVector*     pSolnVec;	/** pressure vector */\
		ForceVector*        fForceVec;	/** forcing term vector */\
		ForceVector*        hForceVec;	/** continuity force vector */\

	struct Stokes_SLE { __Stokes_SLE };	

	Stokes_SLE* Stokes_SLE_New( 		
		Name                                                name,
		SLE_Solver*                                         solver,
		FiniteElementContext*                               context,
		Bool                                                isNonLinear,
		double                                              nonLinearTolerance,
		Iteration_Index                                     nonLinearMaxIterations,
		Bool                                                killNonConvergent,			
		EntryPoint_Register*                                entryPoint_Register,
		MPI_Comm                                            comm,
		StiffnessMatrix*                                    kStiffMat,
		StiffnessMatrix*                                    gStiffMat,
		StiffnessMatrix*                                    dStiffMat,
		StiffnessMatrix*                                    cStiffMat,
		SolutionVector*                                     uSolnVec,
		SolutionVector*                                     pSolnVec,
		ForceVector*                                        fForceVec,
		ForceVector*                                        hForceVec );

	/* Creation implementation / Virtual constructor */
	Stokes_SLE* _Stokes_SLE_New(
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

	void _Stokes_SLE_Init( 		
		void*                                               sle, 
		StiffnessMatrix*                                    kStiffMat,
		StiffnessMatrix*                                    gStiffMat,
		StiffnessMatrix*                                    dStiffMat,
		StiffnessMatrix*                                    cStiffMat,
		SolutionVector*                                     uSolnVec,
		SolutionVector*                                     pSolnVec,
		ForceVector*                                        fForceVec,
		ForceVector*                                        hForceVec ) ;

void Stokes_SLE_InitAll( 
		void*                                               sle, 
		SLE_Solver*                                         solver,
		FiniteElementContext*                               context,
		Bool                                                isNonLinear,
		double                                              nonLinearTolerance,
		Iteration_Index                                     nonLinearMaxIterations,
		Bool                                                killNonConvergent,			
		EntryPoint_Register*                                entryPoint_Register,
		MPI_Comm                                            comm,		
		StiffnessMatrix*                                    kStiffMat,
		StiffnessMatrix*                                    gStiffMat,
		StiffnessMatrix*                                    dStiffMat,
		StiffnessMatrix*                                    cStiffMat,
		SolutionVector*                                     uSolnVec,
		SolutionVector*                                     pSolnVec,
		ForceVector*                                        fForceVec,
		ForceVector*                                        hForceVec );

	void _Stokes_SLE_Delete( void* stokesSleSolver );
	void _Stokes_SLE_Print( void* stokesSleSolver, Stream* stream );

	void* _Stokes_SLE_DefaultNew( Name name ) ;
void _Stokes_SLE_AssignFromXML( void* sle, Stg_ComponentFactory* cf, void* data ) ;
	
	void _Stokes_SLE_MG_SelectStiffMats( void* _sle, unsigned* nSMs, StiffnessMatrix*** sms );
	
#endif
