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


#ifndef __StgFEM_SLE_ProvidedSystems_AdvectionDiffusion_AdvectionDiffusionSLE_h__
#define __StgFEM_SLE_ProvidedSystems_AdvectionDiffusion_AdvectionDiffusionSLE_h__

	extern const Type AdvectionDiffusionSLE_Type;
	
	#define __AdvectionDiffusionSLE \
		__SystemLinearEquations \
		/* Items passed into constructor */ \
		FeVariable*						phiField; \
		ForceVector*					residual; \
		Stg_Component*					massMatrix; \
		Dimension_Index				dim; \
		/* Items Created By solver */ \
		SolutionVector*				phiVector; /* should be passed in */ \
		double*							phiDotArray; \
		FeVariable*						phiDotField; \
		DofLayout*						phiDotDofLayout; \
		SolutionVector*				phiDotVector; \
		/* pointer to force term which must always be given */ \
		AdvDiffResidualForceTerm*	advDiffResidualForceTerm; \
		/* Timestep Stuff */ \
		double							currentDt; \
		double							courantFactor; \
		double							maxDiffusivity; \
		\
		Variable_Register*			variableReg; \
		FieldVariable_Register*		fieldVariableReg;
	
	struct AdvectionDiffusionSLE { __AdvectionDiffusionSLE };
		
	AdvectionDiffusionSLE* AdvectionDiffusionSLE_New( 
		Name							name,
		SLE_Solver*					solver,
		FiniteElementContext*	context,
		Bool							isNonLinear,
		double						nonLinearTolerance,
		Iteration_Index			nonLinearMaxIterations,
		Bool							killNonConvergent,		
		EntryPoint_Register*		entryPoint_Register,
		MPI_Comm						comm,
		FeVariable*					phiField,
		ForceVector*				residual,
		Stg_Component*				massMatrix,
		Dimension_Index			dim,
		double						courantFactor,
		Variable_Register*		variable_Register,
		FieldVariable_Register*	fieldVariable_Register ) ;
	
	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define ADVECTIONDIFFUSIONSLE_DEFARGS \
                SYSTEMLINEAREQUATIONS_DEFARGS

	#define ADVECTIONDIFFUSIONSLE_PASSARGS \
                SYSTEMLINEAREQUATIONS_PASSARGS

	AdvectionDiffusionSLE* _AdvectionDiffusionSLE_New(  ADVECTIONDIFFUSIONSLE_DEFARGS  );

	void _AdvectionDiffusionSLE_Init(
		void*							sle,
		FeVariable*					phiField,
		ForceVector*				residual,
		Stg_Component*				massMatrix,
		Dimension_Index			dim,
		double						courantFactor,
		Variable_Register*		variable_Register,  
		FieldVariable_Register*	fieldVariable_Register );

	/** Virtual Functions from "Class" Class */
	void _AdvectionDiffusionSLE_Delete( void* sle );

	void _AdvectionDiffusionSLE_Print( void* sle, Stream* stream );

	void* _AdvectionDiffusionSLE_Copy( void* sle, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap );

	/** Virtual Functions from "Stg_Component" Class */
	void* _AdvectionDiffusionSLE_DefaultNew( Name name );

	void _AdvectionDiffusionSLE_AssignFromXML( void* sle, Stg_ComponentFactory* cf, void* data );

	void _AdvectionDiffusionSLE_Build( void* sle, void* data );

	void _AdvectionDiffusionSLE_Initialise( void* sle, void* data );

	void _AdvectionDiffusionSLE_Execute( void* sle, void* _context );

	void _AdvectionDiffusionSLE_Destroy( void* sle, void* _context );

	//Vector* _AdvectionDiffusionSLE_GetResidual( void* sle, Index fv_I );
	Vec _AdvectionDiffusionSLE_GetResidual( void* sle, Index fv_I );

	void AdvectionDiffusionSLE_ResetStoredValues( void* sle );

#endif

