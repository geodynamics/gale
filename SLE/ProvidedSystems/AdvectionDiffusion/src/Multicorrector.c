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
** $Id: Multicorrector.c 656 2006-10-18 06:45:50Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgFEM/Discretisation/Discretisation.h>
#include <StgFEM/SLE/LinearAlgebra/LinearAlgebra.h>
#include <StgFEM/SLE/SystemSetup/SystemSetup.h>

#include "types.h"
#include "AdvectionDiffusionSLE.h"
#include "Multicorrector.h"
#include "Residual.h"
#include "MassMatrix_Assembly.h"

#include <assert.h>

/* Textual name of this class */
const Type AdvDiffMulticorrector_Type = "AdvDiffMulticorrector";

AdvDiffMulticorrector* AdvDiffMulticorrector_New( 
		Name                                                name,
		double                                              gamma,
		Iteration_Index                                     multiCorrectorIterations )
{
	AdvDiffMulticorrector* self = (AdvDiffMulticorrector*) _AdvDiffMulticorrector_DefaultNew( name );

	AdvDiffMulticorrector_InitAll( self, gamma, multiCorrectorIterations );
	return self;
}

/* Creation implementation / Virtual constructor */
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
		SLE_Solver_MG_SetupSmootherFunc*                    _mgSetupSmoother,
		Name                                                name )
{
	AdvDiffMulticorrector* self;
	
	/* Allocate memory */
	assert( sizeOfSelf >= sizeof(AdvDiffMulticorrector) );
	self = (AdvDiffMulticorrector*) _SLE_Solver_New( 
		sizeOfSelf, 
		type, 
		_delete, 
		_print, 
		_copy,
		_defaultConstructor,
		_construct,
		_build, 
		_initialise,
		_execute,
		_destroy,
		_solverSetup,
		_solve,
		_getResidual, 
		_mgSetupSmoother, 
		name );
	
	/* Virtual info */
	
	return self;
}

void _AdvDiffMulticorrector_Init( 
		AdvDiffMulticorrector*                              self, 
		double                                              gamma,
		Iteration_Index                                     multiCorrectorIterations )
{
	self->gamma                    = gamma;
	self->multiCorrectorIterations = multiCorrectorIterations;
}

void AdvDiffMulticorrector_InitAll( 
		void*                                               solver,
		double                                              gamma,
		Iteration_Index                                     multiCorrectorIterations )
{
	AdvDiffMulticorrector* self = (AdvDiffMulticorrector*) solver;

	SLE_Solver_InitAll( self, False, 0 );
	_AdvDiffMulticorrector_Init( self, gamma, multiCorrectorIterations );
}

void _AdvDiffMulticorrector_Delete( void* solver ) {
	AdvDiffMulticorrector* self = (AdvDiffMulticorrector*)solver;

	if ( self->matrixSolver )
		MatrixSolver_Destroy( self->matrixSolver );

	_SLE_Solver_Delete( self );
}

void _AdvDiffMulticorrector_Print( void* solver, Stream* stream ) {
	AdvDiffMulticorrector* self = (AdvDiffMulticorrector*)solver;
	
	_SLE_Solver_Print( self, stream );

	Journal_PrintValue( stream, self->gamma );
	Journal_PrintValue( stream, self->multiCorrectorIterations );
}

void* _AdvDiffMulticorrector_DefaultNew( Name name ) {
	return (void*)_AdvDiffMulticorrector_New( 
		sizeof(AdvDiffMulticorrector), 
		AdvDiffMulticorrector_Type,
		_AdvDiffMulticorrector_Delete,
		_AdvDiffMulticorrector_Print,
		NULL,
		_AdvDiffMulticorrector_DefaultNew,
		_AdvDiffMulticorrector_Construct,
		_AdvDiffMulticorrector_Build,
		_AdvDiffMulticorrector_Initialise,
		_AdvDiffMulticorrector_Execute,
		_AdvDiffMulticorrector_Destroy,
		_AdvDiffMulticorrector_SolverSetup,
		_AdvDiffMulticorrector_Solve,
		NULL, /*_AdvDiffMulticorrector_GetResidual, */
		_SLE_Solver_MG_SetupSmoother,
		name );
}

void _AdvDiffMulticorrector_Construct( void* solver, Stg_ComponentFactory* cf, void* data ) {
	AdvDiffMulticorrector*                     self             = (AdvDiffMulticorrector*)solver;
	double                                     gamma;
	Iteration_Index                            multiCorrectorIterations;

	/* Construct Parent */
	_SLE_Solver_Construct( self, cf, data );

	gamma = Stg_ComponentFactory_GetDouble( cf, self->name, "gamma", 0.5 );
	multiCorrectorIterations = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, "multiCorrectorIterations", 2 );

	_AdvDiffMulticorrector_Init( self, gamma, multiCorrectorIterations );
}

void _AdvDiffMulticorrector_Build( void* solver, void* data ) {
	AdvDiffMulticorrector* self   = Stg_CheckType( solver, AdvDiffMulticorrector );

	_SLE_Solver_Build( self, data );
}

void _AdvDiffMulticorrector_Initialise( void* solver, void* data ) {
	AdvDiffMulticorrector*             self             = (AdvDiffMulticorrector*)solver;

	_SLE_Solver_Initialise( self, data );
}

void _AdvDiffMulticorrector_Execute( void* solver, void* data ) {
	_SLE_Solver_Execute( solver, data );
}

void _AdvDiffMulticorrector_Destroy( void* solver, void* data ) {
	_SLE_Solver_Destroy( solver, data );
}

void _AdvDiffMulticorrector_SolverSetup( void* solver, void* data ) {
	AdvDiffMulticorrector* self   = Stg_CheckType( solver, AdvDiffMulticorrector );
	AdvectionDiffusionSLE* sle    = Stg_CheckType( data, AdvectionDiffusionSLE );
	
	if ( self->matrixSolver && Stg_Class_IsInstance( sle->massMatrix, StiffnessMatrix_Type ) ) {
		StiffnessMatrix* massMatrix = Stg_CheckType( sle->massMatrix, StiffnessMatrix );
		MatrixSolver_Setup( self->matrixSolver, massMatrix->matrix );
	}
}

/* See Brooks, Hughes 1982 Section 4.2 
 * All equations refer to this paper if not otherwise indicated */
void _AdvDiffMulticorrector_Solve( void* solver, void* _sle ) {
	AdvDiffMulticorrector* self   = (AdvDiffMulticorrector*) solver;
	AdvectionDiffusionSLE* sle    = (AdvectionDiffusionSLE*) _sle;
	double                 dt     = sle->currentDt;
	Index                  iteration_I;
	Vector*                deltaPhiDot;

	Journal_DPrintf( sle->debug, "In func %s:\n", __func__ );

	/* Put mesh data onto vectors */
	SolutionVector_LoadCurrentFeVariableValuesOntoVector( sle->phiVector );
	SolutionVector_LoadCurrentFeVariableValuesOntoVector( sle->phiDotVector );

	/* Solve for predictor step */
	AdvDiffMulticorrector_Predictors( self, sle, dt );

	/* Allocate Memory For Corrector Step */
	Vector_Duplicate( sle->phiVector->vector, &deltaPhiDot );

	/* Multi-corrector Steps */
	for ( iteration_I = 0 ; iteration_I < self->multiCorrectorIterations ; iteration_I++ ) {
		AdvDiffMulticorrector_Solution( self, sle, deltaPhiDot );
		AdvDiffMulticorrector_Correctors( self, sle, deltaPhiDot, dt );

		/* Put solutions onto meshes */
		SolutionVector_UpdateSolutionOntoNodes( sle->phiVector );
		SolutionVector_UpdateSolutionOntoNodes( sle->phiDotVector );

		SystemLinearEquations_ZeroAllVectors( sle, NULL );
	}

	/* Clean Up */
	Vector_Destroy( deltaPhiDot );
}

/** See Eqns. 4.2.3-4 */
void AdvDiffMulticorrector_Predictors( AdvDiffMulticorrector* self, AdvectionDiffusionSLE* sle, double dt ) {
	double factor       = dt * ( 1.0 - self->gamma );
	Stream* debugStream = sle->debug;

	Journal_DPrintf( debugStream, "In func %s:\n", __func__ );

	#if DEBUG
	if ( Stream_IsPrintableLevel( debugStream, 3 ) ) {
		Journal_DPrintf( debugStream, "At start of %s:\n", __func__ );
		Stream_Indent( debugStream );

		Journal_PrintValue( debugStream, dt );
		Journal_PrintValue( debugStream, self->gamma );
		Journal_PrintValue( debugStream, factor );

		Journal_DPrintf( debugStream, "Phi:\n" );
		Vector_View( sle->phiVector->vector, debugStream );
		Journal_DPrintf( debugStream, "Phi Dot:\n" );
		Vector_View( sle->phiDotVector->vector, debugStream );

		Stream_UnIndent( debugStream );
	}
	#endif

	/* Calculate Predictor for \phi - 
	 * Eq. 4.2.3: \phi_{n+1}^{(0)} = \phi_n + \Delta t(1 - \gamma)\dot \phi_n */
	Vector_AddScaledVector( sle->phiVector->vector, factor, sle->phiDotVector->vector ); 
	
	/* Calculate Predictor for \dot \phi - 
	 * Eq. 4.2.4: \dot \phi_{n+1}^{(0)} = 0 */
	Vector_Zero( sle->phiDotVector->vector );

	#if DEBUG
	if ( Stream_IsPrintableLevel( debugStream, 3 ) ) {
		Journal_DPrintf( debugStream, "At end of %s: Phi is:\n", __func__ );
		Vector_View( sle->phiVector->vector, debugStream );
	}
	#endif
}
	
	
void AdvDiffMulticorrector_Solution( AdvDiffMulticorrector* self, AdvectionDiffusionSLE* sle, Vector* deltaPhiDot ) {
	Journal_DPrintf( sle->debug, "In func %s:\n", __func__ );

	/* Calculate Residual - See Eq. 4.2.6 */
	SystemLinearEquations_VectorSetup( sle, NULL );
	SystemLinearEquations_MatrixSetup( sle, NULL );

	/* Calculate Mass Matrix out of three options - fully explicit, fully implicit, split operators */
	AdvDiffMulticorrector_CalculatePhiDot( self, sle, deltaPhiDot );

	#if DEBUG
	if ( Stream_IsPrintableLevel( self->debug, 3 ) ) {
		Journal_DPrintf( self->debug, "Delta Phi Dot is:\n" );
		Vector_View( deltaPhiDot, self->debug );
	}
	#endif
}


	
/* Correct \phi and \dot \phi - See Eqns. 4.2.7-8 */
void AdvDiffMulticorrector_Correctors( AdvDiffMulticorrector* self, AdvectionDiffusionSLE* sle, Vector* deltaPhiDot, double dt ) {
	double factor = dt * self->gamma;
	
	Journal_DPrintf( sle->debug, "In func %s:\n", __func__ );

	/* Add correction to \phi - Eq. 4.2.7 */
	Vector_AddScaledVector( sle->phiVector->vector, factor, deltaPhiDot );
	
	/* Add correction to \dot \phi - Eq. 4.2.8 */
	Vector_AddScaledVector( sle->phiDotVector->vector, 1.0, deltaPhiDot );
}


void AdvDiffMulticorrector_CalculatePhiDot( AdvDiffMulticorrector* self, AdvectionDiffusionSLE* sle, Vector* deltaPhiDot ) {
	Stg_Component* massMatrix = sle->massMatrix;

	if ( Stg_Class_IsInstance( massMatrix, ForceVector_Type ) ) 
		_AdvDiffMulticorrector_CalculatePhiDot_Explicit( self, sle, deltaPhiDot );
	else if ( Stg_Class_IsInstance( massMatrix, StiffnessMatrix_Type ) )
		_AdvDiffMulticorrector_CalculatePhiDot_Implicit( self, sle, deltaPhiDot );
	else {
		Journal_Firewall( False, Journal_Register( Error_Type, self->type ),
				"Error in func '%s': Cannot understand type '%s' for mass matrix '%s'.\n",
				__func__, massMatrix->name, massMatrix->type );
	}
}

/* Lump all things onto diagonal of matrix - which is stored as a vector - Eq. 4.2.11 */
void _AdvDiffMulticorrector_CalculatePhiDot_Explicit( AdvDiffMulticorrector* self, AdvectionDiffusionSLE* sle, Vector* deltaPhiDot ) {
	ForceVector* massMatrix = Stg_CheckType( sle->massMatrix, ForceVector );

	/* Calculate change in \dot \phi - See Eq. 4.2.5 */
	Vector_PointwiseDivide( sle->residual->vector, massMatrix->vector, deltaPhiDot );
}

void _AdvDiffMulticorrector_CalculatePhiDot_Implicit( AdvDiffMulticorrector* self, AdvectionDiffusionSLE* sle, Vector* deltaPhiDot ) {

	MatrixSolver_Solve( self->matrixSolver, deltaPhiDot, sle->residual->vector );
}
