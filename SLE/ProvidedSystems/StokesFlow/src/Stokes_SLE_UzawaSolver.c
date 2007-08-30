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
** $Id: Stokes_SLE_UzawaSolver.c 948 2007-08-30 07:42:19Z DavidMay $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include "StgFEM/Discretisation/Discretisation.h"
#include "StgFEM/SLE/LinearAlgebra/LinearAlgebra.h"
#include "StgFEM/SLE/SystemSetup/SystemSetup.h"
#include "types.h"
#include "Stokes_SLE_UzawaSolver.h"

#include <assert.h>
#include <string.h>

#include "Stokes_SLE.h"

/* Macro to checking number integrity - i.e. checks if number is infinite or "not a number" */
#define isGoodNumber( number ) \
	( (! isnan( number ) ) && ( ! isinf( number ) ) )

const Type Stokes_SLE_UzawaSolver_Type = "Stokes_SLE_UzawaSolver";

void* _Stokes_SLE_UzawaSolver_DefaultNew( Name name ) {
	return (void*) _Stokes_SLE_UzawaSolver_New( 
		sizeof(Stokes_SLE_UzawaSolver), 
		Stokes_SLE_UzawaSolver_Type, 
		_Stokes_SLE_UzawaSolver_Delete, 
		_Stokes_SLE_UzawaSolver_Print, 
		_Stokes_SLE_UzawaSolver_Copy,
		_Stokes_SLE_UzawaSolver_DefaultNew,
		_Stokes_SLE_UzawaSolver_Construct,
		_Stokes_SLE_UzawaSolver_Build,
		_Stokes_SLE_UzawaSolver_Initialise,
		_SLE_Solver_Execute,
		_SLE_Solver_Destroy,
		_Stokes_SLE_UzawaSolver_SolverSetup, 
		_Stokes_SLE_UzawaSolver_Solve, 
		_Stokes_SLE_UzawaSolver_GetResidual,
		name );
}

Stokes_SLE_UzawaSolver* Stokes_SLE_UzawaSolver_New( 
		Name                                        name,
		Bool                                        useStatSolve, 
		int                                         statReps,
		StiffnessMatrix*                            preconditioner,
		Iteration_Index                             maxUzawaIterations,
		double                                      tolerance,
		Bool                                        useAbsoluteTolerance )
{		
	Stokes_SLE_UzawaSolver* self = _Stokes_SLE_UzawaSolver_DefaultNew( name );

	Stokes_SLE_UzawaSolver_InitAll( self, useStatSolve, statReps, preconditioner, maxUzawaIterations, tolerance, useAbsoluteTolerance );

	return self;
}


/* Creation implementation / Virtual constructor */
Stokes_SLE_UzawaSolver* _Stokes_SLE_UzawaSolver_New( 
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
		Name                                        name )
{
	Stokes_SLE_UzawaSolver* self;

	/* Allocate memory */
	assert( sizeOfSelf >= sizeof(Stokes_SLE_UzawaSolver) );
	self = (Stokes_SLE_UzawaSolver*) _SLE_Solver_New( 
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
		name );
	
	/* Virtual info */
	return self;
}


void _Stokes_SLE_UzawaSolver_Init( 
		Stokes_SLE_UzawaSolver*      self,
		StiffnessMatrix*             preconditioner, 
		Iteration_Index              maxUzawaIterations,
		double                       tolerance,
		Bool                         useAbsoluteTolerance )
{
	self->isConstructed        = True;
	self->tolerance            = tolerance;
	self->maxUzawaIterations   = maxUzawaIterations;
	self->preconditioner       = preconditioner;
	self->useAbsoluteTolerance = useAbsoluteTolerance;
}

void Stokes_SLE_UzawaSolver_InitAll( 
		void*                        solver,
		Bool                         useStatSolve,
		int                          statReps, 
		StiffnessMatrix*             preconditioner, 
		Iteration_Index              maxUzawaIterations,
		double                       tolerance,
		Bool                         useAbsoluteTolerance )
{
	Stokes_SLE_UzawaSolver* self = (Stokes_SLE_UzawaSolver*)solver;

	SLE_Solver_InitAll( self, useStatSolve, statReps );
	_Stokes_SLE_UzawaSolver_Init( self, preconditioner, maxUzawaIterations, tolerance, useAbsoluteTolerance );
}

void _Stokes_SLE_UzawaSolver_Delete( void* solver ) {
	Stokes_SLE_UzawaSolver* self = (Stokes_SLE_UzawaSolver*)solver;
		
	Journal_DPrintf( self->debug, "In: %s \n", __func__);

	Stream_IndentBranch( StgFEM_Debug );
	Journal_DPrintfL( self->debug, 2, "Deleting Solver contexts.\n" );
	FreeObject( self->velSolver );
	FreeObject( self->pcSolver );

	Journal_DPrintfL( self->debug, 2, "Deleting temporary solver vectors.\n" );
	FreeObject( self->pTempVec );
	FreeObject( self->rVec ); 
	FreeObject( self->sVec );
	FreeObject( self->fTempVec );
	FreeObject( self->vStarVec );
	Stream_UnIndentBranch( StgFEM_Debug );
}       


void _Stokes_SLE_UzawaSolver_Print( void* solver, Stream* stream ) {
	Stokes_SLE_UzawaSolver* self = (Stokes_SLE_UzawaSolver*)solver;

	_SLE_Solver_Print( self, stream );

	Journal_PrintValue( stream, self->tolerance );
	Journal_PrintValue( stream, self->maxUzawaIterations );
}


void* _Stokes_SLE_UzawaSolver_Copy( void* stokesSleUzawaSolver, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	Stokes_SLE_UzawaSolver* self = (Stokes_SLE_UzawaSolver*)stokesSleUzawaSolver;
	Stokes_SLE_UzawaSolver*	newStokesSleUzawaSolver;
	
	newStokesSleUzawaSolver = _SLE_Solver_Copy( self, dest, deep, nameExt, ptrMap );
	
	newStokesSleUzawaSolver->velSolver           = self->velSolver;
	newStokesSleUzawaSolver->pcSolver            = self->pcSolver;
	newStokesSleUzawaSolver->preconditioner      = self->preconditioner;
	newStokesSleUzawaSolver->pTempVec            = self->pTempVec;
	newStokesSleUzawaSolver->rVec                = self->rVec;
	newStokesSleUzawaSolver->sVec                = self->sVec;
	newStokesSleUzawaSolver->fTempVec            = self->fTempVec;
	newStokesSleUzawaSolver->vStarVec            = self->vStarVec;
	newStokesSleUzawaSolver->tolerance           = self->tolerance;
	newStokesSleUzawaSolver->maxUzawaIterations  = self->maxUzawaIterations;
	newStokesSleUzawaSolver->useAbsoluteTolerance  = self->useAbsoluteTolerance;
	
	return (void*) newStokesSleUzawaSolver;
}


void _Stokes_SLE_UzawaSolver_Build( void* solver, void* stokesSLE ) {
	Stokes_SLE_UzawaSolver*	self  = (Stokes_SLE_UzawaSolver*)solver;
	Stokes_SLE*             sle   = (Stokes_SLE*)stokesSLE;

 	Journal_DPrintf( self->debug, "In %s\n", __func__ );
	Stream_IndentBranch( StgFEM_Debug );
	
 	Journal_DPrintfL( self->debug, 2, "building a standard solver for the velocity system.\n" );
	MatrixSolver_SetComm( self->velSolver, sle->comm );
	
	/* Build Preconditioner */
	if ( self->preconditioner ) {
		Stg_Component_Build( self->preconditioner, stokesSLE, False );
		SystemLinearEquations_AddStiffnessMatrix( sle, self->preconditioner );

		Journal_DPrintfL( self->debug, 2, "build a standard solver for the preconditioner system.\n" );
		self->pcSolver = (MatrixSolver*)PETScMatrixSolver_New( "" );
		MatrixSolver_SetComm( self->pcSolver, sle->comm );
	}
	else 
		self->pcSolver = NULL;

 	Journal_DPrintfL( self->debug, 2, "Allocate the auxillary vectors pTemp, r, s, fTemp and vStar.\n" ); 
	Vector_Duplicate( sle->pSolnVec->vector, (void**)&self->pTempVec );
	Vector_SetLocalSize( self->pTempVec, Vector_GetLocalSize( sle->pSolnVec->vector ) );
	Vector_Duplicate( sle->pSolnVec->vector, (void**)&self->rVec );
	Vector_SetLocalSize( self->rVec, Vector_GetLocalSize( sle->pSolnVec->vector ) );
	Vector_Duplicate( sle->pSolnVec->vector, (void**)&self->sVec );
	Vector_SetLocalSize( self->sVec, Vector_GetLocalSize( sle->pSolnVec->vector ) );
	
	Vector_Duplicate( sle->fForceVec->vector, (void**)&self->fTempVec );
	Vector_SetLocalSize( self->fTempVec, Vector_GetLocalSize( sle->fForceVec->vector ) );
	Vector_Duplicate( sle->fForceVec->vector, (void**)&self->vStarVec );
	Vector_SetLocalSize( self->vStarVec, Vector_GetLocalSize( sle->fForceVec->vector ) );
	Stream_UnIndentBranch( StgFEM_Debug );
}

void _Stokes_SLE_UzawaSolver_Construct( void* solver, Stg_ComponentFactory* cf, void* data ) {
	Stokes_SLE_UzawaSolver* self         = (Stokes_SLE_UzawaSolver*) solver;
	double                  tolerance;
	Iteration_Index         maxUzawaIterations;
	StiffnessMatrix*        preconditioner;
	Bool                    useAbsoluteTolerance;

	_SLE_Solver_Construct( self, cf, data );

	tolerance            = Stg_ComponentFactory_GetDouble( cf, self->name, "tolerance", 1.0e-5 );
	maxUzawaIterations   = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, "maxIterations", 1000 );
	useAbsoluteTolerance = Stg_ComponentFactory_GetBool( cf, self->name, "useAbsoluteTolerance", False );

	preconditioner = Stg_ComponentFactory_ConstructByKey( cf, self->name, "Preconditioner", StiffnessMatrix, False, data );

	_Stokes_SLE_UzawaSolver_Init( self, preconditioner, maxUzawaIterations, tolerance, useAbsoluteTolerance );

	self->velSolver = Stg_ComponentFactory_ConstructByKey( cf, self->name, "velocitySolver", MatrixSolver, 
							       False, data );
	if( !self->velSolver ) {
#ifdef HAVE_PETSC
		self->velSolver = (MatrixSolver*)PETScMatrixSolver_New( "" );
#else
		self->velSolver = NULL;
#endif
	}
}

void _Stokes_SLE_UzawaSolver_Execute( void* solver, void* data ) {
}

void _Stokes_SLE_UzawaSolver_Destroy( void* solver, void* data ) {
}

void _Stokes_SLE_UzawaSolver_Initialise( void* solver, void* stokesSLE ) {
	Stokes_SLE_UzawaSolver* self = (Stokes_SLE_UzawaSolver*) solver;
	Stokes_SLE*             sle  = (Stokes_SLE*)             stokesSLE;
	
	/* Initialise Parent */
	_SLE_Solver_Initialise( self, sle );
	
	if ( sle->context && (True == sle->context->loadFromCheckPoint) ) {
		/* The previous timestep's velocity solution will be helpful in iterating to a better
		solution faster - and thus make restarting from checkpoint more repeatable compared
		to original non-restart solution */
		SolutionVector_LoadCurrentFeVariableValuesOntoVector( sle->uSolnVec );
		SolutionVector_LoadCurrentFeVariableValuesOntoVector( sle->pSolnVec );
	}

}

/* SolverSetup */


void _Stokes_SLE_UzawaSolver_SolverSetup( void* solver, void* stokesSLE ) {
	Stokes_SLE_UzawaSolver* self = (Stokes_SLE_UzawaSolver*) solver;
	Stokes_SLE*             sle  = (Stokes_SLE*)             stokesSLE;
	
 	Journal_DPrintf( self->debug, "In %s:\n", __func__ );
	Stream_IndentBranch( StgFEM_Debug );

	Journal_DPrintfL( self->debug, 1, "Setting up MatrixSolver for the velocity eqn.\n" );
	MatrixSolver_SetMatrix( self->velSolver, sle->kStiffMat->matrix );

	if( self->pcSolver ) {
		Journal_DPrintfL( self->debug, 1, "Setting up MatrixSolver for the Preconditioner.\n" );
		MatrixSolver_SetMatrix( self->pcSolver, self->preconditioner->matrix );
	}

	Stream_UnIndentBranch( StgFEM_Debug );
}


void _Stokes_SLE_UzawaSolver_Solve( void* solver, void* stokesSLE ) {
	Stokes_SLE_UzawaSolver* self            = (Stokes_SLE_UzawaSolver*)solver;	
	Stokes_SLE*             sle             = (Stokes_SLE*)stokesSLE;
	
	/* Create shortcuts to stuff needed on sle */
	Matrix*                 K_Mat           = sle->kStiffMat->matrix;
	Matrix*                 G_Mat           = sle->gStiffMat->matrix;
	Matrix*                 D_Mat           = NULL;
	Matrix*                 M_Mat           = NULL;
	Vector*                 uVec            = sle->uSolnVec->vector;
	Vector*                 qVec            = sle->pSolnVec->vector;
	Vector*                 fVec            = sle->fForceVec->vector;
	Vector*                 hVec            = sle->hForceVec->vector;
	
	/* Create shortcuts to solver related stuff */
	Vector*                 qTempVec        = self->pTempVec;  
	Vector*                 rVec            = self->rVec;
	Vector*                 sVec            = self->sVec;
	Vector*                 fTempVec        = self->fTempVec;
	Vector*                 vStarVec        = self->vStarVec; 
	
	MatrixSolver*           velSolver       = self->velSolver;	/*  Inner velocity solver */
	MatrixSolver*           pcSolver        = self->pcSolver;   /*  Preconditioner  */

	Iteration_Index         maxIterations   = self->maxUzawaIterations;	
	Iteration_Index         iteration_I     = 0;
	Iteration_Index         outputInterval  = 1;
	
	double                  zdotr_current	= 0.0;
	double                  zdotr_previous 	= 1.0;
	double                  sdotGTrans_v;
	double                  alpha, beta;
	double                  absResidual;  
	double                  relResidual;
	double*                 chosenResidual;	  /* We can opt to use either the absolute or relative residual in termination condition */
    	double                  uzawaRhsScale;      
	double                  divU;
	double                  weightedResidual;
	double                  weightedVelocityScale;
	double                  momentumEquationResidual;
	
	Iteration_Index         innerLoopIterations;
	Stream*                 errorStream     = Journal_Register( Error_Type, Stokes_SLE_UzawaSolver_Type );
	
	double                  qGlobalProblemScale = sqrt( (double)Vector_GetGlobalSize( qTempVec ) );
	double                  qReciprocalGlobalProblemScale = 1.0 / qGlobalProblemScale;
	int			init_info_stream_rank;	

	init_info_stream_rank = Stream_GetPrintingRank( self->info );
	Stream_SetPrintingRank( self->info, 0 ); 


	/*	DEFINITIONS:
					See accompanying documentation
					u - the displacement / velocity solution (to which constraints are applied)
					q - the pressure-like variable which constrains the divergence displacement / velocity	(= pressure for incompressible)	
					F - standard FE force vector
					Fhat - Uzawa RHS = K^{-1} G F  - h 
					K - standard FE stiffness matrix
					Khat - Uzawa transformed stiffness matrix = G^T K^{-1} G
					G matrix - discrete gradient operator
					D matrix - discrete divergence operator = G^T for this particular algorithm
					C matrix - Mass matrix (M) for compressibility 
					
		LM & DAM			
	*/

	/* CHOICE OF RESIDUAL: 
					we may opt to converge on the absolute value (self->useAbsoluteTolerance == True ... default)
					or the relative value of the residual (self->useAbsoluteTolerance == False) 
			 		(another possibility would be always to improve the residual by a given tolerance)
					The Moresi & Solomatov (Phys Fluids, 1995) approach is to use the relative tolerance	
	*/ 

	
	if ( Vector_L2Norm( fVec ) / sqrt( (double)Vector_GetGlobalSize( fVec ) ) <= 1e-99 ) {
		Journal_Printf( errorStream,
			"Error in func %s: The momentum force vector \"%s\" is zero. "
			"The force vector should be non-zero either because of your chosen boundary "
			"conditions, or because of the element force vector assembly. You have %d "
			"element force vectors attached.\n",
			__func__, sle->fForceVec->name, sle->fForceVec->assembleForceVector->hooks->count );
		if ( sle->fForceVec->assembleForceVector->hooks->count > 0 ) {
			Journal_Printf( errorStream, "You used the following force vector assembly terms:\n" );
			EntryPoint_PrintConcise( sle->fForceVec->assembleForceVector, errorStream );
/* 			 TODO : need to print the elementForceVector assembly, not the global guy!! */
		}	
		Journal_Printf( errorStream,
			"Please check values for building the force vector.\n" );
		Journal_Firewall( 0, errorStream, "Exiting.\n" ); 	
	}
	
					
 	Journal_DPrintf( self->debug, "In %s:\n", __func__ );
	Journal_PrintfL( self->debug, 2, "Conjugate Gradient Uzawa solver with:\n");
	
	Stream_IndentBranch( StgFEM_Debug );
	
	Journal_PrintfL( self->debug, 2, "Compressibility %s\n", (sle->cStiffMat)? "on" : "off");
	Journal_PrintfL( self->debug, 2, "Preconditioning %s\n", (pcSolver)? "on" : "off" );   
	
	
	
	if ( sle->cStiffMat ) {
		Journal_DPrintfL( self->debug, 2, "(compressibility active)\n" );
		M_Mat = sle->cStiffMat->matrix;   
	}
	else {
		Journal_DPrintfL( self->debug, 2, "(compressibility inactive)\n" );
	}
	if ( sle->dStiffMat ) {
		Journal_DPrintfL( self->debug, 2, "(asymmetric geometry: handling D Matrix [incorrectly - will be ignored])\n" );
		D_Mat = sle->dStiffMat->matrix;
	}
	else {
		Journal_DPrintfL( self->debug, 2, "(No D -> symmetric geometry: D = Gt)\n" );
	}
	
	#if DEBUG
	if ( Stream_IsPrintableLevel( self->debug, 3 ) ) {
		Journal_DPrintf( self->debug, "Matrices and Vectors to solve are:\n" );
		Journal_DPrintf( self->debug, "K Matrix:\n" );
		/* No nice way of viewing Matrices, so commented out as incompatible with
		 * new 3D decomp at present --Kathleen Humble 30-04-07 
		 * Matrix_View( sle->kStiffMat->matrix, self->debug ); */
		Journal_DPrintf( self->debug, "G Matrix:\n" );
		/*Matrix_View( G_Mat, self->debug ); */
		if ( D_Mat ) {
			Journal_DPrintf( self->debug, "D Matrix:\n" );
			/*Matrix_View( D_Mat, self->debug );*/
		}	
		if ( M_Mat ) {
			Journal_DPrintf( self->debug, "M Matrix:\n" );
			/*Matrix_View( M_Mat, self->debug );*/
		}	
		Journal_DPrintf( self->debug, "Z (preconditioner) Matrix:\n" );
		/*Matrix_View( self->preconditioner->matrix, self->debug ); */
		Journal_DPrintf( self->debug, "f Vector:\n" );
		Vector_View( fVec, self->debug );
		Journal_DPrintf( self->debug, "h Vector:\n" );
		Vector_View( hVec, self->debug );
	}
	#endif
	
	/* STEP 1: Estimate the magnitude of the RHS for the transformed problem
			   we compute (usually to lower accuracy than elsewhere) the RHS (Fhat - h) 
	         and store the result in qTempVec.
		LM & DAM			
	*/
	
	Journal_DPrintfL( self->debug, 2, "Building Fhat - h.\n" );
	
	MatrixSolver_SetRelativeTolerance( velSolver, self->tolerance );
	MatrixSolver_Solve( velSolver, fVec, vStarVec );
	innerLoopIterations = MatrixSolver_GetIterations( velSolver );
	
	Journal_DPrintfL( self->debug, 2, "Fhat inner solution: Number of iterations: %d\n", innerLoopIterations );
	
	Matrix_TransposeMultiply( G_Mat, vStarVec, qTempVec );
 	Vector_AddScaled( qTempVec, -1.0, hVec );
	
	/*  WARNING:
			If D != G^T then the resulting \hat{K} is not likely to be symmetric, positive definite as
			required by this implementation of the Uzawa iteration.  This next piece of code
			is VERY unlikely to work properly so it's in the sin bin for the time being - LM.
			
			if ( D_Mat ) {
				MatrixMultiply( D_Mat, vStarVec, qTempVec );
			}
			else {
				MatrixTransposeMultiply( G_Mat, vStarVec, qTempVec );
			}
		LM & DAM			
	*/	

	
	/* STEP 2: The problem scaling - optionally normalize the uzawa residual by the magnitude of the RHS (use a relative tolerance)
			For the inner velocity solver,  Citcom uses a relative tolerance equal to that used for the Uzawa iteration as a whole
		LM & DAM			
	*/
	
	if (self->useAbsoluteTolerance) {
		chosenResidual = &absResidual;
		Journal_PrintfL( self->info, 2, "Absolute residual < %g for Uzawa stopping condition\n", self->tolerance);
		/* We should calculate the effective relative tolerance and insert that here !! */
		MatrixSolver_SetRelativeTolerance( velSolver, 0.1 * self->tolerance );
	}
	else {  /* The CITCOM compatible choice */
		chosenResidual = &relResidual;
		Journal_PrintfL( self->info, 2, "Relative residual < %g for Uzawa stopping condition\n", self->tolerance);	
		MatrixSolver_SetRelativeTolerance( velSolver, 0.1 * self->tolerance );
	}
	
	Journal_DPrintfL( self->debug, 2, "Determining scaling factor for residual:\n" );
	uzawaRhsScale = Vector_L2Norm( qTempVec ) * qReciprocalGlobalProblemScale;
	
	Journal_DPrintfL( self->debug, 2, "uzawaRhsScale = %f\n", uzawaRhsScale );	
	Journal_Firewall( isGoodNumber( uzawaRhsScale ), errorStream, 
			"Error in func '%s' for %s '%s' - uzawaRhsScale has illegal value '%g'.\n", __func__, self->type, self->name, uzawaRhsScale );
	
	
	/* STEP 3: Calculate initial residual for transformed equation  (\hat{F} - h - \hat{K} q_0)
	    Compute the solution to K u_0 = F - G q_0  (u_0 unknown)
		  Then G^T u* = \hat{F} - \hat{K} q_0 
	    u_0 is also the initial velocity solution to which the constraint is applied by the subsequent iteration
		LM & DAM			
	*/
	
	Journal_DPrintfL( self->debug, 2, "Solving for transformed Uzawa RHS.\n" );
	
	Vector_CopyEntries( fVec, fTempVec );
	Vector_Scale( fTempVec, -1.0 );
	Matrix_MultiplyAdd( G_Mat, qVec, fTempVec, fTempVec );
	Vector_Scale( fTempVec, -1.0 );
	MatrixSolver_Solve( velSolver, fTempVec, uVec ); 
	
	/* Handling for NON-SYMMETRIC: relegated to sin bin (see comment above) 
	 if ( D_Mat ) {
    	MatrixMultiply( D_Mat, uVec, rVec );
	 }
	 else {
		 MatrixTransposeMultiply( G_Mat, uVec, rVec );
	 }	
		LM & DAM			
	*/
	
	Matrix_TransposeMultiply( G_Mat, uVec, rVec );
	divU = Vector_L2Norm( rVec ) / Vector_L2Norm( uVec );  
	
	Journal_PrintfL( self->info, 2, "Initial l2Norm( Div u ) / l2Norm( u ) = %f \n", divU);
	
	Journal_Firewall( isGoodNumber( divU ), errorStream, 
			"Error in func '%s' for %s '%s' - l2Norm( Div u ) has illegal value '%g'.\n",
			__func__, self->type, self->name, divU );
	
	
	Journal_DPrintfL( self->debug, 2, "Adding compressibility and prescribed divergence terms.\n" );
	
	if ( M_Mat ) {
		Matrix_MultiplyAdd( M_Mat, qVec, rVec, rVec );
	}	
	Vector_AddScaled( rVec, -1.0, hVec );
			
			
	/* STEP 4: Preconditioned conjugate gradient iteration loop */	
		
	Journal_DPrintfL( self->debug, 1, "Beginning main Uzawa conjugate gradient loop:\n" );	
	
	iteration_I = 0;
	do{	
		Journal_DPrintfL( self->debug, 2, "Beginning solve '%u'.\n", iteration_I );
		Stream_IndentBranch( StgFEM_Debug );
		
		/* STEP 4.1: Preconditioner
			Solve:
				Q_\hat{K} z_1 =  r_1
				Q_\hat{K} is an approximation to \hat{K} which is simple / trivial / quick to invert
			LM & DAM			
		*/
		
		if ( pcSolver ) 
			MatrixSolver_Solve( pcSolver, rVec, qTempVec );
		else
			Vector_CopyEntries( rVec, qTempVec );
				
		/* STEP 4.2: Calculate s_I, the pressure search direction
				z_{I-1} . r_{I-1}  
				\beta = (z_{I-1} . r_{I-1}) / (z_{I-2} . r_{I-2})  
					\beta = 0 for the first iteration
		      s_I = z_(I-1) + \beta * s_(I-1) 
			LM & DAM			
		*/ 
		
		zdotr_current = Vector_DotProduct( qTempVec, rVec );
		
		Journal_DPrintfL( self->debug, 2, "l2Norm (qTempVec) %g; (rVec) %g \n", 
			Vector_L2Norm( qTempVec ) * qReciprocalGlobalProblemScale, 
			Vector_L2Norm( rVec ) * qReciprocalGlobalProblemScale );
		
		if ( iteration_I == 0 ) {
			Vector_CopyEntries( qTempVec, sVec );
		}
		else {
			beta = zdotr_current/zdotr_previous;
			Vector_ScaleAdd( sVec, beta, qTempVec );
		}
		
		/* STEP 4.3: Velocity search direction corresponding to s_I is found by solving
				K u* = G s_I
			LM & DAM			
		*/
			
		Matrix_Multiply( G_Mat, sVec, fTempVec );
		
		Journal_DPrintfL( self->debug, 2, "Uzawa inner iteration step\n");
		MatrixSolver_Solve( velSolver, fTempVec, vStarVec );
		innerLoopIterations = MatrixSolver_GetIterations( velSolver );
		Journal_DPrintfL( self->debug, 2, "Completed Uzawa inner iteration in '%u' iterations \n", innerLoopIterations );
				
		
		/* STEP 4.4: Calculate the step size ( \alpha = z_{I-1} . r_{I-1} / (s_I . \hat{K} s_I) )
				 \hat{K} s_I = G^T u* - M s_I (u* from step 4.3) 	
			LM & DAM			
		*/ 
		
		Matrix_TransposeMultiply( G_Mat, vStarVec, qTempVec );
		
		/* Handling for NON-SYMMETRIC: relegated to sin bin (see comment above) 
		
			if ( D_Mat ) {
				MatrixMultiply( D_Mat, vStarVec, qTempVec );
			}
			else {
				MatrixTransposeMultiply( G_Mat, vStarVec, qTempVec );
			}
			LM & DAM			
		*/

		if ( M_Mat ) {
			Journal_DPrintfL( self->debug, 2, "Correcting for Compressibility\n" );
			Vector_Scale( qTempVec, -1.0 );
			Matrix_MultiplyAdd( M_Mat, sVec, qTempVec, qTempVec );
			Vector_Scale( qTempVec, -1.0 );
		}

		sdotGTrans_v = Vector_DotProduct( sVec, qTempVec );
		
		alpha = zdotr_current/sdotGTrans_v;
		
		
		/* STEP 4.5: Update pressure, velocity and value of residual
				 by \alpha times corresponding search direction 
			LM & DAM			
		*/
		
		Journal_DPrintfL( self->debug, 2, "zdotr_current = %g \n", zdotr_current);
		Journal_DPrintfL( self->debug, 2, "sdotGTrans_v = %g \n", sdotGTrans_v);
		Journal_DPrintfL( self->debug, 2, "alpha = %g \n", alpha);
	
		Journal_Firewall( 
				isGoodNumber( zdotr_current ) && isGoodNumber( sdotGTrans_v ) && isGoodNumber( alpha ), 
				errorStream, 
				"Error in func '%s' for %s '%s' - zdotr_current, sdotGTrans_v or alpha has an illegal value: '%g','%g' or '%g'\n",
				__func__, self->type, self->name, zdotr_current, sdotGTrans_v, alpha );
		
		Vector_AddScaled( qVec, alpha, sVec );
		Vector_AddScaled( uVec, -alpha, vStarVec );
		Vector_AddScaled( rVec, -alpha, qTempVec );
		
		/* STEP 4.6: store the value of z_{I-1} . r_{I-1} for the next iteration
		 LM & DAM
		*/
		
		zdotr_previous = zdotr_current; 
		
		absResidual = Vector_L2Norm( rVec ) * qReciprocalGlobalProblemScale;  
		relResidual = absResidual / uzawaRhsScale;
		
		Stream_UnIndentBranch( StgFEM_Debug );
		
		if( iteration_I % outputInterval == 0 ) {
			Journal_PrintfL( self->info, 2, "\tLoop = %u, absResidual = %.8e, relResidual = %.8e\n", 
				iteration_I, absResidual, relResidual );
		}
		
		Journal_Firewall( isGoodNumber( absResidual ), errorStream, 
				"Error in func '%s' for %s '%s' - absResidual has an illegal value: '%g'\n",
				__func__, self->type, self->name, absResidual );
		
		Journal_Firewall( iteration_I < maxIterations, 
				errorStream, "In func %s: Reached maximum number of iterations %u without converging; absResidual = %.5g, relResidual = %.5g \n",
				__func__, iteration_I, absResidual, relResidual );

/* 		 TODO: test for small change in 10 iterations and if so restart? */
			
	iteration_I++;  
	}  while (*chosenResidual > self->tolerance );  

	Journal_DPrintfL( self->debug, 1, "Pressure solution converged. Exiting uzawa \n ");
	
	/* STEP 5:  Check all the relevant residuals and report back */
	
	
	if (Stream_IsEnable( self->info ) ) {
	
	/* This information should be in an info stream */
	Journal_PrintfL( self->info, 1, "Summary: Uzawa_its = %04d; Uzawa Residual  = %.8e\n", iteration_I, relResidual );  
	Matrix_TransposeMultiply( G_Mat, uVec, rVec );	
	divU = Vector_L2Norm( rVec ) / Vector_L2Norm( uVec );  
	Journal_PrintfL( self->info, 1, "Summary: || Div ( u ) || / || u ||         = %.8e\n", divU);
	
	/* Residual for the momentum equation 
		Compute r = || F - Ku - Gp || / || F ||
	*/
	
	Matrix_Multiply( G_Mat, qVec, vStarVec );	
	Matrix_MultiplyAdd( K_Mat, uVec, vStarVec, fTempVec );
	Vector_ScaleAdd( fTempVec, -1.0, fVec );
	
	momentumEquationResidual = Vector_L2Norm( fTempVec ) / Vector_L2Norm( fVec ); 
	Journal_PrintfL( self->info, 1, "Summary: || F - Ku - Gp || / || F ||       = %.8e\n", momentumEquationResidual );
	Journal_Firewall( isGoodNumber( momentumEquationResidual ), errorStream, 
			"Bad residual for the momentum equation (|| F - Ku - Gp || / || F || = %g):\n"
			"\tCheck to see if forcing term is zero or nan - \n\t|| F - Ku - Gp || = %g \n\t|| F || = %g.\n", 
			momentumEquationResidual,
			Vector_L2Norm( fTempVec ),
			Vector_L2Norm( fVec ) );
		
	/* "Preconditioned"	residual for the momentum equation 
	 		r_{w} = || Q_{K}(r) || / || Q_{K}(F)
			fTempVec contains the residual but is overwritten once used
			vStarVec is used to hold the diagonal preconditioner Q_{K} 
	*/
			
	Matrix_GetDiagonal( K_Mat, vStarVec );
	Vector_Reciprocal( vStarVec );	

	Vector_PointwiseMultiply( vStarVec, fTempVec, fTempVec );	
	weightedResidual = Vector_L2Norm( fTempVec );
	
	Vector_PointwiseMultiply( vStarVec, fVec, fTempVec );	
	weightedVelocityScale = Vector_L2Norm( fTempVec );
		
	Journal_PrintfL( self->info, 1, "Summary: || F - Ku - Gp ||_w / || F ||_w   = %.8e\n", 
		weightedResidual / weightedVelocityScale );
		
		
	/* Report back on the solution - velocity and pressure 
	 Note - correction for dof in Vrms ??
	*/
	
	Journal_PrintfL( self->info, 1, "Summary: Max velocity component = %.8e; Velocity magnitude = %.8e\n",
		Vector_LInfNorm( uVec ),
		Vector_L2Norm( uVec ) / sqrt( (double)Vector_GetGlobalSize( uVec ) )  );
		
	
	Journal_PrintfL( self->info, 1, "Summary: Max pressure value     = %.8e; Pressure magnitude = %.8e\n",
		Vector_LInfNorm( qVec ),
		Vector_L2Norm( qVec ) / sqrt( (double)Vector_GetGlobalSize( qVec ) ) );
		
	}	


	#if DEBUG
	if ( Stream_IsPrintableLevel( self->debug, 3 ) ) {
		Journal_DPrintf( self->debug, "Velocity solution:\n" );
		Vector_View( uVec, self->debug );
		Journal_DPrintf( self->debug, "Pressure solution:\n" );
		Vector_View( qVec, self->debug );
	}
	#endif
	Stream_UnIndentBranch( StgFEM_Debug );

        Stream_SetPrintingRank( self->info, init_info_stream_rank );	
}


Vector* _Stokes_SLE_UzawaSolver_GetResidual( void* solver, Index fv_I ) {
	return NULL;
}

