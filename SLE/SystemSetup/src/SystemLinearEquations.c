/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street, Melbourne, 3053, Australia.
**
** Authors:
**	Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**	David May, PhD Student Monash University, VPAC. (davidm@vpac.org)
**	Patrick D. Sunter, Software Engineer, VPAC. (pds@vpac.org)
**	Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**	Siew-Ching Tan, Software Engineer, VPAC. (siew@vpac.org)
**	Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
**	Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
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
** $Id: SystemLinearEquations.c 1017 2008-02-01 00:36:32Z DavidLee $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include "StgFEM/Discretisation/Discretisation.h"
#include "StgFEM/SLE/LinearAlgebra/LinearAlgebra.h"
#include "units.h"
#include "types.h"
#include "shortcuts.h"
#include "SystemLinearEquations.h"
#include "SLE_Solver.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "StiffnessMatrix.h"
#include "SolutionVector.h"
#include "ForceVector.h"
#include "Context.h"


/* Textual name of this class */
const Type SystemLinearEquations_Type = "SystemLinearEquations";

/** Constructor */
SystemLinearEquations* SystemLinearEquations_New(
		Name                                               name,
		SLE_Solver*                                        solver,
		NonlinearSolver*				   nlSolver,
		FiniteElementContext*                              context,
		Bool                                               isNonLinear,
		double                                             nonLinearTolerance,
		Iteration_Index                                    nonLinearMaxIterations,
		Bool                                               killNonConvergent,		
		EntryPoint_Register*                               entryPoint_Register,
		MPI_Comm                                           comm )
{
	SystemLinearEquations* self = _SystemLinearEquations_DefaultNew( name );

	SystemLinearEquations_InitAll( 
			self, 
			solver,
			nlSolver,
			context,
			isNonLinear,
			nonLinearTolerance, 
			nonLinearMaxIterations,
			killNonConvergent, 
			entryPoint_Register,
			comm );

	return self;
}

/* Creation implementation / Virtual constructor */
SystemLinearEquations* _SystemLinearEquations_New( 
		SizeT                                              sizeOfSelf,
		Type                                               type,
		Stg_Class_DeleteFunction*                          _delete,
		Stg_Class_PrintFunction*                           _print,
		Stg_Class_CopyFunction*                            _copy, 
		Stg_Component_DefaultConstructorFunction*          _defaultConstructor,
		Stg_Component_ConstructFunction*                   _construct,
		Stg_Component_BuildFunction*                       _build,
		Stg_Component_InitialiseFunction*                  _initialise,
		Stg_Component_ExecuteFunction*                     _execute,
		Stg_Component_DestroyFunction*                     _destroy,
		SystemLinearEquations_LM_SetupFunction*            _LM_Setup,
		SystemLinearEquations_MatrixSetupFunction*         _matrixSetup,
		SystemLinearEquations_VectorSetupFunction*         _vectorSetup,
		SystemLinearEquations_UpdateSolutionOntoNodesFunc* _updateSolutionOntoNodes, 
		SystemLinearEquations_MG_SelectStiffMatsFunc*		_mgSelectStiffMats, 

		Name                                               name ) 
{
	SystemLinearEquations*					self;

	/* Allocate memory */
	assert( sizeOfSelf >= sizeof(SystemLinearEquations) );
	self = (SystemLinearEquations*) _Stg_Component_New( 
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
			name,
			NON_GLOBAL );
	
	/* Virtual info */
	self->_LM_Setup = _LM_Setup;
	self->_matrixSetup = _matrixSetup;
	self->_vectorSetup = _vectorSetup;
	self->_updateSolutionOntoNodes = _updateSolutionOntoNodes;
	self->_mgSelectStiffMats = _mgSelectStiffMats;
	
	return self;
}

void _SystemLinearEquations_Init( 
		void*                                              sle, 
		SLE_Solver*                                        solver, 
		NonlinearSolver*				   nlSolver,
		FiniteElementContext*                              context, 
		Bool                                               makeConvergenceFile,  
		Bool                                               isNonLinear,
		double                                             nonLinearTolerance,
		Iteration_Index                                    nonLinearMaxIterations,
		Bool                                               killNonConvergent,
		Iteration_Index                                    nonLinearMinIterations,
		Name						   nonLinearSolutionType,
		EntryPoint_Register*                               entryPoint_Register,
		MPI_Comm                                           comm ) 
{
	SystemLinearEquations* self = (SystemLinearEquations*)sle;
	char* filename;

	self->isConstructed = True;
	self->extensionManager = ExtensionManager_New_OfExistingObject( self->name, self );
	
	self->debug = Stream_RegisterChild( StgFEM_SLE_SystemSetup_Debug, self->type );
	self->info =  Journal_MyStream( Info_Type, self );
    /* Note: currently we're sending self->info to the master proc only so there's not too much
	   identical timing info printed. May want to fine-tune later so that some info does get 
       printed on all procs. */	
	Stream_SetPrintingRank( self->info, 0 );

	self->makeConvergenceFile = makeConvergenceFile;
	if ( self->makeConvergenceFile ) {
		self->convergenceStream = Journal_Register( InfoStream_Type, "Convergence Info" );
		Stg_asprintf( &filename, "Convergence.dat" );
		Stream_RedirectFile_WithPrependedPath( self->convergenceStream, context->outputPath, filename );
		Stream_SetPrintingRank( self->convergenceStream, 0 );
		Memory_Free( filename );
		Journal_Printf( self->convergenceStream , "Timestep\tIteration\tResidual\tTolerance\n" );
	}
	
	
	self->comm = comm;
	self->solver = solver;
	self->nlSolver = nlSolver;
	self->stiffnessMatrices = Stg_ObjectList_New();
	self->forceVectors = Stg_ObjectList_New();
	self->solutionVectors = Stg_ObjectList_New();
	self->bcRemoveQuery = True;

	/* Init NonLinear Stuff */
	if ( isNonLinear ) {
		self->nonLinearSolutionType	= nonLinearSolutionType;
		SystemLinearEquations_SetToNonLinear( self );
	}
	self->nonLinearTolerance        = nonLinearTolerance;
	self->nonLinearMaxIterations    = nonLinearMaxIterations;
	self->killNonConvergent         = killNonConvergent;
	self->nonLinearMinIterations    = nonLinearMinIterations;

	/* BEGIN LUKE'S FRICTIONAL BCS BIT */
	Stg_asprintf( &self->nlEPName, "%s-nlEP", self->name );
	self->nlEP = EntryPoint_New( self->nlEPName, EntryPoint_2VoidPtr_CastType );
	/* END LUKE'S FRICTIONAL BCS BIT */
	
	/* Initialise MG stuff. */
	self->mgEnabled = False;
	self->mgUpdate = True;
	self->nMGHandles = 0;
	self->mgHandles = NULL;
	
	/* Create Execute Entry Point */
	Stg_asprintf( &self->executeEPName, "%s-execute", self->name );
	self->executeEP = EntryPoint_New( self->executeEPName, EntryPoint_2VoidPtr_CastType );

	/* Add default hooks to Execute E.P. */
	EntryPoint_Append( self->executeEP, "BC_Setup", SystemLinearEquations_BC_Setup, self->type);
	EntryPoint_Append( self->executeEP, "LM_Setup", SystemLinearEquations_LM_Setup, self->type);
	EntryPoint_Append( self->executeEP, "IntegrationSetup", SystemLinearEquations_IntegrationSetup, self->type );
	EntryPoint_Append( self->executeEP, "ZeroAllVectors", SystemLinearEquations_ZeroAllVectors, self->type);
	EntryPoint_Append( self->executeEP, "MatrixSetup", SystemLinearEquations_MatrixSetup, self->type);
	EntryPoint_Append( self->executeEP, "VectorSetup", SystemLinearEquations_VectorSetup, self->type);
	EntryPoint_Append( self->executeEP, "ExecuteSolver", SystemLinearEquations_ExecuteSolver, self->type);
	EntryPoint_Append( self->executeEP, "UpdateSolutionOntoNodes",SystemLinearEquations_UpdateSolutionOntoNodes,self->type);

	/* Create Integration Setup EP */
	Stg_asprintf( &self->integrationSetupEPName, "%s-integrationSetup", self->name );
	self->integrationSetupEP = EntryPoint_New( self->integrationSetupEPName, EntryPoint_Class_VoidPtr_CastType );
	
	if ( entryPoint_Register )
		EntryPoint_Register_Add( entryPoint_Register, self->executeEP );
	self->entryPoint_Register = entryPoint_Register;

	/* Add SLE to Context */
	self->context = context;
	if ( context )
		FiniteElementContext_AddSLE( context, self );
}

void SystemLinearEquations_InitAll( 
		void*                                              sle, 
		SLE_Solver*                                        solver,
		NonlinearSolver*				   nlSolver,
		FiniteElementContext*                              context, 
		Bool                                               isNonLinear,
		double                                             nonLinearTolerance,
		Iteration_Index                                    nonLinearMaxIterations,
		Bool                                               killNonConvergent,
		EntryPoint_Register*                               entryPoint_Register,
		MPI_Comm                                           comm ) 
{
	/* TODO - Init Parent */
	_SystemLinearEquations_Init( 
			sle, 
			solver,
			nlSolver,
			context,
			False, /* TODO: A hack put in place for setting the convergence stream to 'off' if the SLE class is created from within the code, not via an xml */
			isNonLinear,
			nonLinearTolerance, 
			nonLinearMaxIterations,
			killNonConvergent, 
			1,/* TODO : hack for setting the minimum number of iterations to 1- same hack as above */
			"",
			entryPoint_Register,
			comm );
}


void _SystemLinearEquations_Delete( void* sle ) {
	SystemLinearEquations*					self = (SystemLinearEquations*)sle;
	
	Journal_DPrintf( self->debug, "In %s\n", __func__ );
	Stream_IndentBranch( StgFEM_Debug );
	
	/* BEGIN LUKE'S FRICTIONAL BCS BIT */
	Memory_Free( self->nlEPName );
	/* END LUKE'S FRICTIONAL BCS BIT */
	Memory_Free( self->executeEPName );
	
	Stg_Class_Delete( self->extensionManager );

	Stg_Class_Delete( self->stiffnessMatrices ); 
	Stg_Class_Delete( self->forceVectors ); 
	Stg_Class_Delete( self->solutionVectors ); 

/* 	 delete parent */
	_Stg_Component_Delete( self );
	Stream_UnIndentBranch( StgFEM_Debug );
}


void _SystemLinearEquations_Print( void* sle, Stream* stream ) {
	SystemLinearEquations*		self = (SystemLinearEquations*)sle;
	
	/* General info */
	Journal_Printf( stream, "SystemLinearEquations (ptr): %p\n", self );
	_Stg_Component_Print( self, stream );
	
	/* Virtual info */
	Stg_Class_Print( self->stiffnessMatrices, stream );
	Stg_Class_Print( self->forceVectors, stream );
	Stg_Class_Print( self->solutionVectors, stream );

	/* other info */
	Journal_PrintPointer( stream, self->extensionManager );
	Journal_Printf( stream, "\tcomm: %u\n", self->comm );
	Journal_Printf( stream, "\tsolver (ptr): %p\n", self->solver );
	Stg_Class_Print( self->solver, stream );
}


void* _SystemLinearEquations_Copy( void* sle, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	SystemLinearEquations*	self = (SystemLinearEquations*)sle;
	SystemLinearEquations*	newSLE;
	PtrMap*			map = ptrMap;
	Bool			ownMap = False;
	
	if( !map ) {
		map = PtrMap_New( 10 );
		ownMap = True;
	}
	
	newSLE = _Stg_Component_Copy( sle, dest, deep, nameExt, map );
	
	/* Virtual methods */
	newSLE->_LM_Setup = self->_LM_Setup;
	newSLE->_matrixSetup = self->_matrixSetup;
	newSLE->_vectorSetup = self->_vectorSetup;
	newSLE->_mgSelectStiffMats = self->_mgSelectStiffMats;
	
	newSLE->debug = Stream_RegisterChild( StgFEM_SLE_SystemSetup_Debug, newSLE->type );
	newSLE->comm = self->comm;
	
	if( deep ) {
		newSLE->solver = (SLE_Solver*)Stg_Class_Copy( self->solver, NULL, deep, nameExt, map );
		newSLE->stiffnessMatrices = (StiffnessMatrixList*)Stg_Class_Copy( self->stiffnessMatrices, NULL, deep, nameExt, map );
		newSLE->forceVectors = (ForceVectorList*)Stg_Class_Copy( self->forceVectors, NULL, deep, nameExt, map );
		newSLE->solutionVectors = (SolutionVectorList*)Stg_Class_Copy( self->solutionVectors, NULL, deep, nameExt, map );
		if( (newSLE->extensionManager = PtrMap_Find( map, self->extensionManager )) == NULL ) {
			newSLE->extensionManager = Stg_Class_Copy( self->extensionManager, NULL, deep, nameExt, map );
			PtrMap_Append( map, self->extensionManager, newSLE->extensionManager );
		}
	}
	else {
		newSLE->solver = self->solver;
		newSLE->stiffnessMatrices = self->stiffnessMatrices;
		newSLE->forceVectors = self->forceVectors;
		newSLE->solutionVectors = self->solutionVectors;
	}
	
	if( ownMap ) {
		Stg_Class_Delete( map );
	}
	
	return newSLE;
}

void* _SystemLinearEquations_DefaultNew( Name name ) {
	return _SystemLinearEquations_New(
		sizeof(SystemLinearEquations), 
		SystemLinearEquations_Type,
		_SystemLinearEquations_Delete,
		_SystemLinearEquations_Print,
		_SystemLinearEquations_Copy,
		_SystemLinearEquations_DefaultNew,
		_SystemLinearEquations_Construct,
		_SystemLinearEquations_Build,
		_SystemLinearEquations_Initialise,
		_SystemLinearEquations_Execute,
		_SystemLinearEquations_Destroy,
		_SystemLinearEquations_LM_Setup,
		_SystemLinearEquations_MatrixSetup,
		_SystemLinearEquations_VectorSetup,
		_SystemLinearEquations_UpdateSolutionOntoNodes,
		_SystemLinearEquations_MG_SelectStiffMats, 
		name );
}

void _SystemLinearEquations_Construct( void* sle, Stg_ComponentFactory* cf, void* data ){
	SystemLinearEquations*  self                     = (SystemLinearEquations*)sle;
	SLE_Solver*             solver                   = NULL;
	void*                   entryPointRegister       = NULL;
	FiniteElementContext*   context                  = NULL;
	double                  nonLinearTolerance;
	Iteration_Index         nonLinearMaxIterations;
	Bool                    isNonLinear;
	Bool                    killNonConvergent;
	Bool                    makeConvergenceFile;
	Iteration_Index         nonLinearMinIterations;                     
	Name			nonLinearSolutionType;
	NonlinearSolver*	nlSolver		= NULL;
	
	solver = Stg_ComponentFactory_ConstructByKey( cf, self->name, SLE_Solver_Type, SLE_Solver, False, data ) ;

	/* TODO - Construct Parent */

	makeConvergenceFile       = Stg_ComponentFactory_GetBool(   cf, self->name, "makeConvergenceFile",    False );
	isNonLinear               = Stg_ComponentFactory_GetBool(   cf, self->name, "isNonLinear",            False );
	nonLinearTolerance        = Stg_ComponentFactory_GetDouble( cf, self->name, "nonLinearTolerance",     0.01 );
	nonLinearMaxIterations    = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, "nonLinearMaxIterations", 500 );
	killNonConvergent         = Stg_ComponentFactory_GetBool(   cf, self->name, "killNonConvergent",      True );
	nonLinearMinIterations    = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, "nonLinearMinIterations", 1 );
	nonLinearSolutionType	  = Stg_ComponentFactory_GetString( cf, self->name, "nonLinearSolutionType", "" );
	
	entryPointRegister = Stg_ObjectList_Get( cf->registerRegister, "EntryPoint_Register" );
	assert( entryPointRegister );

	context = Stg_ComponentFactory_ConstructByKey( cf, self->name, "Context", FiniteElementContext, False, data );
								/* should this be a PETSc nls?? */
	//nlSolver = Stg_ComponentFactory_ConstructByKey( cf, self->name, NonlinearSolver_Type, NonlinearSolver, False, data );
	if( isNonLinear )
		nlSolver = PETScNonlinearSolver_New( "nonLinearSolver" );

	_SystemLinearEquations_Init( 
			self,
			solver,
			nlSolver,
			context,
			makeConvergenceFile,
			isNonLinear,
			nonLinearTolerance, 
			nonLinearMaxIterations,
			killNonConvergent, 
			nonLinearMinIterations,
			nonLinearSolutionType,
			entryPointRegister,
			MPI_COMM_WORLD );

	self->delta_x   = PETScVector_New( "delta x" );
	self->F		= PETScVector_New( "residual" );
	self->J		= PETScMatrix_New( "Jacobian" ); 
	//PETScMatrix_SetComm( self->J, MPI_COMM_WORLD );
	PETScMatrix_SetComm( self->J, self->comm );
}

/* Build */
void _SystemLinearEquations_Build( void* sle, void* _context ) {
	SystemLinearEquations*		self = (SystemLinearEquations*)sle;
	Index				index;
	
	Journal_DPrintf( self->debug, "In %s\n", __func__ );
	Stream_IndentBranch( StgFEM_Debug );

	/* build the matrices */
	for ( index = 0; index < self->stiffnessMatrices->count; index++ ) {
		/* Update rowSize and colSize if boundary conditions have been applied */		
		Stg_Component_Build( self->stiffnessMatrices->data[index], _context, False );
	}	
	
	/* and the vectors */
	for ( index = 0; index < self->forceVectors->count; index++ ) {
		/* Build the force vectors - includes updateing matrix size based on Dofs */
		Stg_Component_Build( self->forceVectors->data[index], _context, False );
	}
	
	/* and the solutions */
	for ( index = 0; index < self->solutionVectors->count; index++ ) {
		/* Build the force vectors - includes updateing matrix size based on Dofs */
		Stg_Component_Build( self->solutionVectors->data[index], _context, False );
	}
	
	/* lastly, the solver - if required */
	if( self->solver )
		Stg_Component_Build( self->solver, self, False );
	if( self->nlSolver )
		Stg_Component_Build( self->nlSolver, self, False );

	Stream_UnIndentBranch( StgFEM_Debug );
}


void _SystemLinearEquations_Initialise( void* sle, void* _context ) {
	SystemLinearEquations*					self = (SystemLinearEquations*)sle;
	Index						index;
	
	Journal_DPrintf( self->debug, "In %s\n", __func__ );
	Stream_IndentBranch( StgFEM_Debug );

	/* initialise the matrices */
	for ( index = 0; index < self->stiffnessMatrices->count; index++ ) {
		/* Update rowSize and colSize if boundary conditions have been applied */		
		Stg_Component_Initialise( self->stiffnessMatrices->data[index], _context, False );
	}	
	
	/* and the vectors */
	for ( index = 0; index < self->forceVectors->count; index++ ) {
		/* Initialise the force vectors - includes updateing matrix size based on Dofs */
		Stg_Component_Initialise( self->forceVectors->data[index], _context, False );
	}
	
	/* and the solutions */
	for ( index = 0; index < self->solutionVectors->count; index++ ) {
		/* Initialise the force vectors - includes updateing matrix size based on Dofs */
		Stg_Component_Initialise( self->solutionVectors->data[index], _context, False );
	}

	/* Check to see if any of the components need to make the SLE non-linear */
	SystemLinearEquations_CheckIfNonLinear( self );
	
	/* Setup Location Matrix */
	SystemLinearEquations_LM_Setup( self, _context );

	/* lastly, the solver, if required */
	if( self->solver )
		Stg_Component_Initialise( self->solver, self, False );
	if( self->nlSolver )
		Stg_Component_Initialise( self->nlSolver, self, False );
	Stream_UnIndentBranch( StgFEM_Debug );
}


void _SystemLinearEquations_Execute( void* sle, void* _context ) {
	SystemLinearEquations*	self = (SystemLinearEquations*)sle;
	
	Journal_DPrintf( self->debug, "In %s\n", __func__ );
	Stream_IndentBranch( StgFEM_Debug );

	_EntryPoint_Run_2VoidPtr( self->executeEP, sle, _context );
	
	Stream_UnIndentBranch( StgFEM_Debug );
}


void SystemLinearEquations_ExecuteSolver( void* sle, void* _context ) {
	SystemLinearEquations*	self = (SystemLinearEquations*)sle;
	double wallTime;
	/* Actually run the solver to get the new values into the SolutionVectors */
	
	Journal_Printf(self->info,"Linear solver (%s) \n",self->executeEPName);
		
	wallTime = MPI_Wtime();
	if( self->solver )	
		Stg_Component_Execute( self->solver, self, True );
	
	Journal_Printf(self->info,"Linear solver (%s), solution time %6.6e (secs)\n",self->executeEPName, MPI_Wtime() - wallTime);
		
}

void _SystemLinearEquations_Destroy( void* sle, void* _context ) {
	SystemLinearEquations*	self = (SystemLinearEquations*)sle;
	
	/* Free the the MG handles. */
	FreeArray( self->mgHandles );
}

void SystemLinearEquations_BC_Setup( void* sle, void* _context ) {
	SystemLinearEquations*				self = (SystemLinearEquations*)sle;
	Index						index;
	
	Journal_DPrintf( self->debug, "In %s\n", __func__ );
	for ( index = 0; index < self->solutionVectors->count; index++ ) {
		SolutionVector_ApplyBCsToVariables( self->solutionVectors->data[index], _context );
	}
}


void SystemLinearEquations_LM_Setup( void* sle, void* _context ) {
	SystemLinearEquations*					self = (SystemLinearEquations*)sle;
	
	self->_LM_Setup( self, _context );
}

void SystemLinearEquations_IntegrationSetup( void* sle, void* _context ) {
	SystemLinearEquations*					self = (SystemLinearEquations*)sle;
	
	_EntryPoint_Run_Class_VoidPtr( self->integrationSetupEP, _context );
}

void _SystemLinearEquations_LM_Setup( void* sle, void* _context ) {
	SystemLinearEquations*					self = (SystemLinearEquations*)sle;
	Index						index;

	Journal_DPrintf( self->debug, "In %s\n", __func__ );
	Stream_IndentBranch( StgFEM_Debug );
	/* For each feVariable of each stiffness matrix, build the LM  */
	for ( index = 0; index < self->stiffnessMatrices->count; index++ ) {
		StiffnessMatrix*				sm = (StiffnessMatrix*)self->stiffnessMatrices->data[index];
		
		FeEquationNumber_BuildLocationMatrix( sm->rowVariable->eqNum );
		FeEquationNumber_BuildLocationMatrix( sm->columnVariable->eqNum );
	}
	Stream_UnIndentBranch( StgFEM_Debug );
}

void SystemLinearEquations_MatrixSetup( void* sle, void* _context ) {
	SystemLinearEquations*					self = (SystemLinearEquations*)sle;
	
	self->_matrixSetup( self, _context );
}

void _SystemLinearEquations_MatrixSetup( void* sle, void* _context ) {
	SystemLinearEquations*				self = (SystemLinearEquations*)sle;
	FiniteElementContext*				context = (FiniteElementContext*)_context;
	Index						index;
	
	Journal_DPrintf( self->debug, "In %s\n", __func__ );
	Stream_IndentBranch( StgFEM_Debug );
	for ( index = 0; index < self->stiffnessMatrices->count; index++ ) {
		StiffnessMatrix_Assemble( self->stiffnessMatrices->data[index], self->bcRemoveQuery, self, context );
	}
	Stream_UnIndentBranch( StgFEM_Debug );
}


void SystemLinearEquations_VectorSetup( void* sle, void* _context ) {
	SystemLinearEquations*				self = (SystemLinearEquations*)sle;
	
	self->_vectorSetup( self, _context );
}

void _SystemLinearEquations_VectorSetup( void* sle, void* _context ) {
	SystemLinearEquations*				self = (SystemLinearEquations*)sle;
	Index						index;
	
	Journal_DPrintf( self->debug, "In %s\n", __func__ );
	Stream_IndentBranch( StgFEM_Debug );
	for ( index = 0; index < self->forceVectors->count; index++ ) {
		ForceVector_Assemble( self->forceVectors->data[index] );
	}
	Stream_UnIndentBranch( StgFEM_Debug );
}



Index _SystemLinearEquations_AddStiffnessMatrix( void* sle, StiffnessMatrix* stiffnessMatrix ) {
	SystemLinearEquations*					self = (SystemLinearEquations*)sle;
	
	return SystemLinearEquations_AddStiffnessMatrix( self, stiffnessMatrix );
}

StiffnessMatrix* _SystemLinearEquations_GetStiffnessMatrix( void* sle, Name stiffnessMatrixName ) {
	SystemLinearEquations*					self = (SystemLinearEquations*)sle;
	
	return SystemLinearEquations_GetStiffnessMatrix( self, stiffnessMatrixName );
}

Index _SystemLinearEquations_AddForceVector( void* sle, ForceVector* forceVector ) {
	SystemLinearEquations*					self = (SystemLinearEquations*)sle;
	
	return SystemLinearEquations_AddForceVector( self, forceVector );
}

ForceVector* _SystemLinearEquations_GetForceVector( void* sle, Name forceVectorName ) {
	SystemLinearEquations*					self = (SystemLinearEquations*)sle;
	
	return SystemLinearEquations_GetForceVector( self, forceVectorName );
}

Index _SystemLinearEquations_AddSolutionVector( void* sle, SolutionVector* solutionVector ) {
	SystemLinearEquations*					self = (SystemLinearEquations*)sle;
	
	return SystemLinearEquations_AddSolutionVector( self, solutionVector );
}

SolutionVector* _SystemLinearEquations_GetSolutionVector( void* sle, Name solutionVectorName ) {
	SystemLinearEquations*					self = (SystemLinearEquations*)sle;
	
	return SystemLinearEquations_GetSolutionVector( self, solutionVectorName );
}

void SystemLinearEquations_UpdateSolutionOntoNodes( void* sle, void* _context ) {
	SystemLinearEquations*	self = (SystemLinearEquations*)sle;

	self->_updateSolutionOntoNodes( self, _context );
}

void _SystemLinearEquations_UpdateSolutionOntoNodes( void* sle, void* _context ) {
	SystemLinearEquations*	self = (SystemLinearEquations*)sle;
	SolutionVector_Index	solnVec_I;
	SolutionVector*		currentSolnVec;

	for ( solnVec_I=0; solnVec_I < self->solutionVectors->count; solnVec_I++ ) {
		currentSolnVec = (SolutionVector*)self->solutionVectors->data[solnVec_I];
		SolutionVector_UpdateSolutionOntoNodes( currentSolnVec );
	}	
}

void SystemLinearEquations_ZeroAllVectors( void* sle, void* _context ) {
	SystemLinearEquations*      self = (SystemLinearEquations*)sle;
	Index                       index;
	ForceVector*                forceVector;
	
	Journal_DPrintf( self->debug, "In %s\n", __func__ );
	for ( index = 0; index < self->forceVectors->count; index++ ) {
		forceVector = (ForceVector*) self->forceVectors->data[index];
		
		Vector_Zero( forceVector->vector );
	}
}

void SystemLinearEquations_NewtonExecute( void* sle, void* _context ) {
	SystemLinearEquations*	self            = (SystemLinearEquations*) sle;
	//NonlinearSolver*	nlSolver	= self->nlSolver; //need to build this guy...
	PETScVector*		F;
	SNES			snes;

	Vector_Duplicate( SystemLinearEquations_GetSolutionVectorAt( self, 0 )->vector, &F );
	
	SNESCreate( MPI_COMM_WORLD, &snes );
	SNESSetJacobian( snes, ((PETScMatrix*)self->J)->petscMat, ((PETScMatrix*)self->J)->petscMat, self->_buildJ, self->buildJContext );
	SNESSetFunction( snes, ((PETScVector*)self->F)->petscVec, self->_buildF, self->buildFContext );
	
	SNESSetFromOptions( snes );

	SNESSolve( snes, PETSC_NULL, ((PETScVector*)self->delta_x)->petscVec );

	SNESDestroy( snes );

}

void SystemLinearEquations_NewtonMFFDExecute( void* sle, void* _context ) {
	SystemLinearEquations*	self            = (SystemLinearEquations*) sle;
	Vector*			F;
	NonlinearSolver*	nlSolver	= self->nlSolver; //need to build this guy...

	Vector_Duplicate( SystemLinearEquations_GetSolutionVectorAt( self, 0 )->vector, &F );

	/* creates the nonlinear solver */
	NonlinearSolver_SetComm( nlSolver, MPI_COMM_WORLD );
	NonlinearSolver_SetFunction( nlSolver, F, self->_buildF, _context );

	// set J (jacobian)
	
	// set F (residual vector)
	
	// call non linear solver func (SNES wrapper)
}

void SystemLinearEquations_NonLinearExecute( void* sle, void* _context ) {
	SystemLinearEquations*	self            = (SystemLinearEquations*) sle;
	Vector*                 previousVector;
	Vector*                 currentVector;
	double                  residual;
	double                  tolerance       = self->nonLinearTolerance;
	Iteration_Index         maxIterations   = self->nonLinearMaxIterations;
	Bool                    converged;
	Stream*                 errorStream     = Journal_Register( Error_Type, self->type );
	double					wallTime;
	Iteration_Index         minIterations   = self->nonLinearMinIterations;

	Journal_Printf( self->info, "In %s\n", __func__ );
	Stream_IndentBranch( StgFEM_Debug );
		
	wallTime = MPI_Wtime();		
		
	/* First Solve */
	self->nonLinearIteration_I = 0;
	Journal_Printf(self->info,"\nNon linear solver - iteration %d\n", self->nonLinearIteration_I);
	
	self->linearExecute( self, _context );
	self->hasExecuted = True;

	/* TODO - Give option which solution vector to test */
	currentVector   = SystemLinearEquations_GetSolutionVectorAt( self, 0 )->vector; 
	Vector_Duplicate( currentVector, (void**)&previousVector );
	Vector_SetLocalSize( previousVector, Vector_GetLocalSize( currentVector ) );
	
	for ( self->nonLinearIteration_I = 1 ; self->nonLinearIteration_I < maxIterations ; self->nonLinearIteration_I++ ) {
		/*
		** BEGIN LUKE'S FRICTIONAL BCS BIT
		**
		** Adding an interface for allowing other components to add some form of non-linearity to the system.
		** This is with a focus on frictional BCs, where we want to examine the stress field and modify
		** traction BCs to enforce friction rules. - Luke 18/07/2007
		*/

		_EntryPoint_Run_2VoidPtr( self->nlEP, sle, _context );

		/*
		** END LUKE'S FRICTIONAL BCS BIT
		*/


		Vector_CopyEntries( currentVector, previousVector );
	
		Journal_Printf(self->info,"Non linear solver - iteration %d\n", self->nonLinearIteration_I);
			
		self->linearExecute( self, _context );

		/* Calculate Residual */
		Vector_AddScaled( previousVector, -1.0, currentVector );
		residual = Vector_L2Norm( previousVector ) / Vector_L2Norm( currentVector );

		Journal_Printf( self->info, "In func %s: Iteration %u of %u - Residual %.5g - Tolerance = %.5g\n", 
				__func__, self->nonLinearIteration_I, maxIterations, residual, tolerance );
		if ( self->makeConvergenceFile ) {
			Journal_Printf( self->convergenceStream, "%d\t\t%d\t\t%.5g\t\t%.5g\n", 
							 self->context->timeStep, self->nonLinearIteration_I, residual, tolerance );
		}
			

		/* Check if residual is below tolerance */
		converged = (residual < tolerance);
		
		Journal_Printf(self->info,"Non linear solver - Residual %.8e; Tolerance %.4e%s%s - %6.6e (secs)\n\n", residual, tolerance, 
			(converged) ? " - Converged" : " - Not converged",
			(self->nonLinearIteration_I < maxIterations) ? "" : " - Reached iteration limit",
			MPI_Wtime() - wallTime );
		
		if ( (converged) && (self->nonLinearIteration_I>=minIterations) )
			break;
	}

	/* Print Info */
	if ( converged ) {
		Journal_Printf( self->info, "In func %s: Converged after %u iterations.\n",
				__func__, self->nonLinearIteration_I );
	}
	else {
		Journal_Printf( errorStream, "In func %s: Failed to converge after %u iterations.\n", 
				__func__, self->nonLinearIteration_I);
		if ( self->killNonConvergent ) {
			abort();
		}
	}

	Stream_UnIndentBranch( StgFEM_Debug );

	FreeObject( previousVector );
}

void SystemLinearEquations_AddNonLinearEP( void* sle, const char* name, EntryPoint_2VoidPtr_Cast func ) {
	SystemLinearEquations* self = (SystemLinearEquations*)sle;

	SystemLinearEquations_SetToNonLinear( self );
	EntryPoint_Append( self->nlEP, (char*)name, func, self->type );
}

void SystemLinearEquations_SetToNonLinear( void* sle ) {
	SystemLinearEquations*	self            = (SystemLinearEquations*) sle;

	assert( self );
	if ( self->isNonLinear )
		return;

	self->isNonLinear = True;

	self->linearExecute = self->_execute;
	self->_execute = SystemLinearEquations_NonLinearExecute;

	if( !strcmp( self->nonLinearSolutionType, "MatrixFreeNewton" ) )
		self->_execute = SystemLinearEquations_NewtonMFFDExecute;

	if( !strcmp( self->nonLinearSolutionType, "Newton" ) )
		self->_execute = SystemLinearEquations_NewtonExecute;
}

void SystemLinearEquations_CheckIfNonLinear( void* sle ) {
	SystemLinearEquations*	self            = (SystemLinearEquations*) sle;
	Index                   index;
	
	for ( index = 0; index < self->stiffnessMatrices->count; index++ ) {
		StiffnessMatrix* stiffnessMatrix = SystemLinearEquations_GetStiffnessMatrixAt( self, index );
		
		if ( stiffnessMatrix->isNonLinear )
			SystemLinearEquations_SetToNonLinear( self );

		/* TODO CHECK FOR FORCE VECTORS */
	}	
}

/*
** All the MG functions and their general implementations.
*/

void SystemLinearEquations_MG_Enable( void* _sle ) {
	SystemLinearEquations*	self = (SystemLinearEquations*)_sle;
	
	if( !self->isBuilt ) {
		Journal_Printf(self->info, "Warning: SLE has not been built, can't enable multi-grid.\n" );
		return;
	}
	
	self->mgEnabled = True;
}


void SystemLinearEquations_MG_SelectStiffMats( void* _sle, unsigned* nSMs, StiffnessMatrix*** sms ) {
	SystemLinearEquations*	self = (SystemLinearEquations*)_sle;
	
	assert( self->_mgSelectStiffMats );
	self->_mgSelectStiffMats( self, nSMs, sms );
}


void _SystemLinearEquations_MG_SelectStiffMats( void* _sle, unsigned* nSMs, StiffnessMatrix*** sms ) {
	SystemLinearEquations*	self = (SystemLinearEquations*)_sle;
	
	
	/*
	** As we have nothing else to go on, attempt to apply MG to all stiffness matrices in the list.
	*/
	
	{
		unsigned	sm_i;
		
		*nSMs = 0;
		for( sm_i = 0; sm_i < self->stiffnessMatrices->count; sm_i++ ) {
			StiffnessMatrix*	sm = ((StiffnessMatrix**)self->stiffnessMatrices->data)[sm_i];
			
			/* Add this one to the list. */
			*sms = Memory_Realloc_Array( *sms, StiffnessMatrix*, (*nSMs) + 1 );
			(*sms)[*nSMs] = sm;
			(*nSMs)++;
		}
	}
}
