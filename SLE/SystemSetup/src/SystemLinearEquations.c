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
** $Id: SystemLinearEquations.c 1167 2008-06-30 09:43:05Z LukeHodkinson $
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

	self->_sleFormFunction = NULL;
	
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
		Name						   optionsPrefix,
		EntryPoint_Register*                               entryPoint_Register,
		MPI_Comm                                           comm ) 
{
	SystemLinearEquations*  self 		= (SystemLinearEquations*)sle;
	char* 			filename;
	char* 			optionsName;

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
	self->context = context;

	/* Init NonLinear Stuff */
	self->nonLinearSolutionType = nonLinearSolutionType; /* This will never got propogated through to _Initialise->SetToNonLinear if we keep it in the loop */
	if ( isNonLinear ) {
		self->nonLinearSolutionType	= nonLinearSolutionType;
		SystemLinearEquations_SetToNonLinear( self );
	}
	self->nonLinearTolerance        = nonLinearTolerance;
	self->nonLinearMaxIterations    = nonLinearMaxIterations;
	self->killNonConvergent         = killNonConvergent;
	self->nonLinearMinIterations    = nonLinearMinIterations;    
									     /* _  /0 */
	optionsName = Memory_Alloc_Array_Unnamed( char, strlen(optionsPrefix) + 1 + 1 );
	sprintf( optionsName, "%s_", optionsPrefix );
	self->optionsPrefix = optionsName;

	/* BEGIN LUKE'S FRICTIONAL BCS BIT */
	Stg_asprintf( &self->nlSetupEPName, "%s-nlSetupEP", self->name );
	self->nlSetupEP = EntryPoint_New( self->nlSetupEPName, EntryPoint_2VoidPtr_CastType );
	Stg_asprintf( &self->nlEPName, "%s-nlEP", self->name );
	self->nlEP = EntryPoint_New( self->nlEPName, EntryPoint_2VoidPtr_CastType );
	Stg_asprintf( &self->nlEPName, "%s-postNlEP", self->name );
	self->postNlEP = EntryPoint_New( self->postNlEPName, EntryPoint_2VoidPtr_CastType );
	Stg_asprintf( &self->nlConvergedEPName, "%s-nlConvergedEP", self->name );
	self->nlConvergedEP = EntryPoint_New( self->nlConvergedEPName, EntryPoint_2VoidPtr_CastType );
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
			"",
			entryPoint_Register,
			comm );
}


void _SystemLinearEquations_Delete( void* sle ) {
	SystemLinearEquations*					self = (SystemLinearEquations*)sle;
	
	Journal_DPrintf( self->debug, "In %s\n", __func__ );
	Stream_IndentBranch( StgFEM_Debug );
	
	/* BEGIN LUKE'S FRICTIONAL BCS BIT */
        Memory_Free( self->nlSetupEPName );
	Memory_Free( self->nlEPName );
	Memory_Free( self->postNlEPName );
        Memory_Free( self->nlConvergedEPName );
	/* END LUKE'S FRICTIONAL BCS BIT */
	Memory_Free( self->executeEPName );
	
	Stg_Class_Delete( self->extensionManager );

	Stg_Class_Delete( self->stiffnessMatrices ); 
	Stg_Class_Delete( self->forceVectors ); 
	Stg_Class_Delete( self->solutionVectors ); 

/* 	 delete parent */
	_Stg_Component_Delete( self );
	Stream_UnIndentBranch( StgFEM_Debug );

	Memory_Free( self->optionsPrefix );
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
	Name			optionsPrefix;
	double double_from_dict;
	Bool   bool_from_dict;
	int    int_from_dict;
	
	solver = Stg_ComponentFactory_ConstructByKey( cf, self->name, SLE_Solver_Type, SLE_Solver, False, data ) ;

	/* TODO - Construct Parent */

	makeConvergenceFile       = Stg_ComponentFactory_GetBool(   cf, self->name, "makeConvergenceFile",    False );
	isNonLinear               = Stg_ComponentFactory_GetBool(   cf, self->name, "isNonLinear",            False );
	nonLinearTolerance        = Stg_ComponentFactory_GetDouble( cf, self->name, "nonLinearTolerance",     0.01 );
	nonLinearMaxIterations    = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, "nonLinearMaxIterations", 500 );
	killNonConvergent         = Stg_ComponentFactory_GetBool(   cf, self->name, "killNonConvergent",      True );
	nonLinearMinIterations    = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, "nonLinearMinIterations", 1 );
	nonLinearSolutionType	  = Stg_ComponentFactory_GetString( cf, self->name, "nonLinearSolutionType",  "default" );
	optionsPrefix			  = Stg_ComponentFactory_GetString( cf, self->name, "optionsPrefix",          "" );


	/* Read some value for Picard */
	self->picard_form_function_type = Stg_ComponentFactory_GetString( cf, self->name, "picard_FormFunctionType", "PicardFormFunction_KSPResidual" );

	self->alpha  = Stg_ComponentFactory_GetDouble( cf, self->name, "picard_alpha", 1.0 );
	self->rtol   = Stg_ComponentFactory_GetDouble( cf, self->name, "picard_rtol", 1.0e-8 );
	self->abstol = Stg_ComponentFactory_GetDouble( cf, self->name, "picard_atol", 1.0e-50 );        
	self->xtol           = Stg_ComponentFactory_GetDouble( cf, self->name, "picard_xtol", 1.0e-8 );
	self->picard_monitor = Stg_ComponentFactory_GetBool(   cf, self->name, "picard_ActivateMonitor", False );
	
	entryPointRegister = Stg_ObjectList_Get( cf->registerRegister, "EntryPoint_Register" );
	assert( entryPointRegister );

	context = Stg_ComponentFactory_ConstructByKey( cf, self->name, "Context", FiniteElementContext, False, data );

	if( isNonLinear ) {
		nlSolver = PETScNonlinearSolver_New( "nonLinearSolver" );
		self->linearSolveInitGuess = Stg_ComponentFactory_GetBool( cf, self->name, "linearSolveInitialGuess", False );
	}
	
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
			optionsPrefix,
			entryPointRegister,
			MPI_COMM_WORLD );

	self->X  = PETScVector_New( "delta x" );
	self->F	 = PETScVector_New( "residual" );
	self->J	 = PETScMatrix_New( "Jacobian" ); 
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

/* need to do this before the SLE specific function to set up the 
pre conditioners is called (beginning of solve) */
void SystemLinearEquations_NewtonInitialise( void* _context, void* data ) {
	FiniteElementContext*	context		= (FiniteElementContext*)_context;
	SystemLinearEquations*	sle             = (SystemLinearEquations*)context->slEquations->data[0];
	PETScVector*		F;
	SNES			snes;
	SNES			oldSnes		= ((PETScNonlinearSolver*)sle->nlSolver)->snes;

	/* don't assume that a snes is being used for initial guess, check for this!!! */
	if( oldSnes && context->timeStep == 1 && !sle->linearSolveInitGuess )
		SNESDestroy( oldSnes );

	SNESCreate( sle->comm, &snes );

	((PETScNonlinearSolver*)sle->nlSolver)->snes = snes;
	sle->_setFFunc( &sle->F, context );

	SNESSetJacobian( snes, ((PETScMatrix*)sle->J)->petscMat, ((PETScMatrix*)sle->P)->petscMat, sle->_buildJ, sle->buildJContext );
	SNESSetFunction( snes, ((PETScVector*)sle->F)->petscVec, sle->_buildF, sle->buildFContext );

	/* configure the KSP */
	sle->_configureNLSolverFunc( snes, context );
}

/* do this after the pre conditoiners have been set up in the problem specific SLE */
void SystemLinearEquations_NewtonExecute( void* sle, void* _context ) {
	SystemLinearEquations*	self            = (SystemLinearEquations*) sle;
	SNES			snes		= ((PETScNonlinearSolver*)self->nlSolver)->snes;

	SNESSetOptionsPrefix( snes, self->optionsPrefix );
	SNESSetFromOptions( snes );
	SNESSolve( snes, PETSC_NULL, ((PETScVector*)self->X)->petscVec );
}

/* do this at end of solve step */
void SystemLinearEquations_NewtonFinalise( void* _context, void* data ) {
	FiniteElementContext*	context		= (FiniteElementContext*)_context;
	SystemLinearEquations*	sle             = (SystemLinearEquations*)context->slEquations->data[0];
	SNES			snes		= ((PETScNonlinearSolver*)sle->nlSolver)->snes;

	sle->_updateOldFields( &sle->X, context );
	
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

        /* More of Luke's stuff. I need an entry point for a non-linear setup operation. */
        _EntryPoint_Run_2VoidPtr( self->nlSetupEP, sle, _context );
	
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
//		PetscPrintf( PETSC_COMM_WORLD, "|Xn+1| = %12.12e \n", Vector_L2Norm(SystemLinearEquations_GetSolutionVectorAt(self,1)->vector) );

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

                /*
                ** Include an entry point to do some kind of post-non-linear-iteration operation. */
		_EntryPoint_Run_2VoidPtr( self->postNlEP, sle, _context );
                
		
		if ( (converged) && (self->nonLinearIteration_I>=minIterations) ) {
		   int result;

		   /* Adding in another entry point so we can insert out own custom
		      convergeance checks. For example, with frictional boundary
		      conditions we need to ensure envery node was gone from the
		      original searching state to a fixed slipping or sticking state. */
		   _EntryPoint_Run_2VoidPtr( self->nlConvergedEP, _context, &converged );
		   MPI_Allreduce( &converged, &result, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD );
		   if( result )
		      break;
		}
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

void SystemLinearEquations_AddNonLinearSetupEP( void* sle, const char* name, EntryPoint_2VoidPtr_Cast func ) {
	SystemLinearEquations* self = (SystemLinearEquations*)sle;

	SystemLinearEquations_SetToNonLinear( self );
	EntryPoint_Append( self->nlSetupEP, (char*)name, func, self->type );
}

void SystemLinearEquations_AddPostNonLinearEP( void* sle, const char* name, EntryPoint_2VoidPtr_Cast func ) {
	SystemLinearEquations* self = (SystemLinearEquations*)sle;

	EntryPoint_Append( self->postNlEP, (char*)name, func, self->type );
}

void SystemLinearEquations_AddNonLinearConvergedEP( void* sle,
						    const char* name,
						    EntryPoint_2VoidPtr_Cast func )
{
	SystemLinearEquations* self = (SystemLinearEquations*)sle;

	SystemLinearEquations_SetToNonLinear( self );
	EntryPoint_Append( self->nlConvergedEP, (char*)name, func, self->type );
}



/*
Computes
  F1 := A(x) x -b  = -r,
  where r = b - A(x) x
*/
void SystemLinearEquations_SNESPicardFormalResidual( void *someSLE, Vector *stg_X, Vector *stg_F, void *_context )
{
	SystemLinearEquations *sle = (SystemLinearEquations*)someSLE;
    	SLE_Solver            *solver = (SLE_Solver*)sle->solver;
	Vec                   F;
	Stream*                 errorStream     = Journal_Register( Error_Type, sle->type );

	Journal_Printf( errorStream, "    **** SystemLinearEquations_SNESPicardFormalResidual: This option is un-tested and does not yet function correctly. \n");
	Journal_Printf( errorStream, "    **** Use the default form function or specify --components.XXX.picard_FormFunctionType=PicardFormFunction_KSPResidual instead. \n");
	Journal_Printf( errorStream, "    ****     [Dave May - 12 May, 2008] \n");
	abort();

	solver->_formResidual( (void*)sle,  (void*)solver, stg_F );
	F = StgVectorGetPetscVec( stg_F );	
//	VecScale( F, -1.0 );
}

/*
Computes
  F2 := x - A(x)^{-1} b
*/
#if 0
void SystemLinearEquations_SNESPicardKSPResidual( void *someSLE, Vector *stg_X, Vector *stg_F, void *_context )
{
        SystemLinearEquations *sle = (SystemLinearEquations*)someSLE;
        SLE_Solver            *solver = (SLE_Solver*)sle->solver;
        Vec                   F,X,Xcopy;
        PetscReal norm;

        F = StgVectorGetPetscVec( stg_F );
        X = StgVectorGetPetscVec( stg_X );

        VecDuplicate( X, &Xcopy );
        VecCopy( X, Xcopy );

        VecCopy( X, F );                        /* F <- X */
//      VecNorm( X, NORM_2, &norm );
//      PetscPrintf(PETSC_COMM_WORLD,"  |X|_pre = %5.5e \n", norm );
        sle->linearExecute( sle, _context );    /* X = A^{-1} b */
        X = StgVectorGetPetscVec( stg_X );
//      VecNorm( X, NORM_2, &norm );
//      PetscPrintf(PETSC_COMM_WORLD,"  |X|_post = %5.5e \n", norm );

        VecAXPY( F, -1.0, X );                  /* F <- X - F */
//      VecNorm( F, NORM_2, &norm );
//      PetscPrintf(PETSC_COMM_WORLD,"  |F|_post = %5.5e \n", norm );

        VecCopy( Xcopy, X );
        VecDestroy( Xcopy );
}
#endif

/*
stg_X must not get modified by this function !!
*/
void SystemLinearEquations_SNESPicardKSPResidual( void *someSLE, Vector *stg_X, Vector *stg_F, void *_context )
{
        SystemLinearEquations *sle = (SystemLinearEquations*)someSLE;
        SLE_Solver            *solver = (SLE_Solver*)sle->solver;
        Vector                *stg_Xstar;
        Vec                   F,X,Xstar;
	PetscReal norm,norms;

	solver->_getSolution( sle, solver, &stg_Xstar );

        F = StgVectorGetPetscVec( stg_F );
	X = StgVectorGetPetscVec( stg_X );
	Xstar = StgVectorGetPetscVec( stg_Xstar );
	
	/* Map most current solution into stg object, vec->mesh  */
	VecCopy( X, Xstar );  /* X* <- X */
	/* Map onto nodes */
	SystemLinearEquations_UpdateSolutionOntoNodes( someSLE, _context );

//	VecCopy( X, F );      /*  F <- X */

	VecNorm( X, NORM_2, &norm );
	VecNorm( Xstar, NORM_2, &norms );
//	PetscPrintf(PETSC_COMM_WORLD,"  |X| = %12.12e : |x*| = %12.12e <pre>\n", norm, norms );

	sle->linearExecute( sle, _context );    /* X* = A^{-1} b */
	Xstar = StgVectorGetPetscVec( stg_Xstar );

	VecNorm( X, NORM_2, &norm );
	VecNorm( Xstar, NORM_2, &norms );
//	PetscPrintf(PETSC_COMM_WORLD,"  |X| = %12.12e : |x*| = %12.12e <post> \n", norm, norms );

	VecWAXPY( F, -1.0, Xstar, X ); /* F = -X* + X  */
}

void SLEComputeFunction( void *someSLE, Vector *stg_X, Vector *stg_F, void *_context )
{
        SystemLinearEquations *sle = (SystemLinearEquations*)someSLE;	

	if (sle->_sleFormFunction!=NULL) {
		sle->_sleFormFunction( sle, stg_X, stg_F, _context );
	}
	else {
		SETERRQ( PETSC_ERR_SUP, "SLEComputeFunction in not valid" );
	}
}

void SLE_SNESMonitor( void *sle, PetscInt iter, PetscReal fnorm )
{
  PetscPrintf( PETSC_COMM_WORLD, "  %.4d SLE_NS Function norm %12.12e    --------------------------------------------------------------------------------\n", iter, fnorm );
}

void SLE_SNESMonitor2( void *sle, PetscInt iter, PetscReal fnorm0, PetscReal fnorm, PetscReal dX, PetscReal X1 )
{
  if(iter==0) {
    PetscPrintf( PETSC_COMM_WORLD, "  SLE_NS  it       |F|              |F|/|F0|          |X1-X0|        |X1-X0|/|X1| \n" );
  }
  PetscPrintf( PETSC_COMM_WORLD,   "  SLE_NS  %1.4d     %2.4e      %2.4e       %2.4e     %2.4e \n", iter, fnorm, fnorm/fnorm0, dX, dX/X1 );
}

void _monitor_progress( PetscReal initial, PetscReal target, PetscReal current, PetscReal *p )
{
  PetscReal p0;

  p0 = log10(initial) - log10(target);
  *p = 100.0 * ( 1.0 - (log10(current)-log10(target)) / p0 );
}

void SLE_SNESMonitorProgress( void *sle, 
  PetscInt iter, 
  PetscReal fnorm0, PetscReal fnorm, PetscReal dX, PetscReal X1, 
  PetscReal fatol, PetscReal frtol, PetscReal xtol )
{
  PetscReal f_abs_s, f_abs_e, w1, v1, p1;
  PetscReal f_rel_s, f_rel_e, w2, v2, p2;
  PetscReal x_del_s, x_del_e, w3, v3, p3;

  if(iter==0) {
    PetscPrintf( PETSC_COMM_WORLD, "  SLE_NS  it       |F|                 |F|/|F0|             |X1-X0|/|X1| \n" );
  }

  f_abs_s = fnorm0;
  f_abs_e = fatol;
  v1 = fnorm;
  _monitor_progress( f_abs_s, f_abs_e, v1, &p1 );

  f_rel_s = fnorm0;
  f_rel_e = fnorm0 * frtol;
  v2 = fnorm;
  _monitor_progress( f_rel_s, f_rel_e, v2, &p2 );

  x_del_s = 1.0;
  x_del_e = xtol;
  v3 = dX/X1;
  _monitor_progress( x_del_s, x_del_e, v3, &p3 );

  PetscPrintf( PETSC_COMM_WORLD,   "  SLE_NS  %1.4d     %2.4e [%2.0f%%]     %2.4e [%2.0f%%]      %2.4e [%2.0f%%] \n", iter, fnorm,p1, fnorm/fnorm0,p2,  dX/X1,p3 );
}


void SLE_SNESConverged(
	PetscReal snes_abstol, PetscReal snes_rtol, PetscReal snes_ttol, PetscReal snes_xtol,
	PetscInt it,PetscReal xnorm,PetscReal pnorm,PetscReal fnorm,SNESConvergedReason *reason )
{
  PetscErrorCode ierr;

  *reason = SNES_CONVERGED_ITERATING;

  if (!it) {
    /* set parameter for default relative tolerance convergence test */
    snes_ttol = fnorm*snes_rtol;
  }
  if (fnorm != fnorm) {
//    ierr = PetscInfo(snes,"Failed to converged, function norm is NaN\n");CHKERRQ(ierr);
    *reason = SNES_DIVERGED_FNORM_NAN;
  } else if (fnorm < snes_abstol) {
//    ierr = PetscInfo2(snes,"Converged due to function norm %G < %G\n",fnorm,snes_abstol);CHKERRQ(ierr);
    *reason = SNES_CONVERGED_FNORM_ABS;
  }
//  } else if (snes->nfuncs >= snes->max_funcs) {
//    ierr = PetscInfo2(snes,"Exceeded maximum number of function evaluations: %D > %D\n",snes->nfuncs,snes->max_funcs);CHKERRQ(ierr);
//    *reason = SNES_DIVERGED_FUNCTION_COUNT;
//  }

  if (it && !*reason) {
    if (fnorm <= snes_ttol) {
//      ierr = PetscInfo2(snes,"Converged due to function norm %G < %G (relative tolerance)\n",fnorm,snes_ttol);CHKERRQ(ierr);
      *reason = SNES_CONVERGED_FNORM_RELATIVE;
    } else if (pnorm < snes_xtol*xnorm) {
//      ierr = PetscInfo3(snes,"Converged due to small update length: %G < %G * %G\n",pnorm,snes_xtol,xnorm);CHKERRQ(ierr);
      *reason = SNES_CONVERGED_PNORM_RELATIVE;
    }
  }
}


#if defined(PETSC_HAVE_ISINF) && defined(PETSC_HAVE_ISNAN)
#define PetscIsInfOrNanScalar(a) (isinf(PetscAbsScalar(a)) || isnan(PetscAbsScalar(a)))
#define PetscIsInfOrNanReal(a) (isinf(a) || isnan(a))
#elif defined(PETSC_HAVE__FINITE) && defined(PETSC_HAVE__ISNAN)
#if defined(PETSC_HAVE_FLOAT_H)
#include "float.h"  /* windows defines _finite() in float.h */ 
#endif
#define PetscIsInfOrNanScalar(a) (!_finite(PetscAbsScalar(a)) || _isnan(PetscAbsScalar(a)))
#define PetscIsInfOrNanReal(a) (!_finite(a) || _isnan(a))
#else
#define PetscIsInfOrNanScalar(a) ((a - a) != 0.0)
#define PetscIsInfOrNanReal(a) ((a - a) != 0.0)
#endif

/*
This will be replaced by SNESPicard in petsc 2.4.0.
*/

void SystemLinearEquations_PicardExecute( void *sle, void *_context )
{
  SystemLinearEquations *self = (SystemLinearEquations*)sle;
  SLE_Solver            *solver = (SLE_Solver*)self->solver;

  Vec            X, Y, F, Xstar,delta_X;
  PetscReal      alpha = 1.0;
  PetscReal      fnorm,norm_X,pnorm,fnorm0;
  PetscInt       i;
  PetscErrorCode ierr;
  SNESConvergedReason snes_reason;
 
  PetscReal snes_norm;
  PetscInt snes_iter;

  PetscReal snes_ttol, snes_rtol, snes_abstol, snes_xtol;
  PetscInt  snes_maxits;

  Vector *stg_X, *stg_F, *stg_Xstar;
  PetscTruth monitor_flg;

  /* setup temporary some vectors */
  solver->_getSolution( self, solver, &stg_Xstar );
  Xstar = StgVectorGetPetscVec( stg_Xstar );

  Vector_Duplicate( stg_Xstar, (void**)&stg_X );
  Vector_SetLocalSize( stg_X, Vector_GetLocalSize(stg_Xstar) );
  X = StgVectorGetPetscVec( stg_X );

  Vector_Duplicate( stg_X, (void**)&stg_F );
  Vector_SetLocalSize( stg_F, Vector_GetLocalSize(stg_X) );
  F = StgVectorGetPetscVec( stg_F );

  VecDuplicate( F, &Y );
  VecDuplicate( F, &delta_X );

  /* Get some values from dictionary */
  snes_maxits =  (PetscInt)self->nonLinearMaxIterations;
  snes_ttol   = 0.0;
  snes_rtol   = (PetscReal)self->rtol;
  snes_abstol = (PetscReal)self->abstol;
  snes_xtol   = (PetscReal)self->xtol;
  monitor_flg = PETSC_FALSE;
  if (self->picard_monitor==True) { monitor_flg = PETSC_TRUE; }
  alpha       = (PetscReal)self->alpha;

  snes_reason = SNES_CONVERGED_ITERATING;

  /* Map X <- X* */
  VecCopy( Xstar, X );  // Vector_CopyEntries( currentVector, previousVector );

  /* Get an initial guess if |X| ~ 0, by solving the linear problem  */
  VecNorm( X, NORM_2, &norm_X );
  if (norm_X <1.0e-20) {
    if (monitor_flg==PETSC_TRUE)
      PetscPrintf( PETSC_COMM_WORLD, "SLE_Picard: Computing an initial guess for X from the linear problem\n");

    self->linearExecute( sle, _context );    /* X* = A^{-1} b */
    self->hasExecuted = True;

    Xstar = StgVectorGetPetscVec( stg_Xstar );
    
    /* Map X <- X* */
    VecCopy( Xstar, X ); 
    VecNorm( X, NORM_2, &norm_X );
  }


  snes_iter = 0;
  snes_norm = 0;
  SLEComputeFunction( sle, stg_X, stg_F, _context );
  ierr = VecNorm(F, NORM_2, &fnorm);CHKERRQ(ierr); /* fnorm <- ||F||  */
  fnorm0 = fnorm;
  if( PetscIsInfOrNanReal(fnorm) ) SETERRQ(PETSC_ERR_FP,"Infinite or not-a-number generated in norm");

  snes_norm = fnorm;
  if(monitor_flg==PETSC_TRUE) {
  /*  SLE_SNESMonitor(sle,0,fnorm); */
  /*  SLE_SNESMonitor2(sle,0,fnorm0,fnorm, norm_X, norm_X ); */
    SLE_SNESMonitorProgress( sle, snes_iter,  fnorm0, fnorm, norm_X, norm_X, snes_abstol, snes_rtol, snes_xtol );
  }

  /* set parameter for default relative tolerance convergence test */
  snes_ttol = fnorm*snes_rtol;
  /* test convergence */
  SLE_SNESConverged( snes_abstol,snes_rtol,snes_ttol,snes_xtol , 0,norm_X,0.0,fnorm,&snes_reason);
  if (snes_reason) return;

  for(i = 0; i < snes_maxits; i++) {
    PetscTruth lsSuccess = PETSC_TRUE;
    PetscReal  dummyNorm;

    /* Update guess Y = X^n - F(X^n) */
    ierr = VecWAXPY(Y, -1.0, F, X);CHKERRQ(ierr);

    VecCopy( X, delta_X );  /* delta_X <- X */

    /* X^{n+1} = (1 - \alpha) X^n + alpha Y */
    ierr = VecAXPBY(X, alpha, 1 - alpha, Y);CHKERRQ(ierr);
    VecNorm( X, NORM_2, &norm_X ); 
/*    PetscPrintf( PETSC_COMM_WORLD, "  Xn+1 = %12.12e \n", norm_X ); */

    VecAYPX( delta_X, -1.0, X );   /* delta_X <- Xn+1 - delta_X */
    VecNorm( delta_X, NORM_2, &pnorm );

    /* Compute F(X^{new}) */
    SLEComputeFunction( sle, stg_X, stg_F, _context );
    ierr = VecNorm(F, NORM_2, &fnorm);CHKERRQ(ierr);
    if( PetscIsInfOrNanReal(fnorm) ) SETERRQ(PETSC_ERR_FP,"Infinite or not-a-number generated norm");

    /* Monitor convergence */
    snes_iter = i+1;
    snes_norm = fnorm;
    if (monitor_flg==PETSC_TRUE) {
    /*  SLE_SNESMonitor(sle,snes_iter,snes_norm); */
    /*  SLE_SNESMonitor2(sle,snes_iter,fnorm0,fnorm, pnorm, norm_X ); */
      SLE_SNESMonitorProgress( sle, snes_iter,  fnorm0, fnorm, pnorm, norm_X, snes_abstol, snes_rtol, snes_xtol );
    }

    /* Test for convergence */
    SLE_SNESConverged( snes_abstol,snes_rtol,snes_ttol,snes_xtol , snes_iter,norm_X,pnorm,fnorm,&snes_reason);
    if (snes_reason) break;
  }
  if (i == snes_maxits) {
//    ierr = PetscInfo1(snes, "Maximum number of iterations has been reached: %D\n", maxits);CHKERRQ(ierr);
    if (!snes_reason) snes_reason = SNES_DIVERGED_MAX_IT;
  }

  /* If monitoring, report reason converged */
  if (monitor_flg==PETSC_TRUE)
    PetscPrintf( PETSC_COMM_WORLD, "Nonlinear solve converged due to %s \n", SNESConvergedReasons[snes_reason] );

  FreeObject( stg_X );
  FreeObject( stg_F );
  VecDestroy( Y );
  VecDestroy( delta_X );

}


///////////////

void SystemLinearEquations_AddNonLinearEP( void* sle, const char* name, EntryPoint_2VoidPtr_Cast func ) {
	SystemLinearEquations* self = (SystemLinearEquations*)sle;

	SystemLinearEquations_SetToNonLinear( self );
	EntryPoint_Append( self->nlEP, (char*)name, func, self->type );
}



void SystemLinearEquations_SetToNonLinear( void* sle ) {
	SystemLinearEquations*	self            	= (SystemLinearEquations*) sle;
	Hook*			nonLinearInitHook	= NULL;
	Hook*			nonLinearFinaliseHook	= NULL;
	FiniteElementContext*	context			= NULL;

	assert( self );
	if ( self->isNonLinear )
		return;

	self->isNonLinear = True;

	self->linearExecute = self->_execute;
	self->_execute = SystemLinearEquations_NonLinearExecute;

	if( self->nonLinearSolutionType ) {
		if( !strcmp( self->nonLinearSolutionType, "default" ) ) {
			self->_execute = SystemLinearEquations_NonLinearExecute;
		}		

		if( !strcmp( self->nonLinearSolutionType, "MatrixFreeNewton" ) )
			self->_execute = SystemLinearEquations_NewtonMFFDExecute;

		if( !strcmp( self->nonLinearSolutionType, "Newton" ) ) {
			context = self->context;

			nonLinearInitHook = Hook_New( "NewtonInitialise", 
						SystemLinearEquations_NewtonInitialise, self->name );
			_EntryPoint_PrependHook_AlwaysFirst( Context_GetEntryPoint( context, AbstractContext_EP_Solve ), 
								nonLinearInitHook );
			nonLinearFinaliseHook = Hook_New( "NewtonFinalise",
						SystemLinearEquations_NewtonFinalise, self->name );
			_EntryPoint_AppendHook_AlwaysLast( Context_GetEntryPoint( context, AbstractContext_EP_Solve ), 
						nonLinearFinaliseHook );
			self->_execute = SystemLinearEquations_NewtonExecute;
		}
		
		if (!strcmp( self->nonLinearSolutionType, "Picard") ) {
			/* set function pointer for execute */
			self->_execute = SystemLinearEquations_PicardExecute;

			/* set form function */
			if (!strcmp(self->picard_form_function_type,"PicardFormFunction_KSPResidual") ) {
				self->_sleFormFunction = SystemLinearEquations_SNESPicardKSPResidual;
			}
			else if (!strcmp(self->picard_form_function_type,"PicardFormFunction_FormalResidual") ) {
				self->_sleFormFunction = SystemLinearEquations_SNESPicardFormalResidual;
			}
			else {
				 Stream *errorStream = Journal_Register( Error_Type, self->type );

		                Journal_Printf( errorStream, "Unknown the Picard FormFunction type %s is unrecognised. .\n", self->picard_form_function_type );
				Journal_Printf( errorStream, "Supported types include <PicardFormFunction_FormalResidual, PicardFormFunction_KSPResidual> \n" );
                        	abort();

			}
		}
	}
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
