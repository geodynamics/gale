/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003, Victorian Partnership for Advanced Computing (VPAC) Ltd,
** 110 Victoria Street, Melbourne, 3053, Australia.
**
** Authors:
**	Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
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
** $Id: Multigrid.c 2192 2004-10-15 02:45:38Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <mpi.h>
#include <petsc.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscksp.h>
#include <petscsnes.h>

#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include <Underworld/Underworld.h>

#include "Multigrid.h"


const Type Underworld_Multigrid_Type = "Underworld_Multigrid";
Underworld_Multigrid* Underworld_Multigrid_selfPointer = NULL;


void Underworld_Multigrid_SolverSetup( void* _solver, void* _stokesSLE ) {
    Underworld_Multigrid* self = Underworld_Multigrid_selfPointer;
	Stokes_SLE_UzawaSolver* solver = (Stokes_SLE_UzawaSolver*)_solver;
	Stokes_SLE* sle = (Stokes_SLE*)_stokesSLE;
	
 	Journal_DPrintf( solver->debug, "In %s:\n", __func__ );
	Stream_IndentBranch( StgFEM_Debug );

	Journal_DPrintfL( solver->debug, 1, "Setting up MatrixSolver for the velocity eqn.\n" );
   self->mgSolver->mgData->matrixChanged = True;
	//MatrixSolver_SetMatrix( self->velSolver, sle->kStiffMat->matrix );
	//KSPSetOperators( self->velSolver, ((PETScMatrix*)sle->kStiffMat->matrix)->petscMat, 
	//		((PETScMatrix*)sle->kStiffMat->matrix)->petscMat, DIFFERENT_NONZERO_PATTERN );

        self->mgSolver->mgData->matrix = sle->kStiffMat->matrix;
        PETScMGSolver_Setup( self->mgSolver, NULL, NULL );
        solver->velSolver = self->mgSolver->mgData->ksp;
	KSPSetOperators( solver->velSolver, sle->kStiffMat->matrix, sle->kStiffMat->matrix, DIFFERENT_NONZERO_PATTERN );
/*
	KSPSetOperators( solver->velSolver, sle->kStiffMat->matrix, sle->kStiffMat->matrix, DIFFERENT_NONZERO_PATTERN );
        KSPSetFromOptions( solver->velSolver );
*/

	if( solver->pcSolver ) {
		Journal_DPrintfL( solver->debug, 1, "Setting up MatrixSolver for the Preconditioner.\n" );
		//MatrixSolver_SetMatrix( self->pcSolver, self->preconditioner->matrix );
		//KSPSetOperators( self->pcSolver, ((PETScMatrix*)self->preconditioner->matrix)->petscMat, 
		//		((PETScMatrix*)self->preconditioner->matrix)->petscMat, DIFFERENT_NONZERO_PATTERN );
		KSPSetOperators( solver->pcSolver, solver->preconditioner->matrix, solver->preconditioner->matrix, DIFFERENT_NONZERO_PATTERN );
    KSPSetFromOptions( solver->pcSolver );
	}

	Stream_UnIndentBranch( StgFEM_Debug );
}


void Underworld_Multigrid_AssignFromXML( void* _self, Stg_ComponentFactory* cf, void* data ) {
   Underworld_Multigrid* self = (Underworld_Multigrid*)_self;
   UnderworldContext* ctx;

   Underworld_Multigrid_selfPointer = self;

   self->ctx = Stg_ComponentFactory_ConstructByName( cf, "context", UnderworldContext, True, data );
   self->sle = Stg_ComponentFactory_ConstructByName( cf, "stokesEqn", Stokes_SLE, True, data );
   self->mgSolver = Stg_ComponentFactory_ConstructByName( cf, "mgSolver", PETScMGSolver, True, data );

   /* Replace the setup routine... hacky. */
   self->sle->solver->_solverSetup = Underworld_Multigrid_SolverSetup;   
   
}

void Underworld_Multigrid_Build( void* _self, void* data ) {
   Underworld_Multigrid* self = (Underworld_Multigrid*)_self;

   Stg_Component_Build( self->mgSolver, data, False );
   Stg_Component_Build( self->sle, data, False );
}

void Underworld_Multigrid_Initialise( void* _self, void* data ) {
   Underworld_Multigrid* self = (Underworld_Multigrid*)_self;

   Stg_Component_Initialise( self->mgSolver, data, False );
   Stg_Component_Initialise( self->sle, data, False );

   /* Setup the MG solver. */
   PETScMGSolver_SetComm( self->mgSolver, MPI_COMM_WORLD );
}

void* Underworld_Multigrid_New( Name name ) {
	/* Variables set in this function */
	SizeT                                              _sizeOfSelf = sizeof(Underworld_Multigrid);
	Type                                                      type = Underworld_Multigrid_Type;
	Stg_Class_DeleteFunction*                              _delete = _Codelet_Delete;
	Stg_Class_PrintFunction*                                _print = _Codelet_Print;
	Stg_Class_CopyFunction*                                  _copy = _Codelet_Copy;
	Stg_Component_DefaultConstructorFunction*  _defaultConstructor = Underworld_Multigrid_New;
	Stg_Component_ConstructFunction*                    _construct = Underworld_Multigrid_AssignFromXML;
	Stg_Component_BuildFunction*                            _build = Underworld_Multigrid_Build;
	Stg_Component_InitialiseFunction*                  _initialise = Underworld_Multigrid_Initialise;
	Stg_Component_ExecuteFunction*                        _execute = _Codelet_Execute;
	Stg_Component_DestroyFunction*                        _destroy = _Codelet_Destroy;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

   return _Codelet_New(  CODELET_PASSARGS  );
}

Index Underworld_Multigrid_Register( PluginsManager* mgr ) {
   return PluginsManager_Submit( mgr,
                                 Underworld_Multigrid_Type,
                                 "0",
                                 Underworld_Multigrid_New );
}


