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
** $Id: SLE_Solver.c 1125 2008-05-12 14:22:02Z DavidMay $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include "StgFEM/Discretisation/Discretisation.h"
#include "units.h"
#include "types.h"
#include "shortcuts.h"
#include "FiniteElementContext.h"
#include "SLE_Solver.h"

#include <assert.h>

/** Textual name of this class */
const Type SLE_Solver_Type = "SLE_Solver";

SLE_Solver* _SLE_Solver_New(  SLE_SOLVER_DEFARGS  )
{	
	SLE_Solver*		self;

	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(SLE_Solver) );
	/* The following terms are parameters that have been passed into this function but are being set before being passed onto the parent */
	/* This means that any values of these parameters that are passed into this function are not passed onto the parent function
	   and so should be set to ZERO in any children of this class. */
	nameAllocationType = NON_GLOBAL;

	self = (SLE_Solver*) _Stg_Component_New(  STG_COMPONENT_PASSARGS  );
	
	/* General info */
	
	/* Virtual info */
	self->_solverSetup = _solverSetup;
	self->_solve = _solve;
	self->_getResidual = _getResidual;

	self->_formResidual = NULL;
	self->_getRhs       = NULL;
	self->_getSolution  = NULL;

	return self;

}

void _SLE_Solver_Init( SLE_Solver* self, Bool useStatSolve, int statReps ) {
	self->isConstructed = True;
	self->extensionManager = ExtensionManager_New_OfExistingObject( self->name, self );
	
	self->debug         = Stream_RegisterChild( StgFEM_SLE_SystemSetup_Debug, self->type );
	self->info          = Journal_MyStream( Info_Type, self );
	self->maxIterations = 0;

	self->previoustimestep = 0;
	self->currenttimestep = 0;
	self->nonlinearitsinitialtime = 0; 
	self->nonlinearitsendtime = 0; 
	self->totalnonlinearitstime = 0; 
	self->totalnumnonlinearits = 0; 
	self->avgtimenonlinearits = 0; 	
	self->inneritsinitialtime = 0; 
	self->outeritsinitialtime = 0; 
	self->inneritsendtime = 0; 
	self->outeritsendtime = 0; 
	self->totalinneritstime = 0; 
	self->totalouteritstime = 0; 
	self->totalnuminnerits = 0; 
	self->totalnumouterits = 0; 
	self->avgnuminnerits = 0; 
	self->avgnumouterits = 0; 
	self->avgtimeinnerits = 0; 
	self->avgtimeouterits = 0;
	
	self->useStatSolve = useStatSolve;
	self->nStatReps     = statReps;
}

void SLE_Solver_InitAll( void* sleSolver, Bool useStatSolve, int statReps ) {
	_SLE_Solver_Init( (SLE_Solver*) sleSolver, useStatSolve, statReps );
}

void* _SLE_Solver_Copy( const void* sleSolver, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	SLE_Solver*		self = (SLE_Solver*)sleSolver;
	SLE_Solver*		newSleSolver;
	PtrMap*			map = ptrMap;
	Bool			ownMap = False;
	
	if( !map ) {
		map = PtrMap_New( 10 );
		ownMap = True;
	}
	
	newSleSolver = (SLE_Solver*)_Stg_Component_Copy( self, dest, deep, nameExt, map );
	
	/* virtual functions */
	newSleSolver->_solverSetup  = self->_solverSetup;
	newSleSolver->_solve        = self->_solve;
	newSleSolver->maxIterations = self->maxIterations;

	newSleSolver->inneritsinitialtime = self->inneritsinitialtime;
	newSleSolver->outeritsinitialtime = self->outeritsinitialtime;
	newSleSolver->nonlinearitsinitialtime = self->nonlinearitsinitialtime;
	newSleSolver->inneritsendtime = self->inneritsendtime;
	newSleSolver->outeritsendtime = self->outeritsendtime;
	newSleSolver->nonlinearitsendtime = self->nonlinearitsendtime;
	newSleSolver->totalinneritstime = self->totalinneritstime;
	newSleSolver->totalouteritstime = self->totalouteritstime;
	newSleSolver->totalnonlinearitstime = self->totalnonlinearitstime;
	newSleSolver->totalnuminnerits = self->totalnuminnerits; 
	newSleSolver->totalnumouterits = self->totalnumouterits; 
	newSleSolver->totalnumnonlinearits = self->totalnumnonlinearits; 	
	newSleSolver->avgnuminnerits = self->avgnuminnerits;
    newSleSolver->avgnumouterits = self->avgnumouterits;
	newSleSolver->avgtimeinnerits = self->avgtimeinnerits; 
	newSleSolver->avgtimeouterits = self->avgtimeouterits; 
	newSleSolver->avgtimenonlinearits = self->avgtimenonlinearits; 
	newSleSolver->currenttimestep = self->currenttimestep; 
	newSleSolver->previoustimestep = self->previoustimestep;
	
	if( deep ) {
          if( (newSleSolver->debug = (Stream*)PtrMap_Find( map, self->debug )) == NULL ) {
			newSleSolver->debug = (Stream*)Stg_Class_Copy( self->debug, NULL, deep, nameExt, map );
			PtrMap_Append( map, self->debug, newSleSolver->debug );
		}
          if( (newSleSolver->extensionManager = (ExtensionManager*)PtrMap_Find( map, self->extensionManager )) == NULL ) {
			newSleSolver->extensionManager = (ExtensionManager*)Stg_Class_Copy( self->extensionManager, NULL, deep, nameExt, map );
			PtrMap_Append( map, self->extensionManager, newSleSolver->extensionManager );
		}
	}
	else {
		newSleSolver->debug = self->debug;
		newSleSolver->extensionManager = (ExtensionManager*)Stg_Class_Copy( self->extensionManager, NULL, deep, nameExt, map );
	}
	
	if( ownMap ) {
		Stg_Class_Delete( map );
	}
	
	return (void*)newSleSolver;
}


void _SLE_Solver_Delete( void* sleSolver ) {
	SLE_Solver* self = (SLE_Solver*)sleSolver;

	_Stg_Component_Delete( self );	
}

void _SLE_Solver_Print( void* sleSolver, Stream* stream ) {
	SLE_Solver*		self = (SLE_Solver*)sleSolver;

	_Stg_Component_Print( self, stream );

	Journal_PrintPointer( stream, self->extensionManager );
	
	Journal_PrintPointer( stream, self->_solverSetup );
	Journal_PrintPointer( stream, self->_solve );
	Journal_PrintPointer( stream, self->_getResidual );

	Journal_PrintPointer( stream, self->debug );
	Journal_PrintValue( stream, self->maxIterations );
}

void _SLE_Solver_Build( void* sleSolver, void* data ) {
	/* Do nothing by default */
}

void _SLE_Solver_AssignFromXML( void* sleSolver, Stg_ComponentFactory* cf, void* data ) {
	SLE_Solver*		self = (SLE_Solver*)sleSolver;
	Bool            useStatSolve;
	int             nStatReps;

	self->context = Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"Context", FiniteElementContext, False, data );
	if( !self->context  )
		self->context = Stg_ComponentFactory_ConstructByName( cf, (Name)"context", FiniteElementContext, True, data  );

	useStatSolve = Stg_ComponentFactory_GetBool( cf, self->name, (Dictionary_Entry_Key)"statSolve", False  );
	nStatReps = Stg_ComponentFactory_GetInt( cf, self->name, (Dictionary_Entry_Key)"statReps", 0  );

	_SLE_Solver_Init( self, useStatSolve, nStatReps );
}

void _SLE_Solver_Initialise( void* sleSolver, void* data ) {
}


void _SLE_Solver_Execute( void* sleSolver, void* data ) {
	SLE_Solver*		self = (SLE_Solver*)sleSolver;

	Journal_DPrintf( self->debug, "In %s()\n", __func__ );

	Stream_IndentBranch( StgFEM_Debug );
	SLE_Solver_SolverSetup( self, data );
	SLE_Solver_Solve( self, data );
	Stream_UnIndentBranch( StgFEM_Debug );
}

void _SLE_Solver_Destroy( void* sleSolver, void* data ) {
	SLE_Solver* self = (SLE_Solver*)sleSolver;

	Stg_Class_Delete( self->extensionManager );
}

void SLE_Solver_SolverSetup( void* sleSolver, void* sle ) {
	SLE_Solver*		self = (SLE_Solver*)sleSolver;
	
	self->_solverSetup( self, sle );
}


void SLE_Solver_Solve( void* sleSolver, void* sle ) {
	SLE_Solver*		self = (SLE_Solver*)sleSolver;
	
	self->_solve( self, sle );
}


