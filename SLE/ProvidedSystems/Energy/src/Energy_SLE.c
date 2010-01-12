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
** $Id: Energy_SLE.c 999 2008-01-09 04:13:42Z DavidLee $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include "StgFEM/Discretisation/Discretisation.h"
#include "StgFEM/SLE/SystemSetup/SystemSetup.h"
#include "types.h"
#include "Energy_SLE.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

/* Textual name of this class */
const Type Energy_SLE_Type = "Energy_SLE";

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
	ForceVector*				fVector )
{
	Energy_SLE* self = _Energy_SLE_DefaultNew( name );

	self->isConstructed = True;
	_SystemLinearEquations_Init( self, solver, NULL, context, False, isNonLinear, nonLinearTolerance,
		nonLinearMaxIterations, killNonConvergent, 1, "", "", entryPoint_Register, comm );
	_Energy_SLE_Init( self, stiffMat, solutionVec, fVector );

	return self;
}

/* Creation implementation / Virtual constructor */
Energy_SLE* _Energy_SLE_New(  ENERGY_SLE_DEFARGS  )
{
	Energy_SLE* self;
	
	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(Energy_SLE) );
	self = (Energy_SLE*) _SystemLinearEquations_New(  SYSTEMLINEAREQUATIONS_PASSARGS  );
	
	/* Virtual info */
	return self;
}

void _Energy_SLE_Init( void* sle, StiffnessMatrix* stiffMat, SolutionVector* solutionVec, ForceVector* fVector ) {
	Energy_SLE* self = (Energy_SLE*)sle;

	self->stiffMat		= stiffMat;
	self->solutionVec = solutionVec;
	self->fVector     = fVector;

	SystemLinearEquations_AddStiffnessMatrix( self, stiffMat );
	SystemLinearEquations_AddSolutionVector( self, solutionVec );
	SystemLinearEquations_AddForceVector( self, fVector );
}

void _Energy_SLE_Print( void* sle, Stream* stream ) {
	Energy_SLE* self = (Energy_SLE*)sle;
	
	/* General info */
	Journal_Printf( stream, "Energy_SLE (ptr): %p\n", self );
	_SystemLinearEquations_Print( self, stream );

	Journal_PrintString( stream, self->stiffMat->name );
	Journal_PrintString( stream, self->solutionVec->name );
	Journal_PrintString( stream, self->fVector->name );

	Stg_Class_Print( self->solver, stream );
}

void* _Energy_SLE_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                                            _sizeOfSelf = sizeof(Energy_SLE);
	Type                                                                    type = Energy_SLE_Type;
	Stg_Class_DeleteFunction*                                            _delete = _SystemLinearEquations_Delete;
	Stg_Class_PrintFunction*                                              _print = _Energy_SLE_Print;
	Stg_Class_CopyFunction*                                                _copy = _SystemLinearEquations_Copy;
	Stg_Component_DefaultConstructorFunction*                _defaultConstructor = _Energy_SLE_DefaultNew;
	Stg_Component_ConstructFunction*                                  _construct = _Energy_SLE_AssignFromXML;
	Stg_Component_BuildFunction*                                          _build = _SystemLinearEquations_Build;
	Stg_Component_InitialiseFunction*                                _initialise = _SystemLinearEquations_Initialise;
	Stg_Component_ExecuteFunction*                                      _execute = _SystemLinearEquations_Execute;
	Stg_Component_DestroyFunction*                                      _destroy = _SystemLinearEquations_Destroy;
	SystemLinearEquations_LM_SetupFunction*                            _LM_Setup = _SystemLinearEquations_LM_Setup;
	SystemLinearEquations_MatrixSetupFunction*                      _matrixSetup = _SystemLinearEquations_MatrixSetup;
	SystemLinearEquations_VectorSetupFunction*                      _vectorSetup = _SystemLinearEquations_VectorSetup;
	SystemLinearEquations_UpdateSolutionOntoNodesFunc*  _updateSolutionOntoNodes = _SystemLinearEquations_UpdateSolutionOntoNodes;
	SystemLinearEquations_MG_SelectStiffMatsFunc*             _mgSelectStiffMats = _SystemLinearEquations_MG_SelectStiffMats;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

	return (void*)_Energy_SLE_New(  ENERGY_SLE_PASSARGS  );
}

void _Energy_SLE_AssignFromXML( void* sle, Stg_ComponentFactory* cf, void* data ){
	Energy_SLE*      self = (Energy_SLE*)sle;
	StiffnessMatrix* stiffMat;
	SolutionVector*  solutionVec;
	ForceVector*     fVector;

	/* Construct Parent */
	_SystemLinearEquations_AssignFromXML( self, cf, data );

	stiffMat    =  Stg_ComponentFactory_ConstructByKey( cf, self->name, StiffnessMatrix_Type, StiffnessMatrix, True, data  ) ;
	solutionVec =  Stg_ComponentFactory_ConstructByKey( cf, self->name, SolutionVector_Type,  SolutionVector,  True, data  ) ;
	fVector     =  Stg_ComponentFactory_ConstructByKey( cf, self->name, ForceVector_Type,     ForceVector,     True, data  ) ;

	_Energy_SLE_Init( self, stiffMat, solutionVec, fVector );
}


