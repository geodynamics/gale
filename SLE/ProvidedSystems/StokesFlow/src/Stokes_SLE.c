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
** $Id: Stokes_SLE.c 999 2008-01-09 04:13:42Z DavidLee $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include "StgFEM/Discretisation/Discretisation.h"
#include "StgFEM/SLE/SystemSetup/SystemSetup.h"
#include "types.h"
#include "Stokes_SLE.h"
#include "UpdateDt.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

/* Textual name of this class */
const Type Stokes_SLE_Type = "Stokes_SLE";

Stokes_SLE* Stokes_SLE_New( 
	Name							name,
	SLE_Solver*					solver,
	FiniteElementContext*	context,
	Bool							isNonLinear,
	double						nonLinearTolerance,
	Iteration_Index			nonLinearMaxIterations,
	Bool							killNonConvergent,		
	EntryPoint_Register*		entryPoint_Register,
	MPI_Comm						comm,
	StiffnessMatrix*			kStiffMat,
	StiffnessMatrix*			gStiffMat,
	StiffnessMatrix*			dStiffMat,
	StiffnessMatrix*			cStiffMat,
	SolutionVector*			uSolnVec,
	SolutionVector*			pSolnVec,
	ForceVector*				fForceVec,
	ForceVector*				hForceVec )
{
	Stokes_SLE* self = (Stokes_SLE*) _Stokes_SLE_DefaultNew( name );

	self->isConstructed = True;
	_SystemLinearEquations_Init( self, solver, NULL, context, False, isNonLinear, 
      nonLinearTolerance, nonLinearMaxIterations, killNonConvergent, 1,  "", "", entryPoint_Register, comm ); 
	_Stokes_SLE_Init( self, kStiffMat, gStiffMat, dStiffMat, cStiffMat, uSolnVec, pSolnVec, fForceVec, hForceVec );

	return self;
}

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
	Name                                                name )
{
	Stokes_SLE* self;
	
	/* Allocate memory */
	assert( sizeOfSelf >= sizeof(Stokes_SLE) );
	self = (Stokes_SLE*) _SystemLinearEquations_New( 
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
		_LM_Setup,
		_matrixSetup,
		_vectorSetup,
		_updateSolutionOntoNodes, 
		_mgSelectStiffMats, 
		name );
	
	/* Virtual info */
	
	return self;
}

void _Stokes_SLE_Init( 
		void*                                               sle, 
		StiffnessMatrix*                                    kStiffMat,
		StiffnessMatrix*                                    gStiffMat,
		StiffnessMatrix*                                    dStiffMat,
		StiffnessMatrix*                                    cStiffMat,
		SolutionVector*                                     uSolnVec,
		SolutionVector*                                     pSolnVec,
		ForceVector*                                        fForceVec,
		ForceVector*                                        hForceVec ) 
{
	Stokes_SLE* self = (Stokes_SLE*)sle;

	self->kStiffMat = kStiffMat;
	self->gStiffMat = gStiffMat;
	self->dStiffMat = dStiffMat;
	self->cStiffMat = cStiffMat;
	self->uSolnVec  = uSolnVec; 
	self->pSolnVec  = pSolnVec;
	self->fForceVec = fForceVec;   
	self->hForceVec = hForceVec;   

	/* add the vecs and matrices to the Base SLE class's dynamic lists, so they can be 
	initialised and built properly */
	SystemLinearEquations_AddStiffnessMatrix( self, kStiffMat );
	SystemLinearEquations_AddStiffnessMatrix( self, gStiffMat );

	if ( dStiffMat ) 
		SystemLinearEquations_AddStiffnessMatrix( self, dStiffMat );
	if ( cStiffMat ) 
		SystemLinearEquations_AddStiffnessMatrix( self, cStiffMat );

	SystemLinearEquations_AddSolutionVector( self, uSolnVec );
	SystemLinearEquations_AddSolutionVector( self, pSolnVec );
	SystemLinearEquations_AddForceVector( self, fForceVec );
	SystemLinearEquations_AddForceVector( self, hForceVec );

	/* Put timestep function onto entry point */
	if ( self->context )
		EP_AppendClassHook( self->context->calcDtEP, Stokes_SLE_UpdateDt, self );
}

void _Stokes_SLE_Print( void* sle, Stream* stream ) {
	Stokes_SLE* self = (Stokes_SLE*)sle;
	/* Set the Journal for printing informations */
	
	/* General info */
	Journal_Printf( stream, "Stokes_SLE (ptr): %p\n", self );
	_SystemLinearEquations_Print( self, stream );

	Journal_Printf( stream, "Name of discrete stress tensor (K) matrix = \"%s\" \n",self->kStiffMat->name );
	Journal_Printf( stream, "Name of discrete gradient (G) matrix = \"%s\" \n",self->gStiffMat->name );
	if (self->dStiffMat) {
		Journal_Printf( stream, "Name of discrete divergence (D) matrix = \"%s\" \n",self->dStiffMat->name );
	}
	else {
		Journal_Printf( stream, "No discrete divergence (D) matrix set up (Symmetric geometry).\n" );
	}
	if (self->cStiffMat) {
		Journal_Printf( stream, "Name of compressibility (C) matrix = \"%s\" \n",self->cStiffMat->name );
	}
	else {
		Journal_Printf( stream, "No compressibility (C) matrix set up (incompressible fluids)\n" );
	}
	Journal_Printf( stream, "Name of velocity (u) vector = \"%s\" \n",self->uSolnVec->name );
	Journal_Printf( stream, "Name of pressure (p) vector = \"%s\" \n",self->pSolnVec->name );

	Stg_Class_Print( self->solver, stream );
}

void* _Stokes_SLE_DefaultNew( Name name ) {
	return (void*)_Stokes_SLE_New( 
		sizeof(Stokes_SLE), 
		Stokes_SLE_Type,
		_SystemLinearEquations_Delete,
		_Stokes_SLE_Print,
		_SystemLinearEquations_Copy,
		_Stokes_SLE_DefaultNew,
		_Stokes_SLE_AssignFromXML,
		_SystemLinearEquations_Build,
		_SystemLinearEquations_Initialise,
		_SystemLinearEquations_Execute,
		_SystemLinearEquations_Destroy,
		_SystemLinearEquations_LM_Setup,
		_SystemLinearEquations_MatrixSetup,
		_SystemLinearEquations_VectorSetup,
		_SystemLinearEquations_UpdateSolutionOntoNodes,
		_Stokes_SLE_MG_SelectStiffMats, 
		name );
}

void _Stokes_SLE_AssignFromXML( void* sle, Stg_ComponentFactory* cf, void* data ) {
	Stokes_SLE*       self = (Stokes_SLE*)sle;
	StiffnessMatrix*  kStiffMat;
	StiffnessMatrix*  gStiffMat;
	StiffnessMatrix*  dStiffMat;
	StiffnessMatrix*  cStiffMat;
	SolutionVector*   uSolnVec;
	SolutionVector*   pSolnVec;
	ForceVector*      fForceVec;
	ForceVector*      hForceVec;

	/* Construct Parent */
	_SystemLinearEquations_AssignFromXML( self, cf, data );

	kStiffMat =  Stg_ComponentFactory_ConstructByKey( cf, self->name, "StressTensorMatrix",    StiffnessMatrix, True, data );
	gStiffMat =  Stg_ComponentFactory_ConstructByKey( cf, self->name, "GradientMatrix",        StiffnessMatrix, True, data );
	dStiffMat =  Stg_ComponentFactory_ConstructByKey( cf, self->name, "DivergenceMatrix",      StiffnessMatrix, False, data );
	cStiffMat =  Stg_ComponentFactory_ConstructByKey( cf, self->name, "CompressibilityMatrix", StiffnessMatrix, False, data );

	uSolnVec  =  Stg_ComponentFactory_ConstructByKey( cf, self->name, "VelocityVector",        SolutionVector,  True, data );
	pSolnVec  =  Stg_ComponentFactory_ConstructByKey( cf, self->name, "PressureVector",        SolutionVector,  True, data );

	fForceVec =  Stg_ComponentFactory_ConstructByKey( cf, self->name, "ForceVector",           ForceVector,     True, data );
	hForceVec =  Stg_ComponentFactory_ConstructByKey( cf, self->name, "ContinuityForceVector", ForceVector,     True, data );

	_Stokes_SLE_Init( self, kStiffMat, gStiffMat, dStiffMat, cStiffMat, uSolnVec, pSolnVec, fForceVec, hForceVec );
}


void _Stokes_SLE_MG_SelectStiffMats( void* _sle, unsigned* nSMs, StiffnessMatrix*** sms ) {
	Stokes_SLE*	self = (Stokes_SLE*)_sle;
	
	/*
	** In this implementation, only the velocity matrix will have MG applied.
	*/
	
	*nSMs = 1;
	*sms = Memory_Alloc_Array( StiffnessMatrix*, 1, "Stokes_SLE" );
	(*sms)[0] = self->kStiffMat;
}

