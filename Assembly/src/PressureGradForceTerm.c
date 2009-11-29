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
** $Id: PressureGradForceTerm.c 822 2007-04-27 06:20:35Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#include <assert.h>
#include <string.h>
#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/Discretisation/Discretisation.h>
#include <StgFEM/SLE/SLE.h>

#include "types.h"
#include "PressureGradForceTerm.h"

/* Textual name of this class */
const Type PressureGradForceTerm_Type = "PressureGradForceTerm";

PressureGradForceTerm* PressureGradForceTerm_New( 
	Name							name,
	FiniteElementContext*	context,
	ForceVector*				forceVector,
	Swarm*						integrationSwarm,
	FeVariable*					gradField, 
	FeVariable*					pressureField )
{
	PressureGradForceTerm* self = (PressureGradForceTerm*) _PressureGradForceTerm_DefaultNew( name );

	self->isConstructed = True;
	_ForceTerm_Init( self, context, forceVector, integrationSwarm, NULL );
	_PressureGradForceTerm_Init( self, pressureField, gradField );

	return self;
}

/* Creation implementation / Virtual constructor */
PressureGradForceTerm* _PressureGradForceTerm_New(  PRESSUREGRADFORCETERM_DEFARGS  )
{
	PressureGradForceTerm* self;
	
	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(PressureGradForceTerm) );
	/* The following terms are parameters that have been passed into this function but are being set before being passed onto the parent */
	/* This means that any values of these parameters that are passed into this function are not passed onto the parent function
	   and so should be set to ZERO in any children of this class. */
	nameAllocationType = NON_GLOBAL;

	self = (PressureGradForceTerm*) _ForceTerm_New(  FORCETERM_PASSARGS  );
	
	/* Virtual info */
	
	return self;
}

void _PressureGradForceTerm_Init( void* forceTerm, FeVariable* pressureField, FeVariable* gradField ) {
	PressureGradForceTerm* self = (PressureGradForceTerm*)forceTerm;

	self->asmb = Assembler_New();
	Assembler_SetCallbacks( self->asmb, 
		NULL, 
		NULL, 
		(Assembler_CallbackType*)PressureGradForceTerm_RowCB, 
		(Assembler_CallbackType*)PressureGradForceTerm_ColCB, 
		(Assembler_CallbackType*)PressureGradForceTerm_ColCB, 
		self );
	self->pressureField = pressureField;
	self->gradField = gradField;
	self->forceVec = NULL;
	self->elForceVec = NULL;
	self->factor = 0.0;
}

void _PressureGradForceTerm_Delete( void* forceTerm ) {
	PressureGradForceTerm* self = (PressureGradForceTerm*)forceTerm;

	_ForceTerm_Delete( self );
}

void _PressureGradForceTerm_Print( void* forceTerm, Stream* stream ) {
	PressureGradForceTerm* self = (PressureGradForceTerm*)forceTerm;
	
	_ForceTerm_Print( self, stream );

	/* General info */
	Journal_PrintPointer( stream, self->pressureField );
}

void* _PressureGradForceTerm_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                              _sizeOfSelf = sizeof(PressureGradForceTerm);
	Type                                                      type = PressureGradForceTerm_Type;
	Stg_Class_DeleteFunction*                              _delete = _PressureGradForceTerm_Delete;
	Stg_Class_PrintFunction*                                _print = _PressureGradForceTerm_Print;
	Stg_Class_CopyFunction*                                  _copy = NULL;
	Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _PressureGradForceTerm_DefaultNew;
	Stg_Component_ConstructFunction*                    _construct = _PressureGradForceTerm_AssignFromXML;
	Stg_Component_BuildFunction*                            _build = _PressureGradForceTerm_Build;
	Stg_Component_InitialiseFunction*                  _initialise = _PressureGradForceTerm_Initialise;
	Stg_Component_ExecuteFunction*                        _execute = _PressureGradForceTerm_Execute;
	Stg_Component_DestroyFunction*                        _destroy = _PressureGradForceTerm_Destroy;
	ForceTerm_AssembleElementFunction*            _assembleElement = _PressureGradForceTerm_AssembleElement;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = ZERO;

	return (void*)_PressureGradForceTerm_New(  PRESSUREGRADFORCETERM_PASSARGS  );
}

void _PressureGradForceTerm_AssignFromXML( void* forceTerm, Stg_ComponentFactory* cf, void* data ) {
	PressureGradForceTerm*            self             = (PressureGradForceTerm*)forceTerm;
	FeVariable*                 pressureField;
	FeVariable*                 gradField;

	/* Construct Parent */
	_ForceTerm_AssignFromXML( self, cf, data );

	pressureField = Stg_ComponentFactory_ConstructByKey( cf, self->name, "PressureField", FeVariable, True, data ) ;
	gradField = Stg_ComponentFactory_ConstructByKey( cf, self->name, "PressureGradField", FeVariable, True, data ) ;

	_PressureGradForceTerm_Init( self, pressureField, gradField );
}

void _PressureGradForceTerm_Build( void* forceTerm, void* data ) {
	PressureGradForceTerm*             self             = (PressureGradForceTerm*)forceTerm;

	_ForceTerm_Build( self, data );
	Stg_Component_Build( self->pressureField, data, False );
	Stg_Component_Build( self->gradField, data, False );
	Assembler_SetVariables( self->asmb, self->gradField, self->pressureField );
	Assembler_SetIntegrationSwarm( self->asmb, self->integrationSwarm );
}

void _PressureGradForceTerm_Initialise( void* forceTerm, void* data ) {
	PressureGradForceTerm*             self             = (PressureGradForceTerm*)forceTerm;

	_ForceTerm_Initialise( self, data );
	Stg_Component_Initialise( self->pressureField, data, False );
}

void _PressureGradForceTerm_Execute( void* forceTerm, void* data ) {
	_ForceTerm_Execute( forceTerm, data );
}

void _PressureGradForceTerm_Destroy( void* forceTerm, void* data ) {
	_ForceTerm_Destroy( forceTerm, data );
}


void _PressureGradForceTerm_AssembleElement( void* forceTerm, ForceVector* forceVector, Element_LocalIndex lElement_I, double* elForceVec ) {
	PressureGradForceTerm* self = (PressureGradForceTerm*)forceTerm;

	assert( self );

	self->forceVec = forceVector;
	self->elForceVec = elForceVec;
	Assembler_IntegrateMatrixElement( self->asmb, lElement_I );
}

Bool PressureGradForceTerm_RowCB( PressureGradForceTerm* self, Assembler* assm ) {
	self->factor = assm->particle->weight * assm->detJac;
	return True;
}

Bool PressureGradForceTerm_ColCB( PressureGradForceTerm* self, Assembler* assm ) {
	double val;

	FeVariable_GetValueAtNode( self->pressureField, assm->colNodeInd, &val );
	self->elForceVec[assm->rowInd] += 
		val * 
		assm->shapeFuncs[assm->rowElNodeInd] * 
		assm->globalDerivs[assm->rowDofInd][assm->colElNodeInd] * 
		self->factor;
	return True;
}


