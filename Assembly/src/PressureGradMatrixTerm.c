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
** $Id: PressureGradMatrixTerm.c 822 2007-04-27 06:20:35Z LukeHodkinson $
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
#include "PressureGradMatrixTerm.h"

/* Textual name of this class */
const Type PressureGradMatrixTerm_Type = "PressureGradMatrixTerm";

PressureGradMatrixTerm* PressureGradMatrixTerm_New( 
	Name							name,
	FiniteElementContext*	context,
	StiffnessMatrix*			stiffnessMatrix,
	Swarm*						integrationSwarm,
	FeVariable*					gradField )
{
	PressureGradMatrixTerm* self = (PressureGradMatrixTerm*) _PressureGradMatrixTerm_DefaultNew( name );

	self->isConstructed = True;
	_StiffnessMatrixTerm_Init( self, context, stiffnessMatrix, integrationSwarm, NULL );
	_PressureGradMatrixTerm_Init( self, gradField );

	return self;
}

/* Creation implementation / Virtual constructor */
PressureGradMatrixTerm* _PressureGradMatrixTerm_New(  PRESSUREGRADMATRIXTERM_DEFARGS  )
{
	PressureGradMatrixTerm* self;
	
	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(PressureGradMatrixTerm) );
	self = (PressureGradMatrixTerm*) _StiffnessMatrixTerm_New(  STIFFNESSMATRIXTERM_PASSARGS  );
	
	/* Virtual info */
	
	return self;
}

void _PressureGradMatrixTerm_Init( void* matrixTerm, FeVariable* gradField ) {
	PressureGradMatrixTerm* self = (PressureGradMatrixTerm*)matrixTerm;

	self->assm = Assembler_New();
	Assembler_SetCallbacks( self->assm, 
		NULL, 
		NULL, 
		(Assembler_CallbackType*)PressureGradMatrixTerm_RowCB, 
		(Assembler_CallbackType*)PressureGradMatrixTerm_ColCB, 
		(Assembler_CallbackType*)PressureGradMatrixTerm_ColCB, 
		self );
	self->gradField = gradField;
	self->stiffnessMatrix = NULL;
	self->elStiffMat = NULL;
	self->factor = 0.0;
}

void _PressureGradMatrixTerm_Delete( void* matrixTerm ) {
	PressureGradMatrixTerm* self = (PressureGradMatrixTerm*)matrixTerm;

	_StiffnessMatrixTerm_Delete( self );
}

void _PressureGradMatrixTerm_Print( void* matrixTerm, Stream* stream ) {
	PressureGradMatrixTerm* self = (PressureGradMatrixTerm*)matrixTerm;
	
	_StiffnessMatrixTerm_Print( self, stream );

	/* General info */
	Journal_PrintPointer( stream, self->gradField );
}

void* _PressureGradMatrixTerm_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                                 _sizeOfSelf = sizeof(PressureGradMatrixTerm);
	Type                                                         type = PressureGradMatrixTerm_Type;
	Stg_Class_DeleteFunction*                                 _delete = _PressureGradMatrixTerm_Delete;
	Stg_Class_PrintFunction*                                   _print = _PressureGradMatrixTerm_Print;
	Stg_Class_CopyFunction*                                     _copy = NULL;
	Stg_Component_DefaultConstructorFunction*     _defaultConstructor = _PressureGradMatrixTerm_DefaultNew;
	Stg_Component_ConstructFunction*                       _construct = _PressureGradMatrixTerm_AssignFromXML;
	Stg_Component_BuildFunction*                               _build = _PressureGradMatrixTerm_Build;
	Stg_Component_InitialiseFunction*                     _initialise = _PressureGradMatrixTerm_Initialise;
	Stg_Component_ExecuteFunction*                           _execute = _PressureGradMatrixTerm_Execute;
	Stg_Component_DestroyFunction*                           _destroy = _PressureGradMatrixTerm_Destroy;
	StiffnessMatrixTerm_AssembleElementFunction*     _assembleElement = _PressureGradMatrixTerm_AssembleElement;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = ZERO;

	return (void*)_PressureGradMatrixTerm_New(  PRESSUREGRADMATRIXTERM_PASSARGS  );
}

void _PressureGradMatrixTerm_AssignFromXML( void* matrixTerm, Stg_ComponentFactory* cf, void* data ) {
	PressureGradMatrixTerm*	self = (PressureGradMatrixTerm*)matrixTerm;
	FeVariable*					gradField;

	/* Construct Parent */
	_StiffnessMatrixTerm_AssignFromXML( self, cf, data );

	gradField = Stg_ComponentFactory_ConstructByKey( cf, self->name, "PressureGradField", FeVariable, True, data ) ;

	_PressureGradMatrixTerm_Init( self, gradField );
}

void _PressureGradMatrixTerm_Build( void* matrixTerm, void* data ) {
	PressureGradMatrixTerm* self = (PressureGradMatrixTerm*)matrixTerm;

	_StiffnessMatrixTerm_Build( self, data );
	Stg_Component_Build( self->gradField, data, False );
	Assembler_SetVariables( self->assm, self->gradField, self->gradField );
	Assembler_SetIntegrationSwarm( self->assm, self->integrationSwarm );
}

void _PressureGradMatrixTerm_Initialise( void* matrixTerm, void* data ) {
	PressureGradMatrixTerm* self = (PressureGradMatrixTerm*)matrixTerm;

	_StiffnessMatrixTerm_Initialise( self, data );
	Stg_Component_Initialise( self->gradField, data, False );
}

void _PressureGradMatrixTerm_Execute( void* matrixTerm, void* data ) {
	_StiffnessMatrixTerm_Execute( matrixTerm, data );
}

void _PressureGradMatrixTerm_Destroy( void* matrixTerm, void* data ) {
	_StiffnessMatrixTerm_Destroy( matrixTerm, data );
}

void _PressureGradMatrixTerm_AssembleElement(
	 void*						matrixTerm,
	StiffnessMatrix*			stiffnessMatrix, 
	Element_LocalIndex		lElement_I, 
	SystemLinearEquations*	 sle, 
	FiniteElementContext*	context, 
	double**						elStiffMat )
{
	PressureGradMatrixTerm* self = (PressureGradMatrixTerm*)matrixTerm;

	assert( self );

	self->stiffnessMatrix = stiffnessMatrix;
	self->elStiffMat = elStiffMat;
	Assembler_IntegrateMatrixElement( self->assm, lElement_I );
}

Bool PressureGradMatrixTerm_RowCB( PressureGradMatrixTerm* self, Assembler* assm ) {
	self->factor = assm->particle->weight * assm->detJac;
	return True;
}

Bool PressureGradMatrixTerm_ColCB( PressureGradMatrixTerm* self, Assembler* assm ) {
	if( assm->rowDofInd != assm->colDofInd )
		return;

	self->elStiffMat[assm->rowInd][assm->colInd] += 
		assm->shapeFuncs[assm->rowElNodeInd] * 
		assm->shapeFuncs[assm->colElNodeInd] * 
		self->factor;
	return True;
}


