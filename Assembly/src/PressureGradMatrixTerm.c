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
#include <StgFEM/Discretisation/Discretisation.h>
#include <StgFEM/SLE/SLE.h>

#include "types.h"
#include "PressureGradMatrixTerm.h"

/* Textual name of this class */
const Type PressureGradMatrixTerm_Type = "PressureGradMatrixTerm";

PressureGradMatrixTerm* PressureGradMatrixTerm_New( 
		Name                                                name,
		StiffnessMatrix*                                    stiffMat,
		Swarm*                                              integrationSwarm,
		FeVariable*                                         gradField )
{
	PressureGradMatrixTerm* self = (PressureGradMatrixTerm*) _PressureGradMatrixTerm_DefaultNew( name );

	PressureGradMatrixTerm_InitAll( 
			self,
			stiffMat,
			integrationSwarm,
			gradField );

	return self;
}

/* Creation implementation / Virtual constructor */
PressureGradMatrixTerm* _PressureGradMatrixTerm_New( 
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
		StiffnessMatrixTerm_AssembleElementFunction*        _assembleElement, 
		Name                                                name )
{
	PressureGradMatrixTerm* self;
	
	/* Allocate memory */
	assert( sizeOfSelf >= sizeof(PressureGradMatrixTerm) );
	self = (PressureGradMatrixTerm*) _StiffnessMatrixTerm_New( 
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
		_assembleElement,
		name );
	
	/* Virtual info */
	
	return self;
}

void _PressureGradMatrixTerm_Init( 
		PressureGradMatrixTerm*                             self, 
		FeVariable*                                         gradField )
{
	self->assm = Assembler_New();
	Assembler_SetCallbacks( self->assm, 
				NULL, 
				NULL, 
				(Assembler_CallbackType*)PressureGradMatrixTerm_RowCB, 
				(Assembler_CallbackType*)PressureGradMatrixTerm_ColCB, 
				(Assembler_CallbackType*)PressureGradMatrixTerm_ColCB, 
				self );
	self->gradField = gradField;
	self->stiffMat = NULL;
	self->elStiffMat = NULL;
	self->factor = 0.0;
}

void PressureGradMatrixTerm_InitAll( 
		void*                                               matrixTerm,
		StiffnessMatrix*                                        stiffMat,
		Swarm*                                              integrationSwarm,
		FeVariable*                                         gradField )
{
	PressureGradMatrixTerm* self = (PressureGradMatrixTerm*) matrixTerm;

	StiffnessMatrixTerm_InitAll( self, stiffMat, integrationSwarm, NULL );
	_PressureGradMatrixTerm_Init( self, gradField );
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
	return (void*)_PressureGradMatrixTerm_New( 
		sizeof(PressureGradMatrixTerm), 
		PressureGradMatrixTerm_Type,
		_PressureGradMatrixTerm_Delete,
		_PressureGradMatrixTerm_Print,
		NULL,
		_PressureGradMatrixTerm_DefaultNew,
		_PressureGradMatrixTerm_Construct,
		_PressureGradMatrixTerm_Build,
		_PressureGradMatrixTerm_Initialise,
		_PressureGradMatrixTerm_Execute,
		_PressureGradMatrixTerm_Destroy,
		_PressureGradMatrixTerm_AssembleElement,
		name );
}

void _PressureGradMatrixTerm_Construct( void* matrixTerm, Stg_ComponentFactory* cf, void* data ) {
	PressureGradMatrixTerm*            self             = (PressureGradMatrixTerm*)matrixTerm;
	FeVariable*                 gradField;

	/* Construct Parent */
	_StiffnessMatrixTerm_Construct( self, cf, data );

	gradField = Stg_ComponentFactory_ConstructByKey( cf, self->name, "PressureGradField", FeVariable, True, data ) ;

	_PressureGradMatrixTerm_Init( self, gradField );
}

void _PressureGradMatrixTerm_Build( void* matrixTerm, void* data ) {
	PressureGradMatrixTerm*             self             = (PressureGradMatrixTerm*)matrixTerm;

	_StiffnessMatrixTerm_Build( self, data );
	Stg_Component_Build( self->gradField, data, False );
	Assembler_SetVariables( self->assm, self->gradField, self->gradField );
	Assembler_SetIntegrationSwarm( self->assm, self->integrationSwarm );
}

void _PressureGradMatrixTerm_Initialise( void* matrixTerm, void* data ) {
	PressureGradMatrixTerm*             self             = (PressureGradMatrixTerm*)matrixTerm;

	_StiffnessMatrixTerm_Initialise( self, data );
	Stg_Component_Initialise( self->gradField, data, False );
}

void _PressureGradMatrixTerm_Execute( void* matrixTerm, void* data ) {
	_StiffnessMatrixTerm_Execute( matrixTerm, data );
}

void _PressureGradMatrixTerm_Destroy( void* matrixTerm, void* data ) {
	_StiffnessMatrixTerm_Destroy( matrixTerm, data );
}


void _PressureGradMatrixTerm_AssembleElement( void* matrixTerm,
					      StiffnessMatrix* stiffMat, 
					      Element_LocalIndex lElement_I, 
					      SystemLinearEquations* sle, 
					      FiniteElementContext* context, 
					      double** elStiffMat )
{
	PressureGradMatrixTerm* self = (PressureGradMatrixTerm*)matrixTerm;

	assert( self );

	self->stiffMat = stiffMat;
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
