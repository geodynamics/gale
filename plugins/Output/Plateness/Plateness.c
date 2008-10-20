/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
** Copyright (c) 2005, Monash Cluster Computing 
** All rights reserved.
** Redistribution and use in source and binary forms, with or without modification,
** are permitted provided that the following conditions are met:
**
** 		* Redistributions of source code must retain the above copyright notice, 
** 			this list of conditions and the following disclaimer.
** 		* Redistributions in binary form must reproduce the above copyright 
**			notice, this list of conditions and the following disclaimer in the 
**			documentation and/or other materials provided with the distribution.
** 		* Neither the name of the Monash University nor the names of its contributors 
**			may be used to endorse or promote products derived from this software 
**			without specific prior written permission.
**
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
** THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR 
** PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS 
** BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
** CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
** SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) 
** HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT 
** LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT 
** OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
**
**
** Contact:
*%		Louis Moresi - Louis.Moresi@sci.monash.edu.au
*%
** Contributors:
*+		Robert Turnbull
*+		Vincent Lemiale
*+		Louis Moresi
*+		David May
*+		David Stegman
*+		Mirko Velic
*+		Patrick Sunter
*+		Julian Giordani
*+
** $Id:  $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include <Underworld/Underworld.h>
#include <StgFEM/FrequentOutput/FrequentOutput.h>

#include "Plateness.h"

const Type Underworld_Plateness_Type = "Underworld_Plateness";

void _Underworld_Plateness_Construct( void* component, Stg_ComponentFactory* cf, void* data ) {
	UnderworldContext* context;

	context = Stg_ComponentFactory_ConstructByName( cf, "context", UnderworldContext, True, data ); 

	/* Add functions to entry points */
	Underworld_Plateness_Setup( context );
	ContextEP_Append( context, AbstractContext_EP_FrequentOutput, Underworld_Plateness_Output );
}

void _Underworld_Plateness_Build( void* component, void* data ) {
	Underworld_Plateness*	self = (Underworld_Plateness*)component;

	Stg_Component_Build( self->reducedStrainRateField, data, False );
	Stg_Component_Build( self->reducedStrainRateFieldInvariant, data, False );
	Stg_Component_Build( self->reducedStrainRateFieldInvariantRoot, data, False );

}

void* _Underworld_Plateness_DefaultNew( Name name ) {
	return _Codelet_New(
		sizeof(Underworld_Plateness),
		Underworld_Plateness_Type,
		_Codelet_Delete,
		_Codelet_Print,
		_Codelet_Copy,
		_Underworld_Plateness_DefaultNew,
		_Underworld_Plateness_Construct,
		_Underworld_Plateness_Build,
		_Codelet_Initialise,
		_Codelet_Execute,
		_Codelet_Destroy,
		name );
}

Index Underworld_Plateness_Register( PluginsManager* pluginsManager ) {
	return PluginsManager_Submit( pluginsManager, Underworld_Plateness_Type, "0", _Underworld_Plateness_DefaultNew );
}

void Underworld_Plateness_Setup( UnderworldContext* context ) {
	FieldVariable_Register*              fV_Register               = context->fieldVariable_Register;
	FieldVariable*                       strainRateField;
	Func_Ptr                             _carryOut;
	Dof_Index                            resultDofs;
	Index                                numberOfOperands;
	Operator*                            ownOperator;

	Underworld_Plateness* self;

	self = (Underworld_Plateness*)LiveComponentRegister_Get(
					context->CF->LCRegister,
					Underworld_Plateness_Type );
	
	StgFEM_FrequentOutput_PrintString( context, "Plateness" );

	Journal_Firewall( 
			context->gaussSwarm != NULL, 
			Underworld_Error,
			"Cannot find gauss swarm. Cannot use %s.\n", CURR_MODULE_NAME );

	/* Create Some FeVariables to determine plateness */
	strainRateField  = FieldVariable_Register_GetByName( fV_Register, "StrainRateField" );

	/* setup required parameters to create new operate */
	resultDofs       = ( strainRateField->dim == 2 ? 1 : 3 );
	numberOfOperands = 1;
	_carryOut        = ( strainRateField->dim == 2 ? Underworld_Plateness_SymmetricTensor_LowerDimension_2d : Underworld_Plateness_SymmetricTensor_LowerDimension_3d );
	
	ownOperator = Operator_New( "Plateness_SymmetricTensor_LowerDimension", _carryOut, numberOfOperands, strainRateField->fieldComponentCount, resultDofs, strainRateField->dim );

	self->reducedStrainRateField = OperatorFeVariable_NewUnary_OwnOperator(
			"ReducedStrainRateField",
			strainRateField, 
			ownOperator );
	
	self->reducedStrainRateFieldInvariant = OperatorFeVariable_NewUnary(  
			"ReducedStrainRateFieldInvariant",
			self->reducedStrainRateField, 
			"SymmetricTensor_Invariant" );

	self->reducedStrainRateFieldInvariantRoot = OperatorFeVariable_NewUnary(  
			"ReducedStrainRateFieldInvariantRoot",
			self->reducedStrainRateFieldInvariant, 
			"SquareRoot" );


	/* Add the variables to register so we can checkpoint & examine if necessary */
	FieldVariable_Register_Add( fV_Register, self->reducedStrainRateField );
	FieldVariable_Register_Add( fV_Register, self->reducedStrainRateFieldInvariant );
	FieldVariable_Register_Add( fV_Register, self->reducedStrainRateFieldInvariantRoot );
}

void Underworld_Plateness_Output( UnderworldContext* context ) {

	Underworld_Plateness* self;

	self = (Underworld_Plateness*)LiveComponentRegister_Get(
				context->CF->LCRegister,
				Underworld_Plateness_Type );

	self->plateness = FeVariable_AverageTopLayer( self->reducedStrainRateFieldInvariantRoot, context->gaussSwarm, J_AXIS );

	StgFEM_FrequentOutput_PrintValue( context, self->plateness );
}

void Underworld_Plateness_SymmetricTensor_LowerDimension_2d( void* operator, double* operand0, double* result ) {
	Operator* self = (Operator*) operator;
	
	Operator_FirewallUnary( self );
	Operator_FirewallResultDofs( self, self->operandDofs );

	result[0] = operand0[0];
}

void Underworld_Plateness_SymmetricTensor_LowerDimension_3d( void* operator, double* operand0, double* result ) {
	Operator* self = (Operator*) operator;
	
	Operator_FirewallUnary( self );
	Operator_FirewallResultDofs( self, self->operandDofs );

	result[0] = operand0[0];
	result[1] = operand0[1];
	result[2] = operand0[6];
}

