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
** $Id: Nusselt.c 619 2007-10-29 11:43:40Z RobertTurnbull $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include <Underworld/Underworld.h>
#include <StgFEM/FrequentOutput/FrequentOutput.h>

#include "Nusselt.h"

const Type Underworld_Nusselt_Type = "Underworld_Nusselt";

void _Underworld_Nusselt_AssignFromXML( void* component, Stg_ComponentFactory* cf, void* data ) {
	Underworld_Nusselt*		self = (Underworld_Nusselt*)component;
	UnderworldContext*		context;
	FieldVariable_Register*	fV_Register;
	FieldVariable*				temperatureGradientsField;
	FieldVariable*				velocityField;
	FieldVariable*				temperatureField;

	self->context = (AbstractContext*)Stg_ComponentFactory_PluginConstructByKey( cf, self, "Context", UnderworldContext, True, data );
	self->gaussSwarm = Stg_ComponentFactory_PluginConstructByKey( cf, self, "GaussSwarm", Swarm, True, data );

	StgFEM_FrequentOutput_PrintString( self->context, "Nusselt" );

	context = (UnderworldContext*)self->context;
	fV_Register = context->fieldVariable_Register;
	/* Create Some FeVariables to calculate nusselt number */
	temperatureField          = FieldVariable_Register_GetByName( fV_Register, "TemperatureField" );
	velocityField             = FieldVariable_Register_GetByName( fV_Register, "VelocityField" );
	temperatureGradientsField = FieldVariable_Register_GetByName( fV_Register, "TemperatureGradientsField" );
	
	self->advectiveHeatFluxField = OperatorFeVariable_NewBinary(  
		"AdvectiveHeatFluxField",
		(DomainContext*) context,
		temperatureField, 
		velocityField, 
		"VectorScale" );

	self->temperatureTotalDerivField = OperatorFeVariable_NewBinary(  
		"TemperatureTotalDerivField",
		(DomainContext*) context,
		self->advectiveHeatFluxField, 
		temperatureGradientsField, 
		"Subtraction" );
	
	self->temperatureVertDerivField = (FeVariable*) OperatorFeVariable_NewUnary(  
		"VerticalAdvectiveHeatFluxField",
		(DomainContext*) context,
		self->temperatureTotalDerivField, 
		"TakeSecondComponent" );

	/* Add the variables to register so we can checkpoint & examine if necessary */
	FieldVariable_Register_Add( fV_Register, self->advectiveHeatFluxField );
	FieldVariable_Register_Add( fV_Register, self->temperatureTotalDerivField );
	FieldVariable_Register_Add( fV_Register, self->temperatureVertDerivField );

	/* Add functions to entry points */
	ContextEP_Append( self->context, AbstractContext_EP_FrequentOutput, Underworld_Nusselt_Output );
}

void _Underworld_Nusselt_Build( void* component, void* data ) {
	Underworld_Nusselt*	self = (Underworld_Nusselt*)component;

   Stg_Component_Build( self->gaussSwarm, data, False );
	Stg_Component_Build( self->advectiveHeatFluxField, data, False );
	Stg_Component_Build( self->temperatureTotalDerivField, data, False );
	Stg_Component_Build( self->temperatureVertDerivField, data, False );
   
   _Codelet_Build( component, data );
}

void _Underworld_Nusselt_Initialise( void* component, void* data ) {
	Underworld_Nusselt*	self = (Underworld_Nusselt*)component;

   Stg_Component_Initialise( self->gaussSwarm, data, False );
	Stg_Component_Initialise( self->advectiveHeatFluxField, data, False );
	Stg_Component_Initialise( self->temperatureTotalDerivField, data, False );
	Stg_Component_Initialise( self->temperatureVertDerivField, data, False );
   
   _Codelet_Initialise( component, data );
}

void _Underworld_Nusselt_Destroy( void* component, void* data ) {
	Underworld_Nusselt*	self = (Underworld_Nusselt*)component;

   Stg_Component_Destroy( self->gaussSwarm, data, False );
	Stg_Component_Destroy( self->advectiveHeatFluxField, data, False );
	Stg_Component_Destroy( self->temperatureTotalDerivField, data, False );
	Stg_Component_Destroy( self->temperatureVertDerivField, data, False );
   
   _Codelet_Destroy( component, data );
}

void* _Underworld_Nusselt_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                              _sizeOfSelf = sizeof(Underworld_Nusselt);
	Type                                                      type = Underworld_Nusselt_Type;
	Stg_Class_DeleteFunction*                              _delete = _Codelet_Delete;
	Stg_Class_PrintFunction*                                _print = _Codelet_Print;
	Stg_Class_CopyFunction*                                  _copy = _Codelet_Copy;
	Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _Underworld_Nusselt_DefaultNew;
	Stg_Component_ConstructFunction*                    _construct = _Underworld_Nusselt_AssignFromXML;
	Stg_Component_BuildFunction*                            _build = _Underworld_Nusselt_Build;
	Stg_Component_InitialiseFunction*                  _initialise = _Underworld_Nusselt_Initialise;
	Stg_Component_ExecuteFunction*                        _execute = _Codelet_Execute;
	Stg_Component_DestroyFunction*                        _destroy = _Underworld_Nusselt_Destroy;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = ZERO;

	return _Codelet_New(  CODELET_PASSARGS  );
}

Index Underworld_Nusselt_Register( PluginsManager* pluginsManager ) {
	return PluginsManager_Submit( pluginsManager, Underworld_Nusselt_Type, "0", _Underworld_Nusselt_DefaultNew );
}

void Underworld_Nusselt_Output( UnderworldContext* context ) {
	Underworld_Nusselt* self;

	self = (Underworld_Nusselt*)LiveComponentRegister_Get( context->CF->LCRegister, Underworld_Nusselt_Type );

	/* This performs in integral given in the first half of equation 23 in 
	 *  Moresi, L. N. and Solomatov, V. S.,  
	 *  Numerical investigation of 2d convection with extremely large viscosity variations,  
	 *  Phys. Fluids, 1995, Volume 7, pp 2154-2162. */
	self->nusseltNumber = FeVariable_AverageTopLayer( self->temperatureVertDerivField, self->gaussSwarm, J_AXIS );

	StgFEM_FrequentOutput_PrintValue( context, self->nusseltNumber );
}


