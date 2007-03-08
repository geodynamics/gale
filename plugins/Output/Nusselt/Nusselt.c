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
** $Id: Nusselt.c 358 2006-10-18 06:17:30Z SteveQuenette $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include <Underworld/Underworld.h>
#include <StgFEM/FrequentOutput/FrequentOutput.h>

#include "Nusselt.h"

const Type Underworld_Nusselt_Type = "Underworld_Nusselt";

void _Underworld_Nusselt_Construct( void* component, Stg_ComponentFactory* cf, void* data ) {
	UnderworldContext* context;

	context = Stg_ComponentFactory_ConstructByName( cf, "context", UnderworldContext, True, data ); 

	/* Add functions to entry points */
	Underworld_Nusselt_Setup( context );
	ContextEP_Append( context, AbstractContext_EP_FrequentOutput, Underworld_Nusselt_Output );
}

void* _Underworld_Nusselt_DefaultNew( Name name ) {
	return _Codelet_New(
		sizeof(Underworld_Nusselt),
		Underworld_Nusselt_Type,
		_Codelet_Delete,
		_Codelet_Print,
		_Codelet_Copy,
		_Underworld_Nusselt_DefaultNew,
		_Underworld_Nusselt_Construct,
		_Codelet_Build,
		_Codelet_Initialise,
		_Codelet_Execute,
		_Codelet_Destroy,
		name );
}

Index Underworld_Nusselt_Register( PluginsManager* pluginsManager ) {
	return PluginsManager_Submit( pluginsManager, Underworld_Nusselt_Type, "0", _Underworld_Nusselt_DefaultNew );
}

void Underworld_Nusselt_Setup( UnderworldContext* context ) {
	FieldVariable_Register*              fV_Register               = context->fieldVariable_Register;
	FieldVariable*                       temperatureGradientsField;
	FieldVariable*                       velocityField;
	FieldVariable*                       temperatureField;
	OperatorFeVariable*                  advectiveHeatFluxField;
	OperatorFeVariable*                  temperatureTotalDerivField;

	Underworld_Nusselt* self;

	self = (Underworld_Nusselt*)LiveComponentRegister_Get(
					context->CF->LCRegister,
					Underworld_Nusselt_Type );
	
	StgFEM_FrequentOutput_PrintString( context, "Nusselt" );

	Journal_Firewall( 
			context->gaussSwarm != NULL, 
			Underworld_Error,
			"Cannot find gauss swarm. Cannot use %s.\n", CURR_MODULE_NAME );

	/* Create Some FeVariables to calculate nusselt number */
	temperatureField          = FieldVariable_Register_GetByName( fV_Register, "TemperatureField" );
	velocityField             = FieldVariable_Register_GetByName( fV_Register, "VelocityField" );
	temperatureGradientsField = FieldVariable_Register_GetByName( fV_Register, "TemperatureGradientsField" );
	
	advectiveHeatFluxField = OperatorFeVariable_NewBinary(  
			"AdvectiveHeatFluxField",
			temperatureField, 
			velocityField, 
			"VectorScale" );

	temperatureTotalDerivField = OperatorFeVariable_NewBinary(  
			"TemperatureTotalDerivField",
			advectiveHeatFluxField, 
			temperatureGradientsField, 
			"Subtraction" );
	
	self->temperatureVertDerivField = (FeVariable*) OperatorFeVariable_NewUnary(  
			"VerticalAdvectiveHeatFluxField",
			temperatureTotalDerivField, 
			"TakeSecondComponent" );
	self->temperatureVertDerivField->feMesh = ((FeVariable*)velocityField)->feMesh;

	/* Add the variables to register so we can checkpoint & examine if necessary */
	FieldVariable_Register_Add( fV_Register, advectiveHeatFluxField );
	FieldVariable_Register_Add( fV_Register, temperatureTotalDerivField );
	FieldVariable_Register_Add( fV_Register, self->temperatureVertDerivField );

}

void Underworld_Nusselt_Output( UnderworldContext* context ) {

	Underworld_Nusselt* self;

	self = (Underworld_Nusselt*)LiveComponentRegister_Get(
				context->CF->LCRegister,
				Underworld_Nusselt_Type );

	/* This performs in integral given in the first half of equation 23 in 
	 *  Moresi, L. N. and Solomatov, V. S.,  
	 *  Numerical investigation of 2d convection with extremely large viscosity variations,  
	 *  Phys. Fluids, 1995, Volume 7, pp 2154-2162. */
	self->nusseltNumber = FeVariable_AverageTopLayer( self->temperatureVertDerivField, context->gaussSwarm, J_AXIS );

	StgFEM_FrequentOutput_PrintValue( context, self->nusseltNumber );
}
