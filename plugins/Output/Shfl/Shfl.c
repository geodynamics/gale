/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
** Copyright (c) 2005, Your Mom
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
** THIS SOFTWARE IS PROVIDED BY YOUR MOM AND CONTRIBUTORS 
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
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include <Underworld/Underworld.h>
#include <StgFEM/FrequentOutput/FrequentOutput.h>

#include "Shfl.h"

const Type Underworld_Shfl_Type = "Underworld_Shfl";

void _Underworld_Shfl_AssignFromXML( void* component, Stg_ComponentFactory* cf, void* data ) {
	UnderworldContext* context;

	context = Stg_ComponentFactory_ConstructByName( cf, "context", UnderworldContext, True, data ); 

	/* Add functions to entry points */
	Underworld_Shfl_Setup( context );
	ContextEP_Append( context, AbstractContext_EP_FrequentOutput, Underworld_Shfl_Output );
}

void _Underworld_Shfl_Build( void* component, void* data ) {
	Underworld_Shfl*	self = (Underworld_Shfl*)component;

	assert( self );

	Stg_Component_Build( self->dTField, data, False );
   
   _Codelet_Build( self, data );
}

void _Underworld_Shfl_Destroy( void* component, void* data ) {
	Underworld_Shfl*	self = (Underworld_Shfl*)component;

   _Codelet_Destroy( self, data );

	Stg_Component_Destroy( self->dTField, data, False );
}

void _Underworld_Shfl_Initialise( void* component, void* data ) {
	Underworld_Shfl*	self = (Underworld_Shfl*)component;

	assert( self );

	Stg_Component_Initialise( self->dTField, data, False );
   
   _Codelet_Initialise( self, data );
}

void* _Underworld_Shfl_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                              _sizeOfSelf = sizeof(Underworld_Shfl);
	Type                                                      type = Underworld_Shfl_Type;
	Stg_Class_DeleteFunction*                              _delete = _Codelet_Delete;
	Stg_Class_PrintFunction*                                _print = _Codelet_Print;
	Stg_Class_CopyFunction*                                  _copy = _Codelet_Copy;
	Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _Underworld_Shfl_DefaultNew;
	Stg_Component_ConstructFunction*                    _construct = _Underworld_Shfl_AssignFromXML;
	Stg_Component_BuildFunction*                            _build = _Underworld_Shfl_Build;
	Stg_Component_InitialiseFunction*                  _initialise = _Underworld_Shfl_Initialise;
	Stg_Component_ExecuteFunction*                        _execute = _Codelet_Execute;
	Stg_Component_DestroyFunction*                        _destroy = _Underworld_Shfl_Destroy;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = ZERO;

	return _Codelet_New(  CODELET_PASSARGS  );
}

Index Underworld_Shfl_Register( PluginsManager* pluginsManager ) {
	return PluginsManager_Submit( pluginsManager, Underworld_Shfl_Type, "0", _Underworld_Shfl_DefaultNew );
}


void Underworld_Shfl_Setup( UnderworldContext* context ) {
	FieldVariable_Register*  fV_Register               = context->fieldVariable_Register;
	FieldVariable*           temperatureGradientsField;
	FieldVariable*           velocityField;
	FieldVariable*           temperatureField;
	OperatorFeVariable*      advectiveHeatFluxField;
	OperatorFeVariable*      temperatureTotalDerivField;
	Swarm*					    gaussSwarm = (Swarm*)LiveComponentRegister_Get( context->CF->LCRegister, "gaussSwarm" );

	Underworld_Shfl* self;

	self = (Underworld_Shfl*)LiveComponentRegister_Get(
					context->CF->LCRegister,
					Underworld_Shfl_Type );
	
	StgFEM_FrequentOutput_PrintString( context, "Shfl" );

	Journal_Firewall( 
			gaussSwarm != NULL, 
			Underworld_Error,
			"Cannot find gauss swarm. Cannot use %s.\n", CURR_MODULE_NAME );

	/* Create Some FeVariables to calculate shear heat flux  */
	temperatureField          = FieldVariable_Register_GetByName( fV_Register, "TemperatureField" );
	velocityField             = FieldVariable_Register_GetByName( fV_Register, "VelocityField" );
	temperatureGradientsField = FieldVariable_Register_GetByName( fV_Register, "TemperatureGradientsField" );
	
	/* get uT */
	advectiveHeatFluxField = OperatorFeVariable_NewBinary(  
			"AdvectiveHeatFluxField",
         (DomainContext*)	context,
			temperatureField, 
			velocityField, 
			"VectorScale" );
	/* get dTdz - uT */
	temperatureTotalDerivField = OperatorFeVariable_NewBinary(  
			"TemperatureTotalDerivField",
         (DomainContext*)	context,
			advectiveHeatFluxField, 
			temperatureGradientsField, 
			"Addition" );
	
	self->dTField = (FeVariable*) OperatorFeVariable_NewUnary(  
			"dTfieldz",
         (DomainContext*)	context,
			temperatureGradientsField, 
			"TakeSecondComponent" );
	self->dTField->feMesh = ((FeVariable*)velocityField)->feMesh;
	
}

void Underworld_Shfl_Output( UnderworldContext* context ) {

	Underworld_Shfl* self;
	FieldVariable_Register*              fV_Register               = context->fieldVariable_Register;
	FieldVariable*                       temperatureField;
	double								 temp;

	temperatureField          			 = FieldVariable_Register_GetByName( fV_Register, "TemperatureField" );
	
	self = (Underworld_Shfl*)LiveComponentRegister_Get(
				context->CF->LCRegister,
				Underworld_Shfl_Type );

	/* get_layer value */
	self->shfl = FeVariable_AveragePlane( self->dTField, J_AXIS, 1.0 );
	temp =   FeVariable_AveragePlane( temperatureField, J_AXIS, 0.0 );
	
		
	self->shfl/=temp;
		
	StgFEM_FrequentOutput_PrintValue( context, -self->shfl  );
}


