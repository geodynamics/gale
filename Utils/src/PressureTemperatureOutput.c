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
** $Id: PressureTemperatureOutput.c 610 2007-10-11 08:09:29Z SteveQuenette $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>

#include "types.h"
#include "PressureTemperatureOutput.h"

#include <assert.h>
#include <string.h>

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type PressureTemperatureOutput_Type = "PressureTemperatureOutput";

/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/
PressureTemperatureOutput* PressureTemperatureOutput_New(
		Name                                  name,
		TimeIntegrator*                       timeIntegrator,
		FeVariable*                           velocityField )
{
	PressureTemperatureOutput* self = (PressureTemperatureOutput*) _PressureTemperatureOutput_DefaultNew( name );

	/* 	PressureTemperatureOutput_InitAll */
	abort();

	return self;
}

PressureTemperatureOutput* _PressureTemperatureOutput_New(  PRESSURETEMPERATUREOUTPUT_DEFARGS  )
{
	PressureTemperatureOutput* self;
	
	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( _sizeOfSelf >= sizeof(PressureTemperatureOutput) );
	self = (PressureTemperatureOutput*)_SwarmOutput_New(  SWARMOUTPUT_PASSARGS  );

	
	/* General info */

	/* Virtual Info */
	
	return self;
}

void _PressureTemperatureOutput_Init( 
		void*                                 swarmOutput,
		FeVariable*                           pressureField,
		FeVariable*                           temperatureField )
{
	PressureTemperatureOutput*   self                = (PressureTemperatureOutput*)swarmOutput;

	self->pressureField      = pressureField;
	self->temperatureField   = temperatureField;
}


/*------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _PressureTemperatureOutput_Delete( void* swarmOutput ) {
	PressureTemperatureOutput* self = (PressureTemperatureOutput*)swarmOutput;

	/* Delete parent */
	_SwarmOutput_Delete( self );
}


void _PressureTemperatureOutput_Print( void* swarmOutput, Stream* stream ) {
	PressureTemperatureOutput* self = (PressureTemperatureOutput*)swarmOutput;
	
	/* Print parent */
	_SwarmOutput_Print( self, stream );
}

void* _PressureTemperatureOutput_Copy( void* swarmOutput, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	PressureTemperatureOutput*	self = (PressureTemperatureOutput*)swarmOutput;
	PressureTemperatureOutput*	newPressureTemperatureOutput;
	
	newPressureTemperatureOutput = (PressureTemperatureOutput*)_SwarmOutput_Copy( self, dest, deep, nameExt, ptrMap );
	
	return (void*)newPressureTemperatureOutput;
}

void* _PressureTemperatureOutput_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                              _sizeOfSelf = sizeof(PressureTemperatureOutput);
	Type                                                      type = PressureTemperatureOutput_Type;
	Stg_Class_DeleteFunction*                              _delete = _PressureTemperatureOutput_Delete;
	Stg_Class_PrintFunction*                                _print = _PressureTemperatureOutput_Print;
	Stg_Class_CopyFunction*                                  _copy = _PressureTemperatureOutput_Copy;
	Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _PressureTemperatureOutput_DefaultNew;
	Stg_Component_ConstructFunction*                    _construct = _PressureTemperatureOutput_AssignFromXML;
	Stg_Component_BuildFunction*                            _build = _PressureTemperatureOutput_Build;
	Stg_Component_InitialiseFunction*                  _initialise = _PressureTemperatureOutput_Initialise;
	Stg_Component_ExecuteFunction*                        _execute = _PressureTemperatureOutput_Execute;
	Stg_Component_DestroyFunction*                        _destroy = _PressureTemperatureOutput_Destroy;
	SwarmOutput_PrintHeaderFunction*                  _printHeader = _PressureTemperatureOutput_PrintHeader;
	SwarmOutput_PrintDataFunction*                      _printData = _PressureTemperatureOutput_PrintData;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = ZERO;

	return (void*) _PressureTemperatureOutput_New(  PRESSURETEMPERATUREOUTPUT_PASSARGS  );
}


void _PressureTemperatureOutput_AssignFromXML( void* swarmOutput, Stg_ComponentFactory* cf, void* data ) {
	PressureTemperatureOutput*  self          = (PressureTemperatureOutput*) swarmOutput;
	FeVariable*                 pressureField;
	FeVariable*                 temperatureField;

	_SwarmOutput_AssignFromXML( self, cf, data );

	pressureField    = Stg_ComponentFactory_ConstructByKey( cf, self->name, "PressureField",    FeVariable, True, data ) ;
	temperatureField = Stg_ComponentFactory_ConstructByKey( cf, self->name, "TemperatureField", FeVariable, True, data ) ;

	_PressureTemperatureOutput_Init( self, pressureField, temperatureField );

}

void _PressureTemperatureOutput_Build( void* swarmOutput, void* data ) {
	PressureTemperatureOutput*	self = (PressureTemperatureOutput*) swarmOutput;

	Stg_Component_Build( self->pressureField, data, False );
	Stg_Component_Build( self->temperatureField, data, False );

	_SwarmOutput_Build( self, data );
}
void _PressureTemperatureOutput_Initialise( void* swarmOutput, void* data ) {
	PressureTemperatureOutput*	self = (PressureTemperatureOutput*) swarmOutput;

	Stg_Component_Initialise( self->pressureField, data, False );
	Stg_Component_Initialise( self->temperatureField, data, False );
	
	_SwarmOutput_Initialise( self, data );
}
void _PressureTemperatureOutput_Execute( void* swarmOutput, void* data ) {
	PressureTemperatureOutput*	self = (PressureTemperatureOutput*)swarmOutput;
	
	_SwarmOutput_Execute( self, data );
}
void _PressureTemperatureOutput_Destroy( void* swarmOutput, void* data ) {
	PressureTemperatureOutput*	self = (PressureTemperatureOutput*)swarmOutput;

	Stg_Component_Destroy( self->pressureField, data, False );
	Stg_Component_Destroy( self->temperatureField, data, False );
	
	_SwarmOutput_Destroy( self, data );
}

void _PressureTemperatureOutput_PrintHeader( void* swarmOutput, Stream* stream, Particle_Index lParticle_I, void* context ){
	PressureTemperatureOutput*	self = (PressureTemperatureOutput*)swarmOutput;
	
	_SwarmOutput_PrintHeader( self, stream, lParticle_I, context );
	
	SwarmOutput_PrintString( self, stream, "Pressure" );
	SwarmOutput_PrintString( self, stream, "Temperature" );
}

void _PressureTemperatureOutput_PrintData( void* swarmOutput, Stream* stream, Particle_Index lParticle_I, void* context ){
	PressureTemperatureOutput*	self                = (PressureTemperatureOutput*)swarmOutput;
	Swarm*                      swarm               = self->swarm;
	GlobalParticle*           particle            = (GlobalParticle*)Swarm_ParticleAt( swarm, lParticle_I );
	double*                     coord               = particle->coord;
	double                      pressure;
	double                      temperature;

	Journal_Firewall(
		swarm->particleLayout->coordSystem == GlobalCoordSystem,
		Journal_MyStream( Error_Type, self ),
		"Swarm is not using global coord system! Modify this code to use both systems\n" );

	_SwarmOutput_PrintData( self, stream, lParticle_I, context );

	FieldVariable_InterpolateValueAt( self->pressureField,    coord, &pressure );
	FieldVariable_InterpolateValueAt( self->temperatureField, coord, &temperature );
	
	SwarmOutput_PrintValue( self, stream, pressure );
	SwarmOutput_PrintValue( self, stream, temperature );
}

/*-------------------------------------------------------------------------------------------------------------------------
** Private Functions
*/
/*---------------------------------------------------------------------------------------------------------------------
** Entry Point Hooks
*/

/*-------------------------------------------------------------------------------------------------------------------------
** Public Functions
*/



