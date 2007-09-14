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
*+		Mirko Velic
*+		Julian Giordani
*+		Robert Turnbull
*+		Vincent Lemiale
*+		Louis Moresi
*+		David May
*+		David Stegman
*+		Patrick Sunter
** $Id: solHAy.c 636 2006-09-08 03:07:06Z JulianGiordani $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include "Analytic_solHAy.h"

const Type Velic_solHAy_Type = "Underworld_Velic_solHAy";

void Velic_solHAy_PressureFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* pressure ) {
	Velic_solHAy* self = (Velic_solHAy*) analyticSolution;
	
	_Velic_solHAy( coord, self->sigma, self->eta, self->dx, self->dy, self->x0, self->y0, NULL, pressure, NULL, NULL );
}

void Velic_solHAy_VelocityFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* velocity ) {
	Velic_solHAy* self = (Velic_solHAy*) analyticSolution;
	
	_Velic_solHAy( coord, self->sigma, self->eta, self->dx, self->dy, self->x0, self->y0, velocity, NULL, NULL, NULL );
}

void Velic_solHAy_StressFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* stress ) {
	Velic_solHAy* self = (Velic_solHAy*) analyticSolution;
	
	_Velic_solHAy( coord, self->sigma, self->eta, self->dx, self->dy, self->x0, self->y0, NULL, NULL, stress, NULL );
}


void Velic_solHAy_StrainRateFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* strainRate ) {
	Velic_solHAy* self = (Velic_solHAy*) analyticSolution;
	
	_Velic_solHAy( coord, self->sigma, self->eta, self->dx, self->dy, self->x0, self->y0, NULL, NULL, NULL, strainRate );
}

void _Velic_solHAy_Init( Velic_solHAy* self, double sigma, double eta, double dx, double dy, double x0, double y0 ) {
	// TODO	Bool                     correctInput = True;
	// Set the general AnalyticSolution (parent class) field function pointers to the specific analytic sols to use
        self->_getAnalyticVelocity    = Velic_solHAy_VelocityFunction;
	self->_getAnalyticPressure    = Velic_solHAy_PressureFunction;
	self->_getAnalyticTotalStress = Velic_solHAy_StressFunction;
	self->_getAnalyticStrainRate  = Velic_solHAy_StrainRateFunction;

	self->sigma = sigma;
	self->eta = eta;
	self->dx = dx;
	self->dy = dy;
	self->x0 = x0;
	self->y0 = y0;
}

void _Velic_solHAy_Construct( void* analyticSolution, Stg_ComponentFactory* cf, void* data ) {
	Velic_solHAy* self = (Velic_solHAy*) analyticSolution;
	FeVariable*              velocityField;
	FeVariable*              pressureField;
	FeVariable*              stressField;
	FeVariable*              strainRateField;
	FeVariable*              recoverdStrainRateField;
	FeVariable*              recoveredStressField;
	double                   sigma, eta, dx, dy, x0, y0;

	/* Construct Parent */
	_AnalyticSolution_Construct( self, cf, data );

	/* Create Analytic Fields */
	velocityField = Stg_ComponentFactory_ConstructByName( cf, "VelocityField", FeVariable, True, data );
	AnalyticSolution_RegisterFeVariableWithAnalyticFunction( self, velocityField, Velic_solHAy_VelocityFunction );

	pressureField = Stg_ComponentFactory_ConstructByName( cf, "PressureField", FeVariable, True, data );
	AnalyticSolution_RegisterFeVariableWithAnalyticFunction( self, pressureField, Velic_solHAy_PressureFunction );

	stressField = Stg_ComponentFactory_ConstructByName( cf, "StressField", FeVariable, False, data );
	if ( stressField )
		AnalyticSolution_RegisterFeVariableWithAnalyticFunction( self, stressField, Velic_solHAy_StressFunction );

	strainRateField = Stg_ComponentFactory_ConstructByName( cf, "StrainRateField", FeVariable, False, data );
	if ( strainRateField  ) {
		AnalyticSolution_RegisterFeVariableWithAnalyticFunction( self, strainRateField, Velic_solHAy_StrainRateFunction );
	}

	recoverdStrainRateField = Stg_ComponentFactory_ConstructByName( cf, "recoveredStrainRate", FeVariable, False, data );
	if ( recoverdStrainRateField )
		AnalyticSolution_RegisterFeVariableWithAnalyticFunction( self, recoverdStrainRateField, Velic_solHAy_StrainRateFunction );

	recoveredStressField = Stg_ComponentFactory_ConstructByName( cf, "recoveredStressField", FeVariable, False, data );
	if ( recoveredStressField )
		AnalyticSolution_RegisterFeVariableWithAnalyticFunction( self, recoveredStressField, Velic_solHAy_StressFunction );

	sigma = Stg_ComponentFactory_GetRootDictDouble( cf, "solHA_sigma", 1.0 );
	eta = Stg_ComponentFactory_GetRootDictDouble( cf, "solHA_eta", 1.0 );
	dx = Stg_ComponentFactory_GetRootDictDouble( cf, "solHA_dx", 0.4 );
	dy = Stg_ComponentFactory_GetRootDictDouble( cf, "solHA_dy", 0.35 );
	x0 = Stg_ComponentFactory_GetRootDictDouble( cf, "solHA_x0", 0.3 );
	y0 = Stg_ComponentFactory_GetRootDictDouble( cf, "solHA_y0", 0.6 );

	_Velic_solHAy_Init( self, sigma, eta, dx, dy, x0, y0 );
}

void* _Velic_solHAy_DefaultNew( Name name ) {
	return _AnalyticSolution_New(
			sizeof(Velic_solHAy),
			Velic_solHAy_Type,
			_AnalyticSolution_Delete,
			_AnalyticSolution_Print,
			_AnalyticSolution_Copy,
			_Velic_solHAy_DefaultNew,
			_Velic_solHAy_Construct,
			_AnalyticSolution_Build,
			_AnalyticSolution_Initialise,
			_AnalyticSolution_Execute,
			_AnalyticSolution_Destroy,
			name );
}

Index Underworld_Velic_solHAy_Register( PluginsManager* pluginsManager ) {
	return PluginsManager_Submit( pluginsManager, Velic_solHAy_Type, "0", _Velic_solHAy_DefaultNew );
}
