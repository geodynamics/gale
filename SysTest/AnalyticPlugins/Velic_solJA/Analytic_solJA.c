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
** $Id: solJA.c 636 2006-09-08 03:07:06Z JulianGiordani $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include "Analytic_solJA.h"

const Type Velic_solJA_Type = "Underworld_Velic_solJA";

void Velic_solJA_PressureFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* pressure ) {
	Velic_solJA* self = (Velic_solJA*) analyticSolution;
	
	_Velic_solJA( coord, self->sigma, self->etaA, self->etaB, self->dx, self->x0, self->zc,
		       	NULL, pressure, NULL, NULL );
}

void Velic_solJA_VelocityFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* velocity ) {
	Velic_solJA* self = (Velic_solJA*) analyticSolution;
	
	_Velic_solJA( coord, self->sigma, self->etaA, self->etaB, self->dx, self->x0, self->zc,
		       	velocity, NULL, NULL, NULL );
}

void Velic_solJA_StressFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* stress ) {
	Velic_solJA* self = (Velic_solJA*) analyticSolution;
	
	_Velic_solJA( coord, self->sigma, self->etaA, self->etaB, self->dx, self->x0, self->zc,
		       	NULL, NULL, stress, NULL );
}


void Velic_solJA_StrainRateFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* strainRate ) {
	Velic_solJA* self = (Velic_solJA*) analyticSolution;
	
	_Velic_solJA( coord, self->sigma, self->etaA, self->etaB, self->dx, self->x0, self->zc,
		       	NULL, NULL, NULL, strainRate );
}

void _Velic_solJA_Init( Velic_solJA* self, double sigma, double etaA, double etaB, double dx, double x0, double zc ) {
	// TODO	Bool                     correctInput = True;
	// Set the general AnalyticSolution (parent class) field function pointers to the specific analytic sols to use
        self->_getAnalyticVelocity    = Velic_solJA_VelocityFunction;
	self->_getAnalyticPressure    = Velic_solJA_PressureFunction;
	self->_getAnalyticTotalStress = Velic_solJA_StressFunction;
	self->_getAnalyticStrainRate  = Velic_solJA_StrainRateFunction;

	self->sigma = sigma;
	self->etaA = etaA;
	self->etaB = etaB;
	self->dx = dx;
	self->x0 = x0;
	self->zc = zc;
}

void _Velic_solJA_AssignFromXML( void* analyticSolution, Stg_ComponentFactory* cf, void* data ) {
	Velic_solJA* self = (Velic_solJA*) analyticSolution;
	FeVariable*              velocityField;
	FeVariable*              pressureField;
	FeVariable*              stressField;
	FeVariable*              strainRateField;
	FeVariable*              recoveredStrainRateField;
	FeVariable*              recoveredStressField;
	double                   sigma, etaA, etaB, dx, x0, zc;

	/* Construct Parent */
	_AnalyticSolution_AssignFromXML( self, cf, data );

	/* Create Analytic Fields */
	velocityField = Stg_ComponentFactory_ConstructByName( cf, "VelocityField", FeVariable, True, data );
	AnalyticSolution_RegisterFeVariableWithAnalyticFunction( self, velocityField, Velic_solJA_VelocityFunction );

	pressureField = Stg_ComponentFactory_ConstructByName( cf, "PressureField", FeVariable, True, data );
	AnalyticSolution_RegisterFeVariableWithAnalyticFunction( self, pressureField, Velic_solJA_PressureFunction );

	stressField = Stg_ComponentFactory_ConstructByName( cf, "StressField", FeVariable, False, data );
	if ( stressField )
		AnalyticSolution_RegisterFeVariableWithAnalyticFunction( self, stressField, Velic_solJA_StressFunction );

	strainRateField = Stg_ComponentFactory_ConstructByName( cf, "StrainRateField", FeVariable, False, data );
	if ( strainRateField  ) {
		AnalyticSolution_RegisterFeVariableWithAnalyticFunction( self, strainRateField, Velic_solJA_StrainRateFunction );
	}

	recoveredStrainRateField = Stg_ComponentFactory_ConstructByName( cf, "recoveredStrainRateField", FeVariable, False, data );
	if ( recoveredStrainRateField )
		AnalyticSolution_RegisterFeVariableWithAnalyticFunction( self, recoveredStrainRateField, Velic_solJA_StrainRateFunction );

	recoveredStressField = Stg_ComponentFactory_ConstructByName( cf, "recoveredStressField", FeVariable, False, data );
	if ( recoveredStressField )
		AnalyticSolution_RegisterFeVariableWithAnalyticFunction( self, recoveredStressField, Velic_solJA_StressFunction );
	
	sigma = Stg_ComponentFactory_GetRootDictDouble( cf, "solJA_sigma", 1.0 );
	etaA = Stg_ComponentFactory_GetRootDictDouble( cf, "solJA_etaA", 100.0 );
	etaB = Stg_ComponentFactory_GetRootDictDouble( cf, "solJA_etaB", 1.0 );
	dx = Stg_ComponentFactory_GetRootDictDouble( cf, "solJA_dx", 0.4 );
	x0 = Stg_ComponentFactory_GetRootDictDouble( cf, "solJA_x0", 0.2 );
	zc = Stg_ComponentFactory_GetRootDictDouble( cf, "solJA_zc", 0.8 );

	_Velic_solJA_Init( self, sigma, etaA, etaB, dx, x0, zc );

}

void* _Velic_solJA_DefaultNew( Name name ) {
	return _AnalyticSolution_New(
			sizeof(Velic_solJA),
			Velic_solJA_Type,
			_AnalyticSolution_Delete,
			_AnalyticSolution_Print,
			_AnalyticSolution_Copy,
			_Velic_solJA_DefaultNew,
			_Velic_solJA_AssignFromXML,
			_AnalyticSolution_Build,
			_AnalyticSolution_Initialise,
			_AnalyticSolution_Execute,
			_AnalyticSolution_Destroy,
			name );
}

Index Underworld_Velic_solJA_Register( PluginsManager* pluginsManager ) {
	return PluginsManager_Submit( pluginsManager, Velic_solJA_Type, "0", _Velic_solJA_DefaultNew );
}
