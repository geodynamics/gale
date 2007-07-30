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
** $Id: solKz.c 636 2006-09-08 03:07:06Z JulianGiordani $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include "Analytic_solKz.h"

const Type Velic_solKz_Type = "ExperimentalUnderworld_Velic_solKz";

void Velic_solKz_PressureFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* pressure ) {
	Velic_solKz* self = (Velic_solKz*) analyticSolution;
	
	_Velic_solKz( coord, self->sigma, self->km, self->n, self->B, NULL, pressure, NULL, NULL );
}

void Velic_solKz_VelocityFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* velocity ) {
	Velic_solKz* self = (Velic_solKz*) analyticSolution;
	
	_Velic_solKz( coord, self->sigma, self->km, self->n, self->B, velocity, NULL, NULL, NULL );
}

void Velic_solKz_StressFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* stress ) {
	Velic_solKz* self = (Velic_solKz*) analyticSolution;
	
	_Velic_solKz( coord, self->sigma, self->km, self->n, self->B, NULL, NULL, stress, NULL );
}


void Velic_solKz_StrainRateFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* strainRate ) {
	Velic_solKz* self = (Velic_solKz*) analyticSolution;
	
	_Velic_solKz( coord, self->sigma, self->km, self->n, self->B, NULL, NULL, NULL, strainRate );
}

	
void _Velic_solKz_Init( Velic_solKz* self, double sigma, double km, double B, int n ) {
	// TODO	Bool                     correctInput = True;
	// Set the general AnalyticSolution (parent class) field function pointers to the specific analytic sols to use
        self->_getAnalyticVelocity    = Velic_solKz_VelocityFunction;
	self->_getAnalyticPressure    = Velic_solKz_PressureFunction;
	self->_getAnalyticTotalStress = Velic_solKz_StressFunction;
	self->_getAnalyticStrainRate  = Velic_solKz_StrainRateFunction;

	self->sigma = sigma;
	self->km = km;
	self->B = B;
	self->n = n;
}

void _Velic_solKz_Build( void* analyticSolution, void* data ) {
	Velic_solKz*          self  = (Velic_solKz*)analyticSolution;

	AnalyticSolution_BuildAllAnalyticFields( self, data );
	
	_AnalyticSolution_Build( self, data );
}

void _Velic_solKz_Construct( void* analyticSolution, Stg_ComponentFactory* cf, void* data ) {
	Velic_solKz* self = (Velic_solKz*) analyticSolution;
	FeVariable*              velocityField;
	FeVariable*              pressureField, *recoveredPressureField;
	FeVariable*              stressField;
	FeVariable*              strainRateField;
	FeVariable*              recoverdStrainRateField;
	FeVariable*              recoveredStressField;
	double                   sigma, km, B;
	int                      n;

	/* Construct Parent */
	_AnalyticSolution_Construct( self, cf, data );

	/* Create Analytic Fields */
	velocityField = Stg_ComponentFactory_ConstructByName( cf, "VelocityField", FeVariable, True, data );
	AnalyticSolution_RegisterFeVariableWithAnalyticFunction( self, velocityField, Velic_solKz_VelocityFunction );

	pressureField = Stg_ComponentFactory_ConstructByName( cf, "PressureField", FeVariable, True, data );
	AnalyticSolution_RegisterFeVariableWithAnalyticFunction( self, pressureField, Velic_solKz_PressureFunction );

	stressField = Stg_ComponentFactory_ConstructByName( cf, "StressField", FeVariable, False, data );
	if ( stressField ) 
		AnalyticSolution_RegisterFeVariableWithAnalyticFunction( self, stressField, Velic_solKz_StressFunction );

	strainRateField = Stg_ComponentFactory_ConstructByName( cf, "StrainRateField", FeVariable, False, data );
	if ( strainRateField  )
		AnalyticSolution_RegisterFeVariableWithAnalyticFunction( self, strainRateField, Velic_solKz_StrainRateFunction );

	recoverdStrainRateField = Stg_ComponentFactory_ConstructByName( cf, "recoveredStrainRateField", FeVariable, False, data );
	if ( recoverdStrainRateField )
		AnalyticSolution_RegisterFeVariableWithAnalyticFunction( self, recoverdStrainRateField, Velic_solKz_StrainRateFunction );

	recoveredStressField = Stg_ComponentFactory_ConstructByName( cf, "recoveredStressField", FeVariable, False, data );
	if ( recoveredStressField )
		AnalyticSolution_RegisterFeVariableWithAnalyticFunction( self, recoveredStressField, Velic_solKz_StressFunction );

	recoveredPressureField = Stg_ComponentFactory_ConstructByName( cf, "recoveredPressureField", FeVariable, False, data );
	if ( recoveredPressureField )
		AnalyticSolution_RegisterFeVariableWithAnalyticFunction( self, recoveredPressureField, Velic_solKz_PressureFunction );

	sigma = Stg_ComponentFactory_GetRootDictDouble( cf, "solKz_sigma", 1.0 );
	km = Stg_ComponentFactory_GetRootDictDouble( cf, "solKz_km", M_PI );
	B = Stg_ComponentFactory_GetRootDictDouble( cf, "solKz_B", 1.0 );
	n = Stg_ComponentFactory_GetRootDictInt( cf, "solKz_n", 1 );
	
	_Velic_solKz_Init( self, sigma, km, B, n );
}

void* _Velic_solKz_DefaultNew( Name name ) {
	return _AnalyticSolution_New(
			sizeof(Velic_solKz),
			Velic_solKz_Type,
			_AnalyticSolution_Delete,
			_AnalyticSolution_Print,
			_AnalyticSolution_Copy,
			_Velic_solKz_DefaultNew,
			_Velic_solKz_Construct,
			_Velic_solKz_Build,
			_AnalyticSolution_Initialise,
			_AnalyticSolution_Execute,
			_AnalyticSolution_Destroy,
			name );
}

Index ExperimentalUnderworld_Velic_solKz_Register( PluginsManager* pluginsManager ) {
	return PluginsManager_Submit( pluginsManager, Velic_solKz_Type, "0", _Velic_solKz_DefaultNew );
}
