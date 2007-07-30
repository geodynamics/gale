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
** $Id: solC.c 636 2006-09-08 03:07:06Z JulianGiordani $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include "Analytic_solC.h"

const Type Velic_solC_Type = "Underworld_Velic_solC";

void Velic_solC_PressureFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* pressure ) {
	Velic_solC* self = (Velic_solC*) analyticSolution;
	
	_Velic_solC( coord, self->sigma, self->eta, self->x_c, NULL, pressure, NULL, NULL );
}

void Velic_solC_VelocityFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* velocity ) {
	Velic_solC* self = (Velic_solC*) analyticSolution;
	
	_Velic_solC( coord, self->sigma, self->eta, self->x_c, velocity, NULL, NULL, NULL );
}

void Velic_solC_StressFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* stress ) {
	Velic_solC* self = (Velic_solC*) analyticSolution;
	
	_Velic_solC( coord, self->sigma, self->eta, self->x_c, NULL, NULL, stress, NULL );
}


void Velic_solC_StrainRateFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* strainRate ) {
	Velic_solC* self = (Velic_solC*) analyticSolution;
	
	_Velic_solC( coord, self->sigma, self->eta, self->x_c, NULL, NULL, NULL, strainRate );
}

void _Velic_solC_Init( Velic_solC* self, double sigma, double eta, double x_c ) {
	Bool                     isCorrectInput = True;

	// Set the general AnalyticSolution (parent class) field function pointers to the specific analytic sols to use
        self->_getAnalyticVelocity    = Velic_solC_VelocityFunction;
	self->_getAnalyticPressure    = Velic_solC_PressureFunction;
	self->_getAnalyticTotalStress = Velic_solC_StressFunction;
	self->_getAnalyticStrainRate  = Velic_solC_StrainRateFunction;

	self->sigma = sigma;
	self->eta = eta;
	self->x_c = x_c;

	isCorrectInput = _checkInputParams( self );
	Journal_Firewall( isCorrectInput , Journal_Register( Error_Type, "Velic_solA" ),
			"Error in function %s: Bad Input parameters, solution check valid values in .tex documentation\n",
			__func__ );
}

void _Velic_solC_Construct( void* analyticSolution, Stg_ComponentFactory* cf, void* data ) {
	Velic_solC* self = (Velic_solC*) analyticSolution;
	FeVariable*              velocityField;
	FeVariable*              pressureField;
	FeVariable*              stressField;
	FeVariable*              strainRateField;
	FeVariable*              recoverdStrainRateField;
	FeVariable*              recoveredStressField;
	double                   sigma, eta, x_c;

	/* Construct Parent */
	_AnalyticSolution_Construct( self, cf, data );

	/* Create Analytic Fields */
	velocityField = Stg_ComponentFactory_ConstructByName( cf, "VelocityField", FeVariable, True, data );
	AnalyticSolution_CreateAnalyticVectorField( self, velocityField, Velic_solC_VelocityFunction );

	pressureField = Stg_ComponentFactory_ConstructByName( cf, "PressureField", FeVariable, True, data );
	AnalyticSolution_CreateAnalyticField( self, pressureField, Velic_solC_PressureFunction );

	stressField = Stg_ComponentFactory_ConstructByName( cf, "StressField", FeVariable, False, data );
	if ( stressField )
		AnalyticSolution_CreateAnalyticSymmetricTensorField( self, stressField, Velic_solC_StressFunction );

	strainRateField = Stg_ComponentFactory_ConstructByName( cf, "StrainRateField", FeVariable, False, data );
	if ( strainRateField  ) {
		AnalyticSolution_CreateAnalyticSymmetricTensorField( self, strainRateField, Velic_solC_StrainRateFunction );
	}

	recoverdStrainRateField = Stg_ComponentFactory_ConstructByName( cf, "recoveredStrainRate", FeVariable, False, data );
	if ( recoverdStrainRateField )
		AnalyticSolution_CreateAnalyticSymmetricTensorField( self, recoverdStrainRateField, Velic_solC_StrainRateFunction );

	recoveredStressField = Stg_ComponentFactory_ConstructByName( cf, "recoveredStressField", FeVariable, False, data );
	if ( recoveredStressField )
		AnalyticSolution_CreateAnalyticSymmetricTensorField( self, recoveredStressField, Velic_solC_StressFunction );
	
	self->sigma = Stg_ComponentFactory_GetRootDictDouble( cf, "solC_sigma", 1.0 );
	self->eta = Stg_ComponentFactory_GetRootDictDouble( cf, "solC_eta", 1.0 );
	self->x_c = Stg_ComponentFactory_GetRootDictDouble( cf, "solC_xc", 0.4 );

	_Velic_solC_Init( self, sigma, eta, x_c );
}

Bool _checkInputParams( Velic_solC* self ) {
	return ( 
			( self->sigma > 0.0 ) && ( self->eta > 0.0 ) &&
			( self->x_c >= 0.0 ) && (self->x_c <= 1.0 ) 
		);
}

void* _Velic_solC_DefaultNew( Name name ) {
	return _AnalyticSolution_New(
			sizeof(Velic_solC),
			Velic_solC_Type,
			_AnalyticSolution_Delete,
			_AnalyticSolution_Print,
			_AnalyticSolution_Copy,
			_Velic_solC_DefaultNew,
			_Velic_solC_Construct,
			_AnalyticSolution_Build,
			_AnalyticSolution_Initialise,
			_AnalyticSolution_Execute,
			_AnalyticSolution_Destroy,
			name );
}

Index Underworld_Velic_solC_Register( PluginsManager* pluginsManager ) {
	return PluginsManager_Submit( pluginsManager, Velic_solC_Type, "0", _Velic_solC_DefaultNew );
}
