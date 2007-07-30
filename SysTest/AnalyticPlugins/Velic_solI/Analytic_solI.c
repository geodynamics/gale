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
** $Id: solI.c 636 2006-09-08 03:07:06Z JulianGiordani $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include "Analytic_solI.h"

const Type Velic_solI_Type = "ExperimentalUnderworld_Velic_solI";

void Velic_solI_PressureFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* pressure ) {
	Velic_solI* self = (Velic_solI*) analyticSolution;
	
	_Velic_solI( coord, self->sigma, self->B, self->xc, NULL, pressure, NULL, NULL );
}

void Velic_solI_VelocityFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* velocity ) {
	Velic_solI* self = (Velic_solI*) analyticSolution;
	
	_Velic_solI( coord, self->sigma, self->B, self->xc, velocity, NULL, NULL, NULL );
}

void Velic_solI_StressFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* stress ) {
	Velic_solI* self = (Velic_solI*) analyticSolution;
	
	_Velic_solI( coord, self->sigma, self->B, self->xc, NULL, NULL, stress, NULL );
}


void Velic_solI_StrainRateFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* strainRate ) {
	Velic_solI* self = (Velic_solI*) analyticSolution;
	
	_Velic_solI( coord, self->sigma, self->B, self->xc, NULL, NULL, NULL, strainRate );
}

void _Velic_solI_Init( Velic_solI* self, double sigma, double B, double xc ) {
	//	Bool                     correctInput = True;
	// Set the general AnalyticSolution (parent class) field function pointers to the specific analytic sols to use
        self->_getAnalyticVelocity    = Velic_solI_VelocityFunction;
	self->_getAnalyticPressure    = Velic_solI_PressureFunction;
	self->_getAnalyticTotalStress = Velic_solI_StressFunction;
	self->_getAnalyticStrainRate  = Velic_solI_StrainRateFunction;

	self->sigma = sigma;
	self->B = B;
	self->xc = xc;
}

void _Velic_solI_Construct( void* analyticSolution, Stg_ComponentFactory* cf, void* data ) {
	Velic_solI* self = (Velic_solI*) analyticSolution;
	FeVariable*              velocityField;
	FeVariable*              pressureField;
	FeVariable*              stressField;
	FeVariable*              strainRateField;
	FeVariable*              recoverdStrainRateField;
	FeVariable*              recoveredStressField;
	double                   sigma, B, xc;

	/* Construct Parent */
	_AnalyticSolution_Construct( self, cf, data );

	/* Create Analytic Fields */
	velocityField = Stg_ComponentFactory_ConstructByName( cf, "VelocityField", FeVariable, True, data );
	AnalyticSolution_CreateAnalyticVectorField( self, velocityField, Velic_solI_VelocityFunction );

	pressureField = Stg_ComponentFactory_ConstructByName( cf, "PressureField", FeVariable, True, data );
	AnalyticSolution_CreateAnalyticField( self, pressureField, Velic_solI_PressureFunction );

	stressField = Stg_ComponentFactory_ConstructByName( cf, "StressField", FeVariable, False, data );
	if ( stressField )
		AnalyticSolution_CreateAnalyticSymmetricTensorField( self, stressField, Velic_solI_StressFunction );

	strainRateField = Stg_ComponentFactory_ConstructByName( cf, "StrainRateField", FeVariable, False, data );
	if ( strainRateField  ) {
		AnalyticSolution_CreateAnalyticSymmetricTensorField( self, strainRateField, Velic_solI_StrainRateFunction );
	}

	recoverdStrainRateField = Stg_ComponentFactory_ConstructByName( cf, "recoveredStrainRate", FeVariable, False, data );
	if ( recoverdStrainRateField )
		AnalyticSolution_CreateAnalyticSymmetricTensorField( self, recoverdStrainRateField, Velic_solI_StrainRateFunction );

	recoveredStressField = Stg_ComponentFactory_ConstructByName( cf, "recoveredStressField", FeVariable, False, data );
	if ( recoveredStressField )
		AnalyticSolution_CreateAnalyticSymmetricTensorField( self, recoveredStressField, Velic_solI_StressFunction );
	
	sigma = Stg_ComponentFactory_GetRootDictDouble( cf, "solI_sigma", 1.0 );
	B = Stg_ComponentFactory_GetRootDictDouble( cf, "solI_B", 1.5 );
	xc = Stg_ComponentFactory_GetRootDictDouble( cf, "solI_xc", 0.3 );

	_Velic_solI_Init( self, sigma, B, xc );

}

void* _Velic_solI_DefaultNew( Name name ) {
	return _AnalyticSolution_New(
			sizeof(Velic_solI),
			Velic_solI_Type,
			_AnalyticSolution_Delete,
			_AnalyticSolution_Print,
			_AnalyticSolution_Copy,
			_Velic_solI_DefaultNew,
			_Velic_solI_Construct,
			_AnalyticSolution_Build,
			_AnalyticSolution_Initialise,
			_AnalyticSolution_Execute,
			_AnalyticSolution_Destroy,
			name );
}

Index ExperimentalUnderworld_Velic_solI_Register( PluginsManager* pluginsManager ) {
	return PluginsManager_Submit( pluginsManager, Velic_solI_Type, "0", _Velic_solI_DefaultNew );
}
