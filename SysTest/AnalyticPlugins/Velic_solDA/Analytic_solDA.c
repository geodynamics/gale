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
** $Id: solDA.c 636 2006-09-08 03:07:06Z JulianGiordani $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include "Analytic_solDA.h"

const Type Velic_solDA_Type = "Underworld_Velic_solDA";

void Velic_solDA_PressureFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* pressure ) {
	Velic_solDA* self = (Velic_solDA*) analyticSolution;
	
	_Velic_solDA( coord, self->sigma, self->etaA, self->etaB, self->zc, self->dx, self->x0,
		       	NULL, pressure, NULL, NULL );
}

void Velic_solDA_VelocityFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* velocity ) {
	Velic_solDA* self = (Velic_solDA*) analyticSolution;
	
	_Velic_solDA( coord, self->sigma, self->etaA, self->etaB, self->zc, self->dx, self->x0,
		       	velocity, NULL, NULL, NULL );
}

void Velic_solDA_StressFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* stress ) {
	Velic_solDA* self = (Velic_solDA*) analyticSolution;
	
	_Velic_solDA( coord, self->sigma, self->etaA, self->etaB, self->zc, self->dx, self->x0,
		       	NULL, NULL, stress, NULL );
}


void Velic_solDA_StrainRateFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* strainRate ) {
	Velic_solDA* self = (Velic_solDA*) analyticSolution;
	
	_Velic_solDA( coord, self->sigma, self->etaA, self->etaB, self->zc, self->dx, self->x0,
		       	NULL, NULL, NULL, strainRate );
}

void _Velic_solDA_Init( Velic_solDA* self, double sigma, double etaA, double etaB, double zc, double dx, double x0 ) {
	// Set the general AnalyticSolution (parent class) field function pointers to the specific analytic sols to use
        self->_getAnalyticVelocity    = Velic_solDA_VelocityFunction;
	self->_getAnalyticPressure    = Velic_solDA_PressureFunction;
	self->_getAnalyticTotalStress = Velic_solDA_StressFunction;
	self->_getAnalyticStrainRate  = Velic_solDA_StrainRateFunction;
	
	self->sigma = sigma;
	self->etaA = etaA;
	self->etaB = etaB;
	self->zc = zc;
	self->dx = dx;
	self->x0 = x0;
}

void _Velic_solDA_Construct( void* analyticSolution, Stg_ComponentFactory* cf, void* data ) {
	Velic_solDA* self = (Velic_solDA*) analyticSolution;
	FeVariable*              velocityField;
	FeVariable*              pressureField;
	FeVariable*              stressField;
	FeVariable*              strainRateField;
	FeVariable*              recoverdStrainRateField;
	FeVariable*              recoveredStressField;
	double                   sigma, etaA, etaB, zc, dx, x0;
	double                   startX, endX;
	
	/* Construct Parent */
	_AnalyticSolution_Construct( self, cf, data );

	/* Create Analytic Fields */
	velocityField = Stg_ComponentFactory_ConstructByName( cf, "VelocityField", FeVariable, True, data );
	AnalyticSolution_RegisterFeVariableWithAnalyticFunction( self, velocityField, Velic_solDA_VelocityFunction );

	pressureField = Stg_ComponentFactory_ConstructByName( cf, "PressureField", FeVariable, True, data );
	AnalyticSolution_RegisterFeVariableWithAnalyticFunction( self, pressureField, Velic_solDA_PressureFunction );

	stressField = Stg_ComponentFactory_ConstructByName( cf, "StressField", FeVariable, False, data );
	if ( stressField )
		AnalyticSolution_RegisterFeVariableWithAnalyticFunction( self, stressField, Velic_solDA_StressFunction );

	strainRateField = Stg_ComponentFactory_ConstructByName( cf, "StrainRateField", FeVariable, False, data );
	if ( strainRateField  ) {
		AnalyticSolution_RegisterFeVariableWithAnalyticFunction( self, strainRateField, Velic_solDA_StrainRateFunction );
	}

	recoverdStrainRateField = Stg_ComponentFactory_ConstructByName( cf, "recoveredStrainRate", FeVariable, False, data );
	if ( recoverdStrainRateField )
		AnalyticSolution_RegisterFeVariableWithAnalyticFunction( self, recoverdStrainRateField, Velic_solDA_StrainRateFunction );

	recoveredStressField = Stg_ComponentFactory_ConstructByName( cf, "recoveredStressField", FeVariable, False, data );
	if ( recoveredStressField )
		AnalyticSolution_RegisterFeVariableWithAnalyticFunction( self, recoveredStressField, Velic_solDA_StressFunction );

	sigma = Stg_ComponentFactory_GetRootDictDouble( cf, "solDA_sigma", 3.0 );
	etaA = Stg_ComponentFactory_GetRootDictDouble( cf, "solDA_etaA", 1.0 );
	etaB = Stg_ComponentFactory_GetRootDictDouble( cf, "solDA_etaB", 2.0 );
	zc = Stg_ComponentFactory_GetRootDictDouble( cf, "solDA_zc", 0.8 );
	startX = Stg_ComponentFactory_GetRootDictDouble( cf, "solCA_startX", 0.4 );
	endX   = Stg_ComponentFactory_GetRootDictDouble( cf, "solCA_endX", 0.8 );

	dx = endX - startX;
	x0 = 0.5 * (endX + startX);
	
	_Velic_solDA_Init( self, sigma, etaA, etaB, zc, dx, x0 );
	

}

void* _Velic_solDA_DefaultNew( Name name ) {
	return _AnalyticSolution_New(
			sizeof(Velic_solDA),
			Velic_solDA_Type,
			_AnalyticSolution_Delete,
			_AnalyticSolution_Print,
			_AnalyticSolution_Copy,
			_Velic_solDA_DefaultNew,
			_Velic_solDA_Construct,
			_AnalyticSolution_Build,
			_AnalyticSolution_Initialise,
			_AnalyticSolution_Execute,
			_AnalyticSolution_Destroy,
			name );
}

Index Underworld_Velic_solDA_Register( PluginsManager* pluginsManager ) {
	return PluginsManager_Submit( pluginsManager, Velic_solDA_Type, "0", _Velic_solDA_DefaultNew );
}
