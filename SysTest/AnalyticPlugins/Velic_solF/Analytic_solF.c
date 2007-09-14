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
** $Id: solF.c 636 2006-09-08 03:07:06Z JulianGiordani $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include "Analytic_solF.h"

const Type Velic_solF_Type = "Underworld_Velic_solF";

void Velic_solF_PressureFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* pressure ) {
	Velic_solF* self = (Velic_solF*) analyticSolution;
	
	_Velic_solF( coord, self->sigma, self->etaA, self->etaB, self->xc, self->zc,
		       	NULL, pressure, NULL, NULL );
}

void Velic_solF_VelocityFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* velocity ) {
	Velic_solF* self = (Velic_solF*) analyticSolution;
	
	_Velic_solF( coord, self->sigma, self->etaA, self->etaB, self->xc, self->zc,
		       	velocity, NULL, NULL, NULL );
}

void Velic_solF_StressFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* stress ) {
	Velic_solF* self = (Velic_solF*) analyticSolution;
	
	_Velic_solF( coord, self->sigma, self->etaA, self->etaB, self->xc, self->zc,
		       	NULL, NULL, stress, NULL );
}


void Velic_solF_StrainRateFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* strainRate ) {
	Velic_solF* self = (Velic_solF*) analyticSolution;
	
	_Velic_solF( coord, self->sigma, self->etaA, self->etaB, self->xc, self->zc,
			NULL, NULL, NULL, strainRate );
}

void _Velic_solF_Init( Velic_solF* self, double sigma, double etaA, double etaB, double xc, double zc ) {
	// TODO	Bool                     correctInput = True;
	// Set the general AnalyticSolution (parent class) field function pointers to the specific analytic sols to use
        self->_getAnalyticVelocity    = Velic_solF_VelocityFunction;
	self->_getAnalyticPressure    = Velic_solF_PressureFunction;
	self->_getAnalyticTotalStress = Velic_solF_StressFunction;
	self->_getAnalyticStrainRate  = Velic_solF_StrainRateFunction;

	self->sigma = sigma;
	self->etaA = etaA;
	self->etaB = etaB;
	self->xc = xc;
	self->zc = zc;
}

void _Velic_solF_Construct( void* analyticSolution, Stg_ComponentFactory* cf, void* data ) {
	Velic_solF* self = (Velic_solF*) analyticSolution;
	FeVariable*              velocityField;
	FeVariable*              pressureField;
	FeVariable*              stressField;
	FeVariable*              strainRateField;
	FeVariable*              recoverdStrainRateField;
	FeVariable*              recoveredStressField;
 	double                   sigma, etaA, etaB, xc, zc;

	/* Construct Parent */
	_AnalyticSolution_Construct( self, cf, data );

	/* Create Analytic Fields */
	velocityField = Stg_ComponentFactory_ConstructByName( cf, "VelocityField", FeVariable, True, data );
	AnalyticSolution_RegisterFeVariableWithAnalyticFunction( self, velocityField, Velic_solF_VelocityFunction );

	pressureField = Stg_ComponentFactory_ConstructByName( cf, "PressureField", FeVariable, True, data );
	AnalyticSolution_RegisterFeVariableWithAnalyticFunction( self, pressureField, Velic_solF_PressureFunction );

	stressField = Stg_ComponentFactory_ConstructByName( cf, "StressField", FeVariable, False, data );
	if ( stressField )
		AnalyticSolution_RegisterFeVariableWithAnalyticFunction( self, stressField, Velic_solF_StressFunction );

	strainRateField = Stg_ComponentFactory_ConstructByName( cf, "StrainRateField", FeVariable, False, data );
	if ( strainRateField  ) {
		AnalyticSolution_RegisterFeVariableWithAnalyticFunction( self, strainRateField, Velic_solF_StrainRateFunction );
	}

	recoverdStrainRateField = Stg_ComponentFactory_ConstructByName( cf, "recoveredStrainRate", FeVariable, False, data );
	if ( recoverdStrainRateField )
		AnalyticSolution_RegisterFeVariableWithAnalyticFunction( self, recoverdStrainRateField, Velic_solF_StrainRateFunction );

	recoveredStressField = Stg_ComponentFactory_ConstructByName( cf, "recoveredStressField", FeVariable, False, data );
	if ( recoveredStressField )
		AnalyticSolution_RegisterFeVariableWithAnalyticFunction( self, recoveredStressField, Velic_solF_StressFunction );

	sigma = Stg_ComponentFactory_GetRootDictDouble( cf, "solF_sigma", 1.0 );
	etaA = Stg_ComponentFactory_GetRootDictDouble( cf, "solF_etaA", 2.0 );
	etaB = Stg_ComponentFactory_GetRootDictDouble( cf, "solF_etaB", 1.0 );
	xc = Stg_ComponentFactory_GetRootDictDouble( cf, "solF_xc", 0.25 );
	zc = Stg_ComponentFactory_GetRootDictDouble( cf, "solF_zc", 0.7 );

	_Velic_solF_Init( self, sigma, etaA, etaB, xc, zc );
	

}

void* _Velic_solF_DefaultNew( Name name ) {
	return _AnalyticSolution_New(
			sizeof(Velic_solF),
			Velic_solF_Type,
			_AnalyticSolution_Delete,
			_AnalyticSolution_Print,
			_AnalyticSolution_Copy,
			_Velic_solF_DefaultNew,
			_Velic_solF_Construct,
			_AnalyticSolution_Build,
			_AnalyticSolution_Initialise,
			_AnalyticSolution_Execute,
			_AnalyticSolution_Destroy,
			name );
}

Index Underworld_Velic_solF_Register( PluginsManager* pluginsManager ) {
	return PluginsManager_Submit( pluginsManager, Velic_solF_Type, "0", _Velic_solF_DefaultNew );
}
