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
** $Id: solHy.c 636 2006-09-08 03:07:06Z JulianGiordani $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include "Analytic_solHy.h"

const Type Velic_solHy_Type = "Underworld_Velic_solHy";

void Velic_solHy_PressureFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* pressure ) {
	Velic_solHy* self = (Velic_solHy*) analyticSolution;
	
	_Velic_solHy( coord, self->etaA, self->etaB, self->xc, self->n, NULL, pressure, NULL, NULL );
}

void Velic_solHy_VelocityFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* velocity ) {
	Velic_solHy* self = (Velic_solHy*) analyticSolution;
	
	_Velic_solHy( coord, self->etaA, self->etaB, self->xc, self->n, velocity, NULL, NULL, NULL );
}

void Velic_solHy_StressFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* stress ) {
	Velic_solHy* self = (Velic_solHy*) analyticSolution;
	
	_Velic_solHy( coord, self->etaA, self->etaB, self->xc, self->n, NULL, NULL, stress, NULL );
}


void Velic_solHy_StrainRateFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* strainRate ) {
	Velic_solHy* self = (Velic_solHy*) analyticSolution;
	
	_Velic_solHy( coord, self->etaA, self->etaB, self->xc, self->n, NULL, NULL, NULL, strainRate );
}

void _Velic_solHy_Init( Velic_solHy* self, double etaA, double etaB, double xc, int n ) {
	Bool                     correctInput = True;
	// Set the general AnalyticSolution (parent class) field function pointers to the specific analytic sols to use
        self->_getAnalyticVelocity    = Velic_solHy_VelocityFunction;
	self->_getAnalyticPressure    = Velic_solHy_PressureFunction;
	self->_getAnalyticTotalStress = Velic_solHy_StressFunction;
	self->_getAnalyticStrainRate  = Velic_solHy_StrainRateFunction;

	self->etaA = etaA;
	self->etaB = etaB;
	self->xc = xc;
	self->n = n;

	if( (self->xc < 0.0) || self->xc > 1.0 ) 
	      correctInput = False;	

	Journal_Firewall( correctInput,  Journal_Register( Error_Type, "solHy_ErrorStream") ,
			"Error: In func %s. The input parameters you supplied to the analytic Solution are incorrect.\nValid range of xc values is [0,1]. Currently it is %f\n", __func__, self->xc );
}

void _Velic_solHy_AssignFromXML( void* analyticSolution, Stg_ComponentFactory* cf, void* data ) {
	Velic_solHy* self = (Velic_solHy*) analyticSolution;
	FeVariable*              velocityField;
	FeVariable*              pressureField;
	FeVariable*              stressField;
	FeVariable*              strainRateField;
	FeVariable*              recoveredStrainRateField;
	FeVariable*              recoveredStressField;
	double                   etaA, etaB, xc;
	int                      n;

	/* Construct Parent */
	_AnalyticSolution_AssignFromXML( self, cf, data );

	/* Create Analytic Fields */
	velocityField = Stg_ComponentFactory_ConstructByName( cf, "VelocityField", FeVariable, True, data );
	AnalyticSolution_RegisterFeVariableWithAnalyticFunction( self, velocityField, Velic_solHy_VelocityFunction );

	pressureField = Stg_ComponentFactory_ConstructByName( cf, "PressureField", FeVariable, True, data );
	AnalyticSolution_RegisterFeVariableWithAnalyticFunction( self, pressureField, Velic_solHy_PressureFunction );

	stressField = Stg_ComponentFactory_ConstructByName( cf, "StressField", FeVariable, False, data );
	if ( stressField )
		AnalyticSolution_RegisterFeVariableWithAnalyticFunction( self, stressField, Velic_solHy_StressFunction );

	strainRateField = Stg_ComponentFactory_ConstructByName( cf, "StrainRateField", FeVariable, False, data );
	if ( strainRateField  ) {
		AnalyticSolution_RegisterFeVariableWithAnalyticFunction( self, strainRateField, Velic_solHy_StrainRateFunction );
	}

	recoveredStrainRateField = Stg_ComponentFactory_ConstructByName( cf, "recoveredStrainRateField", FeVariable, False, data );
	if ( recoveredStrainRateField )
		AnalyticSolution_RegisterFeVariableWithAnalyticFunction( self, recoveredStrainRateField, Velic_solHy_StrainRateFunction );

	recoveredStressField = Stg_ComponentFactory_ConstructByName( cf, "recoveredStressField", FeVariable, False, data );
	if ( recoveredStressField )
		AnalyticSolution_RegisterFeVariableWithAnalyticFunction( self, recoveredStressField, Velic_solHy_StressFunction );
	
	etaA = Stg_ComponentFactory_GetRootDictDouble( cf, "solHy_etaA", 1.0 );
	etaB = Stg_ComponentFactory_GetRootDictDouble( cf, "solHy_etaB", 2.0 );
	xc = Stg_ComponentFactory_GetRootDictDouble( cf, "solHy_xc", 0.25 );
	n = Stg_ComponentFactory_GetRootDictInt( cf, "solHy_n", 1 );

	_Velic_solHy_Init( self, etaA, etaB, xc, n );
	
}

void* _Velic_solHy_DefaultNew( Name name ) {
	return _AnalyticSolution_New(
			sizeof(Velic_solHy),
			Velic_solHy_Type,
			_AnalyticSolution_Delete,
			_AnalyticSolution_Print,
			_AnalyticSolution_Copy,
			_Velic_solHy_DefaultNew,
			_Velic_solHy_AssignFromXML,
			_AnalyticSolution_Build,
			_AnalyticSolution_Initialise,
			_AnalyticSolution_Execute,
			_AnalyticSolution_Destroy,
			name );
}

Index Underworld_Velic_solHy_Register( PluginsManager* pluginsManager ) {
	return PluginsManager_Submit( pluginsManager, Velic_solHy_Type, "0", _Velic_solHy_DefaultNew );
}
