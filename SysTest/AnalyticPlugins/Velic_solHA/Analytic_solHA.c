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
** $Id: solHA.c 636 2006-09-08 03:07:06Z JulianGiordani $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <stdio.h>
#include <string.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include "Analytic_solHA.h"

const Type Velic_solHA_Type = "Underworld_Velic_solHA";

void Velic_solHA_PressureFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* pressure ) {
	Velic_solHA* self = (Velic_solHA*) analyticSolution;
	
	_Velic_solHA( coord, self->sigma, self->eta, self->dx, self->dy, self->x0, self->y0, NULL, pressure, NULL, NULL );
}

void Velic_solHA_VelocityFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* velocity ) {
	Velic_solHA* self = (Velic_solHA*) analyticSolution;
	
	_Velic_solHA( coord, self->sigma, self->eta, self->dx, self->dy, self->x0, self->y0, velocity, NULL, NULL, NULL );
}

void Velic_solHA_StressFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* stress ) {
	Velic_solHA* self = (Velic_solHA*) analyticSolution;
	
	_Velic_solHA( coord, self->sigma, self->eta, self->dx, self->dy, self->x0, self->y0, NULL, NULL, stress, NULL );
}


void Velic_solHA_StrainRateFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* strainRate ) {
	Velic_solHA* self = (Velic_solHA*) analyticSolution;
	
	_Velic_solHA( coord, self->sigma, self->eta, self->dx, self->dy, self->x0, self->y0, NULL, NULL, NULL, strainRate );
}

void _Velic_solHA_Init( Velic_solHA* self, double sigma, double eta, double dx, double dy, double x0, double y0 ) {
	// TODO	Bool                     correctInput = True;
	// Set the general AnalyticSolution (parent class) field function pointers to the specific analytic sols to use
        self->_getAnalyticVelocity    = Velic_solHA_VelocityFunction;
	self->_getAnalyticPressure    = Velic_solHA_PressureFunction;
	self->_getAnalyticTotalStress = Velic_solHA_StressFunction;
	self->_getAnalyticStrainRate  = Velic_solHA_StrainRateFunction;

	self->sigma = sigma;
	self->eta = eta;
	self->dx = dx;
	self->dy = dy;
	self->x0 = x0;
	self->y0 = y0;
}

void _Velic_solHA_Build( void* analyticSolution, void* data ) {
	Velic_solHA*          self  = (Velic_solHA*)analyticSolution;
	
	_AnalyticSolution_Build( self, data );
}

void _Velic_solHA_Construct( void* analyticSolution, Stg_ComponentFactory* cf, void* data ) {
	Velic_solHA* self = (Velic_solHA*) analyticSolution;
	FeVariable*              velocityField;
	FeVariable*              pressureField;
	FeVariable*              recoveredPressureField;
	FeVariable*              stressField;
	FeVariable*              strainRateField;
	FeVariable*              recoverdStrainRateField;
	FeVariable*              recoveredStressField;
	double                   sigma, eta, dx, dy, x0, y0;

	/* Construct Parent */
	_AnalyticSolution_Construct( self, cf, data );

	/* Create Analytic Fields */
	velocityField = Stg_ComponentFactory_ConstructByName( cf, "VelocityField", FeVariable, True, data );
	AnalyticSolution_RegisterFeVariableWithAnalyticFunction( self, velocityField, Velic_solHA_VelocityFunction );

	pressureField = Stg_ComponentFactory_ConstructByName( cf, "PressureField", FeVariable, True, data );
	AnalyticSolution_RegisterFeVariableWithAnalyticFunction( self, pressureField, Velic_solHA_PressureFunction );

	stressField = Stg_ComponentFactory_ConstructByName( cf, "StressField", FeVariable, False, data );
	if ( stressField ) 
		AnalyticSolution_RegisterFeVariableWithAnalyticFunction( self, stressField, Velic_solHA_StressFunction );

	strainRateField = Stg_ComponentFactory_ConstructByName( cf, "StrainRateField", FeVariable, False, data );
	if ( strainRateField  )
		AnalyticSolution_RegisterFeVariableWithAnalyticFunction( self, strainRateField, Velic_solHA_StrainRateFunction );

	recoverdStrainRateField = Stg_ComponentFactory_ConstructByName( cf, "recoveredStrainRateField", FeVariable, False, data );
	if ( recoverdStrainRateField )
		AnalyticSolution_RegisterFeVariableWithAnalyticFunction( self, recoverdStrainRateField, Velic_solHA_StrainRateFunction );

	recoveredStressField = Stg_ComponentFactory_ConstructByName( cf, "recoveredStressField", FeVariable, False, data );
	if ( recoveredStressField )
		AnalyticSolution_RegisterFeVariableWithAnalyticFunction( self, recoveredStressField, Velic_solHA_StressFunction );

	recoveredPressureField = Stg_ComponentFactory_ConstructByName( cf, "recoveredPressureField", FeVariable, False, data );
	if ( recoveredPressureField )
		AnalyticSolution_RegisterFeVariableWithAnalyticFunction( self, recoveredPressureField, Velic_solHA_PressureFunction );

	sigma = Stg_ComponentFactory_GetRootDictDouble( cf, "solHA_sigma", 1.0 );
	eta = Stg_ComponentFactory_GetRootDictDouble( cf, "solHA_eta", 1.0 );
	dx = Stg_ComponentFactory_GetRootDictDouble( cf, "solHA_dx", 0.4 );
	dy = Stg_ComponentFactory_GetRootDictDouble( cf, "solHA_dy", 0.35 );
	x0 = Stg_ComponentFactory_GetRootDictDouble( cf, "solHA_x0", 0.3 );
	y0 = Stg_ComponentFactory_GetRootDictDouble( cf, "solHA_y0", 0.6 );

	_Velic_solHA_Init( self, sigma, eta, dx, dy, x0, y0 );

}

void* _Velic_solHA_DefaultNew( Name name ) {
	return _AnalyticSolution_New(
			sizeof(Velic_solHA),
			Velic_solHA_Type,
			_AnalyticSolution_Delete,
			_AnalyticSolution_Print,
			_AnalyticSolution_Copy,
			_Velic_solHA_DefaultNew,
			_Velic_solHA_Construct,
			_Velic_solHA_Build,
			_AnalyticSolution_Initialise,
			_AnalyticSolution_Execute,
			_AnalyticSolution_Destroy,
			name );
}

Index Underworld_Velic_solHA_Register( PluginsManager* pluginsManager ) {
	return PluginsManager_Submit( pluginsManager, Velic_solHA_Type, "0", _Velic_solHA_DefaultNew );
}
