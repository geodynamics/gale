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
** $Id: solB.c 636 2006-09-08 03:07:06Z JulianGiordani $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include<mpi.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<StGermain/StGermain.h>
#include<StgFEM/StgFEM.h>
#include "Analytic_solB.h"

const Type Velic_solB_Type = "Underworld_Velic_solB";

void Velic_solB_PressureFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* pressure ) {
	Velic_solB* self = (Velic_solB*) analyticSolution;
	
	_Velic_solB( coord, self->sigma, self->Z, self->n, self->km, NULL, pressure, NULL, NULL );
}

void Velic_solB_VelocityFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* velocity ) {
	Velic_solB* self = (Velic_solB*) analyticSolution;
	
	_Velic_solB( coord, self->sigma, self->Z, self->n, self->km, velocity, NULL, NULL, NULL );
}

void Velic_solB_StressFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* stress ) {
	Velic_solB* self = (Velic_solB*) analyticSolution;
	
	_Velic_solB( coord, self->sigma, self->Z, self->n, self->km, NULL, NULL, stress, NULL );
}


void Velic_solB_StrainRateFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* strainRate ) {
	Velic_solB* self = (Velic_solB*) analyticSolution;
	
	_Velic_solB( coord, self->sigma, self->Z, self->n, self->km, NULL, NULL, NULL, strainRate );
}

void _Velic_solB_Init( Velic_solB* self, double sigma, double Z, double km, int n ) {
	Bool                     isCorrectInput = True;
	
	// Set the general AnalyticSolution (parent class) field function pointers to the specific analytic sols to use
        self->_getAnalyticVelocity    = Velic_solB_VelocityFunction;
	self->_getAnalyticPressure    = Velic_solB_PressureFunction;
	self->_getAnalyticTotalStress = Velic_solB_StressFunction;
	self->_getAnalyticStrainRate  = Velic_solB_StrainRateFunction;
	
	self->sigma = sigma;
	self->Z = Z;
	self->km = km;
	self->n = n;

	isCorrectInput = _checkInputParams( self );
	Journal_Firewall( isCorrectInput , Journal_Register( Error_Type, "Velic_solB" ),
			"Error in function %s: Bad Input parameters, solution check valid values in .tex documentation\n",
			__func__ );
}

void _Velic_solB_Build( void* analyticSolution, void* data ) {
	Velic_solB*           self  = (Velic_solB*)analyticSolution;

	AnalyticSolution_BuildAllAnalyticFields( self, data );

	_AnalyticSolution_Build( self, data );
}
	
void _Velic_solB_Construct( void* analyticSolution, Stg_ComponentFactory* cf, void* data ) {
	Velic_solB* self = (Velic_solB*) analyticSolution;
	FeVariable*              velocityField;
	FeVariable*              pressureField;
	FeVariable*              stressField;
	FeVariable*              strainRateField;
	FeVariable*              recoverdStrainRateField;
	FeVariable*              recoveredStressField;
	FeVariable*              recoveredPressureField;
	double                   sigma, Z, km, n;

	/* Construct Parent */
	_AnalyticSolution_Construct( self, cf, data );

	/* Create Analytic Fields */
	velocityField = Stg_ComponentFactory_ConstructByName( cf, "VelocityField", FeVariable, True, data );
	AnalyticSolution_RegisterFeVariableWithAnalyticFunction( self, velocityField, Velic_solB_VelocityFunction );

	pressureField = Stg_ComponentFactory_ConstructByName( cf, "PressureField", FeVariable, True, data );
	AnalyticSolution_RegisterFeVariableWithAnalyticFunction( self, pressureField, Velic_solB_PressureFunction );

	stressField = Stg_ComponentFactory_ConstructByName( cf, "StressField", FeVariable, False, data );
	if ( stressField ) 
		AnalyticSolution_RegisterFeVariableWithAnalyticFunction( self, stressField, Velic_solB_StressFunction );

	strainRateField = Stg_ComponentFactory_ConstructByName( cf, "StrainRateField", FeVariable, False, data );
	if ( strainRateField  )
		AnalyticSolution_RegisterFeVariableWithAnalyticFunction( self, strainRateField, Velic_solB_StrainRateFunction );

	recoverdStrainRateField = Stg_ComponentFactory_ConstructByName( cf, "recoveredStrainRateField", FeVariable, False, data );
	if ( recoverdStrainRateField )
		AnalyticSolution_RegisterFeVariableWithAnalyticFunction( self, recoverdStrainRateField, Velic_solB_StrainRateFunction );

	recoveredStressField = Stg_ComponentFactory_ConstructByName( cf, "recoveredStressField", FeVariable, False, data );
	if ( recoveredStressField )
		AnalyticSolution_RegisterFeVariableWithAnalyticFunction( self, recoveredStressField, Velic_solB_StressFunction );

	recoveredPressureField = Stg_ComponentFactory_ConstructByName( cf, "recoveredPressureField", FeVariable, False, data );
	if ( recoveredPressureField )
		AnalyticSolution_RegisterFeVariableWithAnalyticFunction( self, recoveredPressureField, Velic_solB_PressureFunction );

	sigma = Stg_ComponentFactory_GetRootDictDouble( cf, "solB_sigma", 1.0 );
	Z = Stg_ComponentFactory_GetRootDictDouble( cf, "solB_Z", 1.0 );
	km = Stg_ComponentFactory_GetRootDictDouble( cf, "wavenumberY", 2.0 );
	n = Stg_ComponentFactory_GetRootDictInt( cf, "wavenumberX", 1 );

	_Velic_solB_Init( self, sigma, Z, km, n );
}

Bool _checkInputParams( Velic_solB* self ) {
	return ( 
			( self->sigma > 0.0 ) && ( self->Z > 0.0 ) &&
			( self->km > 0.0 )    && ( self->n > 0 )  
		);
}
void* _Velic_solB_DefaultNew( Name name ) {
	return _AnalyticSolution_New(
			sizeof(Velic_solB),
			Velic_solB_Type,
			_AnalyticSolution_Delete,
			_AnalyticSolution_Print,
			_AnalyticSolution_Copy,
			_Velic_solB_DefaultNew,
			_Velic_solB_Construct,
			_Velic_solB_Build,
			_AnalyticSolution_Initialise,
			_AnalyticSolution_Execute,
			_AnalyticSolution_Destroy,
			name );
}

Index Underworld_Velic_solB_Register( PluginsManager* pluginsManager ) {
	return PluginsManager_Submit( pluginsManager, Velic_solB_Type, "0", _Velic_solB_DefaultNew );
}
