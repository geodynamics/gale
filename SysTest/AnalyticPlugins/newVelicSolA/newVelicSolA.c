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
** $Id: solA.c 636 2006-09-08 03:07:06Z JulianGiordani $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include "newVelicSolA.h"

const Type newVelicSolA_Type = "Analytic_newVelicSolA";

void newVelicSolA_PressureFunction( void* analyticSolution, double* coord, double* pressure ) {
	newVelicSolA* self = (newVelicSolA*) analyticSolution;
	
	_Velic_solA( coord, self->sigma, self->Z, self->n, self->km, NULL, pressure, NULL, NULL );
}

void newVelicSolA_VelocityFunction( void* analyticSolution, double* coord, double* velocity ) {
	newVelicSolA* self = (newVelicSolA*) analyticSolution;
	
	_Velic_solA( coord, self->sigma, self->Z, self->n, self->km, velocity, NULL, NULL, NULL );
}

void newVelicSolA_StressFunction( void* analyticSolution, double* coord, double* stress ) {
	newVelicSolA* self = (newVelicSolA*) analyticSolution;
	
	_Velic_solA( coord, self->sigma, self->Z, self->n, self->km, NULL, NULL, stress, NULL );
}


void newVelicSolA_StrainRateFunction( void* analyticSolution, double* coord, double* strainRate ) {
	newVelicSolA* self = (newVelicSolA*) analyticSolution;
	
	_Velic_solA( coord, self->sigma, self->Z, self->n, self->km, NULL, NULL, NULL, strainRate );
}

void _newVelicSolA_Init( newVelicSolA* self, double sigma, double Z, double wavenumberY, int n ) {
	self->sigma = sigma;
	self->Z     = Z;
	self-> km   = M_PI * wavenumberY;
	self->n     = n;	
}

void _newVelicSolA_Construct( void* analyticSolution, Stg_ComponentFactory* cf, void* data ) {
	newVelicSolA* 		self = (newVelicSolA*) analyticSolution;
	/*
	FeVariable*              velocityField;
	FeVariable*              pressureField;
	FeVariable*              stressField;
	FeVariable*              strainRateField;
	FeVariable*              recoverdStrainRateField;
	FeVariable*              recoveredStressField;
	*/
	Bool                     isCorrectInput = True;
	double                   sigma, Z, wavenumberY, n;
	

	/* Construct Parent */
	_FieldTest_Construct( self, cf, data );

	/* Create Analytic Fields */
	/*
	velocityField = Stg_ComponentFactory_ConstructByName( cf, "VelocityField", FeVariable, True, data );
	FieldTest_RegisterFeVariableWithAnalyticFunction( self, velocityField, newVelicSolA_VelocityFunction );

	pressureField = Stg_ComponentFactory_ConstructByName( cf, "PressureField", FeVariable, True, data );
	FieldTest_RegisterFeVariableWithAnalyticFunction( self, pressureField, newVelicSolA_PressureFunction );

	stressField = Stg_ComponentFactory_ConstructByName( cf, "StressField", FeVariable, False, data );
	if ( stressField )
		FieldTest_RegisterFeVariableWithAnalyticFunction( self, stressField, newVelicSolA_StressFunction );

	strainRateField = Stg_ComponentFactory_ConstructByName( cf, "StrainRateField", FeVariable, False, data );
	if ( strainRateField  ) {
		FieldTest_RegisterFeVariableWithAnalyticFunction( self, strainRateField, newVelicSolA_StrainRateFunction );
	}

	recoverdStrainRateField = Stg_ComponentFactory_ConstructByName( cf, "recoveredStrainRate", FeVariable, False, data );
	if ( recoverdStrainRateField )
		FieldTest_RegisterFeVariableWithAnalyticFunction( self, recoverdStrainRateField, newVelicSolA_StrainRateFunction );

	recoveredStressField = Stg_ComponentFactory_ConstructByName( cf, "recoveredStressField", FeVariable, False, data );
	if ( recoveredStressField )
		FieldTest_RegisterFeVariableWithAnalyticFunction( self, recoveredStressField, newVelicSolA_StressFunction );
	*/

	sigma = Stg_ComponentFactory_GetRootDictDouble( cf, "solA_sigma", 1.0 );
	Z = Stg_ComponentFactory_GetRootDictDouble( cf, "solA_Z", 1.0 );
	wavenumberY = Stg_ComponentFactory_GetRootDictDouble( cf, "wavenumberY", 1.0 );
	n = Stg_ComponentFactory_GetRootDictInt( cf, "wavenumberX", 1 );
	
	_newVelicSolA_Init( self, sigma, Z, wavenumberY, n );

	isCorrectInput = _checkInputParams( self );
	Journal_Firewall( isCorrectInput , Journal_Register( Error_Type, "newVelicSolA" ),
			"Error in function %s: Bad Input parameters, solution check valid values in .tex documentation\n",
			__func__ );
}

void _newVelicSolA_Build( void* analyticSolution, void* data ) {
	newVelicSolA* 		self = (newVelicSolA*) analyticSolution;

	_FieldTest_Build( self, data );

	FieldTest_AddAnalyticSolutionFuncToListAtIndex( self, 0, newVelicSolA_VelocityFunction, 0 );
	FieldTest_AddAnalyticSolutionFuncToListAtIndex( self, 1, newVelicSolA_PressureFunction, 1 );
}

Bool _checkInputParams( newVelicSolA* self ) {
	return ( 
			( self->sigma > 0.0 ) && ( self->Z > 0.0 ) &&
			( self->km > 0.0 )    && ( self->n > 0 )  
		);
}
void* _newVelicSolA_DefaultNew( Name name ) {
	return _FieldTest_New(
			sizeof(newVelicSolA),
			newVelicSolA_Type,
			_FieldTest_Delete,
			_FieldTest_Print,
			_FieldTest_Copy,
			_newVelicSolA_DefaultNew,
			_newVelicSolA_Construct,
			_newVelicSolA_Build,
			_FieldTest_Initialise,
			_FieldTest_Execute,
			_FieldTest_Destroy,
			name );
}

Index Underworld_newVelicSolA_Register( PluginsManager* pluginsManager ) {
	return PluginsManager_Submit( pluginsManager, newVelicSolA_Type, "0", _newVelicSolA_DefaultNew );
}
