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
** $Id: solCx.c 636 2006-09-08 03:07:06Z JulianGiordani $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
#include<mpi.h>
#include<stdio.h>
#include<string.h>
#include<StGermain/StGermain.h>
#include<StgFEM/StgFEM.h>

#include "Analytic_solCx.h"

const Type ExperimentalUnderworld_Velic_solCx_Type = "ExperimentalUnderworld_Velic_solCx";

void Velic_solCx_PressureFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* pressure ) {
	Velic_solCx* self = (Velic_solCx*) analyticSolution;
	
	_Velic_solCx( coord, self->etaA, self->etaB, self->xc, self->n, NULL, pressure, NULL, NULL );
}

void Velic_solCx_VelocityFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* velocity ) {
	Velic_solCx* self = (Velic_solCx*) analyticSolution;
	
	_Velic_solCx( coord, self->etaA, self->etaB, self->xc, self->n, velocity, NULL, NULL, NULL );
}

void Velic_solCx_StressFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* stress ) {
	Velic_solCx* self = (Velic_solCx*) analyticSolution;
	
	_Velic_solCx( coord, self->etaA, self->etaB, self->xc, self->n, NULL, NULL, stress, NULL );
}


void Velic_solCx_StrainRateFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* strainRate ) {
	Velic_solCx* self = (Velic_solCx*) analyticSolution;
	
	_Velic_solCx( coord, self->etaA, self->etaB, self->xc, self->n, NULL, NULL, NULL, strainRate );
}

void _Velic_solCx_Init( Velic_solCx* self, double etaA, double etaB, double xc, int n ) {
	Bool                     correctInput = True;

	self->_getAnalyticVelocity = Velic_solCx_VelocityFunction;
	self->_getAnalyticPressure = Velic_solCx_PressureFunction;
	self->_getAnalyticTotalStress = Velic_solCx_StressFunction;
	self->_getAnalyticStrainRate = Velic_solCx_StrainRateFunction;

	self->etaA = etaA;
	self->etaB = etaB;
	self->xc   = xc;
	self->n    = n;

	if( (self->xc < 0.0) || self->xc > 1.0 ) 
	      correctInput = False;	

	Journal_Firewall( correctInput,  Journal_Register( Error_Type, "solCx_ErrorStream") ,
			"Error: In func %s. The input parameters you supplied to the analytic Solution are incorrect.\nValid range of xc values is [0,1]. Currently it is %f\n", __func__, self->xc );
}

void _Velic_solCx_Build( void* analyticSolution, void* data ) {
	Velic_solCx*          self  = (Velic_solCx*)analyticSolution;

	AnalyticSolution_BuildAllAnalyticFields( self, data );

	_AnalyticSolution_Build( self, data );
}

void _Velic_solCx_Construct( void* analyticSolution, Stg_ComponentFactory* cf, void* data ) {
	Velic_solCx* self = (Velic_solCx*) analyticSolution;
	FeVariable*              velocityField;
	FeVariable*              pressureField;
	FeVariable*              recoveredPressureField;
	FeVariable*              stressField;
	FeVariable*              strainRateField;
	FeVariable*              recoverdStrainRateField;
	FeVariable*              recoveredStressField;
	double                   etaA,  etaB,  xc;
	int                      n;

	/* Construct Parent */
	_AnalyticSolution_Construct( self, cf, data );

	/* Create Analytic Fields */
	velocityField = Stg_ComponentFactory_ConstructByName( cf, "VelocityField", FeVariable, True, data );
	AnalyticSolution_RegisterFeVariableWithAnalyticFunction( self, velocityField, Velic_solCx_VelocityFunction );

	pressureField = Stg_ComponentFactory_ConstructByName( cf, "PressureField", FeVariable, True, data );
	AnalyticSolution_RegisterFeVariableWithAnalyticFunction( self, pressureField, Velic_solCx_PressureFunction );

	stressField = Stg_ComponentFactory_ConstructByName( cf, "StressField", FeVariable, False, data );
	if ( stressField ) 
		AnalyticSolution_RegisterFeVariableWithAnalyticFunction( self, stressField, Velic_solCx_StressFunction );

	strainRateField = Stg_ComponentFactory_ConstructByName( cf, "StrainRateField", FeVariable, False, data );
	if ( strainRateField  )
		AnalyticSolution_RegisterFeVariableWithAnalyticFunction( self, strainRateField, Velic_solCx_StrainRateFunction );

	recoverdStrainRateField = Stg_ComponentFactory_ConstructByName( cf, "recoveredStrainRateField", FeVariable, False, data );
	if ( recoverdStrainRateField )
		AnalyticSolution_RegisterFeVariableWithAnalyticFunction( self, recoverdStrainRateField, Velic_solCx_StrainRateFunction );

	recoveredStressField = Stg_ComponentFactory_ConstructByName( cf, "recoveredStressField", FeVariable, False, data );
	if ( recoveredStressField )
		AnalyticSolution_RegisterFeVariableWithAnalyticFunction( self, recoveredStressField, Velic_solCx_StressFunction );

	recoveredPressureField = Stg_ComponentFactory_ConstructByName( cf, "recoveredPressureField", FeVariable, False, data );
	if ( recoveredPressureField )
		AnalyticSolution_RegisterFeVariableWithAnalyticFunction( self, recoveredPressureField, Velic_solCx_PressureFunction );
	
	etaA = Stg_ComponentFactory_GetRootDictDouble( cf, "solCx_etaA", 1.0 );
	etaB = Stg_ComponentFactory_GetRootDictDouble( cf, "solCx_etaB", 2.0 );
	xc   = Stg_ComponentFactory_GetRootDictDouble( cf, "solCx_xc", 0.25 );
	n    = Stg_ComponentFactory_GetRootDictInt( cf, "solCx_n", 1 );

	_Velic_solCx_Init( self,  etaA,  etaB,  xc, n );
}

void* _Velic_solCx_DefaultNew( Name name ) {
	return _AnalyticSolution_New(
			sizeof(Velic_solCx),
			ExperimentalUnderworld_Velic_solCx_Type,
			_AnalyticSolution_Delete,
			_AnalyticSolution_Print,
			_AnalyticSolution_Copy,
			_Velic_solCx_DefaultNew,
			_Velic_solCx_Construct,
			_Velic_solCx_Build,
			_AnalyticSolution_Initialise,
			_AnalyticSolution_Execute,
			_AnalyticSolution_Destroy,
			name );
}

Index ExperimentalUnderworld_Velic_solCx_Register( PluginsManager* pluginsManager ) {
	return PluginsManager_Submit( pluginsManager, ExperimentalUnderworld_Velic_solCx_Type, "0", _Velic_solCx_DefaultNew );
}
