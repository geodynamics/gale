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

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include "Analytic_solKz.h"

const Type Velic_solKz_Type = "Underworld_Velic_solKz";

void Velic_solKz_PressureFunction( void* analyticSolution, double* coord, double* pressure ) {
	Velic_solKz* self = (Velic_solKz*) analyticSolution;
	
	_Velic_solKz( coord, self->sigma, self->km, self->n, self->B, NULL, pressure, NULL, NULL );
}

void Velic_solKz_VelocityFunction( void* analyticSolution, double* coord, double* velocity ) {
	Velic_solKz* self = (Velic_solKz*) analyticSolution;
	
	_Velic_solKz( coord, self->sigma, self->km, self->n, self->B, velocity, NULL, NULL, NULL );
}

void Velic_solKz_StressFunction( void* analyticSolution, double* coord, double* stress ) {
	Velic_solKz* self = (Velic_solKz*) analyticSolution;
	
	_Velic_solKz( coord, self->sigma, self->km, self->n, self->B, NULL, NULL, stress, NULL );
}


void Velic_solKz_StrainRateFunction( void* analyticSolution, double* coord, double* strainRate ) {
	Velic_solKz* self = (Velic_solKz*) analyticSolution;
	
	_Velic_solKz( coord, self->sigma, self->km, self->n, self->B, NULL, NULL, NULL, strainRate );
}

void Velic_solKz_ViscosityFunction( void* analyticSolution, double* coord, double* viscosity ) {
	Velic_solKz* self = (Velic_solKz*) analyticSolution;
	
	*viscosity = exp( 2.0 * self->B * coord[ J_AXIS ] );
}

	
void _Velic_solKz_Init( Velic_solKz* self, double sigma, double km, double B, int n ) {
	self->sigma = sigma;
	self->km = km;
	self->B = B;
	self->n = n;
}

void _Velic_solKz_Build( void* analyticSolution, void* data ) {
	Velic_solKz*          self  = (Velic_solKz*)analyticSolution;

	_FieldTest_Build( self, data );

	/* here we assign the memory and the func ptr for analytic sols */
	self->_analyticSolutionList = Memory_Alloc_Array_Unnamed( FieldTest_AnalyticSolutionFunc*, 4 );
	/* this order MUST be consistent with the xml file definition */
	self->_analyticSolutionList[0] = Velic_solKz_VelocityFunction;
	self->_analyticSolutionList[1] = Velic_solKz_PressureFunction;
	self->_analyticSolutionList[2] = Velic_solKz_StrainRateFunction;
	self->_analyticSolutionList[3] = Velic_solKz_StressFunction;
}

void _Velic_solKz_Construct( void* analyticSolution, Stg_ComponentFactory* cf, void* data ) {
	Velic_solKz* self = (Velic_solKz*) analyticSolution;
	double                   sigma, km, B, twiceB;
	int                      n;

	/* Construct Parent */
	_FieldTest_Construct( self, cf, data );

	sigma = Stg_ComponentFactory_GetRootDictDouble( cf, "solKz_sigma", 1.0 );
	km = Stg_ComponentFactory_GetRootDictDouble( cf, "solKz_km", M_PI );
	twiceB = Stg_ComponentFactory_GetRootDictDouble( cf, "solKz_twiceB", 2.0 );
	B = Stg_ComponentFactory_GetRootDictDouble( cf, "solKz_B", 0.5 * twiceB );
	n = Stg_ComponentFactory_GetRootDictInt( cf, "solKz_n", 1 );
	
	_Velic_solKz_Init( self, sigma, km, B, n );
}

void* _Velic_solKz_DefaultNew( Name name ) {
	return _FieldTest_New(
			sizeof(Velic_solKz),
			Velic_solKz_Type,
			_FieldTest_Delete,
			_FieldTest_Print,
			_FieldTest_Copy,
			_Velic_solKz_DefaultNew,
			_Velic_solKz_Construct,
			_Velic_solKz_Build,
			_FieldTest_Initialise,
			_FieldTest_Execute,
			_FieldTest_Destroy,
			name );
}

Index Underworld_Velic_solKz_Register( PluginsManager* pluginsManager ) {
	return PluginsManager_Submit( pluginsManager, Velic_solKz_Type, "0", _Velic_solKz_DefaultNew );
}
