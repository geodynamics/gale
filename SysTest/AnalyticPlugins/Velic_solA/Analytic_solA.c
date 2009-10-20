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

#include "Analytic_solA.h"

const Type Underworld_solA_Type = "Underworld_Velic_solA";

void Underworld_solA_PressureFunction( void* analyticSolution, double* coord, double* pressure ) {
	Underworld_solA* self = (Underworld_solA*) analyticSolution;
	
	_Velic_solutionA( coord, self->sigma, self->Z, self->n, self->km, NULL, pressure, NULL, NULL );
}

void Underworld_solA_VelocityFunction( void* analyticSolution, double* coord, double* velocity ) {
	Underworld_solA* self = (Underworld_solA*) analyticSolution;
	
	_Velic_solutionA( coord, self->sigma, self->Z, self->n, self->km, velocity, NULL, NULL, NULL );
}

void Underworld_solA_StressFunction( void* analyticSolution, double* coord, double* stress ) {
	Underworld_solA* self = (Underworld_solA*) analyticSolution;
	
	_Velic_solutionA( coord, self->sigma, self->Z, self->n, self->km, NULL, NULL, stress, NULL );
}


void Underworld_solA_StrainRateFunction( void* analyticSolution, double* coord, double* strainRate ) {
	Underworld_solA* self = (Underworld_solA*) analyticSolution;
	
	_Velic_solutionA( coord, self->sigma, self->Z, self->n, self->km, NULL, NULL, NULL, strainRate );
}

void _Underworld_solA_Init( Underworld_solA* self, double sigma, double Z, double wavenumberY, int n ) {
	self->sigma = sigma;
	self->Z     = Z;
	self-> km   = M_PI * wavenumberY;
	self->n     = n;	
}

void _Underworld_solA_Construct( void* analyticSolution, Stg_ComponentFactory* cf, void* data ) {
	Underworld_solA*	self = (Underworld_solA*) analyticSolution;
	Bool					isCorrectInput = True;
	double				sigma, Z, wavenumberY, n;
	
	/* Construct Parent */
	_FieldTest_Construct( self, cf, data );

	sigma = Stg_ComponentFactory_GetRootDictDouble( cf, "solA_sigma", 1.0 );
	Z = Stg_ComponentFactory_GetRootDictDouble( cf, "solA_Z", 1.0 );
	wavenumberY = Stg_ComponentFactory_GetRootDictDouble( cf, "wavenumberY", 1.5 );
	n = Stg_ComponentFactory_GetRootDictInt( cf, "wavenumberX", 1 );
	
	_Underworld_solA_Init( self, sigma, Z, wavenumberY, n );

	isCorrectInput = solA_checkInputParams( self );
	Journal_Firewall( isCorrectInput , Journal_Register( Error_Type, "solA" ),
			"Error in function %s: Bad Input parameters, solution check valid values in .tex documentation\n",
			__func__ );
}

void _Underworld_solA_Build( void* analyticSolution, void* data ) {
	Underworld_solA* 		self = (Underworld_solA*) analyticSolution;

	_FieldTest_Build( self, data );

	/* here we assign the memory and the func ptr for analytic sols */
	self->_analyticSolutionList = Memory_Alloc_Array_Unnamed( FieldTest_AnalyticSolutionFunc*, 4 );
	/* this order MUST be consistent with the xml file definition */
	self->_analyticSolutionList[0] = Underworld_solA_VelocityFunction;
	self->_analyticSolutionList[1] = Underworld_solA_PressureFunction;
	self->_analyticSolutionList[2] = Underworld_solA_StrainRateFunction;
	self->_analyticSolutionList[3] = Underworld_solA_StressFunction;

}

Bool solA_checkInputParams( Underworld_solA* self ) {
	return ( 
		( self->sigma > 0.0 ) && ( self->Z > 0.0 ) &&
		( self->km > 0.0 )    && ( self->n > 0 )  
	);
}

void* _Underworld_solA_DefaultNew( Name name ) {
	return _FieldTest_New(
			sizeof(Underworld_solA),
			Underworld_solA_Type,
			_FieldTest_Delete,
			_FieldTest_Print,
			_FieldTest_Copy,
			_Underworld_solA_DefaultNew,
			_Underworld_solA_Construct,
			_Underworld_solA_Build,
			_FieldTest_Initialise,
			_FieldTest_Execute,
			_FieldTest_Destroy,
			name );
}

Index Underworld_Velic_solA_Register( PluginsManager* pluginsManager ) {
	return PluginsManager_Submit( pluginsManager, Underworld_solA_Type, "0", _Underworld_solA_DefaultNew );
}
