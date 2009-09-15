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
** $Id: solS.c 636 2006-09-08 03:07:06Z JulianGiordani $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <stdio.h>
#include <string.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include "Analytic_solS.h"

const Type Velic_solS_Type = "Underworld_Velic_solS";

void Velic_solS_PressureFunction( void* analyticSolution, double* coord, double* pressure ) {
	Velic_solS* self = (Velic_solS*) analyticSolution;
	
	_Velic_solS( coord, self->_n, self->_eta, NULL, pressure, NULL, NULL );
}

void Velic_solS_VelocityFunction( void* analyticSolution, double* coord, double* velocity ) {
	Velic_solS* self = (Velic_solS*) analyticSolution;
	
	_Velic_solS( coord, self->_n, self->_eta, velocity, NULL, NULL, NULL );
}

void Velic_solS_StressFunction( void* analyticSolution, double* coord, double* stress ) {
	Velic_solS* self = (Velic_solS*) analyticSolution;
	
	_Velic_solS( coord, self->_n, self->_eta, NULL, NULL, stress, NULL );
}


void Velic_solS_StrainRateFunction( void* analyticSolution, double* coord, double* strainRate ) {
	Velic_solS* self = (Velic_solS*) analyticSolution;
	
	_Velic_solS( coord, self->_n, self->_eta, NULL, NULL, NULL, strainRate );
}

void _Velic_solS_Init( Velic_solS* self, double _eta, int _n ) {
	Bool                     isCorrectInput = True;

	self->_eta = _eta;
	self->_n = _n;
	
	isCorrectInput = solS_checkInputParams( self );
	Journal_Firewall( isCorrectInput , Journal_Register( Error_Type, "Velic_solS" ),
			"Error in function %s: Bad Input parameters, solution check valid values in .tex documentation\n",
			__func__ );
}

void _Velic_solS_Build( void* analyticSolution, void* data ) {
	Velic_solS*          self  = (Velic_solS*)analyticSolution;
	
	_FieldTest_Build( self, data );

	/* here we assign the memory and the func ptr for analytic sols */
	self->_analyticSolutionList = Memory_Alloc_Array_Unnamed( FieldTest_AnalyticSolutionFunc*, 4 );
	/* this order MUST be consistent with the xml file definition */
	self->_analyticSolutionList[0] = Velic_solS_VelocityFunction;
	self->_analyticSolutionList[1] = Velic_solS_PressureFunction;
	self->_analyticSolutionList[2] = Velic_solS_StrainRateFunction;
	self->_analyticSolutionList[3] = Velic_solS_StressFunction;
}

void _Velic_solS_Construct( void* analyticSolution, Stg_ComponentFactory* cf, void* data ) {
	Velic_solS* self = (Velic_solS*) analyticSolution;
	double                   _eta;
	int                      _n;

	/* Construct Parent */
	_FieldTest_Construct( self, cf, data );
	
	_eta = Stg_ComponentFactory_GetRootDictDouble( cf, "solS_eta", 1.0 );
	_n = Stg_ComponentFactory_GetRootDictInt( cf, "sinusoidalLidWavenumber", 1 );

        _Velic_solS_Init( self, _eta, _n );
}

Bool solS_checkInputParams( Velic_solS* self ) {
	return ( 
			( self->_eta > 0.0 ) && ( self->_n > 0.0 )
		);
}
void* _Velic_solS_DefaultNew( Name name ) {
	return  _FieldTest_New(
			sizeof(Velic_solS),
			Velic_solS_Type,
			_FieldTest_Delete,
			_FieldTest_Print,
			_FieldTest_Copy,
			_Velic_solS_DefaultNew,
			_Velic_solS_Construct,
			_Velic_solS_Build,
			_FieldTest_Initialise,
			_FieldTest_Execute,
			_FieldTest_Destroy,
			name );
}

Index Underworld_Velic_solS_Register( PluginsManager* pluginsManager ) {
	return PluginsManager_Submit( pluginsManager, Velic_solS_Type, "0", _Velic_solS_DefaultNew );
}
