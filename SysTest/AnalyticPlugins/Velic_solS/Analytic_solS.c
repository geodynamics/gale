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

#include "Analytic_solS.h"

const Type Underworld_solS_Type = "Underworld_Velic_solS";

void Underworld_solS_PressureFunction( void* analyticSolution, double* coord, double* pressure ) {
	Underworld_solS* self = (Underworld_solS*) analyticSolution;
	
	_Velic_solS( coord, self->_n, self->_eta, NULL, pressure, NULL, NULL );
}

void Underworld_solS_VelocityFunction( void* analyticSolution, double* coord, double* velocity ) {
	Underworld_solS* self = (Underworld_solS*) analyticSolution;
	
	_Velic_solS( coord, self->_n, self->_eta, velocity, NULL, NULL, NULL );
}

void Underworld_solS_StressFunction( void* analyticSolution, double* coord, double* stress ) {
	Underworld_solS* self = (Underworld_solS*) analyticSolution;
	
	_Velic_solS( coord, self->_n, self->_eta, NULL, NULL, stress, NULL );
}

void Underworld_solS_StrainRateFunction( void* analyticSolution, double* coord, double* strainRate ) {
	Underworld_solS* self = (Underworld_solS*) analyticSolution;
	
	_Velic_solS( coord, self->_n, self->_eta, NULL, NULL, NULL, strainRate );
}

void _Underworld_solS_Init( Underworld_solS* self, double _eta, int _n ) {
	Bool isCorrectInput = True;

	self->_eta = _eta;
	self->_n = _n;
	
	isCorrectInput = _checkInputParams( self );
	Journal_Firewall( isCorrectInput , Journal_Register( Error_Type, "Underworld_solS" ),
			"Error in function %s: Bad Input parameters, solution check valid values in .tex documentation\n",
			__func__ );
}

void _Underworld_solS_Build( void* analyticSolution, void* data ) {
	Underworld_solS* self = (Underworld_solS*)analyticSolution;
	
	_FieldTest_Build( self, data );

	/* here we assign the memory and the func ptr for analytic sols */
	self->_analyticSolutionList = Memory_Alloc_Array_Unnamed( FieldTest_AnalyticSolutionFunc*, 4 );
	/* this order MUST be consistent with the xml file definition */
	self->_analyticSolutionList[0] = Underworld_solS_VelocityFunction;
	self->_analyticSolutionList[1] = Underworld_solS_PressureFunction;
	self->_analyticSolutionList[2] = Underworld_solS_StrainRateFunction;
	self->_analyticSolutionList[3] = Underworld_solS_StressFunction;
}

void _Underworld_solS_AssignFromXML( void* analyticSolution, Stg_ComponentFactory* cf, void* data ) {
	Underworld_solS*	self = (Underworld_solS*) analyticSolution;
	Bool					isCorrectInput = True;
	double				_eta;
	int					_n;

	/* Construct Parent */
	_FieldTest_AssignFromXML( self, cf, data );
	
	_eta = Stg_ComponentFactory_GetRootDictDouble( cf, "solS_eta", 1.0 );
	_n = Stg_ComponentFactory_GetRootDictInt( cf, "sinusoidalLidWavenumber", 1 );

	_Underworld_solS_Init( self, _eta, _n );

	isCorrectInput = _checkInputParams( self );
	Journal_Firewall( isCorrectInput , Journal_Register( Error_Type, "solS" ),
			"Error in function %s: Bad Input parameters, solution check valid values in .tex documentation\n",
			__func__ );
}

Bool _checkInputParams( Underworld_solS* self ) {
	return ( 
		( self->_eta > 0.0 ) && ( self->_n > 0.0 )
	);
}

void* _Underworld_solS_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                              _sizeOfSelf = sizeof(Underworld_solS);
	Type                                                      type = Underworld_solS_Type;
	Stg_Class_DeleteFunction*                              _delete = _FieldTest_Delete;
	Stg_Class_PrintFunction*                                _print = _FieldTest_Print;
	Stg_Class_CopyFunction*                                  _copy = _FieldTest_Copy;
	Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _Underworld_solS_DefaultNew;
	Stg_Component_ConstructFunction*                    _construct = _Underworld_solS_AssignFromXML;
	Stg_Component_BuildFunction*                            _build = _Underworld_solS_Build;
	Stg_Component_InitialiseFunction*                  _initialise = _FieldTest_Initialise;
	Stg_Component_ExecuteFunction*                        _execute = _FieldTest_Execute;
	Stg_Component_DestroyFunction*                        _destroy = _FieldTest_Destroy;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = ZERO;

	return  _FieldTest_New(  FIELDTEST_PASSARGS  );
}

Index Underworld_Velic_solS_Register( PluginsManager* pluginsManager ) {
	return PluginsManager_Submit( pluginsManager, Underworld_solS_Type, "0", _Underworld_solS_DefaultNew );
}


