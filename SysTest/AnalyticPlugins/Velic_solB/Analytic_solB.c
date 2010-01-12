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

#include "Analytic_solB.h"

const Type Underworld_solB_Type = "Underworld_Velic_solB";

void Underworld_solB_PressureFunction( void* analyticSolution, double* coord, double* pressure ) {
	Underworld_solB* self = (Underworld_solB*) analyticSolution;
	
	_Velic_solB( coord, self->sigma, self->Z, self->n, self->km, NULL, pressure, NULL, NULL );
}

void Underworld_solB_VelocityFunction( void* analyticSolution, double* coord, double* velocity ) {
	Underworld_solB* self = (Underworld_solB*) analyticSolution;
	
	_Velic_solB( coord, self->sigma, self->Z, self->n, self->km, velocity, NULL, NULL, NULL );
}

void Underworld_solB_StressFunction( void* analyticSolution, double* coord, double* stress ) {
	Underworld_solB* self = (Underworld_solB*) analyticSolution;
	
	_Velic_solB( coord, self->sigma, self->Z, self->n, self->km, NULL, NULL, stress, NULL );
}


void Underworld_solB_StrainRateFunction( void* analyticSolution, double* coord, double* strainRate ) {
	Underworld_solB* self = (Underworld_solB*) analyticSolution;
	
	_Velic_solB( coord, self->sigma, self->Z, self->n, self->km, NULL, NULL, NULL, strainRate );
}

void _Underworld_solB_Init( Underworld_solB* self, double sigma, double Z, double wavenumberY, int n ) {
	self->sigma = sigma;
	self->Z     = Z;
	self-> km   = wavenumberY;
	self->n     = n;	
}

void _Underworld_solB_AssignFromXML( void* analyticSolution, Stg_ComponentFactory* cf, void* data ) {
	Underworld_solB* 		self = (Underworld_solB*) analyticSolution;
	Bool                     isCorrectInput = True;
	double                   sigma, Z, wavenumberY, n;

	/* Construct Parent */
	_FieldTest_AssignFromXML( self, cf, data );

	sigma = Stg_ComponentFactory_GetRootDictDouble( cf, "solB_sigma", 1.0 );
	Z = Stg_ComponentFactory_GetRootDictDouble( cf, "solB_Z", 1.0 );
	wavenumberY = Stg_ComponentFactory_GetRootDictDouble( cf, "wavenumberY", 1.0 );
	n = Stg_ComponentFactory_GetRootDictInt( cf, "wavenumberX", 1 );
	
	_Underworld_solB_Init( self, sigma, Z, wavenumberY, n );

	isCorrectInput = solB_checkInputParams( self );
	Journal_Firewall( isCorrectInput , Journal_Register( Error_Type, "Analytic_solB" ),
			"Error in function %s: Bad Input parameters, solution check valid values in .tex documentation\n",
			__func__ );
}

void _Underworld_solB_Build( void* analyticSolution, void* data ) {
	Underworld_solB* 		self = (Underworld_solB*) analyticSolution;

	_FieldTest_Build( self, data );

	/* here we assign the memory and the func ptr for analytic sols */
	self->_analyticSolutionList = Memory_Alloc_Array_Unnamed( FieldTest_AnalyticSolutionFunc*, 4 );
	/* this order MUST be consistent with the xml file definition */
	self->_analyticSolutionList[0] = Underworld_solB_VelocityFunction;
	self->_analyticSolutionList[1] = Underworld_solB_PressureFunction;
	self->_analyticSolutionList[2] = Underworld_solB_StrainRateFunction;
	self->_analyticSolutionList[3] = Underworld_solB_StressFunction;
}

Bool solB_checkInputParams( Underworld_solB* self ) {
	return ( 
			( self->sigma > 0.0 ) && ( self->Z > 0.0 ) &&
			( self->km > 0.0 )    && ( self->n > 0 )  
		);
}
void* _Underworld_solB_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                              _sizeOfSelf = sizeof(Underworld_solB);
	Type                                                      type = Underworld_solB_Type;
	Stg_Class_DeleteFunction*                              _delete = _FieldTest_Delete;
	Stg_Class_PrintFunction*                                _print = _FieldTest_Print;
	Stg_Class_CopyFunction*                                  _copy = _FieldTest_Copy;
	Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _Underworld_solB_DefaultNew;
	Stg_Component_ConstructFunction*                    _construct = _Underworld_solB_AssignFromXML;
	Stg_Component_BuildFunction*                            _build = _Underworld_solB_Build;
	Stg_Component_InitialiseFunction*                  _initialise = _FieldTest_Initialise;
	Stg_Component_ExecuteFunction*                        _execute = _FieldTest_Execute;
	Stg_Component_DestroyFunction*                        _destroy = _FieldTest_Destroy;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

	return _FieldTest_New(  FIELDTEST_PASSARGS  );
}

Index Underworld_Velic_solB_Register( PluginsManager* pluginsManager ) {
	return PluginsManager_Submit( pluginsManager, Underworld_solB_Type, "0", _Underworld_solB_DefaultNew );
}


