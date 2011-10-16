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
** $Id: solCA.c 636 2006-09-08 03:07:06Z JulianGiordani $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include "Analytic_solCA.h"

const Type Velic_solCA_Type = "Underworld_Velic_solCA";

void Velic_solCA_PressureFunction( void* analyticSolution, FeVariable* analyticFeVariable, const double *coord, double* pressure ) {
	Velic_solCA* self = (Velic_solCA*) analyticSolution;
	
	_Velic_solCA( coord, self->sigma, self->eta, self->dx, self->x0, NULL, pressure, NULL, NULL );
}

void Velic_solCA_VelocityFunction( void* analyticSolution, FeVariable* analyticFeVariable, const double *coord, double* velocity ) {
	Velic_solCA* self = (Velic_solCA*) analyticSolution;
	
	_Velic_solCA( coord, self->sigma, self->eta, self->dx, self->x0, velocity, NULL, NULL, NULL );
}

void Velic_solCA_StressFunction( void* analyticSolution, FeVariable* analyticFeVariable, const double *coord, double* stress ) {
	Velic_solCA* self = (Velic_solCA*) analyticSolution;
	
	_Velic_solCA( coord, self->sigma, self->eta, self->dx, self->x0, NULL, NULL, stress, NULL );
}


void Velic_solCA_StrainRateFunction( void* analyticSolution, FeVariable* analyticFeVariable, const double *coord, double* strainRate ) {
	Velic_solCA* self = (Velic_solCA*) analyticSolution;
	
	_Velic_solCA( coord, self->sigma, self->eta, self->dx, self->x0, NULL, NULL, NULL, strainRate );
}

void _Velic_solCA_Init( Velic_solCA* self, double sigma, double eta, double dx, double x0 ) {
	Bool                     isCorrectInput = True;
	// Set the general AnalyticSolution (parent class) field function pointers to the specific analytic sols to use
        self->_getAnalyticVelocity    = Velic_solCA_VelocityFunction;
	self->_getAnalyticPressure    = Velic_solCA_PressureFunction;
	self->_getAnalyticTotalStress = Velic_solCA_StressFunction;
	self->_getAnalyticStrainRate  = Velic_solCA_StrainRateFunction;

	self->sigma = sigma;
	self->eta   = eta;
	self->dx    = dx;
	self->x0    = x0;

	isCorrectInput = _checkInputParams( self );
	Journal_Firewall( isCorrectInput , Journal_Register( Error_Type, (Name)"Velic_solCA"  ),
			"Error in function %s: Bad Input parameters, solution check valid values in .tex documentation\n",
			__func__ );
}

void _Velic_solCA_AssignFromXML( void* analyticSolution, Stg_ComponentFactory* cf, void* data ) {
	Velic_solCA* self = (Velic_solCA*) analyticSolution;
	FeVariable*              velocityField;
	FeVariable*              pressureField;
	FeVariable*              stressField;
	FeVariable*              strainRateField;
	FeVariable*              recoveredStrainRateField;
	FeVariable*              recoveredStressField;
	double                   sigma, eta, dx, x0, startX, endX;

	/* Construct Parent */
	_AnalyticSolution_AssignFromXML( self, cf, data );

	/* Create Analytic Fields */
	velocityField = Stg_ComponentFactory_ConstructByName( cf, (Name)"VelocityField", FeVariable, True, data  );
	AnalyticSolution_RegisterFeVariableWithAnalyticFunction( self, velocityField, Velic_solCA_VelocityFunction );

	pressureField = Stg_ComponentFactory_ConstructByName( cf, (Name)"PressureField", FeVariable, True, data  );
	AnalyticSolution_RegisterFeVariableWithAnalyticFunction( self, pressureField, Velic_solCA_PressureFunction );

	stressField = Stg_ComponentFactory_ConstructByName( cf, (Name)"StressField", FeVariable, False, data );
	if ( stressField  )
		AnalyticSolution_RegisterFeVariableWithAnalyticFunction( self, stressField, Velic_solCA_StressFunction );

	strainRateField = Stg_ComponentFactory_ConstructByName( cf, (Name)"StrainRateField", FeVariable, False, data );
	if ( strainRateField   ) {
		AnalyticSolution_RegisterFeVariableWithAnalyticFunction( self, strainRateField, Velic_solCA_StrainRateFunction );
	}

	recoveredStrainRateField = Stg_ComponentFactory_ConstructByName( cf, (Name)"recoveredStrainRateField", FeVariable, False, data );
	if ( recoveredStrainRateField  )
		AnalyticSolution_RegisterFeVariableWithAnalyticFunction( self, recoveredStrainRateField, Velic_solCA_StrainRateFunction );

	recoveredStressField = Stg_ComponentFactory_ConstructByName( cf, (Name)"recoveredStressField", FeVariable, False, data );
	if ( recoveredStressField  )
		AnalyticSolution_RegisterFeVariableWithAnalyticFunction( self, recoveredStressField, Velic_solCA_StressFunction );

	sigma = Stg_ComponentFactory_GetRootDictDouble( cf, (Dictionary_Entry_Key)"solCA_sigma", 1.0  );
	eta = Stg_ComponentFactory_GetRootDictDouble( cf, (Dictionary_Entry_Key)"solCA_eta", 1.0  );
	startX = Stg_ComponentFactory_GetRootDictDouble( cf, (Dictionary_Entry_Key)"solCA_startX", 0.4  );
	endX   = Stg_ComponentFactory_GetRootDictDouble( cf, (Dictionary_Entry_Key)"solCA_endX", 0.8 );

	dx = endX - startX;
	x0 = 0.5 * (endX + startX );
	
	_Velic_solCA_Init( self, sigma, eta, dx, x0 );

}

Bool _checkInputParams( Velic_solCA* self ) {
	return ( 
			( self->sigma > 0.0 ) && ( self->eta > 0.0 ) &&
			( self->dx >= 0.0 )    && ( self->x0 > 0 ) &&
			( self->dx < 1.0 )    && ( self->x0 < 1.0 )
		);
}

void* _Velic_solCA_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                              _sizeOfSelf = sizeof(Velic_solCA);
	Type                                                      type = Velic_solCA_Type;
	Stg_Class_DeleteFunction*                              _delete = _AnalyticSolution_Delete;
	Stg_Class_PrintFunction*                                _print = _AnalyticSolution_Print;
	Stg_Class_CopyFunction*                                  _copy = _AnalyticSolution_Copy;
	Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _Velic_solCA_DefaultNew;
	Stg_Component_ConstructFunction*                    _construct = _Velic_solCA_AssignFromXML;
	Stg_Component_BuildFunction*                            _build = _AnalyticSolution_Build;
	Stg_Component_InitialiseFunction*                  _initialise = _AnalyticSolution_Initialise;
	Stg_Component_ExecuteFunction*                        _execute = _AnalyticSolution_Execute;
	Stg_Component_DestroyFunction*                        _destroy = _AnalyticSolution_Destroy;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

	return _AnalyticSolution_New(  ANALYTICSOLUTION_PASSARGS  );
}

Index Underworld_Velic_solCA_Register( PluginsManager* pluginsManager ) {
	return PluginsManager_Submit( pluginsManager, Velic_solCA_Type, (Name)"0", _Velic_solCA_DefaultNew  );
}


