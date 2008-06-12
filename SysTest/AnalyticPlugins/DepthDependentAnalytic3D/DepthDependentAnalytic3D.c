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
*+		Robert Turnbull
*+		Vincent Lemiale
*+		Louis Moresi
*+		David May
*+		David Stegman
*+		Mirko Velic
*+		Patrick Sunter
*+		Julian Giordani
*+
** $Id: DepthDependentAnalytic3D.c 739 2008-06-12 06:05:14Z RobertTurnbull $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include <Underworld/Rheology/Rheology.h>

#include <assert.h>
#include <string.h>

const Type DepthDependentAnalytic3D_Type = "DepthDependentAnalytic3D";

typedef struct {
	__AnalyticSolution
	AnalyticSolution_SolutionFunction* viscosityFunc;
	AnalyticSolution_SolutionFunction* viscosityDerivativeFunc;
	AnalyticSolution_SolutionFunction* viscosity2ndDerivativeFunc;
	FeVariable*                        velocityField;
	FeVariable*                        pressureField;
	FeVariable*                        temperatureField;
	double                             Ra;
	double                             V0;
} DepthDependentAnalytic3D;

void DepthDependentAnalytic3D_ViscosityFunc_Isoviscous( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* eta ) {
	*eta = 1.0;
}
void DepthDependentAnalytic3D_ViscosityDerivativeFunc_Isoviscous( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* d_eta_dy ) {
	*d_eta_dy = 0.0;
}
void DepthDependentAnalytic3D_Viscosity2ndDerivativeFunc_Isoviscous( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* d2_eta_dy2 ) {
	*d2_eta_dy2 = 0.0;
}

void DepthDependentAnalytic3D_ViscosityFunc_Linear( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* eta ) {
	double y = coord[ J_AXIS ];
	*eta = 2.0 - y;
}
void DepthDependentAnalytic3D_ViscosityDerivativeFunc_Linear( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* d_eta_dy ) {
	*d_eta_dy = -1.0;
}
void DepthDependentAnalytic3D_Viscosity2ndDerivativeFunc_Linear( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* d2_eta_dy2 ) {
	*d2_eta_dy2 = 0.0;
}
void DepthDependentAnalytic3D_ViscosityFunc_Exponential( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* eta ) {
	double y = coord[ J_AXIS ];
	*eta = exp( 1.0 - y );
}
void DepthDependentAnalytic3D_ViscosityDerivativeFunc_Exponential( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* d_eta_dy ) {
	double y = coord[ J_AXIS ];
	*d_eta_dy = - exp( 1.0 - y );
}
void DepthDependentAnalytic3D_Viscosity2ndDerivativeFunc_Exponential( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* d2_eta_dy2 ) {
	double y = coord[ J_AXIS ];
	*d2_eta_dy2 = exp( 1.0 - y );
}
void DepthDependentAnalytic3D_ViscosityFunc_Exponential2( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* eta ) {
	double y = coord[ J_AXIS ];
	double viscosityContrast = 1e6;
	double gamma = log(viscosityContrast);
	*eta = viscosityContrast * exp( - gamma *( 1.0 - y ) );
}
void DepthDependentAnalytic3D_ViscosityDerivativeFunc_Exponential2( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* d_eta_dy ) {
	double y = coord[ J_AXIS ];
	double viscosityContrast = 1e6;
	double gamma = log(viscosityContrast);
	*d_eta_dy = viscosityContrast * gamma * exp( - gamma *( 1.0 - y ) );
}
void DepthDependentAnalytic3D_Viscosity2ndDerivativeFunc_Exponential2( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* d2_eta_dy2 ) {
	double y = coord[ J_AXIS ];
	double viscosityContrast = 1e6;
	double gamma = log(viscosityContrast);
	*d2_eta_dy2 = viscosityContrast * gamma * gamma * exp( - gamma *( 1.0 - y ) );
}


void _DepthDependentAnalytic3D_TemperatureFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* temperature ) {
	DepthDependentAnalytic3D* self            = (DepthDependentAnalytic3D*)analyticSolution;
	double                 x; 
	double                 y;
	double                 z;
	double                 V0                 = self->V0;
	double                 Ra                 = self->Ra;
	double                 perturbation;
	double                 eta, d_eta_dy, d2_eta_dy2;
	
	/* Find coordinate of node */
	x = coord[ I_AXIS ];
	y = coord[ J_AXIS ];
	z = coord[ K_AXIS ];
	
	/* Get Viscositiy and derivatives */
	self->viscosityFunc( self, NULL, coord, &eta );
	self->viscosityDerivativeFunc( self, NULL, coord, &d_eta_dy );
	self->viscosity2ndDerivativeFunc( self, NULL, coord, &d2_eta_dy2 );

	perturbation = V0 / Ra * cos( M_PI * x ) * cos( M_PI * z ) * 
		( 6.0 * M_PI * d_eta_dy * cos( M_PI * y ) - ( 9.0 * M_PI * M_PI * eta + d2_eta_dy2 ) * sin( M_PI * y ) );

	*temperature = 1 - y + perturbation;
}

void DepthDependentAnalytic3D_TemperatureIC( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	DomainContext*         context            = (DomainContext*)_context;
	FeVariable*            temperatureField   = (FeVariable*) FieldVariable_Register_GetByName( context->fieldVariable_Register, "TemperatureField" );
	FeMesh*                mesh               = temperatureField->feMesh;
	double*                temperature        = (double*) _result;
	double*                coord;
	DepthDependentAnalytic3D*  self           = Stg_ComponentFactory_ConstructByName( context->CF, DepthDependentAnalytic3D_Type, DepthDependentAnalytic3D, True, context );
	
	/* Find coordinate of node */
	coord = Mesh_GetVertex( mesh, node_lI );
	
	_DepthDependentAnalytic3D_TemperatureFunction( self, NULL, coord, temperature );
}

void _DepthDependentAnalytic3D_VelocityFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* velocity ) {
	DepthDependentAnalytic3D* self            = (DepthDependentAnalytic3D*)analyticSolution;
	double                 V0                 = self->V0;
	double                 x; 
	double                 y;
	double                 z;
	XYZ                    min, max;

	Mesh_GetGlobalCoordRange( self->velocityField->feMesh, min, max );
	x = coord[ I_AXIS ] - min[ I_AXIS ];
	y = coord[ J_AXIS ] - min[ J_AXIS ];
	z = coord[ K_AXIS ] - min[ K_AXIS ];

	velocity[ I_AXIS ] =          V0 * sin( M_PI * x ) * cos( M_PI * y ) * cos( M_PI * z ) ;
	velocity[ J_AXIS ] =  - 2.0 * V0 * cos( M_PI * x ) * sin( M_PI * y ) * cos( M_PI * z ) ;
	velocity[ K_AXIS ] =          V0 * cos( M_PI * x ) * cos( M_PI * y ) * sin( M_PI * z ) ;
}

void _DepthDependentAnalytic3D_ViscosityFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* viscosity ) {
	DepthDependentAnalytic3D* self            = (DepthDependentAnalytic3D*)analyticSolution;

	self->viscosityFunc( self, analyticFeVariable, coord, viscosity );
}

void _DepthDependentAnalytic3D_PressureFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* pressure ) {
	DepthDependentAnalytic3D* self            = (DepthDependentAnalytic3D*)analyticSolution;
	double                 x; 
	double                 y;
	double                 z;
	XYZ                    min, max;
	double                 V0                 = self->V0;
	double                 Ra                 = self->Ra;
	double                 C                  = 0;
	double                 f                  = 0;
	double                 eta, d_eta_dy;

	Mesh_GetGlobalCoordRange( self->velocityField->feMesh, min, max );
	x = coord[ I_AXIS ] - min[ I_AXIS ];
	y = coord[ J_AXIS ] - min[ J_AXIS ];
	z = coord[ K_AXIS ] - min[ K_AXIS ];

	f = - Ra*y*y*0.5 + Ra*y + C;
	
	/* Get Viscositiy and derivatives */
	self->viscosityFunc( self, NULL, coord, &eta );
	self->viscosityDerivativeFunc( self, NULL, coord, &d_eta_dy );

	*pressure = V0 * cos( M_PI * x ) * cos( M_PI * z ) * ( 3.0 * eta * M_PI * cos( M_PI * y ) - d_eta_dy * sin( M_PI * y ) ) + f;
}


void _DepthDependentAnalytic3D_Construct( void* analyticSolution, Stg_ComponentFactory* cf, void* data ) {
	DepthDependentAnalytic3D*         self           = (DepthDependentAnalytic3D*)analyticSolution;
	AbstractContext*        context;
	ConditionFunction*      condFunc;
	char*                   viscosityType;
	
	/* Construct Parent */
	_AnalyticSolution_Construct( self, cf, data );

	context = Stg_ComponentFactory_ConstructByName( cf, "context", AbstractContext, True, data ); 
	
	/* Add temperature initial condition */
	condFunc = ConditionFunction_New( DepthDependentAnalytic3D_TemperatureIC, "DepthDependentAnalytic3D_TemperatureIC" );
	ConditionFunction_Register_Add( context->condFunc_Register, condFunc );
	
	/* Create Analytic Fields */
	self->velocityField = AnalyticSolution_RegisterFeVariableFromCF( self, "VelocityField", _DepthDependentAnalytic3D_VelocityFunction, cf, True, data );
	self->pressureField = AnalyticSolution_RegisterFeVariableFromCF( self, "PressureField", _DepthDependentAnalytic3D_PressureFunction, cf, True, data );
	self->temperatureField = AnalyticSolution_RegisterFeVariableFromCF( self, "TemperatureField", _DepthDependentAnalytic3D_TemperatureFunction, cf, True, data );
	AnalyticSolution_RegisterFeVariableFromCF( self, "ViscosityField", _DepthDependentAnalytic3D_ViscosityFunction, cf, False, data );

	/* Setup Viscosity Functions */
	viscosityType = Stg_ComponentFactory_GetRootDictString( cf, "ViscosityType", "Isoviscous" );
	if ( strcasecmp( viscosityType, "Isoviscous" ) == 0 ) {
		self->viscosityFunc              = DepthDependentAnalytic3D_ViscosityFunc_Isoviscous;
		self->viscosityDerivativeFunc    = DepthDependentAnalytic3D_ViscosityDerivativeFunc_Isoviscous;
		self->viscosity2ndDerivativeFunc = DepthDependentAnalytic3D_Viscosity2ndDerivativeFunc_Isoviscous;
	}
	else if ( strcasecmp( viscosityType, "Linear" ) == 0 ) {
		self->viscosityFunc              = DepthDependentAnalytic3D_ViscosityFunc_Linear;
		self->viscosityDerivativeFunc    = DepthDependentAnalytic3D_ViscosityDerivativeFunc_Linear;
		self->viscosity2ndDerivativeFunc = DepthDependentAnalytic3D_Viscosity2ndDerivativeFunc_Linear;
	}
	else if ( strcasecmp( viscosityType, "Exponential" ) == 0 ) {
		self->viscosityFunc              = DepthDependentAnalytic3D_ViscosityFunc_Exponential;
		self->viscosityDerivativeFunc    = DepthDependentAnalytic3D_ViscosityDerivativeFunc_Exponential;
		self->viscosity2ndDerivativeFunc = DepthDependentAnalytic3D_Viscosity2ndDerivativeFunc_Exponential;
	}
	else if ( strcasecmp( viscosityType, "Exponential2" ) == 0 ) {
		self->viscosityFunc              = DepthDependentAnalytic3D_ViscosityFunc_Exponential2;
		self->viscosityDerivativeFunc    = DepthDependentAnalytic3D_ViscosityDerivativeFunc_Exponential2;
		self->viscosity2ndDerivativeFunc = DepthDependentAnalytic3D_Viscosity2ndDerivativeFunc_Exponential2;
	}

	self->Ra = Stg_ComponentFactory_GetRootDictDouble( cf, "Ra", 0.0 );
	self->V0 = Stg_ComponentFactory_GetRootDictDouble( cf, "V0", 0.0 );
}

void _DepthDependentAnalytic3D_Build( void* analyticSolution, void* data ) {
	DepthDependentAnalytic3D*         self = (DepthDependentAnalytic3D*)analyticSolution;

	_AnalyticSolution_Build( self, data );
}


void* _DepthDependentAnalytic3D_DefaultNew( Name name ) {
	return _AnalyticSolution_New(
			sizeof(DepthDependentAnalytic3D),
			DepthDependentAnalytic3D_Type,
			_AnalyticSolution_Delete,
			_AnalyticSolution_Print,
			_AnalyticSolution_Copy,
			_DepthDependentAnalytic3D_DefaultNew,
			_DepthDependentAnalytic3D_Construct,
			_DepthDependentAnalytic3D_Build,
			_AnalyticSolution_Initialise,
			_AnalyticSolution_Execute,
			_AnalyticSolution_Destroy,
			name );
}

Index Underworld_DepthDependentAnalytic3D_Register( PluginsManager* pluginsManager ) {
	return PluginsManager_Submit( pluginsManager, DepthDependentAnalytic3D_Type, "0", _DepthDependentAnalytic3D_DefaultNew );
}
