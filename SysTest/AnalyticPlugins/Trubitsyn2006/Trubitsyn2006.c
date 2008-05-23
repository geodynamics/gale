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
** $Id: Trubitsyn2006.c 737 2008-05-23 06:08:42Z RobertTurnbull $
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

#define SMALL 1.0e-5
#define IS_ODD(A)   ( (A) % 2 == 1 )

const Type Trubitsyn2006_Type = "Trubitsyn2006";

typedef struct {
	__AnalyticSolution
	AnalyticSolution_SolutionFunction* viscosityFunc;
	AnalyticSolution_SolutionFunction* viscosityDerivativeFunc;
	FeVariable*                        velocityField;
	FeVariable*                        pressureField;
	double                             Ra;
	double                             T0;
	int                                wavenumberX;
	int                                wavenumberY;
} Trubitsyn2006;

/** Analytic Solution taken from 
 * V. P. Trubitsyn, A. A. Baranov, A. N. Eyseev, and A. P. Trubitsyn. 
 * 	Exact analytical solutions of the stokes equation for testing the equations of mantle convection with a variable viscosity. 
 *	Izvestiya Physics of the Solid Earth, 42(7):537–545, 2006. 
 *  
 *  All equations refer to this paper */

double Trubitsyn_G( double x, double x0, double a ) {
	return 1/(1 + exp( -2.0 * a * ( x - x0 ) ) );
}
double Trubitsyn_GDeriv( double x, double x0, double a ) {
	return 2.0 * a * pow( exp( a * ( x - x0 ) ) + exp( -a * ( x - x0 ) ), -2.0 );
}

void _Trubitsyn2006_ViscosityFunc_Isoviscous( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* eta ) {
	*eta = 1.0;
}
void _Trubitsyn2006_ViscosityDerivativeFunc_Isoviscous( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* viscDeriv ) {
	viscDeriv[0] = 0.0;
	viscDeriv[1] = 0.0;
}
void _Trubitsyn2006_ViscosityFunc_Model1( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* eta ) {
	*eta = 1.0 + 100.0 * coord[ I_AXIS ] * coord[ J_AXIS ];
}
void _Trubitsyn2006_ViscosityDerivativeFunc_Model1( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* viscDeriv ) {
	viscDeriv[0] = 100.0 * coord[ J_AXIS ];
	viscDeriv[1] = 100.0 * coord[ I_AXIS ];
}
void _Trubitsyn2006_ViscosityFunc_Model2( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* eta ) {
//	*eta = 999.0 * ( 1.0 - Trubitsyn_G( coord[ I_AXIS ], 0.5, 50 ) ) + 1;
	*eta = ( coord[ I_AXIS ] <= 0.5 ? 1000.0 : 1 );
}
void _Trubitsyn2006_ViscosityDerivativeFunc_Model2( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* viscDeriv ) {
//	viscDeriv[0] = -999.0 * Trubitsyn_GDeriv( coord[ I_AXIS ], 0.5, 50 );
	viscDeriv[0] = 0.0;
	viscDeriv[1] = 0.0;
}
void _Trubitsyn2006_ViscosityFunc_Model3( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* eta ) {
	//*eta = 999.0 * ( 1.0 - Trubitsyn_G( coord[ J_AXIS ], 0.5, 50 ) ) + 1;
	*eta = ( coord[ J_AXIS ] <= 0.5 ? 1000.0 : 1 );
}
void _Trubitsyn2006_ViscosityDerivativeFunc_Model3( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* viscDeriv ) {
	viscDeriv[0] = 0.0;
	viscDeriv[1] = 0.0;
//	viscDeriv[1] = -999.0 * Trubitsyn_GDeriv( coord[ J_AXIS ], 0.5, 50 );
}


double Trubitsyn2006_V0( void* analyticSolution ) {
	Trubitsyn2006*         self               = (Trubitsyn2006*)analyticSolution;
	double                 T0                 = self->T0;
	double                 Ra                 = self->Ra;

	/* NB - Paper Trubitsyn [2006] is in error for this substitution - they leave out the square on pi */
	return 0.25 * Ra * T0 * M_1_PI * M_1_PI;
}

void Trubitsyn2006_TemperatureIC( Node_LocalIndex node_lI, Variable_Index var_I, void* _context, void* _result ) {
	DomainContext*         context            = (DomainContext*)_context;
	FeVariable*            temperatureField   = (FeVariable*) FieldVariable_Register_GetByName( context->fieldVariable_Register, "TemperatureField" );
	FeMesh*                mesh               = temperatureField->feMesh;
	double*                temperature        = (double*) _result;
	double*                coord;
	double                 x; 
	double                 y;
	Trubitsyn2006*         self               = Stg_ComponentFactory_ConstructByName( context->CF, Trubitsyn2006_Type, Trubitsyn2006, True, context );
	double                 T0                 = self->T0;
	double                 Ra                 = self->Ra;
	double                 v0                 = Trubitsyn2006_V0( self );
	XYZ                    min, max;
	double                 eta;
	double                 d_eta_dy;
	XYZ                    viscDeriv;
	
	/* Find coordinate of node */
	coord = Mesh_GetVertex( mesh, node_lI );
	Mesh_GetGlobalCoordRange( mesh, min, max );
	x = coord[ I_AXIS ] - min[ I_AXIS ];
	y = coord[ J_AXIS ] - min[ J_AXIS ];

	self->viscosityFunc( self, NULL, coord, &eta );
	self->viscosityDerivativeFunc( self, NULL, coord, viscDeriv );
	d_eta_dy = viscDeriv[ J_AXIS ];

	*temperature = 1 - y + T0 * cos( M_PI * x ) * ( eta * sin( M_PI * y ) - d_eta_dy / M_PI * cos( M_PI * y ) ) ;    /* Equation 47 */ 
	
	*temperature = 1 - y + 4.0 * M_PI * v0 / Ra * cos( M_PI * x ) * 
		( M_PI * eta * sin( M_PI * y ) - d_eta_dy * cos( M_PI * y ) ) ;    /* Equation 47 */ 

}

void _Trubitsyn2006_VelocityFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* velocity ) {
	Trubitsyn2006*         self               = (Trubitsyn2006*)analyticSolution;
	double                 v0                 = Trubitsyn2006_V0( self );
	double                 x; 
	double                 y;
	XYZ                    min, max;

	Mesh_GetGlobalCoordRange( self->velocityField->feMesh, min, max );
	x = coord[ I_AXIS ] - min[ I_AXIS ];
	y = coord[ J_AXIS ] - min[ J_AXIS ];

	velocity[ J_AXIS ] =   v0 * sin( M_PI * y ) * cos( M_PI * x );         /* Equation 31 */
	velocity[ I_AXIS ] = - v0 * cos( M_PI * y ) * sin( M_PI * x );         /* Equation 32 */
}

void _Trubitsyn2006_ViscosityFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* viscosity ) {
	Trubitsyn2006*         self               = (Trubitsyn2006*)analyticSolution;

	self->viscosityFunc( self, analyticFeVariable, coord, viscosity );
}


void _Trubitsyn2006_PressureFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* pressure ) {
	Trubitsyn2006*         self               = (Trubitsyn2006*)analyticSolution;
	double                 T0                 = self->T0;
	double                 Ra                 = self->Ra;
	double                 x; 
	double                 y;
	XYZ                    min, max;
	double                 eta;
	double                 eta01;
	double                 upperLeftCorner[]  = { 0.0, 1.0, 0.0 };
	
	Mesh_GetGlobalCoordRange( self->velocityField->feMesh, min, max );
	x = coord[ I_AXIS ] - min[ I_AXIS ];
	y = coord[ J_AXIS ] - min[ J_AXIS ];

	/* Get Viscosities */
	self->viscosityFunc( self, NULL, coord, &eta );
	self->viscosityFunc( self, NULL, upperLeftCorner, &eta01 );

	*pressure = - Ra * ( ( y * y * 0.5 - y + 0.5 )
			+ 0.5 * T0 / M_PI * ( eta * cos( M_PI * y ) * cos( M_PI * x ) + eta01 ) );        /* Equation 46 */
	//printf("pressure from t0 = %g\n", *pressure );
	//*pressure = - 2.0 * M_PI * v0 * ( eta * cos( M_PI * y ) * cos( M_PI * x ) + eta01 ) 
	//		- Ra * ( y*y * 0.5 - y + 0.5 );
	//printf("pressure from v0 = %g\n\n", *pressure );
}

void _Trubitsyn2006_StreamFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* psi ) {
	Trubitsyn2006*         self               = (Trubitsyn2006*)analyticSolution;
	double                 v0                 = Trubitsyn2006_V0( self );
	double                 x; 
	double                 y;
	XYZ                    min, max;

	Mesh_GetGlobalCoordRange( self->velocityField->feMesh, min, max );
	x = coord[ I_AXIS ] - min[ I_AXIS ];
	y = coord[ J_AXIS ] - min[ J_AXIS ];

	*psi = - v0 / M_PI * sin( M_PI * y ) * sin( M_PI * x ) ;                                          /* Equation 40 */
}


void _Trubitsyn2006_Construct( void* analyticSolution, Stg_ComponentFactory* cf, void* data ) {
	Trubitsyn2006*         self           = (Trubitsyn2006*)analyticSolution;
	AbstractContext*        context;
	ConditionFunction*      condFunc;
	char*                   viscosityType;
	
	/* Construct Parent */
	_AnalyticSolution_Construct( self, cf, data );

	context = Stg_ComponentFactory_ConstructByName( cf, "context", AbstractContext, True, data ); 
	
	/* Add temperature initial condition */
	condFunc = ConditionFunction_New( Trubitsyn2006_TemperatureIC, "Trubitsyn2006_TemperatureIC" );
	ConditionFunction_Register_Add( context->condFunc_Register, condFunc );
	
	/* Setup Viscosity Functions */
	viscosityType = Stg_ComponentFactory_GetRootDictString( cf, "ViscosityType", "Isoviscous" );
	if ( strcasecmp( viscosityType, "Isoviscous" ) == 0 ) {
		self->viscosityFunc           = _Trubitsyn2006_ViscosityFunc_Isoviscous;
		self->viscosityDerivativeFunc = _Trubitsyn2006_ViscosityDerivativeFunc_Isoviscous;
	}
	else if ( strcasecmp( viscosityType, "Model1" ) == 0 ) {
		self->viscosityFunc           = _Trubitsyn2006_ViscosityFunc_Model1;
		self->viscosityDerivativeFunc = _Trubitsyn2006_ViscosityDerivativeFunc_Model1;
	}
	else if ( strcasecmp( viscosityType, "Model2" ) == 0 ) {
		self->viscosityFunc           = _Trubitsyn2006_ViscosityFunc_Model2;
		self->viscosityDerivativeFunc = _Trubitsyn2006_ViscosityDerivativeFunc_Model2;
	}
	else if ( strcasecmp( viscosityType, "Model3" ) == 0 ) {
		self->viscosityFunc           = _Trubitsyn2006_ViscosityFunc_Model3;
		self->viscosityDerivativeFunc = _Trubitsyn2006_ViscosityDerivativeFunc_Model3;
	}
	else {
		Journal_Printf( Journal_Register( Error_Type, self->type ), "Cannot understand viscosity type = '%s'\n", viscosityType );
		abort();
	}
	
	/* Create Analytic Fields */
	self->velocityField = AnalyticSolution_RegisterFeVariableFromCF( self, "VelocityField",  _Trubitsyn2006_VelocityFunction, cf, True, data );
	self->pressureField = AnalyticSolution_RegisterFeVariableFromCF( self, "PressureField",  _Trubitsyn2006_PressureFunction, cf, True, data );
	AnalyticSolution_RegisterFeVariableFromCF( self, "ViscosityField", _Trubitsyn2006_ViscosityFunction, cf, False, data );

	self->Ra = Stg_ComponentFactory_GetRootDictDouble( cf, "Ra", 0.0 );
	self->T0 = Stg_ComponentFactory_GetRootDictDouble( cf, "T0", 0.0 );
	self->wavenumberX = Stg_ComponentFactory_GetRootDictInt( cf, "wavenumberX", 1 );
	self->wavenumberY = Stg_ComponentFactory_GetRootDictInt( cf, "wavenumberY", 1 );
}

void _Trubitsyn2006_Build( void* analyticSolution, void* data ) {
	Trubitsyn2006*         self = (Trubitsyn2006*)analyticSolution;

	_AnalyticSolution_Build( self, data );
}


void* _Trubitsyn2006_DefaultNew( Name name ) {
	return _AnalyticSolution_New(
			sizeof(Trubitsyn2006),
			Trubitsyn2006_Type,
			_AnalyticSolution_Delete,
			_AnalyticSolution_Print,
			_AnalyticSolution_Copy,
			_Trubitsyn2006_DefaultNew,
			_Trubitsyn2006_Construct,
			_Trubitsyn2006_Build,
			_AnalyticSolution_Initialise,
			_AnalyticSolution_Execute,
			_AnalyticSolution_Destroy,
			name );
}

Index Underworld_Trubitsyn2006_Register( PluginsManager* pluginsManager ) {
	return PluginsManager_Submit( pluginsManager, Trubitsyn2006_Type, "0", _Trubitsyn2006_DefaultNew );
}
