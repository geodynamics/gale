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
** $Id: solA.c 822 2007-04-27 06:20:35Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgFEM/StgFEM.h>
#include "solA.h"

const Type Velic_solA_Type = "Velic_solA";
/*
typedef struct {
	__AnalyticSolution
	double sigma;
	double Z;
} Velic_solA;
*/

double Calculate__C1 ( Velic_solA* self, double x, double z, double y, double n ) {
	double t1, t3, t4, t5, t10, t17, t21, t24;
	double sigma = self->sigma;
	double Z = self->Z;
	double km = M_PI;
	double kn = n * M_PI;

	t1 = sin(km);
	t3 = kn * kn;
	t4 = exp(-kn);
	t5 = t4 * t4;
	t10 = km * km;
	t17 = pow(t10 + t3, 0.2e1);
	t21 = pow(t4 - 0.1e1, 0.2e1);
	t24 = pow(t4 + 0.1e1, 0.2e1);

	return t1 * sigma * (t3 * t5 + t3 - 0.2e1 * kn * t5 + 0.2e1 * kn + t10 * t5 + t10) * t4 / Z / t17 / t21 / t24 / 0.2e1;
}

double Calculate__C2 ( Velic_solA* self, double x, double z, double y, double n ) {
	double t1, t3, t4, t5, t10, t16, t20, t23;
	double sigma = self->sigma;
	double Z = self->Z;
	double km = M_PI;

	double kn = n * M_PI;

	t1 = sin(km);
	t3 = kn * kn;
	t4 = exp(-kn);
	t5 = t4 * t4;
	t10 = km * km;
	t16 = pow(t10 + t3, 0.2e1);
	t20 = pow(t4 - 0.1e1, 0.2e1);
	t23 = pow(t4 + 0.1e1, 0.2e1);

	return -t1 * sigma * (t3 * t5 + t3 - 0.2e1 * kn * t5 + 0.2e1 * kn + t10 * t5 + t10) / Z / t16 / t20 / t23 / 0.2e1;
}

double Calculate__C3 ( Velic_solA* self, double x, double z, double y, double n ) {
	double t1, t3, t6, t7;
	double sigma = self->sigma;
	double Z = self->Z;
	double km = M_PI;

	double kn = n * M_PI;

	t1 = sin(km);
	t3 = exp(-kn);
	t6 = km * km;
	t7 = kn * kn;

	return -t1 * sigma * t3 / Z / (t6 + t7) / (t3 - 0.1e1) / (t3 + 0.1e1) / 0.2e1;
}

double Calculate__C4 ( Velic_solA* self, double x, double z, double y, double n ) {
	double t1, t5, t6, t9;
	double sigma = self->sigma;
	double Z = self->Z;
	double km = M_PI;

	double kn = n * M_PI;

	t1 = sin(km);
	t5 = km * km;
	t6 = kn * kn;
	t9 = exp(-kn);

	return -t1 * sigma / Z / (t5 + t6) / (t9 - 0.1e1) / (t9 + 0.1e1) / 0.2e1;
}

double Calculate_ss ( Velic_solA* self, double x, double z, double y, double n ) {
	double t4, t10, t14, t15, t16, t18;
	double sigma = self->sigma;
	double Z = self->Z;
	double km = M_PI;

	double kn = n * M_PI;

	t4 = exp(-kn * z);
	t10 = exp(kn * (z - 0.1e1));
	t14 = sin(km * z);
	t15 = km * km;
	t16 = kn * kn;
	t18 = pow(t15 + t16, 0.2e1);

	return (Calculate__C1( self, x, z, y, n ) + z * Calculate__C3( self, x, z, y, n )) * t4 + (Calculate__C2( self, x, z, y, n ) + z * Calculate__C4( self, x, z, y, n )) * t10 + kn * sigma * t14 / t18 / Z;
}

double Calculate_ss_z ( Velic_solA* self, double x, double z, double y, double n ) {
	double t6, t14, t18, t20, t21, t23;
	double sigma = self->sigma;
	double Z = self->Z;
	double km = M_PI;

	double kn = n * M_PI;

	t6 = exp(-kn * z);
	t14 = exp(kn * (z - 0.1e1));
	t18 = cos(km * z);
	t20 = km * km;
	t21 = kn * kn;
	t23 = pow(t20 + t21, 0.2e1);

	return (-(Calculate__C1( self, x, z, y, n ) + z * Calculate__C3( self, x, z, y, n )) * kn + Calculate__C3( self, x, z, y, n )) * t6 + (Calculate__C4( self, x, z, y, n ) + (Calculate__C2( self, x, z, y, n ) + z * Calculate__C4( self, x, z, y, n )) * kn) * t14 + kn * sigma * t18 * km / t23 / Z;
}

double Calculate_ss_zz ( Velic_solA* self, double x, double z, double y, double n ) {
	double t3, t9, t19, t23, t25, t27;
	double sigma = self->sigma;
	double Z = self->Z;
	double km = M_PI;

	double kn = n * M_PI;

	t3 = kn * kn;
	t9 = exp(-kn * z);
	t19 = exp(kn * (z - 0.1e1));
	t23 = sin(km * z);
	t25 = km * km;
	t27 = pow(t25 + t3, 0.2e1);

	return ((Calculate__C1( self, x, z, y, n ) + z * Calculate__C3( self, x, z, y, n )) * t3 - 0.2e1 * Calculate__C3( self, x, z, y, n ) * kn) * t9 + (0.2e1 * Calculate__C4( self, x, z, y, n ) * kn + (Calculate__C2( self, x, z, y, n ) + z * Calculate__C4( self, x, z, y, n )) * t3) * t19 - kn * sigma * t23 * t25 / t27 / Z;
}

double Calculate_ss_zzz ( Velic_solA* self, double x, double z, double y, double n ) {
	double t3, t4, t10, t20, t24, t26, t29;
	double sigma = self->sigma;
	double Z = self->Z;
	double km = M_PI;

	double kn = n * M_PI;

	t3 = kn * kn;
	t4 = t3 * kn;
	t10 = exp(-kn * z);
	t20 = exp(kn * (z - 0.1e1));
	t24 = cos(km * z);
	t26 = km * km;
	t29 = pow(t26 + t3, 0.2e1);

	return (-(Calculate__C1( self, x, z, y, n ) + z * Calculate__C3( self, x, z, y, n )) * t4 + 0.3e1 * Calculate__C3( self, x, z, y, n ) * t3) * t10 + (0.3e1 * Calculate__C4( self, x, z, y, n ) * t3 + (Calculate__C2( self, x, z, y, n ) + z * Calculate__C4( self, x, z, y, n )) * t4) * t20 - kn * sigma * t24 * t26 * km / t29 / Z;
}




void Velic_solA_PressureFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* pressure ) {
	Velic_solA* self = (Velic_solA*) analyticSolution;
	double x    = coord[ I_AXIS ];
	double z    = coord[ J_AXIS ];
	double y    = coord[ K_AXIS ];
	double n    = 1;
	double kn   = n * M_PI;
	double ss_zzz, ss_z, pp;

	ss_zzz = Calculate_ss_zzz( self, x, z, y, n );
	ss_z   = Calculate_ss_z( self, x, z, y, n );
	pp     = self->Z*(ss_zzz-kn*kn*ss_z)/kn;
	
	*pressure = pp * cos(kn * x);
}

void Velic_solA_VelocityFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* velocity ) {
	Velic_solA* self = (Velic_solA*) analyticSolution;
	double x    = coord[ I_AXIS ];
	double z    = coord[ J_AXIS ];
	double y    = coord[ K_AXIS ];
	double n    = 1;
	double u1, u2;


	u1 = n * M_PI * Calculate_ss( self, x, z, y, n );
	u2 = -Calculate_ss_z( self, x, z, y, n );

	u1 *= cos( n * M_PI * x );
	u2 *= sin( n * M_PI * x );

	velocity[ I_AXIS ] = u2;
	velocity[ J_AXIS ] = u1;
	if ( analyticFeVariable->dim == 3 ) 
		velocity[ K_AXIS ] = 0.0;
}

void Velic_solA_StressFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* stress ) {
	Velic_solA* self = (Velic_solA*) analyticSolution;
	double x    = coord[ I_AXIS ];
	double z    = coord[ J_AXIS ];
	double y    = coord[ K_AXIS ];
	double n    = 1;
	double kn   = n * M_PI;
	double pp, txx, txz, tzz;
	double ss_zz, ss_z, ss;
	double mean;

	ss_zz = Calculate_ss_zz( self, x, z, y, n );
	ss_z   = Calculate_ss_z( self, x, z, y, n );
	ss     = Calculate_ss( self, x, z, y, n );
	
	pp  = self->Z*(ss_zz-kn*kn*ss_z)/kn;
	tzz = 2.0*kn*ss_z - pp;
	txz = -self->Z*(ss_zz + kn*kn*ss);
	txx = -2.0*self->Z*kn*ss_z - pp;

	tzz *= cos(n*M_PI*x);
	txz *= sin(n*M_PI*x);
	txx *= cos(n*M_PI*x);
/*  only deviatoric stress is calculated */
	mean = 0.5 * ( txx + tzz );
	stress[ 0 ] = txx - mean;
	stress[ 1 ] = tzz - mean;
	if ( analyticFeVariable->dim == 2 ) {
		stress[ 2 ] = txz;
	}
	else { 
		stress[ 2 ] = 0.0;
		stress[ 3 ] = txz;
		stress[ 4 ] = 0.0;
		stress[ 5 ] = 0.0;
	}
}


void Velic_solA_StrainRateFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* strainRate ) {
	Velic_solA* self = (Velic_solA*) analyticSolution;
	double x    = coord[ I_AXIS ];
	double z    = coord[ J_AXIS ];
	double y    = coord[ K_AXIS ];
	double n    = 1;
	double kn   = n * M_PI;
	double e_xx, e_xz, e_zz;
	double ss_zz, ss_z, ss;

	ss_zz = Calculate_ss_zz( self, x, z, y, n );
	ss_z   = Calculate_ss_z( self, x, z, y, n );
	ss     = Calculate_ss( self, x, z, y, n );
	
	e_zz = kn*ss_z*cos(kn*x);  /* zz rate of strain */
	e_xx = -e_zz;  /* xx rate of strain */
	e_xz = -0.5*(ss_zz+kn*kn*ss)*sin(kn*x); /* xz rate of strain */

	strainRate[ 0 ] = e_xx;
	strainRate[ 1 ] = e_zz;
	if ( analyticFeVariable->dim == 2 ) {
		strainRate[ 2 ] = e_xz;
	}
	else { 
		strainRate[ 2 ] = 0.0;
		strainRate[ 3 ] = e_xz;
		strainRate[ 4 ] = 0.0;
		strainRate[ 5 ] = 0.0;
	}
}

void _Velic_solA_Construct( void* analyticSolution, Stg_ComponentFactory* cf, void* data ) {
	Velic_solA* self = (Velic_solA*) analyticSolution;

	/* Construct Parent */
	_AnalyticSolution_Construct( self, cf, data );

	/* Create Analytic Fields */
	self->velocityField = Stg_ComponentFactory_ConstructByName( cf, "VelocityField", FeVariable, True, data ); 
	self->pressureField = Stg_ComponentFactory_ConstructByName( cf, "PressureField", FeVariable, True, data ); 
	self->stressField = Stg_ComponentFactory_ConstructByName( cf, "StressField", FeVariable, False, data ); 
	self->strainRateField = Stg_ComponentFactory_ConstructByName( cf, "StrainRateField", FeVariable, False, data ); 
	self->recoveredStrainRateField = Stg_ComponentFactory_ConstructByName( cf, "recoveredStrainRate", FeVariable, False, data ); 
	self->recoveredStressField = Stg_ComponentFactory_ConstructByName( cf, "recoveredStress", FeVariable, False, data );

	self->sigma = Stg_ComponentFactory_GetRootDictDouble( cf, "sigma", 1.0 );
	self->Z = Stg_ComponentFactory_GetRootDictDouble( cf, "Z", 1.0 );
}

void _Velic_solA_Build( void* analyticSolution, void* data ) {
	Velic_solA* self = (Velic_solA*) analyticSolution;

	Build( self->velocityField, data, False );
	Build( self->pressureField, data, False );

	AnalyticSolution_CreateAnalyticVectorField( self, self->velocityField, Velic_solA_VelocityFunction );
	AnalyticSolution_CreateAnalyticField( self, self->pressureField, Velic_solA_PressureFunction );
	if ( self->stressField ) {
		Build( self->stressField, data, False );
		AnalyticSolution_CreateAnalyticSymmetricTensorField( self, self->stressField, Velic_solA_StressFunction );
	}
	if ( self->strainRateField  ) {
		Build( self->strainRateField, data, False );
		AnalyticSolution_CreateAnalyticSymmetricTensorField( self, self->strainRateField, Velic_solA_StrainRateFunction );
	}
	if ( self->recoveredStrainRateField ) {
		Build( self->recoveredStrainRateField, data, False );
		AnalyticSolution_CreateAnalyticSymmetricTensorField( self, self->recoveredStrainRateField, Velic_solA_StrainRateFunction );
	}
	if ( self->recoveredStressField ) {
		Build( self->recoveredStressField, data, False );
		AnalyticSolution_CreateAnalyticSymmetricTensorField( self, self->recoveredStressField, Velic_solA_StressFunction );
	}

	_AnalyticSolution_Build( self, data );
}

void* _Velic_solA_DefaultNew( Name name ) {
	return _AnalyticSolution_New(
			sizeof(Velic_solA),
			Velic_solA_Type,
			_AnalyticSolution_Delete,
			_AnalyticSolution_Print,
			_AnalyticSolution_Copy,
			_Velic_solA_DefaultNew,
			_Velic_solA_Construct,
			_Velic_solA_Build, 
			_AnalyticSolution_Initialise,
			_AnalyticSolution_Execute,
			_AnalyticSolution_Destroy,
			name );
}

Index _Velic_solA_Register( PluginsManager* pluginsManager ) {
	return PluginsManager_Submit( pluginsManager, Velic_solA_Type, "0", _Velic_solA_DefaultNew );
}
