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
** $Id: solA.c 565 2006-05-19 02:33:01Z RobertTurnbull $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgFEM/StgFEM.h>
#include <assert.h>

const Type Velic_solB_Type = "Velic_solB";

typedef struct {
	__AnalyticSolution
	double sigma;
	double Z;
	int waveNum_n;
	double waveNum_m;
} Velic_solB;


double Calculate__C1 ( Velic_solB* self, double x, double z, double y ) {
	double t1, t5, t6, t7, t12, t22, t26, t29, t33;
	double sigma = self->sigma;
	double Z = self->Z;
	double km = self->waveNum_m * M_PI;
	double kn = self->waveNum_n * M_PI;

	t1 = exp(km);
	t5 = kn * kn;
	t6 = exp(-kn);
	t7 = t6 * t6;
	t12 = km * km;
	t22 = pow(-0.1e1 + t6, 0.2e1);
	t26 = pow(t6 + 0.1e1, 0.2e1);
	t29 = pow(kn - km, 0.2e1);
	t33 = pow(kn + km, 0.2e1);
	return sigma * (-0.1e1 + t1) * (t1 + 0.1e1) * (t5 * t7 + t5 + 0.2e1 * kn - 0.2e1 * kn * t7 - t12 * t7 - t12) * t6 / Z / t1 / t22 / t26 / t29 / t33 / 0.4e1;
}

double Calculate__C2 ( Velic_solB* self, double x, double z, double y ) {
	double t1, t2, t3, t8, t11, t21, t25, t28, t32; 
	double sigma = self->sigma;
	double Z = self->Z;
	double km = self->waveNum_m * M_PI;
	double kn = self->waveNum_n * M_PI;

	t1 = kn * kn;
	t2 = exp(-kn);
	t3 = t2 * t2;
	t8 = km * km;
	t11 = exp(km);
	t21 = pow(-0.1e1 + t2, 0.2e1);
	t25 = pow(t2 + 0.1e1, 0.2e1);
	t28 = pow(kn - km, 0.2e1);
	t32 = pow(kn + km, 0.2e1);
	return -(t1 * t3 + t1 + 0.2e1 * kn - 0.2e1 * kn * t3 - t8 * t3 - t8) * (t11 + 0.1e1) * (-0.1e1 + t11) * sigma / Z / t11 / t21 / t25 / t28 / t32 / 0.4e1;

}

double Calculate__C3 ( Velic_solB* self, double x, double z, double y ) {
	double t1, t5;
	double sigma = self->sigma;
	double Z = self->Z;
	double km = self->waveNum_m * M_PI;
	double kn = self->waveNum_n * M_PI;

	t1 = exp(km);
	t5 = exp(-kn);
	return -sigma * (-0.1e1 + t1) * (t1 + 0.1e1) * t5 / Z / t1 / (-0.1e1 + t5) / (t5 + 0.1e1) / (kn - km) / (kn + km) / 0.4e1;
}

double Calculate__C4 ( Velic_solB* self, double x, double z, double y ) {
	double t1, t9;
	double sigma = self->sigma;
	double Z = self->Z;
	double km = self->waveNum_m * M_PI;
	double kn = self->waveNum_n * M_PI;

	t1 = exp(km);
	t9 = exp(-kn);
	return -(t1 + 0.1e1) * (-0.1e1 + t1) * sigma / Z / t1 / (-0.1e1 + t9) / (t9 + 0.1e1) / (kn - km) / (kn + km) / 0.4e1;
}

double Calculate_ss ( Velic_solB* self, double x, double z, double y ) {
	double t4, t10, t14, t15, t16, t18;
	double sigma = self->sigma;
	double Z = self->Z;
	double km = self->waveNum_m * M_PI;
	double kn = self->waveNum_n * M_PI;

	t4 = exp(-kn * z);
	t10 = exp(kn * (z - 0.1e1));
	t14 = sinh(km * z);
	t15 = km * km;
	t16 = kn * kn;
	t18 = pow(t15 - t16, 0.2e1);
	return (Calculate__C1( self, x, z, y) + z * Calculate__C3( self, x, z, y)) * t4 + (Calculate__C2( self, x, z, y) + z * Calculate__C4( self, x, z, y )) * t10 + kn * sigma * t14 / t18 / Z;
	
}

double Calculate_ss_z ( Velic_solB* self, double x, double z, double y ) {
	double t6, t14, t18, t20, t21, t23;
	double sigma = self->sigma;
	double Z = self->Z;
	double km = self->waveNum_m * M_PI;
	double kn = self->waveNum_n * M_PI;

	t6 = exp(-kn * z);
	t14 = exp(kn * (z - 0.1e1));
	t18 = cosh(km * z);
	t20 = km * km;
	t21 = kn * kn;
	t23 = pow(t20 - t21, 0.2e1);
	return (-(Calculate__C1( self, x, z, y ) + z * Calculate__C3( self, x, z, y )) * kn + Calculate__C3( self, x, z, y )) * t6 + (Calculate__C4( self, x, z, y ) + (Calculate__C2( self, x, z, y ) + z * Calculate__C4( self, x, z, y )) * kn) * t14 + kn * sigma * t18 * km / t23 / Z;
	
}

double Calculate_ss_zz ( Velic_solB* self, double x, double z, double y ) {
	double t3, t9, t19, t23, t25, t27;
	double sigma = self->sigma;
	double Z = self->Z;
	double km = self->waveNum_m * M_PI;
	double kn = self->waveNum_n * M_PI;

	t3 = kn * kn;
	t9 = exp(-kn * z);
	t19 = exp(kn * (z - 0.1e1));
	t23 = sinh(km * z);
	t25 = km * km;
	t27 = pow(t25 - t3, 0.2e1);
	return ((Calculate__C1( self, x, z, y ) + z * Calculate__C3( self, x, z, y )) * t3 - 0.2e1 * Calculate__C3( self, x, z, y ) * kn) * t9 + (0.2e1 * Calculate__C4( self, x, z, y ) * kn + (Calculate__C2( self, x, z, y ) + z * Calculate__C4( self, x, z, y )) * t3) * t19 + kn * sigma * t23 * t25 / t27 / Z;

}

double Calculate_ss_zzz ( Velic_solB* self, double x, double z, double y ) {
	double t3, t4, t10, t20, t24, t26, t29;
	double sigma = self->sigma;
	double Z = self->Z;
	double km = self->waveNum_m * M_PI;
	double kn = self->waveNum_n * M_PI;

	t3 = kn * kn;
	t4 = t3 * kn;
	t10 = exp(-kn * z);
	t20 = exp(kn * (z - 0.1e1));
	t24 = cosh(km * z);
	t26 = km * km;
	t29 = pow(t26 - t3, 0.2e1);
	return (-(Calculate__C1( self, x, z, y ) + z * Calculate__C3( self, x, z, y )) * t4 + 0.3e1 * Calculate__C3( self, x, z, y ) * t3) * t10 + (0.3e1 * Calculate__C4( self, x, z, y ) * t3 + (Calculate__C2( self, x, z, y ) + z * Calculate__C4( self, x, z, y )) * t4) * t20 + kn * sigma * t24 * t26 * km / t29 / Z;
}

void Velic_solB_PressureFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* pressure ) {
	Velic_solB* self = (Velic_solB*) analyticSolution;
	double x    = coord[ I_AXIS ];
	double z    = coord[ J_AXIS ];
	double y    = coord[ K_AXIS ];
	double n    = self->waveNum_n;
	double kn   = n * M_PI;
	double ss_zzz, ss_z, pp;

	ss_zzz = Calculate_ss_zzz( self, x, z, y );
	ss_z   = Calculate_ss_z( self, x, z, y );
	pp     = self->Z*(ss_zzz-kn*kn*ss_z)/kn;

	*pressure = pp * cos(kn * x);
}

void Velic_solB_VelocityFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* velocity ) {
	Velic_solB* self = (Velic_solB*) analyticSolution;
	double x    = coord[ I_AXIS ];
	double z    = coord[ J_AXIS ];
	double y    = coord[ K_AXIS ];
	double n    = self->waveNum_n;
	double u1, u2;


	u1 = n * M_PI * Calculate_ss( self, x, z, y );
	u2 = -Calculate_ss_z( self, x, z, y );

	u1 *= cos( n * M_PI * x );
	u2 *= sin( n * M_PI * x );

	velocity[ I_AXIS ] = u2;
	velocity[ J_AXIS ] = u1;
	if ( analyticFeVariable->dim == 3 ) 
		velocity[ K_AXIS ] = 0.0;
}

void Velic_solB_StressFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* stress ) {
	Velic_solB* self = (Velic_solB*) analyticSolution;
	double x    = coord[ I_AXIS ];
	double z    = coord[ J_AXIS ];
	double y    = coord[ K_AXIS ];
	double n    = self->waveNum_n;
	double Z    = self->Z;
	double kn   = n * M_PI;
	double pp, u3, u4, txx;
	double ss_zzz, ss_z, ss;

	ss_zzz      = Calculate_ss_zzz( self, x, z, y );
	ss_z        = Calculate_ss_z( self, x, z, y );
	ss          = Calculate_ss( self, x, z, y );

	pp = Z*(ss_zzz-kn*kn*ss_z)/kn;
	u3 = 2.0*Z*kn*ss_z - pp;
	u4 = -Z*(ss_zzz + kn*kn*ss);
	txx = -2.0*Z*kn*ss_z - pp;
	 
	stress[1] = u3 * cos(kn*x); /* zz stress */
	stress[2] = u4 * sin(kn*x); /* zx stress */
	stress[0] = txx * cos(kn*x); /* xx stress */
	
}
void Velic_solB_StrainRateFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* strainRate ) {
	Velic_solB* self = (Velic_solB*) analyticSolution;
	double x    = coord[ I_AXIS ];
	double z    = coord[ J_AXIS ];
	double y    = coord[ K_AXIS ];
	double n    = self->waveNum_n;
	double kn   = n * M_PI;
	double ss_zz, ss_z, ss;

	ss_zz       = Calculate_ss_zz( self, x, z, y );
	ss_z        = Calculate_ss_z( self, x, z, y );
	ss          = Calculate_ss( self, x, z, y );

	if( analyticFeVariable->dim == 2 ) {
		strainRate[1] = kn*ss_z*cos(kn*x);                // zz 
		strainRate[0] = -strainRate[1];                   // xx 
		strainRate[2] = -0.5*(ss_zz+kn*kn*ss)*sin(kn*x); // xz 
	}

}
	
void _Velic_solB_Construct( void* analyticSolution, Stg_ComponentFactory* cf, void* data ) {
	Velic_solB* self = (Velic_solB*) analyticSolution;
	FeVariable*              velocityField;
	FeVariable*              pressureField;
	FeVariable*              stressField;
	FeVariable*              strainRateField;
	FeVariable*              recoverdStrainRateField;

	/* Construct Parent */
	_AnalyticSolution_Construct( self, cf, data );

	/* Create Analytic Fields */
	velocityField = Stg_ComponentFactory_ConstructByName( cf, "VelocityField", FeVariable, True, data ); 
	AnalyticSolution_CreateAnalyticVectorField( self, velocityField, Velic_solB_VelocityFunction );

	pressureField = Stg_ComponentFactory_ConstructByName( cf, "PressureField", FeVariable, True, data ); 
	AnalyticSolution_CreateAnalyticField( self, pressureField, Velic_solB_PressureFunction );

	stressField = Stg_ComponentFactory_ConstructByName( cf, "StressField", FeVariable, False, data ); 
	if ( stressField )
		AnalyticSolution_CreateAnalyticSymmetricTensorField( self, stressField, Velic_solB_StressFunction );

	strainRateField = Stg_ComponentFactory_ConstructByName( cf, "StrainRateField", FeVariable, False, data ); 
	if ( strainRateField  ) {
		AnalyticSolution_CreateAnalyticSymmetricTensorField( self, strainRateField, Velic_solB_StrainRateFunction );
	}

	recoverdStrainRateField = Stg_ComponentFactory_ConstructByName( cf, "recoveredStrainRate", FeVariable, False, data ); 
	if( recoverdStrainRateField ) {
		AnalyticSolution_CreateAnalyticSymmetricTensorField( self, recoverdStrainRateField, Velic_solB_StrainRateFunction );
	}

	self->sigma = Stg_ComponentFactory_GetRootDictDouble( cf, "sigma", 1.0 );
	self->Z = Stg_ComponentFactory_GetRootDictDouble( cf, "Z", 1.0 );
	self->waveNum_n = Stg_ComponentFactory_GetRootDictInt( cf, "wavenumberX", 1 );
	self->waveNum_m = Stg_ComponentFactory_GetRootDictDouble( cf, "wavenumberY", 2 );
	assert( self->waveNum_n != self->waveNum_m );
}

void* _Velic_solB_DefaultNew( Name name ) {
	return _AnalyticSolution_New(
			sizeof(Velic_solB),
			Velic_solB_Type,
			_AnalyticSolution_Delete,
			_AnalyticSolution_Print,
			_AnalyticSolution_Copy,
			_Velic_solB_DefaultNew,
			_Velic_solB_Construct,
			_AnalyticSolution_Build,
			_AnalyticSolution_Initialise,
			_AnalyticSolution_Execute,
			_AnalyticSolution_Destroy,
			name );
}

Index _Velic_solB_Register( PluginsManager* pluginsManager ) {
	return PluginsManager_Submit( pluginsManager, Velic_solB_Type, "0", _Velic_solB_DefaultNew );
}
