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
#include "solCA.h"

const Type Velic_solCA_Type = "Velic_solCA";

#define PI 3.14159265358979323846264338328
void Velic_solCA_MirkoFunction( void* analyticSolution, double* coord, int whatWeWant, double* result ) {
	Velic_solCA* self = (Velic_solCA*)analyticSolution;
	
	double Z,u1,u2,u3,u4,u5,u6;
   double _C1,_C2,_C3,_C4;
   double sum1,sum2,sum3,sum4,sum5,sum6,sum7,mag,x,z;
   double sigma,del_rho,k,x0,dx;
   int n;
   double t1,t2,t3,t4,t5,t6,t7,t8,t10,t11;
   double t12,t14,t16,t21;

   /*  XML input and coord */
   /*************************************************************************/
   /*************************************************************************/
   sigma = self->sigma;
   Z = self->Z;
   x0 = self->x0;
   dx = self->dx;
   x = coord[ I_AXIS ];
   z = coord[ J_AXIS ];
   /*************************************************************************/
   /*************************************************************************/
   sum1=0.0;
	 sum2=0.0;
	 sum3=0.0;
	 sum4=0.0;
	 sum5=0.0;
	 sum6=0.0;
	 sum7=0.0;
	 


	 for(n=1;n<155;n++){

	    k = (double)n*M_PI;
	    del_rho = 4.0*sigma*cos(k*x0)*sin(k*dx/2.0)/k;
	    
	    t1 = exp(-k);
	    t7 = k * k;
	    t12 = pow(t1 + 0.1e1, 0.2e1);
	    _C1 = -del_rho * (-0.2e1 + k * t1 - 0.2e1 * t1) / Z / t7 / k / t12 / 0.2e1;

	    t1 = exp(-k);
	    t6 = k * k;
	    t11 = pow(t1 + 0.1e1, 0.2e1);
	    _C2 = del_rho * (k + 0.2e1 + 0.2e1 * t1) / Z / t6 / k / t11 / 0.2e1;

	    t3 = k * k;
	    t5 = exp(-k);
	    _C3 = del_rho / Z / t3 / (t5 + 0.1e1) / 0.2e1;

	    t3 = k * k;
	    t5 = exp(-k);
	    _C4 = -del_rho / Z / t3 / (t5 + 0.1e1) / 0.2e1;

	    /* Vz */
	    t4 = exp(-k * z);
	    t10 = exp(k * (z - 0.1e1));
	    t14 = k * k;
	    u1 = k * ((_C1 + z * _C3) * t4 + (_C2 + z * _C4) * t10 - del_rho / Z / t14 / k);
	    /* Vx */
	    t7 = exp(k * (z - 0.1e1));
	    t14 = exp(-k * z);
	    u2 = (-_C4 - (_C2 + z * _C4) * k) * t7 + (-_C3 + (_C1 + z * _C3) * k) * t14;
	    /* tzz */
	    t1 = Z * k;
	    t8 = exp(k * (z - 0.1e1));
	    t16 = exp(-k * z);
	    u3 = 0.2e1 * t1 * (_C4 + (_C2 + z * _C4) * k) * t8 + 0.2e1 * t1 * (_C3 - (_C1 + z * _C3) * k) * t16;
	    /* txz */
	    t2 = k * k;
	    t11 = exp(k * (z - 0.1e1));
	    t21 = exp(-k * z);
	    u4 = -Z * (0.2e1 * _C4 * k + 0.2e1 * t2 * (_C2 + z * _C4)) * t11 - Z * (-0.2e1 * _C3 * k + 0.2e1 * t2 * (_C1 + z * _C3)) * t21 + 0.1e1 / k * del_rho;
	    /* txx */
	    t1 = Z * k;
	    t8 = exp(k * (z - 0.1e1));
	    t16 = exp(-k * z);
	    u6 = -0.2e1 * t1 * (_C4 + (_C2 + z * _C4) * k) * t8 - 0.2e1 * t1 * (_C3 - (_C1 + z * _C3) * k) * t16;
	    /* pressure */
	    t1 = Z * k;
	    t4 = exp(k * (z - 0.1e1));
	    t8 = exp(-k * z);
	    u5 = 0.2e1 * t1 * _C4 * t4 + 0.2e1 * t1 * _C3 * t8;


	    

	    u5 = u5*cos(n*PI*x); /* pressure */
	    u6 = u6*cos(n*PI*x); /* xx stress */	    
	    sum5 +=u5;
	    sum6 +=u6;

	    u1 *= cos(n*PI*x); /* z velocity */
	    sum1 += u1;
	    u2 *= sin(n*PI*x); /* x velocity */
	    sum2 += u2;
	    u3 *= 2*n*PI*cos(n*PI*x); /* zz stress */
	    sum3 += u3;
	    u4 *= 2*n*PI*sin(n*PI*x); /* zx stress */
	    sum4 += u4;
	    /* density */
	    sum7 += del_rho*cos(k*x);

	 }/* n */
	 /* n=0 term*/
	 sum7 += sigma*dx;
	 /* pressure 0th term integration constant is arbitrarily chosen so that this term is 0 at z=0.5 */
	 sum5 += sigma*dx*(0.5-z); /* now have total pressure */
	 sum3 += -sigma*dx*(0.5-z); /* now have total zz stress */
	 sum6 += -sigma*dx*(0.5-z); /* now have total xx stress */
	 mag=sqrt(sum1*sum1+sum2*sum2);
	 /* printf("%0.7f %0.7f %0.7f %0.7f %0.7f %0.7f %0.7f %0.7f %0.7f %0.7f\n",x,z,sum1,sum2,sum3,sum4,sum5,sum6,mag,sum7); */
	 
	 if ( whatWeWant == 1 ) { 
		 result[0] = sum2; /* vx */
		 result[1] = sum1; /* vz */
	 }else if( whatWeWant == 2 ) {
		 result[0] = sum6; /* txx */
		 result[1] = sum3; /* tzz */
		 result[2] = sum4; /* txz */
	 } else if (whatWeWant == 3 ) {
		 *result = sum5;
	 } else if( whatWeWant == 4 ) {
		 result[0] = ( sum6 - 0.5*(sum6 + sum3) ) / ( 2 * Z );
		 result[1] = ( sum3 - 0.5*(sum6 + sum3) ) / ( 2 * Z );
		 result[2] = ( sum4 ) / ( 2 * Z ); 
	 }
}

void Velic_solCA_PressureFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* pressure ) {
	Velic_solCA_MirkoFunction( analyticSolution, coord, 3, pressure );
}

void Velic_solCA_VelocityFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* velocity ){
	Velic_solCA_MirkoFunction( analyticSolution, coord, 1, velocity );
}

void Velic_solCA_StressFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* stress ) {
	Velic_solCA_MirkoFunction( analyticSolution, coord, 2, stress );
}

void Velic_solCA_StrainRateFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* strainRate ) {
	Velic_solCA_MirkoFunction( analyticSolution, coord, 4, strainRate );
}
	
void _Velic_solCA_Construct( void* analyticSolution, Stg_ComponentFactory* cf, void* data ) {
	Velic_solCA* self = (Velic_solCA*) analyticSolution;
	FeVariable*              velocityField;
	FeVariable*              pressureField;
	FeVariable*              stressField;
	FeVariable*              strainRateField;
	/* Construct Parent */
	_AnalyticSolution_Construct( self, cf, data );

	/* Create Analytic Fields */
	velocityField = Stg_ComponentFactory_ConstructByName( cf, "VelocityField", FeVariable, True, data );
	AnalyticSolution_CreateAnalyticVectorField( self, velocityField, Velic_solCA_VelocityFunction );

	pressureField = Stg_ComponentFactory_ConstructByName( cf, "PressureField", FeVariable, True, data );
	AnalyticSolution_CreateAnalyticField( self, pressureField, Velic_solCA_PressureFunction );

	strainRateField = Stg_ComponentFactory_ConstructByName( cf, "StrainRateField", FeVariable, False, data );
	AnalyticSolution_CreateAnalyticField( self, strainRateField, Velic_solCA_StrainRateFunction );

	stressField = Stg_ComponentFactory_ConstructByName( cf, "StressField", FeVariable, False, data );
	if ( stressField )
		AnalyticSolution_CreateAnalyticSymmetricTensorField( self, stressField, Velic_solCA_StressFunction );

	
	self->sigma = Stg_ComponentFactory_GetRootDictDouble( cf, "solCA_sigma", 1.0 );
	self->Z = Stg_ComponentFactory_GetRootDictDouble( cf, "solCA_Z", 1.0);
	self->x0 = Stg_ComponentFactory_GetRootDictDouble( cf, "solCA_x0", 0.4 );
	self->dx = Stg_ComponentFactory_GetRootDictDouble( cf, "solCA_dx", 0.2 );
}

void* _Velic_solCA_DefaultNew( Name name ) {
	return _AnalyticSolution_New(
			sizeof(Velic_solCA),
			Velic_solCA_Type,
			_AnalyticSolution_Delete,
			_AnalyticSolution_Print,
			_AnalyticSolution_Copy,
			_Velic_solCA_DefaultNew,
			_Velic_solCA_Construct,
			_AnalyticSolution_Build,
			_AnalyticSolution_Initialise,
			_AnalyticSolution_Execute,
			_AnalyticSolution_Destroy,
			name );
}

Index _Velic_solCA_Register( PluginsManager* pluginsManager ) {
	return PluginsManager_Submit( pluginsManager, Velic_solCA_Type, "0", _Velic_solCA_DefaultNew );
}
