#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>

#ifndef __Underworld_Velic_solC_h__
#define __Underworld_Velic_solC_h__

	extern const Type Velic_solC_Type;

	typedef struct {
		__AnalyticSolution
		double sigma;
		double eta;  /* Input parameters: density, viscosity A */
		double x_c;  
	} Velic_solC;

	void* _Velic_solC_DefaultNew( Name name );
	void _Velic_solC_Init( Velic_solC* self, double sigma, double eta, double x_c );
	void _Velic_solC_AssignFromXML( void* analyticSolution, Stg_ComponentFactory* cf, void* data );

	void Velic_solC_PressureFunction( void* analyticSolution, FeVariable* analyticFeVariable, const double *coord, double* pressure );
	void Velic_solC_VelocityFunction( void* analyticSolution, FeVariable* analyticFeVariable, const double *coord, double* velocity );
	void Velic_solC_StressFunction( void* analyticSolution, FeVariable* analyticFeVariable, const double *coord, double* stress );
	void Velic_solC_StrainRateFunction( void* analyticSolution, FeVariable* analyticFeVariable, const double *coord, double* strainRate );

	void _Velic_solC( 
		const double pos[], 
		double sigma, double eta, double x_c,
		double vel[], double* presssure, 
		double total_stress[], double strain_rate[] );

	Bool _checkInputParams( Velic_solC* self );

#endif
