#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgFEM/StgFEM.h>

#ifndef __Underworld_Velic_solA_h__
#define __Underworld_Velic_solA_h__

	extern const Type Velic_solA_Type;

	typedef struct {
		__AnalyticSolution
		double sigma;
		double Z;
		double km;
		int n;
	} Velic_solA;

	void Velic_solA_PressureFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* pressure );
	void Velic_solA_VelocityFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* velocity );
	void Velic_solA_StressFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* stress );
	void Velic_solA_StrainRateFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* strainRate );

	void _Velic_solA( double* pos, 
		double sigma, double Z, int n, double km,
		double* velocity, double* pressure, double* Tstress, double* strainRate );

	Bool _checkInputParams( Velic_solA* self );
#endif
