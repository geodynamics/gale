#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>

#ifndef __Underworld_solA_h__
#define __Underworld_solA_h__

	extern const Type Underworld_solA_Type;

	typedef struct {
		__FieldTest
		double sigma;
		double Z;
		double km;
		int n;
	} Underworld_solA;


	void Underworld_solA_PressureFunction( void* analyticSolution, const double* coord, double* pressure );
	void Underworld_solA_VelocityFunction( void* analyticSolution, const double* coord, double* velocity );
	void Underworld_solA_StressFunction( void* analyticSolution, const double* coord, double* stress );
	void Underworld_solA_StrainRateFunction( void* analyticSolution, const double* coord, double* strainRate );

	void _Velic_solA( double* pos, 
		double sigma, double Z, int n, double km,
		double* velocity, double* pressure, double* Tstress, double* strainRate );

	Bool solA_checkInputParams( Underworld_solA* self );
#endif
