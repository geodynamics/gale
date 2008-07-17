#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>

#ifndef __Underworld_newVelicSolA_h__
#define __Underworld_newVelicSolA_h__

	extern const Type newVelicSolA_Type;

	typedef struct {
		__FieldTest
		double sigma;
		double Z;
		double km;
		int n;
	} newVelicSolA;

	void newVelicSolA_PressureFunction( void* analyticSolution, double* coord, double* pressure );
	void newVelicSolA_VelocityFunction( void* analyticSolution, double* coord, double* velocity );
	void newVelicSolA_StressFunction( void* analyticSolution, double* coord, double* stress );
	void newVelicSolA_StrainRateFunction( void* analyticSolution, double* coord, double* strainRate );

	void _Velic_solA( double* pos, 
		double sigma, double Z, int n, double km,
		double* velocity, double* pressure, double* Tstress, double* strainRate );

	Bool _checkInputParams( newVelicSolA* self );
#endif
