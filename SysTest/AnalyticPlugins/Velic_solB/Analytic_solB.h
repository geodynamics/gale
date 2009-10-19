#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>

#ifndef __Underworld_solB_h__
#define __Underworld_solB_h__

	extern const Type Underworld_solB_Type;

	typedef struct {
		__FieldTest
		double sigma;
		double Z;
		double km;
		int n;
	} Underworld_solB;

	void Underworld_solB_PressureFunction( void* analyticSolution, double* coord, double* pressure );
	void Underworld_solB_VelocityFunction( void* analyticSolution, double* coord, double* velocity );
	void Underworld_solB_StressFunction( void* analyticSolution, double* coord, double* stress );
	void Underworld_solB_StrainRateFunction( void* analyticSolution, double* coord, double* strainRate );

	void _Velic_solB( double* pos, 
		double sigma, double Z, int n, double km,
		double* velocity, double* pressure, double* Tstress, double* strainRate );

	Bool solB_checkInputParams( Underworld_solB* self );
#endif
