#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>

#ifndef __Underworld_Velic_solHy_h__
#define __Underworld_Velic_solHy_h__

	extern const Type Velic_solHy_Type;

	typedef struct {
		__AnalyticSolution
	//	double sigma;
		double etaA;  /* Input parameters: density, viscosity A */
		double etaB;  /* Input parameters: density, viscosity B */
		double xc;    /* Input parameters: viscosity jump location */
		int n;
	} Velic_solHy;

	Index Underworld_Velic_solHy_Register( PluginsManager* pluginsManager );
	void* _Velic_solHy_DefaultNew( Name name );
	void _Velic_solHy_Construct( void* analyticSolution, Stg_ComponentFactory* cf, void* data );
	void _Velic_solHy_Init( Velic_solHy* self, double etaA, double etaB, double xc, int n );

	void Velic_solHy_PressureFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* pressure );
	void Velic_solHy_VelocityFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* velocity );
	void Velic_solHy_StressFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* stress );
	void Velic_solHy_StrainRateFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* strainRate );

#endif
