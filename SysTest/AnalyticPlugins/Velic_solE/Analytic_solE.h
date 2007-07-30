#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgFEM/StgFEM.h>

#ifndef __Underworld_Velic_solE_h__
#define __Underworld_Velic_solE_h__

	extern const Type Velic_solE_Type;

	typedef struct {
		__AnalyticSolution
		double sigma;
		double etaA;
		double etaB;
		double zc;
		double km;
		int n;
	} Velic_solE;

	Index ExperimentalUnderworld_Velic_solE_Register( PluginsManager* pluginsManager );
	void* _Velic_solE_DefaultNew( Name name );
	void _Velic_solE_Construct( void* analyticSolution, Stg_ComponentFactory* cf, void* data );
	void _Velic_solE_Init( Velic_solE* self, double sigma, double etaA, double etaB, double zc, double km, double n );
	void Velic_solE_PressureFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* pressure );
	void Velic_solE_VelocityFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* velocity );
	void Velic_solE_StressFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* stress );
	void Velic_solE_StrainRateFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* strainRate );

	void _Velic_solE(
		double pos[],
		double _sigma,
		double _eta_A, double _eta_B, 
		double _z_c, double _km, int _n,
		double vel[], double* presssure, 
		double total_stress[], double strain_rate[] );
#endif
