#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgFEM/StgFEM.h>

#ifndef __Underworld_Velic_solG_h__
#define __Underworld_Velic_solG_h__

	extern const Type Velic_solG_Type;

	typedef struct {
		__AnalyticSolution
		double sigma;
		double etaA;
	       	double etaB; 
		double dx;
	       	double x0;
	       	double zc;
	} Velic_solG;

	Index Underworld_Velic_solG_Register( PluginsManager* pluginsManager );
	void* _Velic_solG_DefaultNew( Name name );
	void _Velic_solG_Construct( void* analyticSolution, Stg_ComponentFactory* cf, void* data );
	void _Velic_solG_Init( Velic_solG* self, double sigma, double etaA, double etaB, double dx, double x0, double zc );

	void Velic_solG_PressureFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* pressure );
	void Velic_solG_VelocityFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* velocity );
	void Velic_solG_StressFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* stress );
	void Velic_solG_StrainRateFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* strainRate );

	void _Velic_solG(
		double pos[],
		double _sigma,
		double _eta_A, double _eta_B, 
		double _dx, double _x_0, double _z_c,
		double vel[], double* presssure, 
		double total_stress[], double strain_rate[] );

#endif
