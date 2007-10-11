#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>

#ifndef __Underworld_Velic_solH_h__
#define __Underworld_Velic_solH_h__

	extern const Type Velic_solH_Type;

	typedef struct {
		__AnalyticSolution
		double sigma;
		double eta;
		double dx;
	       	double dy;
	} Velic_solH;

	Index Underworld_Velic_solH_Register( PluginsManager* pluginsManager );
	void* _Velic_solH_DefaultNew( Name name );
	void _Velic_solH_Construct( void* analyticSolution, Stg_ComponentFactory* cf, void* data );
	void _Velic_solH_Init( Velic_solH* self, double sigma, double eta, double dx, double dy );
	void Velic_solH_PressureFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* pressure );
	void Velic_solH_VelocityFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* velocity );
	void Velic_solH_StressFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* stress );
	void Velic_solH_StrainRateFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* strainRate );

	void _Velic_solH(
		double pos[],
		double _sigma,
		double _eta,
		double _dx, double _dy,
		double vel[], double* presssure, 
		double total_stress[], double strain_rate[] );
#endif
