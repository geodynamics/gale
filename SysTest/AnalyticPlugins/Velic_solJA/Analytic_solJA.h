#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>

#ifndef __Underworld_Velic_solJA_h__
#define __Underworld_Velic_solJA_h__

	extern const Type Velic_solJA_Type;

	typedef struct {
		__AnalyticSolution
		double sigma;
		double etaA;
	       	double etaB;
		double dx;
		double x0;
	       	double zc;
	} Velic_solJA;

	Index Underworld_Velic_solJA_Register( PluginsManager* pluginsManager );
	void* _Velic_solJA_DefaultNew( Name name );
	void _Velic_solJA_Construct( void* analyticSolution, Stg_ComponentFactory* cf, void* data );
	void _Velic_solJA_Init( Velic_solJA* self, double sigma, double etaA, double etaB, double dx, double x0, double zc );
	void Velic_solJA_PressureFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* pressure );
	void Velic_solJA_VelocityFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* velocity );
	void Velic_solJA_StressFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* stress );
	void Velic_solJA_StrainRateFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* strainRate );

	void _Velic_solJA( 
		double pos[],
		double _sigma, /* density */
		double _eta_A, double _eta_B, /* viscosity A, viscosity B */
		double _dx, /* width of upper dense block */
		double _x_0, double _z_c, /* centre of upper dense block, bottom of upper dense block */
		double vel[], double* presssure, 
		double total_stress[], double strain_rate[] );

#endif
