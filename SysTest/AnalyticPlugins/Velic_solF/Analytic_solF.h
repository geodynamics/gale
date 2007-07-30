#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgFEM/StgFEM.h>

#ifndef __Underworld_Velic_solF_h__
#define __Underworld_Velic_solF_h__

	extern const Type Velic_solF_Type;

	typedef struct {
		__AnalyticSolution
		double sigma; /* density */
		double etaA;
	       	double etaB; 
		double xc;
	       	double zc;
	} Velic_solF;

	Index ExperimentalUnderworld_Velic_solF_Register( PluginsManager* pluginsManager );
	void* _Velic_solF_DefaultNew( Name name );
	void _Velic_solF_Construct( void* analyticSolution, Stg_ComponentFactory* cf, void* data );
	void _Velic_solF_Init( Velic_solF* self, double sigma, double etaA, double etaB, double xc, double zc );

	void Velic_solF_PressureFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* pressure );
	void Velic_solF_VelocityFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* velocity );
	void Velic_solF_StressFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* stress );
	void Velic_solF_StrainRateFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* strainRate );

	void _Velic_solF(
		double pos[],
		double _sigma, /* density */
		double _eta_A, double _eta_B, /* viscosity A, viscosity B */ 
		double _x_c, double _z_c, /* width of dense block, bottom of dense block */
		double vel[], double* presssure, 
		double total_stress[], double strain_rate[] );

#endif
