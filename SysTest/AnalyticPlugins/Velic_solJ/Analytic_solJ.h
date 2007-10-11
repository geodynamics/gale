#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>

#ifndef __Underworld_Velic_solJ_h__
#define __Underworld_Velic_solJ_h__

	extern const Type Velic_solJ_Type;

	typedef struct {
		__AnalyticSolution
		double sigmaB;
	       	double sigmaA;
		double etaB;
		double etaA;
		double dxB;
		double dxA;
		double x0B;
		double x0A; 
		double zc;
	} Velic_solJ;

	Index Underworld_Velic_solJ_Register( PluginsManager* pluginsManager );
	void* _Velic_solJ_DefaultNew( Name name );
	void _Velic_solJ_Construct( void* analyticSolution, Stg_ComponentFactory* cf, void* data );
	void _Velic_solJ_Init( Velic_solJ* self, double sigmaA, double sigmaB, double etaB, double etaA, double dxB, double dxA, double x0B, double x0A, double zc );

	void Velic_solJ_PressureFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* pressure );
	void Velic_solJ_VelocityFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* velocity );
	void Velic_solJ_StressFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* stress );
	void Velic_solJ_StrainRateFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* strainRate );

	void _Velic_solJ( 
		double pos[],
		double _sigma_B ,double _sigma_A, /* density B, density A */
		double _eta_B   ,double _eta_A  , /* viscosity B, viscosity A */
		double _dx_B    ,double _dx_A   , /* width of the upper dense block, width of the lower dense block */
		double _x_0_B   ,double _x_0_A  , /* centre of the upper dense block, centre of lower dense block */
		double _z_c, /* bottom of the upper dense block */
		double vel[], double* presssure, 
		double total_stress[], double strain_rate[] );
#endif
