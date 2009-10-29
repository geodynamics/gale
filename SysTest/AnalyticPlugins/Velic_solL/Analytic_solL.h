#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>

#ifndef __Underworld_Velic_solL_h__
#define __Underworld_Velic_solL_h__

	extern const Type Velic_solL_Type;

	typedef struct {
		__AnalyticSolution
		double sigmaB;
	       	double sigmaA;
		double eta; 
	} Velic_solL;

	Index Underworld_Velic_solL_Register( PluginsManager* pluginsManager );
	void* _Velic_solL_DefaultNew( Name name );
	void _Velic_solL_AssignFromXML( void* analyticSolution, Stg_ComponentFactory* cf, void* data );
	void _Velic_solL_Init( Velic_solL* self, double sigmaB, double sigmaA, double eta );

	void Velic_solL_PressureFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* pressure );
	void Velic_solL_VelocityFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* velocity );
	void Velic_solL_StressFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* stress );
	void Velic_solL_StrainRateFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* strainRate );

	void _Velic_solL( 
		double pos[],
		double _sigma_B, double _sigma_A, /* density B, density A */
		double _eta, /* viscosity */
		double vel[], double* presssure, 
		double total_stress[], double strain_rate[] );
#endif
