#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>

#ifndef __Underworld_Velic_solIA_h__
#define __Underworld_Velic_solIA_h__

	extern const Type Velic_solIA_Type;

	typedef struct {
		__AnalyticSolution
		double sigma;
		double B;
		double dx;
		double x0;
	} Velic_solIA;

	Index Underworld_Velic_solIA_Register( PluginsManager* pluginsManager );
	void* _Velic_solIA_DefaultNew( Name name );
	void _Velic_solIA_AssignFromXML( void* analyticSolution, Stg_ComponentFactory* cf, void* data );
	void _Velic_solIA_Init( Velic_solIA* self, double sigma, double B, double dx, double x0 );

	void Velic_solIA_PressureFunction( void* analyticSolution, FeVariable* analyticFeVariable, const double *coord, double* pressure );
	void Velic_solIA_VelocityFunction( void* analyticSolution, FeVariable* analyticFeVariable, const double *coord, double* velocity );
	void Velic_solIA_StressFunction( void* analyticSolution, FeVariable* analyticFeVariable, const double *coord, double* stress );
	void Velic_solIA_StrainRateFunction( void* analyticSolution, FeVariable* analyticFeVariable, const double *coord, double* strainRate );

	void _Velic_solIA( 
		const double pos[],
		double _sigma, double BB, /* density, viscosity parameter */
		double _dx, double _x_0, /* width of dense column, centre of dense column */
		double vel[], double* presssure, 
		double total_stress[], double strain_rate[] );
#endif
