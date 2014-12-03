#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>

#ifndef __Underworld_Velic_solCA_h__
#define __Underworld_Velic_solCA_h__

	extern const Type Velic_solCA_Type;

	typedef struct {
		__AnalyticSolution
		double sigma;
		double eta;  /* Input parameters: density, viscosity A */
		double dx;  
		double x0;  
	} Velic_solCA;

	Index Underworld_Velic_solCA_Register( PluginsManager* pluginsManager );
	void* _Velic_solCA_DefaultNew( Name name );
	void _Velic_solCA_Init( Velic_solCA* self, double sigma, double eta, double dx, double x0 ); 
	void _Velic_solCA_AssignFromXML( void* analyticSolution, Stg_ComponentFactory* cf, void* data );

	void Velic_solCA_PressureFunction( void* analyticSolution, FeVariable* analyticFeVariable, const double *coord, double* pressure );
	void Velic_solCA_VelocityFunction( void* analyticSolution, FeVariable* analyticFeVariable, const double *coord, double* velocity );
	void Velic_solCA_StressFunction( void* analyticSolution, FeVariable* analyticFeVariable, const double *coord, double* stress );
	void Velic_solCA_StrainRateFunction( void* analyticSolution, FeVariable* analyticFeVariable, const double *coord, double* strainRate );

	void _Velic_solCA( 
		const double pos[], 
		double _sigma, double _eta, double _dx, double _x_0,
		double vel[], double* presssure, 
		double total_stress[], double strain_rate[] );

	Bool _checkInputParams( Velic_solCA* self );	

#endif
