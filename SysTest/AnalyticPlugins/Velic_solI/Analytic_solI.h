#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>

#ifndef __Underworld_Velic_solI_h__
#define __Underworld_Velic_solI_h__

	extern const Type Velic_solI_Type;

	typedef struct {
		__AnalyticSolution
		double sigma; /* density */
		double B;
	       	double xc; /* width of dense column */
	
	} Velic_solI;

	Index Underworld_Velic_solI_Register( PluginsManager* pluginsManager );
	void* _Velic_solI_DefaultNew( Name name );
	void _Velic_solI_AssignFromXML( void* analyticSolution, Stg_ComponentFactory* cf, void* data );
	void _Velic_solI_Init( Velic_solI* self, double sigma, double B, double xc );

	void Velic_solI_PressureFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* pressure );
	void Velic_solI_VelocityFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* velocity );
	void Velic_solI_StressFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* stress );
	void Velic_solI_StrainRateFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* strainRate );

	void _Velic_solI( 
		double pos[],
		double _sigma, /* density */
		double _B, double _x_c, /* viscosity parameter, width of dense column */
		double vel[], double* presssure, 
		double total_stress[], double strain_rate[] );
#endif
