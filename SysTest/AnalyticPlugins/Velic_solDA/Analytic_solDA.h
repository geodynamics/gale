#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>

#ifndef __Underworld_Velic_solDA_h__
#define __Underworld_Velic_solDA_h__

	extern const Type Velic_solDA_Type;

	typedef struct {
		__AnalyticSolution
		double sigma;
		double etaA;  
		double etaB; 
		double zc;    
		double dx;    
		double x0;    
	} Velic_solDA;

	Index Underworld_Velic_solDA_Register( PluginsManager* pluginsManager );
	void* _Velic_solDA_DefaultNew( Name name );
	void _Velic_solDA_AssignFromXML( void* analyticSolution, Stg_ComponentFactory* cf, void* data );
	void _Velic_solDA_Init( Velic_solDA* self, double sigma, double etaA, double etaB, double zc, double dx, double x0 );

	void Velic_solDA_PressureFunction( void* analyticSolution, FeVariable* analyticFeVariable, const double *coord, double* pressure );
	void Velic_solDA_VelocityFunction( void* analyticSolution, FeVariable* analyticFeVariable, const double *coord, double* velocity );
	void Velic_solDA_StressFunction( void* analyticSolution, FeVariable* analyticFeVariable, const double *coord, double* stress );
	void Velic_solDA_StrainRateFunction( void* analyticSolution, FeVariable* analyticFeVariable, const double *coord, double* strainRate );

	void _Velic_solDA(		
		const double pos[],
		double _sigma, /* density */
		double _eta_A, double _eta_B, /* viscosity A, viscosity B */ 
		double _z_c, double _dx, double _x_0, /* viscosity jump location, width of dense column, centre of dense column */
		double vel[], double* presssure, 
		double total_stress[], double strain_rate[] );	
#endif
