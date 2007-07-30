#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgFEM/StgFEM.h>

#ifndef __Underworld_Velic_solKz_h__
#define __Underworld_Velic_solKz_h__

	extern const Type Velic_solKz_Type;

	typedef struct {
		__AnalyticSolution
		double sigma;
		double km;
	       	int n;
		double B;
	} Velic_solKz;

	Index Underworld_Velic_solKz_Register( PluginsManager* pluginsManager );
	void* _Velic_solKz_DefaultNew( Name name );
	void _Velic_solKz_Construct( void* analyticSolution, Stg_ComponentFactory* cf, void* data );
	void _Velic_solKz_Init( Velic_solKz* self, double sigma, double km, double B, int n );
	void Velic_solKz_PressureFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* pressure );
	void Velic_solKz_VelocityFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* velocity );
	void Velic_solKz_StressFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* stress );
	void Velic_solKz_StrainRateFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* strainRate );

	void _Velic_solKz( 
		double pos[],
		double _sigma, /* density */
		double _km, int _n, /* wavelength in z, wavenumber in x */
		double _B, /* viscosity parameter */
		double vel[], double* presssure, 
		double total_stress[], double strain_rate[] );
#endif
