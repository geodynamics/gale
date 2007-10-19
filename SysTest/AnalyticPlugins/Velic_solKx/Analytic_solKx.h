#ifndef __Underworld_Velic_solKx_h__
#define __Underworld_Velic_solKx_h__

	extern const Type Velic_solKx_Type;

	typedef struct {
		__AnalyticSolution
		double sigma;
		double _m;
	       	int n;
		double B;
	} Velic_solKx;

	Index Underworld_Velic_solKx_Register( PluginsManager* pluginsManager );
	void* _Velic_solKx_DefaultNew( Name name );
	void _Velic_solKx_Construct( void* analyticSolution, Stg_ComponentFactory* cf, void* data );
	void _Velic_solKx_Init( Velic_solKx* self, double sigma, double _m, double B, int n );

	void Velic_solKx_PressureFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* pressure );
	void Velic_solKx_VelocityFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* velocity );
	void Velic_solKx_StressFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* stress );
	void Velic_solKx_StrainRateFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* strainRate );

	void _Velic_solKx( 
		double pos[],
		double _sigma, /* density */
		double _m, int _n, /* wavelength in z, wavenumber in x */
		double _B, /* viscosity parameter */
		double vel[], double* presssure, 
		double total_stress[], double strain_rate[], double* viscosity );
#endif
