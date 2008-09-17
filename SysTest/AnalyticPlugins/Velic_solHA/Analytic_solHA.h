#ifndef __Underworld_Velic_solHA_h__
#define __Underworld_Velic_solHA_h__

	extern const Type Velic_solHA_Type;

	typedef struct {
		__AnalyticSolution
		double sigma;
		double eta;
		double dx;
		double dy;
		double x0;
		double y0;
	} Velic_solHA;

	Index Underworld_Velic_solHA_Register( PluginsManager* pluginsManager );
	void* _Velic_solHA_DefaultNew( Name name );
	void _Velic_solHA_Construct( void* analyticSolution, Stg_ComponentFactory* cf, void* data );
	void _Velic_solHA_Init( Velic_solHA* self, double sigma, double eta, double dx, double dy, double x0, double y0 );

	void Velic_solHA_PressureFunction( void* analyticSolution, double* coord, double* pressure );
	void Velic_solHA_VelocityFunction( void* analyticSolution, double* coord, double* velocity );
	void Velic_solHA_StressFunction( void* analyticSolution, double* coord, double* stress );
	void Velic_solHA_StrainRateFunction( void* analyticSolution, double* coord, double* strainRate );

	void _Velic_solHA( 
		double pos[],
		double _sigma, double _eta,
		double _dx, double _dy,
		double _x_0, double _y_0,
		double vel[], double* presssure, 
		double total_stress[], double strain_rate[] );

#endif
