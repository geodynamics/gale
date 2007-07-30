#ifndef __Underworld_Velic_solS_h__
#define __Underworld_Velic_solS_h__

	extern const Type Velic_solS_Type;

	typedef struct {
		__AnalyticSolution
		int _n;
		double _eta;
	} Velic_solS;

	Index Underworld_Velic_solS_Register( PluginsManager* pluginsManager );
	void* _Velic_solS_DefaultNew( Name name );
	void _Velic_solS_Construct( void* analyticSolution, Stg_ComponentFactory* cf, void* data );
	void _Velic_solS_Init( Velic_solS* self, double eta, int _n );

	void Velic_solS_PressureFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* pressure );
	void Velic_solS_VelocityFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* velocity );
	void Velic_solS_StressFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* stress );
	void Velic_solS_StrainRateFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* strainRate );

	void _Velic_solS( 
		double pos[],
		int _n, double _eta,
		double vel[], double* presssure, 
		double total_stress[], double strain_rate[] );

	Bool _checkInputParams( Velic_solS* self );
#endif
