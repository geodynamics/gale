#ifndef __Underworld_Velic_solB_h__
#define __Underworld_Velic_solB_h__

	extern const Type Velic_solB_Type;

	typedef struct {
		__AnalyticSolution
		double sigma;
		double Z;
		double km;
		int n;
	} Velic_solB;

	Index Underworld_Velic_solB_Register( PluginsManager* pluginsManager );
	void* _Velic_solB_DefaultNew( Name name );
	void _Velic_solB_Construct( void* analyticSolution, Stg_ComponentFactory* cf, void* data );
	void _Velic_solB_Init( Velic_solB* self, double sigma, double Z, double km, int n );

	void Velic_solB_PressureFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* pressure );
	void Velic_solB_VelocityFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* velocity );
	void Velic_solB_StressFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* stress );
	void Velic_solB_StrainRateFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* strainRate );

	void _Velic_solB( double* pos, 
		double sigma, double Z, int n, double km,
		double* velocity, double* pressure, double* Tstress, double* strainRate );

	Bool _checkInputParams( Velic_solB* self );
#endif
