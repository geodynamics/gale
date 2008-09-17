#ifndef __Underworld_Velic_solCx_h__
#define __Underworld_Velic_solCx_h__

	extern const Type Underworld_solCx_Type;

	typedef struct {
		__FieldTest
		double etaA;  /* Input parameters: density, viscosity A */
		double etaB;  /* Input parameters: density, viscosity B */
		double xc;    /* Input parameters: viscosity jump location */
		int n;
	} Velic_solCx;

	Index Underworld_Velic_solCx_Register( PluginsManager* pluginsManager );
	void* _Underworld_solCx_DefaultNew( Name name );
	void _Underworld_solCx_Construct( void* analyticSolution, Stg_ComponentFactory* cf, void* data );
	void _Underworld_solCx_Init( Velic_solCx* self, double etaA, double etaB, double xc, int n );

	void Velic_solCx_PressureFunction( void* analyticSolution, double* coord, double* pressure );
	void Velic_solCx_VelocityFunction( void* analyticSolution, double* coord, double* velocity );
	void Velic_solCx_StressFunction( void* analyticSolution, double* coord, double* stress );
	void Velic_solCx_StrainRateFunction( void* analyticSolution, double* coord, double* strainRate );

	void _Velic_solCx(
		double pos[], 
		double _eta_A, double _eta_B, 	/* Input parameters: density, viscosity A, viscosity B */
		double _x_c, int _n, 			/* Input parameters: viscosity jump location, wavenumber in x */
		double vel[], double* presssure, 
		double total_stress[], double strain_rate[] );

#endif
