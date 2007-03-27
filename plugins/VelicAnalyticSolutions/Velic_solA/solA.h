#ifndef __StgFEM_Velic_solA_h__
#define __StgFEM_Velic_solA_h__

	extern const Type Velic_solA_Type;

	typedef struct {
		__AnalyticSolution
		double sigma;
		double Z;
	} Velic_solA;

	void Velic_solA_PressureFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* pressure );
	void Velic_solA_VelocityFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* velocity );
	void Velic_solA_StressFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* stress );
	void Velic_solA_StrainRateFunction( void* analyticSolution, FeVariable* analyticFeVariable, double* coord, double* strainRate );

#endif
