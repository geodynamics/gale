#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>

#ifndef __AnalyticPressure_h__
#define __AnalyticPressure_h__

	extern const Type AnalyticPressure_Type;

	typedef struct {
		__FieldTest
		double density;
		double gravity;
		double maxY;
		double minY;
	} AnalyticPressure;

	void AnalyticPressure_PressureFunction( void* _self, double* coord, double* pressure );

#endif
