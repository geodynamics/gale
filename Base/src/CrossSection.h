#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>

#ifndef __lucCrossSection_h__
#define __lucCrossSection_h__

/* Cross section utility */
typedef struct {
  		double   value;
		Axis     axis;
      Bool     interpolate;
} lucCrossSection;

/* Returns a cross section struct parsed from XML string */
lucCrossSection* lucCrossSection_Read( Stg_ComponentFactory* cf, Name component);

/* Setup cross section values from passed parameters 
 * Returns pointer to passed in cross section so can be used in function calls 
 * If input object is NULL a new one is created */
lucCrossSection* lucCrossSection_Set(lucCrossSection* self, double value, Axis axis, Bool interpolate);

/* Returns the cross section value, interpolating where necessary */
double lucCrossSection_GetValue(lucCrossSection* self, double min, double max);

/* Free cross-section memory */
void lucCrossSection_Delete(lucCrossSection* self);

#endif
