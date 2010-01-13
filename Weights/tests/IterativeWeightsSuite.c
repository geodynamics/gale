/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street, Melbourne, 3053, Australia.
**
** Authors:
**   Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**   Patrick D. Sunter, Software Engineer, VPAC. (pds@vpac.org)
**   Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**   Siew-Ching Tan, Software Engineer, VPAC. (siew@vpac.org)
**   Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
**   Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
**
**  This library is free software; you can redistribute it and/or
**  modify it under the terms of the GNU Lesser General Public
**  License as published by the Free Software Foundation; either
**  version 2.1 of the License, or (at your option) any later version.
**
**  This library is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
**  Lesser General Public License for more details.
**
**  You should have received a copy of the GNU Lesser General Public
**  License along with this library; if not, write to the Free Software
**  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
**
** $Id: testList.c 2136 2004-09-30 02:47:13Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pcu/pcu.h"
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>

#include "ConstantWeightsSuite.h"
#include "IterativeWeightsSuite.h"

typedef struct {
   PICelleratorContext*       context;
   Stg_ComponentFactory* cf;
} IterativeWeightsSuiteData;


void IterativeWeightsSuite_Setup( IterativeWeightsSuiteData* data ) {
   char              xmlInputFilename[PCU_PATH_MAX];

   pcu_filename_input( "testIterativeWeights.xml", xmlInputFilename );
   data->cf = stgMainInitFromXML( xmlInputFilename, MPI_COMM_WORLD, NULL );
   data->context = (PICelleratorContext*)LiveComponentRegister_Get( data->cf->LCRegister, (Name)"context" );
   stgMainBuildAndInitialise( data->cf );
} 


void IterativeWeightsSuite_Teardown( IterativeWeightsSuiteData* data ) {
   stgMainDestroy( data->cf );
}

void IterativeWeightsSuite_TestConstantFunction( IterativeWeightsSuiteData* data  ) {
   WeightsSuite_TestElementIntegral( data->context, "ConstantFunction", 1000,
      1e-10, /* --mean-tolerance */
      1e-10, /* --standardDeviation-tolerance */
      0.0, /* --mean-expectedValue */
      0.0 /* --standardDeviation-expectedValue */ );
}
void IterativeWeightsSuite_TestLinearFunction ( IterativeWeightsSuiteData* data ) {
   WeightsSuite_TestElementIntegral( data->context, "LinearFunction", 1000,
      1e-4, /* --mean-tolerance */
      1e-4, /* --standardDeviation-tolerance */
      0.0, /* --mean-expectedValue */
      0.0 /* --standardDeviation-expectedValue */ );
}
void IterativeWeightsSuite_TestQuadraticFunction ( IterativeWeightsSuiteData* data ) {
   WeightsSuite_TestElementIntegral( data->context, "QuadraticFunction", 1000,
      0.000001, /* --mean-tolerance */
      0.000001, /* --standardDeviation-tolerance */
      0.0430721, /* --mean-expectedValue */
      0.0326016 /* --standardDeviation-expectedValue */ );
}

void IterativeWeightsSuite_TestPolynomialFunction( IterativeWeightsSuiteData* data ) {
   WeightsSuite_TestElementIntegral( data->context, "PolynomialFunction", 1000,
      0.000001, /* --mean-tolerance */
      0.000001, /* --standardDeviation-tolerance */
      0.0175259, /* --mean-expectedValue */
      0.013522 /* --standardDeviation-expectedValue */ );
}

void IterativeWeightsSuite_TestCircleInterface( IterativeWeightsSuiteData* data ) {
   WeightsSuite_TestElementIntegral( data->context, "CircleInterface", 1000,
      0.000001, /* --mean-tolerance */
      0.000001, /* --standardDeviation-tolerance */
      0.10172, /* --mean-expectedValue */
      0.070065 /* --standardDeviation-expectedValue */ );
}

void IterativeWeightsSuite_TestExponentialInterface( IterativeWeightsSuiteData* data ) {
   WeightsSuite_TestElementIntegral( data->context, "ExponentialInterface", 1000,
      0.000001, /* --mean-tolerance */
      0.000001, /* --standardDeviation-tolerance */
      0.088927, /* --mean-expectedValue */
      0.06681 /* --standardDeviation-expectedValue */ );
}

void IterativeWeightsSuite( pcu_suite_t* suite ) {
   pcu_suite_setData( suite, IterativeWeightsSuiteData );
   pcu_suite_setFixtures( suite, IterativeWeightsSuite_Setup, IterativeWeightsSuite_Teardown );
   pcu_suite_addTest( suite, IterativeWeightsSuite_TestConstantFunction );
   pcu_suite_addTest( suite, IterativeWeightsSuite_TestLinearFunction );
   pcu_suite_addTest( suite, IterativeWeightsSuite_TestQuadraticFunction );
   pcu_suite_addTest( suite, IterativeWeightsSuite_TestPolynomialFunction );
   pcu_suite_addTest( suite, IterativeWeightsSuite_TestCircleInterface );
   pcu_suite_addTest( suite, IterativeWeightsSuite_TestExponentialInterface );
}


