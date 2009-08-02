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
#include "PICellerator/PopulationControl/PopulationControl.h"
#include "PICellerator/Weights/Weights.h"

#include "ConstantWeightsSuite.h"
#include "MomentBalanceWeightsSuite.h"

typedef struct {
   DomainContext*       context;
} MomentBalanceWeightsSuiteData;


void MomentBalanceWeightsSuite_Setup( MomentBalanceWeightsSuiteData* data ) {
   char              xmlInputFilename[PCU_PATH_MAX];

   pcu_filename_input( "testMomentBalanceWeights.xml", xmlInputFilename );
   data->context = (DomainContext*)stgMainInitFromXML( xmlInputFilename, MPI_COMM_WORLD );
} 


void MomentBalanceWeightsSuite_Teardown( MomentBalanceWeightsSuiteData* data ) {
   stgMainDestroy( (AbstractContext*)data->context );
}


void MomentBalanceWeightsSuite_TestElementIntegral_Circle( MomentBalanceWeightsSuiteData* data ) {
   WeightsSuite_TestElementIntegral( data->context, "CircleInterface", 1000,
      0.000001, /* --mean-tolerance */
      0.000001, /* --standardDeviation-tolerance */
      0.097814, /* --mean-expectedValue */
      0.068607 /* --standardDeviation-expectedValue */ );
}


void MomentBalanceWeightsSuite_TestElementIntegral_Exponential( MomentBalanceWeightsSuiteData* data ) {
   WeightsSuite_TestElementIntegral( data->context, "ExponentialInterface", 1000,
      0.000001, /* --mean-tolerance */
      0.000001, /* --standardDeviation-tolerance */
      0.094671, /* --mean-expectedValue */
      0.075287 /* --standardDeviation-expectedValue */ );
}


void MomentBalanceWeightsSuite_TestElementIntegral_Polynomial( MomentBalanceWeightsSuiteData* data ) {
   WeightsSuite_TestElementIntegral( data->context, "PolynomialFunction", 1000,
      0.000001, /* --mean-tolerance */
      0.000001, /* --standardDeviation-tolerance */
      0.016697, /* --mean-expectedValue */
      0.013041 /* --standardDeviation-expectedValue */ );
}


void MomentBalanceWeightsSuite( pcu_suite_t* suite ) {
   pcu_suite_setData( suite, MomentBalanceWeightsSuiteData );
   pcu_suite_setFixtures( suite, MomentBalanceWeightsSuite_Setup, MomentBalanceWeightsSuite_Teardown );
   pcu_suite_addTest( suite, MomentBalanceWeightsSuite_TestElementIntegral_Circle );
   pcu_suite_addTest( suite, MomentBalanceWeightsSuite_TestElementIntegral_Exponential );
   pcu_suite_addTest( suite, MomentBalanceWeightsSuite_TestElementIntegral_Polynomial );
}
