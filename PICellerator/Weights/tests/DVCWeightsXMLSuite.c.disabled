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
#include "DVCWeightsXMLSuite.h"

typedef struct {
   DomainContext*       context;
   Stg_ComponentFactory* cf;
} DVCWeightsXMLSuiteData;


void DVCWeightsXMLSuite_Setup( DVCWeightsXMLSuiteData* data ) {
   char              xmlInputFilename[PCU_PATH_MAX];

   pcu_filename_input( "testDVCWeights.xml", xmlInputFilename );
   data->cf = stgMainInitFromXML( xmlInputFilename, MPI_COMM_WORLD, NULL );
   data->context = LiveComponentRegister_Get( data->cf, "context" );
} 


void DVCWeightsXMLSuite_Teardown( DVCWeightsXMLSuiteData* data ) {
   stgMainDestroy( data->cf );
}


void DVCWeightsXMLSuite_TestElementIntegral_Circle( DVCWeightsXMLSuiteData* data ) {
   WeightsSuite_TestElementIntegral( data->context, "CircleInterface", 1000,
      0.005, /* --mean-tolerance */
      0.001, /* --standardDeviation-tolerance */
      0.090837, /* --mean-expectedValue */
      0.055628 /* --standardDeviation-expectedValue */ );
}


void DVCWeightsXMLSuite_TestElementIntegral_Exponential( DVCWeightsXMLSuiteData* data ) {
   WeightsSuite_TestElementIntegral( data->context, "ExponentialInterface", 1000,
      0.001, /* --mean-tolerance */
      0.001, /* --standardDeviation-tolerance */
      0.047945, /* --mean-expectedValue */
      0.035636 /* --standardDeviation-expectedValue */ );
}


void DVCWeightsXMLSuite_TestElementIntegral_Polynomial( DVCWeightsXMLSuiteData* data ) {
   WeightsSuite_TestElementIntegral( data->context, "PolynomialFunction", 1000,
      0.001, /* --mean-tolerance */
      0.001, /* --standardDeviation-tolerance */
      0.0077474, /* --mean-expectedValue */
      0.0014641 /* --standardDeviation-expectedValue */ );
}


void DVCWeightsXMLSuite( pcu_suite_t* suite ) {
   pcu_suite_setData( suite, DVCWeightsXMLSuiteData );
   pcu_suite_setFixtures( suite, DVCWeightsXMLSuite_Setup, DVCWeightsXMLSuite_Teardown );
   pcu_suite_addTest( suite, DVCWeightsXMLSuite_TestElementIntegral_Circle );
/* TEMPORARILY disable multiple tests until Toolbox/Context issue sorted out (See tickets #70,#71
 in StgNumerics trac) -- PatrickSunter, 19 Aug 2009 */
#if 0
   pcu_suite_addTest( suite, DVCWeightsXMLSuite_TestElementIntegral_Exponential );
   pcu_suite_addTest( suite, DVCWeightsXMLSuite_TestElementIntegral_Polynomial );
#endif
}
