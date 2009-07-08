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
#include "PICellerator/PICellerator.h"

typedef struct {
   PICelleratorContext* context;
} MaterialFeVariableSuiteData;

void MaterialFeVariableSuite_Setup( MaterialFeVariableSuiteData* data ) {
   char              xmlInputFilename[PCU_PATH_MAX];
   Dictionary*       dictionary = NULL;
   Bool              result;
   XML_IO_Handler*   ioHandler;

   dictionary = Dictionary_New();
   ioHandler = XML_IO_Handler_New();
   pcu_filename_input( "testMaterialFeVariable.xml", xmlInputFilename );
   result = IO_Handler_ReadAllFromFile( ioHandler, xmlInputFilename, dictionary );
   Journal_ReadFromDictionary( dictionary );

   /* Magic happens here! */
   data->context = (PICelleratorContext*)stgMainInit( dictionary, MPI_COMM_WORLD );
   stgMainLoop( (AbstractContext*)data->context );
} 


void MaterialFeVariableSuite_Teardown( MaterialFeVariableSuiteData* data ) {
   stgMainDestroy( (AbstractContext*)data->context );
}


void MaterialFeVariableSuite_TestVolume( MaterialFeVariableSuiteData* data ) {
/** Test Definition: This test compares the volume of a swarm as calculated based on the material points swarm and as calculate from a MaterialFeVariable */
   IntegrationPointsSwarm* swarm;
   MaterialFeVariable* materialFeVariable;
   Material*           material;
   Swarm*              gaussSwarm;
   double              volumePIC;
   double              volumeFEM;
   Coord               centroid;

   char   testStr[] = "This test compares the volume of a swarm as calculated based on the material points swarm and as calculate from a MaterialFeVariable\n";

   materialFeVariable = (MaterialFeVariable*) LiveComponentRegister_Get( data->context->CF->LCRegister, "materialFeVariable" );
   gaussSwarm = (Swarm*) LiveComponentRegister_Get( data->context->CF->LCRegister, "gaussSwarm" );

   assert(materialFeVariable);
   assert(gaussSwarm);

   material = materialFeVariable->material;
   swarm = materialFeVariable->picIntegrationPoints;

   volumePIC = Material_Volume( material, swarm, centroid );
   volumeFEM = FeVariable_Integrate( materialFeVariable, gaussSwarm );

   pcu_check_true(  fabs(volumePIC-volumeFEM) < 1e-7 );
}

void MaterialFeVariableSuite( pcu_suite_t* suite ) {
   pcu_suite_setData( suite, MaterialFeVariableSuiteData );
   pcu_suite_setFixtures( suite, MaterialFeVariableSuite_Setup, MaterialFeVariableSuite_Teardown );
   pcu_suite_addTest( suite, MaterialFeVariableSuite_TestVolume );
}
