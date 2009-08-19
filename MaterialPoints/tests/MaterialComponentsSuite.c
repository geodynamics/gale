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
#include "PICellerator/MaterialPoints/MaterialPoints.h"

struct _Node {
   Coord            coord;
};

struct _Element {
   Coord            coord;
};

struct _Particle {
   __MaterialPoint
};

typedef struct {
   FeMesh*                       feMesh;
   CellLayout*                   cellLayout;
   ParticleLayout*               particleLayout;
   MaterialPointsSwarm*          mpSwarm;
   Variable_Register*            svRegister;
   ExtensionManager_Register*    eRegister;
   Materials_Register*           mRegister;
   Dictionary*                   matDict1;
   Dictionary*                   matDict2;
   Stg_Shape*                    shape1;
   Stg_Shape*                    shape2;
   Material*                     mat1;
   Material*                     mat2;
} MaterialComponentsSuiteData;


FeMesh* buildFeMesh( unsigned nDims, unsigned* size, 
           double* minCrds, double* maxCrds, 
           ExtensionManager_Register* emReg )
{
   CartesianGenerator*   gen;
   FeMesh*         feMesh;

   gen = CartesianGenerator_New( "" );
   CartesianGenerator_SetDimSize( gen, nDims );
   CartesianGenerator_SetTopologyParams( gen, size, 0, NULL, NULL );
   CartesianGenerator_SetGeometryParams( gen, minCrds, maxCrds );
   MeshGenerator_SetIncidenceState( gen, nDims, nDims, True );

   feMesh = FeMesh_New( "" );
   Mesh_SetExtensionManagerRegister( feMesh, emReg );
   Mesh_SetGenerator( feMesh, gen );
   FeMesh_SetElementFamily( feMesh, "linear" );

   Stg_Component_Build( feMesh, NULL, False );
   Stg_Component_Initialise( feMesh, NULL, False );

   KillObject( feMesh->generator );

   return feMesh;
}


void MaterialComponentsSuite_BuildInitSwarmBasics( MaterialComponentsSuiteData* data ) {
   Index    var_I;

   /* Mesh will already be initialised */
   /* Manually do the basic Swarm building and initialising */
   _Swarm_Build( data->mpSwarm, NULL );
   for( var_I = 0 ; var_I < data->mpSwarm->nSwarmVars ; var_I++ ) {
		Stg_Component_Build( data->mpSwarm->swarmVars[var_I], NULL , False );
	}
   _Swarm_Initialise( data->mpSwarm, NULL );
   for( var_I = 0 ; var_I < data->mpSwarm->nSwarmVars ; var_I++ ) {
		Stg_Component_Initialise( data->mpSwarm->swarmVars[var_I], NULL , False );
	}
}


void MaterialComponentsSuite_Setup( MaterialComponentsSuiteData* data ) {
   unsigned    dim = 3;
   unsigned    meshSize[3] = {3, 3, 3};
   double      minCrds[3] = {0.0, 0.0, 0.0};
   double      maxCrds[3] = {1.0, 1.0, 1.0};
   XYZ         boxCentre = {0.25,0.25,0.25};
   XYZ         boxWidth = {0.1,0.1,0.1};
   XYZ         sphereCentre = {0.75,0.75,0.75};
   double      sphereRadius = 0.1;

   data->svRegister = Variable_Register_New();
   data->eRegister = ExtensionManager_Register_New();
   data->mRegister = Materials_Register_New();

   data->feMesh = buildFeMesh( dim, meshSize, minCrds, maxCrds, data->eRegister );

   data->cellLayout = (CellLayout*)ElementCellLayout_New( "elementCellLayout", data->feMesh );
   data->particleLayout = (ParticleLayout*)RandomParticleLayout_New( "randomParticleCellLayout", 4, 13 );

   data->mpSwarm = MaterialPointsSwarm_New(
      "testSwarm",
      data->cellLayout,
      data->particleLayout,
      dim,                      /* dim */
      sizeof(Particle),
      data->feMesh,
      NULL,                   /* escapedRoutine*/
      NULL,                   /* material */
      data->svRegister,
      data->eRegister,
      data->mRegister,
      MPI_COMM_WORLD );

   data->shape1 = (Stg_Shape*)Box_New( "boxShape", dim, boxCentre, 0, 0, 0, boxWidth );
   data->shape2 = (Stg_Shape*)Sphere_New( "sphereShape", dim, sphereCentre, 0, 0, 0, sphereRadius );
   data->matDict1 = Dictionary_New();
   data->matDict2 = Dictionary_New();
   data->mat1 = Material_New( "mat1", data->shape1, data->matDict1, data->mRegister );
   data->mat2 = Material_New( "mat2", data->shape2, data->matDict2, data->mRegister );

   MaterialComponentsSuite_BuildInitSwarmBasics( data );
}


void MaterialComponentsSuite_Teardown( MaterialComponentsSuiteData* data ) {
   Stg_Class_Delete( data->mpSwarm );
   Stg_Class_Delete( data->cellLayout );
   Stg_Class_Delete( data->particleLayout );
   Stg_Class_Delete( data->feMesh );
   Stg_Class_Delete( data->eRegister );
   Stg_Class_Delete( data->svRegister );
   Stg_Class_Delete( data->mRegister );
   Stg_Class_Delete( data->mat1 );
   Stg_Class_Delete( data->mat2 );
   Stg_Class_Delete( data->matDict1 );
   Stg_Class_Delete( data->matDict2 );
   Stg_Class_Delete( data->shape1 );
   Stg_Class_Delete( data->shape2 );
}


/*Tests*/
void MaterialComponentsSuite_TestRegisterSetup( MaterialComponentsSuiteData* data ) {

   pcu_check_true( data->mRegister != NULL );
   
}


// Test a basic MP swarm can be created, with an MP register pointing to a set of materials
// Test the basic access functions: eg the MPS_GetMaterialOn
// Material: Test the M_Layout function works correctly
// M_Register: Test the multiple layout and extension functions work correctly. 


void MaterialComponentsSuite( pcu_suite_t* suite ) {
   pcu_suite_setData( suite, MaterialComponentsSuiteData );
   pcu_suite_setFixtures( suite, MaterialComponentsSuite_Setup, MaterialComponentsSuite_Teardown );
   pcu_suite_addTest( suite, MaterialComponentsSuite_TestRegisterSetup );
}
