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

/*******************
** Note: this test is currently designed to test the 3 components MaterialPointsSwarm, Material, and Materials_Register.
**  The reason is that their design and implementation is very tightly coupled (especially during Build and Init phases).
** 
********************/


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
   unsigned int   matProp1;
   double         matProp2;
   Bool           matProp3;
   Bool           matProp4;
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
   SwarmVariable*                matPropVar1;
   SwarmVariable*                matPropVar2;
   SwarmVariable*                matPropVar3;
   SwarmVariable*                matPropVar4;
} MaterialComponentsSuiteData;


FeMesh* buildFeMesh( unsigned nDims, unsigned* size, 
           double* minCrds, double* maxCrds, 
           ExtensionManager_Register* emReg )
{
   CartesianGenerator*   gen;
   FeMesh*         feMesh;

   gen = CartesianGenerator_New( "", NULL );
   CartesianGenerator_SetDimSize( gen, nDims );
   CartesianGenerator_SetTopologyParams( gen, size, 0, NULL, NULL );
   CartesianGenerator_SetGeometryParams( gen, minCrds, maxCrds );
   MeshGenerator_SetIncidenceState( gen, nDims, nDims, True );

   feMesh = FeMesh_New( "", NULL );
   Mesh_SetExtensionManagerRegister( feMesh, emReg );
   Mesh_SetGenerator( feMesh, gen );
   FeMesh_SetElementFamily( feMesh, "linear" );

   Stg_Component_Build( feMesh, NULL, False );
   Stg_Component_Initialise( feMesh, NULL, False );

   KillObject( feMesh->generator );

   return feMesh;
}


void MaterialComponentsSuite_BuildInitSwarmBasics( MaterialComponentsSuiteData* data ) {
   Index             var_I;
   Particle_Index    pI;
   MaterialPoint*    matPoint;

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
   /* Pre-set the swarm particle material indexes to undefined, as in the build function */
   for ( pI = 0; pI < data->mpSwarm->particleLocalCount; pI++ ) {
      matPoint = (MaterialPoint*)Swarm_ParticleAt( data->mpSwarm, pI );
      matPoint->materialIndex = UNDEFINED_MATERIAL;
   }
}


void MaterialComponentsSuite_Setup( MaterialComponentsSuiteData* data ) {
   unsigned    dim = 3;
   unsigned    meshSize[3] = {2, 2, 2};
   double      minCrds[3] = {0.0, 0.0, 0.0};
   double      maxCrds[3] = {1.0, 1.0, 1.0};
   XYZ         boxCentre = {0.25,0.25,0.25};
   XYZ         boxWidth = {0.4,0.4,0.4};
   XYZ         sphereCentre = {0.75,0.75,0.75};
   double      sphereRadius = 0.2;
   Particle_Index    pI;
   Particle*         matPoint;
   Particle          particle;

   data->svRegister = Variable_Register_New();
   data->eRegister = ExtensionManager_Register_New();
   data->mRegister = Materials_Register_New();

   /* Turn off the journal to avoid debug messages about parsing */
   Journal_Enable_TypedStream( Debug_Type, False );
   Journal_Enable_TypedStream( Info_Type, False );

   data->feMesh = buildFeMesh( dim, meshSize, minCrds, maxCrds, data->eRegister );

   data->cellLayout = (CellLayout*)ElementCellLayout_New( "elementCellLayout", NULL, data->feMesh );
   data->particleLayout = (ParticleLayout*)RandomParticleLayout_New( "randomParticleCellLayout", NULL, 
         GlobalCoordSystem, False, 
         20, 13 );

   data->mpSwarm = MaterialPointsSwarm_New(
         "testSwarm", NULL,
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
         MPI_COMM_WORLD,
         NULL );

   data->shape1 = (Stg_Shape*)Box_New( "boxShape", dim, boxCentre, 0, 0, 0, boxWidth );
   data->shape2 = (Stg_Shape*)Sphere_New( "sphereShape", dim, sphereCentre, 0, 0, 0, sphereRadius );

   /* Set up the dictionaries, and set some sample properties for testing */
   data->matDict1 = Dictionary_New();
   Dictionary_Add( data->matDict1, "testSwarm-matProp1", Dictionary_Entry_Value_FromUnsignedInt( 1 ) );
   Dictionary_Add( data->matDict1, "testSwarm-matProp2", Dictionary_Entry_Value_FromDouble( 1.1 ) );
   Dictionary_Add( data->matDict1, "testSwarm-matProp3", Dictionary_Entry_Value_FromBool( False ) );
   data->matDict2 = Dictionary_New();
   Dictionary_Add( data->matDict2, "testSwarm-matProp1", Dictionary_Entry_Value_FromUnsignedInt( 2 ) );
   Dictionary_Add( data->matDict2, "testSwarm-matProp2", Dictionary_Entry_Value_FromDouble( 2.2 ) );
   Dictionary_Add( data->matDict2, "testSwarm-matProp3", Dictionary_Entry_Value_FromBool( True ) );
   Dictionary_Add( data->matDict2, "testSwarm-matProp4", Dictionary_Entry_Value_FromBool( True ) );

   /* Now update the svRegister to match the material properties */
	data->matPropVar1 = Swarm_NewScalarVariable( 
			data->mpSwarm,
			"matProp1",
			GetOffsetOfMember( particle , matProp1 ), 
			Variable_DataType_Int );
	data->matPropVar2 = Swarm_NewScalarVariable( 
			data->mpSwarm,
			"matProp2",
			GetOffsetOfMember( particle , matProp2 ), 
			Variable_DataType_Double );
	data->matPropVar3 = Swarm_NewScalarVariable( 
			data->mpSwarm,
			"matProp3",
			GetOffsetOfMember( particle , matProp3 ), 
			Variable_DataType_Int );
	data->matPropVar4 = Swarm_NewScalarVariable( 
			data->mpSwarm,
			"matProp4",
			GetOffsetOfMember( particle , matProp4 ), 
			Variable_DataType_Int );

   data->mat1 = Material_New( "mat1", NULL, data->shape1, data->matDict1, data->mRegister );
   data->mat2 = Material_New( "mat2", NULL, data->shape2, data->matDict2, data->mRegister );

   MaterialComponentsSuite_BuildInitSwarmBasics( data );
   /* Now set them all the swarm var properties to 0 / False initially */
   for ( pI = 0; pI < data->mpSwarm->particleLocalCount; pI++ ) {
      matPoint = (Particle*)Swarm_ParticleAt( data->mpSwarm, pI );
      matPoint->matProp1 = 0;
      matPoint->matProp2 = 0.0;
      matPoint->matProp3 = False;
      matPoint->matProp4 = False;
   }
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

   pcu_docstring( "This test simply checks that the Material_Register associated to a test MP_Swarm is correctly set-up as part of the "
      "Swarm's Build phase (performed already in test suite setup." );

   pcu_check_true( data->mRegister != NULL );

   pcu_check_true( Materials_Register_GetCount( data->mRegister ) == 2 );
   pcu_check_true( Materials_Register_GetByIndex( data->mRegister, 0 ) == data->mat1 );
   pcu_check_true( Materials_Register_GetByIndex( data->mRegister, 1 ) == data->mat2 );


}


void MaterialComponentsSuite_TestMaterialLayout( MaterialComponentsSuiteData* data ) {
   Particle_Index    pI;
   MaterialPoint*    matPoint;

   pcu_docstring( "Tests that 2 given materials can be successfully 'layed out' - that is, the particles of a swarm have their materialIndex"
      "correctly updated to index the right material in the material register, thus setting up an ownership relationship." );
   
   for ( pI = 0; pI < data->mpSwarm->particleLocalCount; pI++ ) {
      matPoint = (MaterialPoint*)Swarm_ParticleAt( data->mpSwarm, pI );
      pcu_check_true( matPoint->materialIndex == UNDEFINED_MATERIAL );
   }

   Material_Layout( data->mat1, data->mpSwarm );
   Material_Layout( data->mat2, data->mpSwarm );

   for ( pI = 0; pI < data->mpSwarm->particleLocalCount; pI++ ) {
      matPoint = (MaterialPoint*)Swarm_ParticleAt( data->mpSwarm, pI );
      if ( Stg_Shape_IsCoordInside( data->shape1, matPoint->coord ) ) {
         pcu_check_true( matPoint->materialIndex == 0 );
      }
      else if ( Stg_Shape_IsCoordInside( data->shape2, matPoint->coord ) ) {
         pcu_check_true( matPoint->materialIndex == 1 );
      }
      else {
         pcu_check_true( matPoint->materialIndex == UNDEFINED_MATERIAL );
      }
   }
}


void MaterialComponentsSuite_TestLayoutGeometry( MaterialComponentsSuiteData* data ) {
   Particle_Index    pI;
   MaterialPoint*    matPoint;

   pcu_docstring( "Similar to TestMaterialLayout, checks that all materials in a register can be 'layed out' correctly - ie the individual"
      "Particles in a MaterialPointSwarm updated to index the correct material in the register." );
   _Materials_Register_LayoutGeometry( data->mRegister, data->mpSwarm );
   
   for ( pI = 0; pI < data->mpSwarm->particleLocalCount; pI++ ) {
      matPoint = (MaterialPoint*)Swarm_ParticleAt( data->mpSwarm, pI );
      if ( Stg_Shape_IsCoordInside( data->shape1, matPoint->coord ) ) {
         pcu_check_true( matPoint->materialIndex == 0 );
      }
      else if ( Stg_Shape_IsCoordInside( data->shape2, matPoint->coord ) ) {
         pcu_check_true( matPoint->materialIndex == 1 );
      }
      else {
         pcu_check_true( matPoint->materialIndex == UNDEFINED_MATERIAL );
      }
   }
}

void MaterialComponentsSuite_TestGetMaterial( MaterialComponentsSuiteData* data ) {
   Particle_Index    pI;
   MaterialPoint*    matPoint;
   Material*         material = NULL;
   Material*         materialOn = NULL;
   Material_Index    matIndex;

   pcu_docstring( "Tests that after a layout has been performed, the material that a particle belongs to can be successfully looked up."
      "Thus requires TestLayoutGeometry to pass beforehand." );

   _Materials_Register_LayoutGeometry( data->mRegister, data->mpSwarm );

   for ( pI = 0; pI < data->mpSwarm->particleLocalCount; pI++ ) {
      matPoint = (MaterialPoint*)Swarm_ParticleAt( data->mpSwarm, pI );
      material = MaterialPointsSwarm_GetMaterialAt( data->mpSwarm, pI );
      materialOn = MaterialPointsSwarm_GetMaterialOn( data->mpSwarm, matPoint );
      matIndex = MaterialPointsSwarm_GetMaterialIndexAt( data->mpSwarm, pI );

      if ( Stg_Shape_IsCoordInside( data->shape1, matPoint->coord ) ) {
         pcu_check_true( material == data->mat1 );
         pcu_check_true( materialOn == data->mat1 );
         pcu_check_true( matIndex == 0 );
      }
      else if ( Stg_Shape_IsCoordInside( data->shape2, matPoint->coord ) ) {
         pcu_check_true( material == data->mat2 );
         pcu_check_true( materialOn == data->mat2 );
         pcu_check_true( matIndex == 1 );
      }
      else {
         pcu_check_true( NULL == material );
         pcu_check_true( NULL == materialOn );
         pcu_check_true( UNDEFINED_MATERIAL == matIndex );
      }
   }
}


void MaterialComponentsSuite_TestAssignParticleProperties( MaterialComponentsSuiteData* data ) {
   Particle_Index    pI;
   Particle*         matPoint;

   pcu_docstring( "This tests the function Materials_Register_AssignParticleProperties: essentially, material "
      "properties that should be allowed to vary on a per-particle basis during simulation are set up as "
      "SwarmVariables, and this test ensures they are applied correctly at startup." );

   _Materials_Register_LayoutGeometry( data->mRegister, data->mpSwarm );
   Materials_Register_AssignParticleProperties( data->mRegister, data->mpSwarm, data->svRegister );

   /* check the properties have been assigned as expected - see Setup phase */
   for ( pI = 0; pI < data->mpSwarm->particleLocalCount; pI++ ) {
      matPoint = (Particle*)Swarm_ParticleAt( data->mpSwarm, pI );
      if ( Stg_Shape_IsCoordInside( data->shape1, matPoint->coord ) ) {
         pcu_check_true( matPoint->matProp1 == 1 );
         pcu_check_true( matPoint->matProp2 == 1.1 );
         pcu_check_true( matPoint->matProp3 == False );
         pcu_check_true( matPoint->matProp4 == False );
      }
      else if ( Stg_Shape_IsCoordInside( data->shape2, matPoint->coord ) ) {
         pcu_check_true( matPoint->matProp1 == 2 );
         pcu_check_true( matPoint->matProp2 == 2.2 );
         pcu_check_true( matPoint->matProp3 == True );
         pcu_check_true( matPoint->matProp4 == True );
      }
      else {
         pcu_check_true( matPoint->matProp1 == 0 );
         pcu_check_true( matPoint->matProp2 == 0.0 );
         pcu_check_true( matPoint->matProp3 == False );
         pcu_check_true( matPoint->matProp4 == False );
      }
   }
}

/* TODO M_Register: Test the extension functions work correctly */


void MaterialComponentsSuite( pcu_suite_t* suite ) {
   pcu_suite_setData( suite, MaterialComponentsSuiteData );
   pcu_suite_setFixtures( suite, MaterialComponentsSuite_Setup, MaterialComponentsSuite_Teardown );
   pcu_suite_addTest( suite, MaterialComponentsSuite_TestRegisterSetup );
   pcu_suite_addTest( suite, MaterialComponentsSuite_TestMaterialLayout );
   pcu_suite_addTest( suite, MaterialComponentsSuite_TestLayoutGeometry );
   pcu_suite_addTest( suite, MaterialComponentsSuite_TestGetMaterial );
   pcu_suite_addTest( suite, MaterialComponentsSuite_TestAssignParticleProperties );
}


