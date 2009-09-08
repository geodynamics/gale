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
**
** $Id: testJournal-Dictionary.c 2745 2005-03-05 08:12:18Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pcu/pcu.h"
#include "StGermain/Base/Foundation/Foundation.h"
#include "StGermain/Base/IO/IO.h"
#include "StGermain/Base/Container/Container.h"
#include "StGermain/Base/Automation/Automation.h"
#include "SetVC_Suite.h"

/* Note that these params below must match up with setVC.xml input file for the test to work */
#define NUM_VARS 7
#define X_INDEX 0
#define Y_INDEX 1
#define Z_INDEX 2
#define VX_INDEX 3
#define VY_INDEX 4
#define VZ_INDEX 5
#define TEMP_INDEX 6

#define NUM_NODES 64
#define LAST_VC_NODE_INDEX 26
#define TEMP_VAR_ARRAYSIZE 5

#define VX_CONDVALUE 2
#define VY_CONDVALUE 20.0
#define VZ_CONDVALUE 1
#define TEMP_CONDVALUE_BASE 5


typedef struct {
   Variable_Register*            vr;
   ConditionFunction_Register*   conFunc_Register;
   Dictionary*                   dict;
   ConditionFunction*            quadCF;
   double*                       array[NUM_VARS];
   Index                         arraySize;
   char*                         varNames[NUM_VARS];
   Variable*                     var[NUM_VARS];
   VariableCondition*            vc;
} SetVC_SuiteData;


void quadratic(Index index, Variable_Index var_I, void* context, void* result)
{
	*(double *)result = 20.0;
}


void SetVC_Suite_Setup( SetVC_SuiteData* data ) {
   Index             ii=0;
   char              setVC_XMLFilename[PCU_PATH_MAX];
   XML_IO_Handler*   io_handler = XML_IO_Handler_New();

   data->vr = Variable_Register_New();
   data->conFunc_Register = ConditionFunction_Register_New();
   data->dict = Dictionary_New();
   data->arraySize = NUM_NODES;

   Stg_asprintf( &data->varNames[X_INDEX], "x" );
   Stg_asprintf( &data->varNames[Y_INDEX], "y" );
   Stg_asprintf( &data->varNames[Z_INDEX], "z" );
   Stg_asprintf( &data->varNames[VX_INDEX], "vx" );
   Stg_asprintf( &data->varNames[VY_INDEX], "vy" );
   Stg_asprintf( &data->varNames[VZ_INDEX], "vz" );
   Stg_asprintf( &data->varNames[TEMP_INDEX], "temp" );

   /* Load the definition of the SetVC from XML */
   pcu_filename_input( "setVC.xml", setVC_XMLFilename );
   IO_Handler_ReadAllFromFile( io_handler, setVC_XMLFilename, data->dict );
   
   /* Create CF stuff */
   data->quadCF = ConditionFunction_New( quadratic, "quadratic" );
   ConditionFunction_Register_Add( data->conFunc_Register, data->quadCF );
   
   /* Create variables */
   for (ii = 0; ii < (TEMP_INDEX); ii++) {
      data->array[ii] = Memory_Alloc_Array( double, data->arraySize, "data->array[]" );
      data->var[ii] = Variable_NewScalar( data->varNames[ii], Variable_DataType_Double, &data->arraySize, NULL,
         (void**)&data->array[ii], 0 ); 
      Variable_Register_Add(data->vr, data->var[ii]);
   }

   data->array[TEMP_INDEX] = Memory_Alloc_Array( double , data->arraySize*TEMP_VAR_ARRAYSIZE, "data->array[TEMP_INDEX]" );
   data->var[TEMP_INDEX] = Variable_NewVector( data->varNames[TEMP_INDEX], Variable_DataType_Double, TEMP_VAR_ARRAYSIZE,
       &data->arraySize, NULL, (void**)&data->array[TEMP_INDEX], 0, "a", "b", "c", "d", "e" );
   Variable_Register_Add( data->vr, data->var[TEMP_INDEX] );

   Variable_Register_BuildAll( data->vr );
   
   /* Create the VC */
   data->vc = (VariableCondition*)SetVC_New( "setVC", "setVC", data->vr, data->conFunc_Register, data->dict );
   Stg_Component_Build( data->vc, 0, False );

   /* Blank the memory to be applied to before the test */
   for (ii = 0; ii < (TEMP_INDEX); ii++) {
      memset(data->array[ii], 0, sizeof(double)*NUM_NODES);
   }
   memset(data->array[TEMP_INDEX], 0, sizeof(double)*NUM_NODES*TEMP_VAR_ARRAYSIZE);

   Memory_Free( io_handler );
}


void SetVC_Suite_Teardown( SetVC_SuiteData* data ) {
   Index    ii=0;

   Stg_Class_Delete(data->vr);
   Stg_Class_Delete(data->conFunc_Register);
   Stg_Class_Delete(data->dict);

   Stg_Class_Delete(data->quadCF);

   for (ii = 0; ii < NUM_VARS; ii++) {
      if (data->array[ii]) Memory_Free(data->array[ii]);
      Stg_Class_Delete( data->var[ii] );
   }

   Stg_Class_Delete(data->vc);
}


void SetVC_Suite_TestIsCondition( SetVC_SuiteData* data ) {
   Index      var_I=0, node_I=0;

   /* remember varName[] = {"x", "y", "z", "vx", "vy", "vz", "temp"}*/

   /* X, Y, and Z shouldn't have Conditions applied to them */
   for( var_I = 0; var_I < VX_INDEX; var_I++ ) {
      for (node_I = 0; node_I < NUM_NODES; node_I++) {
         pcu_check_true( False == VariableCondition_IsCondition(data->vc, node_I, var_I) );
      }
   }
   /* VX, VY, VZ, and temp should have Conditions applied to them, at even indices */
   for( var_I = VX_INDEX; var_I < NUM_VARS; var_I++ ) {
      for (node_I = 0; node_I < NUM_NODES; node_I++) {
         if ( (node_I % 2 == 0) && node_I <= LAST_VC_NODE_INDEX ) {
            pcu_check_true( True == VariableCondition_IsCondition(data->vc, node_I, var_I) );
         }
         else {
            pcu_check_true( False == VariableCondition_IsCondition(data->vc, node_I, var_I) );
         }
      }
   }
}


void SetVC_Suite_TestGetValueIndex( SetVC_SuiteData* data ) {
   Index      var_I=0, node_I=0;

   /* remember varName[] = {"x", "y", "z", "vx", "vy", "vz", "temp"}*/
   /* These results are dependent on the order the conditions are added in the dictionary */

   for( var_I = 0; var_I < VX_INDEX; var_I++ ) {
      for (node_I = 0; node_I < NUM_NODES; node_I++) {
         pcu_check_true( (unsigned)-1 == VariableCondition_GetValueIndex(data->vc, node_I, var_I) );
      }
   }
   for( var_I = VX_INDEX; var_I < NUM_VARS; var_I++ ) {
      for (node_I = 0; node_I < NUM_NODES; node_I++) {
         if ( (node_I % 2 == 0) && node_I <= LAST_VC_NODE_INDEX ) {
            /* Condition indices will start at 0, starting at condition 3 "vx" */
            pcu_check_true( (var_I-VX_INDEX) == VariableCondition_GetValueIndex(data->vc, node_I, var_I) );
         }
         else {
            pcu_check_true( (unsigned)-1 == VariableCondition_GetValueIndex(data->vc, node_I, var_I) );
         }
      }
   }
}


void SetVC_Suite_TestApply( SetVC_SuiteData* data ) {
   Index      var_I=0, node_I=0, array_I=0;
   
   /* remember varName[] = {"x", "y", "z", "vx", "vy", "vz", "temp"}; */

   VariableCondition_Apply(data->vc, NULL);
   
   for( var_I = 0; var_I < VX_INDEX; var_I++ ) {
      for (node_I = 0; node_I < NUM_NODES; node_I++) {
         pcu_check_true( 0 == data->array[var_I][node_I] );
      }
   }
   for( var_I = VX_INDEX; var_I < TEMP_INDEX; var_I++ ) {
      for (node_I = 0; node_I < NUM_NODES; node_I++) {
         if ( (node_I % 2 == 0) && node_I <= LAST_VC_NODE_INDEX ) {
            if (var_I == VX_INDEX) {
               pcu_check_true( VX_CONDVALUE == data->array[var_I][node_I] );
            }
            else if (var_I == VY_INDEX) {
               pcu_check_true( VY_CONDVALUE == data->array[var_I][node_I] );
            }
            else if (var_I == VZ_INDEX) {
               pcu_check_true( VZ_CONDVALUE == data->array[var_I][node_I] );
            }
         }
         else {
            pcu_check_true( 0 == data->array[var_I][node_I] );
         }
      }
   }
   /* Now for the "temp" variable */
   for (node_I = 0; node_I < NUM_NODES; node_I++) {
      for( array_I = 0; array_I < TEMP_VAR_ARRAYSIZE; array_I++ ) {
         if ( (node_I % 2 == 0) && node_I <= LAST_VC_NODE_INDEX ) {
            pcu_check_true( (TEMP_CONDVALUE_BASE+array_I) == data->array[TEMP_INDEX][node_I*TEMP_VAR_ARRAYSIZE + array_I] );
         }
         else {
            pcu_check_true( 0 == data->array[TEMP_INDEX][node_I*TEMP_VAR_ARRAYSIZE + array_I] );
         }
      }
   }
}


void SetVC_Suite( pcu_suite_t* suite ) {
   pcu_suite_setData( suite, SetVC_SuiteData );
   pcu_suite_setFixtures( suite, SetVC_Suite_Setup, SetVC_Suite_Teardown );
   pcu_suite_addTest( suite, SetVC_Suite_TestIsCondition );
   pcu_suite_addTest( suite, SetVC_Suite_TestGetValueIndex );
   pcu_suite_addTest( suite, SetVC_Suite_TestApply );
}
