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

#include "pcu/pcu.h"
#include "StGermain/Base/Foundation/Foundation.h"
#include "StGermain/Base/IO/IO.h"
#include "StGermain/Base/Container/Container.h"
#include "StGermain/Base/Automation/Automation.h"
#include "VariableSuite.h"

typedef double VectorD[3];
typedef double VectorF[3];

typedef struct {
   int      mass;
   VectorF  force;
   short    num;
   char*    info;
} Particle;

typedef struct {
   Variable_Register*      vr;
   Index                   aSize[3];
   double*                 temperature;
   VectorD*                velocity;
   Particle*               particle;
} VariableSuiteData;



void VariableSuite_Setup( VariableSuiteData* data ) {
   Particle                tmpParticle;
   Name                    pNames[] = { "mass", "force", "info" };
   SizeT                   pOffsets[] = { 0, 0, 0 };   /* Init later... */
   Variable_DataType       pDataTypes[] = {
                              Variable_DataType_Int,
                              Variable_DataType_Float,
                              Variable_DataType_Pointer, };
   Index                   pDtCounts[] = { 1, 3, 1 };
   static SizeT            pSize = sizeof(Particle);
   
   data->aSize[0] = 16;
   data->aSize[1] = 16;
   data->aSize[2] = 16;

   pOffsets[0] = (ArithPointer)&tmpParticle.mass - (ArithPointer)&tmpParticle;
   pOffsets[1] = (ArithPointer)&tmpParticle.force - (ArithPointer)&tmpParticle;
   pOffsets[2] = (ArithPointer)&tmpParticle.info - (ArithPointer)&tmpParticle;

   /* Construction phase --------------------------------------------------------------------------------------------*/
   data->vr = Variable_Register_New();
   Variable_NewScalar( "temperature", Variable_DataType_Double, &data->aSize[0], NULL, (void**)&data->temperature, data->vr );
   Variable_NewVector( "velocity", Variable_DataType_Double, 3, &data->aSize[1], NULL, (void**)&data->velocity, data->vr, "vx", "vy", "vz" );
   Variable_New( "particle", 3, pOffsets, pDataTypes, pDtCounts, pNames, &pSize, &data->aSize[2], NULL, (void**)&data->particle, data->vr );
   
   /* Build phase ---------------------------------------------------------------------------------------------------*/
   data->temperature = Memory_Alloc_Array( double, data->aSize[0], "temperature" );
   data->velocity = Memory_Alloc_Array( VectorD, data->aSize[1], "velocity" );
   data->particle = Memory_Alloc_Array( Particle, data->aSize[2], "array" );
   
   Variable_Register_BuildAll( data->vr );
}


void VariableSuite_Teardown( VariableSuiteData* data ) {
   Variable_Index          var_I;

   /* manually delete all the created Variables */
   for( var_I = 0; var_I < data->vr->count; var_I++ ) {
      Stg_Class_Delete( data->vr->_variable[var_I] );
   }
   
   Memory_Free( data->particle );
   Memory_Free( data->velocity );
   Memory_Free( data->temperature );
}


void VariableSuite_TestGetValueDouble( VariableSuiteData* data ) {
   Index                   ii;

   /* Test the Get and Set of a scalar double....................................................................... */
   /* Fill the temperature array with a known pattern of kinda random (bit filling) numbers. */
   for( ii = 0; ii < data->aSize[0]; ii++ ) {
      data->temperature[ii] = 1.0f / (data->aSize[0]+2) * (ii+1); 
   }
   
   /* Check that Variable_GetValueDouble on the temperature Variable returns the right numbers */
   for( ii = 0; ii < data->aSize[0]; ii++ ) {
      Variable*      var = Variable_Register_GetByName( data->vr, "temperature" );
      const double      tmp = 1.0f / (data->aSize[0]+2) * (ii+1);
      
      pcu_check_true( Variable_GetValueDouble( var, ii ) == tmp );
   }
}
   

void VariableSuite_TestSetValueDouble( VariableSuiteData* data ) {
   Index                   ii;

   /* Fill the temperature Variable with another known pattern of kinda random (bit filling) numbers */
   for( ii = 0; ii < data->aSize[0]; ii++ ) {
      Variable*      var = Variable_Register_GetByName( data->vr, "temperature" );
      
      Variable_SetValueDouble( var, ii, 1.0f - ( 1.0f / (data->aSize[0]+2) * (ii+1) ) );
   }
   
   /* Check that Variable_SetValueDouble on the temperature Variable set the right numbers */
   for( ii = 0; ii < data->aSize[0]; ii++ ) {
      const double      tmp = 1.0f - 1.0f / (data->aSize[0]+2) * (ii+1);
      
      pcu_check_true( data->temperature[ii] == tmp );
   }
}
   

/* Test the Get and Set of a vector double....................................................................... */
void VariableSuite_TestGetValueAtDouble( VariableSuiteData* data ) {
   Index                   ii;

/* Fill the velocity array with a known pattern of kinda random (bit filling) numbers. */
   for( ii = 0; ii < data->aSize[1]; ii++ ) {
      int         d;
      
      for( d = 0; d < 3; d++ ) {
         data->velocity[ii][d] = 1.0f / ((data->aSize[1]*3)+2) * (ii*3+d+1); 
      }
   }
   
   /* Check that Variable_GetPtrDouble on the velocity Variable returns the right numbers */
   for( ii = 0; ii < data->aSize[1]; ii++ ) {
      Variable*      var = Variable_Register_GetByName( data->vr, "velocity" );
      int         d;
      
      for( d = 0; d < 3; d++ ) {
         const double       tmp = 1.0f / ((data->aSize[1]*3)+2) * (ii*3+d+1);
         
         pcu_check_true( Variable_GetValueAtDouble( var, ii, d ) == tmp );
      }
   }
}


void VariableSuite_TestSetValueAtDouble( VariableSuiteData* data ) {
   Index                   ii;

   /* Fill the variable Variable with another known pattern of kinda random (bit filling) numbers */
   for( ii = 0; ii < data->aSize[1]; ii++ ) {
      Variable*      var = Variable_Register_GetByName( data->vr, "velocity" );
      int         d;
      
      for( d = 0; d < 3; d++ ) {
         Variable_SetValueAtDouble( var, ii, d, 1.0f - ( 1.0f / ((data->aSize[1]*3)+2) * (ii*3+d+1) ) );
      }
   }
   
   /* Check that Variable_SetValueDouble on the temperature Variable set the right numbers */
   for( ii = 0; ii < data->aSize[1]; ii++ ) {
      int         d;
      
      for( d = 0; d < 3; d++ ) {
         const double      tmp = 1.0f - ( 1.0f / ((data->aSize[1]*3)+2) * (ii*3+d+1) );
         
         pcu_check_true( data->velocity[ii][d] == tmp );
      }
   }
}

   
/* TODO: try out vx, vy, vz, complex tests */

void VariableSuite_TestVariable_Char( VariableSuiteData* data ) {
   typedef char Triple[3];
   char* array;
   Triple* structArray;
   Index length = 10;
   /* List of values to test the variable with.
    * Values to test are hex 5's and a's because they are a series of 0101 and 1010 respectively so they test
    * each bit in memory to read/set.
    */
   long int testValues[] = { 0x55, 0xaa };
   Index testValueCount = 2;
   Index test_I;
   long int testValue;
   Variable* var;
   Variable* vec;
   Variable* vecVar[3];
   int i, j;

   array = Memory_Alloc_Array( char, length, "test" );
   structArray = Memory_Alloc_Array( Triple, length, "test" );

   var = Variable_NewScalar(
      "Scalar",
      Variable_DataType_Char,
      &length,
      NULL,
      (void**)&array,
      data->vr );

   vec = Variable_NewVector(
      "Three",
      Variable_DataType_Char,
      3,
      &length,
      NULL,
      (void**)&structArray,
      data->vr,
      "a",
      "b",
      "c" );
   vecVar[0] = Variable_Register_GetByName( data->vr, "a" );
   vecVar[1] = Variable_Register_GetByName( data->vr, "b" );
   vecVar[2] = Variable_Register_GetByName( data->vr, "c" );

   Variable_Register_BuildAll( data->vr );

   for ( test_I = 0; test_I < testValueCount; ++test_I ) {

      testValue = testValues[test_I];

      for ( i = 0; i < length; ++i ) {
         Variable_SetValueChar( var, i, testValue );
         Variable_SetValueAtChar( vec, i, 0, testValue );
         Variable_SetValueAtChar( vec, i, 1, testValue );
         Variable_SetValueAtChar( vec, i, 2, testValue );
      }

      /* ~~~Scalar~~~*/
      for ( i = 0; i < length; ++i ) {
         pcu_check_true( Variable_GetValueChar( var, i ) == (char)(char)testValue );
         pcu_check_true( Variable_GetValueCharAsShort( var, i ) == (short)(char)testValue );
         pcu_check_true( Variable_GetValueCharAsInt( var, i ) == (int)(char)testValue );
         pcu_check_true( Variable_GetValueCharAsFloat( var, i ) == (float)(char)testValue );
         pcu_check_true( Variable_GetValueCharAsDouble( var, i ) == (double)(char)testValue );
      }

      /*~~~Vector~~~*/
      for ( i = 0; i < length; ++i ) {
         pcu_check_true( Variable_GetValueAtChar( vec, i, 0 ) == (char)(char)testValue );
         pcu_check_true( Variable_GetValueAtCharAsShort( vec, i, 0 ) == (short)(char)testValue );
         pcu_check_true( Variable_GetValueAtCharAsInt( vec, i, 0 ) == (int)(char)testValue );
         pcu_check_true( Variable_GetValueAtCharAsFloat( vec, i, 0 ) == (float)(char)testValue );
         pcu_check_true( Variable_GetValueAtCharAsDouble( vec, i, 0 ) == (double)(char)testValue );

         pcu_check_true( Variable_GetPtrAtChar( vec, i, 0 ) == &structArray[i][0] );

         pcu_check_true( Variable_GetValueAtChar( vec, i, 1 ) == (char)(char)testValue );
         pcu_check_true( Variable_GetValueAtCharAsShort( vec, i, 1 ) == (short)(char)testValue );
         pcu_check_true( Variable_GetValueAtCharAsInt( vec, i, 1 ) == (int)(char)testValue );
         pcu_check_true( Variable_GetValueAtCharAsFloat( vec, i, 1 ) == (float)(char)testValue );
         pcu_check_true( Variable_GetValueAtCharAsDouble( vec, i, 1 ) == (double)(char)testValue );
         pcu_check_true( Variable_GetPtrAtChar( vec, i, 1 ) == &structArray[i][1] );

         pcu_check_true( Variable_GetValueAtChar( vec, i, 2 ) == (char)(char)testValue );
         pcu_check_true( Variable_GetValueAtCharAsShort( vec, i, 2 ) == (short)(char)testValue );
         pcu_check_true( Variable_GetValueAtCharAsInt( vec, i, 2 ) == (int)(char)testValue );
         pcu_check_true( Variable_GetValueAtCharAsFloat( vec, i, 2 ) == (float)(char)testValue );
         pcu_check_true( Variable_GetValueAtCharAsDouble( vec, i, 2 ) == (double)(char)testValue );
         pcu_check_true( Variable_GetPtrAtChar( vec, i, 2 ) == &structArray[i][2] );
      }

      /*~~~Vector: Sub-Variable~~~*/
      for ( i = 0; i < length; ++i ) {
         for ( j = 0; j < 3; ++j ) {
            pcu_check_true( _Variable_GetStructPtr( vecVar[j], i ) == &structArray[i][j] );
         }
      }
   }
}


void VariableSuite( pcu_suite_t* suite ) {
   pcu_suite_setData( suite, VariableSuiteData );
   pcu_suite_setFixtures( suite, VariableSuite_Setup, VariableSuite_Teardown );
   pcu_suite_addTest( suite, VariableSuite_TestGetValueDouble );
   pcu_suite_addTest( suite, VariableSuite_TestSetValueDouble );
   pcu_suite_addTest( suite, VariableSuite_TestGetValueAtDouble );
   pcu_suite_addTest( suite, VariableSuite_TestSetValueAtDouble );
   pcu_suite_addTest( suite, VariableSuite_TestVariable_Char );
}
