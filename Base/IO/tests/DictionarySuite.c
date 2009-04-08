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
** Role:
**   Tests the data->dict functionality
**
** $Id: testDictionary.c 3462 2006-02-19 06:53:24Z WalterLandry $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "pcu/pcu.h"
#include "StGermain/Base/Foundation/Foundation.h"
#include "StGermain/Base/IO/IO.h"
#include "DictionarySuite.h"

typedef struct {
   double startx;
   double starty;
   double startz;
} Geom;

typedef struct {
   double   height;
   Bool     anisotropic;
   char*    person;
   Geom     geom;
} TestStruct;

typedef struct {
   Dictionary*                dict;
   Index                      testEntriesCount;
   char**                     testKeys;
   Dictionary_Entry_Value**   testValues;
   char*                      testString;
   double                     testDouble;
   unsigned int               testUint;
   int                        testInt;
   unsigned long              testUnsignedlong;
   Bool                       testBool;
   double*                    testList;
   Index                      testListCount;
   TestStruct*                testStruct;
} DictionarySuiteData;


void DictionarySuite_Setup( DictionarySuiteData* data ) {
   Index                   ii=0;
   Index                   iter=0;
   Dictionary_Entry_Value* testStruct; 
   Dictionary_Entry_Value* testStruct2; 
   Dictionary_Entry_Value* testList; 

   data->dict     = Dictionary_New();
   data->testEntriesCount = 8;
   data->testKeys = Memory_Alloc_Array_Unnamed( char*, data->testEntriesCount );
   data->testValues = Memory_Alloc_Array_Unnamed( Dictionary_Entry_Value*,
      data->testEntriesCount );
   for ( ii=0; ii< data->testEntriesCount; ii++ ) {
      data->testKeys[ii] = NULL;
      data->testValues[ii] = NULL;
   }

   Stg_asprintf( &data->testString, "hello" );
   data->testDouble=45.567;
   data->testUint = 5;
   data->testInt = -5;
   data->testUnsignedlong = 52342423;
   data->testBool = True;
   iter = 0;
   Stg_asprintf( &data->testKeys[iter], "test_cstring" );
   data->testValues[iter] = Dictionary_Entry_Value_FromString( data->testString );
   Stg_asprintf( &data->testKeys[++iter], "test_double" );
   data->testValues[iter] = Dictionary_Entry_Value_FromDouble( data->testDouble );
   Stg_asprintf( &data->testKeys[++iter], "test_uint" );
   data->testValues[iter] = Dictionary_Entry_Value_FromUnsignedInt( data->testUint );
   Stg_asprintf( &data->testKeys[++iter], "test_int" );
   data->testValues[iter] = Dictionary_Entry_Value_FromInt( data->testInt );
   Stg_asprintf( &data->testKeys[++iter], "test_unsignedlong" );
   data->testValues[iter] = Dictionary_Entry_Value_FromUnsignedLong( data->testUnsignedlong );
   Stg_asprintf( &data->testKeys[++iter], "test_bool" );
   data->testValues[iter] = Dictionary_Entry_Value_FromUnsignedInt( data->testBool );

   /* adding a list */
   data->testListCount = 5;
   data->testList = Memory_Alloc_Array_Unnamed( double, data->testListCount );
   for (ii=0; ii<data->testListCount; ii++ ) {
      data->testList[ii] = 10.0 * ii;
   }

   testList = Dictionary_Entry_Value_NewList();
   Stg_asprintf( &data->testKeys[++iter], "test_list" );
   data->testValues[iter] = testList;
   for (ii=0; ii < data->testListCount; ii++ ) {
      Dictionary_Entry_Value_AddElement( testList, Dictionary_Entry_Value_FromDouble(data->testList[ii]) );
   }

   /* Adding members to a struct */
   data->testStruct = Memory_Alloc_Unnamed( TestStruct );
   data->testStruct->height = 37;
   data->testStruct->anisotropic = True;
   Stg_asprintf( &data->testStruct->person, "Patrick" );
   data->testStruct->geom.startx = 45;
   data->testStruct->geom.starty = 60;
   data->testStruct->geom.startz = 70;

   testStruct = Dictionary_Entry_Value_NewStruct();
   Stg_asprintf( &data->testKeys[++iter], "test_struct" );
   data->testValues[iter] = testStruct;
   Dictionary_Entry_Value_AddMember( testStruct, "height",
      Dictionary_Entry_Value_FromDouble( 37 ) );
   Dictionary_Entry_Value_AddMember( testStruct, "anisotropic",
      Dictionary_Entry_Value_FromBool( True ) );
   Dictionary_Entry_Value_AddMember( testStruct, "person",
      Dictionary_Entry_Value_FromString( "Patrick" ) );

   /* Adding a 2nd struct within the first struct */
   testStruct2 = Dictionary_Entry_Value_NewStruct();
   Dictionary_Entry_Value_AddMember( testStruct, "geom", testStruct2 );
   Dictionary_Entry_Value_AddMember( testStruct2, "startx",
      Dictionary_Entry_Value_FromUnsignedInt( data->testStruct->geom.startx ) );
   Dictionary_Entry_Value_AddMember( testStruct2, "starty",
      Dictionary_Entry_Value_FromUnsignedInt( data->testStruct->geom.starty ) );
   Dictionary_Entry_Value_AddMember( testStruct2, "startz",
      Dictionary_Entry_Value_FromUnsignedInt( data->testStruct->geom.startz ) );
}

void DictionarySuite_Teardown( DictionarySuiteData* data ) {
   Index ii;

   Stg_Class_Delete( data->dict );
   for ( ii=0; ii< data->testEntriesCount; ii++ ) {
      Memory_Free( data->testKeys[ii] );
      /* Note: we don't free the testValues here, as expect that deleting the
       * dictionary has already done this */
   }
   Memory_Free( data->testKeys );
   Memory_Free( data->testValues );
   Memory_Free( data->testStruct->person );
   Memory_Free( data->testStruct );
   Memory_Free( data->testList );
}


void DictionarySuite_TestCreateValues( DictionarySuiteData* data ) {
   Dictionary_Entry_Value*   dev;

   /* Don't use the pre-created test values. Want to do a very fundamental test here */
   dev = Dictionary_Entry_Value_FromString( "hello" );
   pcu_check_true( Dictionary_Entry_Value_Type_String == dev->type );
   pcu_check_true( 0 == strcmp( "hello", dev->as.typeString ) );
   Dictionary_Entry_Value_Delete( dev );
   dev = Dictionary_Entry_Value_FromDouble( 45.567 );
   pcu_check_true( Dictionary_Entry_Value_Type_Double == dev->type );
   pcu_check_true( 45.567 == dev->as.typeDouble );
   Dictionary_Entry_Value_Delete( dev );
   dev = Dictionary_Entry_Value_FromUnsignedInt( 5 );
   pcu_check_true( Dictionary_Entry_Value_Type_UnsignedInt == dev->type );
   pcu_check_true( 5 == dev->as.typeUnsignedInt );
   Dictionary_Entry_Value_Delete( dev );
   dev = Dictionary_Entry_Value_FromInt( -5 );
   pcu_check_true( Dictionary_Entry_Value_Type_Int == dev->type );
   pcu_check_true( -5 == dev->as.typeInt );
   Dictionary_Entry_Value_Delete( dev );
   dev = Dictionary_Entry_Value_FromUnsignedLong( 52342423 );
   pcu_check_true( Dictionary_Entry_Value_Type_UnsignedLong == dev->type );
   pcu_check_true( 52342423 == dev->as.typeUnsignedLong );
   Dictionary_Entry_Value_Delete( dev );
   dev = Dictionary_Entry_Value_FromBool( True );
   pcu_check_true( Dictionary_Entry_Value_Type_Bool == dev->type );
   pcu_check_true( True == dev->as.typeBool );
   Dictionary_Entry_Value_Delete( dev );

   /* Since we know the DEV Struct is basically a Dictionary, won't test that one
    *  until after we've tested Dictionary_Add works */   
}


DictionarySuite_PopulateDictWithTestValues( DictionarySuiteData* data ) {
   Index ii;

   for ( ii=0; ii< data->testEntriesCount; ii++ ) {
      Dictionary_Add( data->dict, data->testKeys[ii],\
         Dictionary_Entry_Value_Copy( data->testValues[ii], True ) );
   }
}


DictionarySuite_TestCopyCompare( DictionarySuiteData* data ) {
   Index                    ii=0, jj=0;
   Dictionary_Entry_Value*  copiedDev;
   
   for( ii = 0; ii < data->dict->count; ii++ ) {
      copiedDev = Dictionary_Entry_Value_Copy( data->testValues[ii], True );
      pcu_check_true( Dictionary_Entry_Value_Compare( data->testValues[ii],
         copiedDev ) ); 
      Dictionary_Entry_Value_Delete( copiedDev );

      for( jj = 0; ii < data->dict->count; ii++ ) {
         if ( ii == jj ) continue;
         pcu_check_true( False == Dictionary_Entry_Value_Compare( data->testValues[ii],
            data->testValues[jj] ) ); 
      }
   }
}

/* Add a set of values to a dictionary, and check they were added as desired
 * using the Compare functions.
 */
void DictionarySuite_TestAddEmpty( DictionarySuiteData* data ) {
   Dictionary_Index        ii;
   Dictionary_Entry*       entryPtr;

   DictionarySuite_PopulateDictWithTestValues( data );

   pcu_check_true( data->testEntriesCount == data->dict->count );

   for( ii = 0; ii < data->dict->count; ii++ ) {
      entryPtr = data->dict->entryPtr[ii];
      pcu_check_true( Dictionary_Entry_Compare( entryPtr, data->testKeys[ii] ) ); 
      pcu_check_true( Dictionary_Entry_Value_Compare( entryPtr->value, data->testValues[ii] ) ); 
   }

   /* Empty the dict for future tests, but don't call Delete on the test values */
   Dictionary_Empty( data->dict );
   pcu_check_true( 0 == data->dict->count );
}


void DictionarySuite_TestGet( DictionarySuiteData* data ) {
   Dictionary_Entry_Value* yValue;
   Dictionary_Entry_Value* testStruct;

   DictionarySuite_PopulateDictWithTestValues( data );

   testStruct = Dictionary_Get( data->dict, "test_struct" );
   yValue = Dictionary_Entry_Value_GetMember(
      Dictionary_Entry_Value_GetMember(testStruct, "geom"), "starty");
   pcu_check_true( data->testStruct->geom.starty == Dictionary_Entry_Value_AsDouble( yValue ) );

   Dictionary_Empty( data->dict );
}


void DictionarySuite_TestSet( DictionarySuiteData* data ) {
   Dictionary_Entry_Value* currValue;
   Dictionary_Entry_Value* listValue;
   double                  newVal1 = 34.3, newVal2 = 38.9;
   Dictionary_Entry_Value* yValue;
   Dictionary_Entry_Value* testStruct;
   Index                   ii=0;

   DictionarySuite_PopulateDictWithTestValues( data );

   listValue = Dictionary_Get( data->dict, "test_list" );
   /* getting dictionary out of a list */
   currValue = Dictionary_Entry_Value_GetFirstElement( listValue );
   /* do something to this value */
   Dictionary_Entry_Value_SetFromDouble( currValue, newVal1 );
   currValue = currValue->next;
   /* do something to this value */
   Dictionary_Entry_Value_SetFromDouble( currValue, newVal2 );
   
   pcu_check_true( 5 == Dictionary_Entry_Value_GetCount( listValue ) );
   currValue = Dictionary_Entry_Value_GetFirstElement( listValue );
   pcu_check_le( fabs(newVal1 - Dictionary_Entry_Value_AsDouble( currValue )), 0.01 );
   currValue = currValue->next;
   pcu_check_le( fabs(newVal2 - Dictionary_Entry_Value_AsDouble( currValue )), 0.01 );
   currValue = currValue->next;
   for ( ii=2; ii<data->testListCount; ii++ ) {
      pcu_check_le( fabs(data->testList[ii]
         - Dictionary_Entry_Value_AsDouble( currValue )), 0.01 );
      currValue = currValue->next;
   }

   Dictionary_Empty( data->dict );
}


void DictionarySuite_TestAddElement( DictionarySuiteData* data ) {
   Dictionary_Entry_Value* yValue;
   Dictionary_Entry_Value* testStruct;
   Dictionary_Entry_Value* currValue;
   double                  newVal = -45.0;

   DictionarySuite_PopulateDictWithTestValues( data );

   /* turning the starty value into a list using add element */
   testStruct = Dictionary_Get( data->dict, "test_struct" );
   yValue = Dictionary_Entry_Value_GetMember(
      Dictionary_Entry_Value_GetMember(testStruct, "geom"), "starty");
   Dictionary_Entry_Value_AddElement( yValue, Dictionary_Entry_Value_FromDouble(newVal) );

   pcu_check_true( Dictionary_Entry_Value_Type_List == yValue->type );
   pcu_check_true( 2 == Dictionary_Entry_Value_GetCount( yValue ) );

   currValue = Dictionary_Entry_Value_GetFirstElement( yValue );
   pcu_check_le( fabs( data->testStruct->geom.starty - Dictionary_Entry_Value_AsDouble( currValue )), 0.01 );
   currValue = currValue->next;
   pcu_check_le( fabs( newVal - Dictionary_Entry_Value_AsDouble( currValue )), 0.01 );

   Dictionary_Empty( data->dict );
}


void DictionarySuite_TestShortcuts( DictionarySuiteData* data ) {

   DictionarySuite_PopulateDictWithTestValues( data );

   /* Testing GetString_WithDefault. If an entry with that key already exists, then
    *  the value of the existing key should be returned, and the default passed in
    *  ignored. However if the given key _doesn't_ exist, the default should be 
    *  returned, and a new entry with the given key added to the dict. */
   pcu_check_true( 0 == strcmp( data->testString,
      Dictionary_GetString_WithDefault( data->dict, "test_cstring", "heya" ) ) );
   pcu_check_true( 0 == strcmp( "heya",
      Dictionary_GetString_WithDefault( data->dict, "test_cstring2", "heya" ) ) );
   pcu_check_true( NULL != Dictionary_Get( data->dict, "test_cstring2" ) );
   pcu_check_true( data->testDouble =
      Dictionary_GetDouble_WithDefault( data->dict, "test_double", 2.8 ) );
   pcu_check_true( 2.8 ==
      Dictionary_GetDouble_WithDefault( data->dict, "test_double2", 2.8 ) );
   pcu_check_true( NULL != Dictionary_Get( data->dict, "test_double2" ) );
   pcu_check_true( data->testUint =
      Dictionary_GetUnsignedInt_WithDefault( data->dict, "test_uint", 33 ) );
   pcu_check_true( 33 ==
      Dictionary_GetUnsignedInt_WithDefault( data->dict, "test_uint2", 33 ) );
   pcu_check_true( NULL != Dictionary_Get( data->dict, "test_uint2" ) );
   pcu_check_true( data->testInt =
      Dictionary_GetInt_WithDefault( data->dict, "test_int", -24 ) );
   pcu_check_true( -24 ==
      Dictionary_GetInt_WithDefault( data->dict, "test_int2", -24 ) );
   pcu_check_true( NULL != Dictionary_Get( data->dict, "test_int2" ) );
   pcu_check_true( data->testUnsignedlong =
      Dictionary_GetUnsignedLong_WithDefault( data->dict, "test_unsignedlong", 32433 ) );
   pcu_check_true( 32433 ==
      Dictionary_GetUnsignedLong_WithDefault( data->dict, "test_unsignedlong2", 32433 ) );
   pcu_check_true( NULL != Dictionary_Get( data->dict, "test_unsignedlong2" ) );
   pcu_check_true( data->testBool =
      Dictionary_GetBool_WithDefault( data->dict, "test_bool", False ) );
   pcu_check_true( False ==
      Dictionary_GetBool_WithDefault( data->dict, "test_bool2", False ) );
   pcu_check_true( NULL != Dictionary_Get( data->dict, "test_bool2" ) );

   pcu_check_true( 0 == strcmp( data->testString,
      Dictionary_GetString_WithPrintfDefault( data->dict, "test_cstring",
         "heya%s%u", "hey", 3 ) ) );
   pcu_check_true( 0 == strcmp( "heyahey3",
      Dictionary_GetString_WithPrintfDefault( data->dict, "test_cstring3",
         "heya%s%u", "hey", 3 ) ) );
   pcu_check_true( NULL != Dictionary_Get( data->dict, "test_cstring3" ) );

   Dictionary_Empty( data->dict );
}


void DictionarySuite( pcu_suite_t* suite ) {
   pcu_suite_setData( suite, DictionarySuiteData );
   pcu_suite_setFixtures( suite, DictionarySuite_Setup, DictionarySuite_Teardown );
   pcu_suite_addTest( suite, DictionarySuite_TestCreateValues );
   pcu_suite_addTest( suite, DictionarySuite_TestCopyCompare );
   pcu_suite_addTest( suite, DictionarySuite_TestAddEmpty );
   pcu_suite_addTest( suite, DictionarySuite_TestGet );
   pcu_suite_addTest( suite, DictionarySuite_TestSet );
   pcu_suite_addTest( suite, DictionarySuite_TestAddElement );
   pcu_suite_addTest( suite, DictionarySuite_TestShortcuts );
}
