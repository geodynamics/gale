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
**   Tests the IO handler
**
** $Id: testIO_Handler-normal.c 3743 2006-08-03 03:14:38Z KentHumphries $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "pcu/pcu.h"
#include "StGermain/Base/Foundation/Foundation.h"
#include "StGermain/Base/IO/IO.h"
#include "DictionarySuite.h"  /* Because want to re-use sample dictionary structs/funcs */
#include "IO_HandlerSuite.h"


const char* IO_HandlerSuite_XMLStartString1 = "<?xml version=\"1.0\"?>\n";
const char* IO_HandlerSuite_XMLStartString2 = "<StGermainData xmlns=\"http://www.vpac.org/StGermain/XML_IO_Handler/Jun2003\">\n";
const char* IO_HandlerSuite_XMLEndString = "</StGermainData>\n";
const char* IO_HandlerSuite_XMLEmptyDataString = "<StGermainData xmlns=\"http://www.vpac.org/StGermain/XML_IO_Handler/Jun2003\"/>\n";

typedef struct {
   XML_IO_Handler*                io_handler;
   Dictionary*                    dict1;
   Dictionary*                    dict2;
   DictionarySuite_TestDictData*  testDD;
} IO_HandlerSuiteData;


void _IO_HandlerSuite_CreateTestXMLFile( const char* testXMLFileName,
     const char* entriesString )
{
   FILE*         testFile = NULL;
   testFile = fopen(testXMLFileName, "w");
   fwrite( IO_HandlerSuite_XMLStartString1, sizeof(char),
      strlen( IO_HandlerSuite_XMLStartString1 ), testFile );
   fwrite( IO_HandlerSuite_XMLStartString2, sizeof(char),
      strlen( IO_HandlerSuite_XMLStartString2 ), testFile );
   fwrite( entriesString, sizeof(char), strlen( entriesString ), testFile );
   fwrite( IO_HandlerSuite_XMLEndString, sizeof(char),
      strlen( IO_HandlerSuite_XMLEndString ), testFile );
   fclose( testFile );
}


void IO_HandlerSuite_Setup( IO_HandlerSuiteData* data ) {
   data->io_handler = XML_IO_Handler_New();
   /* We don't want output in the tests by default */
   Stream_Enable( Journal_Register( Debug_Type, XML_IO_Handler_Type ), False );
   Stream_Enable( Journal_Register( Info_Type, XML_IO_Handler_Type ), False );
   data->dict1 = Dictionary_New();
   data->dict2 = Dictionary_New();
   data->testDD   = Memory_Alloc_Unnamed( DictionarySuite_TestDictData );
   DictionarySuite_SetupTestDictData( data->testDD );
}


void IO_HandlerSuite_Teardown( IO_HandlerSuiteData* data ) {
   Stg_Class_Delete( data->io_handler );
   Stg_Class_Delete( data->dict1 );
   Stg_Class_Delete( data->dict2 );
   DictionarySuite_DictionaryData_Free( data->testDD );
   Memory_Free( data->testDD );

}


/* Just populate a test dictionary, write it out to a file, read it back in again to a different dict, and check all the values are the same */
void IO_HandlerSuite_TestWriteReadNormalEntries( IO_HandlerSuiteData* data ) {
   Index         ii;
   const char*   xmlTestFileName = "xmlTest.xml";

   Dictionary_Empty( data->dict1 );
   Dictionary_Empty( data->dict2 );
   DictionarySuite_PopulateDictWithTestValues( data->dict1, data->testDD );

   IO_Handler_WriteAllToFile( data->io_handler, xmlTestFileName, data->dict1 );
   IO_Handler_ReadAllFromFile( data->io_handler, xmlTestFileName, data->dict2 ); 

   pcu_check_true( data->dict1->count == data->dict2->count );
   for (ii=0; ii<data->dict1->count; ii++) {
      pcu_check_true( Dictionary_Entry_Compare( data->dict1->entryPtr[ii],
         data->dict2->entryPtr[ii]->key) );
      pcu_check_true( Dictionary_Entry_Value_Compare( data->dict1->entryPtr[ii]->value,
         data->dict2->entryPtr[ii]->value) );
   }

   remove(xmlTestFileName);
   Dictionary_Empty( data->dict1 );
   Dictionary_Empty( data->dict2 );
}


/* Similar to above test, except using the function to write just one entry at a time */
void IO_HandlerSuite_TestWriteReadNormalSingleEntry( IO_HandlerSuiteData* data ) {
   Index         ii;
   const char*   fileName = "singleEntry.xml";

   Dictionary_Empty( data->dict1 );
   Dictionary_Empty( data->dict2 );
   DictionarySuite_PopulateDictWithTestValues( data->dict1, data->testDD );

   for (ii=0; ii<data->dict1->count; ii++) {
      XML_IO_Handler_WriteEntryToFile( data->io_handler, fileName,
            data->testDD->testKeys[ii],
            data->testDD->testValues[ii], 
            NULL );
      IO_Handler_ReadAllFromFile( data->io_handler, fileName, data->dict2 ); 

      pcu_check_true( 1 == data->dict2->count );
      pcu_check_true( Dictionary_Entry_Compare( data->dict2->entryPtr[0],
         data->testDD->testKeys[ii]) );
      pcu_check_true( Dictionary_Entry_Value_Compare( data->dict2->entryPtr[0]->value,
         data->testDD->testValues[ii] ) );

      Dictionary_Empty( data->dict2 );
      remove(fileName);
   }

   Dictionary_Empty( data->dict1 );
   Dictionary_Empty( data->dict2 );
}


/* Similar to above test, except test we can write out an empty Dictionary, then read in */
void IO_HandlerSuite_TestWriteReadEmpty( IO_HandlerSuiteData* data ) {
   Index         ii;
   const char*   xmlTestFileName = "empty.xml";
   FILE*         testFile = NULL;
   const int     MAXLINE = 1000;
   char*         xmlLine = NULL;

   Dictionary_Empty( data->dict1 );
   Dictionary_Empty( data->dict2 );

   IO_Handler_WriteAllToFile( data->io_handler, xmlTestFileName, data->dict1 );

   testFile = fopen(xmlTestFileName, "r");
   xmlLine = Memory_Alloc_Array_Unnamed( char, MAXLINE );
   pcu_check_true( fgets( xmlLine, sizeof(char*)*MAXLINE, testFile ) );
   pcu_check_true( 0 == strcmp( IO_HandlerSuite_XMLStartString1, xmlLine ) );
   pcu_check_true( fgets( xmlLine, sizeof(char*)*MAXLINE, testFile ) );
   pcu_check_true( 0 == strcmp( IO_HandlerSuite_XMLEmptyDataString, xmlLine ) );
   Memory_Free( xmlLine );
   fclose(testFile);

   IO_Handler_ReadAllFromFile( data->io_handler, xmlTestFileName, data->dict2 ); 
   pcu_check_true( 0 == data->dict2->count );

   remove(xmlTestFileName);
   Dictionary_Empty( data->dict1 );
   Dictionary_Empty( data->dict2 );
}


/* In this case, want to make sure the types are written explicitly into the output, so will
 * check against expected text. */
void IO_HandlerSuite_TestWriteExplicitTypes( IO_HandlerSuiteData* data ) {
   Index         ii=0;
   const char*   testFileName = "xmlTest-explicittypes.xml";
   const int     MAXLINE = 1000;
   FILE*         testFile = NULL;
   char*         xmlLine = NULL;
   const unsigned explicityTypesExpectedLineNum = 23;
   const char* explicitTypesExpected[] = {
      "  <element type=\"param\" name=\"test_cstring\" paramType=\"string\">hello</element>\n",
      "  <element type=\"param\" name=\"test_double\" paramType=\"double\">45.567</element>\n",
      "  <element type=\"param\" name=\"test_uint\" paramType=\"uint\">5</element>\n",
      "  <element type=\"param\" name=\"test_int\" paramType=\"int\">-5</element>\n",
      "  <element type=\"param\" name=\"test_unsignedlong\" paramType=\"ulong\">52342423</element>\n",
      "  <element type=\"param\" name=\"test_bool\" paramType=\"uint\">1</element>\n",
      "  <element type=\"list\" name=\"test_list\">\n",
      "    <element type=\"param\" paramType=\"double\">0</element>\n",
      "    <element type=\"param\" paramType=\"double\">10</element>\n",
      "    <element type=\"param\" paramType=\"double\">20</element>\n",
      "    <element type=\"param\" paramType=\"double\">30</element>\n",
      "    <element type=\"param\" paramType=\"double\">40</element>\n",
      "  </element>\n",
      "  <element type=\"struct\" name=\"test_struct\">\n",
      "    <element type=\"param\" name=\"height\" paramType=\"double\">37</element>\n",
      "    <element type=\"param\" name=\"anisotropic\" paramType=\"bool\">true</element>\n",
      "    <element type=\"param\" name=\"person\" paramType=\"string\">Patrick</element>\n",
      "    <element type=\"struct\" name=\"geom\">\n",
      "      <element type=\"param\" name=\"startx\" paramType=\"uint\">45</element>\n",
      "      <element type=\"param\" name=\"starty\" paramType=\"uint\">60</element>\n",
      "      <element type=\"param\" name=\"startz\" paramType=\"uint\">70</element>\n",
      "    </element>\n",
      "  </element>\n"};

   xmlLine = Memory_Alloc_Array_Unnamed( char, MAXLINE );

   Dictionary_Empty( data->dict1 );
   DictionarySuite_PopulateDictWithTestValues( data->dict1, data->testDD );

   XML_IO_Handler_SetWriteExplicitTypes( data->io_handler, True );
   IO_Handler_WriteAllToFile( data->io_handler, testFileName, data->dict1 );

   testFile = fopen(testFileName, "r");
   pcu_check_true( fgets( xmlLine, sizeof(char*)*MAXLINE, testFile ) );
   pcu_check_true( 0 == strcmp( IO_HandlerSuite_XMLStartString1, xmlLine ) );
   pcu_check_true( fgets( xmlLine, sizeof(char*)*MAXLINE, testFile ) );
   pcu_check_true( 0 == strcmp( IO_HandlerSuite_XMLStartString2, xmlLine ) );
   for ( ii=0; ii< explicityTypesExpectedLineNum; ii++ ) {
      pcu_check_true( fgets( xmlLine, sizeof(char*)*MAXLINE, testFile ) );
      pcu_check_true( 0 == strcmp( explicitTypesExpected[ii], xmlLine ) );
   }
   pcu_check_true( fgets( xmlLine, sizeof(char*)*MAXLINE, testFile ) );
   pcu_check_true( 0 == strcmp( IO_HandlerSuite_XMLEndString, xmlLine ) );
   fclose(testFile);

   remove(testFileName);
   Memory_Free( xmlLine );
   Dictionary_Empty( data->dict1 );
}


void IO_HandlerSuite_TestReadWhitespaceEntries( IO_HandlerSuiteData* data ) {
   Index             ii;
   const char*       testFileName = "xmlTest-whitespaces.xml";
   char*             whiteSpacesEntry = NULL;
   const char*       testKey = "spacedKey";
   const char*       testValString = "spacedVal";

   Dictionary_Empty( data->dict2 );

   Stg_asprintf( &whiteSpacesEntry, "<param name=\"    %s   \"> \t %s \n\t</param>\n",
      testKey, testValString );
   _IO_HandlerSuite_CreateTestXMLFile( testFileName, whiteSpacesEntry );
   Memory_Free( whiteSpacesEntry );

   IO_Handler_ReadAllFromFile( data->io_handler, testFileName, data->dict2 ); 

   pcu_check_true( 1 == data->dict2->count );
   pcu_check_true( Dictionary_Entry_Compare( data->dict2->entryPtr[0],
      (Dictionary_Entry_Key)testKey) );
   pcu_check_true( 0 == strcmp(
      Dictionary_Entry_Value_AsString( data->dict2->entryPtr[0]->value ), testValString ) );

   remove( testFileName );
   Dictionary_Empty( data->dict2 );
}


/* Testing the functionality of using included files. Including specifying a search path */
void IO_HandlerSuite_TestReadIncludedFile( IO_HandlerSuiteData* data ) {
   Index             ii;
   const char*       testFileName = "xmlTest-include.xml";
   const char*       testIncludedFileName = "xmlTest-included.xml";
   const char*       testSearchPathSubdir = "./testXML-subdir";
   const char*       testIncludedFileNameSP = "xmlTest-includedSP.xml";
   char*             subdirIncludedFileNameSP = NULL;
   char*             xmlEntry = NULL;
   char*             includeLine = NULL;
   char*             searchPathLine = NULL;
   char*             includeLineSP = NULL;
   char*             xmlTestEntries = NULL;
   const char*       testKey = "regularKey";
   const char*       testValString = "regularVal";
   const char*       testKeyInc = "keyInc";
   const char*       testValStringInc = "valInc";
   const char*       testKeyIncSP = "keyIncSP";
   const char*       testValStringIncSP = "valIncSP";

   Dictionary_Empty( data->dict2 );

   Stg_asprintf( &subdirIncludedFileNameSP, "%s/%s", testSearchPathSubdir, testIncludedFileNameSP );

   Stg_asprintf( &xmlEntry, "<param name=\"%s\">%s</param>\n",
      testKey, testValString );
   Stg_asprintf( &includeLine, "<include>%s</include>\n", testIncludedFileName );
   Stg_asprintf( &searchPathLine, "<searchPath>%s</searchPath>\n", testSearchPathSubdir );
   Stg_asprintf( &includeLineSP, "<include>%s</include>\n", testIncludedFileNameSP );
   Stg_asprintf( &xmlTestEntries, "%s%s%s%s", xmlEntry, includeLine, searchPathLine,
      includeLineSP );
   _IO_HandlerSuite_CreateTestXMLFile( testFileName, xmlTestEntries );
   Memory_Free( xmlEntry );
   Memory_Free( includeLine );
   Memory_Free( searchPathLine );
   Memory_Free( includeLineSP );
   Memory_Free( xmlTestEntries );

   Stg_asprintf( &xmlEntry, "<param name=\"%s\">%s</param>\n",
      testKeyInc, testValStringInc );
   _IO_HandlerSuite_CreateTestXMLFile( testIncludedFileName, xmlEntry );
   Memory_Free( xmlEntry );

   mkdir( testSearchPathSubdir, 0755 );
   Stg_asprintf( &xmlEntry, "<param name=\"%s\">%s</param>\n",
      testKeyIncSP, testValStringIncSP );
   _IO_HandlerSuite_CreateTestXMLFile( subdirIncludedFileNameSP, xmlEntry );
   Memory_Free( xmlEntry );

   IO_Handler_ReadAllFromFile( data->io_handler, testFileName, data->dict2 ); 

   pcu_check_true( 3 == data->dict2->count );
   pcu_check_true( Dictionary_Entry_Compare( data->dict2->entryPtr[0],
      (Dictionary_Entry_Key)testKey) );
   pcu_check_true( 0 == strcmp(
      Dictionary_Entry_Value_AsString( data->dict2->entryPtr[0]->value ), testValString ) );
   pcu_check_true( Dictionary_Entry_Compare( data->dict2->entryPtr[1],
      (Dictionary_Entry_Key)testKeyInc) );
   pcu_check_true( 0 == strcmp(
      Dictionary_Entry_Value_AsString( data->dict2->entryPtr[1]->value ), testValStringInc ) );
   pcu_check_true( Dictionary_Entry_Compare( data->dict2->entryPtr[2],
      (Dictionary_Entry_Key)testKeyIncSP) );
   pcu_check_true( 0 == strcmp(
      Dictionary_Entry_Value_AsString( data->dict2->entryPtr[2]->value ), testValStringIncSP ) );

   remove( testFileName );
   remove( testIncludedFileName );
   remove( subdirIncludedFileNameSP );
   Memory_Free( subdirIncludedFileNameSP );
   rmdir( testSearchPathSubdir );
   Dictionary_Empty( data->dict2 );
}


void IO_HandlerSuite_TestReadRawDataEntries( IO_HandlerSuiteData* data ) {
   Index             ii;
   const char*       testFileName = "xmlTest-rawData.xml";
   char*             xmlEntries = NULL;
   char*             rawDataEntry1 = NULL;
   char*             rawDataEntry2 = NULL;
   char*             entryLine = NULL;
   const char*       list1Name = "bcs";
   const int         list1EntryCount = 2;
   const int         list1Vals[2][3] = { {1, 3, 6}, {2, 9, 14} };
   const char*       list2Name = "boundary_conditions2";
   const int         list2CompCount = 5;
   const int         list2EntryCount = 3;
   const char*       list2CompNames[5] = {"side", "xval", "yval", "zval", "active"};
   const char*       list2CompTypes[5] = {"string", "int", "int", "int", "bool"};
   const char*       list2StringVals[3] = {"top", "bottom", "left"};
   const int         list2CoordVals[3][3] = { {4,5,8}, {3,5,9}, {9,3,4} };
   const Bool        list2BoolVals[3] = { True, False, True };
   const char*       list2BoolValStrings[3] = { "True", "False", "1" };

   Dictionary_Empty( data->dict2 );

   Stg_asprintf( &rawDataEntry1, "<list name=\"%s\">\n<asciidata>\n%d %d %d\n%d %d %d\n"
      "</asciidata>\n</list>\n",
      list1Name, list1Vals[0][0], list1Vals[0][1], list1Vals[0][2], 
      list1Vals[1][0], list1Vals[1][1], list1Vals[1][2] );

   rawDataEntry2 = Memory_Alloc_Array_Unnamed( char, 10000 );
   entryLine = Memory_Alloc_Array_Unnamed( char, 1000 );
   sprintf( rawDataEntry2, "<list name=\"%s\">\n<asciidata>\n", list2Name );
   for (ii=0; ii < list2CompCount; ii++ ) {
      sprintf( entryLine, "<columnDefinition name=\"%s\" type=\"%s\"/>\n",
         list2CompNames[ii], list2CompTypes[ii] );
      strcat( rawDataEntry2, entryLine );
   }
   for (ii=0; ii < list2EntryCount; ii++ ) {
      sprintf( entryLine, "%s %i %i %i %s\n", list2StringVals[ii],
         list2CoordVals[ii][0], list2CoordVals[ii][1], list2CoordVals[ii][2],
         list2BoolValStrings[ii] );
      strcat( rawDataEntry2, entryLine );
   }
   sprintf( entryLine, "</asciidata>\n</list>\n" );
   strcat( rawDataEntry2, entryLine );

   Stg_asprintf( &xmlEntries, "%s%s", rawDataEntry1, rawDataEntry2 );
   _IO_HandlerSuite_CreateTestXMLFile( testFileName, xmlEntries );
   Memory_Free( xmlEntries );
   Memory_Free( rawDataEntry1 );
   Memory_Free( rawDataEntry2 );
   Memory_Free( entryLine );

   IO_Handler_ReadAllFromFile( data->io_handler, testFileName, data->dict2 ); 

   {
      Dictionary_Entry_Value* dev = NULL;
      int                     intVal = 0;
      char*                   strVal = 0;
      Bool                    boolVal = False;

      pcu_check_true( 2 == data->dict2->count );
      pcu_check_true( Dictionary_Entry_Compare( data->dict2->entryPtr[0],
         (Dictionary_Entry_Key)list1Name) );
      pcu_check_true( Dictionary_Entry_Value_Type_List ==
         data->dict2->entryPtr[0]->value->type );
      for (ii=0; ii < list1EntryCount; ii++ ) {
         dev = Dictionary_Entry_Value_GetElement( data->dict2->entryPtr[0]->value, ii );
         intVal = Dictionary_Entry_Value_AsInt( Dictionary_Entry_Value_GetMember( dev, "0" ) );
         pcu_check_true( intVal == list1Vals[ii][0] );
         intVal = Dictionary_Entry_Value_AsInt( Dictionary_Entry_Value_GetMember( dev, "1" ) );
         pcu_check_true( intVal == list1Vals[ii][1] );
         intVal = Dictionary_Entry_Value_AsInt( Dictionary_Entry_Value_GetMember( dev, "2" ) );
         pcu_check_true( intVal == list1Vals[ii][2] );
      }
      pcu_check_true( Dictionary_Entry_Compare( data->dict2->entryPtr[1],
         (Dictionary_Entry_Key)list2Name) );
      pcu_check_true( Dictionary_Entry_Value_Type_List ==
         data->dict2->entryPtr[1]->value->type );
      for (ii=0; ii < list2EntryCount; ii++ ) {
         dev = Dictionary_Entry_Value_GetElement( data->dict2->entryPtr[1]->value, ii );
         strVal = Dictionary_Entry_Value_AsString( Dictionary_Entry_Value_GetMember(
            dev, (Dictionary_Entry_Key)list2CompNames[0] ) );
         pcu_check_true( 0 == strcmp( list2StringVals[ii], strVal ) );
         intVal = Dictionary_Entry_Value_AsInt( Dictionary_Entry_Value_GetMember(
            dev, (Dictionary_Entry_Key)list2CompNames[1] ) );
         pcu_check_true( intVal == list2CoordVals[ii][0] );
         intVal = Dictionary_Entry_Value_AsInt( Dictionary_Entry_Value_GetMember(
            dev, (Dictionary_Entry_Key)list2CompNames[2] ) );
         pcu_check_true( intVal == list2CoordVals[ii][1] );
         intVal = Dictionary_Entry_Value_AsInt( Dictionary_Entry_Value_GetMember(
            dev, (Dictionary_Entry_Key)list2CompNames[3] ) );
         pcu_check_true( intVal == list2CoordVals[ii][2] );
         boolVal = Dictionary_Entry_Value_AsBool( Dictionary_Entry_Value_GetMember(
            dev, (Dictionary_Entry_Key)list2CompNames[4] ) );
         pcu_check_true( boolVal == list2BoolVals[ii] );
      }
   }

   remove( testFileName );
   Dictionary_Empty( data->dict2 );
}


// TODO read: mock cmd line, with multiple files
// TODO read: various bad entries

void IO_HandlerSuite( pcu_suite_t* suite ) {
   pcu_suite_setData( suite, IO_HandlerSuiteData );
   pcu_suite_setFixtures( suite, IO_HandlerSuite_Setup, IO_HandlerSuite_Teardown );
   pcu_suite_addTest( suite, IO_HandlerSuite_TestWriteReadNormalEntries );
   pcu_suite_addTest( suite, IO_HandlerSuite_TestWriteReadNormalSingleEntry );
   pcu_suite_addTest( suite, IO_HandlerSuite_TestWriteReadEmpty );
   pcu_suite_addTest( suite, IO_HandlerSuite_TestWriteExplicitTypes );
   pcu_suite_addTest( suite, IO_HandlerSuite_TestReadWhitespaceEntries );
   pcu_suite_addTest( suite, IO_HandlerSuite_TestReadIncludedFile );
   pcu_suite_addTest( suite, IO_HandlerSuite_TestReadRawDataEntries );
}
