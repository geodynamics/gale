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
   if ( data->dict1->count == data->dict2->count ) {
      for (ii=0; ii<data->dict1->count; ii++) {
         pcu_check_true( Dictionary_Entry_Compare( data->dict1->entryPtr[ii],
            data->dict2->entryPtr[ii]->key) );
         pcu_check_true( Dictionary_Entry_Value_Compare( data->dict1->entryPtr[ii]->value,
            data->dict2->entryPtr[ii]->value) );
      }
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
      if ( 1 == data->dict2->count ) {
         pcu_check_true( Dictionary_Entry_Compare( data->dict2->entryPtr[0],
            data->testDD->testKeys[ii]) );
         pcu_check_true( Dictionary_Entry_Value_Compare( data->dict2->entryPtr[0]->value,
            data->testDD->testValues[ii] ) );
      }

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
   pcu_check_true( fgets( xmlLine, MAXLINE, testFile ) );
   pcu_check_true( 0 == strcmp( IO_HandlerSuite_XMLStartString1, xmlLine ) );
   pcu_check_true( fgets( xmlLine, MAXLINE, testFile ) );
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
   pcu_check_true( fgets( xmlLine, MAXLINE, testFile ) );
   pcu_check_true( 0 == strcmp( IO_HandlerSuite_XMLStartString1, xmlLine ) );
   pcu_check_true( fgets( xmlLine, MAXLINE, testFile ) );
   pcu_check_true( 0 == strcmp( IO_HandlerSuite_XMLStartString2, xmlLine ) );
   for ( ii=0; ii< explicityTypesExpectedLineNum; ii++ ) {
      pcu_check_true( fgets( xmlLine, MAXLINE, testFile ) );
      pcu_check_true( 0 == strcmp( explicitTypesExpected[ii], xmlLine ) );
   }
   pcu_check_true( fgets( xmlLine, MAXLINE, testFile ) );
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
   if ( 1 == data->dict2->count ) {
      pcu_check_true( Dictionary_Entry_Compare( data->dict2->entryPtr[0],
         (Dictionary_Entry_Key)testKey) );
      pcu_check_true( 0 == strcmp(
         Dictionary_Entry_Value_AsString( data->dict2->entryPtr[0]->value ), testValString ) );
   }

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
   if ( 3 == data->dict2->count ) {
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
   }

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

void IO_HandlerSuite_TestReadAllFromCommandLine( IO_HandlerSuiteData* data ) {
   Index          ii;
   char**         xmlTestFileNames;
   int            argc;
   char**         argv;
   int            fakeParamArgsCount = 2;
   
   Dictionary_Empty( data->dict1 );
   Dictionary_Empty( data->dict2 );

   DictionarySuite_PopulateDictWithTestValues( data->dict1, data->testDD );

   xmlTestFileNames = Memory_Alloc_Array_Unnamed( char*, data->testDD->testEntriesCount );
   argc = data->testDD->testEntriesCount + 1 + fakeParamArgsCount;
   argv = Memory_Alloc_Array_Unnamed( char*, argc );

   for ( ii=0; ii < data->testDD->testEntriesCount; ii++ ) {
      Stg_asprintf( &xmlTestFileNames[ii], "readFromCommandLineTest%u.xml", ii );
      XML_IO_Handler_WriteEntryToFile( data->io_handler, xmlTestFileNames[ii],
            data->testDD->testKeys[ii],
            data->testDD->testValues[ii], 
            NULL );
   }

   /* Create the argv command line */
   Stg_asprintf( &argv[0], "./testStGermain");
   for ( ii=0; ii < data->testDD->testEntriesCount; ii++ ) {
      Stg_asprintf( &argv[1+ii], "%s", xmlTestFileNames[ii] );
   }
   /* Now just add a couple of extra cmd line entries, to simulate user passing other 
    *  parameters, which should be ignored by the XML IO Handler */
   for ( ii=0; ii < fakeParamArgsCount; ii++ ) {
      Stg_asprintf( &argv[1+data->testDD->testEntriesCount+ii], "simParam%u=test", ii );
   }

   IO_Handler_ReadAllFilesFromCommandLine( data->io_handler, argc, argv, data->dict2 );

   /* Now, dict2 should correspond to dict1, having read in and combined all the
    *  separate files. */
   pcu_check_true( data->dict1->count == data->dict2->count );
   for (ii=0; ii<data->dict1->count; ii++) {
      pcu_check_true( Dictionary_Entry_Compare( data->dict1->entryPtr[ii],
         data->dict2->entryPtr[ii]->key) );
      pcu_check_true( Dictionary_Entry_Value_Compare( data->dict1->entryPtr[ii]->value,
         data->dict2->entryPtr[ii]->value) );
   }


   for ( ii=0; ii < data->testDD->testEntriesCount; ii++ ) {
      remove(xmlTestFileNames[ii]);
      Memory_Free( xmlTestFileNames[ii] );
   }
   Memory_Free( xmlTestFileNames );
   for ( ii=0; ii < argc; ii++ ) {
      Memory_Free( argv[ii] );
   }
   Memory_Free( argv );
   Dictionary_Empty( data->dict1 );
   Dictionary_Empty( data->dict2 );
}


/* It's only worthwhile to test the different mergeType operations for one type (eg struct)
 * in this test - the thorough testing of all the merge operations themselves should be 
 * done in DictionarySuite.c */
void IO_HandlerSuite_TestReadDuplicateEntryKeys( IO_HandlerSuiteData* data ) {
   Index                   ii=0;
   const char*             xmlTestFileName = "testXML-dupKeys.xml";
   const char*             struct1Name = "structOne";
   const int               struct1_OrigParamCount = 2;
   const char*             paramNames[2] = { "paramOne", "paramTwo" };
   const char*             paramNames2[2] = { "2nd-paramOne", "2nd-paramTwo" };
   const unsigned int      paramVals[2] = { 1, 2 };
   const unsigned int      paramVals2[2] = { 3, 4 };
   char                    struct1Entry[10000];
   char                    struct2Entry[10000];
   char                    xmlLine[1000];
   char*                   testEntries = NULL;
   Dictionary_Entry_Value* structDev = NULL;
   Dictionary_Entry_Value* elementDev = NULL;
   Dictionary*             structDict = NULL;

   Dictionary_Empty( data->dict2 );

   sprintf( xmlLine, "<struct name=\"%s\">\n", struct1Name );
   strcat( struct1Entry, xmlLine );
   sprintf( xmlLine, "<param name=\"%s\">%u</param>\n", paramNames[0], paramVals[0] );
   strcat( struct1Entry, xmlLine );
   sprintf( xmlLine, "<param name=\"%s\">%u</param>\n", paramNames[1], paramVals[1] );
   strcat( struct1Entry, xmlLine );
   sprintf( xmlLine, "</struct>\n" );
   strcat( struct1Entry, xmlLine );

   /* Sub-test 1: we expect default behaviour is "replace", therefore the 2nd struct
    *  should be the only entry */
   pcu_check_true( IO_Handler_DefaultMergeType == Dictionary_MergeType_Replace );

   sprintf( xmlLine, "<struct name=\"%s\">\n", struct1Name );
   strcat( struct2Entry, xmlLine );
   sprintf( xmlLine, "<param name=\"%s\">%u</param>\n", paramNames2[0], paramVals2[0] );
   strcat( struct2Entry, xmlLine );
   sprintf( xmlLine, "<param name=\"%s\">%u</param>\n", paramNames2[1], paramVals2[1] );
   strcat( struct2Entry, xmlLine );
   sprintf( xmlLine, "</struct>\n" );
   strcat( struct2Entry, xmlLine );

   Stg_asprintf( &testEntries, "%s%s", struct1Entry, struct2Entry );
   _IO_HandlerSuite_CreateTestXMLFile( xmlTestFileName, testEntries );
   Memory_Free( testEntries );
   IO_Handler_ReadAllFromFile( data->io_handler, xmlTestFileName, data->dict2 );

   pcu_check_true( 1 == data->dict2->count );
   pcu_check_true( Dictionary_Entry_Compare( data->dict2->entryPtr[0],
      (Dictionary_Entry_Key)struct1Name) );
   structDev = data->dict2->entryPtr[0]->value;
   pcu_check_true( Dictionary_Entry_Value_Type_Struct == structDev->type );
   pcu_check_true( struct1_OrigParamCount == Dictionary_Entry_Value_GetCount( structDev ) );
   for (ii=0; ii < struct1_OrigParamCount; ii++ ) {
      elementDev = Dictionary_Entry_Value_GetMember( structDev,
         (Dictionary_Entry_Key)paramNames2[ii] );
      pcu_check_true( paramVals2[ii] == Dictionary_Entry_Value_AsUnsignedInt( elementDev ) );
   }
   remove( xmlTestFileName );
   Dictionary_Empty( data->dict2 );

   /* Sub-test 2: with mergeType as "append", the 2 structs should be 2 separate entries */
   sprintf( struct2Entry, "" );
   sprintf( xmlLine, "<struct name=\"%s\" mergeType=\"append\">\n", struct1Name );
   strcat( struct2Entry, xmlLine );
   sprintf( xmlLine, "<param name=\"%s\">%u</param>\n", paramNames[0], paramVals2[0] );
   strcat( struct2Entry, xmlLine );
   sprintf( xmlLine, "<param name=\"%s\">%u</param>\n", paramNames[1], paramVals2[1] );
   strcat( struct2Entry, xmlLine );
   sprintf( xmlLine, "</struct>\n" );
   strcat( struct2Entry, xmlLine );

   Stg_asprintf( &testEntries, "%s%s", struct1Entry, struct2Entry );
   _IO_HandlerSuite_CreateTestXMLFile( xmlTestFileName, testEntries );
   Memory_Free( testEntries );
   IO_Handler_ReadAllFromFile( data->io_handler, xmlTestFileName, data->dict2 );

   pcu_check_true( 2 == data->dict2->count );
   /* First entry should be unchanged */
   pcu_check_true( Dictionary_Entry_Compare( data->dict2->entryPtr[0],
      (Dictionary_Entry_Key)struct1Name) );
   structDev = data->dict2->entryPtr[0]->value;
   pcu_check_true( Dictionary_Entry_Value_Type_Struct == structDev->type );
   pcu_check_true( struct1_OrigParamCount == Dictionary_Entry_Value_GetCount( structDev ) );
   for (ii=0; ii < struct1_OrigParamCount; ii++ ) {
      elementDev = Dictionary_Entry_Value_GetMember( structDev,
         (Dictionary_Entry_Key)paramNames[ii] );
      pcu_check_true( paramVals[ii] == Dictionary_Entry_Value_AsUnsignedInt( elementDev ) );
   }
   /* Second entry should be struct2 */
   pcu_check_true( Dictionary_Entry_Compare( data->dict2->entryPtr[1],
      (Dictionary_Entry_Key)struct1Name) );
   structDev = data->dict2->entryPtr[1]->value;
   pcu_check_true( Dictionary_Entry_Value_Type_Struct == structDev->type );
   pcu_check_true( struct1_OrigParamCount == Dictionary_Entry_Value_GetCount( structDev ) );
   for (ii=0; ii < struct1_OrigParamCount; ii++ ) {
      elementDev = Dictionary_Entry_Value_GetMember( structDev,
         (Dictionary_Entry_Key)paramNames[ii] );
      pcu_check_true( paramVals2[ii] == Dictionary_Entry_Value_AsUnsignedInt( elementDev ) );
   }
   //remove( xmlTestFileName );
   Dictionary_Empty( data->dict2 );

   /* Sub-test 3.1: with mergeType as "merge", structs to be merged.
    * However, default childrenMergeType is "append", so all entries added */
   sprintf( struct2Entry, "" );
   sprintf( xmlLine, "<struct name=\"%s\" mergeType=\"merge\">\n", struct1Name );
   strcat( struct2Entry, xmlLine );
   sprintf( xmlLine, "<param name=\"%s\">%u</param>\n", paramNames[1], paramVals2[1] );
   strcat( struct2Entry, xmlLine );
   sprintf( xmlLine, "<param name=\"%s\">%u</param>\n", paramNames2[0], paramVals2[0] );
   strcat( struct2Entry, xmlLine );
   sprintf( xmlLine, "</struct>\n" );
   strcat( struct2Entry, xmlLine );

   Stg_asprintf( &testEntries, "%s%s", struct1Entry, struct2Entry );
   _IO_HandlerSuite_CreateTestXMLFile( xmlTestFileName, testEntries );
   Memory_Free( testEntries );
   IO_Handler_ReadAllFromFile( data->io_handler, xmlTestFileName, data->dict2 );

   pcu_check_true( 1 == data->dict2->count );
   pcu_check_true( Dictionary_Entry_Compare( data->dict2->entryPtr[0],
      (Dictionary_Entry_Key)struct1Name) );
   structDev = data->dict2->entryPtr[0]->value;
   structDict = structDev->as.typeStruct;
   pcu_check_true( Dictionary_Entry_Value_Type_Struct == structDev->type );
   pcu_check_true( struct1_OrigParamCount*2 == Dictionary_Entry_Value_GetCount( structDev ) );
   pcu_check_true( 0 == strcmp( structDict->entryPtr[0]->key, paramNames[0] ) );
   elementDev = structDict->entryPtr[0]->value;
   pcu_check_true( paramVals[0] == Dictionary_Entry_Value_AsUnsignedInt( elementDev ) );
   pcu_check_true( 0 == strcmp( structDict->entryPtr[1]->key, paramNames[1] ) );
   elementDev = structDict->entryPtr[1]->value;
   pcu_check_true( paramVals[1] == Dictionary_Entry_Value_AsUnsignedInt( elementDev ) );
   pcu_check_true( 0 == strcmp( structDict->entryPtr[2]->key, paramNames[1] ) );
   elementDev = structDict->entryPtr[2]->value;
   pcu_check_true( paramVals2[1] == Dictionary_Entry_Value_AsUnsignedInt( elementDev ) );
   pcu_check_true( 0 == strcmp( structDict->entryPtr[3]->key, paramNames2[0] ) );
   elementDev = structDict->entryPtr[3]->value;
   pcu_check_true( paramVals2[0] == Dictionary_Entry_Value_AsUnsignedInt( elementDev ) );
   remove( xmlTestFileName );
   Dictionary_Empty( data->dict2 );

   /* Sub-test 3.2: with mergeType as "merge", structs to be merged.
    * childrenMergeType set to merge also */
   sprintf( struct2Entry, "" );
   sprintf( xmlLine, "<struct name=\"%s\" mergeType=\"merge\" childrenMergeType=\"merge\">\n", struct1Name );
   strcat( struct2Entry, xmlLine );
   sprintf( xmlLine, "<param name=\"%s\">%u</param>\n", paramNames[1], paramVals2[1] );
   strcat( struct2Entry, xmlLine );
   sprintf( xmlLine, "<param name=\"%s\">%u</param>\n", paramNames2[0], paramVals2[0] );
   strcat( struct2Entry, xmlLine );
   sprintf( xmlLine, "</struct>\n" );
   strcat( struct2Entry, xmlLine );

   Stg_asprintf( &testEntries, "%s%s", struct1Entry, struct2Entry );
   _IO_HandlerSuite_CreateTestXMLFile( xmlTestFileName, testEntries );
   Memory_Free( testEntries );
   IO_Handler_ReadAllFromFile( data->io_handler, xmlTestFileName, data->dict2 );

   pcu_check_true( 1 == data->dict2->count );
   pcu_check_true( Dictionary_Entry_Compare( data->dict2->entryPtr[0],
      (Dictionary_Entry_Key)struct1Name) );
   structDev = data->dict2->entryPtr[0]->value;
   structDict = structDev->as.typeStruct;
   pcu_check_true( Dictionary_Entry_Value_Type_Struct == structDev->type );
   pcu_check_true( struct1_OrigParamCount+1 == Dictionary_Entry_Value_GetCount( structDev ) );
   pcu_check_true( 0 == strcmp( structDict->entryPtr[0]->key, paramNames[0] ) );
   elementDev = structDict->entryPtr[0]->value;
   pcu_check_true( paramVals[0] == Dictionary_Entry_Value_AsUnsignedInt( elementDev ) );
   pcu_check_true( 0 == strcmp( structDict->entryPtr[1]->key, paramNames[1] ) );
   elementDev = structDict->entryPtr[1]->value;
   pcu_check_true( paramVals2[1] == Dictionary_Entry_Value_AsUnsignedInt( elementDev ) );
   pcu_check_true( 0 == strcmp( structDict->entryPtr[2]->key, paramNames2[0] ) );
   elementDev = structDict->entryPtr[2]->value;
   pcu_check_true( paramVals2[0] == Dictionary_Entry_Value_AsUnsignedInt( elementDev ) );
   remove( xmlTestFileName );
   Dictionary_Empty( data->dict2 );
}


void IO_HandlerSuite_TestReadNonExistent( IO_HandlerSuiteData* data ) {
   const char*    errorFileName = "./errorMsg-NonExist.txt";
   const char*    notExistFilename = "I_Dont_Exist.xml";
   FILE*          errorFile;
   #define        MAXLINE 1000
   char           errorLine[MAXLINE];
   char           expectedErrorMsg[MAXLINE];

   Stream_RedirectFile( Journal_Register( Error_Type, XML_IO_Handler_Type ), errorFileName );
   pcu_check_assert( IO_Handler_ReadAllFromFile( data->io_handler, notExistFilename, data->dict2 ) );
   pcu_check_true(errorFile = fopen( errorFileName, "r" ));
   if (errorFile) {
      pcu_check_true( fgets( errorLine, MAXLINE, errorFile ) );
      sprintf( expectedErrorMsg, "Error: File %s doesn't exist, not readable, or not valid.\n",
         notExistFilename );
      pcu_check_true( 0 == strcmp( errorLine, expectedErrorMsg ) );
      remove( errorFileName );
   }
   
   Dictionary_Empty( data->dict2 );
}


void errorHandler( void* ctx, const char* msg, ... )
{
   va_list ap;
   va_start( ap, msg );
   Stream_Printf(
      Journal_Register( Error_Type, XML_IO_Handler_Type ),
      msg, ap );
   va_end(ap);
}


void IO_HandlerSuite_TestReadInvalid( IO_HandlerSuiteData* data ) {
   const char*    errorFileName = "./errorMsg-Invalid.txt";
   const char*    invalidXMLFilename = "Invalid.xml";
   FILE*          errorFile;
   #define        MAXLINE 1000
   char           errorLine[MAXLINE];
   char           expectedErrorMsg[MAXLINE];

   Stream_RedirectFile( Journal_Register( Error_Type, XML_IO_Handler_Type ), errorFileName );

   /* Create an invalid XML file */
   _IO_HandlerSuite_CreateTestXMLFile( invalidXMLFilename, "<invalid></param>\n" );

   pcu_check_assert( IO_Handler_ReadAllFromFile( data->io_handler, invalidXMLFilename, data->dict2 ) );
   pcu_check_true( (errorFile = fopen( errorFileName, "r" )) );
   if ( errorFile ) { 
      sprintf( expectedErrorMsg, ".//%s:3: parser error : Opening and ending tag mismatch: "
         "invalid line 3 and param\n",
         invalidXMLFilename );
      pcu_check_true( fgets( errorLine, MAXLINE, errorFile ) );
      pcu_check_true( 0 == strcmp( errorLine, expectedErrorMsg ) );
      sprintf( expectedErrorMsg, "<invalid></param>\n" );
      pcu_check_true( fgets( errorLine, MAXLINE, errorFile ) );
      pcu_check_true( 0 == strcmp( errorLine, expectedErrorMsg ) );
      sprintf( expectedErrorMsg, "                 ^\n" );
      pcu_check_true( fgets( errorLine, MAXLINE, errorFile ) );
      pcu_check_true( 0 == strcmp( errorLine, expectedErrorMsg ) );
      sprintf( expectedErrorMsg, "Error: File %s doesn't exist, not readable, or not valid.\n",
         invalidXMLFilename );
      pcu_check_true( fgets( errorLine, MAXLINE, errorFile ) );
      pcu_check_true( 0 == strcmp( errorLine, expectedErrorMsg ) );
      remove( errorFileName );
   }
   remove( invalidXMLFilename );

   Dictionary_Empty( data->dict2 );
}


void IO_HandlerSuite_TestReadWrongNS( IO_HandlerSuiteData* data ) {
   const char*    wrongNS_XMLFilename = "WrongNS.xml";
   FILE*          wrongNS_XMLFile = NULL;
   const char*    wrongNS_Line = "<StGermainData xmlns=\"http://www.wrong.com/StGermain/XML_IO_Handler/Jun2003\">";
   const char*    errorFileName = "./errorMsg-wrongNS.txt";
   FILE*          errorFile;
   #define        MAXLINE 1000
   char           errorLine[MAXLINE];
   char           expectedErrorMsg[MAXLINE];

   Stream_RedirectFile( Journal_Register( Error_Type, XML_IO_Handler_Type ), errorFileName );

   /* Create an XML file with wrong NS - thus can't use our standard test XML creator */
   wrongNS_XMLFile = fopen(wrongNS_XMLFilename, "w");
   fwrite( IO_HandlerSuite_XMLStartString1, sizeof(char),
      strlen( IO_HandlerSuite_XMLStartString1 ), wrongNS_XMLFile );
   fwrite( wrongNS_Line, sizeof(char), strlen( wrongNS_Line ), wrongNS_XMLFile );
   fwrite( IO_HandlerSuite_XMLEndString, sizeof(char),
      strlen( IO_HandlerSuite_XMLEndString ), wrongNS_XMLFile );
   fclose( wrongNS_XMLFile );

   pcu_check_assert( IO_Handler_ReadAllFromFile( data->io_handler, wrongNS_XMLFilename, data->dict2 ) );
   pcu_check_true( errorFile = fopen( errorFileName, "r" ));
   if (errorFile) {
      sprintf( expectedErrorMsg, "Error: resource .//%s of the wrong type, unknown namespace "
         "wasn't expected value of http://www.vpac.org/StGermain/XML_IO_Handler/Jun2003.\n",
         wrongNS_XMLFilename );
      pcu_check_true( fgets( errorLine, MAXLINE, errorFile ) );
      pcu_check_true( 0 == strcmp( errorLine, expectedErrorMsg ) );
      sprintf( expectedErrorMsg, "Error: File %s not valid/readable.\n",
         wrongNS_XMLFilename );
      pcu_check_true( fgets( errorLine, MAXLINE, errorFile ) );
      pcu_check_true( 0 == strcmp( errorLine, expectedErrorMsg ) );
      remove( errorFileName );
   }
   remove( wrongNS_XMLFilename );
   Dictionary_Empty( data->dict2 );
}


void IO_HandlerSuite_TestReadWrongRootNode( IO_HandlerSuiteData* data ) {
   const char*    wrongRootNode_XMLFilename = "WrongRootNode.xml";
   FILE*          wrongRootNode_XMLFile = NULL;
   const char*    wrongRootNode_Line = "<Wrong xmlns=\"http://www.vpac.org/StGermain/XML_IO_Handler/Jun2003\">\n";
   const char*    wrongRootNode_EndLine = "</Wrong>";
   const char*    errorFileName = "./errorMsg-wrongRootNode.txt";
   FILE*          errorFile;
   #define        MAXLINE 1000
   char           errorLine[MAXLINE];
   char           expectedErrorMsg[MAXLINE];

   Stream_RedirectFile( Journal_Register( Error_Type, XML_IO_Handler_Type ), errorFileName );

   /* Create an XML file with wrong NS - thus can't use our standard test XML creator */
   wrongRootNode_XMLFile = fopen(wrongRootNode_XMLFilename, "w");
   fwrite( IO_HandlerSuite_XMLStartString1, sizeof(char),
      strlen( IO_HandlerSuite_XMLStartString1 ), wrongRootNode_XMLFile );
   fwrite( wrongRootNode_Line, sizeof(char), strlen( wrongRootNode_Line ), wrongRootNode_XMLFile );
   fwrite( wrongRootNode_EndLine, sizeof(char), strlen( wrongRootNode_EndLine ), wrongRootNode_XMLFile );
   fclose( wrongRootNode_XMLFile );

   pcu_check_assert( IO_Handler_ReadAllFromFile( data->io_handler, wrongRootNode_XMLFilename, data->dict2 ) );
   pcu_check_true( (errorFile = fopen( errorFileName, "r" )) );
   if ( errorFile ) {
      sprintf( expectedErrorMsg, "resource .//%s of wrong type, root node "
         "=%s, should be <StGermainData>.\n", wrongRootNode_XMLFilename, "<Wrong>" );
      pcu_check_true( fgets( errorLine, MAXLINE, errorFile ) );
      pcu_check_true( 0 == strcmp( errorLine, expectedErrorMsg ) );
      sprintf( expectedErrorMsg, "Not parsing.\n" );
      pcu_check_true( fgets( errorLine, MAXLINE, errorFile ) );
      pcu_check_true( 0 == strcmp( errorLine, expectedErrorMsg ) );
      sprintf( expectedErrorMsg, "Error: File %s not valid/readable.\n",
         wrongRootNode_XMLFilename );
      pcu_check_true( fgets( errorLine, MAXLINE, errorFile ) );
      pcu_check_true( 0 == strcmp( errorLine, expectedErrorMsg ) );
      remove( errorFileName );
   }
   remove( wrongRootNode_XMLFilename );
   Dictionary_Empty( data->dict2 );
}


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
   pcu_suite_addTest( suite, IO_HandlerSuite_TestReadAllFromCommandLine );
   pcu_suite_addTest( suite, IO_HandlerSuite_TestReadDuplicateEntryKeys );
   pcu_suite_addTest( suite, IO_HandlerSuite_TestReadNonExistent );
   pcu_suite_addTest( suite, IO_HandlerSuite_TestReadInvalid );
   pcu_suite_addTest( suite, IO_HandlerSuite_TestReadWrongNS );
   pcu_suite_addTest( suite, IO_HandlerSuite_TestReadWrongRootNode );
}
