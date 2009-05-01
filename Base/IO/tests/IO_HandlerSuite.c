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

#include "pcu/pcu.h"
#include "StGermain/Base/Foundation/Foundation.h"
#include "StGermain/Base/IO/IO.h"
#include "DictionarySuite.h"  /* Because want to re-use sample dictionary structs/funcs */
#include "IO_HandlerSuite.h"


const char* IO_HandlerSuite_XMLStartString1 = "<?xml version=\"1.0\"?>\n";
const char* IO_HandlerSuite_XMLStartString2 = "<StGermainData xmlns=\"http://www.vpac.org/StGermain/XML_IO_Handler/Jun2003\">\n";
const char* IO_HandlerSuite_XMLEndString = "</StGermainData>\n";

typedef struct {
   XML_IO_Handler*                io_handler;
   Dictionary*                    dict1;
   Dictionary*                    dict2;
   DictionarySuite_TestDictData*  testDD;
} IO_HandlerSuiteData;


void IO_HandlerSuite_Setup( IO_HandlerSuiteData* data ) {
   data->io_handler = XML_IO_Handler_New();
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

   Dictionary_Empty( data->dict1 );
   Dictionary_Empty( data->dict2 );
   DictionarySuite_PopulateDictWithTestValues( data->dict1, data->testDD );

   IO_Handler_WriteAllToFile( data->io_handler, "dictTest.xml", data->dict1 );
   IO_Handler_ReadAllFromFile( data->io_handler, "dictTest.xml", data->dict2 ); 

   pcu_check_true( data->dict1->count == data->dict2->count );
   for (ii=0; ii<data->dict1->count; ii++) {
      pcu_check_true( Dictionary_Entry_Compare( data->dict1->entryPtr[ii],
         data->dict2->entryPtr[ii]->key) );
      pcu_check_true( Dictionary_Entry_Value_Compare( data->dict1->entryPtr[ii]->value,
         data->dict2->entryPtr[ii]->value) );
   }

   remove("dictTest.xml");
   Dictionary_Empty( data->dict1 );
   Dictionary_Empty( data->dict2 );
}


/* Similar to above test, except using the function to write just one entry at a time */
void IO_HandlerSuite_TestWriteNormalSingleEntry( IO_HandlerSuiteData* data ) {
   Index         ii;
   const char*   fileName = "singleEntry.xml";

   Dictionary_Empty( data->dict2 );

   for (ii=0; ii<data->dict1->count; ii++) {
      XML_IO_Handler_WriteEntryToFile( data->io_handler, fileName,
            data->testDD->testKeys[ii],
            data->testDD->testValues[ii], 
            NULL );
      IO_Handler_ReadAllFromFile( data->io_handler, fileName, data->dict2 ); 

      pcu_check_true( 1 == data->dict2->count );
      pcu_check_true( Dictionary_Entry_Compare( data->dict2->entryPtr[1],
         data->testDD->testKeys[ii]) );
      pcu_check_true( Dictionary_Entry_Value_Compare( data->dict2->entryPtr[1]->value,
         data->testDD->testValues[ii] ) );

      Dictionary_Empty( data->dict2 );
      remove(fileName);
   }

   Dictionary_Empty( data->dict2 );
}


/* In this case, want to make sure the types are written explicitly into the output, so will
 * check against expected text. */
void IO_HandlerSuite_TestWriteExplicitTypes( IO_HandlerSuiteData* data ) {
   Index         ii=0;
   const char*   testFileName = "dictTest-explicittypes.xml";
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
   IO_Handler_WriteAllToFile( data->io_handler, "dictTest-explicittypes.xml", data->dict1 );

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

//TODO: tests for reading input with various whitespaces, bad values etc

void IO_HandlerSuite( pcu_suite_t* suite ) {
   pcu_suite_setData( suite, IO_HandlerSuiteData );
   pcu_suite_setFixtures( suite, IO_HandlerSuite_Setup, IO_HandlerSuite_Teardown );
   pcu_suite_addTest( suite, IO_HandlerSuite_TestWriteReadNormalEntries );
   pcu_suite_addTest( suite, IO_HandlerSuite_TestWriteNormalSingleEntry );
   pcu_suite_addTest( suite, IO_HandlerSuite_TestWriteExplicitTypes );
}
