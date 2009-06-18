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
#include "StGermain/Base/Extensibility/Extensibility.h"
#include "StGermain/Base/Context/Context.h"
#include "DictionaryCheckSuite.h"

typedef struct {
   unsigned int   rank;
} DictionaryCheckSuiteData;

void DictionaryCheckSuite_Setup( DictionaryCheckSuiteData* data ) {
   MPI_Comm_rank( MPI_COMM_WORLD, &data->rank );
}

void DictionaryCheckSuite_Teardown( DictionaryCheckSuiteData* data ) {
}

  
void DictionaryCheckSuite_TestCheckKeys( DictionaryCheckSuiteData* data ) {
   Dictionary*       dictionary = Dictionary_New();
   Dictionary*       dictionary2 = Dictionary_New();
   const char*       testFilename = "./testDictionaryCheck.txt";
   FILE*             testFile = NULL;
   #define           MAXLINE 1000
   char              inputLine[MAXLINE];
   Dictionary_Index   index;
   char*             errMessage = "Component dictionary must have unique names\n";
   
   Stream_RedirectFile(Journal_Register( Error_Type, "DictionaryCheck"), testFilename );
   Stream_SetPrintingRank(Journal_Register( Error_Type, "DictionaryCheck"), 0 );
   Stream_ClearCustomFormatters( Journal_Register( Error_Type, "DictionaryCheck") );

   /* Create a set of Dictionary entries */
   /* For dictionary */
   Dictionary_Add( dictionary, "test_dict_string",
      Dictionary_Entry_Value_FromString( "hello" ) );
   Dictionary_Add( dictionary, "test_dict_double",
      Dictionary_Entry_Value_FromDouble( 45.567 ) );
   Dictionary_Add( dictionary, "test_dict_string",
      Dictionary_Entry_Value_FromString( "goodbye" ) );   
   Dictionary_Add( dictionary, "test_dict_string",
      Dictionary_Entry_Value_FromString( "hello" ) );
   Dictionary_Add( dictionary, "test_dict_string2",
      Dictionary_Entry_Value_FromString( "hello" ) );
   
   CheckDictionaryKeys(dictionary,  errMessage);

   if ( data->rank==0 ) {
      testFile = fopen( testFilename, "r" );
      pcu_check_true( fgets( inputLine, MAXLINE, testFile ));
      pcu_check_streq( inputLine, errMessage );
      pcu_check_true( fgets( inputLine, MAXLINE, testFile ));
      /* Ignore the actual 2nd line, since it includes a memory ptr print and hard to compare */
      pcu_check_true( fgets( inputLine, MAXLINE, testFile ));
      pcu_check_streq( inputLine, "The following keys were repeated:\n" );
      pcu_check_true( fgets( inputLine, MAXLINE, testFile ));
      pcu_check_streq( inputLine, "\t\"test_dict_string\"\n" );
      pcu_check_true( fgets( inputLine, MAXLINE, testFile ));
      pcu_check_streq( inputLine, "Error in CheckDictionaryKeys with 1 entries in dictionary keys\n" );
   }

   /* For dictionary2 */
   Dictionary_Add( dictionary2, "test_dict_string",
      Dictionary_Entry_Value_FromString( "hello" ) );
   Dictionary_Add( dictionary2, "test_dict_double",
      Dictionary_Entry_Value_FromDouble( 45.567 ) );
   Dictionary_Add( dictionary2, "test_dict_stuff",
      Dictionary_Entry_Value_FromString( "hello") );

   /* Call DictionaryCheck function */
   CheckDictionaryKeys(dictionary2, errMessage);

   /* there shouldn't be any more error lines printed about the 2nd dictionary */
   if ( data->rank==0 ) {
      pcu_check_true( NULL == fgets( inputLine, MAXLINE, testFile ) );
   }
   
   Stg_Class_Delete( dictionary );
   Stg_Class_Delete( dictionary2 );

   if ( data->rank==0 ) {
      fclose( testFile );
      remove( testFilename );
   }
}


void DictionaryCheckSuite( pcu_suite_t* suite ) {
   pcu_suite_setData( suite, DictionaryCheckSuiteData );
   pcu_suite_setFixtures( suite, DictionaryCheckSuite_Setup, DictionaryCheckSuite_Teardown );
   pcu_suite_addTest( suite, DictionaryCheckSuite_TestCheckKeys );
}
