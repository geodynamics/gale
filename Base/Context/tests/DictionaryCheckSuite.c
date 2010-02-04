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
   int rank;
} DictionaryCheckSuiteData;

void DictionaryCheckSuite_Setup( DictionaryCheckSuiteData* data ) {
   MPI_Comm_rank( MPI_COMM_WORLD, &data->rank );
}

void DictionaryCheckSuite_Teardown( DictionaryCheckSuiteData* data ) {
}

  
void DictionaryCheckSuite_TestCheckKeys( DictionaryCheckSuiteData* data ) {
   Dictionary*       dictionary = Dictionary_New();
   Dictionary*       dictionary2 = Dictionary_New();
   const Name       testFilename1 = "testDictionaryCheck-1.txt";
   const Name       testFilename2 = "testDictionaryCheck-2.txt";
   char              expectedFilename[PCU_PATH_MAX];
   const char*       errMessage = "Component dictionary must have unique names\n";
   
   Stream_RedirectFile(Journal_Register( Error_Type, (Name)"DictionaryCheck" ), testFilename1 );
   Stream_SetPrintingRank(Journal_Register( Error_Type, (Name)"DictionaryCheck" ), 0 );
   Stream_ClearCustomFormatters( Journal_Register( Error_Type, (Name)"DictionaryCheck")  );

   /* Create a set of Dictionary entries */
   /* For dictionary */
   Dictionary_Add( dictionary, (Dictionary_Entry_Key)"test_dict_string", Dictionary_Entry_Value_FromString( "hello" )  );
   Dictionary_Add( dictionary, (Dictionary_Entry_Key)"test_dict_double", Dictionary_Entry_Value_FromDouble( 45.567 )  );
   Dictionary_Add( dictionary, (Dictionary_Entry_Key)"test_dict_string", Dictionary_Entry_Value_FromString( "goodbye" )  );   
   Dictionary_Add( dictionary, (Dictionary_Entry_Key)"test_dict_string", Dictionary_Entry_Value_FromString( "hello" )  );
   Dictionary_Add( dictionary, (Dictionary_Entry_Key)"test_dict_string2", Dictionary_Entry_Value_FromString( "hello" )  );
   
   CheckDictionaryKeys(dictionary,  errMessage);

   if ( data->rank==0 ) {
      pcu_filename_expected( testFilename1, expectedFilename );
      pcu_check_fileEq( testFilename1, expectedFilename );
      remove( testFilename1 );
   }

   /* For dictionary2 */
   Dictionary_Add( dictionary2, (Dictionary_Entry_Key)"test_dict_string", Dictionary_Entry_Value_FromString( "hello" )  );
   Dictionary_Add( dictionary2, (Dictionary_Entry_Key)"test_dict_double", Dictionary_Entry_Value_FromDouble( 45.567 )  );
   Dictionary_Add( dictionary2, (Dictionary_Entry_Key)"test_dict_stuff", Dictionary_Entry_Value_FromString( "hello")  );

   /* Call DictionaryCheck function */
   Stream_RedirectFile(Journal_Register( Error_Type, (Name)"DictionaryCheck" ), testFilename2 );
   CheckDictionaryKeys(dictionary2, errMessage);

   /* This file expected to be empty */
   if ( data->rank==0 ) {
      pcu_filename_expected( testFilename2, expectedFilename );
      pcu_check_fileEq( testFilename2, expectedFilename );
      remove( testFilename2 );
   }
   
   Stg_Class_Delete( dictionary );
   Stg_Class_Delete( dictionary2 );

   if ( data->rank==0 ) {
      remove( testFilename1 );
      remove( testFilename2 );
   }
}


void DictionaryCheckSuite( pcu_suite_t* suite ) {
   pcu_suite_setData( suite, DictionaryCheckSuiteData );
   pcu_suite_setFixtures( suite, DictionaryCheckSuite_Setup, DictionaryCheckSuite_Teardown );
   pcu_suite_addTest( suite, DictionaryCheckSuite_TestCheckKeys );
}


