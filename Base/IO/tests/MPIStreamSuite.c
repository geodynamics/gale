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
#include "MPIStreamSuite.h"

typedef struct {
   MPI_Comm       CommWorld;
   int            rank;
   int            numProcessors;
   Dictionary*    dict;
} MPIStreamSuiteData;


void MPIStreamSuite_Setup( MPIStreamSuiteData* data ) {
   MPI_Comm_dup( MPI_COMM_WORLD, &data->CommWorld );
   MPI_Comm_size( data->CommWorld, &data->numProcessors );
   MPI_Comm_rank( data->CommWorld, &data->rank );
   data->dict = Dictionary_New();
}

void MPIStreamSuite_Teardown( MPIStreamSuiteData* data ) {
   Stg_Class_Delete( data->dict );
}
   

void MPIStreamSuite_TestWriteAllProcessors( MPIStreamSuiteData* data ) {
   char        dataArray[13];
   Index       ii;
   Stream*     stream1;
   const char* testMPIFilename = "./test-mpi1.txt";
   FILE*       testMPIFile;
   #define     MAXLINE 1000
   char        outLine[MAXLINE];

   Journal_Printf( Journal_Register( Error_Type, "MPIStreamSuite" ),
      "Warning-MPIStreamSuite tests disabled since they currently produce segfaults,\n"
      "seems to be a problem in the MPIStream or MPIFile code\n" );
   pcu_check_true( 0 );
   return;
   /* Will not run below......................................................*/

   Dictionary_AddFromString( data->dict, "journal-file.MPIStream.one", testMPIFilename );
   Dictionary_AddFromUnsignedInt( data->dict, "journal-mpi-offset.MPIStream.one", 100 );
   Journal_ReadFromDictionary( data->dict );

   stream1 = Journal_Register( MPIStream_Type, "one" );
   pcu_check_true( STREAM_ALL_RANKS == Stream_GetPrintingRank( stream1 ) );

   /* Write half the alphabet to each process */
   for ( ii = 0; ii < 13; ++ii ) {
      dataArray[ii] = 'a' + 13 * data->rank + ii;
   }
   
   /* Print the alphabet */
   MPIStream_WriteAllProcessors( stream1, dataArray, sizeof(char), 13, data->CommWorld );
   Stream_Flush( stream1 );

   testMPIFile = fopen( testMPIFilename, "r" );
   pcu_check_true( fgets( outLine, MAXLINE, testMPIFile ));
   pcu_check_true( 0 == strcmp( outLine, "abcdefghijklmnopqrstuvwxyz" ));
   fclose( testMPIFile );
   remove( testMPIFilename );
}

   
void MPIStreamSuite_TestPrintWithOffset( MPIStreamSuiteData* data ) {
   Stream* stream2;
   const char* testMPIFilename = "./test-mpi2.txt";
   FILE*       testMPIFile;
   #define     MAXLINE 1000
   char        outLine[MAXLINE];

   Journal_Printf( Journal_Register( Error_Type, "MPIStreamSuite" ),
      "Warning-MPIStreamSuite tests disabled since they currently produce segfaults,\n"
      "seems to be a problem in the MPIStream or MPIFile code\n" );
   pcu_check_true( 0 );
   return;
   /* Will not run below......................................................*/

   Dictionary_AddFromString( data->dict, "journal-file.MPIStream.two", testMPIFilename );
   Journal_ReadFromDictionary( data->dict );

   stream2 = Journal_Register( MPIStream_Type, "two" );
   pcu_check_true( STREAM_ALL_RANKS == Stream_GetPrintingRank( stream2 ) );

   switch ( data->rank ) {
      case 0:
         MPIStream_SetOffset( stream2, 8, data->CommWorld );
         Journal_Printf( stream2, "Hello %d", 10 );
         break;
      case 1:
         MPIStream_SetOffset( stream2, 12, data->CommWorld );
         Journal_Printf( stream2, "world %.3f\n", 0.123f );
         break;
   }
   Stream_Flush( stream2 );

   testMPIFile = fopen( testMPIFilename, "r" );
   pcu_check_true( fgets( outLine, MAXLINE, testMPIFile ));
   pcu_check_true( 0 == strcmp( outLine, "Hello 10world 0.123" ));
   fclose( testMPIFile );
   remove( testMPIFilename );
}


void MPIStreamSuite( pcu_suite_t* suite ) {
   pcu_suite_setData( suite, MPIStreamSuiteData );
   pcu_suite_setFixtures( suite, MPIStreamSuite_Setup, MPIStreamSuite_Teardown );
   pcu_suite_addTest( suite, MPIStreamSuite_TestWriteAllProcessors );
   pcu_suite_addTest( suite, MPIStreamSuite_TestPrintWithOffset );
}
