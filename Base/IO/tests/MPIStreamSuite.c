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
#include <mpi.h>

#include "pcu/pcu.h"
#include "StGermain/Base/Foundation/Foundation.h"
#include "StGermain/Base/IO/IO.h"
#include "MPIStreamSuite.h"

#define MAXLINE 1000

typedef struct {
   MPI_Comm       comm;
   int            rank;
   int            numProcs;
   Dictionary*    dict;
} MPIStreamSuiteData;


void MPIStreamSuite_Setup( MPIStreamSuiteData* data ) {
   MPI_Comm_dup( MPI_COMM_WORLD, &data->comm );
   MPI_Comm_size( data->comm, &data->numProcs );
   MPI_Comm_rank( data->comm, &data->rank );
   data->dict = Dictionary_New();
}

void MPIStreamSuite_Teardown( MPIStreamSuiteData* data ) {
   Stg_Class_Delete( data->dict );
}
   

void MPIStreamSuite_TestWriteAllProcessors( MPIStreamSuiteData* data ) {
   Index          ii;
   Stream*        stream1;
   const char*    testMPIFilename = "./test-mpi1.txt";
   FILE*          testMPIFile;
   char           outLine[MAXLINE];
   char           compString[MAXLINE];
   #define        PER_RANK_COUNT 2
   char           dataArray[PER_RANK_COUNT];
   Index          rank_I;

   Dictionary_AddFromString( data->dict, "journal-file.MPIStream.one", testMPIFilename );
   Dictionary_AddFromUnsignedInt( data->dict, "journal-mpi-offset.MPIStream.one", 100 );
   Journal_ReadFromDictionary( data->dict );

   stream1 = Journal_Register( MPIStream_Type, "one" );
   pcu_check_true( STREAM_ALL_RANKS == Stream_GetPrintingRank( stream1 ) );

   /* Write half the alphabet to each process */
   for ( ii = 0; ii < PER_RANK_COUNT; ++ii ) {
      dataArray[ii] = 'a' + PER_RANK_COUNT * data->rank + ii;
   }
   
   /* Print the data */
   MPIStream_WriteAllProcessors( stream1, dataArray, sizeof(char), PER_RANK_COUNT, data->comm );
   Stream_Flush( stream1 );

   /* Now build up the comparison string. Depends on how many processes are running. */
   for (rank_I=0; rank_I < data->numProcs; rank_I++ ) {
      for ( ii = 0; ii < PER_RANK_COUNT; ++ii ) {
         compString[PER_RANK_COUNT*rank_I+ii] = 'a' + PER_RANK_COUNT*rank_I + ii;
      }
   }
   compString[data->numProcs*PER_RANK_COUNT] = '\0';

   /* Do the following since in parallel on some systems, the file
    * doesn't get re-opened at the start automatically. */
   for ( rank_I = 0; rank_I < data->numProcs; rank_I++ ) {
      MPI_Barrier( data->comm );
      if ( rank_I == data->rank ) {
         testMPIFile = fopen( testMPIFilename, "r" );
         rewind( testMPIFile );
      }
   }

   pcu_check_true( fgets( outLine, MAXLINE, testMPIFile ));
   pcu_check_streq( outLine, compString );
   fclose( testMPIFile );
   /* Make sure not to delete the file until all processors have read from it */
   MPI_Barrier( data->comm );
   if ( data->rank == 0 ) {
      remove( testMPIFilename );
   }
}


void MPIStreamSuite_TestPrintWithOffset( MPIStreamSuiteData* data ) {
   Stream*     stream2;
   const char* testMPIFilename = "./test-mpi2.txt";
   FILE*       testMPIFile;
   char        outLine[MAXLINE];
   char        rankPrintString[MAXLINE];
   char        compString[MAXLINE];
   Index       rank_I;
   Index       ii;
   Index       startPoint = 0;
   unsigned    stringLength=0;

   Dictionary_AddFromString( data->dict, "journal-file.MPIStream.two", testMPIFilename );
   Journal_ReadFromDictionary( data->dict );

   stream2 = Journal_Register( MPIStream_Type, "two" );
   pcu_check_true( STREAM_ALL_RANKS == Stream_GetPrintingRank( stream2 ) );

   for( ii=0; ii <10*(data->rank+1); ii++ ) {
      rankPrintString[ii] = 'a' + data->rank;
   }
   rankPrintString[10*(data->rank+1)] = '\0';

   stringLength = 10*(data->rank+1);

   MPIStream_SetOffset( stream2, stringLength*sizeof(char), data->comm );

   /* now print the data */
   Journal_Printf( stream2, "%s", rankPrintString );
   Stream_Flush( stream2 );

   /* Now build up the comparison string. Depends on how many processes are running. */
   startPoint = 0;
   for (rank_I=0; rank_I < data->numProcs; rank_I++ ) {
      for ( ii = 0; ii < 10*(rank_I+1); ii++ ) {
         compString[startPoint+ii] = 'a' + rank_I;
      }
      startPoint += 10*(rank_I+1);
   }
   compString[startPoint] = '\0';

   /* Do the following since in parallel on some systems, the file
    * doesn't get re-opened at the start automatically. */
   for ( rank_I = 0; rank_I < data->numProcs; rank_I++ ) {
      MPI_Barrier( data->comm );
      if ( rank_I == data->rank ) {
         testMPIFile = fopen( testMPIFilename, "r" );
         rewind( testMPIFile );
      }
   }

   pcu_check_true( fgets( outLine, MAXLINE, testMPIFile ));
   pcu_check_streq( outLine, compString );
   fclose( testMPIFile );
   /* Make sure not to delete the file until all processors have read from it */
   MPI_Barrier( data->comm );
   if ( data->rank == 0 ) {
      remove( testMPIFilename );
   }
}


void MPIStreamSuite( pcu_suite_t* suite ) {
   pcu_suite_setData( suite, MPIStreamSuiteData );
   pcu_suite_setFixtures( suite, MPIStreamSuite_Setup, MPIStreamSuite_Teardown );
   pcu_suite_addTest( suite, MPIStreamSuite_TestWriteAllProcessors );
   pcu_suite_addTest( suite, MPIStreamSuite_TestPrintWithOffset );
}


