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
**   Tests accuracy of memory statistics generation.
**
** $Id: testMemory2.c 3462 2006-02-19 06:53:24Z WalterLandry $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#include "pcu/pcu.h"
#include "StGermain/Base/Foundation/Foundation.h"
#include "StGermain/Base/IO/IO.h" /* For Journal stuff */
#include "MemoryReportSuite.h"

struct StructA
{
   int x;
   float y;
   char z;
};
typedef struct StructA StructA;

struct StructB
{
   double x;
};
typedef struct StructB StructB;

struct StructC
{
   char* x;
   StructA a;
};
typedef struct StructC StructC;

typedef struct {
   Memory*                 savedStgMemory;
   MemoryReport*           report;
   int							rank;   
   void*                   bytesObj;
   void*                   bytesArray;
   StructA*                object;
   StructB*                array1d;
   StructC**               array2d;
   StructA***              array3d;
   StructB****             array4d;
   StructC*                one2d;
   StructA*                one3d;
   StructB*                one4d;
   StructC**               complex2d;
   Index**                 setup;
   StructA***              complex3d;
   unsigned                strA_alloc;
   unsigned                strB_alloc;
   unsigned                strC_alloc;
   unsigned                bytes_alloc;
   unsigned                index_alloc;
   unsigned                strA_total;
   unsigned                strB_total;
   unsigned                strC_total;
   unsigned                bytes_total;
   unsigned                index_total;
   unsigned                groupTwo_strA_alloc;
   unsigned                groupTwo_strB_alloc;
   unsigned                groupTwo_strC_alloc;
   unsigned                groupTwo_strA_total;
   unsigned                groupTwo_strB_total;
   unsigned                groupTwo_strC_total;
   unsigned                groupOneFunc_alloc;
   unsigned                groupTwoFunc_alloc;
   unsigned                groupTwoName_alloc;
   unsigned                bytesFunc_alloc;
   unsigned                groupOneFunc_total;
   unsigned                groupTwoFunc_total;
   unsigned                groupTwoName_total;
   unsigned                bytesFunc_total;
   unsigned                file_alloc;
   unsigned                file_total;
} MemoryReportSuiteData;


void MemoryReportSuite_AllocGroupOne( MemoryReportSuiteData* data ) {
   data->array1d = Memory_Alloc_Array( StructB, 3, "GroupOne" );
   data->strB_alloc++;
   data->strB_total += sizeof(StructB)*3;

   data->array2d = Memory_Alloc_2DArray( StructC, 4, 5, "GroupOne" );
   data->strC_alloc++;
   data->strC_total += Memory_Length_2DArray( sizeof(StructC), 4, 5 );

   data->array3d = Memory_Alloc_3DArray( StructA, 2, 3, 4, "GroupOne" );
   data->strA_alloc++;
   data->strA_total += Memory_Length_3DArray( sizeof(StructA), 2, 3, 4 );

   data->array4d = Memory_Alloc_4DArray_Unnamed( StructB, 5, 4, 3, 2 );
   data->strB_alloc++;
   data->strB_total += Memory_Length_4DArray( sizeof(StructB), 5, 4, 3, 2 );

   data->one2d = Memory_Alloc_2DArrayAs1D_Unnamed( StructC, 4, 2 );
   data->strC_alloc++;
   data->strC_total += Memory_Length_2DAs1D( sizeof(StructC), 4, 2 );

   data->one3d = Memory_Alloc_3DArrayAs1D_Unnamed( StructA, 2, 2, 3 );
   data->strA_alloc++;
   data->strA_total += Memory_Length_3DAs1D( sizeof(StructA), 2, 2, 3 );

   data->groupOneFunc_alloc = data->strA_alloc + data->strB_alloc + data->strC_alloc;
   data->groupOneFunc_total = data->strA_total + data->strB_total + data->strC_total;
}


void MemoryReportSuite_AllocGroupTwo( MemoryReportSuiteData* data ) {
   Index          x1 = 4;
   Index          y1[] = { 1, 2, 3, 4 };
   Index          x2 = 2;
   Index          y2[] = { 1, 1 };

   data->object = Memory_Alloc( StructA, "GroupTwo" );
   data->strA_alloc++;
   data->strA_total += sizeof(StructA);
   data->groupTwo_strA_alloc++;
   data->groupTwo_strA_total += sizeof(StructA);

   data->one4d = Memory_Alloc_4DArrayAs1D( StructB, 4, 2, 3, 5, "GroupTwo" );
   data->strB_alloc++;
   data->strB_total += Memory_Length_4DAs1D( sizeof(StructB), 4, 2, 3, 5 );
   data->groupTwo_strB_alloc++;
   data->groupTwo_strB_total += Memory_Length_4DAs1D( sizeof(StructB), 4, 2, 3, 5 );

   data->complex2d = Memory_Alloc_2DComplex( StructC, x1, y1, "GroupTwo" );
   data->strC_alloc++;
   data->strC_total += Memory_Length_2DComplex( sizeof(StructC), x1, y1 );
   data->groupTwo_strC_alloc++;
   data->groupTwo_strC_total += Memory_Length_2DComplex( sizeof(StructC), x1, y1 );

   data->setup = Memory_Alloc_3DSetup( x2, y2 );
   data->index_alloc++;
   data->index_total += Memory_Length_2DComplex( sizeof(Index), x2, y2 );
   /* The Index allocation won't be classed as coming from this name, as it's called
    *  within Memory_Alloc_3DSetup */   

   data->setup[0][0] = 2;
   data->setup[1][0] = 3;
   data->complex3d = Memory_Alloc_3DComplex( StructA, x2, y2, data->setup, "GroupTwo" );
   data->strA_alloc++;
   data->strA_total += Memory_Length_3DComplex( sizeof(StructA), x2, y2, data->setup );
   data->groupTwo_strA_alloc++;
   data->groupTwo_strA_total += Memory_Length_3DComplex( sizeof(StructA), x2, y2, data->setup );

   data->groupTwoName_alloc = data->groupTwo_strA_alloc + data->groupTwo_strB_alloc + data->groupTwo_strC_alloc;
   data->groupTwoName_total = data->groupTwo_strA_total + data->groupTwo_strB_total + data->groupTwo_strC_total;
   data->groupTwoFunc_alloc = data->groupTwoName_alloc + data->index_alloc;
   data->groupTwoFunc_total = data->groupTwoName_total + data->index_total;
}


void MemoryReportSuite_AllocBytes( MemoryReportSuiteData* data ) {
   data->bytesObj = Memory_Alloc_Bytes( 5, "Bytes", "BytesGroup" );
   data->bytes_alloc++;
   data->bytes_total += 5;

   data->bytesArray = Memory_Alloc_Array_Bytes( 3, 10, "Bytes", "BytesGroup" );
   data->bytes_alloc++;
   data->bytes_total += 3 * 10;

   data->bytesFunc_alloc = data->bytes_alloc;
   data->bytesFunc_total = data->bytes_total;
}


Memory* MemoryReportSuite_SaveStgMemoryAndCreateTemp( MemoryReportSuiteData* data ) {
   /* Save the main stgMemory struct, and create a special one for this test */
   data->savedStgMemory = stgMemory;
   stgMemory = Memory_Init();

   stgMemory->infoStream = Stg_Class_Copy( (Stream*)Journal_GetTypedStream( Info_Type ), NULL, True, NULL, NULL );
   stgMemory->debugStream = Stg_Class_Copy( (Stream*)Journal_GetTypedStream( Debug_Type ), NULL, True, NULL, NULL );
   stgMemory->errorStream = Stg_Class_Copy( (Stream*)Journal_GetTypedStream( Error_Type ), NULL, True, NULL, NULL );
   Journal_Enable_TypedStream( Info_Type, True );
   return stgMemory;
}


void MemoryReportSuite_Setup( MemoryReportSuiteData* data ) {
   MPI_Comm_rank( MPI_COMM_WORLD, &data->rank );

   data->strA_alloc=0;
   data->strB_alloc=0;
   data->strC_alloc=0;
   data->bytes_alloc=0;
   data->index_alloc=0;
   data->strA_total=0;
   data->strB_total=0;
   data->strC_total=0;
   data->bytes_total=0;
   data->index_total=0;
   data->groupTwo_strA_alloc=0;
   data->groupTwo_strB_alloc=0;
   data->groupTwo_strC_alloc=0;
   data->groupTwo_strA_total=0;
   data->groupTwo_strB_total=0;
   data->groupTwo_strC_total=0;
   data->groupOneFunc_alloc=0;
   data->groupTwoFunc_alloc=0;
   data->groupTwoName_alloc=0;
   data->bytesFunc_alloc=0;
   data->groupOneFunc_total=0;
   data->groupTwoFunc_total=0;
   data->groupTwoName_total=0;
   data->bytesFunc_total=0;
   data->file_alloc=0;
   data->file_total=0;
}


void MemoryReportSuite_Teardown( MemoryReportSuiteData* data ) {
}


void MemoryReportSuite_AllocTestMemoryObjects( MemoryReportSuiteData* data ) {
   MemoryReportSuite_AllocGroupOne( data );
   MemoryReportSuite_AllocGroupTwo( data );
   MemoryReportSuite_AllocBytes( data );
   data->file_alloc = data->groupOneFunc_alloc + data->groupTwoFunc_alloc
      + data->bytesFunc_alloc;
   data->file_total = data->groupOneFunc_total + data->groupTwoFunc_total
      + data->bytesFunc_total;
}


void MemoryReportSuite_FreeTestMemoryObjects( MemoryReportSuiteData* data ) {
   Memory_Free( data->bytesObj );
   Memory_Free( data->bytesArray );
   Memory_Free_Type( StructA );
   Memory_Free_Type( StructB );
   Memory_Free_Type( StructC );
}


/*Test 1: MemoryReport: (Type), where file is this test*/
void MemoryReportSuite_TestReportPrintsOne( MemoryReportSuiteData* data ) {
   char*          memoryReportOutputFilename = "./MemoryReportSuite_TestOutput-1.txt";
   Memory*        tempMemoryManager;   

   Stream_RedirectFile( stgMemory->infoStream, memoryReportOutputFilename );

   tempMemoryManager = MemoryReportSuite_SaveStgMemoryAndCreateTemp( data );
   MemoryReportSuite_AllocTestMemoryObjects( data );
   stgMemory = data->savedStgMemory;

   data->report = MemoryReport_New();
   MemoryReport_SetCustomMemoryManager( data->report, tempMemoryManager );
   MemoryReport_AddGroup( data->report, MEMORYREPORT_TYPE );
   MemoryReport_AddCondition( data->report, MEMORYREPORT_FILE, "StGermain/Base/Foundation/tests/MemoryReportSuite.c" );
   MemoryReport_Print( data->report );

   if (data->rank==0) {
      #define        MAXLINE 1000
      FILE*          memoryReportOutputFile = NULL;
      char           memoryReportString[MAXLINE];
      unsigned       timesAlloc, timesFree, currBytes, totalBytes;
      char           valString[1000];

      valString[0] = '\0';

      memoryReportOutputFile = fopen(memoryReportOutputFilename, "r");
      /* Just skip first 2 lines: headings */
      pcu_check_true( fgets( memoryReportString, MAXLINE, memoryReportOutputFile ) );
      pcu_check_true( fgets( memoryReportString, MAXLINE, memoryReportOutputFile ) );

      fscanf( memoryReportOutputFile, "%u %u %u %u\n", &timesAlloc, &timesFree, &currBytes, &totalBytes ); 
      pcu_check_true( timesAlloc == data->file_alloc );
      pcu_check_true( timesFree == 0 );
      pcu_check_true( currBytes == data->file_total );
      pcu_check_true( totalBytes == data->file_total );
      
      /* skip a heading line */
      pcu_check_true( fgets( memoryReportString, MAXLINE, memoryReportOutputFile ) );

      /* Due to sorting, order should be: bytes, index, strA, strB, strC */
      fscanf( memoryReportOutputFile, "%s %u %u %u %u\n", valString, &timesAlloc, &timesFree, &currBytes, &totalBytes ); 
      pcu_check_streq( valString, "Bytes" );
      pcu_check_true( timesAlloc == data->bytes_alloc );
      pcu_check_true( timesFree == 0 );
      pcu_check_true( currBytes == data->bytes_total );
      pcu_check_true( totalBytes == data->bytes_total );
      fscanf( memoryReportOutputFile, "%s %u %u %u %u\n", valString, &timesAlloc, &timesFree, &currBytes, &totalBytes ); 
      pcu_check_streq( valString, "Index" );
      pcu_check_true( timesAlloc == data->index_alloc );
      pcu_check_true( timesFree == 0 );
      pcu_check_true( currBytes == data->index_total );
      pcu_check_true( totalBytes == data->index_total );
      fscanf( memoryReportOutputFile, "%s %u %u %u %u\n", valString, &timesAlloc, &timesFree, &currBytes, &totalBytes ); 
      pcu_check_streq( valString, "StructA" );
      pcu_check_true( timesAlloc == data->strA_alloc );
      pcu_check_true( timesFree == 0 );
      pcu_check_true( currBytes == data->strA_total );
      pcu_check_true( totalBytes == data->strA_total );
      fscanf( memoryReportOutputFile, "%s %u %u %u %u\n", valString, &timesAlloc, &timesFree, &currBytes, &totalBytes ); 
      pcu_check_streq( valString, "StructB" );
      pcu_check_true( timesAlloc == data->strB_alloc );
      pcu_check_true( timesFree == 0 );
      pcu_check_true( currBytes == data->strB_total );
      pcu_check_true( totalBytes == data->strB_total );
      fscanf( memoryReportOutputFile, "%s %u %u %u %u\n", valString, &timesAlloc, &timesFree, &currBytes, &totalBytes ); 
      pcu_check_streq( valString, "StructC" );
      pcu_check_true( timesAlloc == data->strC_alloc );
      pcu_check_true( timesFree == 0 );
      pcu_check_true( currBytes == data->strC_total );
      pcu_check_true( totalBytes == data->strC_total );

      fclose( memoryReportOutputFile );
      remove( memoryReportOutputFilename );
   }

   stgMemory = tempMemoryManager;
   MemoryReportSuite_FreeTestMemoryObjects( data );
   Memory_Delete();
   stgMemory = data->savedStgMemory;

   MemoryReport_Delete( data->report );
}


/*Test 2: MemoryReport: (Type), where name=Test1*/
void MemoryReportSuite_TestReportPrintsTwo( MemoryReportSuiteData* data ) {
   char*          memoryReportOutputFilename = "./MemoryReportSuite_TestOutput-2.txt";
   Memory*        tempMemoryManager;   
   
   Stream_RedirectFile( stgMemory->infoStream, memoryReportOutputFilename );

   tempMemoryManager = MemoryReportSuite_SaveStgMemoryAndCreateTemp( data );
   MemoryReportSuite_AllocTestMemoryObjects( data );
   stgMemory = data->savedStgMemory;

   data->report = MemoryReport_New();
   MemoryReport_SetCustomMemoryManager( data->report, tempMemoryManager );
   MemoryReport_AddGroup( data->report, MEMORYREPORT_TYPE );
   MemoryReport_AddCondition( data->report, MEMORYREPORT_NAME, "GroupTwo" );
   MemoryReport_Print( data->report );

   if (data->rank==0) {
      #define        MAXLINE 1000
      FILE*          memoryReportOutputFile = NULL;
      char           memoryReportString[MAXLINE];
      unsigned       timesAlloc, timesFree, currBytes, totalBytes;
      char           valString[1000];

      valString[0] = '\0';

      memoryReportOutputFile = fopen(memoryReportOutputFilename, "r");
      /* Just skip first 2 lines: headings */
      pcu_check_true( fgets( memoryReportString, MAXLINE, memoryReportOutputFile ) );
      pcu_check_true( fgets( memoryReportString, MAXLINE, memoryReportOutputFile ) );

      fscanf( memoryReportOutputFile, "%u %u %u %u\n", &timesAlloc, &timesFree, &currBytes, &totalBytes ); 
      pcu_check_true( timesAlloc == data->groupTwoName_alloc );
      pcu_check_true( timesFree == 0 );
      pcu_check_true( currBytes == data->groupTwoName_total );
      pcu_check_true( totalBytes == data->groupTwoName_total );
      
      /* skip a heading line */
      pcu_check_true( fgets( memoryReportString, MAXLINE, memoryReportOutputFile ) );

      /* Due to sorting, order should be: index, strA, strB, strC */
      fscanf( memoryReportOutputFile, "%s %u %u %u %u\n", valString, &timesAlloc, &timesFree, &currBytes, &totalBytes ); 
      pcu_check_streq( valString, "StructA" );
      pcu_check_true( timesAlloc == data->groupTwo_strA_alloc );
      pcu_check_true( timesFree == 0 );
      pcu_check_true( currBytes == data->groupTwo_strA_total );
      pcu_check_true( totalBytes == data->groupTwo_strA_total );
      fscanf( memoryReportOutputFile, "%s %u %u %u %u\n", valString, &timesAlloc, &timesFree, &currBytes, &totalBytes ); 
      pcu_check_streq( valString, "StructB" );
      pcu_check_true( timesAlloc == data->groupTwo_strB_alloc );
      pcu_check_true( timesFree == 0 );
      pcu_check_true( currBytes == data->groupTwo_strB_total );
      pcu_check_true( totalBytes == data->groupTwo_strB_total );
      fscanf( memoryReportOutputFile, "%s %u %u %u %u\n", valString, &timesAlloc, &timesFree, &currBytes, &totalBytes ); 
      pcu_check_streq( valString, "StructC" );
      pcu_check_true( timesAlloc == data->groupTwo_strC_alloc );
      pcu_check_true( timesFree == 0 );
      pcu_check_true( currBytes == data->groupTwo_strC_total );
      pcu_check_true( totalBytes == data->groupTwo_strC_total );

      fclose( memoryReportOutputFile );
      remove( memoryReportOutputFilename );
   }

   stgMemory = tempMemoryManager;
   MemoryReportSuite_FreeTestMemoryObjects( data );
   Memory_Delete();
   stgMemory = data->savedStgMemory;

   MemoryReport_Delete( data->report );
}


/*Test 3: MemoryReport: (Func), where file= this file */
void MemoryReportSuite_TestReportPrintsThree( MemoryReportSuiteData* data ) {
   char*          memoryReportOutputFilename = "./MemoryReportSuite_TestOutput-3.txt";
   Memory*        tempMemoryManager;   
   
   Stream_RedirectFile( stgMemory->infoStream, memoryReportOutputFilename );

   tempMemoryManager = MemoryReportSuite_SaveStgMemoryAndCreateTemp( data );
   MemoryReportSuite_AllocTestMemoryObjects( data );
   stgMemory = data->savedStgMemory;

   data->report = MemoryReport_New();
   MemoryReport_SetCustomMemoryManager( data->report, tempMemoryManager );
   MemoryReport_AddGroup( data->report, MEMORYREPORT_FUNC );
   MemoryReport_AddGroup( data->report, MEMORYREPORT_FILE );
   MemoryReport_AddGroup( data->report, MEMORYREPORT_TYPE );
   MemoryReport_AddGroup( data->report, MEMORYREPORT_NAME );
   MemoryReport_AddCondition( data->report, MEMORYREPORT_FILE, "StGermain/Base/Foundation/tests/MemoryReportSuite.c" );
   MemoryReport_Print( data->report );

   if (data->rank==0) {
      #define        MAXLINE 1000
      FILE*          memoryReportOutputFile = NULL;
      char           memoryReportString[MAXLINE];
      unsigned       timesAlloc, timesFree, currBytes, totalBytes;
      char           valString[1000];

      valString[0] = '\0';

      memoryReportOutputFile = fopen(memoryReportOutputFilename, "r");
      /* Just skip first 2 lines: headings */
      pcu_check_true( fgets( memoryReportString, MAXLINE, memoryReportOutputFile ) );
      pcu_check_true( fgets( memoryReportString, MAXLINE, memoryReportOutputFile ) );

      fscanf( memoryReportOutputFile, "%u %u %u %u\n", &timesAlloc, &timesFree, &currBytes, &totalBytes ); 
      pcu_check_true( timesAlloc == data->file_alloc );
      pcu_check_true( timesFree == 0 );
      pcu_check_true( currBytes == data->file_total );
      pcu_check_true( totalBytes == data->file_total );
      
      /* skip a heading line */
      pcu_check_true( fgets( memoryReportString, MAXLINE, memoryReportOutputFile ) );

      /* Due to sorting, order should be: allocBytes, allocGroupOne, allocGroupTwo */
      fscanf( memoryReportOutputFile, "%s %u %u %u %u\n", valString, &timesAlloc, &timesFree, &currBytes, &totalBytes ); 
      pcu_check_streq( valString, "MemoryReportSuite_AllocBytes" );
      pcu_check_true( timesAlloc == data->bytesFunc_alloc );
      pcu_check_true( timesFree == 0 );
      pcu_check_true( currBytes == data->bytesFunc_total );
      pcu_check_true( totalBytes == data->bytesFunc_total );
      fscanf( memoryReportOutputFile, "%s %u %u %u %u\n", valString, &timesAlloc, &timesFree, &currBytes, &totalBytes ); 
      pcu_check_streq( valString, "MemoryReportSuite_AllocGroupOne" );
      pcu_check_true( timesAlloc == data->groupOneFunc_alloc );
      pcu_check_true( timesFree == 0 );
      pcu_check_true( currBytes == data->groupOneFunc_total );
      pcu_check_true( totalBytes == data->groupOneFunc_total );
      fscanf( memoryReportOutputFile, "%s %u %u %u %u\n", valString, &timesAlloc, &timesFree, &currBytes, &totalBytes ); 
      pcu_check_streq( valString, "MemoryReportSuite_AllocGroupTwo" );
      pcu_check_true( timesAlloc == data->groupTwoFunc_alloc );
      pcu_check_true( timesFree == 0 );
      pcu_check_true( currBytes == data->groupTwoFunc_total );
      pcu_check_true( totalBytes == data->groupTwoFunc_total );

      fclose( memoryReportOutputFile );
      remove( memoryReportOutputFilename );
   }

   stgMemory = tempMemoryManager;
   MemoryReportSuite_FreeTestMemoryObjects( data );
   Memory_Delete();
   stgMemory = data->savedStgMemory;

   MemoryReport_Delete( data->report );
}


/*Test 4: MemoryReport: (File), where type=StructA and name=Test2 */
void MemoryReportSuite_TestReportPrintsFour( MemoryReportSuiteData* data ) {
   char*          memoryReportOutputFilename = "./MemoryReportSuite_TestOutput-4.txt";
   Memory*        tempMemoryManager;   
   
   Stream_RedirectFile( stgMemory->infoStream, memoryReportOutputFilename );

   tempMemoryManager = MemoryReportSuite_SaveStgMemoryAndCreateTemp( data );
   MemoryReportSuite_AllocTestMemoryObjects( data );
   stgMemory = data->savedStgMemory;

   data->report = MemoryReport_New();
   MemoryReport_SetCustomMemoryManager( data->report, tempMemoryManager );
   MemoryReport_AddGroup( data->report, MEMORYREPORT_FILE );
   MemoryReport_AddCondition( data->report, MEMORYREPORT_TYPE, "StructA" );
   MemoryReport_AddCondition( data->report, MEMORYREPORT_NAME, "GroupTwo" );
   MemoryReport_Print( data->report );

   if (data->rank==0) {
      #define        MAXLINE 1000
      FILE*          memoryReportOutputFile = NULL;
      char           memoryReportString[MAXLINE];
      unsigned       timesAlloc, timesFree, currBytes, totalBytes;
      char           valString[1000];

      valString[0] = '\0';

      memoryReportOutputFile = fopen(memoryReportOutputFilename, "r");
      /* Just skip first 2 lines: headings */
      pcu_check_true( fgets( memoryReportString, MAXLINE, memoryReportOutputFile ) );
      pcu_check_true( fgets( memoryReportString, MAXLINE, memoryReportOutputFile ) );

      fscanf( memoryReportOutputFile, "%u %u %u %u\n", &timesAlloc, &timesFree, &currBytes, &totalBytes ); 
      pcu_check_true( timesAlloc == data->groupTwo_strA_alloc );
      pcu_check_true( timesFree == 0 );
      pcu_check_true( currBytes == data->groupTwo_strA_total );
      pcu_check_true( totalBytes == data->groupTwo_strA_total );
      
      /* skip a heading line */
      pcu_check_true( fgets( memoryReportString, MAXLINE, memoryReportOutputFile ) );

      /* Due to sorting, order should be: allocBytes, allocGroupOne, allocGroupTwo */
      fscanf( memoryReportOutputFile, "%s %u %u %u %u\n", valString, &timesAlloc, &timesFree, &currBytes, &totalBytes ); 
      pcu_check_streq( valString, "StGermain/Base/Foundation/tests/MemoryReportSuite.c" );
      pcu_check_true( timesAlloc == data->groupTwo_strA_alloc );
      pcu_check_true( timesFree == 0 );
      pcu_check_true( currBytes == data->groupTwo_strA_total );
      pcu_check_true( totalBytes == data->groupTwo_strA_total );

      fclose( memoryReportOutputFile );
      remove( memoryReportOutputFilename );
   }

   stgMemory = tempMemoryManager;
   MemoryReportSuite_FreeTestMemoryObjects( data );
   Memory_Delete();
   stgMemory = data->savedStgMemory;

   MemoryReport_Delete( data->report );
}


void MemoryReportSuite( pcu_suite_t* suite ) {
   pcu_suite_setData( suite, MemoryReportSuiteData );
   pcu_suite_setFixtures( suite, MemoryReportSuite_Setup, MemoryReportSuite_Teardown );
   #ifdef MEMORY_STATS
   pcu_suite_addTest( suite, MemoryReportSuite_TestReportPrintsOne );
   pcu_suite_addTest( suite, MemoryReportSuite_TestReportPrintsTwo );
   pcu_suite_addTest( suite, MemoryReportSuite_TestReportPrintsThree );
   pcu_suite_addTest( suite, MemoryReportSuite_TestReportPrintsFour );
   #endif
}


