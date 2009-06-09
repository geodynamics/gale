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
#include <math.h>

#include "pcu/pcu.h"
#include "StGermain/Base/Foundation/Foundation.h"
#include "StGermain/Base/Foundation/forwardDecl.h" /* For Journal stuff */
#include "MemMonitorSuite.h"


typedef struct {
   Stg_MemMonitor*   mm;
   unsigned int      rank;   
} MemMonitorSuiteData;


void MemMonitorSuite_Setup( MemMonitorSuiteData* data ) {
   Journal_Enable_TypedStream( Info_Type, True );
   Stream_Enable( Journal_Register( Info_Type, Stg_MemMonitor_InfoStreamName ), True );

   data->mm = NULL;
   MPI_Comm_rank( MPI_COMM_WORLD, &data->rank );
}


void MemMonitorSuite_Teardown( MemMonitorSuiteData* data ) {
   Stg_MemMonitor_Delete( data->mm );
}


void MemMonitorSuite_TestMonitor( MemMonitorSuiteData* data ) {
   char*          memoryReportOutputFilename = "./MemMonitorSuite_TestOutput.txt";
   char*          a;
   char*          b;
   char*          c;
   char*          d;
   char*          e;
   char*          f;
   MemMonitorData mmData;
   int            totalMemAtTestStart;
   int            expMemDiff;
   double         expPercentChange;
   Bool           expCritResult = False;
   
   Stg_MemMonitor_SetMemoryWatchCriteria( 0.2 );
   if (data->rank==0) {
      Stream_RedirectFile( Journal_Register( Info_Type, Stg_MemMonitor_InfoStreamName ), memoryReportOutputFilename );
   }

   /* Don't create the MM until now, so we can control the total memory for testing purposes */
   data->mm = Stg_MemMonitor_New( "test", True, True, MPI_COMM_WORLD );
   a = Memory_Alloc_Array( char, 1024, "a" );

   MemoryField_UpdateAsSumOfSubFields( stgMemory->types );
   totalMemAtTestStart = stgMemory->types->currentAllocation;
   Stg_MemMonitor_Begin( data->mm );
   b = Memory_Alloc_Array( char, 1024*2, "b" );
   Stg_MemMonitor_End( data->mm, &mmData );

   expMemDiff = 1024*2;
   pcu_check_true( mmData.memDiff == expMemDiff );
   pcu_check_true( mmData.memFinal == totalMemAtTestStart + expMemDiff );
   pcu_check_true( mmData.avgProcMemDiff == mmData.memDiff );
   expPercentChange = expMemDiff / (double)totalMemAtTestStart*100; 
   pcu_check_true( fabs( expPercentChange - mmData.percentChange ) < 0.1 );
   expCritResult = fabs(expPercentChange/100) >= 0.2;
   pcu_check_true( mmData.criterionPassed == expCritResult );

   MemoryField_UpdateAsSumOfSubFields( stgMemory->types );
   totalMemAtTestStart = stgMemory->types->currentAllocation;
   Stg_MemMonitor_Begin( data->mm );
   c = Memory_Alloc_Array( char, 100, "c" );
   d = Memory_Alloc_Array( char, 100*1024, "d" );
   Stg_MemMonitor_End( data->mm, &mmData );
   
   expMemDiff = 100 + 100*1024;
   pcu_check_true( mmData.memDiff == expMemDiff );
   pcu_check_true( mmData.memFinal == totalMemAtTestStart + expMemDiff );
   pcu_check_true( mmData.avgProcMemDiff == mmData.memDiff );
   expPercentChange = expMemDiff / (double)totalMemAtTestStart*100; 
   pcu_check_true( fabs( expPercentChange - mmData.percentChange ) < 0.1 );
   expCritResult = fabs(expPercentChange/100) >= 0.2;
   pcu_check_true( mmData.criterionPassed == expCritResult );

   MemoryField_UpdateAsSumOfSubFields( stgMemory->types );
   totalMemAtTestStart = stgMemory->types->currentAllocation;
   Stg_MemMonitor_Begin( data->mm );
   Memory_Free( a );
   Memory_Free( c );
   Memory_Free( d );
   Stg_MemMonitor_End( data->mm, &mmData );

   expMemDiff = -1024 - 100 - 100*1024;
   pcu_check_true( mmData.memDiff == expMemDiff );
   pcu_check_true( mmData.memFinal == totalMemAtTestStart + expMemDiff );
   pcu_check_true( mmData.avgProcMemDiff == mmData.memDiff );
   /* Percent should be negative this time */
   expPercentChange = expMemDiff / (double)totalMemAtTestStart*100; 
   pcu_check_true( fabs( expPercentChange - mmData.percentChange ) < 0.1 );
   expCritResult = fabs(expPercentChange/100) >= 0.2;
   pcu_check_true( mmData.criterionPassed == expCritResult );

   MemoryField_UpdateAsSumOfSubFields( stgMemory->types );
   totalMemAtTestStart = stgMemory->types->currentAllocation;
   Stg_MemMonitor_Begin( data->mm );
   e = Memory_Alloc_Array( char, 10*1024, "e" );
   f = Memory_Alloc_Array( char, 10*1024, "f" );
   Stg_MemMonitor_End( data->mm, &mmData );
   
   expMemDiff = 10*1024 + 10*1024;
   pcu_check_true( mmData.memDiff == expMemDiff );
   pcu_check_true( mmData.memFinal == totalMemAtTestStart + expMemDiff );
   pcu_check_true( mmData.avgProcMemDiff == mmData.memDiff );
   expPercentChange = expMemDiff / (double)totalMemAtTestStart*100; 
   pcu_check_true( fabs( expPercentChange - mmData.percentChange ) < 0.1 );
   expCritResult = fabs(expPercentChange/100) >= 0.2;
   pcu_check_true( mmData.criterionPassed == expCritResult );

   Memory_Free( b );
   Memory_Free( e );
   Memory_Free( f );
   
   if (data->rank==0) {
      remove( memoryReportOutputFilename );
   }
}
 

void MemMonitorSuite( pcu_suite_t* suite ) {
   pcu_suite_setData( suite, MemMonitorSuiteData );
   pcu_suite_setFixtures( suite, MemMonitorSuite_Setup, MemMonitorSuite_Teardown );
   #ifdef MEMORY_STATS
   pcu_suite_addTest( suite, MemMonitorSuite_TestMonitor );
   #endif
}
