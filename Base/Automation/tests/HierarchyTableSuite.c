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
#include "HierarchyTableSuite.h"

typedef struct {
   HierarchyTable*   hTable;
   HierarchyTable*   savedHierarchyTable;
   unsigned int      rank;
} HierarchyTableSuiteData;


const Type A_Type = "A";
const Type B_Type = "B";
const Type C_Type = "C";
const Type D_Type = "D";
const Type AA_Type = "AA";
const Type BB_Type = "BB";


void HierarchyTableSuite_Setup( HierarchyTableSuiteData* data ) {
   /* We need to operate on the stgHierarchyTable ptr, since HierarchyTable_New() even uses this */
   data->savedHierarchyTable = stgHierarchyTable;
   stgHierarchyTable = NULL;
   data->hTable = HierarchyTable_New();
   stgHierarchyTable = data->savedHierarchyTable;
   HierarchyTable_RegisterParent( data->hTable, B_Type, A_Type );
   HierarchyTable_RegisterParent( data->hTable, C_Type, B_Type );
   HierarchyTable_RegisterParent( data->hTable, D_Type, C_Type );
   HierarchyTable_RegisterParent( data->hTable, BB_Type, AA_Type );
   MPI_Comm_rank( MPI_COMM_WORLD, &data->rank );
}

void HierarchyTableSuite_Teardown( HierarchyTableSuiteData* data ) {
   Stg_Class_Delete( data->hTable );
}
   

void HierarchyTableSuite_TestIsChild( HierarchyTableSuiteData* data ) {
   pcu_check_true( True == HierarchyTable_IsChild( data->hTable, A_Type, A_Type ) );
   pcu_check_true( False == HierarchyTable_IsChild( data->hTable, A_Type, B_Type ) );
   pcu_check_true( False == HierarchyTable_IsChild( data->hTable, A_Type, C_Type ) );
   pcu_check_true( False == HierarchyTable_IsChild( data->hTable, A_Type, D_Type ) );
   pcu_check_true( False == HierarchyTable_IsChild( data->hTable, A_Type, AA_Type ) );
   pcu_check_true( False == HierarchyTable_IsChild( data->hTable, A_Type, BB_Type ) );

   pcu_check_true( True == HierarchyTable_IsChild( data->hTable, B_Type, A_Type ) );
   pcu_check_true( True == HierarchyTable_IsChild( data->hTable, B_Type, B_Type ) );
   pcu_check_true( False == HierarchyTable_IsChild( data->hTable, B_Type, C_Type ) );
   pcu_check_true( False == HierarchyTable_IsChild( data->hTable, B_Type, D_Type ) );
   pcu_check_true( False == HierarchyTable_IsChild( data->hTable, B_Type, AA_Type ) );
   pcu_check_true( False == HierarchyTable_IsChild( data->hTable, B_Type, BB_Type ) );
   
   pcu_check_true( True == HierarchyTable_IsChild( data->hTable, C_Type, A_Type ) );
   pcu_check_true( True == HierarchyTable_IsChild( data->hTable, C_Type, B_Type ) );
   pcu_check_true( True == HierarchyTable_IsChild( data->hTable, C_Type, C_Type ) );
   pcu_check_true( False == HierarchyTable_IsChild( data->hTable, C_Type, D_Type ) );
   pcu_check_true( False == HierarchyTable_IsChild( data->hTable, C_Type, AA_Type ) );
   pcu_check_true( False == HierarchyTable_IsChild( data->hTable, C_Type, BB_Type ) ); 

   pcu_check_true( True == HierarchyTable_IsChild( data->hTable, D_Type, A_Type ) );
   pcu_check_true( True == HierarchyTable_IsChild( data->hTable, D_Type, B_Type ) );
   pcu_check_true( True == HierarchyTable_IsChild( data->hTable, D_Type, C_Type ) );
   pcu_check_true( True == HierarchyTable_IsChild( data->hTable, D_Type, D_Type ) );
   pcu_check_true( False == HierarchyTable_IsChild( data->hTable, D_Type, AA_Type ) );
   pcu_check_true( False == HierarchyTable_IsChild( data->hTable, D_Type, BB_Type ) );

   pcu_check_true( False == HierarchyTable_IsChild( data->hTable, AA_Type, A_Type ) );
   pcu_check_true( False == HierarchyTable_IsChild( data->hTable, AA_Type, B_Type ) );
   pcu_check_true( False == HierarchyTable_IsChild( data->hTable, AA_Type, C_Type ) );
   pcu_check_true( False == HierarchyTable_IsChild( data->hTable, AA_Type, D_Type ) );
   pcu_check_true( True == HierarchyTable_IsChild( data->hTable, AA_Type, AA_Type ) );
   pcu_check_true( False == HierarchyTable_IsChild( data->hTable, AA_Type, BB_Type ) ); 

   pcu_check_true( False == HierarchyTable_IsChild( data->hTable, BB_Type, A_Type ) );
   pcu_check_true( False == HierarchyTable_IsChild( data->hTable, BB_Type, B_Type ) );
   pcu_check_true( False == HierarchyTable_IsChild( data->hTable, BB_Type, C_Type ) );
   pcu_check_true( False == HierarchyTable_IsChild( data->hTable, BB_Type, D_Type ) );
   pcu_check_true( True == HierarchyTable_IsChild( data->hTable, BB_Type, AA_Type ) );
   pcu_check_true( True == HierarchyTable_IsChild( data->hTable, BB_Type, BB_Type ) ); 
}


void HierarchyTableSuite_TestPrintParents( HierarchyTableSuiteData* data ) {
   Stream*           stream = Journal_Register( Info_Type, "testStream" );
   const char* const testFilename = "testHTable-PrintParents.txt";

   if (data->rank==0) {
      Stream_RedirectFile( stream, testFilename );
   }

   HierarchyTable_PrintParents( data->hTable, A_Type, stream );
   HierarchyTable_PrintParents( data->hTable, B_Type, stream );
   HierarchyTable_PrintParents( data->hTable, C_Type, stream );
   HierarchyTable_PrintParents( data->hTable, D_Type, stream );
   HierarchyTable_PrintParents( data->hTable, AA_Type, stream );
   HierarchyTable_PrintParents( data->hTable, BB_Type, stream );

   /* Just do these print tests with rank 0, to avoid I/O problems */
   if (data->rank==0) {
      char  expectedFilename[PCU_PATH_MAX];

      pcu_filename_expected( testFilename, expectedFilename );
      pcu_check_fileEq( testFilename, expectedFilename );
      remove( testFilename );
   }
}
   

void HierarchyTableSuite_TestPrintChildren( HierarchyTableSuiteData* data ) {
   Stream*     stream = Journal_Register( Info_Type, "testStream" );
   const char* const testFilename = "testHTable-PrintChildren.txt";

   if (data->rank==0) {
      Stream_RedirectFile( stream, testFilename );
   }

   HierarchyTable_PrintChildren( data->hTable, A_Type, stream );
   HierarchyTable_PrintChildren( data->hTable, B_Type, stream );
   HierarchyTable_PrintChildren( data->hTable, C_Type, stream );
   HierarchyTable_PrintChildren( data->hTable, D_Type, stream );
   HierarchyTable_PrintChildren( data->hTable, AA_Type, stream );
   HierarchyTable_PrintChildren( data->hTable, BB_Type, stream );

   /* Just do these print tests with rank 0, to avoid I/O problems */
   if (data->rank==0) {
      char  expectedFilename[PCU_PATH_MAX];

      pcu_filename_expected( testFilename, expectedFilename );
      pcu_check_fileEq( testFilename, expectedFilename );
      remove( testFilename );
   }
}


void HierarchyTableSuite( pcu_suite_t* suite ) {
   pcu_suite_setData( suite, HierarchyTableSuiteData );
   pcu_suite_setFixtures( suite, HierarchyTableSuite_Setup, HierarchyTableSuite_Teardown );
   pcu_suite_addTest( suite, HierarchyTableSuite_TestIsChild );
   pcu_suite_addTest( suite, HierarchyTableSuite_TestPrintParents );
   pcu_suite_addTest( suite, HierarchyTableSuite_TestPrintChildren );
}


