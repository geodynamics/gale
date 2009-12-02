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
#include <math.h>

#include "pcu/pcu.h"
#include "StGermain/Base/Foundation/Foundation.h"
#include "StGermain/Base/IO/IO.h"
#include "StGermain/Base/Container/Container.h"
#include "StGermain/Base/Automation/Automation.h"
#include "StGermain/Base/Extensibility/Extensibility.h"
#include "StGermain/Base/Context/Context.h"
#include "VariableDumpStreamSuite.h"

typedef struct {
   Variable_Register*      vr;
} VariableDumpStreamSuiteData;


void VariableDumpStreamSuite_Setup( VariableDumpStreamSuiteData* data ) {
   /* Construction phase --------------------------------------------------------------------------------------------*/
   data->vr = Variable_Register_New();
}


void VariableDumpStreamSuite_Teardown( VariableDumpStreamSuiteData* data ) {
   Variable_Index          var_I;

   /* manually delete all the created Variables */
   for( var_I = 0; var_I < data->vr->count; var_I++ ) {
      Stg_Class_Delete( data->vr->_variable[var_I] );
   }

   Stg_Class_Delete( data->vr );
}


void VariableDumpStreamSuite_TestDump( VariableDumpStreamSuiteData* data ) {
   typedef double Triple[3];
   double*        array;
   Triple*        structArray;
   Index          length = 10;
   Variable*      var;
   Variable*      vec;
   Stream*        dumpStream = NULL;
   Stream*        dumpStream2 = NULL;
   int            i;
   const Name     scalarDumpFilename = "./scalardump.dat" ;
   const Name     vectorDumpFilename = "./vectordump.dat" ;
   FILE*          scalarDumpFile = NULL;
   FILE*          vectorDumpFile = NULL;
   float          scalarCheck = 0;
   float          vectorCheck[3] = {0,0,0};

   array = Memory_Alloc_Array( double, length, "test" );
   structArray = Memory_Alloc_Array( Triple, length, "test" );

   data->vr = Variable_Register_New();

   var = Variable_NewScalar(
      "Scalar",
		NULL,
      Variable_DataType_Double,
      &length,
      NULL,
      (void**)&array,
      data->vr );

   vec = Variable_NewVector(
      "Three",
		NULL,
      Variable_DataType_Double,
      3,
      &length,
      NULL,
      (void**)&structArray,
      data->vr,
      "a",
      "b",
      "c" );

   Variable_Register_BuildAll( data->vr );

   for ( i = 0; i < length; ++i ) {
      Variable_SetValueDouble( var, i, 123.456 );

      Variable_SetValueAtDouble( vec, i, 0, 1.2 );
      Variable_SetValueAtDouble( vec, i, 1, 3.4 );
      Variable_SetValueAtDouble( vec, i, 2, 5.6 );
   }

   dumpStream = Journal_Register( VariableDumpStream_Type, "scalar dump" );
   VariableDumpStream_SetVariable( dumpStream, var, 1, 0, scalarDumpFilename );
   pcu_check_true( Journal_Dump( dumpStream, NULL ) );
   Stream_Flush( dumpStream );

   dumpStream2 = Journal_Register( VariableDumpStream_Type, "vector dump" );
   VariableDumpStream_SetVariable( dumpStream2, vec, 1, 0, vectorDumpFilename );
   pcu_check_true( Journal_Dump( dumpStream2, NULL ) );
   Stream_Flush( dumpStream2 );

   /* Really, you'd think there'd be a function to read a dumpStream. Anyway, hard-code for now */
   /* Don't want to use expected files here, since we need to continually check the writing ability
    * of these dump streams */ 
   scalarDumpFile = fopen( scalarDumpFilename, "r" );
   vectorDumpFile = fopen( vectorDumpFilename, "r" );
   for ( i = 0; i < length; ++i ) {
      /* See VariableDumpStream.c:166 */
      pcu_check_true( fread( &scalarCheck, sizeof(float), 1, scalarDumpFile ));
      pcu_check_true( fabs( scalarCheck - 123.456) < 0.0001 );

      pcu_check_true( fread( vectorCheck, sizeof(float), 3, vectorDumpFile ));
      pcu_check_true( fabs( vectorCheck[0] - 1.2 ) < 0.001 );
      pcu_check_true( fabs( vectorCheck[1] - 3.4 ) < 0.001 );
      pcu_check_true( fabs( vectorCheck[2] - 5.6 ) < 0.001 );
   }

   Memory_Free( array );
   Memory_Free( structArray );

   fclose( scalarDumpFile );
   fclose( vectorDumpFile );
   remove( scalarDumpFilename );
   remove( vectorDumpFilename );
}


void VariableDumpStreamSuite( pcu_suite_t* suite ) {
   pcu_suite_setData( suite, VariableDumpStreamSuiteData );
   pcu_suite_setFixtures( suite, VariableDumpStreamSuite_Setup, VariableDumpStreamSuite_Teardown );
   pcu_suite_addTest( suite, VariableDumpStreamSuite_TestDump );
}


