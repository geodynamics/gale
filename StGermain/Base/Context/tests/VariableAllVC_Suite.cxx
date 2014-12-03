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
#include <string.h>

#include "pcu/pcu.h"
#include "StGermain/Base/Foundation/Foundation.h"
#include "StGermain/Base/IO/IO.h"
#include "StGermain/Base/Container/Container.h"
#include "StGermain/Base/Automation/Automation.h"
#include "StGermain/Base/Extensibility/Extensibility.h"
#include "StGermain/Base/Context/Context.h"
#include "VariableAllVC_Suite.h"

#define TEST_CONDVALUE 2.0


typedef struct {
   Variable_Register*            vr;
   ConditionFunction_Register*   conFunc_Register;
   Dictionary*                   dict;
   double*                       testArray;
   Index                         arraySize;
   Variable*                     var;
   char*                         vcKey;
   VariableCondition*            vc;
} VariableAllVC_SuiteData;


void VariableAllVC_Suite_CreateDictEntries( VariableAllVC_SuiteData* data ) {
   Dictionary_Entry_Value* info;
   Dictionary_Entry_Value* varList;
   Dictionary_Entry_Value* varValue;

   info = Dictionary_Entry_Value_NewStruct();
   varList = Dictionary_Entry_Value_NewList();
   varValue = Dictionary_Entry_Value_NewStruct();

   Dictionary_Entry_Value_AddMember( varValue, (Dictionary_Entry_Key)"name", Dictionary_Entry_Value_FromString( "test" )  );
   Dictionary_Entry_Value_AddMember( varValue, (Dictionary_Entry_Key)"type", Dictionary_Entry_Value_FromString( "double" )  );
   Dictionary_Entry_Value_AddMember( varValue, (Dictionary_Entry_Key)"value", Dictionary_Entry_Value_FromDouble( TEST_CONDVALUE )  );

   Dictionary_Entry_Value_AddElement( varList, varValue );
   Dictionary_Entry_Value_AddMember( info, (Dictionary_Entry_Key)"variables", varList  );
   Dictionary_Add( data->dict, (Dictionary_Entry_Key)data->vcKey, info );
}


void VariableAllVC_Suite_Setup( VariableAllVC_SuiteData* data  ) {
   data->arraySize = 10;
   data->testArray = Memory_Alloc_Array( double, data->arraySize, "test" );

   data->dict = Dictionary_New();
   data->conFunc_Register = ConditionFunction_Register_New();
   data->vr = Variable_Register_New();
   
   Stg_asprintf( &data->vcKey, "VariableAllVC" );
   VariableAllVC_Suite_CreateDictEntries( data );

   data->var = Variable_NewScalar( "test", NULL, Variable_DataType_Double, &data->arraySize, NULL, (void**)&data->testArray, data->vr );
      
   Variable_Register_BuildAll(data->vr);
   
   data->vc = (VariableCondition*) VariableAllVC_New( "variableAllVC", NULL, "VariableAllVC", data->vr, data->conFunc_Register, data->dict, NULL );
   
   Stg_Component_Build( data->vc, 0, False );
}


void VariableAllVC_Suite_Teardown( VariableAllVC_SuiteData* data ) {
   Stg_Class_Delete(data->vr);
   Stg_Class_Delete(data->conFunc_Register);
   Stg_Class_Delete(data->dict);

   Memory_Free( data->testArray );
   Stg_Class_Delete( data->var );

   _Stg_Component_Delete(data->vc);
}


void VariableAllVC_Suite_TestIsCondition( VariableAllVC_SuiteData* data ) {
   Index node_I=0;

   for (node_I = 0; node_I < data->arraySize; node_I++) {
      pcu_check_true( True == VariableCondition_IsCondition(data->vc, node_I, 0) );
   }
}


void VariableAllVC_Suite_TestGetValueIndex( VariableAllVC_SuiteData* data ) {
   Index node_I=0;

   for (node_I = 0; node_I < data->arraySize; node_I++) {
      pcu_check_true( 0 == VariableCondition_GetValueIndex(data->vc, node_I, 0) );
   }
}


void VariableAllVC_Suite_TestApply( VariableAllVC_SuiteData* data ) {
   Index node_I=0;
   
   VariableCondition_Apply(data->vc, NULL);

   for (node_I = 0; node_I < data->arraySize; node_I++) {
      pcu_check_true( TEST_CONDVALUE == data->testArray[node_I] );
   }
}


void VariableAllVC_Suite( pcu_suite_t* suite ) {
   pcu_suite_setData( suite, VariableAllVC_SuiteData );
   pcu_suite_setFixtures( suite, VariableAllVC_Suite_Setup, VariableAllVC_Suite_Teardown );
   pcu_suite_addTest( suite, VariableAllVC_Suite_TestIsCondition );
   pcu_suite_addTest( suite, VariableAllVC_Suite_TestGetValueIndex );
   pcu_suite_addTest( suite, VariableAllVC_Suite_TestApply );
}


