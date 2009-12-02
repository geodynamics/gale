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
#include "Variable_RegisterSuite.h"

typedef struct {
   Variable_Register*	reg;
} Variable_RegisterSuiteData;


void Variable_RegisterSuite_Setup( Variable_RegisterSuiteData* data ) {
   data->reg = Variable_Register_New();
}

void Variable_RegisterSuite_Teardown( Variable_RegisterSuiteData* data ) {
   Stg_Class_Delete( data->reg );
}
   

void Variable_RegisterSuite_TestAddGet( Variable_RegisterSuiteData* data ) {
   Variable*		var[10];
   #define ARRAY_SIZE	4
   #define STRUCT_SIZE	4
   double			array[ARRAY_SIZE];
   Index			   arraySize = ARRAY_SIZE;
   char*			   name[10] = {"testVar0", "testVar1", "testVar2", "testVar3",
                  "testVar4", "testVar5", "testVar6", "testVar7",
                  "testVar8", "testVar9"};
   Index		   	i;

   for (i = 0; i < 10; i++) {
      var[i] = Variable_NewVector( name[i], NULL, Variable_DataType_Double, 4, &arraySize, NULL, (void**)&array, 0 );
   }

   for (i = 0; i < 10; i++)
   {
      Variable_Register_Add(data->reg, var[i]);
   }

   for (i = 0; i < 10; i++) {
      pcu_check_true( i == Variable_Register_GetIndex(data->reg, name[i]));
   }

   for (i = 0; i < 10; i++) {
      pcu_check_true( var[i] == Variable_Register_GetByName(data->reg, name[i]));
   }

   for (i = 0; i < 10; i++) {
      Stg_Class_Delete(var[i]);
   }
}


void Variable_RegisterSuite( pcu_suite_t* suite ) {
   pcu_suite_setData( suite, Variable_RegisterSuiteData );
   pcu_suite_setFixtures( suite, Variable_RegisterSuite_Setup, Variable_RegisterSuite_Teardown );
   pcu_suite_addTest( suite, Variable_RegisterSuite_TestAddGet );
}


