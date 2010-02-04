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
#include "ConditionFunctionSuite.h"

typedef struct {
} ConditionFunctionSuiteData;


#define TEST_CF_RESULT 10

void func(Index index, Variable_Index var_I, void* context, void* result)
{
   *((double*)result) = TEST_CF_RESULT;
}


void ConditionFunctionSuite_Setup( ConditionFunctionSuiteData* data ) {
}

void ConditionFunctionSuite_Teardown( ConditionFunctionSuiteData* data ) {
}
   

void ConditionFunctionSuite_TestApply( ConditionFunctionSuiteData* data ) {
   ConditionFunction*	cf;
   double               result;

   cf = ConditionFunction_New( func, (Name)"quadratic" );

   ConditionFunction_Apply(cf, 4, 2, NULL, &result);
   pcu_check_true( TEST_CF_RESULT == result );

   Stg_Class_Delete(cf);
}


void ConditionFunctionSuite( pcu_suite_t* suite ) {
   pcu_suite_setData( suite, ConditionFunctionSuiteData );
   pcu_suite_setFixtures( suite, ConditionFunctionSuite_Setup, ConditionFunctionSuite_Teardown );
   pcu_suite_addTest( suite, ConditionFunctionSuite_TestApply );
}


