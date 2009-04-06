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
** $Id: testNamedStg_ObjectList.c 2432 2005-08-08 23:01:59Z Raquibul Hassan $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <unistd.h>

#include "pcu/pcu.h"
#include "StGermain/Base/Foundation/Foundation.h"
#include "StGermain/Base/Foundation/forwardDecl.h" /* For Journal stuff */
#include "PrimitiveObjectSuite.h"


typedef struct {
} PrimitiveObjectSuiteData;

void PrimitiveObjectSuite_Setup( PrimitiveObjectSuiteData* data ) {
}

void PrimitiveObjectSuite_Teardown( PrimitiveObjectSuiteData* data ) {
}


/* Note that since the key feature of the PrimitiveObject component that was tested was really
 * the print capacity, I've had to use the journal outputting to file approach again */
void PrimitiveObjectSuite_TestPrimObjects( PrimitiveObjectSuiteData* data ) {
   Stg_ObjectList*      list;
   Stg_PrimitiveObject* primObject;

   list = Stg_ObjectList_New();

   Stg_ObjectList_Append( list, Stg_PrimitiveObject_New_UnsignedChar( 'a', "char item" ) );
   Stg_ObjectList_Append( list, Stg_PrimitiveObject_New_UnsignedShort( 123, "short item" ) );
   Stg_ObjectList_Append( list, Stg_PrimitiveObject_New_UnsignedInt( 456, "int item" ) );
   Stg_ObjectList_Append( list, Stg_PrimitiveObject_New_UnsignedLong( 789, "long item" ) );
   Stg_ObjectList_Append( list, Stg_PrimitiveObject_New_Char( 'a', "char item" ) );
   Stg_ObjectList_Append( list, Stg_PrimitiveObject_New_Short( -123, "short item" ) );
   Stg_ObjectList_Append( list, Stg_PrimitiveObject_New_Int( -456, "int item" ) );
   Stg_ObjectList_Append( list, Stg_PrimitiveObject_New_Long( -789, "long item" ) );
   Stg_ObjectList_Append( list, Stg_PrimitiveObject_New_Float( 1.2f, "float item" ) );
   Stg_ObjectList_Append( list, Stg_PrimitiveObject_New_Double( 2.4, "double item" ) );

   primObject = (Stg_PrimitiveObject*)Stg_ObjectList_At( list, 0 );
   pcu_assert_true( primObject->dataType == Stg_C_Primitive_Type_UnsignedChar );
   pcu_assert_true( primObject->value.asUnsignedChar == 'a' );
   primObject = (Stg_PrimitiveObject*)Stg_ObjectList_At( list, 1 );
   pcu_assert_true( primObject->dataType == Stg_C_Primitive_Type_UnsignedShort );
   pcu_assert_true( primObject->value.asUnsignedShort == 123 );
   primObject = (Stg_PrimitiveObject*)Stg_ObjectList_At( list, 2 );
   pcu_assert_true( primObject->dataType == Stg_C_Primitive_Type_UnsignedInt );
   pcu_assert_true( primObject->value.asUnsignedInt == 456 );
   primObject = (Stg_PrimitiveObject*)Stg_ObjectList_At( list, 3 );
   pcu_assert_true( primObject->dataType == Stg_C_Primitive_Type_UnsignedLong );
   pcu_assert_true( primObject->value.asUnsignedLong == 789 );
   primObject = (Stg_PrimitiveObject*)Stg_ObjectList_At( list, 4 );
   pcu_assert_true( primObject->dataType == Stg_C_Primitive_Type_Char );
   pcu_assert_true( primObject->value.asChar == 'a' );
   primObject = (Stg_PrimitiveObject*)Stg_ObjectList_At( list, 5 );
   pcu_assert_true( primObject->dataType == Stg_C_Primitive_Type_Short );
   pcu_assert_true( primObject->value.asShort == -123 );
   primObject = (Stg_PrimitiveObject*)Stg_ObjectList_At( list, 6 );
   pcu_assert_true( primObject->dataType == Stg_C_Primitive_Type_Int );
   pcu_assert_true( primObject->value.asInt == -456 );
   primObject = (Stg_PrimitiveObject*)Stg_ObjectList_At( list, 7 );
   pcu_assert_true( primObject->dataType == Stg_C_Primitive_Type_Long );
   pcu_assert_true( primObject->value.asLong == -789 );
   primObject = (Stg_PrimitiveObject*)Stg_ObjectList_At( list, 8 );
   pcu_assert_true( primObject->dataType == Stg_C_Primitive_Type_Float );
   pcu_assert_true( primObject->value.asFloat == 1.2f );
   primObject = (Stg_PrimitiveObject*)Stg_ObjectList_At( list, 9 );
   pcu_assert_true( primObject->dataType == Stg_C_Primitive_Type_Double );
   pcu_assert_true( primObject->value.asDouble == 2.4 );

   Stg_Class_Delete( list );
}


void PrimitiveObjectSuite( pcu_suite_t* suite ) {
   pcu_suite_setData( suite, PrimitiveObjectSuiteData );
   pcu_suite_setFixtures( suite, PrimitiveObjectSuite_Setup, PrimitiveObjectSuite_Teardown );
   pcu_suite_addTest( suite, PrimitiveObjectSuite_TestPrimObjects );
}
