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
** $Id: testFieldVariable_Register.c 2432 2004-12-16 23:01:59Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>

#include "pcu/pcu.h"
#include "StGermain/Base/Foundation/Foundation.h"
#include "StGermain/Base/Foundation/forwardDecl.h" /* For Journal stuff */
#include "NamedObject_RegisterSuite.h"

typedef struct {
	__Stg_Object
} TestObject;

Stg_Object* TestObject_New( Name name ) {
	return _Stg_Object_New(
		sizeof( TestObject ),
		"TestObject",
		_Stg_Object_Delete,
		_Stg_Object_Print,
		_Stg_Object_Copy, 
		name,
		NON_GLOBAL );
}

typedef struct {
   NamedObject_Register* reg;
   Stg_Object**          testObjects;
   char**                testObjectNames;
   Index                 testObjectsCount;
} NamedObject_RegisterSuiteData;

void NamedObject_RegisterSuite_Setup( NamedObject_RegisterSuiteData* data ) {
   Index  ii=0;
   char   letter='0';

   data->reg = NamedObject_Register_New();
   data->testObjectsCount = 5;
   data->testObjectNames = malloc(sizeof(char*) * data->testObjectsCount);
   data->testObjects = malloc(sizeof(Stg_Object*) * data->testObjectsCount);

   letter='a';
   for (ii=0; ii < data->testObjectsCount; ii++ ) {
      data->testObjectNames[ii] = malloc(sizeof(char) * 2 );
      sprintf( data->testObjectNames[ii], "%c", letter );
      letter++;
      data->testObjects[ii] = TestObject_New( data->testObjectNames[ii] );
   }
}

void NamedObject_RegisterSuite_Teardown( NamedObject_RegisterSuiteData* data ) {
   Index ii;

   _NamedObject_Register_Delete( data->reg );
   for (ii=0; ii < data->testObjectsCount; ii++ ) {
      _Stg_Object_Delete( data->testObjects[ii] );
      free( data->testObjectNames[ii] );
   }
   free( data->testObjects );
   free( data->testObjectNames );
}

void NamedObject_RegisterSuite_TestAdd( NamedObject_RegisterSuiteData* data ) {
   Index ii;

   for (ii=0; ii < data->testObjectsCount; ii++ ) {
      NamedObject_Register_Add( data->reg, data->testObjects[ii] );
   }

   pcu_check_true( data->testObjectsCount == data->reg->objects->count );
   for (ii=0; ii < data->testObjectsCount; ii++ ) {
      pcu_check_true( data->testObjects[ii] == data->reg->objects->data[ii] );
   }
}


void NamedObject_RegisterSuite_TestGetFunctions( NamedObject_RegisterSuiteData* data ) {
   Index ii;

   for (ii=0; ii < data->testObjectsCount; ii++ ) {
      NamedObject_Register_Add( data->reg, data->testObjects[ii] );
   }

   for (ii=0; ii < data->testObjectsCount; ii++ ) {
      pcu_check_true( ii == NamedObject_Register_GetIndex( data->reg,
         data->testObjectNames[ii] ) );
      pcu_check_true( data->testObjects[ii] = NamedObject_Register_GetByIndex( data->reg, ii ) );
      pcu_check_true( data->testObjects[ii] = NamedObject_Register_GetByName( data->reg,
         data->testObjectNames[ii] ) );
   }
}

void NamedObject_RegisterSuite( pcu_suite_t* suite ) {
   pcu_suite_setData( suite, NamedObject_RegisterSuiteData );
   pcu_suite_setFixtures( suite, NamedObject_RegisterSuite_Setup, NamedObject_RegisterSuite_Teardown );
   pcu_suite_addTest( suite, NamedObject_RegisterSuite_TestAdd );
   pcu_suite_addTest( suite, NamedObject_RegisterSuite_TestGetFunctions );
}
