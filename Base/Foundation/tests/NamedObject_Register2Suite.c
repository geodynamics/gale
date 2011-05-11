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
#include "NamedObject_Register2Suite.h"

typedef struct {
	__Stg_Object
} TestObject2;

Stg_Object* TestObject2_New( Name name ) {
	/* Variables set in this function */
	SizeT                             _sizeOfSelf = sizeof( TestObject2 );
	Type                                     type = "TestObject";
	Stg_Class_DeleteFunction*             _delete = _Stg_Object_Delete;
	Stg_Class_PrintFunction*               _print = _Stg_Object_Print;
	Stg_Class_CopyFunction*                 _copy = _Stg_Object_Copy;
	AllocationType             nameAllocationType = NON_GLOBAL;

	return _Stg_Object_New(  STG_OBJECT_PASSARGS  );
}

typedef struct {
   NamedObject_Register* reg;
   Stg_Object**          testObjects;
   char**                testObjectNames;
   Index                 testObjectsCount;
} NamedObject_RegisterSuite2Data;

void NamedObject_RegisterSuite2_Setup( NamedObject_RegisterSuite2Data* data ) {
   Index  ii=0;
   char   letter='0';

   data->reg = NamedObject_Register_New();
   data->testObjectsCount = 5;
   data->testObjectNames = (char**)malloc(sizeof(char*) * data->testObjectsCount);
   data->testObjects = (Stg_Object**)malloc(sizeof(Stg_Object*) * data->testObjectsCount);

   letter='a';
   for (ii=0; ii < data->testObjectsCount; ii++ ) {
     data->testObjectNames[ii] = (char*)malloc(sizeof(char) * 2 );
      sprintf( data->testObjectNames[ii], "%c", letter );
      letter++;
      data->testObjects[ii] = TestObject2_New( data->testObjectNames[ii] );
   }
}

void NamedObject_RegisterSuite2_Teardown( NamedObject_RegisterSuite2Data* data ) {
   Index ii;

   for (ii=0; ii < data->testObjectsCount; ii++ ) {
      free( data->testObjectNames[ii] );
   }
   free( data->testObjects );
   free( data->testObjectNames );
}

void NamedObject_RegisterSuite2_TestDeleteAll( NamedObject_RegisterSuite2Data* data ) {
   Index ii;

   pcu_docstring( "Tests a series of new objects can be added to a register, and the register's data "
      "fields are updated correctly" );

   for (ii=0; ii < data->testObjectsCount; ii++ ) {
      NamedObject_Register_Add( data->reg, data->testObjects[ii] );
   }

   pcu_check_true( data->testObjectsCount == data->reg->objects->count );
   for (ii=0; ii < data->testObjectsCount; ii++ ) {
      pcu_check_true( data->testObjects[ii] == data->reg->objects->data[ii] );
   }

   NamedObject_Register_DeleteAll( data->reg );
}


void NamedObject_Register2Suite( pcu_suite_t* suite ) {
   pcu_suite_setData( suite, NamedObject_RegisterSuite2Data );
   pcu_suite_setFixtures( suite, NamedObject_RegisterSuite2_Setup, NamedObject_RegisterSuite2_Teardown );
   pcu_suite_addTest( suite, NamedObject_RegisterSuite2_TestDeleteAll );
}


