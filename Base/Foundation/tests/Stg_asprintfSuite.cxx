/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street, Melbourne, 3053, Australia.
**
** Authors:
**	Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**	Patrick D. Sunter, Software Engineer, VPAC. (pds@vpac.org)
**	Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**	Siew-Ching Tan, Software Engineer, VPAC. (siew@vpac.org)
**	Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
**	Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
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
#include <string.h>

#include "pcu/pcu.h"
#include "StGermain/Base/Foundation/Foundation.h"
#include "StGermain/Base/Foundation/forwardDecl.h" /* For Journal stuff */
#include "Stg_asprintfSuite.h"

typedef struct {
} Stg_asprintfSuiteData;

void Stg_asprintfSuite_Setup( Stg_asprintfSuiteData* data ) {
}

void Stg_asprintfSuite_Teardown( Stg_asprintfSuiteData* data ) {
}

void Stg_asprintfSuite_TestPrint( Stg_asprintfSuiteData* data ) {
	Name  fiftyBytes = "01234567890123456789012345678901234567890123456789";
	char*        testString;
	char*        testStringPtr;
   unsigned int offset=0;

	/* Stress testing Stg_asprintf beyond the default alloc number of bytes */
	Stg_asprintf( &testString, "%s%s%s%s", fiftyBytes, fiftyBytes, fiftyBytes, fiftyBytes );
	pcu_check_true( 200 == strlen( testString ) );
   for ( offset=0; offset < 200; offset+=50 ) {
      testStringPtr = testString + offset;
      pcu_check_true( 0 == strncmp( fiftyBytes, testStringPtr, 50 ) );
   }
	Memory_Free( testString );
}

void Stg_asprintfSuite( pcu_suite_t* suite ) {
   pcu_suite_setData( suite, Stg_asprintfSuiteData );
   pcu_suite_setFixtures( suite, Stg_asprintfSuite_Setup, Stg_asprintfSuite_Teardown );
   pcu_suite_addTest( suite, Stg_asprintfSuite_TestPrint );
}


