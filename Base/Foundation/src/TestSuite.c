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
** $Id: TestSuite.c 3584 2006-05-16 11:11:07Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#include "types.h"
#include "shortcuts.h"
#include "forwardDecl.h"
#include "debug.h"
#include "MemoryTag.h"
#include "Memory.h"
#include "Class.h"
#include "TestSuite.h"


/* Textual name of this class */
const Type TestSuite_Type = "TestSuite";


/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/

TestSuite* TestSuite_New() {
	/* Variables set in this function */
	SizeT                      _sizeOfSelf = sizeof(TestSuite);
	Type                              type = TestSuite_Type;
	Stg_Class_DeleteFunction*      _delete = _TestSuite_Delete;
	Stg_Class_PrintFunction*        _print = _TestSuite_Print;
	Stg_Class_CopyFunction*          _copy = _TestSuite_Copy;

	return _TestSuite_New(  TESTSUITE_PASSARGS  );
}

TestSuite* _TestSuite_New(  TESTSUITE_DEFARGS  ) {
	TestSuite* self;

	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(TestSuite) );
	self = (TestSuite*)_Stg_Class_New(  STG_CLASS_PASSARGS  );

	/* Virtual info */

	/* TestSuite info */
	_TestSuite_Init( self );

	return self;
}

void _TestSuite_Init( TestSuite* self ) {
	assert( self );

	insist( MPI_Comm_size( MPI_COMM_WORLD, &self->nProcs ), == MPI_SUCCESS );
	insist( MPI_Comm_rank( MPI_COMM_WORLD, &self->rank ), == MPI_SUCCESS );
	self->nTests = 0;
	self->tests = NULL;
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _TestSuite_Delete( void* testSuite ) {
	TestSuite*	self = (TestSuite*)testSuite;

	assert( self );

	TestSuite_Destruct( self );

	/* Delete the parent. */
	_Stg_Class_Delete( self );
}

void _TestSuite_Print( void* testSuite, struct Stream* stream ) {
	TestSuite*	self = (TestSuite*)testSuite;

	/* Print parent */
	Journal_Printf( stream, "TestSuite (ptr): (%p)\n", self );
	_Stg_Class_Print( self, stream );
}

void* _TestSuite_Copy( void* testSuite, void* destProc_I, Bool deep, Name nameExt, struct PtrMap* ptrMap ) {
#if 0
	TestSuite*	self = (TestSuite*)testSuite;
	TestSuite*	newTestSuite;
	PtrMap*	map = ptrMap;
	Bool	ownMap = False;

	/* Damn me for making copying so difficult... what was I thinking? */
	
	/* We need to create a map if it doesn't already exist. */
	if( !map ) {
		map = PtrMap_New( 10 );
		ownMap = True;
	}
	
	newTestSuite = (TestSuite*)_Mesh_Copy( self, destProc_I, deep, nameExt, map );
	
	/* Copy the virtual methods here. */

	/* Deep or shallow? */
	if( deep ) {
	}
	else {
	}
	
	/* If we own the map, get rid of it here. */
	if( ownMap ) Stg_Class_Delete( map );
	
	return (void*)newTestSuite;
#endif

	return NULL;
}


/*--------------------------------------------------------------------------------------------------------------------------
** Public Functions
*/

void TestSuite_SetTests( void* testSuite, unsigned nTests, TestSuite_Test* tests ) {
	TestSuite*	self = (TestSuite*)testSuite;

	assert( self );
	assert( !nTests || tests );

	TestSuite_Destruct( self );

	self->nTests = nTests;
	if( nTests ) {
		self->tests = Memory_Alloc_Array( TestSuite_Test, nTests, "TestSuite::tests" );
		memcpy( self->tests, tests, nTests * sizeof(TestSuite_Test) );
	}
}

void TestSuite_Run( void* testSuite ) {
	TestSuite*	self = (TestSuite*)testSuite;
	unsigned	t_i;

	assert( self );

	for( t_i = 0; t_i < self->nTests; t_i++ ) {
		TestSuite_Test*	test = self->tests + t_i;
		Bool		result;

		assert( test );
		assert( test->name );
		assert( test->func );

		if( !self->rank )
			printf( "   Running test '%s'... ", test->name );
		result = test->func( self );
		if( !self->rank )
			printf( "%s\n", result ? "passed" : "failed" );
		insist( MPI_Barrier( MPI_COMM_WORLD ), == MPI_SUCCESS );
	}
}

/*----------------------------------------------------------------------------------------------------------------------------------
** Private Functions
*/

void TestSuite_Destruct( TestSuite* self ) {
	assert( self );

	self->nTests = 0;
	KillArray( self->tests );
}


