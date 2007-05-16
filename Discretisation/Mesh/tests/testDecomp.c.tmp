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
** $Id: testRangeSet.c 2136 2004-09-30 02:47:13Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <mpi.h>

#include "Base/Base.h"
#include "Discretisation/Geometry/Geometry.h"
#include "Discretisation/Mesh/Mesh.h"


Bool testSetEmpty( unsigned rank, unsigned nProcs, unsigned watch ) {
	Bool		result = True;
	Decomp*		decomp = Decomp_New( "" );
	unsigned	nLocals = 0;
	unsigned*	locals = NULL;;

	decomp =  Decomp_New( "" );
	Decomp_SetLocals( decomp, nLocals, locals );

	if( rank == watch ) {
		if( Decomp_GetGlobalSize( decomp ) != 0 || 
		    Decomp_GetLocalSize( decomp ) != 0 ||
		    decomp->locals != NULL )
		{
			result = False;
			goto done;
		}
	}

done:
	FreeObject( decomp );

	return result;
}

Bool testSetLocals( unsigned rank, unsigned nProcs, unsigned watch ) {
	Bool		result = True;
	Decomp*		decomp = Decomp_New( "" );
	unsigned	nLocals = 10;
	unsigned	nGlobals = nProcs * nLocals;
	unsigned*	locals;
	unsigned	l_i;

	locals = Memory_Alloc_Array_Unnamed( unsigned, nLocals );
	for( l_i = 0; l_i < nLocals; l_i++ )
		locals[l_i] = rank * nLocals + l_i;

	Decomp_SetLocals( decomp, nLocals, locals );

	if( rank == watch ) {
		unsigned	g_i;

		if( Decomp_GetGlobalSize( decomp ) != nGlobals || 
		    Decomp_GetLocalSize( decomp ) != nLocals || 
		    decomp->locals == NULL )
		{
			result = False;
			goto done;
		}

		for( l_i = 0; l_i < nLocals; l_i++ ) {
			if( Decomp_LocalToGlobal( decomp, l_i ) != locals[l_i] ) {
				result = False;
				goto done;
			}
		}

		for( g_i = 0; g_i < nGlobals; g_i++ ) {
			Bool		mapRes;
			unsigned	local;

			mapRes = Decomp_GlobalToLocal( decomp, g_i, &local );
			if( g_i >= rank * nLocals && g_i < (rank + 1) * nLocals && 
			    (!mapRes || local != g_i - rank * nLocals) )
			{
				result = False;
				goto done;
			}
		}
	}

done:
	FreeArray( locals );
	FreeObject( decomp );

	return result;
}


#define nTests	2

TestSuite*	suite;
TestSuite_Test	tests[nTests] = {{"set empty", testSetEmpty, 10}, 
				 {"set locals", testSetLocals, 10}};


int main( int argc, char* argv[] ) {
	/* Initialise MPI, get world info. */
	MPI_Init( &argc, &argv );

	/* Initialise StGermain. */
	BaseFoundation_Init( &argc, &argv );
	BaseIO_Init( &argc, &argv );
	BaseContainer_Init( &argc, &argv );

	/* Create the test suite. */
	suite = TestSuite_New();
	TestSuite_SetProcToWatch( suite, (argc >= 2) ? atoi( argv[1] ) : 0 );
	TestSuite_SetTests( suite, nTests, tests );

	/* Run the tests. */
	TestSuite_Run( suite );

	/* Destroy test suites. */
	FreeObject( suite );

	/* Finalise StGermain. */
	BaseContainer_Finalise();
	BaseIO_Finalise();
	BaseFoundation_Finalise();

	/* Close off MPI */
	MPI_Finalize();

	return MPI_SUCCESS;
}
