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
** $Id: testNumberHasher.c 2136 2004-09-30 02:47:13Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <mpi.h>

#include "Base/Foundation/Foundation.h"
#include "Base/IO/IO.h"
#include "Base/Container/Container.h"


void genSet( unsigned maxItems, unsigned range, unsigned* set ) {
	unsigned	curItem;

	srand( time( NULL ) );

	curItem = 0;
	while( curItem < maxItems ) {
		unsigned	item;
		unsigned	itm_j;

		item = (unsigned)(((double)rand() / (double)RAND_MAX) * (double)range);
		for( itm_j = 0; itm_j < curItem; itm_j++ ) {
			if( set[itm_j] == item )
				break;
		}
		if( itm_j < curItem )
			continue;

		set[curItem++] = item;
	}
}

double genFill( unsigned maxItems, unsigned range ) {
	unsigned	nIts = 100;
	NumberHasher*	hasher;
	unsigned*	items;
	unsigned*	stats;
	double		fill, avgFill;
	unsigned	itm_i, it_i;

	hasher = NumberHasher_New();
	Hasher_SetKeySize( hasher, sizeof(unsigned) );
	Hasher_SetMaxItems( hasher, maxItems );

	items = Memory_Alloc_Array_Unnamed( unsigned, maxItems );
	stats = Memory_Alloc_Array_Unnamed( unsigned, Hasher_GetTableSize( hasher ) );

	avgFill = 0.0;

	for( it_i = 0; it_i < nIts; it_i++ ) {
		genSet( maxItems, range, items );
		memset( stats, 0, Hasher_GetTableSize( hasher ) * sizeof(unsigned) );

		for( itm_i = 0; itm_i < maxItems; itm_i++ ) {
			unsigned	hashInd;

			hashInd = Hasher_Hash( hasher, items + itm_i );
			stats[hashInd]++;
		}

		fill = 0.0;

		for( itm_i = 0; itm_i < Hasher_GetTableSize( hasher ); itm_i++ ) {
			if( stats[itm_i] )
				fill += 1.0;
		}
		fill /= (double)maxItems;

		avgFill += fill;
	}

	avgFill /= (double)nIts;

	FreeArray( items );
	FreeArray( stats );
	FreeObject( hasher );

	return avgFill;
}


Bool testCompleteFill( unsigned rank, unsigned nProcs, unsigned watch ) {
	unsigned	maxItems = 100;
	double		avgFill;

	avgFill = genFill( maxItems, maxItems );

	return avgFill > (1.0 - 1e-8) && avgFill < (1.0 + 1e-8);
}

Bool testCloseFill( unsigned rank, unsigned nProcs, unsigned watch ) {
	unsigned	maxItems = 100;
	double		fillTol = 0.6;
	unsigned	rangeFac = 10;
	double		avgFill;

	avgFill = genFill( maxItems, maxItems * rangeFac );

	return avgFill >= fillTol;
}

Bool testFarFill( unsigned rank, unsigned nProcs, unsigned watch ) {
	unsigned	maxItems = 100;
	double		fillTol = 0.6;
	unsigned	rangeFac = 1000;
	double		avgFill;

	avgFill = genFill( maxItems, maxItems * rangeFac );

	return avgFill >= fillTol;
}


#define nTests	3

TestSuite_Test	tests[nTests] = {{"complete fill ratio", testCompleteFill, 1}, 
				 {"close fill ratio", testCloseFill, 100}, 
				 {"far fill ratio", testFarFill, 100}};


int main( int argc, char* argv[] ) {
	TestSuite*	suite;

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
