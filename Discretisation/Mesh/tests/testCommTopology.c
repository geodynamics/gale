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
** $Id: testCommTopology.c 2136 2004-09-30 02:47:13Z PatrickSunter $
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


CommTopology* buildCompleteIncidence( unsigned rank, unsigned nProcs, 
				      unsigned* nInc, unsigned** inc )
{
	CommTopology*	commTopo;
	unsigned	inc_i;

	commTopo = CommTopology_New();

	*nInc = nProcs - 1;
	*inc = AllocArray( unsigned, *nInc );
	for( inc_i = 0; inc_i < *nInc; inc_i++ )
		(*inc)[inc_i] = (rank + inc_i + 1) % nProcs;

	CommTopology_SetIncidence( commTopo, *nInc, *inc );

	return commTopo;
}

CommTopology* buildPartialIncidence( unsigned rank, unsigned nProcs, 
				     unsigned* nInc, unsigned** inc )
{
	CommTopology*	commTopo;

	commTopo = CommTopology_New();

	if( nProcs > 1 ) {
		*nInc = (rank == 0 || rank == nProcs - 1) ? 1 : 2;
		*inc = Memory_Alloc_Array_Unnamed( unsigned, *nInc );

		if( rank > 0 )
			(*inc)[0] = rank - 1;
		if( rank < nProcs - 1 )
			(*inc)[(rank > 0) ? 1 : 0] = rank + 1;

		CommTopology_SetIncidence( commTopo, *nInc, *inc );
	}
	else {
		*nInc = 0;
		*inc = NULL;
	}

	return commTopo;
}

CommTopology* buildAddIncidence( unsigned rank, unsigned nProcs, 
				 unsigned* nInc, unsigned** inc )
{
	CommTopology*	commTopo;

	commTopo = CommTopology_New();

	if( nProcs > 1 ) {
		*nInc = (rank == 0 || rank == nProcs - 1) ? 1 : 2;
		*inc = AllocArray( unsigned, *nInc );

		if( rank > 0 )
			(*inc)[0] = rank - 1;
		CommTopology_SetIncidence( commTopo, (rank > 0) ? 1 : 0, *inc );

		if( rank < nProcs - 1 )
			(*inc)[(rank > 0) ? 1 : 0] = rank + 1;
		CommTopology_AddIncidence( commTopo, 
					   (rank < nProcs - 1) ? 1 : 0, 
					   (*inc) + ((rank > 0) ? 1 : 0) );
	}
	else {
		*nInc = 0;
		*inc = NULL;
	}

	return commTopo;
}

void buildGatherArrays( unsigned rank, unsigned nProcs, 
			unsigned* srcSize, unsigned** src )
{
	unsigned	p_i;

	*srcSize = nProcs;
	*src = Memory_Alloc_Array_Unnamed( unsigned, *srcSize );
	for( p_i = 0; p_i < *srcSize; p_i++ )
		(*src)[p_i] = rank;
}

void buildAlltoallArrays( unsigned rank, unsigned nProcs, CommTopology* commTopo, 
			  unsigned** srcSizes, unsigned*** srcs )
{
	unsigned	nIncRanks;
	unsigned	p_i, p_j;

	if( nProcs == 1 ) {
		*srcSizes = 0;
		*srcs = NULL;
		return;
	}

	nIncRanks = CommTopology_GetIncidenceSize( commTopo );
	*srcSizes = Memory_Alloc_Array_Unnamed( unsigned, nIncRanks );
	*srcs = Memory_Alloc_2DArray_Unnamed( unsigned, nIncRanks, nProcs );
	for( p_i = 0; p_i < nIncRanks; p_i++ ) {
		(*srcSizes)[p_i] = nProcs;
		for( p_j = 0; p_j < nProcs; p_j++ )
			(*srcs)[p_i][p_j] = rank;
	}
}


Bool testSetComm( unsigned rank, unsigned nProcs, unsigned watch ) {
	CommTopology*	commTopo;
	Bool		result = True;

	commTopo = CommTopology_New();

	CommTopology_SetComm( commTopo, MPI_COMM_WORLD );
	if( rank == watch ) {
		if( commTopo->comm != MPI_COMM_WORLD || 
		    commTopo->nIncRanks != 0 || 
		    commTopo->incRanks != NULL || 
		    commTopo->order != NULL || 
		    UIntMap_GetSize( commTopo->glMap ) )
		{
			result = False;
			goto done;
		}
	}

done:
	FreeObject( commTopo );

	return result;
}

Bool testEmptyInc( unsigned rank, unsigned nProcs, unsigned watch ) {
	CommTopology*	commTopo;
	Bool		result = True;

	commTopo = CommTopology_New();
	CommTopology_SetIncidence( commTopo, 0, NULL );

	if( rank == watch ) {
		if( commTopo->comm != MPI_COMM_WORLD || 
		    commTopo->nIncRanks != 0 || 
		    commTopo->incRanks != NULL || 
		    commTopo->order != NULL || 
		    UIntMap_GetSize( commTopo->glMap ) )
		{
			result = False;
			goto done;
		}
	}

done:
	FreeObject( commTopo );

	return result;
}

Bool testCompInc( unsigned rank, unsigned nProcs, unsigned watch ) {
	CommTopology*	commTopo;
	Bool		result = True;
	unsigned	nInc;
	unsigned*	inc;
	unsigned	inc_i;

	commTopo = buildCompleteIncidence( rank, nProcs, &nInc, &inc );

	if( rank == watch ) {
		if( commTopo->comm != MPI_COMM_WORLD || 
		    commTopo->nIncRanks != nInc )
		{
			result = False;
			goto done;
		}

		if( nProcs > 1 && 
		    (!commTopo->incRanks || 
		     !commTopo->order || 
		     !UIntMap_GetSize( commTopo->glMap )) )
		{
			result = False;
			goto done;
		}

		for( inc_i = 0; inc_i < nInc; inc_i++ ) {
			unsigned	local;

			if( CommTopology_LocalToGlobal( commTopo, inc_i ) != inc[inc_i] || 
			    !CommTopology_GlobalToLocal( commTopo, inc[inc_i], &local ) || 
			    local != inc_i )
			{
				break;
			}
		}
		if( inc_i < nInc ) {
			result = False;
			goto done;
		}
	}

done:
	FreeArray( inc );
	FreeObject( commTopo );

	return result;
}

Bool testPartInc( unsigned rank, unsigned nProcs, unsigned watch ) {
	CommTopology*	commTopo;
	Bool		result = True;
	unsigned	nInc;
	unsigned*	inc;
	unsigned	inc_i;

	commTopo = buildPartialIncidence( rank, nProcs, &nInc, &inc );

	if( rank == watch ) {
		if( commTopo->comm != MPI_COMM_WORLD || 
		    commTopo->nIncRanks != nInc )
		{
			result = False;
			goto done;
		}

		if( nProcs > 1 && 
		    (!commTopo->incRanks || 
		     !commTopo->order || 
		     !UIntMap_GetSize( commTopo->glMap )) )
		{
			result = False;
			goto done;
		}

		for( inc_i = 0; inc_i < nInc; inc_i++ ) {
			unsigned	local;

			if( CommTopology_LocalToGlobal( commTopo, inc_i ) != inc[inc_i] || 
			    !CommTopology_GlobalToLocal( commTopo, inc[inc_i], &local ) || 
			    local != inc_i )
			{
				break;
			}
		}
		if( inc_i < nInc ) {
			result = False;
			goto done;
		}
	}

done:
	FreeArray( inc );
	FreeObject( commTopo );

	return result;
}

Bool testAddInc( unsigned rank, unsigned nProcs, unsigned watch ) {
	CommTopology*	commTopo;
	Bool		result = True;
	unsigned	nInc;
	unsigned*	inc;
	unsigned	inc_i;

	commTopo = buildAddIncidence( rank, nProcs, &nInc, &inc );

	if( rank == watch ) {
		if( commTopo->comm != MPI_COMM_WORLD || 
		    commTopo->nIncRanks != nInc )
		{
			result = False;
			goto done;
		}

		if( nProcs > 1 && 
		    (!commTopo->incRanks || 
		     !commTopo->order || 
		     !UIntMap_GetSize( commTopo->glMap )) )
		{
			result = False;
			goto done;
		}

		for( inc_i = 0; inc_i < nInc; inc_i++ ) {
			unsigned	local;

			if( CommTopology_LocalToGlobal( commTopo, inc_i ) != inc[inc_i] || 
			    !CommTopology_GlobalToLocal( commTopo, inc[inc_i], &local ) || 
			    local != inc_i )
			{
				break;
			}
		}
		if( inc_i < nInc ) {
			result = False;
			goto done;
		}
	}

done:
	FreeArray( inc );
	FreeObject( commTopo );

	return result;
}

Bool testAllgather( unsigned rank, unsigned nProcs, unsigned watch ) {
	CommTopology*	commTopo;
	Bool		result = True;
	unsigned	nInc;
	unsigned*	inc;
	unsigned	srcSize;
	unsigned*	src;
	unsigned*	dstSizes;
	unsigned**	dst;
	unsigned	p_i, e_i;

	commTopo = buildPartialIncidence( rank, nProcs, &nInc, &inc );
	buildGatherArrays( rank, nProcs, &srcSize, &src );

	CommTopology_Allgather( commTopo, 
				srcSize, src, 
				&dstSizes, (void***)&dst, 
				sizeof(unsigned) );

	if( rank == watch ) {
		for( p_i = 0; p_i < nInc; p_i++ ) {
			for( e_i = 0; e_i < srcSize; e_i++ ) {
				if( dst[p_i][e_i] != inc[p_i] ) {
					result = False;
					goto done;
				}
			}
		}
	}

done:
	FreeObject( commTopo );
	FreeArray( inc );
	FreeArray( src );
	FreeArray( dst );

	return result;
}

Bool testAlltoall( unsigned rank, unsigned nProcs, unsigned watch ) {
	CommTopology*	commTopo;
	Bool		result = True;
	unsigned	nInc;
	unsigned*	inc;
	unsigned*	srcSizes;
	unsigned**	srcs;
	unsigned*	dstSizes;
	unsigned**	dst;
	unsigned	p_i, e_i;

	commTopo = buildPartialIncidence( rank, nProcs, &nInc, &inc );
	buildAlltoallArrays( rank, nProcs, commTopo, &srcSizes, &srcs );

	CommTopology_Alltoall( commTopo, 
			       srcSizes, (void**)srcs, 
			       &dstSizes, (void***)&dst, 
			       sizeof(unsigned) );

	if( rank == watch ) {
		for( p_i = 0; p_i < nInc; p_i++ ) {
			for( e_i = 0; e_i < srcSizes[p_i]; e_i++ ) {
				if( dst[p_i][e_i] != inc[p_i] ) {
					result = False;
					goto done;
				}
			}
		}
	}

done:
	FreeObject( commTopo );
	FreeArray( inc );
	FreeArray( srcSizes );
	FreeArray( srcs );
	FreeArray( dst );

	return result;
}


#define nTests	7

TestSuite_Test	tests[nTests] = {{"set communicator", testSetComm, 10}, 
				 {"empty incidence", testEmptyInc, 10 }, 
				 {"partial incidence", testPartInc, 10}, 
				 {"complete incidence", testCompInc, 10}, 
				 {"adding incidence", testAddInc, 10}, 
				 {"all gather", testAllgather, 10}, 
				 {"all to all", testAlltoall, 10}};


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
