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
** $Id: testDecomp_Sync.c 2136 2004-09-30 02:47:13Z PatrickSunter $
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


Decomp* buildDecomp( unsigned rank, unsigned nProcs, unsigned* nLocals, unsigned** locals ) {
	Decomp*		decomp = Decomp_New();
	unsigned	l_i;

	*nLocals = 100;
	*locals = AllocArray( unsigned, *nLocals );
	for( l_i = 0; l_i < *nLocals; l_i++ )
		(*locals)[l_i] = rank * *nLocals + l_i;
	Decomp_SetLocals( decomp, *nLocals, *locals );

	return decomp;
}

void buildRequired( unsigned rank, unsigned nProcs, 
		    unsigned* nRequired, unsigned** required )
{
	unsigned	start;
	unsigned	r_i;

	*nRequired = 100;
	*required = AllocArray( unsigned, *nRequired );
	start = rank * (*nRequired - 10);

	for( r_i = 0; r_i < *nRequired; r_i++ )
		(*required)[r_i] = start + r_i;
}

CommTopology* buildCommTopo( unsigned rank, unsigned nProcs ) {
	CommTopology*	commTopo;
	unsigned	nInc, inc[2];

	commTopo = CommTopology_New();

	if( nProcs > 2 ) {
		nInc = (rank == 0 || rank == nProcs - 1) ? 1 : 2;
		if( rank > 0 )
			inc[0] = rank - 1;
		if( rank < nProcs - 1 )
			inc[(rank > 0) ? 1 : 0] = rank + 1;
	}
	else if( nProcs == 2 ) {
		nInc = 1;
		inc[0] = (rank == 0) ? 1 : 0;
	}
	else
		nInc = 0;

	CommTopology_SetIncidence( commTopo, nInc, inc );

	return commTopo;
}

Decomp_Sync* buildSync( unsigned rank, unsigned nProcs, 
			Decomp** decomp, unsigned* nLocals, unsigned** locals, 
			unsigned* remPerSide, unsigned* nRemotes, unsigned** remotes )
{
	Decomp_Sync*	sync;
	unsigned	r_i;

	sync = Decomp_Sync_New();
	*decomp = buildDecomp( rank, nProcs, nLocals, locals );
	Decomp_Sync_SetDecomp( sync, *decomp );
	Decomp_Sync_SetCommTopology( sync, buildCommTopo( rank, nProcs ) );

	*remPerSide = (*nLocals) / 10;
	*nRemotes = *remPerSide * ((rank > 0 && rank < nProcs - 1) ? 2 : 1) * ((nProcs > 1) ? 1 : 0);
	if( *nRemotes ) {
		*remotes = AllocArray( unsigned, *nRemotes );
		if( rank > 0 ) {
			for( r_i = 0; r_i < *remPerSide; r_i++ )
				(*remotes)[r_i] = rank * *nLocals - *remPerSide + r_i;
		}
		if( rank < nProcs - 1 ) {
			for( r_i = 0; r_i < *remPerSide; r_i++ ) {
				unsigned	ind = r_i + ((rank > 0) ? *remPerSide : 0);

				(*remotes)[ind] = (rank + 1) * *nLocals + r_i;
			}
		}
	}
	else
		*remotes = NULL;

	Decomp_Sync_SetRemotes( sync, *nRemotes, *remotes );

	return sync;
}


Bool testRemotes( unsigned rank, unsigned nProcs, unsigned watch ) {
	Bool		result = True;
	Decomp*		decomp;
	Decomp_Sync*	sync;
	unsigned	remPerSide;
	unsigned	nRemotes, *remotes;
	unsigned	nLocals, *locals;
	unsigned	r_i;

	sync = buildSync( rank, nProcs, 
			  &decomp, &nLocals, &locals, 
			  &remPerSide, &nRemotes, &remotes );

	if( rank == watch ) {
		if( sync->nRemotes != nRemotes ) {
			result = False;
			goto done;
		}

		if( rank > 0 ) {
			for( r_i = 0; r_i < remPerSide; r_i++ ) {
				if( sync->remotes[r_i] != remotes[r_i] ) {
					result = False;
					goto done;
				}
			}
		}
		if( rank < nProcs - 1 ) {
			for( r_i = 0; r_i < remPerSide; r_i++ ) {
				unsigned	ind = r_i + ((rank > 0) ? remPerSide : 0);

				if( sync->remotes[ind] != remotes[ind] ) {
					result = False;
					goto done;
				}
			}
		}
	}

done:
	FreeArray( locals );
	FreeArray( remotes );
	FreeObject( sync );
	FreeObject( decomp );

	return result;
}

Bool testSnkSrc( unsigned rank, unsigned nProcs, unsigned watch ) {
	Bool		result = True;
	Decomp*		decomp;
	Decomp_Sync*	sync;
	unsigned	remPerSide;
	unsigned	nRemotes, *remotes;
	unsigned	nLocals, *locals;

	sync = buildSync( rank, nProcs, 
			  &decomp, &nLocals, &locals, 
			  &remPerSide, &nRemotes, &remotes );

	if( rank == watch ) {
		if( sync->netSnks != nRemotes || sync->netSrcs != nRemotes ) {
			result = False;
			goto done;
		}
	}

done:
	FreeArray( locals );
	FreeArray( remotes );
	FreeObject( sync );
	FreeObject( decomp );

	return result;
}

typedef struct {
	int	one;
	int	two;
	int	three;
} theStruct;

Bool testArrays( unsigned rank, unsigned nProcs, unsigned watch ) {
	Bool			result = True;
	Decomp*			decomp;
	Decomp_Sync*		sync;
	unsigned		remPerSide;
	unsigned		nRemotes, *remotes;
	unsigned		nLocals, *locals;
	Decomp_Sync_Array*	arrays[2];
	int			*intLocals, *intRemotes;
	theStruct		*structLocals, *structRemotes;
	unsigned		r_i;

	sync = buildSync( rank, nProcs, 
			  &decomp, &nLocals, &locals, 
			  &remPerSide, &nRemotes, &remotes );

	if( nLocals ) {
		intLocals = AllocArray( int, nLocals );
		structLocals = AllocArray( theStruct, nLocals );
	}
	else {
		intLocals = NULL;
		structLocals = NULL;
	}
	if( nRemotes ) {
		intRemotes = AllocArray( int, nRemotes );
		structRemotes = AllocArray( theStruct, nRemotes );
	}
	else {
		intRemotes = NULL;
		structRemotes = NULL;
	}
	for( r_i = 0; r_i < nLocals; r_i++ ) {
		intLocals[r_i] = rank;
		structLocals[r_i].one = -rank;
		structLocals[r_i].two = rank;
		structLocals[r_i].three = -rank;
	}
	for( r_i = 0; r_i < nRemotes; r_i++ ) {
		intRemotes[r_i] = rank;
		structRemotes[r_i].one = -rank;
		structRemotes[r_i].two = rank;
		structRemotes[r_i].three = -rank;
	}

	arrays[0] = Decomp_Sync_Array_New();
	Decomp_Sync_Array_SetSync( arrays[0], sync );
	Decomp_Sync_Array_SetMemory( arrays[0], intLocals, intRemotes, sizeof(int), sizeof(int), sizeof(int) );

	arrays[1] = Decomp_Sync_Array_New();
	Decomp_Sync_Array_SetSync( arrays[1], sync );
	Decomp_Sync_Array_SetMemory( arrays[1], &structLocals[0].two, &structRemotes[0].two, 
				     sizeof(theStruct), sizeof(theStruct), sizeof(int) );

	Decomp_Sync_Array_Sync( arrays[0] );
	Decomp_Sync_Array_Sync( arrays[1] );

	if( rank == watch ) {
		for( r_i = 0; r_i < nLocals; r_i++ ) {
			if( intLocals[r_i] != rank || structLocals[r_i].two != rank || 
			    structLocals[r_i].one != -rank || structLocals[r_i].three != -rank )
			{
				break;
			}
		}
		if( r_i < nLocals ) {
			result = False;
			goto done;
		}

		if( rank > 0 ) {
			for( r_i = 0; r_i < remPerSide; r_i++ ) {
				if( intRemotes[r_i] != rank - 1 || structRemotes[r_i].two != rank - 1 || 
				    structRemotes[r_i].one != -rank || structRemotes[r_i].three != -rank )
				{
					break;
				}
			}
			if( r_i < remPerSide ) {
				result = False;
				goto done;
			}
		}
		if( rank < nProcs - 1 ) {
			for( r_i = 0; r_i < remPerSide; r_i++ ) {
				unsigned	ind = r_i + ((rank > 0) ? remPerSide : 0);

				if( intRemotes[ind] != rank + 1 || structRemotes[ind].two != rank + 1 || 
				    structRemotes[ind].one != -rank || structRemotes[ind].three != -rank )
				{
					break;
				}
			}
			if( r_i < remPerSide ) {
				result = False;
				goto done;
			}
		}
	}

done:
	FreeArray( intLocals );
	FreeArray( intRemotes );
	FreeArray( structLocals );
	FreeArray( structRemotes );
	FreeArray( locals );
	FreeArray( remotes );
	FreeObject( arrays[0] );
	FreeObject( arrays[1] );
	FreeObject( sync );
	FreeObject( decomp );

	return result;
}

Bool testClaim( unsigned rank, unsigned nProcs, unsigned watch ) {
	Bool		result = True;
	Decomposer*	decomposer;
	Decomp*		decomp;
	Decomp_Sync*	sync;
	unsigned	nRequired;
	unsigned*	required;

	buildRequired( rank, nProcs, &nRequired, &required );
	decomposer = Decomposer_New();
	decomp = NULL;
	sync = NULL;
	Decomposer_Decompose( decomposer, nRequired, required, 
			      NULL, &decomp, &sync );

	if( rank == watch ) {
		if( Decomp_Sync_GetGlobalSize( sync ) != nProcs * nRequired - (nProcs - 1) * 10 ) {
			result = False;
			goto done;
		}

		if( rank == 0 ) {
			if( Decomp_Sync_GetLocalSize( sync ) != nRequired || 
			    Decomp_Sync_GetRemoteSize( sync ) != 0 || 
			    Decomp_Sync_GetDomainSize( sync ) != nRequired || 
			    (nProcs > 1 && Decomp_Sync_GetSharedSize( sync ) != 10) )
			{
				result = False;
				goto done;
			}
		}
		else if( rank == nProcs - 1 ) {
			if( Decomp_Sync_GetLocalSize( sync ) != nRequired - 10 || 
			    Decomp_Sync_GetRemoteSize( sync ) != 10 || 
			    Decomp_Sync_GetDomainSize( sync ) != nRequired || 
			    Decomp_Sync_GetSharedSize( sync ) != 0 )
			{
				result = False;
				goto done;
			}
		}
		else {
			if( Decomp_Sync_GetLocalSize( sync ) != nRequired - 10 || 
			    Decomp_Sync_GetRemoteSize( sync ) != 10 || 
			    Decomp_Sync_GetDomainSize( sync ) != nRequired || 
			    Decomp_Sync_GetSharedSize( sync ) != 10 )
			{
				result = False;
				goto done;
			}
		}
	}

done:
	FreeArray( required );
	FreeObject( sync );
	FreeObject( decomp );
	FreeObject( decomposer );

	return result;
}


#define nTests	4

TestSuite*	suite;
TestSuite_Test	tests[nTests] = {{"set remotes", testRemotes, 10}, 
				 {"sink/sources", testSnkSrc, 10}, 
				 {"arrays", testArrays, 10}, 
				 {"claim", testClaim, 10}};


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
