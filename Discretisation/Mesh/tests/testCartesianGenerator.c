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
** $Id: testMeshTopology.c 2136 2004-09-30 02:47:13Z PatrickSunter $
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


Bool testGen1D( unsigned rank, unsigned nProcs, unsigned watch ) {
	Bool			result = True;
	CartesianGenerator*	gen;
	Mesh*			mesh;
	unsigned		sizes[3];
	double			minCrd[3];
	double			maxCrd[3];

	/*
	** Generate a single element per CPU.
	*/

	sizes[0] = nProcs;
	minCrd[0] = 0.0;
	maxCrd[0] = (double)nProcs;

	gen = CartesianGenerator_New( "" );
	MeshGenerator_SetDimSize( gen, 1 );
	CartesianGenerator_SetShadowDepth( gen, 0 );
	CartesianGenerator_SetTopologyParams( gen, sizes, 0, NULL, NULL );
	CartesianGenerator_SetGeometryParams( gen, minCrd, maxCrd );

	mesh = Mesh_New( "" );
	CartesianGenerator_Generate( gen, mesh );

	if( rank == watch ) {
		if( Mesh_GetDimSize( mesh ) != 1 || 
		    Mesh_GetGlobalSize( mesh, MT_VERTEX ) != nProcs + 1 || 
		    Mesh_GetGlobalSize( mesh, MT_EDGE ) != nProcs )
		{
			result = False;
			goto done;
		}

		if( rank == 0 ) {
			if( nProcs > 1 ) {
				if( Mesh_GetLocalSize( mesh, MT_VERTEX ) != 2 || 
				    Mesh_GetLocalSize( mesh, MT_EDGE ) != 1 || 
				    Mesh_GetRemoteSize( mesh, MT_VERTEX ) != 0 || 
				    Mesh_GetRemoteSize( mesh, MT_EDGE ) != 0 || 
				    Mesh_GetDomainSize( mesh, MT_VERTEX ) != 2 || 
				    Mesh_GetDomainSize( mesh, MT_EDGE ) != 1 || 
				    Mesh_GetSharedSize( mesh, MT_VERTEX ) != 1 || 
				    Mesh_GetSharedSize( mesh, MT_EDGE ) != 0 )
				{
					result = False;
					goto done;
				}
			}
			else {
				if( Mesh_GetLocalSize( mesh, MT_VERTEX ) != 2 || 
				    Mesh_GetLocalSize( mesh, MT_EDGE ) != 1 || 
				    Mesh_GetRemoteSize( mesh, MT_VERTEX ) != 0 || 
				    Mesh_GetRemoteSize( mesh, MT_EDGE ) != 0 || 
				    Mesh_GetDomainSize( mesh, MT_VERTEX ) != 2 || 
				    Mesh_GetDomainSize( mesh, MT_EDGE ) != 1 || 
				    Mesh_GetSharedSize( mesh, MT_VERTEX ) != 0 || 
				    Mesh_GetSharedSize( mesh, MT_EDGE ) != 0 )
				{
					result = False;
					goto done;
				}
			}
		}
		else if( rank == nProcs - 1 ) {
			if( Mesh_GetLocalSize( mesh, MT_VERTEX ) != 1 || 
			    Mesh_GetLocalSize( mesh, MT_EDGE ) != 1 || 
			    Mesh_GetRemoteSize( mesh, MT_VERTEX ) != 1 || 
			    Mesh_GetRemoteSize( mesh, MT_EDGE ) != 0 || 
			    Mesh_GetDomainSize( mesh, MT_VERTEX ) != 2 || 
			    Mesh_GetDomainSize( mesh, MT_EDGE ) != 1 || 
			    Mesh_GetSharedSize( mesh, MT_VERTEX ) != 0 || 
			    Mesh_GetSharedSize( mesh, MT_EDGE ) != 0 )
			{
				result = False;
				goto done;
			}
		}
		else {
			if( Mesh_GetLocalSize( mesh, MT_VERTEX ) != 1 || 
			    Mesh_GetLocalSize( mesh, MT_EDGE ) != 1 || 
			    Mesh_GetRemoteSize( mesh, MT_VERTEX ) != 1 || 
			    Mesh_GetRemoteSize( mesh, MT_EDGE ) != 0 || 
			    Mesh_GetDomainSize( mesh, MT_VERTEX ) != 2 || 
			    Mesh_GetDomainSize( mesh, MT_EDGE ) != 1 || 
			    Mesh_GetSharedSize( mesh, MT_VERTEX ) != 1 || 
			    Mesh_GetSharedSize( mesh, MT_EDGE ) != 0 )
			{
				result = False;
				goto done;
			}
		}
	}

done:
	FreeObject( gen );
	FreeObject( mesh );

	return result;
}

Bool testShadow1D( unsigned rank, unsigned nProcs, unsigned watch ) {
	Bool			result = True;
	CartesianGenerator*	gen;
	Mesh*			mesh;
	unsigned		sizes[3];
	double			minCrd[3];
	double			maxCrd[3];

	/*
	** Generate a single element per CPU.
	*/

	sizes[0] = nProcs;
	minCrd[0] = 0.0;
	maxCrd[0] = (double)nProcs;

	gen = CartesianGenerator_New( "" );
	MeshGenerator_SetDimSize( gen, 1 );
	CartesianGenerator_SetShadowDepth( gen, 1 );
	CartesianGenerator_SetTopologyParams( gen, sizes, 0, NULL, NULL );
	CartesianGenerator_SetGeometryParams( gen, minCrd, maxCrd );

	mesh = Mesh_New( "" );
	CartesianGenerator_Generate( gen, mesh );

	if( rank == watch ) {
		if( Mesh_GetDimSize( mesh ) != 1 || 
		    Mesh_GetGlobalSize( mesh, MT_VERTEX ) != nProcs + 1 || 
		    Mesh_GetGlobalSize( mesh, MT_EDGE ) != nProcs )
		{
			result = False;
			goto done;
		}

		if( rank == 0 ) {
			if( nProcs > 1 ) {
				if( Mesh_GetLocalSize( mesh, MT_VERTEX ) != 2 || 
				    Mesh_GetLocalSize( mesh, MT_EDGE ) != 1 || 
				    Mesh_GetRemoteSize( mesh, MT_VERTEX ) != 1 || 
				    Mesh_GetRemoteSize( mesh, MT_EDGE ) != 1 || 
				    Mesh_GetDomainSize( mesh, MT_VERTEX ) != 3 || 
				    Mesh_GetDomainSize( mesh, MT_EDGE ) != 2 || 
				    Mesh_GetSharedSize( mesh, MT_VERTEX ) != 2 || 
				    Mesh_GetSharedSize( mesh, MT_EDGE ) != 1 )
				{
					result = False;
					goto done;
				}
			}
			else {
				if( Mesh_GetLocalSize( mesh, MT_VERTEX ) != 2 || 
				    Mesh_GetLocalSize( mesh, MT_EDGE ) != 1 || 
				    Mesh_GetRemoteSize( mesh, MT_VERTEX ) != 0 || 
				    Mesh_GetRemoteSize( mesh, MT_EDGE ) != 0 || 
				    Mesh_GetDomainSize( mesh, MT_VERTEX ) != 2 || 
				    Mesh_GetDomainSize( mesh, MT_EDGE ) != 1 || 
				    Mesh_GetSharedSize( mesh, MT_VERTEX ) != 0 || 
				    Mesh_GetSharedSize( mesh, MT_EDGE ) != 0 )
				{
					result = False;
					goto done;
				}
			}
		}
		else if( rank == nProcs - 1 ) {
			if( Mesh_GetLocalSize( mesh, MT_VERTEX ) != 1 || 
			    Mesh_GetLocalSize( mesh, MT_EDGE ) != 1 || 
			    Mesh_GetRemoteSize( mesh, MT_VERTEX ) != 2 || 
			    Mesh_GetRemoteSize( mesh, MT_EDGE ) != 1 || 
			    Mesh_GetDomainSize( mesh, MT_VERTEX ) != 3 || 
			    Mesh_GetDomainSize( mesh, MT_EDGE ) != 2 || 
			    Mesh_GetSharedSize( mesh, MT_VERTEX ) != 1 || 
			    Mesh_GetSharedSize( mesh, MT_EDGE ) != 1 )
			{
				result = False;
				goto done;
			}
		}
		else {
			if( Mesh_GetLocalSize( mesh, MT_VERTEX ) != 1 || 
			    Mesh_GetLocalSize( mesh, MT_EDGE ) != 1 || 
			    Mesh_GetRemoteSize( mesh, MT_VERTEX ) != 3 || 
			    Mesh_GetRemoteSize( mesh, MT_EDGE ) != 2 || 
			    Mesh_GetDomainSize( mesh, MT_VERTEX ) != 4 || 
			    Mesh_GetDomainSize( mesh, MT_EDGE ) != 3 || 
			    Mesh_GetSharedSize( mesh, MT_VERTEX ) != 1 || 
			    Mesh_GetSharedSize( mesh, MT_EDGE ) != 1 )
			{
				result = False;
				goto done;
			}
		}
	}

done:
	FreeObject( gen );
	FreeObject( mesh );

	return result;
}


#define nTests	2

TestSuite_Test	tests[nTests] = {{"generate 1D mesh", testGen1D, 10}, 
				 {"generate shadowed 1D mesh", testShadow1D, 10}};


int main( int argc, char* argv[] ) {
	TestSuite*	suite;

	/* Initialise MPI, get world info. */
	MPI_Init( &argc, &argv );

	/* Initialise StGermain. */
	Base_Init( &argc, &argv );
	DiscretisationMesh_Init( &argc, &argv );

	/* Create the test suite. */
	suite = TestSuite_New();
	TestSuite_SetProcToWatch( suite, (argc >= 2) ? atoi( argv[1] ) : 0 );
	TestSuite_SetTests( suite, nTests, tests );

	/* Run the tests. */
	TestSuite_Run( suite );

	/* Destroy test suites. */
	FreeObject( suite );

	/* Finalise StGermain. */
	DiscretisationMesh_Finalise();
	Base_Finalise();

	/* Close off MPI */
	MPI_Finalize();

	return MPI_SUCCESS;
}
