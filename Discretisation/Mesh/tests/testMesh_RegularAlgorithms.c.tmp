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
** $Id: testMesh.c 3911 2006-12-22 01:37:13Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <mpi.h>
#include <mpi.h>

#include "Base/Base.h"
#include "Discretisation/Geometry/Geometry.h"
#include "Discretisation/Shape/Shape.h"
#include "Discretisation/Mesh/Mesh.h"


Mesh* buildMesh( unsigned nProcs, unsigned nDims ) {
	CartesianGenerator*	gen;
	Mesh*			mesh;
	unsigned*		sizes;
	double			*minCrd, *maxCrd;
	unsigned		d_i;

	sizes = AllocArray( unsigned, nDims );
	minCrd = AllocArray( double, nDims );
	maxCrd = AllocArray( double, nDims );
	for( d_i = 0; d_i < nDims; d_i++ ) {
		sizes[d_i] = nProcs;
		minCrd[d_i] = 0.0;
		maxCrd[d_i] = (double)nProcs;
	}

	gen = CartesianGenerator_New( "" );
	MeshGenerator_SetDimSize( gen, nDims );
	CartesianGenerator_SetTopologyParams( gen, sizes, 0, NULL, NULL );
	CartesianGenerator_SetGeometryParams( gen, minCrd, maxCrd );

	mesh = Mesh_New( "" );
	Mesh_SetGenerator( mesh, gen );
	Mesh_SetAlgorithms( mesh, Mesh_RegularAlgorithms_New( "" ) );
	Build( mesh, NULL, False );

	FreeObject( gen );
	FreeArray( sizes );
	FreeArray( minCrd );
	FreeArray( maxCrd );

	return mesh;
}


Bool testElSearch1D( unsigned rank, unsigned nProcs, unsigned watch ) {
	Bool	result = True;
	Mesh*	mesh;

	mesh = buildMesh( nProcs, 1 );

	if( rank == watch ) {
		unsigned	e_i;

		for( e_i = 0; e_i < Mesh_GetDomainSize( mesh, MT_EDGE ); e_i++ ) {
			unsigned	nInc, *inc;
			double		point[1];
			unsigned	elDim, elInd;
			unsigned	inc_i;

			Mesh_GetIncidence( mesh, MT_EDGE, e_i, MT_VERTEX, &nInc, &inc );
			for( inc_i = 0; inc_i < nInc; inc_i++ ) {
				double*	vert;

				vert = Mesh_GetVertex( mesh, inc[inc_i] );
				if( !Mesh_Search( mesh, vert, &elDim, &elInd ) || 
				    elDim != MT_VERTEX || 
				    elInd != inc[inc_i] )
				{
					result = False;
					goto done;
				}
			}

			Mesh_GetIncidence( mesh, MT_EDGE, e_i, MT_VERTEX, &nInc, &inc );
			point[0] = 0.0;
			for( inc_i = 0; inc_i < nInc; inc_i++ )
				point[0] += Mesh_GetVertex( mesh, inc[inc_i] )[0];
			point[0] /= (double)nInc;
			if( !Mesh_Search( mesh, point, &elDim, &elInd ) || 
			    elDim != MT_EDGE || 
			    elInd != e_i )
			{
				result = False;
				goto done;
			}
		}
	}

done:
	FreeObject( mesh );

	return result;
}

Bool testElSearch2D( unsigned rank, unsigned nProcs, unsigned watch ) {
	Bool	result = True;
	Mesh*	mesh;

	mesh = buildMesh( nProcs, 2 );

	if( rank == watch ) {
		unsigned	e_i;

		for( e_i = 0; e_i < Mesh_GetDomainSize( mesh, MT_FACE ); e_i++ ) {
			unsigned	nInc, *inc;
			double		point[2];
			unsigned	elDim, elInd;
			unsigned	inc_i;

			Mesh_GetIncidence( mesh, MT_FACE, e_i, MT_VERTEX, &nInc, &inc );
			for( inc_i = 0; inc_i < nInc; inc_i++ ) {
				double*	vert;

				vert = Mesh_GetVertex( mesh, inc[inc_i] );
				if( !Mesh_Search( mesh, vert, &elDim, &elInd ) || 
				    elDim != MT_VERTEX || 
				    elInd != inc[inc_i] )
				{
					result = False;
					goto done;
				}
			}

			Mesh_GetIncidence( mesh, MT_FACE, e_i, MT_EDGE, &nInc, &inc );
			for( inc_i = 0; inc_i < nInc; inc_i++ ) {
				unsigned	nEdgeInc, *edgeInc;
				unsigned	inc_j;

				Mesh_GetIncidence( mesh, MT_EDGE, inc[inc_i], MT_VERTEX, &nEdgeInc, &edgeInc );
				point[0] = point[1] = 0.0;
				for( inc_j = 0; inc_j < nEdgeInc; inc_j++ ) {
					point[0] += Mesh_GetVertex( mesh, edgeInc[inc_j] )[0];
					point[1] += Mesh_GetVertex( mesh, edgeInc[inc_j] )[1];
				}
				point[0] /= (double)nEdgeInc;
				point[1] /= (double)nEdgeInc;
				if( !Mesh_Search( mesh, point, &elDim, &elInd ) || 
				    elDim != MT_EDGE || 
				    elInd != inc[inc_i] )
				{
					result = False;
					goto done;
				}
			}

			Mesh_GetIncidence( mesh, MT_FACE, e_i, MT_VERTEX, &nInc, &inc );
			point[0] = point[1] = 0.0;
			for( inc_i = 0; inc_i < nInc; inc_i++ ) {
				point[0] += Mesh_GetVertex( mesh, inc[inc_i] )[0];
				point[1] += Mesh_GetVertex( mesh, inc[inc_i] )[1];
			}
			point[0] /= (double)nInc;
			point[1] /= (double)nInc;
			if( !Mesh_Search( mesh, point, &elDim, &elInd ) || 
			    elDim != MT_FACE || 
			    elInd != e_i )
			{
				result = False;
				goto done;
			}
		}
	}

done:
	FreeObject( mesh );

	return result;
}

Bool testElSearch3D( unsigned rank, unsigned nProcs, unsigned watch ) {
	Bool	result = True;
	Mesh*	mesh;

	mesh = buildMesh( nProcs, 3 );

	if( rank == watch ) {
		unsigned	e_i;

		for( e_i = 0; e_i < Mesh_GetDomainSize( mesh, MT_VOLUME ); e_i++ ) {
			unsigned	nInc, *inc;
			double		point[3];
			unsigned	elDim, elInd;
			unsigned	inc_i;

			Mesh_GetIncidence( mesh, MT_VOLUME, e_i, MT_VERTEX, &nInc, &inc );
			for( inc_i = 0; inc_i < nInc; inc_i++ ) {
				double*	vert;

				vert = Mesh_GetVertex( mesh, inc[inc_i] );
				if( !Mesh_Search( mesh, vert, &elDim, &elInd ) || 
				    elDim != MT_VERTEX || 
				    elInd != inc[inc_i] )
				{
					result = False;
					goto done;
				}
			}

			Mesh_GetIncidence( mesh, MT_VOLUME, e_i, MT_EDGE, &nInc, &inc );
			for( inc_i = 0; inc_i < nInc; inc_i++ ) {
				unsigned	nEdgeInc, *edgeInc;
				unsigned	inc_j;

				Mesh_GetIncidence( mesh, MT_EDGE, inc[inc_i], MT_VERTEX, &nEdgeInc, &edgeInc );
				point[0] = point[1] = point[2] = 0.0;
				for( inc_j = 0; inc_j < nEdgeInc; inc_j++ ) {
					point[0] += Mesh_GetVertex( mesh, edgeInc[inc_j] )[0];
					point[1] += Mesh_GetVertex( mesh, edgeInc[inc_j] )[1];
					point[2] += Mesh_GetVertex( mesh, edgeInc[inc_j] )[2];
				}
				point[0] /= (double)nEdgeInc;
				point[1] /= (double)nEdgeInc;
				point[2] /= (double)nEdgeInc;
				if( !Mesh_Search( mesh, point, &elDim, &elInd ) || 
				    elDim != MT_EDGE || 
				    elInd != inc[inc_i] )
				{
					result = False;
					goto done;
				}
			}

			Mesh_GetIncidence( mesh, MT_VOLUME, e_i, MT_FACE, &nInc, &inc );
			for( inc_i = 0; inc_i < nInc; inc_i++ ) {
				unsigned	nFaceInc, *faceInc;
				unsigned	inc_j;

				Mesh_GetIncidence( mesh, MT_FACE, inc[inc_i], MT_VERTEX, &nFaceInc, &faceInc );
				point[0] = point[1] = point[2] = 0.0;
				for( inc_j = 0; inc_j < nFaceInc; inc_j++ ) {
					point[0] += Mesh_GetVertex( mesh, faceInc[inc_j] )[0];
					point[1] += Mesh_GetVertex( mesh, faceInc[inc_j] )[1];
					point[2] += Mesh_GetVertex( mesh, faceInc[inc_j] )[2];
				}
				point[0] /= (double)nFaceInc;
				point[1] /= (double)nFaceInc;
				point[2] /= (double)nFaceInc;
				if( !Mesh_Search( mesh, point, &elDim, &elInd ) || 
				    elDim != MT_FACE || 
				    elInd != inc[inc_i] )
				{
					result = False;
					goto done;
				}
			}

			Mesh_GetIncidence( mesh, MT_VOLUME, e_i, MT_VERTEX, &nInc, &inc );
			point[0] = point[1] = point[2] = 0.0;
			for( inc_i = 0; inc_i < nInc; inc_i++ ) {
				point[0] += Mesh_GetVertex( mesh, inc[inc_i] )[0];
				point[1] += Mesh_GetVertex( mesh, inc[inc_i] )[1];
				point[2] += Mesh_GetVertex( mesh, inc[inc_i] )[2];
			}
			point[0] /= (double)nInc;
			point[1] /= (double)nInc;
			point[2] /= (double)nInc;
			if( !Mesh_Search( mesh, point, &elDim, &elInd ) || 
			    elDim != MT_VOLUME || 
			    elInd != e_i )
			{
				result = False;
				goto done;
			}
		}
	}

done:
	FreeObject( mesh );

	return result;
}


#define nTests	3

TestSuite_Test	tests[nTests] = {{"test element search (1D)", testElSearch1D, 5}, 
				 {"test element search (2D)", testElSearch2D, 5}, 
				 {"test element search (3D)", testElSearch3D, 5}};


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
