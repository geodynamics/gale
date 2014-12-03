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
** Role:
**   Tests the ParallelDelaunaySuite
**
** $Id: testParallelDelaunay.c 3462 2006-02-19 06:53:24Z WalterLandry $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>

#include "pcu/pcu.h"
#include <StGermain/StGermain.h>
#include <StgDomain/Geometry/Geometry.h>
#include "ParallelDelaunaySuite.h"

typedef enum pointsType_t{
	Irregular,
	Regular,
	Polygon
}pointsType;

#define PI 3.1415926535897932384626
#define EPS 0.0001
#define MAX_NEIGHBOURS 10000

typedef struct {
	DelaunayAttributes	attr;
	Dictionary*		dict;
} ParallelDelaunaySuiteData;

void _GeneratePoints( CoordF* sites, int numSites, pointsType *p ) {
	int i, j, count;
	int num = numSites;

	if( *p == Irregular ) {
		for( i = 0; i < num; i++ ) {
			sites[i][0] = drand48();
			sites[i][1] = drand48();
		}
	}
	else if( *p == Regular ) {
		count = 0;
		num = sqrt((float)num);
		for( i = 0; i < num; i++ ) {
			for( j = 0; j < num; j++ ) {
				sites[count][0] = ((float)i)*(1.0/((float)num));
				sites[count][1] = ((float)j)*(1.0/((float)num));

				count++;
			}
		}
	}
	else if( *p == Polygon ) {
		for( i = 0; i < num-1; i++ ) {
			sites[i][0] = cos( ((float)i)*2.0*PI/((float)(num-1)) );
			sites[i][1] = sin( ((float)i)*2.0*PI/((float)(num-1)) );
		}
	}
}

int CompareFunction( const void* a, const void* b ) {
	int *p, *q;

	p = (int*)a;
	q = (int*)b;

	if( *p > *q )
		return 1;
	else if( *p < *q )
		return -1;
	else
		return 0;
}

void ParallelDelaunaySuite_Setup( ParallelDelaunaySuiteData* data ) {
	data->attr.BuildBoundingTriangle = 0;
	data->attr.BuildTriangleIndices = 1;
	data->attr.BuildTriangleNeighbours = 1;
	data->attr.CreateVoronoiVertices = 1;
	data->attr.CalculateVoronoiSides = 1;
	data->attr.CalculateVoronoiSurfaceArea = 1;
	data->attr.FindNeighbours = 1;

	data->dict = Dictionary_New();
}

void ParallelDelaunaySuite_Teardown( ParallelDelaunaySuiteData* data ) {
	Stg_Class_Delete( data->dict );
}

/* NOTE: neither this test, nor the next one, works for a triangulation of 900 sites
 *       or greater... note sure why
 */
void ParallelDelaunaySuite_TestIrregular( ParallelDelaunaySuiteData* data ) {
	ParallelDelaunay*	pd;
	Delaunay*		d;
	CoordF*			sites;
	pointsType		p 		= Irregular;
	int			numSites[4]	= { 100, 400, 900, 1600 };
	int			nsites;
	MPI_Comm		mpiComm;
	int			rank, nRanks;
	int			i, j, k = 0, m;
	int			dNeighbours[MAX_NEIGHBOURS];
	int			pdNeighbours[MAX_NEIGHBOURS];
	int			voronoiAreaTest	= 1;
	int			voronoiSideTest	= 1;
	int			triangleCountTest = 1;

	MPI_Comm_dup( MPI_COMM_WORLD, &mpiComm );
	MPI_Comm_rank( mpiComm, &rank );
	MPI_Comm_size( mpiComm, &nRanks );

	for( i = 0; i < 4; i++ ) {
		nsites = numSites[i];

		if( rank == 0 ) {
			sites = Memory_Alloc_Array( CoordF, nsites, "testDelaunayParallel_sites" );
			memset( sites, 0, sizeof( CoordF ) * nsites );
			_GeneratePoints( sites, nsites, &p );
		}

		pd = ParallelDelaunay_New( "Delaunay-Parallel", data->dict, sites, nsites, rank, nRanks, &mpiComm, &data->attr );
		Stg_Component_Build( pd, NULL, True );
		Stg_Component_Execute( pd, NULL, True );

		ParallelDelaunay_GatherTriangulation( pd );
		MPI_Barrier( mpiComm );

		if( rank == 0 ) {
			d = Delaunay_New( "Delaunay-Serial", data->dict, sites, nsites, 0, &data->attr );
			Stg_Component_Build( d, NULL, True );

			for( j = 0; j < nsites; j++ ) {
				if( fabs( d->voronoiArea[j] - pd->voronoiArea[j] ) > EPS )
					voronoiAreaTest = 0;

				if( dNeighbours[k] != pdNeighbours[k] )
					continue;

				memset( dNeighbours, 0, sizeof( dNeighbours ) );
				memset( pdNeighbours, 0, sizeof( pdNeighbours ) );

				memcpy( dNeighbours, d->neighbours[j], sizeof( int ) * d->numNeighbours[j] );
				memcpy( pdNeighbours, pd->neighbours[j], sizeof( int ) * pd->numNeighbours[j] );

				qsort( dNeighbours, d->numNeighbours[j], sizeof( int ), CompareFunction );
				qsort( pdNeighbours, pd->numNeighbours[j], sizeof( int ), CompareFunction );

				for( k = 0; k < (int)(d->numNeighbours[j]); k++ ) {
					if( dNeighbours[k] != pdNeighbours[k] )
						continue;

					for( m = 0; m < (int)(d->numNeighbours[j]); m++ ) {
						if( d->neighbours[j][k] == pd->neighbours[j][m] )
							if( fabs( d->voronoiSides[j][k] - pd->voronoiSides[j][m] ) > EPS )
								voronoiSideTest = 0;
					}
				}
			}

			if( d->numTriangles != pd->numTriangles )
				triangleCountTest = 0;

			Stg_Class_Delete( d );
		}

		Stg_Class_Delete( pd );
		
		if( rank == 0 ) {
			Memory_Free( sites );

			pcu_check_true( voronoiAreaTest );
			pcu_check_true( voronoiSideTest );
			pcu_check_true( triangleCountTest );
		}
	}
}

void ParallelDelaunaySuite_TestRegular( ParallelDelaunaySuiteData* data ) {
	ParallelDelaunay*	pd;
	Delaunay*		d;
	CoordF*			sites;
	pointsType		p 		= Regular;
	int			numSites[4]	= { 100, 400, 900, 1600 };
	MPI_Comm		mpiComm;
	int			rank, nRanks;
	int			i, j, k, m;
	int			dNeighbours[MAX_NEIGHBOURS];
	int			pdNeighbours[MAX_NEIGHBOURS];
	int			voronoiAreaTest	= 1;
	int			voronoiSideTest	= 1;
	int			triangleCountTest = 1;

	MPI_Comm_dup( MPI_COMM_WORLD, &mpiComm );
	MPI_Comm_rank( mpiComm, &rank );
	MPI_Comm_size( mpiComm, &nRanks );

	for( i = 0; i < 4; i++ ) {
		if( rank == 0 ) {
			sites = Memory_Alloc_Array( CoordF, numSites[i], "testDelaunayParallel_sites" );
			memset( sites, 0, sizeof( CoordF ) * numSites[i] );
			_GeneratePoints( sites, numSites[i], &p );
		}

		pd = ParallelDelaunay_New( "Delaunay-Parallel", data->dict, sites, numSites[i], rank, nRanks, &mpiComm, &data->attr );
		Stg_Component_Build( pd, NULL, True );
		Stg_Component_Execute( pd, NULL, True );

		ParallelDelaunay_GatherTriangulation( pd );
		MPI_Barrier( mpiComm );

		if( rank == 0 ) {
			d = Delaunay_New( "Delaunay-Serial", data->dict, sites, numSites[i], 0, &data->attr );
			Stg_Component_Build( d, NULL, True );

			for( j = 0; j < numSites[i]; j++ ) {
				if( fabs( d->voronoiArea[j] - pd->voronoiArea[j] ) > EPS )
					voronoiAreaTest = 0;

				memset( dNeighbours, 0, sizeof( dNeighbours ) );
				memset( pdNeighbours, 0, sizeof( pdNeighbours ) );

				memcpy( dNeighbours, d->neighbours[j], sizeof( int ) * d->numNeighbours[j] );
				memcpy( pdNeighbours, pd->neighbours[j], sizeof( int ) * pd->numNeighbours[j] );

				qsort( dNeighbours, d->numNeighbours[j], sizeof( int ), CompareFunction );
				qsort( pdNeighbours, pd->numNeighbours[j], sizeof( int ), CompareFunction );

				for( k = 0; k < (int)(d->numNeighbours[j]); k++ ) {
					if( dNeighbours[k] != pdNeighbours[k] )
						continue;

					for( m = 0; m < (int)(d->numNeighbours[j]); m++ ) {
						if( d->neighbours[j][k] == pd->neighbours[j][m] )
							if( fabs( d->voronoiSides[j][k] - pd->voronoiSides[j][m] ) > EPS )
								voronoiSideTest = 0;
					}
				}
			}

			if( d->numTriangles != pd->numTriangles )
				triangleCountTest = 0;

			Stg_Class_Delete( d );
		}

		Stg_Class_Delete( pd );
		
		if( rank == 0 ) {
			Memory_Free( sites );

			pcu_check_true( voronoiAreaTest );
			pcu_check_true( voronoiSideTest );
			pcu_check_true( triangleCountTest );
		}
	}
}

void ParallelDelaunaySuite( pcu_suite_t* suite ) {
   pcu_suite_setData( suite, ParallelDelaunaySuiteData );
   pcu_suite_setFixtures( suite, ParallelDelaunaySuite_Setup, ParallelDelaunaySuite_Teardown );
   pcu_suite_addTest( suite, ParallelDelaunaySuite_TestIrregular );
   pcu_suite_addTest( suite, ParallelDelaunaySuite_TestRegular );
}


