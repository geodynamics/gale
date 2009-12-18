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
**   Tests the DelaunaySuite
**
** $Id: testDelaunay.c 3462 2006-02-19 06:53:24Z WalterLandry $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>

#include "pcu/pcu.h"
#include <StGermain/StGermain.h>
#include "StgDomain/Geometry/Geometry.h"
#include "DelaunaySuite.h"

#define PI 3.1415926535897932384626
#define epsilon 0.001

typedef enum pointsType_t {
	Irregular,
	Regular,
	Polygon
} pointsType;

typedef struct {
	DelaunayAttributes 	attr;
	Dictionary*		dict;
} DelaunaySuiteData;

void GeneratePoints( CoordF* sites, int numSites, pointsType *p ) {
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

void DelaunaySuite_Setup( DelaunaySuiteData* data ) {
	data->attr.BuildBoundingTriangle = 0;
	data->attr.BuildTriangleIndices = 1;
	data->attr.BuildTriangleNeighbours = 1;
	data->attr.CreateVoronoiVertices = 1;
	data->attr.CalculateVoronoiSides = 1;
	data->attr.CalculateVoronoiSurfaceArea = 1;
	data->attr.FindNeighbours = 1;

	data->dict = Dictionary_New();
}

void DelaunaySuite_Teardown( DelaunaySuiteData* data ) {
	Stg_Class_Delete( data->dict );
}

void DelaunaySuite_TriangulateRegularMesh( DelaunaySuiteData* data ) {
	Delaunay*	delaunay;
	pointsType	p		= Regular;
	CoordF*		sites;
	int		array[10]	= { 16, 25, 36, 49, 64, 81, 100, 121, 144, 269 };
	float		sqroot, area;
	int 		i, j;
	Bool		pass		= True;

	for( i = 0; i < 10; i++ ) {
		sites = Memory_Alloc_Array_Unnamed( CoordF, array[i] );
		memset( sites, 0, sizeof(CoordF)*array[i] );

		GeneratePoints( sites, array[i], &p );

		delaunay = Delaunay_New( "Delaunay-Regular", data->dict, sites, array[i], 0, &data->attr );
		Stg_Component_Build( delaunay, NULL, True );

		sqroot = sqrt((float)array[i]);
		area = (1.0/sqroot)*(1.0/sqroot);
		
		for( j = 0; j < array[i]; j++ ) {
			if( delaunay->hull[j] )
				continue;

			if( j > 0 && j < sqrt((float)j) )
				if( FABS( area - delaunay->voronoiArea[j] ) > epsilon )
					pass = False;
		}

		Stg_Class_Delete( delaunay );
		Memory_Free( sites );
	}

	pcu_check_true( pass );
}

void DelaunaySuite_TriangulateIrregularMesh( DelaunaySuiteData* data ) {
	Delaunay*	delaunay;
	pointsType	p		= Irregular;
	CoordF*		sites;
	int 		pass		= True;
	int 		i;
		
	for( i = 10; i < 210; i++ ) {
		sites = Memory_Alloc_Array_Unnamed( CoordF, i );
		memset( sites, 0, sizeof( CoordF )*i );
		
		GeneratePoints( sites, i, &p );
	
		delaunay = Delaunay_New( "Delaunay-Regular", data->dict, sites, i, 0, &data->attr );
		Stg_Component_Build( delaunay, NULL, True );

		if( delaunay->numFaces != delaunay->numEdges - delaunay->numSites + 2 )
			pass = False;

		Stg_Class_Delete( delaunay );
		Memory_Free( sites );
	}

	pcu_check_true( pass );
}

void DelaunaySuite_TriangulatePolygonMesh( DelaunaySuiteData* data ) {
	Delaunay*	delaunay;
	CoordF*		sites;
	char		filename[PCU_PATH_MAX];
	FILE*		fp;
	float		theta;
	float		voronoiArea;
	float		a, area, side;
	int 		i, j;
	int		pass		= True;
	
	pcu_filename_input( "small.txt", filename );
	fp = fopen( filename, "r+" );

	fscanf( fp, "%d", &i );

	theta = 2.0*PI/((float)(i-1));

	sites = Memory_Alloc_Array_Unnamed( CoordF, i );
	memset( sites, 0, sizeof(CoordF)*i );

	j = 0;
	while( fscanf( fp, "%f %f", sites[j], sites[j] + 1 ) != EOF ) {
		j++;
	}

	delaunay = Delaunay_New( "Delaunay-Polygon", data->dict, sites, i, 0, &data->attr );
	Stg_Component_Build( delaunay, NULL, True );

	for( j = 0; j < i; j++ ) {
		if( fabs( sites[j][0] ) < 0.00001 && fabs( sites[j][1] ) < 0.00001 ) {
			side = delaunay->voronoiSides[j][0];
			voronoiArea = delaunay->voronoiArea[j];
			break;
		}
	}

	a = ( (0.5*side)/sin(0.5*theta) );
	area = (0.5*a*a*sin(theta))*((float)(i-1));

	if( FABS( area - voronoiArea ) > epsilon )
		pass = False;

	fclose( fp );

	Stg_Class_Delete( delaunay );
	Memory_Free( sites );

	pcu_check_true( pass );
}

void DelaunaySuite( pcu_suite_t* suite ) {
   pcu_suite_setData( suite, DelaunaySuiteData );
   pcu_suite_setFixtures( suite, DelaunaySuite_Setup, DelaunaySuite_Teardown );
   pcu_suite_addTest( suite, DelaunaySuite_TriangulateRegularMesh );
   pcu_suite_addTest( suite, DelaunaySuite_TriangulateIrregularMesh );
   pcu_suite_addTest( suite, DelaunaySuite_TriangulatePolygonMesh );
}


