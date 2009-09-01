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
**   Tests the TrigMathSuite
**
** $Id: testTemplate.c 3462 2006-02-19 06:53:24Z WalterLandry $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>

#include "pcu/pcu.h"
#include <StGermain/StGermain.h> 
#include "StgDomain/Geometry/Geometry.h"
#include "StgDomain/Shape/Shape.h"
#include "StgDomain/Mesh/Mesh.h" 
#include "StgDomain/Utils/Utils.h"
#include "StgDomain/Swarm/Swarm.h"

#include "TrigMathSuite.h"

typedef struct {
	MPI_Comm	comm;
	unsigned rank;
	unsigned nProcs;
	double	angle;
	double	rectOriginal[3];
	double	spherical[3];
	double	rectangular[3];
	Index		dim;
} TrigMathSuiteData;

void TrigMathSuite_DomainFindingFunction( TrigMathSuiteData* data, double angle ) {
	char base[24];
	char stg[24];	

	data->angle = angle;
	sprintf( base, "%lf", sin(data->angle) );
	sprintf( stg, "%lf", sin( StGermain_TrigDomain(data->angle) ) );
	pcu_check_streq( base, stg );
}

void TrigMathSuite_Setup( TrigMathSuiteData* data ) {
	/* MPI Initializations */
	data->comm = MPI_COMM_WORLD;
	MPI_Comm_rank( data->comm, &data->rank );
	MPI_Comm_size( data->comm, &data->nProcs );

	data->rectOriginal[0] = 2.4; data->rectOriginal[1] = 5; data->rectOriginal[2] = -10;
}

void TrigMathSuite_Teardown( TrigMathSuiteData* data ) {
}

void TrigMathSuite_TestAngleConversionMacros( TrigMathSuiteData* data ) {
	char buffer[24];
	
	sprintf( buffer, "%2.3f", StGermain_DegreeToRadian( 37.0 + 3 ));
	pcu_check_streq( buffer, "0.698" );
	sprintf( buffer, "%2.3f", StGermain_RadianToDegree( 1.234 * 2 ));
	pcu_check_streq( buffer, "141.406" );
}

void TrigMathSuite_TestDomainFindingFunction( TrigMathSuiteData* data ) {
	TrigMathSuite_DomainFindingFunction( data, 1.5 );
	TrigMathSuite_DomainFindingFunction( data, -1.2 );
	TrigMathSuite_DomainFindingFunction( data, 20.0 );
}

void TrigMathSuite_TestQuadrantFindingFunction( TrigMathSuiteData* data ) {
	data->angle = 45;
	pcu_check_true( StGermain_TrigQuadrant( StGermain_DegreeToRadian(data->angle) ) == 0 );
	data->angle = 120;
	pcu_check_true( StGermain_TrigQuadrant( StGermain_DegreeToRadian(data->angle) ) == 1 );
	data->angle = 195;
	pcu_check_true( StGermain_TrigQuadrant( StGermain_DegreeToRadian(data->angle) ) == 2 );
	data->angle = 340;
	pcu_check_true( StGermain_TrigQuadrant( StGermain_DegreeToRadian(data->angle) ) == 3 );
	data->angle = 730;
	pcu_check_true( StGermain_TrigQuadrant( StGermain_DegreeToRadian(data->angle) ) == 0 );
	data->angle = -135;
	pcu_check_true( StGermain_TrigQuadrant( StGermain_DegreeToRadian(data->angle) ) == 2 );
}

void TrigMathSuite_TestCoordinateConversionFunction2D( TrigMathSuiteData* data ) {
	char x[24], y[24];

	data->dim = 2;
	sprintf( x, "%lf", data->rectOriginal[0] );
	sprintf( y, "%lf", data->rectOriginal[1] );
	pcu_check_streq( x, "2.400000" );
	pcu_check_streq( y, "5.000000" );
	StGermain_RectangularToSpherical( data->spherical, data->rectOriginal, data->dim );
	sprintf( x, "%lf", data->spherical[0] );
	sprintf( y, "%lf", data->spherical[1] );
	pcu_check_streq( x, "5.546170" );
	pcu_check_streq( y, "1.123276" );
	StGermain_SphericalToRectangular( data->rectangular, data->spherical,data-> dim );
	sprintf( x, "%lf", data->rectangular[0] );
	sprintf( y, "%lf", data->rectangular[1] );
	pcu_check_streq( x, "2.400000" );
	pcu_check_streq( y, "5.000000" );
}

void TrigMathSuite_TestCoordinateConversionFunction3D( TrigMathSuiteData* data ) {
	char x[24], y[24], z[24];

	data->dim = 3;
	sprintf( x, "%lf", data->rectOriginal[0] );
	sprintf( y, "%lf", data->rectOriginal[1] );
	sprintf( z, "%lf", data->rectOriginal[2] );
	pcu_check_streq( x, "2.400000" );
	pcu_check_streq( y, "5.000000" );
	pcu_check_streq( z, "-10.000000" );
	StGermain_RectangularToSpherical( data->spherical, data->rectOriginal, data->dim );
	sprintf( x, "%lf", data->spherical[0] );
	sprintf( y, "%lf", data->spherical[1] );
	sprintf( z, "%lf", data->spherical[2] );
	pcu_check_streq( x, "11.435034" );
	pcu_check_streq( y, "1.123276" );
	pcu_check_streq( z, "2.635212" );
	StGermain_SphericalToRectangular( data->rectangular, data->spherical,data-> dim );
	sprintf( x, "%lf", data->rectangular[0] );
	sprintf( y, "%lf", data->rectangular[1] );
	sprintf( z, "%lf", data->rectangular[2] );
	pcu_check_streq( x, "2.400000" );
	pcu_check_streq( y, "5.000000" );
	pcu_check_streq( z, "-10.000000" );
}

void TrigMathSuite( pcu_suite_t* suite ) {
	pcu_suite_setData( suite, TrigMathSuiteData );
   pcu_suite_setFixtures( suite, TrigMathSuite_Setup, TrigMathSuite_Teardown );
   pcu_suite_addTest( suite, TrigMathSuite_TestAngleConversionMacros );
   pcu_suite_addTest( suite, TrigMathSuite_TestDomainFindingFunction );
   pcu_suite_addTest( suite, TrigMathSuite_TestQuadrantFindingFunction );
   pcu_suite_addTest( suite, TrigMathSuite_TestCoordinateConversionFunction2D );
   pcu_suite_addTest( suite, TrigMathSuite_TestCoordinateConversionFunction3D );
}
