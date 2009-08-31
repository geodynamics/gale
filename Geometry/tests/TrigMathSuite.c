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
	MPI_Comm			comm;
	unsigned int	rank;
	unsigned int	nProcs;
} TrigMathSuiteData;

void TrigMathSuite_Setup( TrigMathSuiteData* data ) {
	/* MPI Initializations */
	data->comm = MPI_COMM_WORLD;
	MPI_Comm_rank( data->comm, &data->rank );
	MPI_Comm_size( data->comm, &data->nProcs );
}

void TrigMathSuite_Teardown( TrigMathSuiteData* data ) {
}

void TrigMathSuite_TestTrigMath( TrigMathSuiteData* data ) {
	unsigned	procToWatch;
	procToWatch = data->nProcs >=2 ? 1 : 0;
	
	if (data->rank == procToWatch) {
		double	angle;
		double	rectOriginal[] = {2.4,5,-10};
		double	spherical[3];
		double	rectangular[3];
		Index		dim;
		char		expected_file[PCU_PATH_MAX];
		Stream*	stream = Journal_Register( Info_Type, "TrigMathStream" );
		Stream_RedirectFile( stream, "testTrigMath.dat" );

		Journal_Printf( stream, "****************************\n");
		Journal_Printf( stream, "Test angle conversion macros\n");
		Journal_Printf( stream, "39 degrees in radians = %2.3f\n", StGermain_DegreeToRadian( 37.0 + 3 ) );
		Journal_Printf( stream, "2.468 radians in degrees = %2.3f\n", StGermain_RadianToDegree( 1.234 * 2 ) );

		Journal_Printf( stream, "\n****************************\n");
		Journal_Printf( stream, "Test domain finding function\n");
		angle = 1.5;
		Journal_Printf( stream, "Angle %lf is equivalent to %lf\n", angle, StGermain_TrigDomain(angle) );
		Journal_Printf( stream, "Sine test: %lf = %lf\n", sin(angle), sin( StGermain_TrigDomain(angle) ) );
		angle = -1.2;
		Journal_Printf( stream, "Angle %lf is equivalent to %lf\n", angle, StGermain_TrigDomain(angle) );
		Journal_Printf( stream, "Sine test: %lf = %lf\n", sin(angle), sin( StGermain_TrigDomain(angle) ) );
		angle = 20.0;
		Journal_Printf( stream, "Angle %lf is equivalent to %lf\n", angle, StGermain_TrigDomain(angle) );
		Journal_Printf( stream, "Sine test: %lf = %lf\n", sin(angle), sin( StGermain_TrigDomain(angle) ) );

		Journal_Printf( stream, "\n****************************\n");
		Journal_Printf( stream, "Test quadrant finding function\n");
		angle = 45;
		Journal_Printf( stream, "Angle %lf degrees is in quadrant %d\n", angle, StGermain_TrigQuadrant( StGermain_DegreeToRadian(angle) ));
		angle = 120;
		Journal_Printf( stream, "Angle %lf degrees is in quadrant %d\n", angle, StGermain_TrigQuadrant( StGermain_DegreeToRadian(angle) ));
		angle = 195;
		Journal_Printf( stream, "Angle %lf degrees is in quadrant %d\n", angle, StGermain_TrigQuadrant( StGermain_DegreeToRadian(angle) ));
		angle = 340;
		Journal_Printf( stream, "Angle %lf degrees is in quadrant %d\n", angle, StGermain_TrigQuadrant( StGermain_DegreeToRadian(angle) ));
		angle = 730;
		Journal_Printf( stream, "Angle %lf degrees is in quadrant %d\n", angle, StGermain_TrigQuadrant( StGermain_DegreeToRadian(angle) ));
		angle = -135;
		Journal_Printf( stream, "Angle %lf degrees is in quadrant %d\n", angle, StGermain_TrigQuadrant( StGermain_DegreeToRadian(angle) ));

		Journal_Printf( stream, "\n****************************\n");
		Journal_Printf( stream, "Test coordinate conversion functions 2D\n");
		dim = 2;
		StGermain_PrintNamedVector( stream, rectOriginal, dim );
		StGermain_RectangularToSpherical( spherical, rectOriginal, dim );
		StGermain_PrintNamedVector( stream, spherical, dim );
		StGermain_SphericalToRectangular( rectangular, spherical, dim );
		StGermain_PrintNamedVector( stream, rectangular, dim );

		Journal_Printf( stream, "\n****************************\n");
		Journal_Printf( stream, "Test coordinate conversion functions 3D\n");
		dim = 3;
		StGermain_PrintNamedVector( stream, rectOriginal, dim );
		StGermain_RectangularToSpherical( spherical, rectOriginal, dim );
		StGermain_PrintNamedVector( stream, spherical, dim );
		StGermain_SphericalToRectangular( rectangular, spherical, dim );
		StGermain_PrintNamedVector( stream, rectangular, dim );

		pcu_filename_expected( "testTrigMath.expected", expected_file );
		pcu_check_fileEq( "testTrigMath.dat", expected_file );
		remove( "testTrigMath.dat" );
	}
}

void TrigMathSuite( pcu_suite_t* suite ) {
	pcu_suite_setData( suite, TrigMathSuiteData );
   pcu_suite_setFixtures( suite, TrigMathSuite_Setup, TrigMathSuite_Teardown );
   pcu_suite_addTest( suite, TrigMathSuite_TestTrigMath );
}
