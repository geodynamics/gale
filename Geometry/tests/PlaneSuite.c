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
**   Tests the PlaneSuite
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

#include "PlaneSuite.h"

typedef struct {
	MPI_Comm comm;
	unsigned rank;
	unsigned nProcs;
} PlaneSuiteData;

void PlaneSuite_Setup( PlaneSuiteData* data ) {
	/* MPI Initializations */
	data->comm = MPI_COMM_WORLD;
	MPI_Comm_rank( data->comm, &data->rank );
	MPI_Comm_size( data->comm, &data->nProcs );
}

void PlaneSuite_Teardown( PlaneSuiteData* data ) {
}

void PlaneSuite_TestDistanceToPoint( PlaneSuiteData* data ) {
	unsigned procToWatch = data->nProcs >=2 ? 1 : 0;
	
   if( data->rank == procToWatch ) {
		Coord axisA = { 1.0, 0.3, 0.0 };
		Coord axisB = { -0.3, 1.0, 0.0 };
		Coord point1 = { 0.0, 0.0, 2.0 };
		Coord point2 = { 20.0, 100.0, 6.0 };
		Plane a;

		Plane_CalcFromVec( a, axisA, axisB, point1 );
		pcu_check_true( a[0] == 0 && a[1] == -0 && a[2] == 1 && a[3] == 2 );
		pcu_check_true( Plane_DistanceToPoint( a, point2 ) == 4 );
	}
}

void PlaneSuite_TestPointIsInFront( PlaneSuiteData* data ) {
	unsigned procToWatch = data->nProcs >=2 ? 1 : 0;
	
   if( data->rank == procToWatch ) {
		Coord axisA = { 1.0, 0.3, 0.0 };
		Coord axisB = { -0.3, 1.0, 0.0 };
		Coord point1 = { 0.0, 0.0, 2.0 };
		Coord point2 = { 20.0, 100.0, 6.0 };
		Coord point3 = { 20.0, 100.0, -60.0 };
		Plane a;

		Plane_CalcFromVec( a, axisA, axisB, point1 );
		pcu_check_true( a[0] == 0 && a[1] == -0 && a[2] == 1 && a[3] == 2 );
		pcu_check_true( Plane_PointIsInFront( a, point2 ) );
		pcu_check_true( !Plane_PointIsInFront( a, point3 ) );
	}
}

void PlaneSuite( pcu_suite_t* suite ) {
	pcu_suite_setData( suite, PlaneSuiteData );
   pcu_suite_setFixtures( suite, PlaneSuite_Setup, PlaneSuite_Teardown );
   pcu_suite_addTest( suite, PlaneSuite_TestDistanceToPoint );
   pcu_suite_addTest( suite, PlaneSuite_TestPointIsInFront );
}


