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
**   Tests the SingleCellLayoutSuite
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

#include "SingleCellLayoutSuite.h"

typedef struct {
	SingleCellLayout*	singleCellLayout;
	unsigned int		dimExists[3];
	double***			cellPoints;
	MPI_Comm				comm;
	unsigned int		rank;
	unsigned int		nProcs;
} SingleCellLayoutSuiteData;


void SingleCellLayoutSuite_Setup( SingleCellLayoutSuiteData* data ) {
	/* MPI Initializations */
	data->comm = MPI_COMM_WORLD;  
	MPI_Comm_rank( data->comm, &data->rank );
	MPI_Comm_size( data->comm, &data->nProcs );

	data->dimExists[0] = True; data->dimExists[1] = True; data->dimExists[2] = True;
	/* Configure the single-cell-layout */
	data->singleCellLayout = SingleCellLayout_New( "singleCellLayout", NULL, data->dimExists, NULL, NULL );
}

void SingleCellLayoutSuite_Teardown( SingleCellLayoutSuiteData* data ) {
}

void SingleCellLayoutSuite_Driver( SingleCellLayoutSuiteData* data ) {
	Cell_Index  cell;

	for( cell = 0; cell < data->singleCellLayout->_cellLocalCount( data->singleCellLayout ); cell++ ) {
		Cell_PointIndex point;
		Cell_PointIndex count;

		count = data->singleCellLayout->_pointCount( data->singleCellLayout, cell );
		data->cellPoints = Memory_Alloc_Array( double**, count, "cellsPoints" );
		data->singleCellLayout->_initialisePoints( data->singleCellLayout, cell, count, data->cellPoints );
	}
}

void SingleCellLayoutSuite_TestMapElement( SingleCellLayoutSuiteData* data ) {
	int procToWatch;

	procToWatch = data->nProcs >=2 ? 1 : 0;

	if( data->rank == procToWatch ) {
		SingleCellLayoutSuite_Driver( data );
		
		pcu_check_true( CellLayout_MapElementIdToCellId( data->singleCellLayout, 0 ) == 0 );
      pcu_check_true( CellLayout_MapElementIdToCellId( data->singleCellLayout, 5 ) == 0 );
      pcu_check_true( CellLayout_MapElementIdToCellId( data->singleCellLayout, 100 ) == 0 );
      pcu_check_true( CellLayout_CellOf( data->singleCellLayout, data->cellPoints[0] ) == 0 );

		Memory_Free( data->cellPoints );
		Stg_Class_Delete( data->singleCellLayout );	
	}
}

void SingleCellLayoutSuite_TestIsInCell( SingleCellLayoutSuiteData* data ) {
	int		procToWatch;
	double*	testCoord;		

	procToWatch = data->nProcs >=2 ? 1 : 0;

	if( data->rank == procToWatch ) {
		SingleCellLayoutSuite_Driver( data );
		testCoord = Memory_Alloc_Array_Unnamed( double, 3 );
		
		testCoord[0] = testCoord[1] = testCoord[2] = 0;
		pcu_check_true( CellLayout_IsInCell( data->singleCellLayout, 0, &testCoord ) );
		testCoord[0] = testCoord[1] = testCoord[2] = 1;
		pcu_check_true( CellLayout_IsInCell( data->singleCellLayout, 0, &testCoord ) );
		testCoord[0] = testCoord[1] = testCoord[2] = 2;
		pcu_check_true( !CellLayout_IsInCell( data->singleCellLayout, 0, &testCoord ) );
	
		Memory_Free( testCoord );
		Memory_Free( data->cellPoints );
		Stg_Class_Delete( data->singleCellLayout );	
	}
}

void SingleCellLayoutSuite( pcu_suite_t* suite ) {
	pcu_suite_setData( suite, SingleCellLayoutSuiteData );
	pcu_suite_setFixtures( suite, SingleCellLayoutSuite_Setup, SingleCellLayoutSuite_Teardown );
	pcu_suite_addTest( suite, SingleCellLayoutSuite_TestMapElement );
	pcu_suite_addTest( suite, SingleCellLayoutSuite_TestIsInCell );
}
