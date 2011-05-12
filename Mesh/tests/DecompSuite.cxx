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
**   Tests the DecompSuite
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

#include "DecompSuite.h"

typedef struct {
	MPI_Comm	comm;
	int		rank;
	int		nProcs;
} DecompSuiteData;

void DecompSuite_Setup( DecompSuiteData* data ) {
	/* MPI Initializations */
	data->comm = MPI_COMM_WORLD;
	MPI_Comm_rank( data->comm, &data->rank );
	MPI_Comm_size( data->comm, &data->nProcs );
}

void DecompSuite_Teardown( DecompSuiteData* data ) {
}

void DecompSuite_TestDecomp( DecompSuiteData* data ) {
	Decomp*	decomp;
	int		nLocs, *locs, *ranks;
   int		l_i, g_i;

	nLocs = 10;
	locs = MemArray( int, nLocs, "testDecomp" );

	decomp = Decomp_New();
   for( l_i = 0; l_i < nLocs; l_i++ )
		locs[l_i] = data->rank * nLocs + l_i;
	pcu_check_noassert( Decomp_SetLocals( decomp, nLocs, locs ) );
	for( g_i = 0; g_i < data->nProcs * nLocs; g_i++ ) {
		if( g_i >= data->rank * nLocs && g_i < (data->rank + 1) * nLocs ) {
			pcu_check_true( IMap_Map( decomp->owners, g_i ) == data->rank );
		}
		else {
			pcu_check_true( !IMap_Has( decomp->owners, g_i ) );
		}
	}

	for( l_i = 0; l_i < nLocs; l_i++ ) {
		locs[l_i] = (data->rank * nLocs + nLocs / 2 + l_i) % (data->nProcs * nLocs);
   }
	pcu_check_noassert( Decomp_SetLocals( decomp, nLocs, locs ) );
	for( g_i = 0; g_i < data->nProcs * nLocs; g_i++ ) {
		if( g_i >= data->rank * nLocs && g_i < (data->rank + 1) * nLocs ) {
			if( g_i < data->rank * nLocs + nLocs / 2 ) {
				if( data->rank > 0 ) {
					pcu_check_true( IMap_Map( decomp->owners, g_i ) == data->rank - 1 );
				}
				else {
					pcu_check_true( IMap_Map( decomp->owners, g_i ) == data->nProcs - 1 );
				}
			}
			else {
				pcu_check_true( IMap_Map( decomp->owners, g_i ) == data->rank );
			}
		}
		else {
			pcu_check_true( !IMap_Has( decomp->owners, g_i ) );
		}
	}

	locs = MemRearray( locs, int, data->nProcs * nLocs, "testDecomp" );
	ranks = MemArray( int, data->nProcs * nLocs, "testDecomp" );
	for( g_i = 0; g_i < data->nProcs * nLocs; g_i++ )
		locs[g_i] = g_i;
	pcu_check_noassert( Decomp_FindOwners( decomp, data->nProcs * nLocs, locs, ranks ) );

	NewClass_Delete( decomp );
	MemFree( locs );
	MemFree( ranks );
}

void DecompSuite( pcu_suite_t* suite ) {
	pcu_suite_setData( suite, DecompSuiteData );
   pcu_suite_setFixtures( suite, DecompSuite_Setup, DecompSuite_Teardown );
   pcu_suite_addTest( suite, DecompSuite_TestDecomp );
}


