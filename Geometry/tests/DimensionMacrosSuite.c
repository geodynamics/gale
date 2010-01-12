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
**   Tests the DimensionMacrosSuite
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

#include "DimensionMacrosSuite.h"

typedef struct {
	MPI_Comm	comm;
	int		rank;
	int		nProcs;
} DimensionMacrosSuiteData;

void DimensionMacrosSuite_testDimensionMacros_DoOneTest( IJK coord, IJK meshSize, Stream* stream ) {
   IJK newCoord;
   Index element_I;

   Journal_Printf( stream,  "(%d,%d,%d) in mesh sized %d*%d*%d -> ",
      coord[0], coord[1], coord[2],
      meshSize[0], meshSize[1], meshSize[2] );

   Dimension_3DTo1D( coord, meshSize, &element_I );
   Journal_Printf( stream,  "index %d\n", element_I );
   Dimension_1DTo3D( element_I, meshSize, newCoord );
   Journal_Printf( stream,  "which maps back to (%d,%d,%d)\n",
      newCoord[0], newCoord[1], newCoord[2] );
}

void DimensionMacrosSuite_Setup( DimensionMacrosSuiteData* data ) {
	/* MPI Initializations */
	data->comm = MPI_COMM_WORLD;
	MPI_Comm_rank( data->comm, &data->rank );
	MPI_Comm_size( data->comm, &data->nProcs );
}

void DimensionMacrosSuite_Teardown( DimensionMacrosSuiteData* data ) {
}

void DimensionMacrosSuite_TestDimensionMacros( DimensionMacrosSuiteData* data ) {
	char		expected_file[PCU_PATH_MAX];
	Stream*  stream = Journal_Register( Info_Type, "DimensionMacrosStream" );
	IJK		coord;
   IJK		meshSize;
   Index		i, j, k;

	Stream_RedirectFile( stream, "testDimensionMacros.dat" );

	Journal_Printf( stream,  "+++ 1D Tests +++\n\n" );
	coord[I_AXIS] = 3; coord[J_AXIS] = 0; coord[K_AXIS] = 0;
	meshSize[I_AXIS] = 8; meshSize[J_AXIS] = 0; meshSize[K_AXIS] = 0;
	DimensionMacrosSuite_testDimensionMacros_DoOneTest( coord, meshSize, stream );
	
	coord[I_AXIS] = 0; coord[J_AXIS] = 3; coord[K_AXIS] = 0;
	meshSize[I_AXIS] = 0; meshSize[J_AXIS] = 8; meshSize[K_AXIS] = 0;
	DimensionMacrosSuite_testDimensionMacros_DoOneTest( coord, meshSize, stream );
	
	coord[I_AXIS] = 0; coord[J_AXIS] = 0; coord[K_AXIS] = 3;
	meshSize[I_AXIS] = 0; meshSize[J_AXIS] = 0; meshSize[K_AXIS] = 8;
	DimensionMacrosSuite_testDimensionMacros_DoOneTest( coord, meshSize, stream );

	Journal_Printf( stream,  "\n+++ 2D Tests +++\n\n" );
	coord[I_AXIS] = 3; coord[J_AXIS] = 4; coord[K_AXIS] = 0;
	meshSize[I_AXIS] = 8; meshSize[J_AXIS] = 8; meshSize[K_AXIS] = 0;
	DimensionMacrosSuite_testDimensionMacros_DoOneTest( coord, meshSize, stream );
	
	coord[I_AXIS] = 3; coord[J_AXIS] = 0; coord[K_AXIS] = 4;
	meshSize[I_AXIS] = 8; meshSize[J_AXIS] = 0; meshSize[K_AXIS] = 8;
	DimensionMacrosSuite_testDimensionMacros_DoOneTest( coord, meshSize, stream );
	
	coord[I_AXIS] = 0; coord[J_AXIS] = 3; coord[K_AXIS] = 4;
	meshSize[I_AXIS] = 0; meshSize[J_AXIS] = 8; meshSize[K_AXIS] = 8;
	DimensionMacrosSuite_testDimensionMacros_DoOneTest( coord, meshSize, stream );
	
	Journal_Printf( stream,  "\n+++ 3D Tests +++\n\n" );
	meshSize[I_AXIS] = 3; meshSize[J_AXIS] = 4; meshSize[K_AXIS] = 5;
	for ( k=0; k < meshSize[K_AXIS]; k++ ) {
		for ( j=0; j < meshSize[J_AXIS]; j++ ) {
			for ( i=0; i < meshSize[I_AXIS]; i++ ) {
				coord[I_AXIS] = i; coord[J_AXIS] = j; coord[K_AXIS] = k;
				DimensionMacrosSuite_testDimensionMacros_DoOneTest( coord, meshSize, stream );
			}
		}
	}	

	pcu_filename_expected( "testDimensionMacros.expected", expected_file );
	pcu_check_fileEq( "testDimensionMacros.dat", expected_file );
	remove( "testDimensionMacros.dat" );

	Stream_CloseAndFreeFile( stream );
}

void DimensionMacrosSuite( pcu_suite_t* suite ) {
	pcu_suite_setData( suite, DimensionMacrosSuiteData );
   pcu_suite_setFixtures( suite, DimensionMacrosSuite_Setup, DimensionMacrosSuite_Teardown );
   pcu_suite_addTest( suite, DimensionMacrosSuite_TestDimensionMacros );
}


