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
**   Tests the SobolGeneratorSuite
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

#include "SobolGeneratorSuite.h"

typedef struct {
	MPI_Comm			comm;
	unsigned			rank;
	unsigned			nProcs;
} SobolGeneratorSuiteData;

void SobolGeneratorSuite_Setup( SobolGeneratorSuiteData* data ) {
	/* MPI Initializations */
	data->comm = MPI_COMM_WORLD;
	MPI_Comm_rank( data->comm, &data->rank );
	MPI_Comm_size( data->comm, &data->nProcs );
}

void SobolGeneratorSuite_Teardown( SobolGeneratorSuiteData* data ) {

	/* Purge output files */
	remove( "RightmostBit.dat" );
	remove( "testGenerator0" );
	remove( "testGenerator1" );
	remove( "testGenerator2" );
	remove( "testGenerator3" );
	remove( "testGenerator4" );
	remove( "testGenerator5" );
	
}

void SobolGeneratorSuite_TestSobolGenerator( SobolGeneratorSuiteData* data ) {
	Stream*		stream;
	SobolGenerator*	sobolGenerator;
	Index		index;
	Index		sobol_I;
	int		bit_I;
	char		name[50];
	double		result;
	char		rightmostBit_expected[PCU_PATH_MAX];
	char		generator0_expected[PCU_PATH_MAX];
	char		generator1_expected[PCU_PATH_MAX];
	char		generator2_expected[PCU_PATH_MAX];
	char		generator3_expected[PCU_PATH_MAX];
	char		generator4_expected[PCU_PATH_MAX];
	char		generator5_expected[PCU_PATH_MAX];

	stream = Journal_Register (Info_Type, "SobolGeneratorStream");
	Stream_RedirectFile( stream, "RightmostBit.dat" );
	
        Journal_Printf( stream, " *********************** Testing _SobolGenerator_FindRightmostZeroBit *******************\n" );
        for ( index = 0 ; index < 30 ; index++ ) {
                for ( bit_I = sizeof( Index ) * 4 - 1 ; bit_I >= 0 ; bit_I-- )
                        Journal_Printf( stream, "%u", index & 1 << bit_I ? 1 : 0 );
                Journal_Printf( stream, " number %u: %u\n", index, _SobolGenerator_FindRightmostZeroBit( index ) );
        }

        /* constructor  */
        for ( sobol_I = 0 ; sobol_I < 6 ; sobol_I++ ) {
                sprintf( name, "testGenerator%u", sobol_I );
                Stream_RedirectFile( stream, name );
                sobolGenerator = SobolGenerator_NewFromTable( name );

                Journal_Printf( stream," ****************** Testing SobolGenerator_GetDirectionalNumber ***************\n" );
                for ( index = 0 ; index < 30 ; index++ )
                        SobolGenerator_GetDirectionalNumber( sobolGenerator, index );

                /* Checking up to 200000 numbers - this number is arbitary - 
                 * it's only limited because we don't want file size to be huge
                 * This number is intentionally over 25535 = 2^16 - 1 because there was a time when numbers repeated after this */
                for ( index = 0 ; index < 200000 ; index++ ) {
                        result = SobolGenerator_GetNextNumber(sobolGenerator);

                        assert( fabs( result - SobolGenerator_GetNumberByIndex(sobolGenerator, index)) < 1e-8 );

                        /* Only dump subset of data - this output criterion is completely arbitary */
                        if ( index % 773 == 3 )
                                Journal_Printf( stream, "%.4g\n", result );
                }
                Stg_Class_Delete( sobolGenerator );
	}
	pcu_filename_expected( "testSobolGeneratorRightmostBit.expected", rightmostBit_expected );
	pcu_check_fileEq( "RightmostBit.dat", rightmostBit_expected );
	pcu_filename_expected( "testSobolGenerator0.expected", generator0_expected );
	pcu_check_fileEq( "testGenerator0", generator0_expected );
	pcu_filename_expected( "testSobolGenerator1.expected", generator1_expected );
	pcu_check_fileEq( "testGenerator1", generator1_expected );	
	pcu_filename_expected( "testSobolGenerator2.expected", generator2_expected );
	pcu_check_fileEq( "testGenerator2", generator2_expected );
	pcu_filename_expected( "testSobolGenerator3.expected", generator3_expected );
	pcu_check_fileEq( "testGenerator3", generator3_expected );
	pcu_filename_expected( "testSobolGenerator4.expected", generator4_expected );
	pcu_check_fileEq( "testGenerator4", generator4_expected );
	pcu_filename_expected( "testSobolGenerator5.expected", generator5_expected );
	pcu_check_fileEq( "testGenerator5", generator5_expected );
}

void SobolGeneratorSuite( pcu_suite_t* suite ) {
	pcu_suite_setData( suite, SobolGeneratorSuiteData );
	pcu_suite_setFixtures( suite, SobolGeneratorSuite_Setup, SobolGeneratorSuite_Teardown );
	pcu_suite_addTest( suite, SobolGeneratorSuite_TestSobolGenerator );
}
