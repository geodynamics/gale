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
**   Tests the UpwindXiSuite
**
** $Id: testUpwindXi.c 3462 2006-02-19 06:53:24Z WalterLandry $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>

#include "pcu/pcu.h"
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include "StgFEM/Discretisation/Discretisation.h"
#include "StgFEM/SLE/SystemSetup/SystemSetup.h"
#include "StgFEM/SLE/ProvidedSystems/AdvectionDiffusion/AdvectionDiffusion.h"
#include "UpwindXiSuite.h"

typedef struct {
} UpwindXiSuiteData;


void UpwindXiSuite_Setup( UpwindXiSuiteData* data ) {
}

void UpwindXiSuite_Teardown( UpwindXiSuiteData* data ) {
}


void UpwindXiSuite_Test( UpwindXiSuiteData* data ) {
	Stream* dataStream;
	double  minPercletNumber = -15.0;
	double  maxPercletNumber = 15.0;
	double  dPerceltNumber = 0.5;
	double  perceltNumber;
	char	expectedFile[PCU_PATH_MAX];
	int	rank;
	char	outputFilename[100];

	MPI_Comm_rank( MPI_COMM_WORLD, &rank );
	sprintf( outputFilename, "output_%2d.dat", rank );

	dataStream    = Journal_Register( Info_Type, "DataStream" );
	Stream_RedirectFile( dataStream, outputFilename );
		
	Journal_Printf( dataStream, "# File to compare code with Brooks, Hughes 1982 - Fig 3.3\n");
	Journal_Printf( dataStream, "# Integration rule for the optimal upwind scheme, doubly asymptotic approximation and critical approximation.\n");
	Journal_Printf( dataStream, "# Plot using line:\n");
	Journal_Printf( dataStream, "# plot \"%s\" using 1:2 title \"Exact\" with line, ", "output.dat" );
	Journal_Printf( dataStream, "\"%s\" using 1:3 title \"DoublyAsymptoticAssumption\" with line, ", "output.dat" );
	Journal_Printf( dataStream, "\"%s\" using 1:4 title \"CriticalAssumption\" with line\n", "output.dat" );

	Journal_Printf( dataStream, "# Perclet Number \t Exact \t DoublyAsymptoticAssumption \t CriticalAssumption\n" );
	for ( perceltNumber = minPercletNumber ; perceltNumber < maxPercletNumber ; perceltNumber += dPerceltNumber )
		Journal_Printf( dataStream, "%0.3g \t\t %0.3g \t\t %0.3g \t\t %0.3g\n", 
				perceltNumber, 
				AdvDiffResidualForceTerm_UpwindXiExact( NULL, perceltNumber),
				AdvDiffResidualForceTerm_UpwindXiDoublyAsymptoticAssumption( NULL, perceltNumber),
				AdvDiffResidualForceTerm_UpwindXiCriticalAssumption( NULL, perceltNumber) );

	pcu_filename_expected( "UpwindXi.expected", expectedFile );
	pcu_check_fileEq( outputFilename, expectedFile );
	remove( outputFilename );
	/*Journal_Purge( );*/
}

void UpwindXiSuite( pcu_suite_t* suite ) {
   pcu_suite_setData( suite, UpwindXiSuiteData );
   pcu_suite_setFixtures( suite, UpwindXiSuite_Setup, UpwindXiSuite_Teardown );
   pcu_suite_addTest( suite, UpwindXiSuite_Test );
}



