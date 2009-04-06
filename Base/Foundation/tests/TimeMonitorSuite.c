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
** $Id: testNamedStg_TimeMonitor.c 2432 2005-08-08 23:01:59Z Raquibul Hassan $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <unistd.h>
#include <stdio.h>
#include <mpi.h>

#include "pcu/pcu.h"
#include "StGermain/Base/Foundation/Foundation.h"
#include "StGermain/Base/Foundation/forwardDecl.h" /* For Journal stuff */
#include "TimeMonitorSuite.h"

typedef struct {
} TimeMonitorSuiteData;

void TimeMonitorSuite_Setup( TimeMonitorSuiteData* data ) {
}

void TimeMonitorSuite_Teardown( TimeMonitorSuiteData* data ) {
}

void TimeMonitorSuite_TestTimingPeriod( TimeMonitorSuiteData* data ) {
   Stream*         timeMonitorStream = NULL;
   char*           timeMonitorOutputFilename = "./TimeMonitorSuite_TestOutput.txt";
   FILE*           timeMonitorOutputFile = NULL;
   char*           infoString = NULL;

   infoString = malloc( sizeof(char*) * 1000 );
   
   Stg_TimeMonitor_SetTimerWatchCriteria( 0.5 );

   /* Calling the TimeMonitor functions below will cause info about the
    * time period in question to
    * be printed to the info stream called Stg_TimeMonitor_InfoStreamName.
    * We need to redirect it so we can compare the output.
    * Note: ideally would redirect to a char* buffer, but the StG journal doesn't seem
    * to have this functionality yet */
   timeMonitorStream = Journal_Register( Info_Type, Stg_TimeMonitor_InfoStreamName );
	Stream_RedirectFile( timeMonitorStream, timeMonitorOutputFilename );
   Journal_Enable_TypedStream( Info_Type, True );
	Stream_Enable( timeMonitorStream, True );

   Stg_TimeMonitor* tm = Stg_TimeMonitor_New( "test", True, True, MPI_COMM_WORLD );
   sleep( 1 );
   Stg_TimeMonitor_Begin( tm );
   sleep( 2 );
   Stg_TimeMonitor_End( tm );

   /* Now we have to do some good ole' C file manipulation to get the output back */
   timeMonitorOutputFile = fopen(timeMonitorOutputFilename, "r");
   fgets( infoString, sizeof(char*)*1000, timeMonitorOutputFile );
   fclose( timeMonitorOutputFile );
   remove( timeMonitorOutputFilename );

   /* Note: I don't like this test approach, but it seems the only way without changing the
    * Component significantly to return TimeMonitor info by methods, rather than just printing
    * -- PatrickSunter, 3 Apr 2009 */
   pcu_check_true( 0 == strcmp( "\tTimeMonitor(test):  ts: 3 (secs), dt(67%): 2s\n", infoString ) );
   Stg_TimeMonitor_Delete( tm );
   free( infoString );
}


void TimeMonitorSuite( pcu_suite_t* suite ) {
   pcu_suite_setData( suite, TimeMonitorSuiteData );
   pcu_suite_setFixtures( suite, TimeMonitorSuite_Setup, TimeMonitorSuite_Teardown );
   pcu_suite_addTest( suite, TimeMonitorSuite_TestTimingPeriod );
}
