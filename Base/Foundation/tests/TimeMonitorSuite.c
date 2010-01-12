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
#include <math.h>
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
   #define           MAXLINE 1000
   TimeMonitorData   tmData;
   Stg_TimeMonitor*  tm=NULL;
   double            percentOfTotalCalc;

   Stg_TimeMonitor_SetTimerWatchCriteria( 0.5 );

   tm = Stg_TimeMonitor_New( "test", True, False /*Don't print*/, MPI_COMM_WORLD );
   sleep( 1 );
   Stg_TimeMonitor_Begin( tm );
   sleep( 2 );
   Stg_TimeMonitor_End( tm, &tmData );

   pcu_check_true( ( 2.95 < tmData.totalSinceInit ) && ( tmData.totalSinceInit < 3.05 ) );
   pcu_check_true( ( 1.95 < tmData.dt ) && ( tmData.dt < 2.05 ) );
   pcu_check_true( ( 1.95 < tmData.aveProcDt ) && ( tmData.aveProcDt < 2.05 ) );
   percentOfTotalCalc = tmData.aveProcDt / tmData.totalSinceInit * 100;
   pcu_check_true( fabs( percentOfTotalCalc - tmData.percentTM_ofTotal ) < 0.01 );
   pcu_check_true( tmData.criterionPassed == True );

   Stg_TimeMonitor_Delete( tm );
}

void TimeMonitorSuite( pcu_suite_t* suite ) {
   pcu_suite_setData( suite, TimeMonitorSuiteData );
   pcu_suite_setFixtures( suite, TimeMonitorSuite_Setup, TimeMonitorSuite_Teardown );
   pcu_suite_addTest( suite, TimeMonitorSuite_TestTimingPeriod );
}


