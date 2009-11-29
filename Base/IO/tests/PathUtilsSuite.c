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
**   Tests the PathUtilsSuite
**
** $Id: testPathUtils.c 3462 2006-02-19 06:53:24Z WalterLandry $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#include "pcu/pcu.h"
#include "StGermain/Base/Foundation/Foundation.h"
#include "StGermain/Base/IO/IO.h"
#include "PathUtilsSuite.h"

typedef struct {
   unsigned int   rank;
} PathUtilsSuiteData;


void PathUtilsSuite_Setup( PathUtilsSuiteData* data ) {
   MPI_Comm_rank( MPI_COMM_WORLD, &data->rank );
}

void PathUtilsSuite_Teardown( PathUtilsSuiteData* data ) {
}


void PathUtilsSuite_TestPathJoin( PathUtilsSuiteData* data ) {
   char      searchPaths[1024];
   unsigned  len;
   
   PathJoin( searchPaths, 1, "." );
   len = strlen( searchPaths );
   searchPaths[len] = ':';
   
   PathJoin( &searchPaths[len + 1], 2,  ".", "data" );
   len = strlen( searchPaths );
   searchPaths[len] = ':';
   
   PathJoin( &searchPaths[len + 1], 6, "..", "..", "..", "does", "not", "exist" );
   
   pcu_check_streq( searchPaths, ".:./data:../../../does/not/exist" );
}


void PathUtilsSuite_TestFindFile( PathUtilsSuiteData* data ) {
   char*       searchPaths = NULL;
   char        fullPath[1024];
   const char* subDir = "./testSubDir";
   const char* subDirFilename = "./testSubDir/subDirTest.xml";
   const char* currDirFilename = "./currDirTest.xml";


   Stg_asprintf( &searchPaths, ".:%s:/does/not/exist", subDir );
   /* Create necessary test files/dirs */
   if (data->rank==0) {
      FILE*       subDirFile = NULL;
      FILE*       currDirFile = NULL;

      currDirFile = fopen( currDirFilename, "w" );
      fputs( "test.\n", currDirFile );
      fclose( currDirFile );
      mkdir( subDir, 0755 );
      subDirFile = fopen( subDirFilename, "w" );
      fputs( "test.\n", subDirFile );
      fclose( subDirFile );
   }
   MPI_Barrier(MPI_COMM_WORLD);

   /* try and open some files using the search path */
   /* Only do this using proc 0 - for why, see warning in Doxygen comment for the function. */
   if (data->rank==0) {
      /* This first test is to make sure it can handle files preceded with ./ */
      FindFile( fullPath, currDirFilename, searchPaths );
      pcu_check_streq( fullPath, currDirFilename );

      FindFile( fullPath, "currDirTest.xml", searchPaths );
      pcu_check_streq( fullPath, currDirFilename );
      
      FindFile( fullPath, "subDirTest.xml", searchPaths );
      pcu_check_streq( fullPath, subDirFilename );
      
      FindFile( fullPath, "nofile.man", searchPaths );
      pcu_check_streq( fullPath, "" );
      
      FindFile( fullPath, "/Users/luke/Projects/StGermain/env_vars", searchPaths );
      pcu_check_streq( fullPath, "" );
   }

   if (data->rank==0) {
      remove( currDirFilename );
      remove( subDirFilename );
      rmdir( subDir );
   }
}


void PathUtilsSuite( pcu_suite_t* suite ) {
   pcu_suite_setData( suite, PathUtilsSuiteData );
   pcu_suite_setFixtures( suite, PathUtilsSuite_Setup, PathUtilsSuite_Teardown );
   pcu_suite_addTest( suite, PathUtilsSuite_TestPathJoin );
   pcu_suite_addTest( suite, PathUtilsSuite_TestFindFile );
}


