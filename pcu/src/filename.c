/*
**  Copyright 2008 Luke Hodkinson
**
**  This file is part of pcu.
**
**  pcu is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License as published by
**  the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  Foobar is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
**  GNU General Public License for more details.
**
**  You should have received a copy of the GNU General Public License
**  along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <mpi.h>
#include "types.h"

char* PCU_FILENAME_PROJECT = NULL;
char* PCU_FILENAME_MODULE = NULL;

const char* pcu_filename_getProject();
const char* pcu_filename_getModule();

void pcu_filename_setProject( const char* const projectName ) {
   if ( PCU_FILENAME_PROJECT ) {
      free( PCU_FILENAME_PROJECT );
   }

   PCU_FILENAME_PROJECT = strdup( projectName );
}

const char* pcu_filename_getProject() {
   assert(PCU_FILENAME_PROJECT);
   return PCU_FILENAME_PROJECT;
}


void pcu_filename_setModule( const char* const moduleName ) {
   if ( PCU_FILENAME_MODULE ) {
      free( PCU_FILENAME_MODULE );
   }

   PCU_FILENAME_MODULE = strdup( moduleName );
}

const char* pcu_filename_getModule() {
   assert(PCU_FILENAME_MODULE);
   return PCU_FILENAME_MODULE;
}


unsigned pcu_filename_expectedLen( const char* expectedFileName ) {
   const char*          fileType = "expected";
   const char*          projSubDir = pcu_filename_getProject();
   const char*          testSuiteSubDir = pcu_filename_getModule(); 

   assert( expectedFileName );

   return 1+1+strlen(projSubDir)+1+strlen(testSuiteSubDir)+1
      +strlen(fileType)+1+strlen(expectedFileName) + 1;
}


/* Callers of this function should already have allocated the fullPathFileName buffer to the correct size using
 * pcu_filename_expectedLen */
void pcu_filename_expected( const char* const expectedFileName, char* const fullPathFileName ) {
   const char*          fileType = "expected";
   const char*          projSubDir = pcu_filename_getProject();
   const char*          testSuiteSubDir = pcu_filename_getModule(); 

   assert( expectedFileName );
   assert( fullPathFileName );
   sprintf( fullPathFileName, "./%s/%s/%s/%s", projSubDir, testSuiteSubDir, fileType, expectedFileName );
}


unsigned pcu_filename_inputLen( const char* inputFileName ) {
   const char*          fileType = "input";
   const char*          projSubDir = pcu_filename_getProject();
   const char*          testSuiteSubDir = pcu_filename_getModule(); 

   assert( inputFileName );

   return 1+1+strlen(projSubDir)+1+strlen(testSuiteSubDir)+1
      +strlen(fileType)+1+strlen(inputFileName) + 1;
}


void pcu_filename_input( const char* const inputFileName, char* const fullPathFileName ) {
   const char*          fileType = "input";
   const char*          projSubDir = pcu_filename_getProject();
   const char*          testSuiteSubDir = pcu_filename_getModule(); 

   sprintf( fullPathFileName, "./%s/%s/%s/%s", projSubDir, testSuiteSubDir, fileType, inputFileName );
}
