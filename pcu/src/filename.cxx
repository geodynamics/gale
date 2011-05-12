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


unsigned _pcu_filename_expectedLen( const char* expectedFileName, const char* moduleDir ) {
   const char*          fileType = "expected";

   assert( expectedFileName );

   return 1+1+strlen(moduleDir)+1+strlen(fileType)+1+strlen(expectedFileName) + 1;
}


/* Callers of this function should already have allocated the fullPathFileName buffer to the correct size using
 * pcu_filename_expectedLen */
void _pcu_filename_expected( const char* const expectedFileName, char* const fullPathFileName,
      const char* moduleDir )
{
   const char*          fileType = "expected";

   assert( expectedFileName );
   assert( fullPathFileName );
   sprintf( fullPathFileName, "./%s/%s/%s", moduleDir, fileType, expectedFileName );
}


unsigned _pcu_filename_inputLen( const char* inputFileName, const char* moduleDir ) {
   const char*          fileType = "input";

   assert( inputFileName );

   return 1+1+strlen(moduleDir)+1+strlen(fileType)+1+strlen(inputFileName) + 1;
}


void _pcu_filename_input( const char* const inputFileName, char* const fullPathFileName,
      const char* moduleDir )
{
   const char*          fileType = "input";

   sprintf( fullPathFileName, "./%s/%s/%s", moduleDir, fileType, inputFileName );
}


