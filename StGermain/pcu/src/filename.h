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

#ifndef pcu_files_h
#define pcu_files_h

extern pcu_suite_t* pcu_cursuite;

/** Useful for statically allocating character buffers for full path to expected and input files.
 * Note that since PCU paths are always relative, we can be confident that the directoreis will not be
 * unreasonably long. */
#define PCU_PATH_MAX 1024+FILENAME_MAX

#define pcu_filename_expectedLen( expectedFileName ) \
      _pcu_filename_expectedLen( (expectedFileName), pcu_cursuite->moduleDir )
unsigned _pcu_filename_expectedLen( const char* expectedFileName, const char* moduleDir );

#define pcu_filename_inputLen( inputFileName ) \
      _pcu_filename_inputLen( (inputFileName), pcu_cursuite->moduleDir )
unsigned _pcu_filename_inputLen( const char* inputFileName, const char* moduleDir );


/** Get the full path name of a given expected file for use in testing
 * Callers of this function should already have allocated the fullPathFileName buffer to the correct size using
 * pcu_filename_expectedLen() - or else just pass in a static buffer of size PCU_PATH_MAX */
#define pcu_filename_expected( expectedFileName, fullPathFileName ) \
      _pcu_filename_expected( (expectedFileName), (fullPathFileName), pcu_cursuite->moduleDir )
void _pcu_filename_expected( const char* const expectedFileName, char* const fullPathFileName,
      const char* moduleDir );

/** Get the full path name of a given input file for use in testing */
#define pcu_filename_input( inputFileName, fullPathFileName ) \
      _pcu_filename_input( (inputFileName), (fullPathFileName), pcu_cursuite->moduleDir )
void _pcu_filename_input( const char* const inputFileName, char* const fullPathFileName,
      const char* moduleDir );

#endif
