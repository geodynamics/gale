/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street, Melbourne, 3053, Australia.
**
** Authors:
**	Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**	Patrick D. Sunter, Software Engineer, VPAC. (pds@vpac.org)
**	Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**	Siew-Ching Tan, Software Engineer, VPAC. (siew@vpac.org)
**	Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
**	Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
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
** $Id: testLibStGermainDynamic.c 2742 2005-03-05 05:33:43Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef __StGermain_Base_IO_PathUtils_h__
#define __StGermain_Base_IO_PathUtils_h__


/** Note that this function is designed as call-by-reference to modify fullPath to the found
 * full path of the chosen file.
 * WARNING: because this function uses fopen() to check the file's existence, then it
 * shouldn't be called by multiple processors simultaneously */
void FindFile( char* fullPath, Name filename, const char* searchPaths );

Bool FindFileInPathList( char* fullPath, char* filename, char** searchPaths, Index searchPathsSize );

void PathJoin( char* path, unsigned count, ... );

void PathClean( char* outPath, char* inPath );

char* ExpandEnvironmentVariables( char* string );

char* ParentDirectory( Name path );

Bool Stg_CreateDirectory( Name path );

Bool Stg_FileExists( Name path );

Bool Stg_DirectoryExists( Name path );

#endif /* __StGermain_Base_IO_PathUtils_h__ */
