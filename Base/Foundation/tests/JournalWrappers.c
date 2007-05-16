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
** $Id: JournalWrappers.c 4095 2007-05-16 00:53:05Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdarg.h>
#include "Base/Foundation/Foundation.h"
#include "JournalWrappers.h"
#include "stdio.h"
#include "stdlib.h"
#include "mpi.h"


const Type Info_Type = "info";
const Type Error_Type = "error";
const Type Debug_Type = "debug";

struct Stream* Journal_Register( Type type, Name name )
{
	return (void*)1;
}

int Journal_Printf ( void *stream, char *fmt, ... )
{
	va_list ap;

	va_start( ap, fmt );
	
	vfprintf( stdout, fmt, ap );

	va_end( ap );
	
	return 0;
}

int Journal_PrintfL ( void *stream, unsigned int level, char *fmt, ... )
{
	va_list ap;

	va_start( ap, fmt );
	
	vfprintf( stdout, fmt, ap );

	va_end( ap );
	
	return 0;
}

int Journal_Firewall ( int expression, void *stream, char *fmt, ... )
{
	va_list ap;

	va_start( ap, fmt );
	
	if( !expression ) {
		vfprintf( stdout, fmt, ap );
		assert( expression );
	}
	
	va_end( ap );
	
	return 0;
}

void Stream_Indent( void* stream )
{

}
void Stream_UnIndent( void* stream )
{

}
