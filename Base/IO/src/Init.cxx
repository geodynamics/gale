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
** $Id: Init.c 4124 2007-05-27 23:18:25Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdarg.h>
#include <mpi.h>
#include <libxml/xmlerror.h>

#include "Base/Foundation/Foundation.h"

#include "types.h"
#include "shortcuts.h"
#include "Dictionary.h"
#include "Dictionary_Entry.h"
#include "Dictionary_Entry_Value.h"
#include "IO_Handler.h"
#include "XML_IO_Handler.h"
#include "Journal.h"
#include "JournalFile.h"
#include "CFile.h"
#include "MPIFile.h"
#include "Stream.h"
#include "CStream.h"
#include "MPIStream.h"
#include "BinaryStream.h"
#include "StreamFormatter.h"
#include "LineFormatter.h"
#include "IndentFormatter.h"
#include "RankFormatter.h"
#include "PathUtils.h"
#include "CmdLineArgs.h"
#include "Init.h"

const Name     LiveDebugName = "LiveDebug";
Stream*        LiveDebug = NULL;

Stream* stgErrorStream;

Bool BaseIO_Init( int* argc, char** argv[] )
{
	Stream*	general;

	/* Set up a useful map in the XML_IO_Handler */
	XML_IO_Handler_MergeTypeMap[Dictionary_MergeType_Append] = APPEND_TAG;
	XML_IO_Handler_MergeTypeMap[Dictionary_MergeType_Merge] = MERGE_TAG;
	XML_IO_Handler_MergeTypeMap[Dictionary_MergeType_Replace] = REPLACE_TAG;

	stgStreamFormatter_Buffer = StreamFormatter_Buffer_New(); 

	stJournal = Journal_New();
	/* Create default Typed Streams. */
	Journal_SetupDefaultTypedStreams();

	Journal_Printf( Journal_Register( DebugStream_Type, "Context" ), "In: %s\n", __func__ ); /* DO NOT CHANGE OR REMOVE */
	
	stgMemory->infoStream = Stg_Class_Copy( Journal_GetTypedStream( Info_Type ), 0, True, 0, 0 );
	stgMemory->debugStream = Stg_Class_Copy( Journal_GetTypedStream( Debug_Type ), 0, True, 0, 0 );
	stgMemory->errorStream = Stg_Class_Copy( Journal_GetTypedStream( Error_Type ), 0, True, 0, 0 );

	/* more inits from the Foundation level */
	Stream_Enable( Journal_Register( Info_Type, Stg_TimeMonitor_InfoStreamName ), False );
	Stream_Enable( Journal_Register( Info_Type, Stg_MemMonitor_InfoStreamName ), False );

	/* The LiveDebug stream */
	LiveDebug = Journal_Register( Info_Type, LiveDebugName );
	Stream_Enable( LiveDebug, False );
	Stream_SetLevel( LiveDebug, 1 );

	/* General streams. */
	general = Journal_Register( Info_Type, "general" );
	Stream_SetPrintingRank( general, 0 );
	stgErrorStream = Journal_Register( Error_Type, "stgErrorStream" );

	/* Handle the output of libXML properly, by redirecting to the XML_IO_Handler error stream */
	xmlSetGenericErrorFunc( NULL, XML_IO_Handler_LibXMLErrorHandler );
	
	return True;
}


