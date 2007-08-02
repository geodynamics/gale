/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
** Copyright (c) 2005, Monash Cluster Computing 
** All rights reserved.
** Redistribution and use in source and binary forms, with or without modification,
** are permitted provided that the following conditions are met:
**
** 		* Redistributions of source code must retain the above copyright notice, 
** 			this list of conditions and the following disclaimer.
** 		* Redistributions in binary form must reproduce the above copyright 
**			notice, this list of conditions and the following disclaimer in the 
**			documentation and/or other materials provided with the distribution.
** 		* Neither the name of the Monash University nor the names of its contributors 
**			may be used to endorse or promote products derived from this software 
**			without specific prior written permission.
**
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
** THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR 
** PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS 
** BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
** CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
** SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) 
** HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT 
** LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT 
** OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
**
**
** Contact:
*%		Cecile Duboz - Cecile.Duboz@sci.monash.edu.au
*%
** Contributors:
*+		Cecile Duboz
*+		Robert Turnbull
*+		Alan Lo
*+		Louis Moresi
*+		David Stegman
*+		David May
*+		Stevan Quenette
*+		Patrick Sunter
*+		Greg Watson
*+
** $Id: Init.c 728 2007-08-02 08:54:00Z SteveQuenette $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#include <mpi.h>
#include <StGermain/StGermain.h>

#include "glucifer.h"

#include <stdio.h>

Bool glucifer_Init() {
	int tmp;
	char* directory;

	lucBase_Init();
	lucWindowing_Init();
	lucRenderingEngines_Init();
	lucDrawingObjects_Init();
	lucOutputFormats_Init();
	lucInputFormats_Init();
	lucWindowInteractions_Init();
	
	Journal_Printf( Journal_Register( DebugStream_Type, "Context" ), "In: %s\n", __func__ ); /* DO NOT CHANGE OR REMOVE */
	tmp = Stream_GetPrintingRank( Journal_Register( InfoStream_Type, "Context" ) );
	Stream_SetPrintingRank( Journal_Register( InfoStream_Type, "Context" ), 0 );
	Journal_Printf( /* DO NOT CHANGE OR REMOVE */
		Journal_Register( InfoStream_Type, "Context" ), 
		"glucifer (Visualisation framework) revision %s. Copyright (C) 2005 Monash Cluster Computing.\n", VERSION );
	Stream_Flush( Journal_Register( InfoStream_Type, "Context" ) );
	Stream_SetPrintingRank( Journal_Register( InfoStream_Type, "Context" ), tmp );

	/* Add the gLucifer path to the global xml path dictionary */
	directory = Memory_Alloc_Array( char, 200, "xmlDirectory" ) ;
	sprintf(directory, "%s%s", LIB_DIR, "/StGermain" );
	XML_IO_Handler_AddDirectory( "gLucifer", directory );
	Memory_Free(directory);

	/* Add the plugin path to the global plugin list */
	ModulesManager_AddDirectory( "gLucifer", LIB_DIR );

	return True;
}
