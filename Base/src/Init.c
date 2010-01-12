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
** $Id: Init.c 740 2007-10-11 08:05:31Z SteveQuenette $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>

#include "Base.h"

Stream* lucInfo  = NULL;
Stream* lucDebug = NULL;
Stream* lucError = NULL;
lucInputFormat_Register* lucInputFormat_Register_Singleton = NULL;

Bool lucBase_Init() {
	Stg_ComponentRegister* componentRegister = Stg_ComponentRegister_Get_ComponentRegister();

	Journal_Printf( Journal_Register( DebugStream_Type, "Context" ), "In: %s\n", __func__ ); /* DO NOT CHANGE OR REMOVE */

	/* Set up streams */
	lucInfo  = Journal_Register( InfoStream_Type, "lucInfo" );
	lucDebug = Journal_Register( DebugStream_Type, "lucDebug" );
	lucError = Journal_Register( ErrorStream_Type, "lucError" );	
	lucInputFormat_Register_Singleton = lucInputFormat_Register_New();
	
	Stg_ComponentRegister_Add( componentRegister, lucCamera_Type,     "0", _lucCamera_DefaultNew );
	Stg_ComponentRegister_Add( componentRegister, lucColourMap_Type,  "0", _lucColourMap_DefaultNew );
	Stg_ComponentRegister_Add( componentRegister, lucViewport_Type,   "0", _lucViewport_DefaultNew );
	Stg_ComponentRegister_Add( componentRegister, lucWindow_Type,     "0", _lucWindow_DefaultNew );
	Stg_ComponentRegister_Add( componentRegister, lucLight_Type,     "0", _lucLight_DefaultNew );

	/* Register Parents for type checking */
	RegisterParent( lucCamera_Type,            Stg_Component_Type );
	RegisterParent( lucColourMap_Type,         Stg_Component_Type );
	RegisterParent( lucDrawingObject_Type,     Stg_Component_Type );
	RegisterParent( lucViewport_Type,          Stg_Component_Type );
	RegisterParent( lucOutputFormat_Type,      Stg_Component_Type );
	RegisterParent( lucInputFormat_Type,       Stg_Component_Type );
	RegisterParent( lucWindow_Type,            Stg_Component_Type );
	RegisterParent( lucRenderingEngine_Type,   Stg_Component_Type );
	RegisterParent( lucWindowInteraction_Type, Stg_Component_Type );
	RegisterParent( lucLight_Type,             Stg_Component_Type );

	
	
	RegisterParent( lucDrawingObject_Register_Type, NamedObject_Register_Type );
	RegisterParent( lucOutputFormat_Register_Type,  NamedObject_Register_Type );
	RegisterParent( lucWindowInteraction_Register_Type,  NamedObject_Register_Type );
	RegisterParent( lucLight_Register_Type,  NamedObject_Register_Type );

	/* Create MPI Datatypes */
	lucCamera_Create_MPI_Datatype();
	lucViewport_Create_MPI_Datatype();
	lucViewportInfo_Create_MPI_Datatype();
	lucWindow_Create_MPI_Datatype();

	return True;
}



