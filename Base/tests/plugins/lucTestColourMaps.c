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
** $Id: lucTestColourMaps.c 740 2007-10-11 08:05:31Z SteveQuenette $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>

#include <glucifer/Base/Base.h>

void lucTestColourMaps_Function( AbstractContext* context ) {
	Index            component_I;
	Stg_ObjectList*  componentList = context->CF->LCRegister->componentList;
	Stg_Component*   component;
	Stream*          stream = Journal_Register( Info_Type, CURR_MODULE_NAME );

	Stream_RedirectFile_WithPrependedPath( stream, context->outputPath, "colourMap.txt" );

	/* Do Stg_Class_Print for all colour maps */
	for ( component_I = 0 ; component_I < componentList->count ; component_I++ ) {
		component = (Stg_Component*) Stg_ObjectList_At( componentList, component_I );

		if ( Stg_Class_IsInstance( component, lucColourMap_Type ) )
			Stg_Class_Print( component, stream );
	}
}

const Type lucTestColourMaps_Type = "lucTestColourMaps";
typedef struct {
	__Codelet
} lucTestColourMaps;

void _lucTestColourMaps_AssignFromXML( void* component, Stg_ComponentFactory* cf, void* data ) {
	AbstractContext* context;
	context = Stg_ComponentFactory_ConstructByName( cf, "context", AbstractContext, True, data ); 
	ContextEP_Append( context, AbstractContext_EP_AssignFromXMLExtensions, lucTestColourMaps_Function );
}

void* _lucTestColourMaps_DefaultNew( Name name ) {
	return Codelet_New(
		lucTestColourMaps_Type,
		_lucTestColourMaps_DefaultNew,
		_lucTestColourMaps_AssignFromXML,
		_Codelet_Build,
		_Codelet_Initialise,
		_Codelet_Execute,
		_Codelet_Destroy,
		name );
}

Index lucTestColourMaps_Register( PluginsManager* pluginsManager ) {
	return PluginsManager_Submit( pluginsManager, lucTestColourMaps_Type, "0", _lucTestColourMaps_DefaultNew );
}



