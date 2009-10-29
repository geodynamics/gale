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
** $Id: lucTestPlottingObjects.c 740 2007-10-11 08:05:31Z SteveQuenette $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>

#include <glucifer/Base/Base.h>


void lucTestCalculate( AbstractContext* context, lucWindow* window, lucViewport* viewport, lucPlottingObject* object ) {
	Stream* stream = Journal_Register( Info_Type, CURR_MODULE_NAME );

	Journal_Printf( stream, "In func '%s' for object %s\n", __func__, object->name );
}

void lucTestDraw( AbstractContext* context, lucWindow* window, lucViewport* viewport, lucPlottingObject* object ) {
	Stream* stream = Journal_Register( Info_Type, CURR_MODULE_NAME );

	Journal_Printf( stream, "In func '%s' for object %s\n", __func__, object->name );
}

void lucTestPlottingObjects( AbstractContext* context ) {
	lucBaseContextExtension* contextExt = ExtensionManager_Get( context->extensionMgr, context, lucBaseContextExtensionHandle );
	lucPlottingObject*            plottingObject;
	PlottingObject_Index          plottingObject_I;
	Stream*                       stream = Journal_Register( Info_Type, CURR_MODULE_NAME );

	Stream_RedirectFile_WithPrependedPath( stream, context->outputPath, "plottingObject.txt" );
	lucPlottingObject_Register_PrintAllObjects( contextExt->plottingObject_Register, stream );

	for ( plottingObject_I = 0 ; plottingObject_I < lucPlottingObject_Register_GetCount( contextExt->plottingObject_Register ) ; plottingObject_I++ ) {
		plottingObject = lucPlottingObject_Register_GetByIndex( contextExt->plottingObject_Register, plottingObject_I );

		EP_Append( plottingObject->CalculateEP, lucTestCalculate );
		EP_Append( plottingObject->DrawEP, lucTestDraw );
		
	}
	
}

void lucTestPlottingObjects_Draw( AbstractContext* context ) {
	lucBaseContextExtension* contextExt = ExtensionManager_Get( context->extensionMgr, context, lucBaseContextExtensionHandle );
	lucPlottingObject_Register_DrawAll( contextExt->plottingObject_Register, NULL, NULL, context );
}

void lucTestPlottingObjects_Register( void* context ) {
	ContextEP_Append( context, AbstractContext_EP_AssignFromXMLExtensions, lucTestPlottingObjects );
	ContextEP_Append( context, AbstractContext_EP_Dump, lucTestPlottingObjects_Draw );
}

