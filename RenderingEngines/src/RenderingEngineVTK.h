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
** $Id: Arrhenius.c 78 2005-11-29 11:58:21Z RobertTurnbull $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#ifndef __lucRenderingEngineVTK_h__
#define __lucRenderingEngineVTK_h__



	/** Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
	extern const Type lucRenderingEngineVTK_Type;
		
	/** Class contents - this is defined as a macro so that sub-classes of this class can use this macro at the start of the definition of their struct */
	#define __lucRenderingEngineVTK \
		/* Macro defining parent goes here - This means you can cast this class as its parent */ \
		__lucRenderingEngine; \
		/* Virtual functions go here */ \
		/* Other info */\
    int nViewports; \
		lucColour* bgColour; 
		
	struct lucRenderingEngineVTK { __lucRenderingEngineVTK };
	
	/** Private Constructor: This will accept all the virtual functions for this class as arguments. */
	lucRenderingEngineVTK* _lucRenderingEngineVTK_New( 
		SizeT                                              sizeOfSelf,
		Type                                               type,
		Stg_Class_DeleteFunction*                          _delete,
		Stg_Class_PrintFunction*                           _print,
		Stg_Class_CopyFunction*                            _copy, 
		Stg_Component_DefaultConstructorFunction*          _defaultConstructor,
		Stg_Component_ConstructFunction*                   _construct,
		Stg_Component_BuildFunction*                       _build,
		Stg_Component_InitialiseFunction*                  _initialise,
		Stg_Component_ExecuteFunction*                     _execute,
		Stg_Component_DestroyFunction*                     _destroy,
		lucRenderingEngine_RenderFunction*                 _render,
		lucRenderingEngine_ClearFunction*             	   _clear,
		lucRenderingEngine_GetPixelDataFunction*           _getPixelData,
		lucRenderingEngine_CompositeViewportFunction*      _compositeViewport,
		Name                                               name );

	void _lucRenderingEngineVTK_Delete( void* renderingEngine ) ;
	void _lucRenderingEngineVTK_Print( void* renderingEngine, Stream* stream ) ;
	void* _lucRenderingEngineVTK_Copy( void* renderingEngine, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap) ;

	/* 'Stg_Component' implementations */
	void* _lucRenderingEngineVTK_DefaultNew( Name name ) ;
	void _lucRenderingEngineVTK_AssignFromXML( void* renderingEngine, Stg_ComponentFactory* cf, void* data );
	void _lucRenderingEngineVTK_Build( void* renderingEngine, void* data ) ;
	void _lucRenderingEngineVTK_Initialise( void* renderingEngine, void* data ) ;
	void _lucRenderingEngineVTK_Execute( void* renderingEngine, void* data );
	void _lucRenderingEngineVTK_Destroy( void* renderingEngine, void* data ) ;

	void _lucRenderingEngineVTK_Render( void* renderingEngine, lucWindow* window, AbstractContext* context ) ;
	void _lucRenderingEngineVTK_Clear( void* renderingEngine, lucWindow* window, Bool clearAll ) ;
	void _lucRenderingEngineVTK_GetPixelData( void* renderingEngine, lucWindow* window, lucPixel* buffer ) ;

	void lucRenderingEngineVTK_DrawTitle( void* renderingEngine, lucWindow* window, lucViewportInfo* viewportInfo ) ;

	/** Compositing Functions */
	Index lucRenderingEngineVTK_MapBufferIdToRank( void* renderingEngine, Index bufferId, Index mergeCount ) ;
	void lucRenderingEngineVTK_CombineToMaster( 
		void*                                              renderingEngine,
		lucViewportInfo*                                   viewportInfo,
		AbstractContext*                                   context,
		lucPixel*                                          imageBuffer, 
		float*                                             depthBuffer );

	void _lucRenderingEngineVTK_CompositeViewport_Stencil( 
		void*                                              renderingEngine, 
		lucViewportInfo*                                   viewportInfo, 
		AbstractContext*                                   context, 
		Bool                                               broadcast ) ;

	void lucRenderingEngineVTK_CombineToMaster( 
		void*                                              renderingEngine,
		lucViewportInfo*                                   viewportInfo,
		AbstractContext*                                   context,
		lucPixel*                                          imageBuffer, 
		float*                                             depthBuffer );

	void _lucRenderingEngineVTK_CompositeViewport_Manual( 
		void*                                              renderingEngine, 
		lucViewportInfo*                                   viewportInfo, 
		AbstractContext*                                   context, 
		Bool                                               broadcast ) ;

#endif
