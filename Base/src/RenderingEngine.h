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
** $Id: RenderingEngine.h 628 2006-10-12 08:23:07Z SteveQuenette $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#ifndef __lucRenderingEngine_h__
#define __lucRenderingEngine_h__

	extern const Type lucRenderingEngine_Type;

	typedef void (lucRenderingEngine_RenderFunction) ( void* renderingEngine, lucWindow* window, AbstractContext* context);
	typedef void (lucRenderingEngine_ClearFunction) ( void* renderingEngine, lucWindow* window, Bool clearAll );
	typedef void (lucRenderingEngine_GetPixelDataFunction) ( void* renderingEngine, lucWindow* window, lucPixel* pixelData);
	typedef void (lucRenderingEngine_CompositeViewportFunction) (  
		void*                                              renderingEngine, 
		lucViewportInfo*                                   viewportInfo, 
		AbstractContext*                                   context, 
		Bool                                               broadcast );


	#define __lucRenderingEngine                           \
		__Stg_Component                                   \
		AbstractContext*				   context;		    \
		/* Virtual Functions */ \
		lucRenderingEngine_RenderFunction*                 _render;                 \
		lucRenderingEngine_ClearFunction*				   _clear;					\
		lucRenderingEngine_GetPixelDataFunction*           _getPixelData;           \
		lucRenderingEngine_CompositeViewportFunction*      _compositeViewport;	    \

	struct lucRenderingEngine {__lucRenderingEngine};

	lucRenderingEngine* _lucRenderingEngine_New(
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
		lucRenderingEngine_ClearFunction*				   _clear,
		lucRenderingEngine_GetPixelDataFunction*           _getPixelData,
		lucRenderingEngine_CompositeViewportFunction*      _compositeViewport,
		Name                                               name );


	void _lucRenderingEngine_Delete( void* renderingEngine ) ;
	void _lucRenderingEngine_Print( void* renderingEngine, Stream* stream ) ;
	void* _lucRenderingEngine_Copy( void* renderingEngine, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) ;

	void _lucRenderingEngine_Construct( void* renderingEngine, Stg_ComponentFactory* cf, void* data ) ;
	void _lucRenderingEngine_Build( void* renderingEngine, void* data );
	void _lucRenderingEngine_Initialise( void* renderingEngine, void* data );
	void _lucRenderingEngine_Execute( void* renderingEngine, void* data );
	void _lucRenderingEngine_Destroy( void* renderingEngine, void* data );

	/* +++ Public Functions +++ */
	void lucRenderingEngine_Render( void* renderingEngine, lucWindow* window, AbstractContext* context ) ;
	void lucRenderingEngine_Clear( void* renderingEngine, lucWindow* window, Bool clearAll ) ;
	void lucRenderingEngine_GetPixelData( void* renderingEngine, lucWindow* window, lucPixel* pixelData ) ;
	void lucRenderingEngine_CompositeViewport( 
		void*                                              renderingEngine, 
		lucViewportInfo*                                   viewportInfo, 
		AbstractContext*                                   context, 
		Bool                                               broadcast );

#endif
