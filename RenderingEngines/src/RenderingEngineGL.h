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
** $Id: RenderingEngineGL.h 628 2006-10-12 08:23:07Z SteveQuenette $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifdef HAVE_OPENGL_FRAMEWORK
	#include <OpenGL/gl.h>
	#include <OpenGL/glu.h>
#else
	#include <gl.h>
	#include <glu.h>
#endif

#ifndef __lucRenderingEngineGL_h__
#define __lucRenderingEngineGL_h__

	/** Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
	extern const Type lucRenderingEngineGL_Type;
		
	/** Class contents - this is defined as a macro so that sub-classes of this class can use this macro at the start of the definition of their struct */
	#define __lucRenderingEngineGL \
		/* Macro defining parent goes here - This means you can cast this class as its parent */ \
		__lucRenderingEngine \
		/* Virtual functions go here */ \
		/* Other info */\
		GLboolean               doubleBuffered;            \

	struct lucRenderingEngineGL { __lucRenderingEngineGL };
	
	/** Private Constructor: This will accept all the virtual functions for this class as arguments. */
	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define LUCRENDERINGENGINEGL_DEFARGS \
                LUCRENDERINGENGINE_DEFARGS

	#define LUCRENDERINGENGINEGL_PASSARGS \
                LUCRENDERINGENGINE_PASSARGS

	lucRenderingEngineGL* _lucRenderingEngineGL_New(  LUCRENDERINGENGINEGL_DEFARGS  );

	void _lucRenderingEngineGL_Delete( void* renderingEngine ) ;
	void _lucRenderingEngineGL_Print( void* renderingEngine, Stream* stream ) ;
	void* _lucRenderingEngineGL_Copy( void* renderingEngine, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap) ;

	/* 'Stg_Component' implementations */
	void* _lucRenderingEngineGL_DefaultNew( Name name ) ;
	void _lucRenderingEngineGL_AssignFromXML( void* renderingEngine, Stg_ComponentFactory* cf, void* data );
	void _lucRenderingEngineGL_Build( void* renderingEngine, void* data ) ;
	void _lucRenderingEngineGL_Initialise( void* renderingEngine, void* data ) ;
	void _lucRenderingEngineGL_Execute( void* renderingEngine, void* data );
	void _lucRenderingEngineGL_Destroy( void* renderingEngine, void* data ) ;

	void _lucRenderingEngineGL_Render( void* renderingEngine, lucWindow* window, AbstractContext* context ) ;
	void _lucRenderingEngineGL_Clear( void* renderingEngine, lucWindow* window, Bool clearAll ) ;
	void _lucRenderingEngineGL_GetPixelData( void* renderingEngine, lucWindow* window, void* buffer, Bool withAlpha ) ;
	void lucRenderingEngineGL_WriteViewportText( void* renderingEngine, lucWindow* window, lucViewportInfo* viewportInfo, AbstractContext* context ) ;

	/** Compositing Functions */
	Index lucRenderingEngineGL_MapBufferIdToRank( void* renderingEngine, Index bufferId, Index mergeCount ) ;
	void lucRenderingEngineGL_CombineToMaster( 
		void*                                              renderingEngine,
		lucViewportInfo*                                   viewportInfo,
		AbstractContext*                                   context,
		lucPixel*                                          imageBuffer, 
		float*                                             depthBuffer );

	void _lucRenderingEngineGL_CompositeViewport_Stencil( 
		void*                                              renderingEngine, 
		lucViewportInfo*                                   viewportInfo, 
		AbstractContext*                                   context, 
		Bool                                               broadcast ) ;

	void _lucRenderingEngineGL_CompositeViewport_Manual( 
		void*                                              renderingEngine, 
		lucViewportInfo*                                   viewportInfo, 
		AbstractContext*                                   context, 
		Bool                                               broadcast ) ;

#endif

