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
** $Id: RenderingEngineGL.c 791 2008-09-01 02:09:06Z JulianGiordani $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>

#include <glucifer/Base/Base.h>

#include "types.h"
#include "RenderingEngineGL.h"
#include "OpenGlUtil.h"

#include <assert.h>
#ifdef HAVE_OPENGL_FRAMEWORK
	#include <OpenGL/gl.h>
	#include <OpenGL/glu.h>
#else
	#include <gl.h>
	#include <glu.h>
#endif
#include <string.h>

#ifdef HAVE_GL2PS
	#include <gl2ps.h>
#endif

#ifndef MASTER
	#define MASTER 0
#endif

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type lucRenderingEngineGL_Type = "lucRenderingEngineGL";

/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
lucRenderingEngineGL* _lucRenderingEngineGL_New( 
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
		lucRenderingEngine_GetPixelDataFunction*           _getPixelData,
		lucRenderingEngine_CompositeViewportFunction*      _compositeViewport,
		Name                                               name ) 
{
	lucRenderingEngineGL*					self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( sizeOfSelf >= sizeof(lucRenderingEngineGL) );
	self = (lucRenderingEngineGL*) _lucRenderingEngine_New( 
			sizeOfSelf,
			type, 
			_delete,
			_print,
			_copy,
			_defaultConstructor,
			_construct,
			_build,
			_initialise,
			_execute,
			_destroy,
			_render,
			_getPixelData,
			_compositeViewport,
			name );
	
	return self;
}

void _lucRenderingEngineGL_Init( 
		lucRenderingEngineGL*                                      self  ) 
{
}

void _lucRenderingEngineGL_Delete( void* renderingEngine ) {
	lucRenderingEngineGL*  self = (lucRenderingEngineGL*)renderingEngine;

	_lucRenderingEngine_Delete( self );
}

void _lucRenderingEngineGL_Print( void* renderingEngine, Stream* stream ) {
	lucRenderingEngineGL*  self = (lucRenderingEngineGL*)renderingEngine;

	_lucRenderingEngine_Print( self, stream );
}

void* _lucRenderingEngineGL_Copy( void* renderingEngine, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap) {
	lucRenderingEngineGL*  self = (lucRenderingEngineGL*)renderingEngine;
	lucRenderingEngineGL* newRenderingEngine;

	newRenderingEngine = _lucRenderingEngine_Copy( self, dest, deep, nameExt, ptrMap );

	/* TODO */
	abort();

	return (void*) newRenderingEngine;
}


void* _lucRenderingEngineGL_DefaultNew( Name name ) {
	return (void*) _lucRenderingEngineGL_New(
		sizeof(lucRenderingEngineGL),
		lucRenderingEngineGL_Type,
		_lucRenderingEngineGL_Delete,
		_lucRenderingEngineGL_Print,
		NULL,
		_lucRenderingEngineGL_DefaultNew,
		_lucRenderingEngineGL_Construct,
		_lucRenderingEngineGL_Build,
		_lucRenderingEngineGL_Initialise,
		_lucRenderingEngineGL_Execute,
		_lucRenderingEngineGL_Destroy,
		_lucRenderingEngineGL_Render,
		_lucRenderingEngineGL_GetPixelData,
		_lucRenderingEngineGL_CompositeViewport_Manual,
		name );
}

void _lucRenderingEngineGL_Construct( void* renderingEngine, Stg_ComponentFactory* cf, void* data ){
	lucRenderingEngineGL*  self = (lucRenderingEngineGL*)renderingEngine;

	/* Construct Parent */
	_lucRenderingEngine_Construct( self, cf, data );
	
	_lucRenderingEngineGL_Init( self );
}

void _lucRenderingEngineGL_Build( void* renderingEngine, void* data ) {}
void _lucRenderingEngineGL_Initialise( void* renderingEngine, void* data ) {}
void _lucRenderingEngineGL_Execute( void* renderingEngine, void* data ) {}
void _lucRenderingEngineGL_Destroy( void* renderingEngine, void* data ) {}


void _lucRenderingEngineGL_Render( void* renderingEngine, lucWindow* window, AbstractContext* context ) {
	lucRenderingEngineGL* self              = (lucRenderingEngineGL*) renderingEngine;
	lucViewport*          viewport;
	Viewport_Index        viewport_I;
	Viewport_Index        viewportCount     = window->viewportCount;
	lucViewportInfo*      viewportInfo;
	GLint                 viewport_gl2ps[4], state;
	
	Journal_DPrintfL( lucDebug, 2, "In func: %s for %s '%s'\n", __func__, self->type, self->name );
	Stream_Indent( lucDebug );

	/* Set up OpenGl Colour */
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);	
	glEnable(GL_COLOR_MATERIAL);
	glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
 	/*	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);  -- Done in the viewport now	*/	
	glDrawBuffer(GL_BACK_LEFT);

	/* Allow Transparency */
	glEnable (GL_BLEND);
	glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
	glColorMask(GL_TRUE,GL_TRUE,GL_TRUE,GL_TRUE);
	
	glEnable(GL_DEPTH_TEST);

	/* Set up lighting */

	lucWindow_Broadcast( window, 0, MPI_COMM_WORLD );
	lucWindow_CheckCameraFlag( window );
	lucWindow_CheckLightFlag( window );
		
	for ( viewport_I = 0 ; viewport_I < viewportCount ; viewport_I++ ) {
		viewportInfo = &window->viewportInfoList[ viewport_I ];
		viewport = viewportInfo->viewport;

		Journal_DPrintfL( lucDebug, 2, "In loop for viewport '%s'.\n", viewport->name );
		Stream_Indent( lucDebug );
		
		/* Set viewport */
		glViewport( viewportInfo->startx, viewportInfo->starty, viewportInfo->width, viewportInfo->height);
		glScissor( viewportInfo->startx, viewportInfo->starty, viewportInfo->width, viewportInfo->height);
		if ( ! viewportInfo->needsToDraw ) {
			Journal_DPrintfL( lucDebug, 2, "Viewport '%s' doesn't need to be redrawn.\n", viewport->name );
			Stream_UnIndent( lucDebug );
			continue;
		}
		
		#ifdef HAVE_GL2PS
			/* calls to gl2ps so that different viewports are created for vector image outputs */
			glGetIntegerv(GL_VIEWPORT, viewport_gl2ps);
			state = gl2psBeginViewport(viewport_gl2ps);
			if(state == 5)
				Journal_Printf( Journal_MyStream( Error_Type, self ), "\nError. Insufficient GL feedback buffer size. \ 
										       \nConsider increasing the OutputVECTOR buffersize. \
										      \nVector image will not be created correctly.\n\n" );
		#endif
		
		lucRenderingEngineGL_Clear( self, window );

		if (context->rank == MASTER)
			lucRenderingEngineGL_WriteViewportText( self, window, viewportInfo, context );
			
			
		switch ( viewport->camera->stereoType ) {
			case lucMono:
				lucViewportInfo_SetOpenGLCamera( viewportInfo );
				lucViewport_Draw( viewport, window, viewportInfo, context );
				break;
			case lucStereoToeIn: case lucStereoAsymmetric:
				glDrawBuffer(GL_BACK_RIGHT);
				viewport->camera->buffer = lucRight;

				lucViewportInfo_SetOpenGLCamera( viewportInfo );
				lucViewport_Draw( viewport, window, viewportInfo, context );

				glDrawBuffer(GL_BACK_LEFT);
				viewport->camera->buffer = lucLeft;
				
				lucViewportInfo_SetOpenGLCamera( viewportInfo );
				lucViewport_Draw( viewport, window, viewportInfo, context );
		}
		
		viewportInfo->needsToDraw = False;

		Stream_UnIndent( lucDebug );
		Journal_DPrintfL( lucDebug, 2, "Finised loop.\n" );

		#ifdef HAVE_GL2PS		
			state = gl2psEndViewport();
			if(state == 5)
				Journal_Printf( Journal_MyStream( Error_Type, self ), "\nError. Insufficient GL feedback buffer size. \ 
										       \nConsider increasing the OutputVECTOR buffersize. \
										      \nVector image will not be created correctly.\n\n" );
		#endif		
	}
	
	Stream_UnIndent( lucDebug );
	Journal_DPrintfL( lucDebug, 2, "Leaving func %s\n", __func__ );
}

void _lucRenderingEngineGL_GetPixelData( void* renderingEngine, lucWindow* window, lucPixel* buffer ) {
	GLsizei width  = window->width;
	GLsizei height = window->height;
	
	glPixelStorei(GL_PACK_ALIGNMENT,1);
	
	if ( lucWindow_HasStereoCamera( window ) ) {
		if ( window->currStereoBuffer == lucRight )
			glReadBuffer( GL_FRONT_RIGHT );
		else
			glReadBuffer( GL_FRONT_LEFT );
	}

	/* Actually read the pixels. */
	glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, buffer); 
}

void lucRenderingEngineGL_WriteViewportText( void* renderingEngine, lucWindow* window, lucViewportInfo* viewportInfo, AbstractContext* context ) {
	lucViewport* viewport = viewportInfo->viewport;
//	int stringWidth = lucStringWidth( viewport->name );

	/* Set up 2D Viewer the size of the viewport */
	glPushMatrix();
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D((GLfloat) 0.0, (GLfloat) viewportInfo->width, (GLfloat) 0.0, (GLfloat) viewportInfo->height );
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();	

	/* Set the colour so that it'll show up against the background */
	lucColour_SetComplimentaryOpenGLColour( &window->backgroundColour );


	/* Print Time Stamp */
	if (viewport->drawTime) {
		char* timeString;
		Stg_asprintf( &timeString, "%g", context->currentTime );
		glRasterPos2i( 0, viewportInfo->height - 13 );
		lucPrintString( timeString );
		Memory_Free( timeString );
	}

	glPopMatrix();
}

void lucRenderingEngineGL_Clear( void* renderingEngineGL, lucWindow* window ) {
	glEnable (GL_SCISSOR_TEST);
	glClearColor( 
		window->backgroundColour.red, 
		window->backgroundColour.green, 
		window->backgroundColour.blue, 
		window->backgroundColour.opacity );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
}

Index lucRenderingEngineGL_MapBufferIdToRank( void* renderingEngineGL, Index bufferId, Index mergeCount ) {
	Index merge_I;
	Index rank;

	rank = bufferId;
	for ( merge_I = 0 ; merge_I < mergeCount ; merge_I++ )
		rank *= 2;

	return rank;
}

/* This function is quite a bit faster that the _lucRenderingEngineGL_CompositeViewport_Manual one but it will not work if
 * you don't have the stencil buffer - 
 * There's also been some problems with this one on the altix in uq - This shouldn't be the default function yet until this
 * problem is sorted out */
void _lucRenderingEngineGL_CompositeViewport_Stencil( 
		void*                                              renderingEngine, 
		lucViewportInfo*                                   viewportInfo, 
		AbstractContext*                                   context, 
		Bool                                               broadcast ) 
{
	lucRenderingEngineGL* self                 = (lucRenderingEngineGL*) renderingEngine;
	Pixel_Index           width                = viewportInfo->width;
	Pixel_Index           height               = viewportInfo->height;
	Pixel_Index           startx               = viewportInfo->startx;
	Pixel_Index           starty               = viewportInfo->starty;
	Index                 buffersToMerge;
	Index                 neighbourRank;
	Index                 bufferId             = context->rank;
	Pixel_Index           pixelCount           = width*height;
	Index                 mergeCount           = 0;
	float*                depthBuffer;
	lucPixel*             imageBuffer;
	MPI_Status            status;
	MPI_Comm              comm                 = context->communicator;
	
	glEnable(GL_STENCIL_TEST);				
	/* Make sure that we can use the stencil buffer */
	if ( !glIsEnabled( GL_STENCIL_TEST ) ) {
		self->_compositeViewport = _lucRenderingEngineGL_CompositeViewport_Manual;
		_lucRenderingEngineGL_CompositeViewport_Manual( self, viewportInfo, context, broadcast );
		return;
	}

	/* Set matrices up so we are in pixel coordinates */
	glPushMatrix();
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D( 0.0, (GLfloat) width, 0.0, (GLfloat) height );
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();	
	glRasterPos2i( 0, 0 );
	glPixelStorei(GL_PACK_ALIGNMENT,1);
	glPixelStorei(GL_UNPACK_ALIGNMENT,1);
	glClearStencil(0x0);

	imageBuffer = Memory_Alloc_Array( lucPixel, pixelCount, "Image Buffer" );
	depthBuffer = Memory_Alloc_Array( float,    pixelCount, "Depth Buffer" );

	/* Merge Buffers */
	for ( buffersToMerge = context->nproc - 1 ; buffersToMerge > 0 ; buffersToMerge = div( buffersToMerge, 2 ).quot ) {
		/* If my ID is odd - then send the buffer and get out of loop */
		if ( bufferId % 2 == 1 ) {
			/* Send buffer to left */
			neighbourRank = lucRenderingEngineGL_MapBufferIdToRank( self, bufferId - 1, mergeCount );
			Journal_DPrintfL( lucDebug, 2, "Sending buffers to processor '%d'\n", neighbourRank );

			glReadPixels( startx, starty, width, height, GL_DEPTH_COMPONENT, GL_FLOAT, depthBuffer);
			MPI_Send( depthBuffer, pixelCount, MPI_FLOAT, neighbourRank, GL_DEPTH_COMPONENT, comm );

			glReadPixels( startx, starty, width, height, GL_RGB, GL_UNSIGNED_BYTE, imageBuffer);
			MPI_Send( imageBuffer, pixelCount * 3, MPI_UNSIGNED_CHAR, neighbourRank, GL_RGB,             comm );
			
			/* Now that I've sent my info - I can quit this loop */
			break;
		}
		else {
			/* Only merge if you are not the last processor to have a buffer */
			if ( bufferId < buffersToMerge ) {
				/* Receive Buffer from Right */
				neighbourRank = lucRenderingEngineGL_MapBufferIdToRank( self, bufferId + 1, mergeCount );
				Journal_DPrintfL( lucDebug, 2, "Receiving buffers from processor '%d'\n", neighbourRank );

				/**************** Merge these two buffers ***********************/
				/* See http://www.opengl.org/resources/tutorials/advanced/advanced97/notes/node200.html */

				/* Clear the stencil buffer */
				glClear( GL_STENCIL_BUFFER_BIT );

				/* Ensure depth testing is set */
				glEnable(GL_DEPTH_TEST);

				/* Disable the color buffer for writing */
				glColorMask( GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE );

				/* Set stencil values to 1 when the depth test passes */
				glStencilFunc(GL_ALWAYS, 1, 1);
				glStencilOp(GL_KEEP, GL_KEEP, GL_REPLACE);

				/* Draw the depth values to the frame buffer */
				glRasterPos2i( 0, 0 );
				MPI_Recv( depthBuffer, pixelCount, MPI_FLOAT, neighbourRank, GL_DEPTH_COMPONENT, comm, &status );
				glDrawPixels( width, height, GL_DEPTH_COMPONENT, GL_FLOAT, depthBuffer );

				/* Set the stencil buffer to test for stencil values of 1 */
				glStencilFunc(GL_EQUAL, 1, 1);
				glStencilOp(GL_KEEP, GL_KEEP, GL_KEEP);

				/* Disable the depth testing */
				glDisable(GL_DEPTH_TEST);

				/* Draw the color values to the frame buffer */
				glRasterPos2i( 0, 0 );
				glColorMask( GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE );
				MPI_Recv( imageBuffer, pixelCount * 3, MPI_UNSIGNED_CHAR, neighbourRank, GL_RGB, comm, &status );
				glDrawPixels( width, height, GL_RGB, GL_UNSIGNED_BYTE, imageBuffer );
			}

			/* Change buffer id */
			bufferId = div( bufferId, 2 ).quot;
		}
		mergeCount++;
	}

	if ( broadcast ) {
		/* All pixels are composited onto the master's processor - 
		 * we'll grab these pixels */
		if ( context->rank == MASTER ) {
			glReadPixels( startx, starty, width, height, GL_DEPTH_COMPONENT, GL_FLOAT, depthBuffer);
			glReadPixels( startx, starty, width, height, GL_RGB, GL_UNSIGNED_BYTE, imageBuffer);
		}

		/* Send composited pixel data to other processors */
		MPI_Bcast ( imageBuffer, pixelCount * 3, MPI_UNSIGNED_CHAR, 0, comm );
		MPI_Bcast ( depthBuffer, pixelCount,     MPI_FLOAT,         0, comm );

		if ( context->rank != MASTER ) {
			glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
			glDisable( GL_DEPTH_TEST );
			
			/* Apply Master's composited depth buffer to this processor */
			glRasterPos2i( 0, 0 );
			glDrawPixels( width, height, GL_DEPTH_COMPONENT, GL_FLOAT, depthBuffer );
			
			/* Apply Master's composited image buffer to this processor */
			glRasterPos2i( 0, 0 );
			glDrawPixels( width, height, GL_RGB, GL_UNSIGNED_BYTE, imageBuffer );
		}
	}
	
	/* Reset some opengl stuff */
	glDisable( GL_STENCIL_TEST );
	glEnable(GL_DEPTH_TEST);
	glPopMatrix();
	lucViewportInfo_SetOpenGLCamera( viewportInfo );

	/* Clean up allocated memory */
	Memory_Free( imageBuffer );
	Memory_Free( depthBuffer );
}

void lucRenderingEngineGL_CombineToMaster( 
		void*                                              renderingEngine,
		lucViewportInfo*                                   viewportInfo,
		AbstractContext*                                   context,
		lucPixel*                                          imageBuffer, 
		float*                                             depthBuffer )
{
	lucRenderingEngineGL* self                 = (lucRenderingEngineGL*) renderingEngine;
	Pixel_Index           width                = viewportInfo->width;
	Pixel_Index           height               = viewportInfo->height;
	Index                 buffersToMerge;
	Index                 neighbourRank;
	Index                 bufferId             = context->rank;
	Pixel_Index           pixel_I;
	Pixel_Index           pixelCount           = width*height;
	Index                 mergeCount           = 0;
	float*                neighbourDepthBuffer;
	lucPixel*             neighbourImageBuffer;
	MPI_Status            status;
	MPI_Comm              comm                 = context->communicator;
	
	Journal_DPrintfL( lucDebug, 2, "In func: %s\n", __func__ );
	
	neighbourImageBuffer = Memory_Alloc_Array( lucPixel, pixelCount, "Neighbour's Image Buffer" );
	neighbourDepthBuffer = Memory_Alloc_Array( float,    pixelCount, "Neighbour's Depth Buffer" );

	/* Merge Buffers */
	for ( buffersToMerge = context->nproc - 1 ; buffersToMerge > 0 ; buffersToMerge = div( buffersToMerge, 2 ).quot ) {
		/* If my ID is odd - then send the buffer and get out of loop */
		if ( bufferId % 2 == 1 ) {
			/* Send buffer to left */
			neighbourRank = lucRenderingEngineGL_MapBufferIdToRank( self, bufferId - 1, mergeCount );
			MPI_Send( depthBuffer, pixelCount,     MPI_FLOAT,         neighbourRank, GL_DEPTH_COMPONENT, comm );
			MPI_Send( imageBuffer, pixelCount * 3, MPI_UNSIGNED_CHAR, neighbourRank, GL_RGB,             comm );
			
			/* Now that I've sent my info - I can quit this loop */
			break;
		}
		else {
			/* Only merge if you are not the last processor to have a buffer */
			if ( bufferId < buffersToMerge ) {
				/* Receive Buffer from Right */
				neighbourRank = lucRenderingEngineGL_MapBufferIdToRank( self, bufferId + 1, mergeCount );
				MPI_Recv( neighbourDepthBuffer, pixelCount, MPI_FLOAT, neighbourRank, GL_DEPTH_COMPONENT, comm, &status );
				MPI_Recv( neighbourImageBuffer, pixelCount * 3, MPI_UNSIGNED_CHAR, neighbourRank, GL_RGB, comm, &status );

				/* Merge two buffers */
				for ( pixel_I = 0 ; pixel_I < pixelCount ; pixel_I++ ) {
					if ( neighbourDepthBuffer[ pixel_I ] < depthBuffer[ pixel_I ] ) {
						memcpy( &depthBuffer[ pixel_I ], &neighbourDepthBuffer[ pixel_I ], sizeof( float ) );
						memcpy( &imageBuffer[ pixel_I ], &neighbourImageBuffer[ pixel_I ], sizeof( lucPixel ) );
					}
				}
			}

			/* Change buffer id */
			bufferId = div( bufferId, 2 ).quot;
		}
		mergeCount++;
	}
	
	Memory_Free( neighbourImageBuffer );
	Memory_Free( neighbourDepthBuffer );
	
	Journal_DPrintfL( lucDebug, 2, "Leaving: %s\n", __func__ );
}

void _lucRenderingEngineGL_CompositeViewport_Manual( 
		void*                                              renderingEngine, 
		lucViewportInfo*                                   viewportInfo, 
		AbstractContext*                                   context, 
		Bool                                               broadcast ) 
{
	lucRenderingEngineGL* self              = (lucRenderingEngineGL*) renderingEngine;
	Pixel_Index           width             = viewportInfo->width;
	Pixel_Index           height            = viewportInfo->height;
	Pixel_Index           startx            = viewportInfo->startx;
	Pixel_Index           starty            = viewportInfo->starty;
	lucViewport*          viewport          = viewportInfo->viewport;
	Pixel_Index           pixelCount        = width*height;
	float*                depthBuffer;
	lucPixel*             imageBuffer;
	MPI_Comm              comm              = context->communicator;

	Journal_DPrintfL( lucDebug, 2, "In func: %s - Viewport %s\n", __func__, viewport->name );

	/* Allocate Memory */
	imageBuffer          = Memory_Alloc_Array( lucPixel, pixelCount, "Image Buffer" );
	depthBuffer          = Memory_Alloc_Array( float,    pixelCount, "Depth Buffer" );
	
	/* Read depth buffer */
	glPixelStorei(GL_PACK_ALIGNMENT,1);
	glColorMask(GL_TRUE,GL_TRUE,GL_TRUE,GL_TRUE);
	glReadPixels( startx, starty, width, height, GL_RGB,             GL_UNSIGNED_BYTE, imageBuffer);
	glReadPixels( startx, starty, width, height, GL_DEPTH_COMPONENT, GL_FLOAT,         depthBuffer);

	lucRenderingEngineGL_CombineToMaster( self, viewportInfo, context, imageBuffer, depthBuffer );

	/* Broadcast master's pixels info */
	if (broadcast) {
		MPI_Bcast ( imageBuffer, pixelCount * 3, MPI_UNSIGNED_CHAR, 0, comm );
		MPI_Bcast ( depthBuffer, pixelCount,     MPI_FLOAT,         0, comm );
	}
		
	/* Reset Pixels */
	if (context->rank == MASTER || broadcast) {
		glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
		glPushMatrix();
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		gluOrtho2D( 0.0, (GLfloat) width, 0.0, (GLfloat) height );
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();	

		glPixelStorei(GL_UNPACK_ALIGNMENT,1);
		glRasterPos2i( 0, 0 );
		glDrawPixels( width, height, GL_DEPTH_COMPONENT, GL_FLOAT, depthBuffer );
		
		glRasterPos2i( 0, 0 );
		glDisable(GL_DEPTH_TEST);
		glDrawPixels( width, height, GL_RGB, GL_UNSIGNED_BYTE, imageBuffer );
		glEnable(GL_DEPTH_TEST);
		
		glPopMatrix();
		lucViewportInfo_SetOpenGLCamera( viewportInfo );
	}

	/* Clean up */
	Memory_Free( imageBuffer );
	Memory_Free( depthBuffer );
	
	Journal_DPrintfL( lucDebug, 2, "Leaving: %s\n", __func__ );
}
