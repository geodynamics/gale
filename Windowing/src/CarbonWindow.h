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
#ifdef HAVE_CARBON

#include <Carbon/Carbon.h>
#include <AGL/agl.h>
#ifndef __lucCarbonWindow_h__
#define __lucCarbonWindow_h__

	/** Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
	extern const Type lucCarbonWindow_Type;
		
	/** Class contents - this is defined as a macro so that sub-classes of this class can use this macro at the start of the definition of their struct */
	#define __lucCarbonWindow \
		/* Macro defining parent goes here - This means you can cast this class as its parent */ \
		__lucWindow \
		/* Virtual functions go here */ \
		/* Other info */\
		Pixel_Index                                        offsetX;                  \
		Pixel_Index                                        offsetY;                  \
		void*                                              graphicsContext;          \
		Bool                                               windowIsVisible;          \
		/* Pixel buffer for background windows */									 \
		AGLPbuffer										   PixelBuffer;				 \
		/* Stuff for interactive windows */											 \
		EventHandlerUPP                                    handler;                  \
        EventLoopIdleTimerUPP                              timerHandler;             \
        EventLoopTimerRef								   timer;                    \
		WindowPtr                                          window;                   \
		
	struct lucCarbonWindow { __lucCarbonWindow };
	
	/** Private Constructor: This will accept all the virtual functions for this class as arguments. */
	lucCarbonWindow* _lucCarbonWindow_New( 
		SizeT                                           sizeOfSelf,
		Type                                            type,
		Stg_Class_DeleteFunction*                       _delete,
		Stg_Class_PrintFunction*                        _print,
		Stg_Class_CopyFunction*                         _copy, 
		Stg_Component_DefaultConstructorFunction*       _defaultConstructor,
		Stg_Component_ConstructFunction*                _construct,
		Stg_Component_BuildFunction*                    _build,
		Stg_Component_InitialiseFunction*               _initialise,
		Stg_Component_ExecuteFunction*                  _execute,
		Stg_Component_DestroyFunction*                  _destroy,
		lucWindow_DisplayFunction*						_displayWindow,	
		lucWindow_EventsWaitingFunction*				_eventsWaiting,	
		lucWindow_EventProcessorFunction*				_eventProcessor,	
		lucWindow_ResizeFunction*						_resizeWindow,	
		Name                                            name );

	void _lucCarbonWindow__Delete( void* window ) ;
	void _lucCarbonWindow_Print( void* window, Stream* stream ) ;
	void* _lucCarbonWindow_Copy( void* window, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap) ;

	/* 'Stg_Component' implementations */
	void* _lucCarbonWindow_DefaultNew( Name name ) ;
	void _lucCarbonWindow_AssignFromXML( void* window, Stg_ComponentFactory* cf, void* data );
	void _lucCarbonWindow_Build( void* window, void* data ) ;
	void _lucCarbonWindow_Initialise( void* window, void* data ) ;
	void _lucCarbonWindow_Execute( void* window, void* data );
	void _lucCarbonWindow_Destroy( void* window, void* data ) ;

	/* Window Virtuals */
	void _lucCarbonWindow_Display( void* window );
	int _lucCarbonWindow_EventsWaiting( void* window ) ;
	Bool _lucCarbonWindow_EventProcessor( void* window ) ;
	void _lucCarbonWindow_Resize( void* window );

	void lucCarbonWindow_CreateWindow( void* window ) ;
	void lucCarbonWindow_DestroyWindow( void* window ) ;
	
	void lucCarbonWindow_Draw( void* window ) ;
	pascal void lucCarbonWindow_IdleTimer(EventLoopTimerRef inTimer, EventLoopIdleTimerMessage inState, void * inUserData);
#endif

#endif /* HAVE_CARBON */
