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
** $Id: SDLWindow.h 628 2006-10-12 08:23:07Z SteveQuenette $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifdef HAVE_SDL

#include <SDL/SDL.h>

#ifndef __lucSDLWindow_h__
#define __lucSDLWindow_h__

	/** Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
	extern const Type lucSDLWindow_Type;
		
	/** Class contents - this is defined as a macro so that sub-classes of this class can use this macro at the start of the definition of their struct */
	#define __lucSDLWindow \
		/* Macro defining parent goes here - This means you can cast this class as its parent */ \
		__lucWindow \
		/* Virtual functions go here */ \
		/* Other info */\
		int													sdlFlags;				\
		SDL_Surface*										screen;					\
		SDL_Surface*										buffer;					\
		void*                                               osMesaContext;			\
		SDL_TimerID											timer;					\

	struct lucSDLWindow { __lucSDLWindow };
	
	/** Private Constructor: This will accept all the virtual functions for this class as arguments. */
	lucSDLWindow* _lucSDLWindow_New( 
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
		lucWindow_DisplayFunction*						   _displayWindow,	
		lucWindow_EventsWaitingFunction*				   _eventsWaiting,	
		lucWindow_EventProcessorFunction*				   _eventProcessor,	
		Name                                               name );

	void _lucSDLWindow_Delete( void* window ) ;
	void _lucSDLWindow_Print( void* window, Stream* stream ) ;
	void* _lucSDLWindow_Copy( void* window, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap) ;

	/* 'Stg_Component' implementations */
	void* _lucSDLWindow_DefaultNew( Name name ) ;
	void _lucSDLWindow_Construct( void* window, Stg_ComponentFactory* cf, void* data );
	void _lucSDLWindow_Build( void* window, void* data ) ;
	void _lucSDLWindow_Initialise( void* window, void* data ) ;
	void _lucSDLWindow_Execute( void* window, void* data );
	void _lucSDLWindow_Destroy( void* window, void* data ) ;

	/* Window Virtuals */
	void _lucSDLWindow_Display( void* window );
	int _lucSDLWindow_EventsWaiting( void* window ) ;
	Bool _lucSDLWindow_EventProcessor( void* window ) ;

	/* Resize video */	
	void lucSDLWindow_Resize( void* window, Pixel_Index width, Pixel_Index height ) ;

	/* Timer callback */
	Uint32 lucSDLWindow_IdleTimer(Uint32 interval, void* param);

#endif

#endif /* HAVE_SDL */
