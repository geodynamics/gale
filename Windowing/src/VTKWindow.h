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

#ifdef HAVE_VTK

#ifndef __lucVTKWindow_h__
#define __lucVTKWindow_h__


	/** Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
	extern const Type lucVTKWindow_Type;
		
	/** Class contents - this is defined as a macro so that sub-classes of this class can use this macro at the start of the definition of their struct */
	#define __lucVTKWindow \
		/* Macro defining parent goes here - This means you can cast this class as its parent */ \
		__lucWindow \
		/* Virtual functions go here */ \
		/* Other info */\
		/*vtkRenderWindow*  */\
		void*						                                    rendererWindow;          \
		/* vtkRenderer* */\
		void*                                               renderer;                \
		/*vtkXRenderWindowInteractor* */\
		void*                                               interactor;              \
		char*                                               title;                   \
		Bool                                                fullscreen;              \
		int                                                 subframes;               \
		/*Stuff from X11Window*/                                                     \
		double                                              maxIdleTime;             \
		/* Display name stuff */ \
		Name                                                displayName;             \
		Name                                                host;                    \
		unsigned int                                        displayNumber;           \
		unsigned int                                        displayScreen;           \

	struct lucVTKWindow { __lucVTKWindow };
	
	/** Private Constructor: This will accept all the virtual functions for this class as arguments. */
	lucVTKWindow* _lucVTKWindow_New( 
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

	void _lucVTKWindow_Delete( void* window ) ;
	void _lucVTKWindow_Print( void* window, Stream* stream ) ;
	void* _lucVTKWindow_Copy( void* window, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap) ;

	/* 'Stg_Component' implementations */
	void* _lucVTKWindow_DefaultNew( Name name ) ;
	void _lucVTKWindow_Construct( void* window, Stg_ComponentFactory* cf, void* data );
	void _lucVTKWindow_Build( void* window, void* data ) ;
	void _lucVTKWindow_Initialise( void* window, void* data ) ;
	void _lucVTKWindow_Execute( void* window, void* data );
	void _lucVTKWindow_Destroy( void* window, void* data ) ;

	Bool lucVTKWindow_EventLoop( void* window, AbstractContext* context) ;

	void lucVTKWindow_SetupFonts( void* window ) ;
	void lucVTKWindow_SwapBuffers( void* window, AbstractContext* context ) ;
	void lucVTKWindow_MakeCurrent( void* window )  ;

	Bool lucVTKWindow_CreateDisplay( void* window )  ;
	void lucVTKWindow_CreateBackgroundWindow( void* window )  ;
	void lucVTKWindow_CreateInteractiveWindow( void* window ) ;

	//Colormap lucVTKWindow_GetShareableColormap( lucVTKWindow* self ) ;
	void lucVTKWindow_CloseInteractiveWindow( lucVTKWindow* self ) ;
	void lucVTKWindow_CloseBackgroundWindow( lucVTKWindow* self ) ;
	void lucVTKWindow_CloseDisplay( lucVTKWindow* self ) ;
	Bool lucVTKWindow_FindFont( void* window )  ;

#endif

#endif /* HAVE_VTK */
