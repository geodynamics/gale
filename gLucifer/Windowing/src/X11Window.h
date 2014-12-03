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
** $Id: X11Window.h 628 2006-10-12 08:23:07Z SteveQuenette $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifdef HAVE_X11

#ifndef __lucX11Window_h__
#define __lucX11Window_h__


	#include <GL/glx.h>
	#include <X11/Xlib.h>
	#include <X11/Xatom.h>
	#include <X11/Xmu/StdCmap.h>
	#include <X11/keysym.h>

	/** Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
	extern const Type lucX11Window_Type;
		
	/** Class contents - this is defined as a macro so that sub-classes of this class can use this macro at the start of the definition of their struct */
	#define __lucX11Window \
		/* Macro defining parent goes here - This means you can cast this class as its parent */ \
		__lucWindow \
		/* Virtual functions go here */ \
		/* Other info */\
		GLXContext                                          glxcontext;              \
		Display*                                            display;                 \
		Window                                              win;                     \
		Pixmap                                              pmap;                    \
		GLXPixmap                                           glxpmap;                 \
		XVisualInfo*                                        vi;                      \
		Bool                                                doubleBuffer;            \
		Atom             									wmDeleteWindow;			 \
		/* Display name stuff */ \
		Name                                                displayName;             \
		Name                                                host;                    \
		unsigned int                                        displayNumber;           \
		unsigned int                                        displayScreen;           \

	struct lucX11Window { __lucX11Window };
	
	/** Private Constructor: This will accept all the virtual functions for this class as arguments. */
	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define LUCX11WINDOW_DEFARGS \
                LUCWINDOW_DEFARGS

	#define LUCX11WINDOW_PASSARGS \
                LUCWINDOW_PASSARGS

	lucX11Window* _lucX11Window_New(  LUCX11WINDOW_DEFARGS  );

	void _lucX11Window_Delete( void* window ) ;
	void _lucX11Window_Print( void* window, Stream* stream ) ;
	void* _lucX11Window_Copy( void* window, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap) ;

	/* 'Stg_Component' implementations */
	void* _lucX11Window_DefaultNew( Name name ) ;
	void _lucX11Window_AssignFromXML( void* window, Stg_ComponentFactory* cf, void* data );
	void _lucX11Window_Build( void* window, void* data ) ;
	void _lucX11Window_Initialise( void* window, void* data ) ;
	void _lucX11Window_Execute( void* window, void* data );
	void _lucX11Window_Destroy( void* window, void* data ) ;

	/* Window Virtual implementations */
	void _lucX11Window_Display( void* window );
	int _lucX11Window_EventsWaiting( void* window ) ;
	Bool _lucX11Window_EventProcessor( void* window ) ;
	void _lucX11Window_Resize( void* window );
	
	Bool lucX11Window_CreateDisplay( void* window )  ;
	void lucX11Window_CreateBackgroundWindow( void* window )  ;
	void lucX11Window_CreateInteractiveWindow( void* window ) ;

	Colormap lucX11Window_GetShareableColormap( lucX11Window* self ) ;
	void lucX11Window_CloseInteractiveWindow( lucX11Window* self ) ;
	void lucX11Window_CloseBackgroundWindow( lucX11Window* self ) ;

	void lucX11Window_Timer( int ) ;
    int lucX11Window_Error(Display* display, XErrorEvent* error);

#endif

#endif /* HAVE_X11 */

