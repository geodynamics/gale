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


#ifndef __lucWindowInteraction_h__
#define __lucWindowInteraction_h__
	
	typedef void (lucWindowInteraction_MouseMotionFunction) ( void* windowInteractor, lucWindow* window, lucMouseButton button, Pixel_Index xpos, Pixel_Index ypos, Pixel_Index startx, Pixel_Index starty );
	typedef void (lucWindowInteraction_MouseClickFunction) ( void* windowInteractor, lucWindow* window, lucMouseButton button, lucMouseState state, Pixel_Index xpos, Pixel_Index ypos );
	typedef void (lucWindowInteraction_MouseMessageFunction)  ( void* windowInteractor, Stream* stream );
	typedef void (lucWindowInteraction_KeyboardEventFunction) ( void* windowInteractor, lucWindow* window, char key, Pixel_Index xpos, Pixel_Index ypos );
	typedef void (lucWindowInteraction_KeyboardMessageFunction) ( void* windowInteractor, Stream* stream );

	/** Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
	extern const Type lucWindowInteraction_Type;

	#define __lucWindowInteraction                                      \
		/* Parent macro */ \
		__Stg_Component                                                      \
		AbstractContext*				   context;	     \
		/* Virtual functions go here */ \
		lucWindowInteraction_MouseMotionFunction*          _mouseMotion;     \
		lucWindowInteraction_MouseClickFunction*           _mouseClick;      \
		lucWindowInteraction_MouseMessageFunction*         _mouseMessage;    \
		lucWindowInteraction_KeyboardEventFunction*        _keyboardEvent;   \
		lucWindowInteraction_KeyboardMessageFunction*      _keyboardMessage; \
		\
		/* Other Info */ \
			
	struct lucWindowInteraction {__lucWindowInteraction};

	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define LUCWINDOWINTERACTION_DEFARGS \
                STG_COMPONENT_DEFARGS, \
                lucWindowInteraction_MouseMotionFunction*          _mouseMotion, \
                lucWindowInteraction_MouseClickFunction*            _mouseClick, \
                lucWindowInteraction_MouseMessageFunction*        _mouseMessage, \
                lucWindowInteraction_KeyboardEventFunction*      _keyboardEvent, \
                lucWindowInteraction_KeyboardMessageFunction*  _keyboardMessage

	#define LUCWINDOWINTERACTION_PASSARGS \
                STG_COMPONENT_PASSARGS, \
	        _mouseMotion,     \
	        _mouseClick,      \
	        _mouseMessage,    \
	        _keyboardEvent,   \
	        _keyboardMessage

	lucWindowInteraction* lucWindowInteraction_New( Name name ) ;
	lucWindowInteraction* _lucWindowInteraction_New(  LUCWINDOWINTERACTION_DEFARGS  );

	void lucWindowInteraction_InitAll( void* windowInteractor ) ;

	void _lucWindowInteraction_Delete( void* windowInteraction ) ;
	void _lucWindowInteraction_Print( void* windowInteraction, Stream* stream ) ;
	void* _lucWindowInteraction_Copy( void* windowInteraction, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) ;

	/* Stg_Component Virtual Function Implementations */
	void* _lucWindowInteraction_DefaultNew( Name name ) ;
	void _lucWindowInteraction_AssignFromXML( void* windowInteraction, Stg_ComponentFactory* cf, void* data ) ;
	void _lucWindowInteraction_Build( void* windowInteraction, void* data ) ;
	void _lucWindowInteraction_Initialise( void* windowInteraction, void* data ) ;
	void _lucWindowInteraction_Execute( void* windowInteraction, void* data ) ;
	void _lucWindowInteraction_Destroy( void* windowInteraction, void* data ) ;

	/* +++ Public Functions +++ */

	/* Wrappers to virtual functions */
	void lucWindowInteraction_MouseMotion( void* windowInteractor, lucWindow* window, lucMouseButton button, Pixel_Index xpos, Pixel_Index ypos, Pixel_Index startx, Pixel_Index starty) ;
	void lucWindowInteraction_MouseClick( void* windowInteractor, lucWindow* window, lucMouseButton button, lucMouseState state, Pixel_Index xpos, Pixel_Index ypos) ;
	void lucWindowInteraction_MouseMessage( void* windowInteractor, Stream* stream ) ;
	void lucWindowInteraction_KeyboardEvent( void* windowInteractor, lucWindow* window, char key, Pixel_Index xpos, Pixel_Index ypos) ;
	void lucWindowInteraction_KeyboardMessage( void* windowInteraction, Stream* stream ) ;

	/** Default interactivity implementations */
	void _lucWindowInteraction_MouseMotion( void* windowInteraction, lucWindow* window, lucMouseButton button, Pixel_Index xpos, Pixel_Index ypos, Pixel_Index startx, Pixel_Index starty) ;
	void _lucWindowInteraction_MouseClick( void* windowInteraction, lucWindow* window, lucMouseButton button, lucMouseState state, Pixel_Index xpos, Pixel_Index ypos) ;
	void _lucWindowInteraction_MouseMessage( void* windowInteractor, Stream* stream ) ;
	void _lucWindowInteraction_KeyboardEvent( void* windowInteraction, lucWindow* window, char key, Pixel_Index xpos, Pixel_Index ypos) ;
	void _lucWindowInteraction_KeyboardMessage( void* windowInteractor, Stream* stream ) ;

#endif

