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
** $Id: ColourBarInteraction.h 510 2006-02-17 04:33:32Z RobertTurnbull $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#ifndef __lucColourBarInteraction_h__
#define __lucColourBarInteraction_h__

	/** Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
	extern const Type lucColourBarInteraction_Type;
		
	/** Class contents - this is defined as a macro so that sub-classes of this class can use this macro at the start of the definition of their struct */
	#define __lucColourBarInteraction \
		/* Macro defining parent goes here - This means you can cast this class as its parent */ \
		__lucWindowInteraction \
		/* Virtual functions go here */ \
		/* Other info */\

	struct lucColourBarInteraction { __lucColourBarInteraction };
	
	/** Private Constructor: This will accept all the virtual functions for this class as arguments. */
	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define LUCCOLOURBARINTERACTION_DEFARGS \
                LUCWINDOWINTERACTION_DEFARGS

	#define LUCCOLOURBARINTERACTION_PASSARGS \
                LUCWINDOWINTERACTION_PASSARGS

	lucColourBarInteraction* _lucColourBarInteraction_New(  LUCCOLOURBARINTERACTION_DEFARGS  );

	void _lucColourBarInteraction_Delete( void* renderingEngine ) ;
	void _lucColourBarInteraction_Print( void* renderingEngine, Stream* stream ) ;
	void* _lucColourBarInteraction_Copy( void* renderingEngine, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap) ;

	/* 'Stg_Component' implementations */
	void* _lucColourBarInteraction_DefaultNew( Name name ) ;
	void _lucColourBarInteraction_AssignFromXML( void* renderingEngine, Stg_ComponentFactory* cf, void* data );
	void _lucColourBarInteraction_Build( void* renderingEngine, void* data ) ;
	void _lucColourBarInteraction_Initialise( void* renderingEngine, void* data ) ;
	void _lucColourBarInteraction_Execute( void* renderingEngine, void* data );
	void _lucColourBarInteraction_Destroy( void* renderingEngine, void* data ) ;
	
	void _lucColourBarInteraction_MouseMotion( void* windowInteraction, lucWindow* window, lucMouseButton button, Pixel_Index xpos, Pixel_Index ypos, Pixel_Index startx, Pixel_Index starty) ;
	void _lucColourBarInteraction_MouseClick( void* windowInteraction, lucWindow* window, lucMouseButton button, lucMouseState state, Pixel_Index xpos, Pixel_Index ypos) ;
	void _lucColourBarInteraction_MouseMessage( void* windowInteractor, Stream* stream ) ;
	void _lucColourBarInteraction_KeyboardEvent( void* windowInteraction, lucWindow* window, char key, Pixel_Index xpos, Pixel_Index ypos) ;
	void _lucColourBarInteraction_KeyboardMessage( void* windowInteractor, Stream* stream ) ;
	
#endif

