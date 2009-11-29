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


#ifndef __lucSwarmRGBColourViewer_h__
#define __lucSwarmRGBColourViewer_h__

	/** Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
	extern const Type lucSwarmRGBColourViewer_Type;
		
	/** Class contents - this is defined as a macro so that sub-classes of this class can use this macro at the start of the definition of their struct */
	#define __lucSwarmRGBColourViewer \
		/* Macro defining parent goes here - This means you can cast this class as its parent */ \
		__lucSwarmViewer \
		/* Virtual functions go here */ \
		/* Other info */\
		ExtensionInfo_Index                                particleExtHandle;         \
		/* Colour Info */ \
		Name                                               colourRedVariableName;     \
		SwarmVariable*                                     colourRedVariable;         \
		Name                                               colourGreenVariableName;   \
		SwarmVariable*                                     colourGreenVariable;       \
		Name                                               colourBlueVariableName;    \
		SwarmVariable*                                     colourBlueVariable;        \
	
	struct lucSwarmRGBColourViewer { __lucSwarmRGBColourViewer };
	
	/** Private Constructor: This will accept all the virtual functions for this class as arguments. */
	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define LUCSWARMRGBCOLOURVIEWER_DEFARGS \
                LUCSWARMVIEWER_DEFARGS

	#define LUCSWARMRGBCOLOURVIEWER_PASSARGS \
                LUCSWARMVIEWER_PASSARGS

	lucSwarmRGBColourViewer* _lucSwarmRGBColourViewer_New(  LUCSWARMRGBCOLOURVIEWER_DEFARGS  );

	void _lucSwarmRGBColourViewer_Delete( void* drawingObject ) ;
	void _lucSwarmRGBColourViewer_Print( void* drawingObject, Stream* stream ) ;
	void* _lucSwarmRGBColourViewer_Copy( void* drawingObject, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap) ;

	/* 'Stg_Component' implementations */
	void* _lucSwarmRGBColourViewer_DefaultNew( Name name ) ;
	void _lucSwarmRGBColourViewer_AssignFromXML( void* drawingObject, Stg_ComponentFactory* cf, void* data );
	void _lucSwarmRGBColourViewer_Build( void* drawingObject, void* data ) ;
	void _lucSwarmRGBColourViewer_Initialise( void* drawingObject, void* data ) ;
	void _lucSwarmRGBColourViewer_Execute( void* drawingObject, void* data );
	void _lucSwarmRGBColourViewer_Destroy( void* drawingObject, void* data ) ;
	
	void _lucSwarmRGBColourViewer_Setup( void* drawingObject, void* _context ) ;
	void _lucSwarmRGBColourViewer_Draw( void* drawingObject, lucWindow* window, lucViewportInfo* viewportInfo, void* _context ) ;
	void _lucSwarmRGBColourViewer_CleanUp( void* drawingObject, void* context ) ;

	void _lucSwarmRGBColourViewer_BuildDisplayList( void* drawingObject, void* _context ) ;

	void _lucSwarmRGBColourViewer_PlotParticle( void* drawingObject, void* _context, Particle_Index lParticle_I ) ;
	
	/* SwarmViewer Implementations */	
	void _lucSwarmRGBColourViewer_SetParticleColour( void* drawingObject, void* _context, Particle_Index lParticle_I ) ;

	void lucSwarmRGBColourViewer_UpdateVariables( void* drawingObject ) ;
#endif

