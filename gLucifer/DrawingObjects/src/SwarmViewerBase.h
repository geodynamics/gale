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
** $Id: SwarmViewer.h 594 2006-07-13 02:17:37Z CecileDuboz $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#ifndef __lucSwarmViewerBase_h__
#define __lucSwarmViewerBase_h__

	typedef void (lucSwarmViewerBase_PlotParticleFunction) ( void* object, void* context, Particle_Index lParticle_I );
	typedef void (lucSwarmViewerBase_SetParticleColourFunction) ( void* object, void* context, Particle_Index lParticle_I );

	/** Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
	extern const Type lucSwarmViewerBase_Type;
		
	/** Class contents - this is defined as a macro so that sub-classes of this class can use this macro at the start of the definition of their struct */
	#define __lucSwarmViewerBase \
		/* Macro defining parent goes here - This means you can cast this class as its parent */ \
		__lucOpenGLDrawingObject \
		/* Virtual functions go here */ \
		lucSwarmViewerBase_PlotParticleFunction*           _plotParticle;          \
		lucSwarmViewerBase_SetParticleColourFunction*      _setParticleColour;     \
		/* Colour stuff */\
		lucColour                                          colour;                 \
		Name                                               colourVariableName;     \
		SwarmVariable*                                     colourVariable;         \
		lucColourMap*                                      colourMap;			\
		/* Other info */\
		Swarm*                                             swarm;                  \
		/* Opacity Stuff */ \
		Name                                               opacityVariableName;    \
		SwarmVariable*                                     opacityVariable;        \
		/* Mask Info */ \
		Name                                               maskVariableName;       \
		SwarmVariable*                                     maskVariable;           \
		lucDrawingObjectMask                               mask;                   \
		/* Other Stuff */ \
		Bool                                               drawParticleNumber;     \
		Bool                                               sameParticleColour;     \
		int                                                subSetEvery;            \
		Bool                                               positionRange;          \
		Coord                                              minPosition;            \
		Coord                                              maxPosition;

	struct lucSwarmViewerBase { __lucSwarmViewerBase };
	
	/** Private Constructor: This will accept all the virtual functions for this class as arguments. */
	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define LUCSWARMVIEWERBASE_DEFARGS \
                LUCOPENGLDRAWINGOBJECT_DEFARGS, \
                lucSwarmViewerBase_PlotParticleFunction*            _plotParticle, \
                lucSwarmViewerBase_SetParticleColourFunction*  _setParticleColour

	#define LUCSWARMVIEWERBASE_PASSARGS \
                LUCOPENGLDRAWINGOBJECT_PASSARGS, \
	        _plotParticle,      \
	        _setParticleColour

	lucSwarmViewerBase* _lucSwarmViewerBase_New(  LUCSWARMVIEWERBASE_DEFARGS  );

	void _lucSwarmViewerBase_Delete( void* drawingObject ) ;
	void _lucSwarmViewerBase_Print( void* drawingObject, Stream* stream ) ;
	void* _lucSwarmViewerBase_Copy( void* drawingObject, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap) ;

	/* 'Stg_Component' implementations */
	void _lucSwarmViewerBase_AssignFromXML( void* drawingObject, Stg_ComponentFactory* cf, void* data );
	void _lucSwarmViewerBase_Build( void* drawingObject, void* data ) ;
	void _lucSwarmViewerBase_Initialise( void* drawingObject, void* data ) ;
	void _lucSwarmViewerBase_Execute( void* drawingObject, void* data );
	void _lucSwarmViewerBase_Destroy( void* drawingObject, void* data ) ;
	
	void _lucSwarmViewerBase_Setup( void* drawingObject, void* _context ) ;
	void _lucSwarmViewerBase_Draw( void* drawingObject, lucWindow* window, lucViewportInfo* viewportInfo, void* _context ) ;
	void _lucSwarmViewerBase_CleanUp( void* drawingObject, void* context ) ;

	void _lucSwarmViewerBase_BuildDisplayList( void* drawingObject, void* _context ) ;

	void lucSwarmViewBase_DrawParticleNumbers( void* drawingObject, void* _context ) ;

	void lucSwarmViewerBase_UpdateVariables( void* drawingObject ) ;

	void lucSwarmViewerBase_FindParticleLocalIndex(void *drawingObject, Coord coord, Particle_Index  *lParticle_I);

	void _lucSwarmViewerBase_SetParticleColourDefault( void* drawingObject, void* _context, Particle_Index lParticle_I ) ;
#endif

