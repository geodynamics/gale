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
** $Id: SwarmVectors.h 667 2007-03-01 06:43:11Z RobertTurnbull $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#ifndef __lucSwarmVectors_h__
#define __lucSwarmVectors_h__

	/** Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
	extern const Type lucSwarmVectors_Type;
		
	/** Class contents - this is defined as a macro so that sub-classes of this class can use this macro at the start of the definition of their struct */
	#define __lucSwarmVectors \
		/* Macro defining parent goes here - This means you can cast this class as its parent */ \
		__lucSwarmViewerBase \
		/* Virtual functions go here */ \
		/* Other info */\
		SwarmVariable*                                     directionVariable;           \
		Name                                               directionVariableName;       \
		double                                             arrowHeadSize;               \
		SwarmVariable*                                     thicknessVariable;           \
		Name                                               thicknessVariableName;       \
		double                                             thickness;                   \
		SwarmVariable*                                     lengthVariable;              \
		Name                                               lengthVariableName;          \
		double                                             length;                      \

	struct lucSwarmVectors { __lucSwarmVectors };
	
	/** Private Constructor: This will accept all the virtual functions for this class as arguments. */
	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define LUCSWARMVECTORS_DEFARGS \
                LUCSWARMVIEWERBASE_DEFARGS

	#define LUCSWARMVECTORS_PASSARGS \
                LUCSWARMVIEWERBASE_PASSARGS

	lucSwarmVectors* _lucSwarmVectors_New(  LUCSWARMVECTORS_DEFARGS  );

	void _lucSwarmVectors_Delete( void* drawingObject ) ;
	void _lucSwarmVectors_Print( void* drawingObject, Stream* stream ) ;
	void* _lucSwarmVectors_Copy( void* drawingObject, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap) ;

	/* 'Stg_Component' implementations */
	void* _lucSwarmVectors_DefaultNew( Name name ) ;
	void _lucSwarmVectors_AssignFromXML( void* drawingObject, Stg_ComponentFactory* cf, void* data );
	void _lucSwarmVectors_Build( void* drawingObject, void* data ) ;
	void _lucSwarmVectors_Initialise( void* drawingObject, void* data ) ;
	void _lucSwarmVectors_Execute( void* drawingObject, void* data );
	void _lucSwarmVectors_Destroy( void* drawingObject, void* data ) ;
	
	/* Drawing Object Implementations */
	void _lucSwarmVectors_Setup( void* drawingObject, void* _context ) ;
	void _lucSwarmVectors_Draw( void* drawingObject, lucWindow* window, lucViewportInfo* viewportInfo, void* _context ) ;
	void _lucSwarmVectors_CleanUp( void* drawingObject, void* _context ) ;

	/* OpenGL Drawing Object Implementations */
	void _lucSwarmVectors_BuildDisplayList( void* drawingObject, void* _context ) ;

	/* SwarmViewer Base Implementations */
	void _lucSwarmVectors_PlotParticle( void* drawingObject, void* _context, Particle_Index lParticle_I ) ;

	void lucSwarmVectors_UpdateVariables( void* drawingObject ) ;
#endif

