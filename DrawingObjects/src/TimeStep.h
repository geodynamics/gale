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
** $Id: timeStep.h 510 2006-02-17 04:33:32Z RobertTurnbull $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/



#ifndef __lucTimeStep_h__
#define __lucTimeStep_h__

	extern const Type lucTimeStep_Type;

	#define __lucTimeStep                              \
		__lucOpenGLDrawingObject \
		lucColour                           colour; \
		Bool                                frame; \
		Bool                                time;
		
	struct lucTimeStep {__lucTimeStep};

	/** Constructors */
	lucTimeStep* lucTimeStep_New( 
		Name                                               name,
		lucColour                                          colour,
		Bool                                               frame,
		Bool                                               time);

	lucTimeStep* _lucTimeStep_New(
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
		lucDrawingObject_SetupFunction*                    _setup,
		lucDrawingObject_DrawFunction*                     _draw,
		lucDrawingObject_CleanUpFunction*                  _cleanUp,
		lucOpenGLDrawingObject_BuildDisplayListFunction*   _buildDisplayList,
		Name                                               name );

	void lucTimeStep_InitAll( 
		void*                                              timeStep,
		lucColour                                          colour,
		Bool                                               frame,
		Bool                                               time);

	void _lucTimeStep_Init( 
		void*                                              timeStep,
		lucColour                                          colour,
		Bool                                               frame,
		Bool                                               time);


	void _lucTimeStep_Setup( void* drawingObject, void* _context );
		
	/** Virtual Functions */
	void _lucTimeStep_Delete( void* timeStep ) ;
	void _lucTimeStep_Print( void* timeStep, Stream* stream ) ;
	void* _luctTimeStep_Copy( void* timeStep, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) ;
	#define lucTimeStep_Copy( self ) \
		(lucTimeStep*) Stg_Class_Copy( self, NULL, False, NULL, NULL )
	void* _lucTimeStep_DefaultNew( Name name ) ;
void _lucTimeStep_Construct( void* timeStep, Stg_ComponentFactory* cf, void* data ) ;
	void _lucTimeStep_Build( void* timeStep, void* data );
	void _lucTimeStep_Initialise( void* timeStep, void* data );
	void _lucTimeStep_Execute( void* timeStep, void* data );
	void _lucTimeStep_Destroy( void* timeStep, void* data );
	void _lucTimeStep_Draw( void* drawingObject, lucWindow* window, lucViewportInfo* viewportInfo, void* _context );
	void _lucTimeStep_CleanUp( void* drawingObject, void* _context );
	void _lucTimeStep_BuildDisplayList( void* drawingObject, void* _context );

#endif
