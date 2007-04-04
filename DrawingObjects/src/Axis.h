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
** $Id: Axis.h 510 2006-02-17 04:33:32Z RobertTurnbull $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/



#ifndef __lucAxis_h__
#define __lucAxis_h__

	extern const Type lucAxis_Type;

	#define __lucAxis                              \
		__lucOpenGLDrawingObject \
		Coord                               origin;\
		float 				    scale;\
		lucColour                           colourX;\
		lucColour                           colourY;\
		lucColour                           colourZ;
		
	struct lucAxis {__lucAxis};

	/** Constructors */
	lucAxis* lucAxis_New( 
		Name                                               name,
		Coord                                              origin,
		float 						   scale,
		lucColour                                          colourX,
		lucColour                                          colourY,
		lucColour                                          colourZ);

	lucAxis* _lucAxis_New(
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

	void lucAxis_InitAll( 
		void*                                              axis,
		Coord                                              origin,
		float                                              scale, 
		lucColour                                          colourX,
		lucColour                                          colourY,
		lucColour                                          colourZ);

	void _lucAxis_Init( 
		void*                                              axis,
		Coord                                              origin,
		float 						   scale, 
		lucColour                                          colourX,
		lucColour                                          colourY,
		lucColour                                          colourZ);


	void _lucAxis_Setup( void* drawingObject, void* _context );
		
	/** Virtual Functions */
	void _lucAxis_Delete( void* axis ) ;
	void _lucAxis_Print( void* axis, Stream* stream ) ;
	void* _lucAxis_Copy( void* axis, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) ;
	#define lucAxis_Copy( self ) \
		(lucAxis*) Stg_Class_Copy( self, NULL, False, NULL, NULL )
	void* _lucAxis_DefaultNew( Name name ) ;
void _lucAxis_Construct( void* axis, Stg_ComponentFactory* cf, void* data ) ;
	void _lucAxis_Build( void* axis, void* data );
	void _lucAxis_Initialise( void* axis, void* data );
	void _lucAxis_Execute( void* axis, void* data );
	void _lucAxis_Destroy( void* axis, void* data );
	void _lucAxis_Draw( void* drawingObject, lucWindow* window, lucViewportInfo* viewportInfo, void* _context );
	void _lucAxis_CleanUp( void* drawingObject, void* _context );
	void _lucAxis_BuildDisplayList( void* drawingObject, void* _context );

#endif
