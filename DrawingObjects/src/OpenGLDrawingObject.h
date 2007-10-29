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
** $Id: OpenGLDrawingObject.h 746 2007-10-29 04:26:41Z RobertTurnbull $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#ifndef __lucOpenGLDrawingObject_h__
#define __lucOpenGLDrawingObject_h__
	
	typedef void (lucOpenGLDrawingObject_BuildDisplayListFunction) ( void* object, void* context );

	/** Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
	extern const Type lucOpenGLDrawingObject_Type;

	/** Class contents - this is defined as a macro so that sub-classes of this class can use this macro at the start of the definition of their struct */
	#define __lucOpenGLDrawingObject \
		/* Macro defining parent goes here - This means you can cast this class as its parent */ \
		__lucDrawingObject \
		/* Virtual functions go here */ \
		lucOpenGLDrawingObject_BuildDisplayListFunction*   _buildDisplayList;            \
		/* Other info */\
		unsigned int                                       displayList;

	struct lucOpenGLDrawingObject { __lucOpenGLDrawingObject };
	
	/** Private Constructor: This will accept all the virtual functions for this class as arguments. */
	lucOpenGLDrawingObject* _lucOpenGLDrawingObject_New( 
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

	void _lucOpenGLDrawingObject_Delete( void* drawingObject ) ;
	void _lucOpenGLDrawingObject_Print( void* drawingObject, Stream* stream ) ;
	void* _lucOpenGLDrawingObject_Copy( void* drawingObject, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap) ;

	/* 'Stg_Component' implementations */
	void _lucOpenGLDrawingObject_Construct( void* drawingObject, Stg_ComponentFactory* cf, void* data );
	void _lucOpenGLDrawingObject_Build( void* drawingObject, void* data ) ;
	void _lucOpenGLDrawingObject_Initialise( void* drawingObject, void* data ) ;
	void _lucOpenGLDrawingObject_Execute( void* drawingObject, void* data );
	void _lucOpenGLDrawingObject_Destroy( void* drawingObject, void* data ) ;
	
	void _lucOpenGLDrawingObject_Setup( void* drawingObject, void* _context ) ;
	void _lucOpenGLDrawingObject_Draw( void* drawingObject, lucWindow* window, lucViewportInfo* viewportInfo, void* _context ) ;
	void _lucOpenGLDrawingObject_CleanUp( void* drawingObject, void* _context ) ;

	void lucOpenGLDrawingObject_BuildDisplayList( void* drawingObject, void* context ) ;

	/* HACK - a function to check whether a field is an FeVariable or not before it does an FeVariable_SyncShadowValues */
	void lucOpenGLDrawingObject_SyncShadowValues( void* drawingObject, void* field );

#endif
