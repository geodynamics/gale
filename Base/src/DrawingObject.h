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
** $Id: DrawingObject.h 628 2006-10-12 08:23:07Z SteveQuenette $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#ifndef __lucDrawingObject_h__
#define __lucDrawingObject_h__

	extern const Type lucDrawingObject_Type;

	typedef void (lucDrawingObject_SetupFunction) ( void* object, void* context );
	typedef void (lucDrawingObject_DrawFunction)  ( void* object, lucWindow* window, lucViewportInfo* viewportInfo, void* context );
	typedef void (lucDrawingObject_CleanUpFunction)  ( void* object, void* context );

	/** Class contents - this is defined as a macro so that sub-classes of this class can use this macro at the start of the definition of their struct */
	#define __lucDrawingObject                           \
		/* Macro defining parent goes here - This means you can cast this class as its parent */ \
		__Stg_Component                                   \
		AbstractContext*				   context;		       \
		/* Virtual Functions */ \
		lucDrawingObject_SetupFunction*                    _setup;                     \
		lucDrawingObject_DrawFunction*                     _draw;                      \
		lucDrawingObject_CleanUpFunction*                  _cleanUp;                   \
		/* Other Info */ \
		Bool                                               needsToSetup;               \
		Bool                                               needsToCleanUp;             \
		/* Journal Information */ \
		Stream*                                            infoStream;                 \
		Stream*                                            errorStream;                \
		Stream*                                            debugStream;                \
		

	struct lucDrawingObject {__lucDrawingObject};

	lucDrawingObject* _lucDrawingObject_New(
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
		Name                                               name );


	void _lucDrawingObject_Delete( void* drawingObject ) ;
	void _lucDrawingObject_Print( void* drawingObject, Stream* stream ) ;
	void* _lucDrawingObject_Copy( void* drawingObject, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) ;

void _lucDrawingObject_AssignFromXML( void* drawingObject, Stg_ComponentFactory* cf, void* data ) ;
	void _lucDrawingObject_Build( void* camera, void* data );
	void _lucDrawingObject_Initialise( void* camera, void* data );
	void _lucDrawingObject_Execute( void* camera, void* data );
	void _lucDrawingObject_Destroy( void* camera, void* data );

	/* +++ Public Functions +++ */
	void lucDrawingObject_Setup( void* drawingObject, void* context ) ;
	void lucDrawingObject_Draw( void* drawingObject, lucWindow* window, lucViewportInfo* viewportInfo, void* context ) ;
	void lucDrawingObject_CleanUp( void* drawingObject, void* context ) ;

	typedef enum { GreaterThan, LessThan, EqualTo } lucDrawingObjectMask_Type;
	typedef struct { 
		lucDrawingObjectMask_Type type;
		double                    value;
		double                    tolerance;
	} lucDrawingObjectMask;

	void lucDrawingObjectMask_Construct( lucDrawingObjectMask* self, Name drawingObjectName, Stg_ComponentFactory* cf, void* data ) ;
	Bool lucDrawingObjectMask_Test( lucDrawingObjectMask* self, double value ) ;
#endif
