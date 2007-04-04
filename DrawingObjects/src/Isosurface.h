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
** $Id: Isosurface.h 628 2006-10-12 08:23:07Z SteveQuenette $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#ifndef __lucIsosurface_h__
#define __lucIsosurface_h__

	/** Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
	extern const Type lucIsosurface_Type;

	typedef struct {
		double pos1[3], pos2[3], pos3[3];
		double normal1[3], normal2[3], normal3[3];
	} Surface_Triangle;
	
	typedef struct {
		double value;
		double pos[3];
		double normal[3];
	} Vertex;
		
	/** Class contents - this is defined as a macro so that sub-classes of this class can use this macro at the start of the definition of their struct */
	#define __lucIsosurface \
		/* Macro defining parent goes here - This means you can cast this class as its parent */ \
		__lucOpenGLDrawingObject \
		/* Virtual functions go here */ \
		/* Other info */\
		FieldVariable*                                     isosurfaceField;        \
		double                                             isovalue;               \
		IJK                                                resolution;             \
		Bool                                               drawWalls;              \
		Bool                                               wireframe;              \
		Bool                                               cullFrontFace;          \
		Bool                                               cullBackFace;           \
		/* Colour Parameters */ \
		lucColour                                          colour;                 \
		lucColourMap*                                      colourMap;              \
		FieldVariable*                                     colourField;            \
		/* Masking parameters */\
		FieldVariable*                                     maskField;              \
		lucDrawingObjectMask                               mask;                   \
		/* Calculated Values */ \
		Surface_Triangle*                                  triangleList;           \
		Index                                              triangleCount;          \
		Index                                              trianglesAlloced;       \

	struct lucIsosurface { __lucIsosurface };
	
	/** Private Constructor: This will accept all the virtual functions for this class as arguments. */
	lucIsosurface* _lucIsosurface_New( 
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

	void _lucIsosurface_Delete( void* drawingObject ) ;
	void _lucIsosurface_Print( void* drawingObject, Stream* stream ) ;
	void* _lucIsosurface_Copy( void* drawingObject, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap) ;

	/* 'Stg_Component' implementations */
	void* _lucIsosurface_DefaultNew( Name name ) ;
	void _lucIsosurface_Construct( void* drawingObject, Stg_ComponentFactory* cf, void* data );
	void _lucIsosurface_Build( void* drawingObject, void* data ) ;
	void _lucIsosurface_Initialise( void* drawingObject, void* data ) ;
	void _lucIsosurface_Execute( void* drawingObject, void* data );
	void _lucIsosurface_Destroy( void* drawingObject, void* data ) ;
	
	/* Drawing Object Implementations */
	void _lucIsosurface_Setup( void* drawingObject, void* _context ) ;
	void _lucIsosurface_Draw( void* drawingObject, lucWindow* window, lucViewportInfo* viewportInfo, void* _context ) ;
	void _lucIsosurface_CleanUp( void* drawingObject, void* _context ) ;

	/* OpenGL Drawing Object Implementation */
	void _lucIsosurface_BuildDisplayList( void* drawingObject, void* _context ) ;

	Bool lucIsosurface_TestMask( lucIsosurface* self, Coord pos ) ;

	void lucIsosurface_MarchingCubes( lucIsosurface* self, Vertex*** vertex ) ;
	void lucIsosurface_Normals( lucIsosurface* self, Vertex*** vertex ) ;
	void lucIsosurface_DrawWalls( lucIsosurface* self, Vertex*** array ) ;

	void lucIsosurface_CreateIntermediatePoints( lucIsosurface* self, Vertex **points, Dimension_Index axis ) ;
	void lucIsosurface_MarchingRectangles( lucIsosurface* self, Vertex** points, char cubeType, char order) ;
	void lucIsosurface_WallElement( lucIsosurface* self, Vertex** points, char order, Dimension_Index axis ) ;
	void lucIsosurface_AddWallTriangle( lucIsosurface* self, int a, int b, int c, Vertex** points, char order ) ;

	void lucIsosurface_SetupPointsX( Vertex** points, Vertex*** array, Index i, Index j, Index k );
	void lucIsosurface_SetupPointsY( Vertex** points, Vertex*** array, Index i, Index j, Index k );
	void lucIsosurface_SetupPointsZ( Vertex** points, Vertex*** array, Index i, Index j, Index k );
			
#endif
