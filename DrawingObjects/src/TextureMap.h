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
** $Id: TextureMap.h 628 2006-10-12 08:23:07Z SteveQuenette $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#ifndef __lucTextureMap_h__
#define __lucTextureMap_h__

	/** Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
	extern const Type lucTextureMap_Type;
		
	/** Class contents - this is defined as a macro so that sub-classes of this class can use this macro at the start of the definition of their struct */
	#define __lucTextureMap \
		/* Macro defining parent goes here - This means you can cast this class as its parent */ \
		__lucOpenGLDrawingObject \
		/* Virtual functions go here */ \
		/* Other info */\
		Bool                                               iAmMaster;        \
		lucPixel*                                          pixelData;        \
		Pixel_Index                                        imageWidth;       \
		Pixel_Index                                        imageHeight;      \
		Coord                                              bottomLeftCoord;  \
		Coord                                              bottomRightCoord; \
		Coord                                              topRightCoord;    \
		Coord                                              topLeftCoord;     \
		lucInputFormat*                                    inputFormat;      \
		
	struct lucTextureMap { __lucTextureMap };
	
	/** Private Constructor: This will accept all the virtual functions for this class as arguments. */
	
	#ifndef ZERO
	#define ZERO 0
	#endif

	#define LUCTEXTUREMAP_DEFARGS \
                LUCOPENGLDRAWINGOBJECT_DEFARGS

	#define LUCTEXTUREMAP_PASSARGS \
                LUCOPENGLDRAWINGOBJECT_PASSARGS

	lucTextureMap* _lucTextureMap_New(  LUCTEXTUREMAP_DEFARGS  );

	void _lucTextureMap_Delete( void* drawingObject ) ;
	void _lucTextureMap_Print( void* drawingObject, Stream* stream ) ;
	void* _lucTextureMap_Copy( void* drawingObject, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap) ;

	/* 'Stg_Component' implementations */
	void* _lucTextureMap_DefaultNew( Name name ) ;
	void _lucTextureMap_AssignFromXML( void* drawingObject, Stg_ComponentFactory* cf, void* data );
	void _lucTextureMap_Build( void* drawingObject, void* data ) ;
	void _lucTextureMap_Initialise( void* drawingObject, void* data ) ;
	void _lucTextureMap_Execute( void* drawingObject, void* data );
	void _lucTextureMap_Destroy( void* drawingObject, void* data ) ;
	
	void _lucTextureMap_Setup( void* drawingObject, void* _context ) ;
	void _lucTextureMap_Draw( void* drawingObject, lucWindow* window, lucViewportInfo* viewportInfo, void* _context ) ;
	void _lucTextureMap_CleanUp( void* drawingObject, void* _context ) ;

	void _lucTextureMap_BuildDisplayList( void* drawingObject, void* _context ) ;

//	void _lucTextureMap_ReadImage( lucTextureMap* self, Name imageName, PixelIndex* width, PixelIndex* eight ) ;
#endif

