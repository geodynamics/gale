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
** $Id: TextureMap.c 791 2008-09-01 02:09:06Z JulianGiordani $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>

#include <glucifer/Base/Base.h>
#include <glucifer/RenderingEngines/RenderingEngines.h>

#include "types.h"
#include "OpenGLDrawingObject.h"
#include "TextureMap.h"

#include <assert.h>
#include <gl.h>
#include <glu.h>
#include <string.h>

#ifndef MASTER
	#define MASTER 0
#endif


/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type lucTextureMap_Type = "lucTextureMap";

/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
lucTextureMap* _lucTextureMap_New(  LUCTEXTUREMAP_DEFARGS  ) 
{
	lucTextureMap*					self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( _sizeOfSelf >= sizeof(lucTextureMap) );
	self = (lucTextureMap*) _lucOpenGLDrawingObject_New(  LUCOPENGLDRAWINGOBJECT_PASSARGS  );
				
	return self;
}

void _lucTextureMap_Init( 
	lucTextureMap*                                               self,
	Bool                                                         iAmMaster,
	Name                                                         imageName,
	double                                                       bottomLeftX,
	double                                                       bottomLeftY,
	double                                                       bottomLeftZ,
	double                                                       bottomRightX,
	double                                                       bottomRightY,
	double                                                       bottomRightZ,
	double                                                       topRightX,
	double                                                       topRightY,
	double                                                       topRightZ,
	double                                                       topLeftX,
	double                                                       topLeftY,
	double                                                       topLeftZ )
{

 	self->iAmMaster = iAmMaster;

//	if ( ! iAmMaster )
//		return;
	self->bottomLeftCoord[ I_AXIS ] = bottomLeftX;
	self->bottomLeftCoord[ J_AXIS ] = bottomLeftY;
	self->bottomLeftCoord[ K_AXIS ] = bottomLeftZ;
	
	self->bottomRightCoord[ I_AXIS ] = bottomRightX;
	self->bottomRightCoord[ J_AXIS ] = bottomRightY;
	self->bottomRightCoord[ K_AXIS ] = bottomRightZ;
	
	self->topRightCoord[ I_AXIS ] = topRightX;
	self->topRightCoord[ J_AXIS ] = topRightY;
	self->topRightCoord[ K_AXIS ] = topRightZ;
	
	self->topLeftCoord[ I_AXIS ] = topLeftX;
	self->topLeftCoord[ J_AXIS ] = topLeftY;
	self->topLeftCoord[ K_AXIS ] = topLeftZ;

	/* Open and read Image */	
	self->inputFormat = lucInputFormat_Register_CreateFromFileName( lucInputFormat_Register_Singleton, imageName );		
	self->pixelData = lucInputFormat_Input( self->inputFormat, imageName, &self->imageWidth, &self->imageHeight );
}

void _lucTextureMap_Delete( void* drawingObject ) {
	lucTextureMap*  self = (lucTextureMap*)drawingObject;

	if ( self->iAmMaster ) 
		Memory_Free( self->pixelData );

	_lucOpenGLDrawingObject_Delete( self );
}

void _lucTextureMap_Print( void* drawingObject, Stream* stream ) {
	lucTextureMap*  self = (lucTextureMap*)drawingObject;

	_lucOpenGLDrawingObject_Print( self, stream );
	
	Journal_PrintValue( stream, self->iAmMaster );

	Journal_PrintPointer( stream, self->pixelData );

	Journal_PrintValue( stream, self->imageWidth );
	Journal_PrintValue( stream, self->imageHeight );
	
	Journal_PrintArray( stream, self->bottomLeftCoord,  3 );
	Journal_PrintArray( stream, self->bottomRightCoord, 3 );
	Journal_PrintArray( stream, self->topRightCoord,    3 );
	Journal_PrintArray( stream, self->topLeftCoord,     3 );
}

void* _lucTextureMap_Copy( void* drawingObject, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap) {
	lucTextureMap*  self = (lucTextureMap*)drawingObject;
	lucTextureMap* newDrawingObject;

	newDrawingObject = _lucOpenGLDrawingObject_Copy( self, dest, deep, nameExt, ptrMap );
	
	//self->pixelData = _lucInputFormat_Copy(self, ...);

	/* TODO */
	abort();

	return (void*) newDrawingObject;
}


void* _lucTextureMap_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                                     _sizeOfSelf = sizeof(lucTextureMap);
	Type                                                             type = lucTextureMap_Type;
	Stg_Class_DeleteFunction*                                     _delete = _lucTextureMap_Delete;
	Stg_Class_PrintFunction*                                       _print = _lucTextureMap_Print;
	Stg_Class_CopyFunction*                                         _copy = NULL;
	Stg_Component_DefaultConstructorFunction*         _defaultConstructor = _lucTextureMap_DefaultNew;
	Stg_Component_ConstructFunction*                           _construct = _lucTextureMap_AssignFromXML;
	Stg_Component_BuildFunction*                                   _build = _lucTextureMap_Build;
	Stg_Component_InitialiseFunction*                         _initialise = _lucTextureMap_Initialise;
	Stg_Component_ExecuteFunction*                               _execute = _lucTextureMap_Execute;
	Stg_Component_DestroyFunction*                               _destroy = _lucTextureMap_Destroy;
	lucDrawingObject_SetupFunction*                                _setup = _lucTextureMap_Setup;
	lucDrawingObject_DrawFunction*                                  _draw = _lucTextureMap_Draw;
	lucDrawingObject_CleanUpFunction*                            _cleanUp = _lucTextureMap_CleanUp;
	lucOpenGLDrawingObject_BuildDisplayListFunction*    _buildDisplayList = _lucTextureMap_BuildDisplayList;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

	return (void*) _lucTextureMap_New(  LUCTEXTUREMAP_PASSARGS  );
		
		
}

void _lucTextureMap_AssignFromXML( void* drawingObject, Stg_ComponentFactory* cf, void* data ){
	lucTextureMap*  self = (lucTextureMap*)drawingObject;

	/* Construct Parent */
	_lucOpenGLDrawingObject_AssignFromXML( self, cf, data );

	_lucTextureMap_Init( 
			self, 
			Stg_ComponentFactory_GetRootDictUnsignedInt( cf, "rank", (unsigned)-1 ) == MASTER,
			Stg_ComponentFactory_GetString( cf, self->name, "image",        "" ),
			Stg_ComponentFactory_GetDouble( cf, self->name, "bottomLeftX",  0.0 ),
			Stg_ComponentFactory_GetDouble( cf, self->name, "bottomLeftY",  0.0 ),
			Stg_ComponentFactory_GetDouble( cf, self->name, "bottomLeftZ",  0.0 ),
			Stg_ComponentFactory_GetDouble( cf, self->name, "bottomRightX", 0.0 ),
			Stg_ComponentFactory_GetDouble( cf, self->name, "bottomRightY", 0.0 ),
			Stg_ComponentFactory_GetDouble( cf, self->name, "bottomRightZ", 0.0 ),
			Stg_ComponentFactory_GetDouble( cf, self->name, "topRightX",    0.0 ),
			Stg_ComponentFactory_GetDouble( cf, self->name, "topRightY",    0.0 ),
			Stg_ComponentFactory_GetDouble( cf, self->name, "topRightZ",    0.0 ),
			Stg_ComponentFactory_GetDouble( cf, self->name, "topLeftX",     0.0 ),
			Stg_ComponentFactory_GetDouble( cf, self->name, "topLeftY",     0.0 ),
			Stg_ComponentFactory_GetDouble( cf, self->name, "topLeftZ",     0.0 ) );
}

void _lucTextureMap_Build( void* drawingObject, void* data ) {}
void _lucTextureMap_Initialise( void* drawingObject, void* data ) {}
void _lucTextureMap_Execute( void* drawingObject, void* data ) {}
void _lucTextureMap_Destroy( void* drawingObject, void* data ) {}


void _lucTextureMap_Setup( void* drawingObject, void* _context ) {
	lucTextureMap*       self            = (lucTextureMap*)drawingObject;
	_lucOpenGLDrawingObject_Setup( self, _context );
}
	
void _lucTextureMap_Draw( void* drawingObject, lucWindow* window, lucViewportInfo* viewportInfo, void* _context ) {
	lucTextureMap*       self            = (lucTextureMap*)drawingObject;
	_lucOpenGLDrawingObject_Draw( self, window, viewportInfo, _context );
}


void _lucTextureMap_CleanUp( void* drawingObject, void* _context ) {
	lucTextureMap*       self            = (lucTextureMap*)drawingObject;
	_lucOpenGLDrawingObject_CleanUp( self, _context );
}
	
void _lucTextureMap_BuildDisplayList( void* drawingObject, void* _context ) {
	lucTextureMap*          self          = (lucTextureMap*)   drawingObject;

	//if ( !self->iAmMaster )
		//return;
			
	/* Setup stuff to draw */
	
	/* VERY Important, otherwise blending happens between pixels already stored in buffer and the texture Map does not appear !! */
	glDisable(GL_BLEND);
		
	glPixelStorei( GL_UNPACK_ALIGNMENT, 1 );
	
	glTexImage2D( GL_TEXTURE_2D, 0, 3, self->imageWidth, self->imageHeight, 0, GL_RGB, GL_UNSIGNED_BYTE, self->pixelData );
	
	glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP );
	glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP );
	
	glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST );
	glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST );
	
	glEnable( GL_TEXTURE_2D );
	glTexEnvf( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL );

	glShadeModel( GL_FLAT );

	/* Draw Texture Map */
	glBegin( GL_QUADS );
		glTexCoord2f( 0.0, 0.0 ); glVertex3dv( self->bottomLeftCoord );
		glTexCoord2f( 0.0, 1.0 ); glVertex3dv( self->topLeftCoord );
		glTexCoord2f( 1.0, 1.0 ); glVertex3dv( self->topRightCoord );
		glTexCoord2f( 1.0, 0.0 ); glVertex3dv( self->bottomRightCoord );
	glEnd();
	glFlush();

	glDisable(GL_TEXTURE_2D);
	glEnable(GL_BLEND);
}





