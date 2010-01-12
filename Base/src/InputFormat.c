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

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>

#include "types.h"
#include "InputFormat.h"
#include "ColourMap.h"
#include "Window.h"
#include "DrawingObject_Register.h"
#include "DrawingObject.h"
#include "Camera.h"
#include "Init.h"
#include "Window.h"

#include <assert.h>

#ifndef MASTER
	#define MASTER 0
#endif

const Type lucInputFormat_Type = "lucInputFormat";

lucInputFormat* _lucInputFormat_New(  LUCINPUTFORMAT_DEFARGS  )
{
	lucInputFormat*    self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( _sizeOfSelf >= sizeof(lucInputFormat) );
	/* The following terms are parameters that have been passed into this function but are being set before being passed onto the parent */
	/* This means that any values of these parameters that are passed into this function are not passed onto the parent function
	   and so should be set to ZERO in any children of this class. */
	nameAllocationType = NON_GLOBAL;

	self = (lucInputFormat*) _Stg_Component_New(  STG_COMPONENT_PASSARGS  );

	self->_input = _input;

	return self;
}


void _lucInputFormat_Init( 
		lucInputFormat*                                   self, 
		Name                                               imageName)
{
	self->imageName     = StG_Strdup( imageName );
}

void lucInputFormat_InitAll( 
		void*                                              inputFormat,
		Name                                               imageName )
{
	lucInputFormat* self        = inputFormat;

	_lucInputFormat_Init( self, imageName );
}

	
void _lucInputFormat_Delete( void* inputFormat ) {
	lucInputFormat* self        = inputFormat;

	Memory_Free( self->imageName );
	
	_Stg_Component_Delete( self );
}

void _lucInputFormat_Print( void* inputFormat, Stream* stream ) {
	lucInputFormat*          self        = inputFormat;
	
	Journal_Printf( stream, "lucInputFormat: %s\n", self->name );

	/* Print Parent */
	_Stg_Component_Print( self, stream );

	Journal_PrintString( stream, self->imageName );
}

void* _lucInputFormat_Copy( void* inputFormat, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	lucInputFormat* self        = inputFormat;
	lucInputFormat* newInputFormat;

	newInputFormat = (lucInputFormat*)_Stg_Component_Copy( self, dest, deep, nameExt, ptrMap );

	return (void*) newInputFormat;
}

void _lucInputFormat_AssignFromXML( void* inputFormat, Stg_ComponentFactory* cf, void* data ) {
	lucInputFormat*        self               = (lucInputFormat*) inputFormat; 

	self->context = Stg_ComponentFactory_ConstructByKey( cf, self->name, "Context", AbstractContext, False, data );
	if( !self->context ) 
		self->context = Stg_ComponentFactory_ConstructByName( cf, "context", AbstractContext, True, data );
}
void _lucInputFormat_Build( void* inputFormat, void* data ) { }
void _lucInputFormat_Initialise( void* inputFormat, void* data ) { }
void _lucInputFormat_Execute( void* inputFormat, void* data ) { }
void _lucInputFormat_Destroy( void* inputFormat, void* data ) { }


lucPixel* lucInputFormat_Input( void* inputFormat, Name imageName, Pixel_Index *imageWidth, Pixel_Index* imageHeight  ) {
	lucInputFormat*        self               = (lucInputFormat*) inputFormat; 
	lucPixel* pixelData = self->_input(self, imageName, imageWidth, imageHeight );

	return pixelData;
}




