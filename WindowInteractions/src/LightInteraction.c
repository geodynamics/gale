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
** $Id: LightInteraction.c 628 2006-10-12 08:23:07Z SteveQuenette $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>

#include <glucifer/Base/Base.h>
#include <glucifer/Windowing/Windowing.h>
#include <glucifer/RenderingEngines/RenderingEngines.h>
#include <glucifer/DrawingObjects/DrawingObjects.h>

#include "types.h"
#include "LightInteraction.h"

#include <assert.h>
#include <gl.h>
#include <glu.h>

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type lucLightInteraction_Type = "lucLightInteraction";

/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
lucLightInteraction* _lucLightInteraction_New(  LUCLIGHTINTERACTION_DEFARGS  ) 
{
	lucLightInteraction*					self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( _sizeOfSelf >= sizeof(lucLightInteraction) );
	self = (lucLightInteraction*) _lucWindowInteraction_New(  LUCWINDOWINTERACTION_PASSARGS  );
	
	return self;
}

void _lucLightInteraction_Init( 
		lucLightInteraction*                                      self  ) 
{
}

void _lucLightInteraction_Delete( void* LightInteraction ) {
	lucLightInteraction*  self = (lucLightInteraction*)LightInteraction;

	_lucWindowInteraction_Delete( self );
}

void _lucLightInteraction_Print( void* LightInteraction, Stream* stream ) {
	lucLightInteraction*  self = (lucLightInteraction*)LightInteraction;

	_lucWindowInteraction_Print( self, stream );
}

void* _lucLightInteraction_Copy( void* LightInteraction, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap) {
	lucLightInteraction*  self = (lucLightInteraction*) LightInteraction;
	lucLightInteraction* newLightInteraction;

	newLightInteraction = _lucWindowInteraction_Copy( self, dest, deep, nameExt, ptrMap );

	/* TODO */
	abort();

	return (void*) newLightInteraction;
}


void* _lucLightInteraction_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                                  _sizeOfSelf = sizeof(lucLightInteraction);
	Type                                                          type = lucLightInteraction_Type;
	Stg_Class_DeleteFunction*                                  _delete = _lucLightInteraction_Delete;
	Stg_Class_PrintFunction*                                    _print = _lucLightInteraction_Print;
	Stg_Class_CopyFunction*                                      _copy = NULL;
	Stg_Component_DefaultConstructorFunction*      _defaultConstructor = _lucLightInteraction_DefaultNew;
	Stg_Component_ConstructFunction*                        _construct = _lucLightInteraction_AssignFromXML;
	Stg_Component_BuildFunction*                                _build = _lucLightInteraction_Build;
	Stg_Component_InitialiseFunction*                      _initialise = _lucLightInteraction_Initialise;
	Stg_Component_ExecuteFunction*                            _execute = _lucLightInteraction_Execute;
	Stg_Component_DestroyFunction*                            _destroy = _lucLightInteraction_Destroy;
	lucWindowInteraction_MouseMotionFunction*             _mouseMotion = _lucLightInteraction_MouseMotion;
	lucWindowInteraction_MouseClickFunction*               _mouseClick = _lucLightInteraction_MouseClick;
	lucWindowInteraction_MouseMessageFunction*           _mouseMessage = _lucLightInteraction_MouseMessage;
	lucWindowInteraction_KeyboardEventFunction*         _keyboardEvent = _lucLightInteraction_KeyboardEvent;
	lucWindowInteraction_KeyboardMessageFunction*     _keyboardMessage = _lucLightInteraction_KeyboardMessage;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

	return (void*) _lucLightInteraction_New(  LUCLIGHTINTERACTION_PASSARGS  );
}

void _lucLightInteraction_AssignFromXML( void* LightInteraction, Stg_ComponentFactory* cf, void* data ){
	lucLightInteraction*  self = LightInteraction;

	/* Construct Parent */
	_lucWindowInteraction_AssignFromXML(self, cf, data );
	
	_lucLightInteraction_Init( self );
}

void _lucLightInteraction_Build( void* renderingEngine, void* data ) {}
void _lucLightInteraction_Initialise( void* renderingEngine, void* data ) {}
void _lucLightInteraction_Execute( void* renderingEngine, void* data ) {}
void _lucLightInteraction_Destroy( void* renderingEngine, void* data ) {}

/* This component doesn't use any mouse interaction */
void _lucLightInteraction_MouseMotion( void* windowInteraction, lucWindow* window, lucMouseButton button, Pixel_Index xpos, Pixel_Index ypos, Pixel_Index startx, Pixel_Index starty) {}
void _lucLightInteraction_MouseClick( void* windowInteraction, lucWindow* window, lucMouseButton button, lucMouseState state, Pixel_Index xpos, Pixel_Index ypos) { } 
void _lucLightInteraction_MouseMessage( void* windowInteraction, Stream* stream ) { }


void _lucLightInteraction_KeyboardEvent( void* WindowInteraction, lucWindow* window, char key, Pixel_Index xpos, Pixel_Index ypos) {
	lucLightInteraction*   self;
	lucViewportInfo*            viewportInfo;
	lucViewport*                viewport;
	Coord                       coord;
	lucLight*                   light;
	Light_Index                 light_I;
	Light_Index                 lightCount;
	float                       posX = 0;
	float                       posY = 0;
	float                       posZ = 0;
	float                       initialPosition[4];
	int                         i;
	Light_Index                  currentLight_I=0;

	self = (lucLightInteraction*) WindowInteraction;

	/* This function works when the key pressed is one of  'x', 'y, 'z', 'l', 'm', 'n'*/
	if ( ( key != 'x' )&&( key != 'y' ) && ( key != 'z' ) && ( key != 'k' ) && ( key != 'l' ) &&( key != 'm') &&( key != 'p')&&( key !='w') )
		return;
	
	/* Find which viewport the user clicked on */
	viewportInfo = lucWindow_GetViewportInfoByPixel( window, xpos, ypos );
	if ( viewportInfo == NULL )
		return; /* If the user hasn't clicked on a viewport at all, then return */

	viewport = viewportInfo->viewport;

	
	
	/* Get the light */
	/* Loop through lights that are registered on the window */
	lightCount   = lucLight_Register_GetCount( viewport->light_Register );
	
	for ( light_I = 0 ; light_I < lightCount ; light_I++ ) {
		light = lucLight_Register_GetByIndex( viewport->light_Register, light_I );
		//printf(" Initial position SELF is %.2f, %.2f, %.2f, %.2f \n", light->position[0], light->position[1], light->position[2], light->position[3]);

		for(i =0; i<4; i++) initialPosition[i] = 0;
		
	
	       // lucLight_GetGLIndex(light->index, glIndex);
		//glGetLightfv(GL_LIGHT0, GL_POSITION, initialPosition);
	
		//printf(" Initial position is %.2f, %.2f, %.2f, %.2f \n", initialPosition[0], initialPosition[1], initialPosition[2], initialPosition[3]);

	}

	
	/* Implements the light position changes */

	if(key == 'w'){
		/* retrieve the currentLightIndex and increases it if suitable */
		lucLight_Register_ChangeCurrentLightIndex(viewport->light_Register);
		currentLight_I = lucLight_Register_GetCurrentLightIndex( viewport->light_Register );

		printf("\n\n\n Selected light is index %d \n\n\n\n\n", currentLight_I);
	}

	if ( key == 'x' ){
	       // printf(" Increasing light x position by 0.5 \n ");
		posX = 0.5;
	}
	

	if (key == 'y' ){
	       // printf(" Increasing light y position by 0.5 \n ");
		posY = 0.5;
	}

	if (key == 'z' ){
	        //printf(" Increasing light z position by 0.5 \n ");
		posZ = 0.5;
	}
	
	if ( key == 'k' ){
	        //printf(" Decreasing light x position by 0.5 \n ");
		posX = -0.5;
	}
	

	if (key == 'l' ){
	        //printf(" Decreasing light y position by 0.5 \n ");
		posY = -0.5;
	}

	if (key == 'm' ){
	        //printf(" Decreasing light z position by 0.5 \n ");
		posZ = -0.5;
	}
	
	if (key == 'p' ){
		/* prints the light position */
		glGetLightfv(GL_LIGHT0 + currentLight_I, GL_POSITION, initialPosition);
		printf(" position by glGET is %.2f, %.2f, %.2f, %.2f \n", initialPosition[0], initialPosition[1], initialPosition[2], initialPosition[3]);
	}
	/* Retrieves the light corresponding to the currentLightIndex */
	currentLight_I = lucLight_Register_GetCurrentLightIndex( viewport->light_Register );
	light = lucLight_Register_GetByIndex( viewport->light_Register, currentLight_I );
	lucLight_Position(light, currentLight_I, posX, posY, posZ, 0);

	/*float position[4];
	glGetLightfv(GL_LIGHT0, GL_POSITION, position); */
	
	printf(" Position for light index %d is %.2f, %.2f, %.2f, %.2f \n", currentLight_I, light->position[0],  light->position[1],  light->position[2],  light->position[3]);
	printf(" SpotCutOff is %.2f \n", light->spotCutOff);

	/* Get spatial coordinate that the user clicked on */
	lucViewportInfo_GetCoordFromPixel( viewportInfo, xpos, ypos, coord );
}

void _lucLightInteraction_KeyboardMessage( void* windowInteraction, Stream* stream ) {
	Journal_Printf( stream,
			"x:                            Increases the X position of the light by 0.5.\n" );
	Journal_Printf( stream,
			"y:                            Increases the Y position of the light by 0.5.\n" );	
	Journal_Printf( stream,
			"z:                            Increases the Z position of the light by 0.5.\n" );
	Journal_Printf( stream,
			"l:                            Decreases the X position of the light by 0.5.\n" );
	Journal_Printf( stream,
			"m:                            Decreases the Y position of the light by 0.5.\n" );
	Journal_Printf( stream,
			"n:                            Decreases the Z position of the light by 0.5.\n" );
}




