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
** $Id: TimeStep.c 510 2006-02-17 04:33:32Z RobertTurnbull $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>

#include <glucifer/Base/Base.h>
#include <glucifer/RenderingEngines/RenderingEngines.h>

#ifdef HAVE_GL2PS
	#include <gl2ps.h>
#endif

#include "types.h"
#include "OpenGLDrawingObject.h"
#include "TimeStep.h"

#include <assert.h>
#include <gl.h>
#include <glu.h>
#include <string.h>

#ifndef MASTER
	#define MASTER 0
#endif



const Type lucTimeStep_Type = "lucTimeStep";

lucTimeStep* lucTimeStep_New( 
		Name                                               name,
		lucColour                                          colour,
		Bool                                               frame,
		Bool                                               time)
{
	lucTimeStep* self = (lucTimeStep*) _lucTimeStep_DefaultNew( name );

	lucTimeStep_InitAll( self, colour, frame, time);

	return self;
}

lucTimeStep* _lucTimeStep_New(  LUCTIMESTEP_DEFARGS  )
{
	lucTimeStep*    self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( _sizeOfSelf >= sizeof(lucTimeStep) );
	self = (lucTimeStep*)  _lucOpenGLDrawingObject_New(  LUCOPENGLDRAWINGOBJECT_PASSARGS  );
	
	
	return self;
}

void lucTimeStep_Init(		
		lucTimeStep*                                         self,
		lucColour                                         colour,
		Bool                                              frame,
		Bool                                              time)
{
	memcpy( &(self->colour), &colour, sizeof(lucColour) );	
	self->frame = frame;
	self->time = time;
}

void lucTimeStep_InitAll( 
		void*                                              timeStep,
		lucColour                                          colour,
		Bool                                               frame,
		Bool                                               time )
{
	lucTimeStep* self        = timeStep;

	/* TODO Init parent */
	lucTimeStep_Init( self, colour, frame, time );
}

void _lucTimeStep_Delete( void* drawingObject ) {
	lucTimeStep*  self = (lucTimeStep*)drawingObject;

	_lucOpenGLDrawingObject_Delete( self );
}

void _lucTimeStep_Print( void* drawingObject, Stream* stream ) {
	lucTimeStep*  self = (lucTimeStep*)drawingObject;

	_lucOpenGLDrawingObject_Print( self, stream );
}

void* _lucTimeStep_Copy( void* timeStep, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	lucTimeStep* self        = timeStep;
	lucTimeStep* newTimeStep;

	newTimeStep = _Stg_Component_Copy( self, dest, deep, nameExt, ptrMap );

	memcpy( &(newTimeStep->colour),       &(self->colour),       sizeof(lucColour) );
	newTimeStep->frame = self->frame;
	newTimeStep->time = self->time;
	
	return (void*) newTimeStep;
}

void* _lucTimeStep_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                                     _sizeOfSelf = sizeof( lucTimeStep );
	Type                                                             type = lucTimeStep_Type;
	Stg_Class_DeleteFunction*                                     _delete = _lucTimeStep_Delete;
	Stg_Class_PrintFunction*                                       _print = _lucTimeStep_Print;
	Stg_Class_CopyFunction*                                         _copy = _lucTimeStep_Copy;
	Stg_Component_DefaultConstructorFunction*         _defaultConstructor = _lucTimeStep_DefaultNew;
	Stg_Component_ConstructFunction*                           _construct = _lucTimeStep_AssignFromXML;
	Stg_Component_BuildFunction*                                   _build = _lucTimeStep_Build;
	Stg_Component_InitialiseFunction*                         _initialise = _lucTimeStep_Initialise;
	Stg_Component_ExecuteFunction*                               _execute = _lucTimeStep_Execute;
	Stg_Component_DestroyFunction*                               _destroy = _lucTimeStep_Destroy;
	lucDrawingObject_SetupFunction*                                _setup = _lucTimeStep_Setup;
	lucDrawingObject_DrawFunction*                                  _draw = _lucTimeStep_Draw;
	lucDrawingObject_CleanUpFunction*                            _cleanUp = _lucTimeStep_CleanUp;
	lucOpenGLDrawingObject_BuildDisplayListFunction*    _buildDisplayList = _lucTimeStep_BuildDisplayList;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = ZERO;

	return _lucTimeStep_New(  LUCTIMESTEP_PASSARGS  );
}

void _lucTimeStep_AssignFromXML( void* timeStep, Stg_ComponentFactory* cf, void* data ) {
	lucTimeStep*             self               = (lucTimeStep*) timeStep;
        Name colourName;
	Bool frame;
	Bool currentTime;
	
	colourName  = Stg_ComponentFactory_GetString( cf, self->name, "colour", "Black") ;
	frame  = Stg_ComponentFactory_GetBool( cf, self->name, "frame", True) ;
	currentTime = Stg_ComponentFactory_GetBool( cf, self->name, "time", False) ;

	lucColour_FromString( &self->colour, colourName );
	
	lucTimeStep_InitAll( self,
		        self->colour, frame, currentTime);
			
}

void _lucTimeStep_Build( void* TimeStep, void* data ) { }
void _lucTimeStep_Initialise( void* TimeStep, void* data ) { }
void _lucTimeStep_Execute( void* TimeStep, void* data ) { }
void _lucTimeStep_Destroy( void* TimeStep, void* data ) { }

void _lucTimeStep_Setup( void* drawingObject, void* _context ) {
	lucTimeStep*       self            = (lucTimeStep*)drawingObject;
	_lucOpenGLDrawingObject_Setup( self, _context );
}
void _lucTimeStep_Draw( void* drawingObject, lucWindow* window, lucViewportInfo* viewportInfo, void* _context ) {
	lucTimeStep*      self            = (lucTimeStep*)drawingObject;
       	AbstractContext*    context = (AbstractContext*)_context;
	
	int stringWidth = 0;
	int time = context->timeStep;
	double currentTime = context->currentTime;
	/* Allocating Memory */


  /* Luke's weird stuff - as per revision 74
	char* displayString = Memory_Alloc_Array( char, 100, "displayString");
	char* timeStepString = Memory_Alloc_Array( char, 10, "timeStepString");
 	char* currentTimeString = Memory_Alloc_Array( char, 20, "currentTimeString");
  sprintf( displayString, "Shortened distance (cm): %e", currentTime * 6.9444444444444e-4 );
  */
	
  char* displayString = Memory_Alloc_Array( char, 40, "displayString");
	char* timeStepString = Memory_Alloc_Array( char, 10, "timeStepString");
 	char* currentTimeString = Memory_Alloc_Array( char, 20, "currentTimeString");

	sprintf(timeStepString, "%d", time );	
	sprintf(currentTimeString, "%e", currentTime );
	
	if(self->frame){
		sprintf(displayString, "%s", "Frame ");
		strcat(displayString, timeStepString);
	    	if(self->time){
			strcat(displayString, "   Time ");
			strcat(displayString, currentTimeString);
		}
		strcat(displayString, "\0");
	}
	else{	
		sprintf(displayString, "%s", "Time ");
		strcat(displayString, currentTimeString);
		strcat(displayString, "\0");
	}

	
	/* Set up 2D Viewer the size of the viewport */
	lucViewport2d(True, viewportInfo);

	/* Set the colour so that it'll show up against the background */
	lucColour_SetComplimentaryOpenGLColour( &window->backgroundColour );

	/* Print TimeStep */	       
        glColor4f(
		self->colour.red,
		self->colour.green,
		self->colour.blue,
		self->colour.opacity );

	stringWidth = lucStringWidth(displayString );
	lucMoveRaster( - stringWidth/2, -20 );
	lucPrintString(displayString);

	/* Free the memory */
	Memory_Free(timeStepString) ;
	Memory_Free(currentTimeString);
	Memory_Free(displayString);

	/* Restore the viewport */
	lucViewport2d(False, viewportInfo);

}

void _lucTimeStep_CleanUp( void* drawingObject, void* _context ) {
	lucTimeStep*      self            = (lucTimeStep*)drawingObject;
	
	_lucOpenGLDrawingObject_CleanUp( self, _context );

}

void _lucTimeStep_BuildDisplayList( void* drawingObject, void* _context ) {
	}







