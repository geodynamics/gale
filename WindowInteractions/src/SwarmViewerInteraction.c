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
** $Id: SwarmViewerInteraction.c 510 2006-02-17 04:33:32Z RobertTurnbull $
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
#include "SwarmViewerInteraction.h"

#include <assert.h>

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type lucSwarmViewerInteraction_Type = "lucSwarmViewerInteraction";


/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
lucSwarmViewerInteraction* _lucSwarmViewerInteraction_New(  LUCSWARMVIEWERINTERACTION_DEFARGS  ) 
{
	lucSwarmViewerInteraction*					self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( _sizeOfSelf >= sizeof(lucSwarmViewerInteraction) );
	self = (lucSwarmViewerInteraction*) _lucWindowInteraction_New(  LUCWINDOWINTERACTION_PASSARGS  );
	
	return self;
}

void _lucSwarmViewerInteraction_Init( 
		lucSwarmViewerInteraction*                                      self  ) 
{
}

void _lucSwarmViewerInteraction_Delete( void* SwarmViewerInteraction ) {
	lucSwarmViewerInteraction*  self = (lucSwarmViewerInteraction*)SwarmViewerInteraction;

	_lucWindowInteraction_Delete( self );
}

void _lucSwarmViewerInteraction_Print( void* SwarmViewerInteraction, Stream* stream ) {
	lucSwarmViewerInteraction*  self = (lucSwarmViewerInteraction*)SwarmViewerInteraction;

	_lucWindowInteraction_Print( self, stream );
}

void* _lucSwarmViewerInteraction_Copy( void* SwarmViewerInteraction, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap) {
	lucSwarmViewerInteraction*  self = (lucSwarmViewerInteraction*)SwarmViewerInteraction;
	lucSwarmViewerInteraction* newSwarmViewerInteraction;

	newSwarmViewerInteraction = _lucWindowInteraction_Copy( self, dest, deep, nameExt, ptrMap );

	/* TODO */
	abort();

	return (void*) newSwarmViewerInteraction;
}


void* _lucSwarmViewerInteraction_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                                  _sizeOfSelf = sizeof(lucSwarmViewerInteraction);
	Type                                                          type = lucSwarmViewerInteraction_Type;
	Stg_Class_DeleteFunction*                                  _delete = _lucSwarmViewerInteraction_Delete;
	Stg_Class_PrintFunction*                                    _print = _lucSwarmViewerInteraction_Print;
	Stg_Class_CopyFunction*                                      _copy = NULL;
	Stg_Component_DefaultConstructorFunction*      _defaultConstructor = _lucSwarmViewerInteraction_DefaultNew;
	Stg_Component_ConstructFunction*                        _construct = _lucSwarmViewerInteraction_AssignFromXML;
	Stg_Component_BuildFunction*                                _build = _lucSwarmViewerInteraction_Build;
	Stg_Component_InitialiseFunction*                      _initialise = _lucSwarmViewerInteraction_Initialise;
	Stg_Component_ExecuteFunction*                            _execute = _lucSwarmViewerInteraction_Execute;
	Stg_Component_DestroyFunction*                            _destroy = _lucSwarmViewerInteraction_Destroy;
	lucWindowInteraction_MouseMotionFunction*             _mouseMotion = _lucSwarmViewerInteraction_MouseMotion;
	lucWindowInteraction_MouseClickFunction*               _mouseClick = _lucSwarmViewerInteraction_MouseClick;
	lucWindowInteraction_MouseMessageFunction*           _mouseMessage = _lucSwarmViewerInteraction_MouseMessage;
	lucWindowInteraction_KeyboardEventFunction*         _keyboardEvent = _lucSwarmViewerInteraction_KeyboardEvent;
	lucWindowInteraction_KeyboardMessageFunction*     _keyboardMessage = _lucSwarmViewerInteraction_KeyboardMessage;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = ZERO;

	return (void*) _lucSwarmViewerInteraction_New(  LUCSWARMVIEWERINTERACTION_PASSARGS  );
}

void _lucSwarmViewerInteraction_AssignFromXML( void* SwarmViewerInteraction, Stg_ComponentFactory* cf, void* data ){
	lucSwarmViewerInteraction*  self = (lucSwarmViewerInteraction*)SwarmViewerInteraction;

	/* Construct Parent */
	_lucWindowInteraction_AssignFromXML( self, cf, data );
	
	_lucSwarmViewerInteraction_Init( self );
}

void _lucSwarmViewerInteraction_Build( void* renderingEngine, void* data ) {}
void _lucSwarmViewerInteraction_Initialise( void* renderingEngine, void* data ) {}
void _lucSwarmViewerInteraction_Execute( void* renderingEngine, void* data ) {}
void _lucSwarmViewerInteraction_Destroy( void* renderingEngine, void* data ) {}

/* This component doesn't use any mouse interaction */
void _lucSwarmViewerInteraction_MouseMotion( void* windowInteraction, lucWindow* window, lucMouseButton button, Pixel_Index xpos, Pixel_Index ypos, Pixel_Index startx, Pixel_Index starty) {}
void _lucSwarmViewerInteraction_MouseClick( void* windowInteraction, lucWindow* window, lucMouseButton button, lucMouseState state, Pixel_Index xpos, Pixel_Index ypos) { } 
void _lucSwarmViewerInteraction_MouseMessage( void* windowInteraction, Stream* stream ) { }
void _lucSwarmViewerInteraction_KeyboardEvent( void* windowInteraction, lucWindow* window, char key, Pixel_Index xpos, Pixel_Index ypos) {
	lucSwarmViewerInteraction*   self = (lucSwarmViewerInteraction*) windowInteraction;
	lucViewportInfo*            viewportInfo;
	lucViewport*                viewport;
	Coord                       coord;
	lucDrawingObject*           object;
	lucSwarmViewer*	            swarmViewer;
	Stream*                     stream = Journal_MyStream( Info_Type, self );
	DrawingObject_Index         object_I;
	DrawingObject_Index         objectCount;

//	unsigned nodeNumber;
//	unsigned elementNumber;
	
	/* Declare a rank index */
	/* TODO Disabled for now 
	static Partition_Index  rank_I = 0;
	*/

	
	/* Stuff to construct the layout */
	/*MeshDecomp*		decomp;
	MeshLayout*		meshLayout;
	Partition_Index		maxRank = 0;*/
	
	Particle_Index lParticle_I;

	/* This function works when the key pressed is 'n' or 'e'*/
        if ( ( key != 'b' ) )
		return;
	
	/* Find which viewport the user clicked on */
	viewportInfo = lucWindow_GetViewportInfoByPixel( window, xpos, ypos );
	if ( viewportInfo == NULL )
		return; /* If the user hasn't clicked on a viewport at all, then return */

	viewport = viewportInfo->viewport;

	
	/* Loop through lucSwarmViewer that are registered on the viewport */
	objectCount   = lucDrawingObject_Register_GetCount( viewport->drawingObject_Register );
	for ( object_I = 0 ; object_I < objectCount ; object_I++ ) {
		object = lucDrawingObject_Register_GetByIndex( viewport->drawingObject_Register, object_I );

		/* Check if this drawing object is a scalar field */
		if ( Stg_Class_IsInstance( object, lucSwarmViewer_Type ) ) {
			swarmViewer = (lucSwarmViewer*) object;
			/*meshLayout =  swarmViewer->mesh->layout;
			decomp = meshLayout->decomp;
			maxRank = decomp->procsInUse;	*/

		
		        /* TODO  Window interaction do not work well in parallel... As the rank switching feature does not
			make sense in serial , it is disabled for now 
		        */
			if ( key == 'b' ){
			        /* Prints out the particle number */
				/* Get spatial coordinate that the user clicked on */
				lucViewportInfo_GetCoordFromPixel( viewportInfo, xpos, ypos, coord );
				lucSwarmViewerBase_FindParticleLocalIndex(swarmViewer, coord, &lParticle_I);
			//	lucSwarmViewer_ClosestNode( swarmViewer, coord, &nodeNumber);
				Journal_Printf( stream, "Particle number is %d \n", lParticle_I );
			}
		}
	}
	return;
}
void _lucSwarmViewerInteraction_KeyboardMessage( void* windowInteraction, Stream* stream ) {
	Journal_Printf( stream, "b:                            Print the particle number .\n");
}




