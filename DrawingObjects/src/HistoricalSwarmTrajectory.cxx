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
** $Id: HistoricalSwarmTrajectory.c 791 2008-09-01 02:09:06Z JulianGiordani $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>

#include <glucifer/Base/Base.h>
#include <glucifer/RenderingEngines/RenderingEngines.h>

#include "types.h"
#include "OpenGLDrawingObject.h"
#include "HistoricalSwarmTrajectory.h"

#include <assert.h>
#include <gl.h>
#include <glu.h>
#include <string.h>

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type lucHistoricalSwarmTrajectory_Type = "lucHistoricalSwarmTrajectory";

/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
lucHistoricalSwarmTrajectory* _lucHistoricalSwarmTrajectory_New(  LUCHISTORICALSWARMTRAJECTORY_DEFARGS  ) 
{
	lucHistoricalSwarmTrajectory*					self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( _sizeOfSelf >= sizeof(lucHistoricalSwarmTrajectory) );
	self = (lucHistoricalSwarmTrajectory*) _lucOpenGLDrawingObject_New(  LUCOPENGLDRAWINGOBJECT_PASSARGS  );

	return self;
}

void _lucHistoricalSwarmTrajectory_Init( 
		lucHistoricalSwarmTrajectory*                                self,
		Swarm*                                                       swarm,
		lucColourMap*                                                colourMap,
		Name                                                         colourName,
		float                                                        lineWidth,
		unsigned int 						     historySteps,
		double							     historyTime )
{
	self->swarm               = swarm;
	self->colourMap           = colourMap;
	self->lineWidth           = lineWidth;
	self->historySteps	  = historySteps;
	self->historyTime	  = historyTime;

	self->particleExtHandle = ExtensionManager_Add( swarm->particleExtensionMgr, (Name)self->type, sizeof( lucHistoricalSwarmTrajectory_ParticleExt ) + ( historySteps+3 )*sizeof( Coord )  );
	
	lucColour_FromString( &self->colour, colourName );

	self->startTimestepIndex = 0;
}

void _lucHistoricalSwarmTrajectory_Delete( void* drawingObject ) {
	lucHistoricalSwarmTrajectory*  self = (lucHistoricalSwarmTrajectory*)drawingObject;

	_lucOpenGLDrawingObject_Delete( self );
}

void _lucHistoricalSwarmTrajectory_Print( void* drawingObject, Stream* stream ) {
	lucHistoricalSwarmTrajectory*  self = (lucHistoricalSwarmTrajectory*)drawingObject;

	_lucOpenGLDrawingObject_Print( self, stream );
}

void* _lucHistoricalSwarmTrajectory_Copy( void* drawingObject, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap) {
	lucHistoricalSwarmTrajectory*  self = (lucHistoricalSwarmTrajectory*)drawingObject;
	lucHistoricalSwarmTrajectory* newDrawingObject;

	newDrawingObject = _lucOpenGLDrawingObject_Copy( self, dest, deep, nameExt, ptrMap );

	/* TODO */
	abort();

	return (void*) newDrawingObject;
}


void* _lucHistoricalSwarmTrajectory_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                                     _sizeOfSelf = sizeof(lucHistoricalSwarmTrajectory);
	Type                                                             type = lucHistoricalSwarmTrajectory_Type;
	Stg_Class_DeleteFunction*                                     _delete = _lucHistoricalSwarmTrajectory_Delete;
	Stg_Class_PrintFunction*                                       _print = _lucHistoricalSwarmTrajectory_Print;
	Stg_Class_CopyFunction*                                         _copy = NULL;
	Stg_Component_DefaultConstructorFunction*         _defaultConstructor = _lucHistoricalSwarmTrajectory_DefaultNew;
	Stg_Component_ConstructFunction*                           _construct = _lucHistoricalSwarmTrajectory_AssignFromXML;
	Stg_Component_BuildFunction*                                   _build = _lucHistoricalSwarmTrajectory_Build;
	Stg_Component_InitialiseFunction*                         _initialise = _lucHistoricalSwarmTrajectory_Initialise;
	Stg_Component_ExecuteFunction*                               _execute = _lucHistoricalSwarmTrajectory_Execute;
	Stg_Component_DestroyFunction*                               _destroy = _lucHistoricalSwarmTrajectory_Destroy;
	lucDrawingObject_SetupFunction*                                _setup = _lucHistoricalSwarmTrajectory_Setup;
	lucDrawingObject_DrawFunction*                                  _draw = _lucHistoricalSwarmTrajectory_Draw;
	lucDrawingObject_CleanUpFunction*                            _cleanUp = _lucHistoricalSwarmTrajectory_CleanUp;
	lucOpenGLDrawingObject_BuildDisplayListFunction*    _buildDisplayList = _lucHistoricalSwarmTrajectory_BuildDisplayList;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

	return (void*) _lucHistoricalSwarmTrajectory_New(  LUCHISTORICALSWARMTRAJECTORY_PASSARGS  );
}

void _lucHistoricalSwarmTrajectory_AssignFromXML( void* drawingObject, Stg_ComponentFactory* cf, void* data ){
	lucHistoricalSwarmTrajectory*   self	= (lucHistoricalSwarmTrajectory*)drawingObject;
	lucColourMap*           	colourMap;
	Swarm*                  	swarm;
	unsigned int			historySteps;
	double 				historyTime;

	/* Construct Parent */
	_lucOpenGLDrawingObject_AssignFromXML( self, cf, data );

	swarm         =  Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"Swarm", Swarm, True, data  );
	colourMap     =  Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"ColourMap", lucColourMap, False, data  );

	historySteps  =  Stg_ComponentFactory_GetUnsignedInt( cf, self->name, (Dictionary_Entry_Key)"historySteps", DEFAULT_STEPS  );
	historyTime   =  Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"historyTime", 0  );

	Journal_Firewall(
			swarm->particleLayout->coordSystem == GlobalCoordSystem,
			Journal_MyStream( Error_Type, self ),
			"In func %s, unable to visualise swarm %s because it uses a local coord system layout %s of type %s.\n",
			__func__,
			swarm->name,
			swarm->particleLayout->name,
			swarm->particleLayout->type );

	_lucHistoricalSwarmTrajectory_Init( 
			self, 
			swarm,
			colourMap,
			Stg_ComponentFactory_GetString( cf, self->name, (Dictionary_Entry_Key)"colour", "black"  ),
			(float) Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"lineWidth", 1.0  ),
			historySteps,
			historyTime );
}

void _lucHistoricalSwarmTrajectory_Build( void* drawingObject, void* data ) {}
void _lucHistoricalSwarmTrajectory_Initialise( void* drawingObject, void* data ) {}

void _lucHistoricalSwarmTrajectory_Execute( void* drawingObject, void* data ) {}

void _lucHistoricalSwarmTrajectory_Destroy( void* drawingObject, void* data ) {
	lucHistoricalSwarmTrajectory*  self            = (lucHistoricalSwarmTrajectory*)drawingObject;
	
	Memory_Free( self->timeAtStep );
}

void _lucHistoricalSwarmTrajectory_Setup( void* drawingObject, void* _context ) {
	lucHistoricalSwarmTrajectory*  self            = (lucHistoricalSwarmTrajectory*)drawingObject;
	AbstractContext*               context         = (AbstractContext*)             _context;
	Swarm*                         swarm           = self->swarm;
	lucColourMap*                  colourMap       = self->colourMap;
	Particle_Index                 lParticle_I;
	GlobalParticle*                particle;
	int                            currentTimestep = context->timeStep;
	lucHistoricalSwarmTrajectory_ParticleExt* particleExt;
	unsigned int		       historySteps    = self->historySteps;
	
	if ( currentTimestep >= historySteps ) {
		self->startTimestepIndex++;
		self->startTimestepIndex %= historySteps;	
		currentTimestep %= historySteps;		
	}
	
	self->timeAtStep =  Memory_Realloc_Array( self->timeAtStep, double, historySteps );

	self->timeAtStep[ currentTimestep ] = context->currentTime;	

	/* Calibrate Colour Map */
	if ( colourMap )
		lucColourMap_SetMinMax( colourMap, self->timeAtStep[ self->startTimestepIndex ], self->timeAtStep[ currentTimestep ] );

	/* Store Current position for each time step */
	for ( lParticle_I = 0 ; lParticle_I < swarm->particleLocalCount ; lParticle_I++ ) {

		particle    = (GlobalParticle*)Swarm_ParticleAt( swarm, lParticle_I );
		particleExt = ExtensionManager_Get( swarm->particleExtensionMgr, particle, self->particleExtHandle );

		particleExt->historyCoordList = particleExt + sizeof( lucHistoricalSwarmTrajectory_ParticleExt );

		memcpy( particleExt->historyCoordList[ currentTimestep ], particle->coord, sizeof(Coord) );
	}
	
	/* Call parent's 'Setup' function */
	_lucOpenGLDrawingObject_Setup( self, _context );
}

void _lucHistoricalSwarmTrajectory_Draw( void* drawingObject, lucWindow* window, lucViewportInfo* viewportInfo, void* _context ) {
	_lucOpenGLDrawingObject_Draw( drawingObject, window, viewportInfo, _context );
}

void _lucHistoricalSwarmTrajectory_CleanUp( void* drawingObject, void* context ) {
	_lucOpenGLDrawingObject_CleanUp( drawingObject, context );
}

void _lucHistoricalSwarmTrajectory_BuildDisplayList( void* drawingObject, void* _context ) {
	lucHistoricalSwarmTrajectory*  self            = (lucHistoricalSwarmTrajectory*)drawingObject;
	DomainContext*         context         = (DomainContext*)       _context;
	Swarm*                         swarm           = self->swarm;
	lucColourMap*                  colourMap       = self->colourMap;
	unsigned int		       historySteps    = self->historySteps;
	double			       historyTime     = self->historyTime;		
	int                            timestep;
	int                            currentTimestep = context->timeStep % historySteps; 
	Particle_Index                 lParticle_I;
	StandardParticle*              particle;
	lucHistoricalSwarmTrajectory_ParticleExt* particleExt;
	float                          offset          = -1.0;
	double*                        coord;
	Dimension_Index                dim             = context->dim;
	double 			       currentTime     = context->currentTime;	

	lucColour_SetOpenGLColour( &self->colour );
	glLineWidth( self->lineWidth );
	
	glDisable(GL_LIGHTING);

	/* Loop over all particles and draw lines according to where each one has been */
	for ( lParticle_I = 0 ; lParticle_I < swarm->particleLocalCount ; lParticle_I++ ) {
		particle    = Swarm_ParticleAt( swarm, lParticle_I );
		particleExt = ExtensionManager_Get( swarm->particleExtensionMgr, particle, self->particleExtHandle );

		glBegin( GL_LINE_STRIP );
		timestep = self->startTimestepIndex ;
		while (True) {
			if ( (currentTime - self->timeAtStep[ timestep ]) <= historyTime || historyTime == 0 ) {
				coord = particleExt->historyCoordList[ timestep ];

				/* Set the colour from the colour map - if we've passed it in */
				if ( colourMap )
					lucColourMap_SetOpenGLColourFromValue( colourMap, self->timeAtStep[ timestep ] );

				if (dim == 2)
					glVertex3f( (float)coord[0], (float)coord[1], offset );
				else
					glVertex3dv( coord );
			}
				
			/* Stop the loop when we have arrived a the current timestep */
			if ( timestep == currentTimestep )
				break;

			/* Adjust current timestep counter so that the list of stored coordinates loops over itself */
			timestep++;
			timestep %= historySteps;
		}
		glEnd();
	}
		
	glEnable(GL_LIGHTING);

}


