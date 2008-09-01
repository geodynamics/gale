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
** $Id: SwarmViewer.c 595 2006-07-18 06:53:06Z LukeHodkinson $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#ifdef GLUCIFER_USE_PICELLERATOR
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#endif

#include <glucifer/Base/Base.h>
#include <glucifer/RenderingEngines/RenderingEngines.h>

#ifdef HAVE_GL2PS
	#include <gl2ps.h>
#endif

#include "types.h"
#include "OpenGLDrawingObject.h"
#include "SwarmViewerBase.h"
#include "SwarmViewer.h"

#include <assert.h>
#ifdef HAVE_OPENGL_FRAMEWORK
	#include <OpenGL/gl.h>
	#include <OpenGL/glu.h>
#else
	#include <gl.h>
	#include <glu.h>
#endif
#include <string.h>

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type lucSwarmViewerBase_Type = "lucSwarmViewerBase";

/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
lucSwarmViewerBase* _lucSwarmViewerBase_New( 
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
	lucSwarmViewerBase_PlotParticleFunction*           _plotParticle,
	lucSwarmViewerBase_SetParticleColourFunction*      _setParticleColour,
	Name                                               name ) 
{
	
	lucSwarmViewerBase*					self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( sizeOfSelf >= sizeof(lucSwarmViewerBase) );
	self = (lucSwarmViewerBase*) _lucOpenGLDrawingObject_New( 
		sizeOfSelf,
		type, 
		_delete,
		_print,
		_copy,
		_defaultConstructor,
		_construct,
		_build,
		_initialise,
		_execute,
		_destroy,
		_setup,
		_draw,
		_cleanUp,
		_buildDisplayList,
		name );

	self->_plotParticle = _plotParticle;
	self->_setParticleColour = _setParticleColour;

	return self;
}

void _lucSwarmViewerBase_Init( 
	lucSwarmViewerBase*                                          self,
	Swarm*                                                       swarm,
	Name                                                         colourName,
	Name                                                         colourVariableName,
	lucColourMap*                                                colourMap,
	Name                                                         opacityVariableName,
	Name                                                         maskVariableName,
	lucDrawingObjectMask*                                        mask,
	Bool                                                         drawParticleNumber,
	Bool                                                         particleColour,
	int                                                          subSetEvery,
	Bool                                                         positionRange,
	Coord                                                        minPosition,
	Coord                                                        maxPosition
	)
{
	self->swarm               = swarm;	
	self->colourVariableName  = colourVariableName;
	self->colourMap           = colourMap;
	self->opacityVariableName = opacityVariableName;
	self->maskVariableName    = maskVariableName;
	self->drawParticleNumber  = drawParticleNumber;
	self->sameParticleColour  = particleColour;
	self->subSetEvery         = subSetEvery;
	self->positionRange       = positionRange;

	memcpy( &self->mask, mask, sizeof( lucDrawingObjectMask ) );
	memcpy( &self->minPosition, minPosition , sizeof( Coord ) );
	memcpy( &self->maxPosition, maxPosition , sizeof( Coord ) );
	
	lucColour_FromString( &self->colour, colourName );
}

void _lucSwarmViewerBase_Delete( void* drawingObject ) {
	lucSwarmViewerBase*  self = (lucSwarmViewerBase*)drawingObject;
	
	_lucOpenGLDrawingObject_Delete( self );
      
}

void _lucSwarmViewerBase_Print( void* drawingObject, Stream* stream ) {
	lucSwarmViewerBase*  self = (lucSwarmViewerBase*)drawingObject;

	_lucOpenGLDrawingObject_Print( self, stream );
}

void* _lucSwarmViewerBase_Copy( void* drawingObject, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap) {
	lucSwarmViewerBase*  self = (lucSwarmViewerBase*)drawingObject;
	lucSwarmViewerBase* newDrawingObject;

	newDrawingObject = _lucOpenGLDrawingObject_Copy( self, dest, deep, nameExt, ptrMap );

	/* TODO */
	abort();

	return (void*) newDrawingObject;
}

void _lucSwarmViewerBase_Construct( void* drawingObject, Stg_ComponentFactory* cf, void* data ){
	lucSwarmViewerBase*     self = (lucSwarmViewerBase*)drawingObject;
	Swarm*                  swarm;	
	Name                    colourVariableName;
	lucColourMap*           colourMap;
	Name                    opacityVariableName;
	Name                    maskVariableName;
	Bool                    drawParticleNumber;
	Bool                    sameParticleColour;
	int 			subSetEvery;
	Bool                    positionRange;
	Coord                   minPosition;
	Coord                   maxPosition;
	lucDrawingObjectMask    mask;

	/* Construct Parent */
	_lucOpenGLDrawingObject_Construct( self, cf, data );

	swarm = Stg_ComponentFactory_ConstructByKey( cf, self->name, "Swarm", Swarm, True, data ) ;
	
	/* This drawing object will only work for swarms with Global Particle Layouts 
	 *  HACK - Adding in check for Gauss particle Layout here because this can be global too */
	Journal_Firewall(
		swarm->particleLayout->coordSystem == GlobalCoordSystem || Stg_Class_IsInstance( swarm->particleLayout, GaussParticleLayout_Type ),
		Journal_MyStream( Error_Type, self ),
		"In func %s, unable to visualise swarm %s because it uses a local coord system layout %s of type %s.\n",
		__func__,
		swarm->name,
		swarm->particleLayout->name,
		swarm->particleLayout->type );

	colourVariableName  = Stg_ComponentFactory_GetString( cf, self->name, "ColourVariable", "" );
	colourMap           =  Stg_ComponentFactory_ConstructByKey( cf, self->name, "ColourMap", lucColourMap, False, data );
	opacityVariableName = Stg_ComponentFactory_GetString( cf, self->name, "OpacityVariable", "" );
	maskVariableName = Stg_ComponentFactory_GetString( cf, self->name, "MaskVariable", "" );
	
	drawParticleNumber = Stg_ComponentFactory_GetBool( cf, self->name, "drawParticleNumber", False );
	sameParticleColour = Stg_ComponentFactory_GetBool( cf, self->name, "sameParticleColour", False );

	subSetEvery = Stg_ComponentFactory_GetInt( cf, self->name, "subSetEvery", 1 );
	positionRange = Stg_ComponentFactory_GetBool( cf, self->name, "positionRange", False );

	/* Memory allocation */
       	minPosition[I_AXIS]  = Stg_ComponentFactory_GetDouble( cf, self->name, "minPositionX", -100000.0 );
	minPosition[J_AXIS]  = Stg_ComponentFactory_GetDouble( cf, self->name, "minPositionY", -100000.0 );
	minPosition[K_AXIS]  = Stg_ComponentFactory_GetDouble( cf, self->name, "minPositionZ", -100000.0 );

	maxPosition[I_AXIS]  = Stg_ComponentFactory_GetDouble( cf, self->name, "maxPositionX", 100000.0 );
	maxPosition[J_AXIS]  = Stg_ComponentFactory_GetDouble( cf, self->name, "maxPositionY", 100000.0 );
	maxPosition[K_AXIS]  = Stg_ComponentFactory_GetDouble( cf, self->name, "maxPositionZ", 100000.0 );


	lucDrawingObjectMask_Construct( &mask, self->name, cf, data );

	_lucSwarmViewerBase_Init( 
		self, 
		swarm,
		Stg_ComponentFactory_GetString( cf, self->name, "colour", "black" ),
		colourVariableName,
		colourMap,
		opacityVariableName,
		maskVariableName,
		&mask,
		drawParticleNumber,
		sameParticleColour, 
		subSetEvery,
		positionRange,
		minPosition,
		maxPosition);
}

void _lucSwarmViewerBase_Build( void* drawingObject, void* data ) {
	lucSwarmViewerBase*          self                   = (lucSwarmViewerBase*)drawingObject;

	_lucOpenGLDrawingObject_Build( self, data );
}

void _lucSwarmViewerBase_Initialise( void* drawingObject, void* data ) {
	lucSwarmViewerBase*          self                   = (lucSwarmViewerBase*)drawingObject;
	SwarmVariable_Register*  swarmVariable_Register = self->swarm->swarmVariable_Register;
	Stream*                  errorStream               = Journal_Register( Error_Type, self->type );
	
	/* Initialise Parent */
	_lucOpenGLDrawingObject_Initialise( self, data );

	if ( 0 != strcmp( self->colourVariableName, "" ) ) {
		self->colourVariable = SwarmVariable_Register_GetByName( swarmVariable_Register, self->colourVariableName );
		Journal_Firewall( self->colourVariable != NULL, errorStream,
					"Error - for gLucifer drawing object \"%s\" - in %s(): colour Variable name given was \"%s\", "
					"but no corresponding SwarmVariable found in the register for swarm \"%s\".\n",
					self->name, __func__, self->colourVariableName, self->swarm->name );
		
		Stg_Component_Build( self->colourVariable, data, False );
		Stg_Component_Initialise( self->colourVariable, data, False );
	}
	if ( 0 != strcmp( self->opacityVariableName, "" ) ) {
		self->opacityVariable = SwarmVariable_Register_GetByName( swarmVariable_Register, self->opacityVariableName );
		Journal_Firewall( self->opacityVariable != NULL, errorStream,
					"Error - for gLucifer drawing object \"%s\" - in %s(): opacity Variable name given was \"%s\", "
					"but no corresponding SwarmVariable found in the register for swarm \"%s\".\n",
					self->name, __func__, self->opacityVariableName, self->swarm->name );
		
		Stg_Component_Build( self->opacityVariable, data, False );
		Stg_Component_Initialise( self->opacityVariable, data, False );
	}
	if ( 0 != strcmp( self->maskVariableName, "" ) ) {
		self->maskVariable    = SwarmVariable_Register_GetByName( swarmVariable_Register, self->maskVariableName );
		Journal_Firewall( self->maskVariable != NULL, errorStream,
				"Error - for gLucifer drawing object \"%s\" - in %s(): mask Variable name given was \"%s\", "
					"but no corresponding SwarmVariable found in the register for swarm \"%s\".\n",
					self->name, __func__, self->maskVariableName, self->swarm->name );
		
		Stg_Component_Build( self->maskVariable, data, False );
		Stg_Component_Initialise( self->maskVariable, data, False );
	}	
}


void _lucSwarmViewerBase_Execute( void* drawingObject, void* data ) {}
void _lucSwarmViewerBase_Destroy( void* drawingObject, void* data ) {}

void _lucSwarmViewerBase_Setup( void* drawingObject, void* _context ) {
	lucSwarmViewerBase*          self                = (lucSwarmViewerBase*)drawingObject;
	lucColourMap*                colourMap           = self->colourMap;
	SwarmVariable*               colourVariable      = self->colourVariable;

	lucSwarmViewerBase_UpdateVariables( self );
	
	/* Scale Colour Map */
	if ( colourVariable && colourMap ) {
		lucColourMap_CalibrateFromSwarmVariable( colourMap, colourVariable );
	}
	
	_lucOpenGLDrawingObject_Setup( self, _context );
}

void _lucSwarmViewerBase_Draw( void* drawingObject, lucWindow* window, lucViewportInfo* viewportInfo, void* _context ) {
	lucSwarmViewerBase*          self          = (lucSwarmViewerBase*)drawingObject;

	_lucOpenGLDrawingObject_Draw( self, window, viewportInfo, _context );
}


void _lucSwarmViewerBase_CleanUp( void* drawingObject, void* context ) {
	lucSwarmViewerBase*          self          = (lucSwarmViewerBase*)drawingObject;
	
	_lucOpenGLDrawingObject_CleanUp( self, context );
}

void _lucSwarmViewerBase_BuildDisplayList( void* drawingObject, void* _context ) {
	lucSwarmViewerBase*      self                = (lucSwarmViewerBase*)drawingObject;
	Swarm*                   swarm               = self->swarm;
	SwarmVariable*           maskVariable        = self->maskVariable;
	Particle_Index           particleLocalCount  = swarm->particleLocalCount;
	Particle_Index           lParticle_I;
	double                   maskResult;
	int                      subSetEvery         = self->subSetEvery;
	Bool                     positionRange       = self->positionRange;
	GlobalParticle*          particle;
	double*                  coord;
	double*                  minPosition;
	double*                  maxPosition;

	minPosition = self->minPosition;
	maxPosition = self->maxPosition;

	/* Take one of subSetEvery particle */
      	for ( lParticle_I = 0 ; lParticle_I < particleLocalCount ; lParticle_I+=subSetEvery ){

		/* Test to see if this particle should be drawn */
		if ( maskVariable ) {
			SwarmVariable_ValueAt( maskVariable, lParticle_I, &maskResult );
			if ( lucDrawingObjectMask_Test( &self->mask, maskResult ) == False )
				continue;
		}

		/* Sets the colour for the particle */
		self->_setParticleColour( self, _context, lParticle_I );
		
	        /* Check if needed that the particle falls into the right position range */
		if(positionRange){
			particle            = (GlobalParticle*)Swarm_ParticleAt( self->swarm, lParticle_I );
			coord               = particle->coord;

			if( coord[0] < minPosition[I_AXIS] || coord[1] < minPosition[J_AXIS] || 
			    coord[0] > maxPosition[I_AXIS] || coord[1] > maxPosition[J_AXIS] )
			{
				continue;
			}

			if( ((DomainContext*)_context)->dim == 3 ) {
				if( coord[2] < minPosition[K_AXIS] || coord[2] > maxPosition[K_AXIS] )
					continue;
			}
		}

		/* Plot the particle using the function given by the concrete class */
		self->_plotParticle( self, _context, lParticle_I );
	}	

	/* Go through the list of the particles again and write the text of the numbers next to each other */
	if ( self->drawParticleNumber ) {
		lucSwarmViewBase_DrawParticleNumbers( self, _context );
	}
}
void lucSwarmViewBase_DrawParticleNumbers( void* drawingObject, void* _context ) {
	abort();
}

void _lucSwarmViewerBase_PlotParticleNumber( void* drawingObject, void* _context, Particle_Index lParticle_I, lucColour colour ) {
	lucSwarmViewerBase*      self                = (lucSwarmViewerBase*)drawingObject;
	DomainContext*   context             = (DomainContext*) _context;
	GlobalParticle*          particle            = (GlobalParticle*)Swarm_ParticleAt( self->swarm, lParticle_I );
	double*                  coord               = particle->coord;
        Name particle_number;
	Stg_asprintf(&particle_number, "%d", lParticle_I );
	
	glDisable(GL_LIGHTING); /*if the lighting is not disabled, the colour won't appear for the numbers*/

	if (context->dim == 2)
		glRasterPos2f( (float)coord[0] + 0.025, (float)coord[1] );		
	else   
		glRasterPos3f( (float)coord[0] + 0.025, (float)coord[1], (float)coord[2] );

	lucPrintString( particle_number );

	Memory_Free(particle_number);
           
	/* Put back settings */
	glEnable(GL_LIGHTING);

}

void lucSwarmViewerBase_FindParticleLocalIndex(void *drawingObject, Coord coord, Particle_Index *lParticle_I){
	lucSwarmViewerBase*      self  = (lucSwarmViewerBase*) drawingObject;
	Swarm*               swarm = self->swarm;
	Dimension_Index      dim = self->swarm->dim;
	Particle_InCellIndex cParticle_I;
	Cell_LocalIndex      lCell_I;
	GlobalParticle       testParticle;
	double               minDistance;
	
	/* Find cell this coordinate is in */
	memcpy( testParticle.coord, coord, sizeof(Coord) );
	/* First specify the particle doesn't have an owning cell yet, so as
	   not to confuse the search algorithm */
	testParticle.owningCell = swarm->cellDomainCount;
	lCell_I = CellLayout_CellOf( swarm->cellLayout, &testParticle );

	/* Test if this cell is on this processor - if not then bail */
	if (lCell_I >= swarm->cellLocalCount){
		*lParticle_I = (Particle_Index) -1;
		return;
	}

	/* Find Closest Particle in this Cell */
	cParticle_I = Swarm_FindClosestParticleInCell( swarm, lCell_I, dim, coord, &minDistance );

	/* Convert to Local Particle Index */
	*lParticle_I = swarm->cellParticleTbl[ lCell_I ][ cParticle_I ];
}



void lucSwarmViewerBase_UpdateVariables( void* drawingObject ) {
	lucSwarmViewerBase*          self                = (lucSwarmViewerBase*)drawingObject;

	if ( self->opacityVariable && self->opacityVariable->variable ) {
		Variable_Update( self->opacityVariable->variable );
	}
	if ( self->maskVariable && self->maskVariable->variable ) {
		Variable_Update( self->maskVariable->variable );
	}
}

/* Default Swarm Viewer Implementation */
void _lucSwarmViewerBase_SetParticleColourDefault( void* drawingObject, void* context, Particle_Index lParticle_I ) {
	lucSwarmViewer*          self                = (lucSwarmViewer*) drawingObject;
	SwarmVariable*           colourVariable      = self->colourVariable;
	SwarmVariable*           opacityVariable     = self->opacityVariable;
	lucColourMap*            colourMap           = self->colourMap;
	double                   colourValue;
	double                   opacity;
	lucColour                colour;

 	/* Get colour value if there is a colourVariable and a colourMap */
	if ( colourVariable && colourMap ) {
		SwarmVariable_ValueAt( colourVariable, lParticle_I, &colourValue );
		lucColourMap_GetColourFromValue( colourMap, colourValue, &colour );
	}
	else {
		/* Set the default Colour */
		memcpy( &colour, &self->colour, sizeof(lucColour) );
	}
	
	/* Get Opacity Value */
	if ( opacityVariable ){
		SwarmVariable_ValueAt( opacityVariable, lParticle_I, &opacity );
		colour.opacity = (float)opacity;	
	}

	lucColour_SetOpenGLColour( &colour );
}

