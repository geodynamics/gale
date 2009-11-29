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
** $Id: SwarmSquares.c 791 2008-09-01 02:09:06Z JulianGiordani $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>

#include <glucifer/Base/Base.h>
#include <glucifer/RenderingEngines/RenderingEngines.h>

#include "types.h"
#include "OpenGLDrawingObject.h"
#include "SwarmViewerBase.h"
#include "SwarmViewer.h"
#include "SwarmSquares.h"

#include <assert.h>
#include <gl.h>
#include <glu.h>
#include <string.h>

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type lucSwarmSquares_Type = "lucSwarmSquares";

/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
lucSwarmSquares* _lucSwarmSquares_New(  LUCSWARMSQUARES_DEFARGS  ) 
{
	lucSwarmSquares*					self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( _sizeOfSelf >= sizeof(lucSwarmSquares) );
	self = (lucSwarmSquares*) _lucSwarmViewerBase_New(  LUCSWARMVIEWERBASE_PASSARGS  );
	
	return self;
}

void _lucSwarmSquares_Init( 
		lucSwarmSquares*                                   self,
		Name                                               colourVariableName,
		lucColourMap*                                      colourMap,
		Name                                               normalVariableName,
		Name                                               planeVectorVariableName,
		Name                                               lengthVariableName,
		double                                             length )
{
	self->colourMap           = colourMap;
	self->colourVariableName  = colourVariableName;
	self->normalVariableName      = normalVariableName;
	self->planeVectorVariableName = planeVectorVariableName;
	self->lengthVariableName      = lengthVariableName;
	self->length                  = length;
}

void _lucSwarmSquares_Delete( void* drawingObject ) {
	lucSwarmSquares*  self = (lucSwarmSquares*)drawingObject;

	_lucSwarmViewerBase_Delete( self );
}

void _lucSwarmSquares_Print( void* drawingObject, Stream* stream ) {
	lucSwarmSquares*  self = (lucSwarmSquares*)drawingObject;

	_lucSwarmViewerBase_Print( self, stream );
}

void* _lucSwarmSquares_Copy( void* drawingObject, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap) {
	lucSwarmSquares*  self = (lucSwarmSquares*)drawingObject;
	lucSwarmSquares* newDrawingObject;

	newDrawingObject = _lucSwarmViewerBase_Copy( self, dest, deep, nameExt, ptrMap );

	/* TODO */
	abort();

	return (void*) newDrawingObject;
}


void* _lucSwarmSquares_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                                     _sizeOfSelf = sizeof(lucSwarmSquares);
	Type                                                             type = lucSwarmSquares_Type;
	Stg_Class_DeleteFunction*                                     _delete = _lucSwarmSquares_Delete;
	Stg_Class_PrintFunction*                                       _print = _lucSwarmSquares_Print;
	Stg_Class_CopyFunction*                                         _copy = NULL;
	Stg_Component_DefaultConstructorFunction*         _defaultConstructor = _lucSwarmSquares_DefaultNew;
	Stg_Component_ConstructFunction*                           _construct = _lucSwarmSquares_AssignFromXML;
	Stg_Component_BuildFunction*                                   _build = _lucSwarmSquares_Build;
	Stg_Component_InitialiseFunction*                         _initialise = _lucSwarmSquares_Initialise;
	Stg_Component_ExecuteFunction*                               _execute = _lucSwarmSquares_Execute;
	Stg_Component_DestroyFunction*                               _destroy = _lucSwarmSquares_Destroy;
	lucDrawingObject_SetupFunction*                                _setup = _lucSwarmSquares_Setup;
	lucDrawingObject_DrawFunction*                                  _draw = _lucSwarmSquares_Draw;
	lucDrawingObject_CleanUpFunction*                            _cleanUp = _lucSwarmSquares_CleanUp;
	lucOpenGLDrawingObject_BuildDisplayListFunction*    _buildDisplayList = _lucSwarmSquares_BuildDisplayList;
	lucSwarmViewerBase_PlotParticleFunction*                _plotParticle = _lucSwarmSquares_PlotParticle;
	lucSwarmViewerBase_SetParticleColourFunction*      _setParticleColour = _lucSwarmViewerBase_SetParticleColourDefault;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = ZERO;

	return (void*) _lucSwarmSquares_New(  LUCSWARMSQUARES_PASSARGS  );
}

void _lucSwarmSquares_AssignFromXML( void* drawingObject, Stg_ComponentFactory* cf, void* data ){
	lucSwarmSquares*  self = (lucSwarmSquares*)drawingObject;
	lucColourMap*           colourMap;
	Name                    colourVariableName;
	
	/* Construct Parent */
	_lucSwarmViewerBase_AssignFromXML( self, cf, data );
	
	colourMap     =  Stg_ComponentFactory_ConstructByKey(  cf,  self->name,  "ColourMap", lucColourMap,      False, data  ) ;
	colourVariableName = Stg_ComponentFactory_GetString( cf, self->name, "ColourVariable", "" );

	_lucSwarmSquares_Init( 
			self,
			colourVariableName,
			colourMap,
			Stg_ComponentFactory_GetString( cf, self->name, "NormalVariable", "" ),
			Stg_ComponentFactory_GetString( cf, self->name, "PlaneVectorVariable", "" ),
			Stg_ComponentFactory_GetString( cf, self->name, "LengthVariable", "" ),
			Stg_ComponentFactory_GetDouble( cf, self->name, "length", 0.2 ) );
}

void _lucSwarmSquares_Build( void* drawingObject, void* data ) {}
void _lucSwarmSquares_Initialise( void* drawingObject, void* data ) {
	lucSwarmSquares*         self                   = (lucSwarmSquares*)drawingObject;
	SwarmVariable_Register*  swarmVariable_Register = self->swarm->swarmVariable_Register;
	Stream*                  errorStream            = Journal_MyStream( Error_Type, self );

	_lucSwarmViewerBase_Initialise( self, data );

	if ( 0 != strcmp( self->colourVariableName, "" ) ) {
		self->colourVariable  = SwarmVariable_Register_GetByName( swarmVariable_Register, self->colourVariableName );
		Journal_Firewall( self->colourVariable != NULL, errorStream,
				  "Error - for gLucifer drawing object \"%s\" - in %s(): colour Variable name given was \"%s\", "
				  "but no corresponding SwarmVariable found in the register for swarm \"%s\".\n",
				  self->name, __func__, self->colourVariableName, self->swarm->name );
		
		Stg_Component_Build( self->colourVariable, data, False );
		Stg_Component_Initialise( self->colourVariable, data, False );

	}
	
	if ( 0 != strcmp( self->normalVariableName, "" ) ) {
		self->normalVariable = SwarmVariable_Register_GetByName( swarmVariable_Register, self->normalVariableName );
		Journal_Firewall( self->normalVariable != NULL, errorStream,
					"Error - for gLucifer drawing object \"%s\" - in %s(): normal Variable name given was \"%s\", "
					"but no corresponding SwarmVariable found in the register for swarm \"%s\".\n",
					self->name, __func__, self->normalVariableName, self->swarm->name );
		
		Stg_Component_Build( self->normalVariable, data, False );
		Stg_Component_Initialise( self->normalVariable, data, False );
	}

	/*
	self->planeVectorVariable = SwarmVariable_Register_GetByName( swarmVariable_Register, self->planeVectorVariableName );
	Journal_Firewall( self->planeVectorVariable != NULL, errorStream, 
			"Error in func %s for %s '%s' - Cannot find SwarmVariable %s to be variable for plane vector.\n", 
			__func__, self->type, self->name, self->planeVectorVariableName );
	*/		
	
	if ( 0 != strcmp( self->lengthVariableName, "" ) ) {
		self->lengthVariable = SwarmVariable_Register_GetByName( swarmVariable_Register, self->lengthVariableName );
		Journal_Firewall( self->lengthVariable != NULL, errorStream,
					"Error - for gLucifer drawing object \"%s\" - in %s(): normal Variable name given was \"%s\", "
					"but no corresponding SwarmVariable found in the register for swarm \"%s\".\n",
					self->name, __func__, self->lengthVariableName, self->swarm->name );
		
		Stg_Component_Build( self->lengthVariable, data, False );
		Stg_Component_Initialise( self->lengthVariable, data, False );
	}

}
void _lucSwarmSquares_Execute( void* drawingObject, void* data ) {}
void _lucSwarmSquares_Destroy( void* drawingObject, void* data ) {}

void _lucSwarmSquares_Setup( void* drawingObject, void* _context ) {
	lucSwarmSquares_UpdateVariables( drawingObject );
	_lucSwarmViewerBase_Setup( drawingObject, _context );
}
	
void _lucSwarmSquares_Draw( void* drawingObject, lucWindow* window, lucViewportInfo* viewportInfo, void* _context ) {
	_lucSwarmViewerBase_Draw( drawingObject, window, viewportInfo, _context );
}

void _lucSwarmSquares_CleanUp( void* drawingObject, void* _context ) {
	_lucSwarmViewerBase_CleanUp( drawingObject, _context );
}

void _lucSwarmSquares_BuildDisplayList( void* drawingObject, void* _context ) {
	lucSwarmViewer*          self                = (lucSwarmViewer*)drawingObject;

     /* 
		Hack to allow the transparency to work properly 
	 	See : http://www.oreillynet.com/pub/a/network/2000/06/23/magazine/opengl_render.html?page=2 
	  	glBlendFunc(GL_SRC_ALPHA, GL_ONE);
		
		This isn't so great either ... it's hard to overlay darker colours on 
		light with this choice of blending.
	 */
	glEnable(GL_BLEND);
	
	
	_lucSwarmViewerBase_BuildDisplayList( self, _context );
	
	/* glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); */
}

void _lucSwarmSquares_PlotParticle( void* drawingObject, void* _context, Particle_Index lParticle_I ) {
	lucSwarmSquares*         self                = (lucSwarmSquares*)drawingObject;
	DomainContext*   context             = (DomainContext*) _context;
	GlobalParticle*          particle            = (GlobalParticle*)Swarm_ParticleAt( self->swarm, lParticle_I );
	SwarmVariable*           lengthVariable      = self->lengthVariable;
	double*                  coord               = particle->coord;
	double                   length              = self->length;
	XYZ                      normal              = { 0, 0, 0 };
	/*XYZ                      planeVector         = { 0, 0, 0 };*/

	SwarmVariable_ValueAt( self->normalVariable, lParticle_I, normal );
	/* SwarmVariable_ValueAt( self->planeVectorVariable, lParticle_I, planeVector );*/

	if ( lengthVariable )
		SwarmVariable_ValueAt( lengthVariable, lParticle_I, &length );

	
	/*  The fat square has a pizza box shape ... i.e. edges and two faces 
		which have opposite normals */
	
	luc_OpenGlFatSquare( context->dim, coord, normal, NULL, length, length * 0.1);  
	/* luc_OpenGlSquare( context->dim, coord, normal, NULL, length); */

}

void lucSwarmSquares_UpdateVariables( void* drawingObject ) {
	lucSwarmSquares*          self                = (lucSwarmSquares*)drawingObject;
	
	lucSwarmViewerBase_UpdateVariables( drawingObject ) ;

	if ( self->normalVariable && self->normalVariable->variable ) {
		Variable_Update( self->normalVariable->variable );
	}
	if ( self->planeVectorVariable && self->planeVectorVariable->variable ) {
		Variable_Update( self->planeVectorVariable->variable );
	}
	if ( self->lengthVariable && self->lengthVariable->variable ) {
		Variable_Update( self->lengthVariable->variable );
	}
}


