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

#include <glucifer/Base/Base.h>
#include <glucifer/RenderingEngines/RenderingEngines.h>

#include "types.h"
#include "OpenGLDrawingObject.h"
#include "SwarmViewerBase.h"
#include "SwarmViewer.h"
#include "SwarmRGBColourViewer.h"

#include <assert.h>
#include <gl.h>
#include <glu.h>
#include <string.h>

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type lucSwarmRGBColourViewer_Type = "lucSwarmRGBColourViewer";

/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
lucSwarmRGBColourViewer* _lucSwarmRGBColourViewer_New(  LUCSWARMRGBCOLOURVIEWER_DEFARGS  ) 
{
	lucSwarmRGBColourViewer*					self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( _sizeOfSelf >= sizeof(lucSwarmRGBColourViewer) );
	self = (lucSwarmRGBColourViewer*) _lucSwarmViewer_New(  LUCSWARMVIEWER_PASSARGS  );

	return self;
}

void _lucSwarmRGBColourViewer_Init( 
		lucSwarmRGBColourViewer*                                     self,
		Name                                                         colourRedVariableName,
		Name                                                         colourGreenVariableName,
		Name                                                         colourBlueVariableName )
{
	self->colourRedVariableName           = colourRedVariableName;
	self->colourGreenVariableName         = colourGreenVariableName;
	self->colourBlueVariableName          = colourBlueVariableName;
}

void _lucSwarmRGBColourViewer_Delete( void* drawingObject ) {
	lucSwarmRGBColourViewer*  self = (lucSwarmRGBColourViewer*)drawingObject;

	_lucSwarmViewer_Delete( self );
}

void _lucSwarmRGBColourViewer_Print( void* drawingObject, Stream* stream ) {
	lucSwarmRGBColourViewer*  self = (lucSwarmRGBColourViewer*)drawingObject;

	_lucSwarmViewer_Print( self, stream );
}

void* _lucSwarmRGBColourViewer_Copy( void* drawingObject, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap) {
	lucSwarmRGBColourViewer*  self = (lucSwarmRGBColourViewer*)drawingObject;
	lucSwarmRGBColourViewer* newDrawingObject;

	newDrawingObject = _lucSwarmViewer_Copy( self, dest, deep, nameExt, ptrMap );

	/* TODO */
	abort();

	return (void*) newDrawingObject;
}


void* _lucSwarmRGBColourViewer_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                                     _sizeOfSelf = sizeof(lucSwarmRGBColourViewer);
	Type                                                             type = lucSwarmRGBColourViewer_Type;
	Stg_Class_DeleteFunction*                                     _delete = _lucSwarmRGBColourViewer_Delete;
	Stg_Class_PrintFunction*                                       _print = _lucSwarmRGBColourViewer_Print;
	Stg_Class_CopyFunction*                                         _copy = NULL;
	Stg_Component_DefaultConstructorFunction*         _defaultConstructor = _lucSwarmRGBColourViewer_DefaultNew;
	Stg_Component_ConstructFunction*                           _construct = _lucSwarmRGBColourViewer_AssignFromXML;
	Stg_Component_BuildFunction*                                   _build = _lucSwarmRGBColourViewer_Build;
	Stg_Component_InitialiseFunction*                         _initialise = _lucSwarmRGBColourViewer_Initialise;
	Stg_Component_ExecuteFunction*                               _execute = _lucSwarmRGBColourViewer_Execute;
	Stg_Component_DestroyFunction*                               _destroy = _lucSwarmRGBColourViewer_Destroy;
	lucDrawingObject_SetupFunction*                                _setup = _lucSwarmRGBColourViewer_Setup;
	lucDrawingObject_DrawFunction*                                  _draw = _lucSwarmRGBColourViewer_Draw;
	lucDrawingObject_CleanUpFunction*                            _cleanUp = _lucSwarmRGBColourViewer_CleanUp;
	lucOpenGLDrawingObject_BuildDisplayListFunction*    _buildDisplayList = _lucSwarmRGBColourViewer_BuildDisplayList;
	lucSwarmViewerBase_PlotParticleFunction*                _plotParticle = _lucSwarmViewer_PlotParticle;
	lucSwarmViewerBase_SetParticleColourFunction*      _setParticleColour = _lucSwarmRGBColourViewer_SetParticleColour;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

	return (void*) _lucSwarmRGBColourViewer_New(  LUCSWARMRGBCOLOURVIEWER_PASSARGS  );
}

void _lucSwarmRGBColourViewer_AssignFromXML( void* drawingObject, Stg_ComponentFactory* cf, void* data ){
	lucSwarmRGBColourViewer* self = (lucSwarmRGBColourViewer*)drawingObject;
	Name                     colourRedVariableName;
	Name                     colourGreenVariableName;
	Name                     colourBlueVariableName;

	/* Construct Parent */
	_lucSwarmViewer_AssignFromXML( self, cf, data );

	colourRedVariableName   = Stg_ComponentFactory_GetString( cf, self->name, (Dictionary_Entry_Key)"ColourRedVariable", ""  );
	colourGreenVariableName = Stg_ComponentFactory_GetString( cf, self->name, (Dictionary_Entry_Key)"ColourGreenVariable", ""  );
	colourBlueVariableName  = Stg_ComponentFactory_GetString( cf, self->name, (Dictionary_Entry_Key)"ColourBlueVariable", ""  );
	
	_lucSwarmRGBColourViewer_Init( 
			self, 
			colourRedVariableName, 
			colourGreenVariableName, 
			colourBlueVariableName );
}

void _lucSwarmRGBColourViewer_Build( void* drawingObject, void* data ) {}

void _lucSwarmRGBColourViewer_Initialise( void* drawingObject, void* data ) {
	lucSwarmRGBColourViewer*	self                   = (lucSwarmRGBColourViewer*)drawingObject;
	SwarmVariable_Register*  	swarmVariable_Register = self->swarm->swarmVariable_Register;
	Stream*                  	errorStr               = Journal_Register( Error_Type, (Name)self->type  );

	_lucSwarmViewer_Initialise( self, data );

	if ( 0 != strcmp( self->colourRedVariableName, "" ) ) {
		self->colourRedVariable  = SwarmVariable_Register_GetByName( swarmVariable_Register, self->colourRedVariableName );
		Journal_Firewall( self->colourRedVariable != NULL, errorStr,
					"Error - for gLucifer drawing object \"%s\" - in %s(): Colour Variable name given was \"%s\", "
					"but no corresponding SwarmVariable found in the register for swarm \"%s\".\n",
					self->name, __func__, self->colourRedVariableName, self->swarm->name );
		
		Stg_Component_Build( self->colourRedVariable, data, False );
		Stg_Component_Initialise( self->colourRedVariable, data, False );

	}

	if ( 0 != strcmp( self->colourGreenVariableName, "" ) ) {
		self->colourGreenVariable  = SwarmVariable_Register_GetByName( swarmVariable_Register, self->colourGreenVariableName );
		Journal_Firewall( self->colourGreenVariable != NULL, errorStr,
					"Error - for gLucifer drawing object \"%s\" - in %s(): Colour Variable name given was \"%s\", "
					"but no corresponding SwarmVariable found in the register for swarm \"%s\".\n",
					self->name, __func__, self->colourGreenVariableName, self->swarm->name );
	
		Stg_Component_Build( self->colourGreenVariable, data, False );
		Stg_Component_Initialise( self->colourGreenVariable, data, False );
	}
	
	if ( 0 != strcmp( self->colourBlueVariableName, "" ) ) {
		self->colourBlueVariable  = SwarmVariable_Register_GetByName( swarmVariable_Register, self->colourBlueVariableName );
		Journal_Firewall( self->colourBlueVariable != NULL, errorStr,
					"Error - for gLucifer drawing object \"%s\" - in %s(): Colour Variable name given was \"%s\", "
					"but no corresponding SwarmVariable found in the register for swarm \"%s\".\n",
					self->name, __func__, self->colourBlueVariableName, self->swarm->name );
		
		Stg_Component_Build( self->colourBlueVariable, data, False );
		Stg_Component_Initialise( self->colourBlueVariable, data, False );
	}
}


void _lucSwarmRGBColourViewer_Execute( void* drawingObject, void* data ) {}
void _lucSwarmRGBColourViewer_Destroy( void* drawingObject, void* data ) {}

void _lucSwarmRGBColourViewer_Setup( void* drawingObject, void* _context ) {
	lucSwarmRGBColourViewer*          self                = (lucSwarmRGBColourViewer*)drawingObject;
	
	lucSwarmRGBColourViewer_UpdateVariables( self );
		
	_lucSwarmViewer_Setup( self, _context );
}

void _lucSwarmRGBColourViewer_Draw( void* drawingObject, lucWindow* window, lucViewportInfo* viewportInfo, void* _context ) {
	_lucSwarmViewer_Draw( drawingObject, window, viewportInfo, _context );
}


void _lucSwarmRGBColourViewer_CleanUp( void* drawingObject, void* context ) {
	_lucSwarmViewer_CleanUp( drawingObject, context );
}

void _lucSwarmRGBColourViewer_BuildDisplayList( void* drawingObject, void* _context ) {
	lucSwarmRGBColourViewer*          self                = (lucSwarmRGBColourViewer*)drawingObject;
	
	_lucSwarmViewer_BuildDisplayList( self, _context );
}

void _lucSwarmRGBColourViewer_SetParticleColour( void* drawingObject, void* _context, Particle_Index lParticle_I ) {
	lucSwarmRGBColourViewer* self                  = (lucSwarmRGBColourViewer*)drawingObject;
	double                   colourValueRed        = 0.0;
	double                   colourValueGreen      = 0.0;
	double                   colourValueBlue       = 0.0;
	double                   opacity               = 0.0;
	lucColour                colour;

	/* Copy the default colour */
	memcpy( &colour, &self->colour, sizeof( lucColour ) );

	/* Get Red colour */
	if ( self->colourRedVariable ) {
		SwarmVariable_ValueAt( self->colourRedVariable, lParticle_I, &colourValueRed );
		colour.red = (float) colourValueRed;
		/* Other way to do it... */
		/* colour.red = Variable_GetValueFloat(colourRedVariable->variable, lParticle_I);*/
	}
		
	/* Get Green colour */
	if ( self->colourGreenVariable ){
		SwarmVariable_ValueAt( self->colourGreenVariable, lParticle_I, &colourValueGreen );
		colour.green = (float) colourValueGreen;
	}
	
	/* Get Blue colour */
	if ( self->colourBlueVariable ){
		SwarmVariable_ValueAt( self->colourBlueVariable, lParticle_I, &colourValueBlue );
		colour.blue = (float) colourValueBlue;
	}
	
	/* Get Opacity Value */
	if ( self->opacityVariable ){
		SwarmVariable_ValueAt( self->opacityVariable, lParticle_I, &opacity );
		colour.opacity = (float)opacity;	
	}

	lucColour_SetOpenGLColour( &colour );
}




void lucSwarmRGBColourViewer_UpdateVariables( void* drawingObject ) {
	lucSwarmRGBColourViewer*          self                = (lucSwarmRGBColourViewer*)drawingObject;

	lucSwarmViewerBase_UpdateVariables( drawingObject ) ;

	if ( self->colourRedVariable && self->colourRedVariable->variable ) {
		Variable_Update( self->colourRedVariable->variable );
	}	
	if ( self->colourGreenVariable && self->colourGreenVariable->variable ) {
		Variable_Update( self->colourGreenVariable->variable );
	}
	if ( self->colourBlueVariable && self->colourBlueVariable->variable ) {
		Variable_Update( self->colourBlueVariable->variable );
	}
}


