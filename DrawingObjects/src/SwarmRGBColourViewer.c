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
lucSwarmRGBColourViewer* _lucSwarmRGBColourViewer_New( 
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
	lucSwarmRGBColourViewer*					self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( sizeOfSelf >= sizeof(lucSwarmRGBColourViewer) );
	self = (lucSwarmRGBColourViewer*) _lucSwarmViewer_New( 
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
			_plotParticle,
			_setParticleColour,
			name );

	return self;
}

void _lucSwarmRGBColourViewer_Init( 
		lucSwarmRGBColourViewer*                                     self,
		Name                                                         colourRedVariableName,
		Name                                                         colourGreenVariableName,
		Name                                                         colourBlueVariableName)
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
	return (void*) _lucSwarmRGBColourViewer_New(
		sizeof(lucSwarmRGBColourViewer),
		lucSwarmRGBColourViewer_Type,
		_lucSwarmRGBColourViewer_Delete,
		_lucSwarmRGBColourViewer_Print,
		NULL,
		_lucSwarmRGBColourViewer_DefaultNew,
		_lucSwarmRGBColourViewer_Construct,
		_lucSwarmRGBColourViewer_Build,
		_lucSwarmRGBColourViewer_Initialise,
		_lucSwarmRGBColourViewer_Execute,
		_lucSwarmRGBColourViewer_Destroy,
		_lucSwarmRGBColourViewer_Setup,
		_lucSwarmRGBColourViewer_Draw,
		_lucSwarmRGBColourViewer_CleanUp,
		_lucSwarmRGBColourViewer_BuildDisplayList,
		_lucSwarmViewer_PlotParticle,
		_lucSwarmRGBColourViewer_SetParticleColour,
		name );
}

void _lucSwarmRGBColourViewer_Construct( void* drawingObject, Stg_ComponentFactory* cf, void* data ){
	lucSwarmRGBColourViewer*         self = (lucSwarmRGBColourViewer*)drawingObject;
	Name                    colourRedVariableName;
	Name                    colourGreenVariableName;
	Name                    colourBlueVariableName;

	/* Construct Parent */
	_lucSwarmViewer_Construct( self, cf, data );

	colourRedVariableName = Stg_ComponentFactory_GetString( cf, self->name, "ColourRedVariable", "" );
	colourGreenVariableName = Stg_ComponentFactory_GetString( cf, self->name, "ColourGreenVariable", "" );
	colourBlueVariableName = Stg_ComponentFactory_GetString( cf, self->name, "ColourBlueVariable", "" );
	
	_lucSwarmRGBColourViewer_Init( 
			self, 
			colourRedVariableName, 
			colourGreenVariableName, 
			colourBlueVariableName );
}

void _lucSwarmRGBColourViewer_Build( void* drawingObject, void* data ) {}

void _lucSwarmRGBColourViewer_Initialise( void* drawingObject, void* data ) {
	lucSwarmRGBColourViewer*          self                   = (lucSwarmRGBColourViewer*)drawingObject;
	SwarmVariable_Register*  swarmVariable_Register = self->swarm->swarmVariable_Register;
	Stream*                  errorStr               = Journal_Register( Error_Type, self->type );

	_lucSwarmViewer_Initialise( self, data );

	if ( 0 != strcmp( self->colourRedVariableName, "" ) ) {
		self->colourRedVariable  = SwarmVariable_Register_GetByName( swarmVariable_Register, self->colourRedVariableName );
		Journal_Firewall( self->colourRedVariable != NULL, errorStr,
					"Error - for gLucifer drawing object \"%s\" - in %s(): colour Variable name given was \"%s\", "
					"but no corresponding SwarmVariable found in the register for swarm \"%s\".\n",
					self->name, __func__, self->colourRedVariableName, self->swarm->name );
		
		Stg_Component_Build( self->colourRedVariable, data, False );
		Stg_Component_Initialise( self->colourRedVariable, data, False );

	}

	if ( 0 != strcmp( self->colourGreenVariableName, "" ) ) {
		self->colourGreenVariable  = SwarmVariable_Register_GetByName( swarmVariable_Register, self->colourGreenVariableName );
		Journal_Firewall( self->colourGreenVariable != NULL, errorStr,
					"Error - for gLucifer drawing object \"%s\" - in %s(): colour Variable name given was \"%s\", "
					"but no corresponding SwarmVariable found in the register for swarm \"%s\".\n",
					self->name, __func__, self->colourGreenVariableName, self->swarm->name );
	
		Stg_Component_Build( self->colourGreenVariable, data, False );
		Stg_Component_Initialise( self->colourGreenVariable, data, False );
	}
	
	if ( 0 != strcmp( self->colourBlueVariableName, "" ) ) {
		self->colourBlueVariable  = SwarmVariable_Register_GetByName( swarmVariable_Register, self->colourBlueVariableName );
		Journal_Firewall( self->colourBlueVariable != NULL, errorStr,
					"Error - for gLucifer drawing object \"%s\" - in %s(): colour Variable name given was \"%s\", "
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

void _lucSwarmRGBColourViewer_SetParticleColour( void* drawingObject, void* _context, Particle_Index lParticle_I, lucColour* colour ) {
	lucSwarmRGBColourViewer*          self                = (lucSwarmRGBColourViewer*)drawingObject;
	double                   colourValueRed        = 0.0;
	double                   colourValueGreen      = 0.0;
	double                   colourValueBlue       = 0.0;

	/* Get red colour */
	if ( self->colourRedVariable){
		SwarmVariable_ValueAt( self->colourRedVariable, lParticle_I, &colourValueRed );
		(self->colour).red = colourValueRed;
		/* Other way to do it... */
		/* colour.red = Variable_GetValueFloat(colourRedVariable->variable, lParticle_I);*/
	}
	/* Get green colour */
	if ( self->colourGreenVariable){
		SwarmVariable_ValueAt( self->colourGreenVariable, lParticle_I, &colourValueGreen );
		(self->colour).green = colourValueGreen;
	}
	
	/* Get blue colour */
	if ( self->colourBlueVariable){
		SwarmVariable_ValueAt( self->colourBlueVariable, lParticle_I, &colourValueBlue );
		(self->colour).blue = colourValueBlue;
	}

	
	lucColourMap_SetOpenGLColourFromRGB_ExplicitOpacity( (self->colour).red, (self->colour).green, (self->colour).blue, (float)(self->colour).opacity );
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
