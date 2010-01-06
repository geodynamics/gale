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
** $Id: Light.c 510 2006-02-17 04:33:32Z RobertTurnbull $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include "types.h"
#include "Light.h"
#include "ViewportInfo.h"
#include "Init.h"

#include <string.h>
#include <assert.h>
#include <gl.h>
#include <glu.h>

const Type lucLight_Type = "lucLight";

/* The position defaults have been chosen for a directional light source (hence posW=0)
   shining onto the left,top,front corner of a 1x1x1 box, currently the most commonly
   used geometry. -- PatrickSunter, 8 Jun 2006 /
const double LUC_LIGHT_DEFAULT_POS_X = 1.0;
const double LUC_LIGHT_DEFAULT_POS_Y = -2.0;
const double LUC_LIGHT_DEFAULT_POS_Z = -2.0;
const double LUC_LIGHT_DEFAULT_POS_W = 0.0;*/

/* Default light at eye pos shining in all directions without attenuation, ie: sunlight from behind viewer 
 * Setting the light position should be done relative to the eye rather than model, thus if you want the scene lit more from above, 
 * simply increase the y component rather than calculating absolute world coordinates for such a light.
 * this allows the lighting to move with the camera and keep the scene lit in the same way
 * If absolute light coords required in future a flag can easily be implemented to do so */
const double LUC_LIGHT_DEFAULT_POS_X = 0.0;
const double LUC_LIGHT_DEFAULT_POS_Y = 0.0;
const double LUC_LIGHT_DEFAULT_POS_Z = 0.0;
const double LUC_LIGHT_DEFAULT_POS_W = 1.0;

lucLight* _lucLight_New(  LUCLIGHT_DEFARGS  )
{
	lucLight*    self;

	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( _sizeOfSelf >= sizeof(lucLight) );
	/* The following terms are parameters that have been passed into this function but are being set before being passed onto the parent */
	/* This means that any values of these parameters that are passed into this function are not passed onto the parent function
	   and so should be set to ZERO in any children of this class. */
	nameAllocationType = NON_GLOBAL;

	self = (lucLight*) _Stg_Component_New(  STG_COMPONENT_PASSARGS  );
	
	return self;
}

void lucLight_Init(		
		lucLight*                                         self,
		Light_Index index,
		int model,
		int material,
		float position[4],
		float lmodel_ambient[4],
		float spotCutOff,
		float spotDirection[3])
{
	self->index= index;
	self->model         = model;
	self->material          = material;
	self->position[0] = position[0];
	self->position[1] = position[1];
	self->position[2] = position[2];
	self->position[3] = position[3];
	self->lmodel_ambient[0] = lmodel_ambient[0];
	self->lmodel_ambient[1] = lmodel_ambient[1];
	self->lmodel_ambient[2] = lmodel_ambient[2];
	self->lmodel_ambient[3] = lmodel_ambient[3];
	self->spotCutOff = spotCutOff;
	self->spotDirection[0] = spotDirection[0];
	self->spotDirection[1] = spotDirection[1];
	self->spotDirection[2] = spotDirection[2];
	

	
}

lucLight* lucLight_New( 
		Name                                               name,
		Light_Index 				index,
		int                                   model,
		int                                   material,
		float                                 position[4],
		float                                 lmodel_ambient[4],
		float                                 spotCutOff,
		float                                 spotDirection[3]
)
{
	lucLight* self = (lucLight*) _lucLight_DefaultNew( name );

	lucLight_Init( self, index, model, material, position, lmodel_ambient, spotCutOff, spotDirection);

	return self;
}

void _lucLight_Delete( void* light ) {
	lucLight* self        = light;
	
	/*if ( self->originalLight != NULL );
		Stg_Class_Delete( self->originalLight );*/
	
	_Stg_Component_Delete( self );
}

void _lucLight_Print( void* light, Stream* stream ) {
	lucLight* self        = light;
	
	Journal_Printf( stream, "lucLight: %s\n", self->name );

	Stream_Indent( stream );

	/* Print Parent */
	_Stg_Component_Print( self, stream );

	Stream_UnIndent( stream );
}

void* _lucLight_Copy( void* light, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	lucLight* self        = light;
	lucLight* newLight;

	newLight = _Stg_Component_Copy( self, dest, deep, nameExt, ptrMap );

	newLight->index    = self->index;
	newLight->model    = self->model;
	newLight->material = self->material;	
	

	/*TODO*/
       	newLight->position[0] = self->position[0];
	newLight->position[1] = self-> position[1];
	newLight->position[2] = self->position[2];
	newLight->position[3] = self->position[3];
	newLight->lmodel_ambient[0] = self->lmodel_ambient[0];
	newLight->lmodel_ambient[1] = self->lmodel_ambient[1];
	newLight->lmodel_ambient[2] = self->lmodel_ambient[2];
	newLight->lmodel_ambient[3] = self->lmodel_ambient[3];
	newLight->spotCutOff = self->spotCutOff;
	newLight->spotDirection[0] = self->spotDirection[0];
	newLight->spotDirection[1] = self->spotDirection[1];
	newLight->spotDirection[2] = self->spotDirection[2];

	
	return (void*) newLight;
}

void* _lucLight_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                              _sizeOfSelf = sizeof( lucLight );
	Type                                                      type = lucLight_Type;
	Stg_Class_DeleteFunction*                              _delete = _lucLight_Delete;
	Stg_Class_PrintFunction*                                _print = _lucLight_Print;
	Stg_Class_CopyFunction*                                  _copy = _lucLight_Copy;
	Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _lucLight_DefaultNew;
	Stg_Component_ConstructFunction*                    _construct = _lucLight_AssignFromXML;
	Stg_Component_BuildFunction*                            _build = _lucLight_Build;
	Stg_Component_InitialiseFunction*                  _initialise = _lucLight_Initialise;
	Stg_Component_ExecuteFunction*                        _execute = _lucLight_Execute;
	Stg_Component_DestroyFunction*                        _destroy = _lucLight_Destroy;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

	return _lucLight_New(  LUCLIGHT_PASSARGS  );
}

void _lucLight_AssignFromXML( void* light, Stg_ComponentFactory* cf, void* data ) {
	lucLight*             	self               = (lucLight*) light;
	Light_Index             index;
	int 			model;
	int 			material;
	float 	        	position[4];
	float                   spotCutOff;
	float                   spotDirection[3];
	Name                    modelName;
	Name                    materialName;
	float                   lmodel_ambient[4]; 
	
	self->context = Stg_ComponentFactory_ConstructByKey( cf, self->name, "Context", AbstractContext, False, data );
	if( !self->context ) 
		self->context = Stg_ComponentFactory_ConstructByName( cf, "context", AbstractContext, True, data );

	glEnable(GL_LIGHTING);

	/* Spot values */
	spotCutOff = Stg_ComponentFactory_GetDouble( cf, self->name, "spotCutOff", 180.0 );
	spotDirection[1] = Stg_ComponentFactory_GetDouble( cf, self->name, "spotDirectionX", 0.0 );
	spotDirection[2] = Stg_ComponentFactory_GetDouble( cf, self->name, "spotDirectionY", 0.0 );
	spotDirection[3] = Stg_ComponentFactory_GetDouble( cf, self->name, "spotDirectionZ", -1.0 );

        
	/* Ambient values*/
	lmodel_ambient[0] = Stg_ComponentFactory_GetDouble( cf, self->name, "ambR", 0.2 );
	lmodel_ambient[1] = Stg_ComponentFactory_GetDouble( cf, self->name, "ambG", 0.2 );
	lmodel_ambient[2] = Stg_ComponentFactory_GetDouble( cf, self->name, "ambB", 0.2 );
	lmodel_ambient[3] = Stg_ComponentFactory_GetDouble( cf, self->name, "ambA", 1.0 );

	position[0]  = Stg_ComponentFactory_GetDouble( cf, self->name, "posX", LUC_LIGHT_DEFAULT_POS_X );
	position[1]  = Stg_ComponentFactory_GetDouble( cf, self->name, "posY", LUC_LIGHT_DEFAULT_POS_Y );
	position[2]  = Stg_ComponentFactory_GetDouble( cf, self->name, "posZ", LUC_LIGHT_DEFAULT_POS_Z );
        position[3]  = Stg_ComponentFactory_GetDouble( cf, self->name, "posW", LUC_LIGHT_DEFAULT_POS_W );

	/*HACK - Got to retrieve */
	materialName = Stg_ComponentFactory_GetString( cf, self->name, "material", "lucMono" );


	lucLight_Init( self, index, model, material, position, lmodel_ambient, spotCutOff, spotDirection);
	
}

void _lucLight_Build( void* light, void* data ) { }
void _lucLight_Initialise( void* light, void* data ) { }
void _lucLight_Execute( void* light, void* data ) { }
void _lucLight_Destroy( void* light, void* data ) { }



void lucLight_Pickle( void* light, Stream* stream ) {
	lucLight* self                = light;
	       
	Journal_Printf( stream, "<struct name=\"%s\">\n", self->name);
	Stream_Indent( stream );

	Journal_Printf( stream, "<param name=\"Type\">%s</param>\n", self->type );
	Journal_Printf( stream, "<param name=\"posX\">%.5g</param>\n", self->position[ 0 ] );
	Journal_Printf( stream, "<param name=\"posY\">%.5g</param>\n", self->position[ 1 ] );
	Journal_Printf( stream, "<param name=\"posZ\">%.5g</param>\n", self->position[ 2 ] );	
	Journal_Printf( stream, "<param name=\"posW\">%.5g</param>\n", self->position[ 3 ] );

	Journal_Printf( stream, "<param name=\"ambX\">%.5g</param>\n", self->lmodel_ambient[ 0 ] );
	Journal_Printf( stream, "<param name=\"ambY\">%.5g</param>\n", self->lmodel_ambient[ 1 ]  );
	Journal_Printf( stream, "<param name=\"ambZ\">%.5g</param>\n", self->lmodel_ambient[ 2 ]  );
	Journal_Printf( stream, "<param name=\"ambW\">%.5g</param>\n", self->lmodel_ambient[ 3 ]  );

	Journal_Printf( stream, "<param name=\"model\">%s</param>\n", 
			self->model == GL_LIGHT_MODEL_LOCAL_VIEWER ? "Local" :
			self->model == GL_LIGHT_MODEL_AMBIENT ? "Ambient" : 
	                self->model == GL_LIGHT_MODEL_TWO_SIDE ? "TwoSide" : "" );

	Journal_Printf( stream, "<param name=\"spotCutOff\">%.5g</param>\n", self->spotCutOff );
	Journal_Printf( stream, "<param name=\"spotDirectionX\">%.5g</param>\n", self->spotDirection[ 0 ]  );
	Journal_Printf( stream, "<param name=\"spotDirectionY\">%.5g</param>\n", self->spotDirection[ 1 ]  );
	Journal_Printf( stream, "<param name=\"spotDirectionZ\">%.5g</param>\n", self->spotDirection[ 2 ]  );



       	Stream_UnIndent( stream );
	Journal_Printf( stream, "</struct>\n");
}

/* functions to change the lights paramters - unused as yet, only called from light interactions which don't work */
void lucLight_Position( void * light, int lightIndex, float posX, float posY, float posZ, float posW) {
	lucLight*             	self               = (lucLight*) light;
	Light_Index             index              = (Light_Index) lightIndex;

	/* Sets the position of the light index = index */	
	self->position[0]  += posX;
	self->position[1]  += posY;
	self->position[2]  += posZ;
   self->position[3]  += posW;

   /* Light position now relative to eye, not model! */
   glPushMatrix();
   glLoadIdentity();
	glLightfv(GL_LIGHT0 + lightIndex, GL_POSITION, self->position);
   glPopMatrix();

	self->needsToDraw = True;
	
}
void lucLight_Material( int material) {
	
	
}
void lucLight_SetNeedsToDraw( void * light ){
	lucLight* self = (lucLight*) light;
	self->needsToDraw = True;
}


