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

#include "Light_Register.h"
#include "Light.h"
#include <gl.h>
#include <glu.h>

const Type lucLight_Register_Type = "lucLight_Register";

lucLight_Register*	lucLight_Register_New( void ) {
	/* Variables set in this function */
	SizeT                      _sizeOfSelf = sizeof(lucLight_Register);
	Type                              type = lucLight_Register_Type;
	Stg_Class_DeleteFunction*      _delete = _NamedObject_Register_Delete;
	Stg_Class_PrintFunction*        _print = _NamedObject_Register_Print;
	Stg_Class_CopyFunction*          _copy = _NamedObject_Register_Copy;

	lucLight_Register* self;
	
	self = (lucLight_Register*) _NamedObject_Register_New(  NAMEDOBJECT_REGISTER_PASSARGS  );
	self->currentLightIndex = 0;

	return self;
}

void    lucLight_Register_EnableAll( void * lightRegister ) {
	lucLight_Register* self = (lucLight_Register*) lightRegister;
	
	lucLight* light;
	Light_Index lightCount = 0;
	Light_Index light_I = 0;
	
	/* Enabling the lights */
	glEnable(GL_LIGHTING);
	lightCount = lucLight_Register_GetCount( self );
	
	for (light_I = 0; light_I < lightCount; light_I++){
            light = lucLight_Register_GetByIndex(self, light_I);	    
	    glLightfv(GL_LIGHT0 + light_I, GL_POSITION, light->position);
	    glLightf(GL_LIGHT0 + light_I, GL_SPOT_CUTOFF, light->spotCutOff);
	    glLightfv(GL_LIGHT0 + light_I, GL_SPOT_DIRECTION, light->spotDirection);	

            glEnable(GL_LIGHT0 + light_I);	    
	}
	
	if(lightCount == 0) 	/* No light is listed in the xml file, but we still need to enable the default light */
	   glEnable(GL_LIGHT0);

}

Light_Index   lucLight_Register_GetCurrentLightIndex( void * lightRegister ) {
	lucLight_Register* self = (lucLight_Register*) lightRegister;
	return self->currentLightIndex;
}

void lucLight_Register_SetCurrentLightIndex( void * lightRegister, Light_Index index ) {
	lucLight_Register* self = (lucLight_Register*) lightRegister;
	self->currentLightIndex = index;
}
void lucLight_Register_ChangeCurrentLightIndex( void * lightRegister ) {
	lucLight_Register* self = (lucLight_Register*) lightRegister;
	Light_Index lightCount   = lucLight_Register_GetCount( self );
	self->currentLightIndex ++;
	if (self->currentLightIndex == lightCount)  self->currentLightIndex = 0;

}


