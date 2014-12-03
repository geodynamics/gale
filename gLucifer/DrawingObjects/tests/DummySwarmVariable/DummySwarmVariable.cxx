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

#ifdef HAVE_PYTHON
#include <Python.h>
#endif

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <glucifer/Base/Base.h>
#include <glucifer/Windowing/Windowing.h>
#include <glucifer/RenderingEngines/RenderingEngines.h>
#include <glucifer/DrawingObjects/DrawingObjects.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

const Type DummySwarmVariable_Type = "DummySwarmVariable";
			
void DummySwarmVariable_ValueAt( void* swarmVariable, Particle_Index lParticle_I, double* value ) {
	SwarmVariable*    self = (SwarmVariable*) swarmVariable;
	GlobalParticle* particle = (GlobalParticle*)Swarm_ParticleAt( self->swarm, lParticle_I );
	Dof_Index dof_I =0;

        /*	if the dof count is not greater than the number of dimensions, we can just use the coordinate of the particle */				
        for ( dof_I = 0 ; dof_I < self->dofCount ; dof_I++ ){
	    if ( self->dofCount <= self->dim ){
                if( strcmp(self->name, "normalVariable") == 0 ){
			value[ dof_I ] = particle->coord[ dof_I ];
	        }
	        else if( strcmp(self->name, "vectorVariable") == 0 ){
	            value[ dof_I ] = 0 ;
                    if(dof_I == 0) 
			value[ dof_I ] = particle->coord[ dof_I ];
                }
	        else if( strcmp(self->name, "lengthVariable") == 0 )
	            *value = 0.2;
                else 
	            value[ dof_I ] = particle->coord[ dof_I ];
	    }    
            else
                    abort();  /*worry about this later */
        }
}

double DummySwarmVariable_GetMinGlobalSwarmMagnitude( void* SwarmVariable ) {
	return  -1.0;
}
double DummySwarmVariable_GetMaxGlobalSwarmMagnitude( void* SwarmVariable ) {
	return  1.0;
}

void DummySwarmVariable_Init( SwarmVariable* self, Swarm* swarm, Variable* variable, Index dofCount ) {
	/* Add ourselves to the register for later retrieval by clients */
	self->isConstructed = True;

	self->swarm                  = swarm;
	self->variable               = variable;
	self->dofCount               = dofCount;
	self->swarmVariable_Register = swarm->swarmVariable_Register;
	self->dim                    = swarm->dim;
	
	if ( swarm->swarmVariable_Register != NULL )	
		SwarmVariable_Register_Add( swarm->swarmVariable_Register, self );
}


void _DummySwarmVariable_AssignFromXML( void* swarmVariable, Stg_ComponentFactory* cf, void* data ) {
	SwarmVariable*	        self         = (SwarmVariable*)swarmVariable;
	Swarm*                  swarm;
	Variable*               variable;
	Index                   dofCount;

	swarm = Stg_ComponentFactory_ConstructByName( cf, (Name)"swarm", Swarm, True, data  ) ;
	variable = Stg_ComponentFactory_ConstructByName( cf, (Name)"variable", Variable, False, data  ) ;

	/* Check if this component has it's own component dictionary */
	if ( Dictionary_GetDictionary( cf->componentDict, self->name ) ) {
		dofCount = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, (Dictionary_Entry_Key)"dofCount", 0  );
	}
	/* if it doesn't, then just get the value from the root dictionary */
	else {
		dofCount = Stg_ComponentFactory_GetRootDictUnsignedInt( cf, (Dictionary_Entry_Key)"dofCount", 0  );
	}

 	DummySwarmVariable_Init( self, swarm, variable, dofCount );
}

void* _DummySwarmVariable_DefaultNew( Name name_renamed ) {
	/* Variables set in this function */
	SizeT                                                 _sizeOfSelf = sizeof(SwarmVariable);
	Type                                                         type = DummySwarmVariable_Type;
	Stg_Class_DeleteFunction*                                 _delete = _SwarmVariable_Delete;
	Stg_Class_PrintFunction*                                   _print = _SwarmVariable_Print;
	Stg_Class_CopyFunction*                                     _copy = _SwarmVariable_Copy;
	Stg_Component_DefaultConstructorFunction*     _defaultConstructor = _DummySwarmVariable_DefaultNew;
	Stg_Component_ConstructFunction*                       _construct = _DummySwarmVariable_AssignFromXML;
	Stg_Component_BuildFunction*                               _build = _SwarmVariable_Build;
	Stg_Component_InitialiseFunction*                     _initialise = _SwarmVariable_Initialise;
	Stg_Component_ExecuteFunction*                           _execute = _SwarmVariable_Execute;
	Stg_Component_DestroyFunction*                           _destroy = _SwarmVariable_Destroy;
	Name                                                         name = DummySwarmVariable_ValueAt;
	AllocationType                                 nameAllocationType = DummySwarmVariable_GetMinGlobalSwarmMagnitude;
	SwarmVariable_ValueAtFunction*                           _valueAt = DummySwarmVariable_GetMaxGlobalSwarmMagnitude;
	SwarmVariable_GetGlobalValueFunction*      _getMinGlobalMagnitude = name;

	return _SwarmVariable_New(  SWARMVARIABLE_PASSARGS  );
}

Index DummySwarmVariable_Register( PluginsManager* pluginsManager ) {
	RegisterParent( DummySwarmVariable_Type, SwarmVariable_Type );
	return PluginsManager_Submit( pluginsManager, DummySwarmVariable_Type, (Name)"0", _DummySwarmVariable_DefaultNew  );
}


