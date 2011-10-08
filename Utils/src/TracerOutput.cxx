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
*%		Louis Moresi - Louis.Moresi@sci.monash.edu.au
*%
** Contributors:
*+		Robert Turnbull
*+		Vincent Lemiale
*+		Louis Moresi
*+		David May
*+		David Stegman
*+		Mirko Velic
*+		Patrick Sunter
*+		Julian Giordani
*+
** $Id: TracerOutput.c 610 2007-10-11 08:09:29Z SteveQuenette $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>

#include "types.h"
#include "TracerOutput.h"

#include <assert.h>
#include <string.h>

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type TracerOutput_Type = "TracerOutput";

/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/
TracerOutput* TracerOutput_New(Name name,
                               TimeIntegrator* timeIntegrator,
                               FeVariable* velocityField )
{
  TracerOutput* self = (TracerOutput*) _TracerOutput_DefaultNew( name );

  /* 	TracerOutput_InitAll */
  abort();

  return self;
}

TracerOutput* _TracerOutput_New(  TRACEROUTPUT_DEFARGS  )
{
  TracerOutput* self;
	
  /* Call private constructor of parent - this will set virtual
     functions of parent and continue up the hierarchy tree. At the
     beginning of the tree it will allocate memory of the size of
     object and initialise all the memory to zero. */
  assert( _sizeOfSelf >= sizeof(TracerOutput) );
  self = (TracerOutput*)_SwarmOutput_New(  SWARMOUTPUT_PASSARGS  );

  /* General info */

  /* Virtual Info */
	
  return self;
}

void _TracerOutput_Init(void*                                 swarmOutput,
                        FeVariable*                           pressureField,
                        FeVariable**                           fields,
                        unsigned int num_fields)
{
  TracerOutput*   self                = (TracerOutput*)swarmOutput;

  self->pressureField      = pressureField;
  self->num_fields=num_fields;
  self->fields=Memory_Alloc_Array( FeVariable*, num_fields, "fields" );
  memcpy( self->fields, fields, num_fields * sizeof(FeVariable*) );
}


/*------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _TracerOutput_Delete( void* swarmOutput ) {
  TracerOutput* self = (TracerOutput*)swarmOutput;

  Memory_Free(self->fields);
  _SwarmOutput_Delete( self );
}


void _TracerOutput_Print( void* swarmOutput, Stream* stream ) {
  TracerOutput* self = (TracerOutput*)swarmOutput;
	
  /* Print parent */
  _SwarmOutput_Print( self, stream );
}

void* _TracerOutput_Copy( const void* swarmOutput, void* dest, Bool deep,
                          Name nameExt, PtrMap* ptrMap ) {
  TracerOutput*	self = (TracerOutput*)swarmOutput;
  TracerOutput*	newTracerOutput;
	
  newTracerOutput = (TracerOutput*)_SwarmOutput_Copy( self, dest, deep, nameExt, ptrMap );
  memcpy(newTracerOutput->fields, self->fields,
         self->num_fields * sizeof(FeVariable*) );
  newTracerOutput->num_fields = self->num_fields;
	
  return (void*)newTracerOutput;
}

void* _TracerOutput_DefaultNew( Name name ) {
  /* Variables set in this function */
  SizeT _sizeOfSelf = sizeof(TracerOutput);
  Type type = TracerOutput_Type;
  Stg_Class_DeleteFunction* _delete = _TracerOutput_Delete;
  Stg_Class_PrintFunction* _print = _TracerOutput_Print;
  Stg_Class_CopyFunction* _copy = _TracerOutput_Copy;
  Stg_Component_DefaultConstructorFunction* _defaultConstructor
    = _TracerOutput_DefaultNew;
  Stg_Component_ConstructFunction* _construct = _TracerOutput_AssignFromXML;
  Stg_Component_BuildFunction* _build = _TracerOutput_Build;
  Stg_Component_InitialiseFunction* _initialise = _TracerOutput_Initialise;
  Stg_Component_ExecuteFunction* _execute = _TracerOutput_Execute;
  Stg_Component_DestroyFunction* _destroy = _TracerOutput_Destroy;
  SwarmOutput_PrintHeaderFunction* _printHeader = _TracerOutput_PrintHeader;
  SwarmOutput_PrintDataFunction* _printData = _TracerOutput_PrintData;

  /* Variables that are set to ZERO are variables that will be set
     either by the current _New function or another parent _New
     function further up the hierachy */
  AllocationType nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

  return (void*) _TracerOutput_New(  TRACEROUTPUT_PASSARGS  );
}


void _TracerOutput_AssignFromXML( void* swarmOutput,
                                  Stg_ComponentFactory* cf, void* data ) {
  TracerOutput*  self          = (TracerOutput*) swarmOutput;
  FeVariable*                 pressureField;
  FeVariable**                 fields;
  unsigned int num_fields;

  _SwarmOutput_AssignFromXML( self, cf, data );

  pressureField =
    Stg_ComponentFactory_ConstructByKey(cf, self->name,
                                        (Dictionary_Entry_Key)"PressureField",
                                        FeVariable, False, data) ;
  
  fields=Stg_ComponentFactory_ConstructByList(cf,self->name,
                                              (Dictionary_Entry_Key)"Fields",
                                              Stg_ComponentFactory_Unlimited,
                                              FeVariable,False,&num_fields,data);

  _TracerOutput_Init(self,pressureField,fields,num_fields);

  Memory_Free(fields);
}

void _TracerOutput_Build( void* swarmOutput, void* data ) {
  int i;
  TracerOutput*	self = (TracerOutput*) swarmOutput;

  if(self->pressureField)
    Stg_Component_Build( self->pressureField, data, False );
  for(i=0;i<self->num_fields;++i)
    Stg_Component_Build( self->fields[i], data, False );

  _SwarmOutput_Build( self, data );
}
void _TracerOutput_Initialise( void* swarmOutput, void* data ) {
  int i;
  TracerOutput*	self = (TracerOutput*) swarmOutput;

  if(self->pressureField)
    Stg_Component_Initialise( self->pressureField, data, False );
  for(i=0;i<self->num_fields;++i)
    Stg_Component_Initialise( self->fields[i], data, False );
	
  _SwarmOutput_Initialise( self, data );
}
void _TracerOutput_Execute( void* swarmOutput, void* data ) {
  TracerOutput*	self = (TracerOutput*)swarmOutput;
	
  _SwarmOutput_Execute( self, data );
}
void _TracerOutput_Destroy( void* swarmOutput, void* data ) {
  int i;
  TracerOutput*	self = (TracerOutput*)swarmOutput;

  if(self->pressureField)
    Stg_Component_Destroy( self->pressureField, data, False );
  for(i=0;i<self->num_fields;++i)
    Stg_Component_Destroy( self->fields[i], data, False );
	
  _SwarmOutput_Destroy( self, data );
}

void _TracerOutput_PrintHeader( void* swarmOutput, Stream* stream,
                                Particle_Index lParticle_I, void* context ){
  int i;
  char name[32];
  TracerOutput*	self = (TracerOutput*)swarmOutput;
	
  Stream_Enable(stream,True);
  _SwarmOutput_PrintHeader( self, stream, lParticle_I, context );
	
  SwarmOutput_PrintString( self, stream, "Pressure" );
  for(i=0;i<self->num_fields;++i)
    {
      sprintf(name,"Field%d",i);
      SwarmOutput_PrintString( self, stream, name );
    }
}

void _TracerOutput_PrintData( void* swarmOutput, Stream* stream,
                              Particle_Index lParticle_I, void* context ){
  TracerOutput*	self = (TracerOutput*)swarmOutput;
  FiniteElementContext* fe_context = (FiniteElementContext*) context;
  Swarm* swarm = self->swarm;
  GlobalParticle* particle = (GlobalParticle*)Swarm_ParticleAt(swarm,
                                                               lParticle_I);
  double* coord = particle->coord;
  HydrostaticTerm *hydrostaticTerm;
  double pressure, field;
  int i;

  hydrostaticTerm =
    (HydrostaticTerm*)LiveComponentRegister_Get(fe_context->CF->LCRegister,
                                                "hydrostaticTerm" );

  Journal_Firewall(swarm->particleLayout->coordSystem == GlobalCoordSystem,
                   Journal_MyStream( Error_Type, self ),
                   "Swarm is not using global coord system! Modify this code to use both systems\n" );

  _SwarmOutput_PrintData( self, stream, lParticle_I, context );

  if(self->pressureField)
    {
      FieldVariable_InterpolateValueAt(self->pressureField,coord,&pressure );
      if(hydrostaticTerm){
        pressure+=HydrostaticTerm_Pressure(hydrostaticTerm,coord);
      }
      SwarmOutput_PrintValue( self, stream, pressure );
    }
  for(i=0;i<self->num_fields;++i)
    {
      FieldVariable_InterpolateValueAt(self->fields[i],coord,&field);
      SwarmOutput_PrintValue( self, stream, field );
    }
}

/*-------------------------------------------------------------------------------------------------------------------------
** Private Functions
*/
/*---------------------------------------------------------------------------------------------------------------------
** Entry Point Hooks
*/

/*-------------------------------------------------------------------------------------------------------------------------
** Public Functions
*/



