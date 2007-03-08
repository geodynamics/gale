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
** $Id: Underworld_SwarmOutput.c 358 2006-10-18 06:17:30Z SteveQuenette $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>

#include "types.h"
#include "Underworld_SwarmOutput.h"

#include <assert.h>
#include <string.h>

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type Underworld_SwarmOutput_Type = "Underworld_SwarmOutput";

/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/
Underworld_SwarmOutput* Underworld_SwarmOutput_New(
		Name                                  name,
		TimeIntegrator*                       timeIntegrator,
		FeVariable*                           velocityField )
{
	Underworld_SwarmOutput* self = (Underworld_SwarmOutput*) _Underworld_SwarmOutput_DefaultNew( name );

	/* 	Underworld_SwarmOutput_InitAll */
	abort();

	return self;
}

Underworld_SwarmOutput* _Underworld_SwarmOutput_New(
		SizeT                                              _sizeOfSelf, 
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
		SwarmOutput_PrintHeaderFunction*                   _printHeader,		
		SwarmOutput_PrintDataFunction*                     _printData,		
		Name                                               name )
{
	Underworld_SwarmOutput* self;
	
	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( _sizeOfSelf >= sizeof(Underworld_SwarmOutput) );
	self = (Underworld_SwarmOutput*)_SwarmOutput_New( 
			_sizeOfSelf,
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
			_printHeader,
			_printData,
			name );

	
	/* General info */

	/* Virtual Info */
	
	return self;
}

void _Underworld_SwarmOutput_Init( 
		void*                                 swarmOutput,
		Stg_ComponentFactory*                 cf,
	        void*                                 data )
{
	Underworld_SwarmOutput*       self       = (Underworld_SwarmOutput*)swarmOutput;
	Dictionary*             dictionary;
	Dictionary_Entry_Value* list;
	Index                   listCount, list_I;
	char*                   varName;

	self->infoStream = Journal_Register( Info_Type, "Underworld_SwarmOutput Info" );
	self->errorStream = Journal_Register( Error_Type, "Underworld_SwarmOutput Error" );
  
	dictionary = Dictionary_GetDictionary( cf->componentDict, self->name );
	/* Find how many FeVariables to Display */
	list = Dictionary_Get( dictionary, "FieldVariables" );
	Journal_Firewall( list > 0, self->errorStream, "Error in %s:\nNo list of FieldVariables could be found\n", __FILE__ ); 

	listCount = Dictionary_Entry_Value_GetCount( list );
	if( listCount < 1 ) Journal_Printf( self->infoStream, "Warning in %s: 0 fields have been specified to be used for swarm ouput, is this what you want ???\n");
	self->sizeList = listCount;

	/* Allocate the memory to store pointers to them */
	self->feVariableList = Memory_Alloc_Array( FieldVariable*, listCount, "List of FeVariables to output from" );

	/* Assign pointers */
	for( list_I = 0 ; list_I < listCount ; list_I++ ) {
		varName = Dictionary_Entry_Value_AsString( Dictionary_Entry_Value_GetElement( list, list_I ) );
		self->feVariableList[ list_I ] = 
			Stg_ComponentFactory_ConstructByName( cf, varName, FieldVariable, True, data );
	}
}


/*------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _Underworld_SwarmOutput_Delete( void* swarmOutput ) {
	Underworld_SwarmOutput* self = (Underworld_SwarmOutput*)swarmOutput;

	/* Delete parent */
	_SwarmOutput_Delete( self );
}


void _Underworld_SwarmOutput_Print( void* swarmOutput, Stream* stream ) {
	Underworld_SwarmOutput* self = (Underworld_SwarmOutput*)swarmOutput;
	
	/* Print parent */
	_SwarmOutput_Print( self, stream );
}

void* _Underworld_SwarmOutput_Copy( void* swarmOutput, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	Underworld_SwarmOutput*	self = (Underworld_SwarmOutput*)swarmOutput;
	Underworld_SwarmOutput*	newUnderworld_SwarmOutput;
	
	newUnderworld_SwarmOutput = (Underworld_SwarmOutput*)_SwarmOutput_Copy( self, dest, deep, nameExt, ptrMap );
	
	return (void*)newUnderworld_SwarmOutput;
}

void* _Underworld_SwarmOutput_DefaultNew( Name name ) {
	return (void*) _Underworld_SwarmOutput_New(
			sizeof(Underworld_SwarmOutput),
			Underworld_SwarmOutput_Type,
			_Underworld_SwarmOutput_Delete,
			_Underworld_SwarmOutput_Print,
			_Underworld_SwarmOutput_Copy,
			_Underworld_SwarmOutput_DefaultNew,
			_Underworld_SwarmOutput_Construct,
			_Underworld_SwarmOutput_Build,
			_Underworld_SwarmOutput_Initialise,
			_Underworld_SwarmOutput_Execute,
			_Underworld_SwarmOutput_Destroy,
			_Underworld_SwarmOutput_PrintHeader,
			_Underworld_SwarmOutput_PrintData,
			name );
}


void _Underworld_SwarmOutput_Construct( void* swarmOutput, Stg_ComponentFactory* cf, void* data ) {
	Underworld_SwarmOutput*  self          = (Underworld_SwarmOutput*) swarmOutput;

	_SwarmOutput_Construct( self, cf, data );
	
	_Underworld_SwarmOutput_Init( self, cf, data );

}

void _Underworld_SwarmOutput_Build( void* swarmOutput, void* data ) {
	Underworld_SwarmOutput*	self = (Underworld_SwarmOutput*) swarmOutput;

	_SwarmOutput_Build( self, data );
}
void _Underworld_SwarmOutput_Initialise( void* swarmOutput, void* data ) {
	Underworld_SwarmOutput*	self = (Underworld_SwarmOutput*) swarmOutput;
	
	_SwarmOutput_Initialise( self, data );
}
void _Underworld_SwarmOutput_Execute( void* swarmOutput, void* data ) {
	Underworld_SwarmOutput*	self = (Underworld_SwarmOutput*)swarmOutput;
	
	_SwarmOutput_Execute( self, data );
}
void _Underworld_SwarmOutput_Destroy( void* swarmOutput, void* data ) {
	Underworld_SwarmOutput*	self = (Underworld_SwarmOutput*)swarmOutput;
	
	_SwarmOutput_Destroy( self, data );

	Memory_Free( self->feVariableList );
}

void _Underworld_SwarmOutput_PrintHeader( void* swarmOutput, Stream* stream, Particle_Index lParticle_I, void* context ){
	Underworld_SwarmOutput*	self = (Underworld_SwarmOutput*)swarmOutput;
	Index                   list_I;
	
	_SwarmOutput_PrintHeader( self, stream, lParticle_I, context );
	
	for( list_I = 0; list_I < self->sizeList ; list_I++ ) 
		SwarmOutput_PrintString( self, stream, self->feVariableList[list_I]->name );
	
}

void _Underworld_SwarmOutput_PrintData( void* swarmOutput, Stream* stream, Particle_Index lParticle_I, void* context ){
	Underworld_SwarmOutput*	self                = (Underworld_SwarmOutput*)swarmOutput;
	Swarm*                  swarm               = self->swarm;
	GlobalParticle*         particle            = (GlobalParticle*)Swarm_ParticleAt( swarm, lParticle_I );
	double*                 coord               = particle->coord;
	double                  value[50];  /* Assume no field has more than 50 components */
	InterpolationResult     result;
	Index                   list_I;

	Journal_Firewall(
		swarm->particleLayout->coordSystem == GlobalCoordSystem,
		Journal_MyStream( Error_Type, self ),
		"Swarm is not using global coord system! Modify this code to use both systems\n" );

	_SwarmOutput_PrintData( self, stream, lParticle_I, context );

	for( list_I = 0 ; list_I < self->sizeList ; list_I++ ) {
		result = FieldVariable_InterpolateValueAt( (FieldVariable*)self->feVariableList[list_I], coord, value );
		if( result == OTHER_PROC || result == OUTSIDE_GLOBAL )
			assert( 0 );

		SwarmOutput_PrintTuple( self, stream, value, self->feVariableList[list_I]->fieldComponentCount );
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

