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


/*------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _Underworld_SwarmOutput_Delete( void* uwSwarmOutput ) {
	Underworld_SwarmOutput* self = (Underworld_SwarmOutput*)uwSwarmOutput;

	/* Delete parent */
	_SwarmOutput_Delete( self );
}


void _Underworld_SwarmOutput_Print( void* uwSwarmOutput, Stream* stream ) {
	Underworld_SwarmOutput* self = (Underworld_SwarmOutput*)uwSwarmOutput;
	
	/* Print parent */
	_SwarmOutput_Print( self, stream );
}

void* _Underworld_SwarmOutput_Copy( void* uwSwarmOutput, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	Underworld_SwarmOutput*	self = (Underworld_SwarmOutput*)uwSwarmOutput;
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


void _Underworld_SwarmOutput_Construct( void* uwSwarmOutput, Stg_ComponentFactory* cf, void* data ) {
	Underworld_SwarmOutput*  self          = (Underworld_SwarmOutput*) uwSwarmOutput;

	MaterialPointsSwarm*    materialSwarm;
	Dictionary*             dictionary = Dictionary_GetDictionary( cf->componentDict, self->name );
	Dictionary_Entry_Value* list;
	unsigned int            listCount, feVar_I;
	char*                   varName;

	materialSwarm = (MaterialPointsSwarm*)Stg_ComponentFactory_ConstructByKey( cf, self->name, "Swarm", MaterialPointsSwarm, True, data );
	/* Get all Swarms specified in input file, under swarms */
	list = Dictionary_Get( dictionary, "FeVariables" );
	listCount = Dictionary_Entry_Value_GetCount( list );

	/* Allocate the memory to store pointers to them */
	self->feVariableList = Memory_Alloc_Array( FeVariable*, listCount, "List FeVariables" );

	for( feVar_I = 0 ; feVar_I < listCount ; feVar_I++ ) {
		varName = Dictionary_Entry_Value_AsString( Dictionary_Entry_Value_GetElement( list, feVar_I ) );
		self->feVariableList[ feVar_I ] = Stg_ComponentFactory_ConstructByName( cf, varName, FeVariable, True, data );
	}

	self->materialSwarm = materialSwarm;
	self->sizeList = listCount;
	self->_getFeValuesFunc = _Underworld_SwarmOutput_GetFeVariableValues;
	self->_printFunc = _Underworld_SwarmOutput_PrintStandardFormat;
}

void _Underworld_SwarmOutput_Build( void* uwSwarmOutput, void* data ) {
	Underworld_SwarmOutput*	self = (Underworld_SwarmOutput*) uwSwarmOutput;
	unsigned int feVar_I, feVarNum;
	
	/* Build all StGermain based data structure this component uses */
	feVarNum = self->sizeList;
	for( feVar_I = 0 ; feVar_I < feVarNum; feVar_I++ ) 
		Stg_Component_Build( self->feVariableList[feVar_I], data, False );

	Stg_Component_Build( self->materialSwarm, data, False );
}
void _Underworld_SwarmOutput_Initialise( void* uwSwarmOutput, void* data ) {
	Underworld_SwarmOutput*	self = (Underworld_SwarmOutput*) uwSwarmOutput;
	unsigned int feVar_I, feVarNum;
	
	/* Initialise all StGermain based data structure this component uses */
	feVarNum = self->sizeList;
	for( feVar_I = 0 ; feVar_I < feVarNum; feVar_I++ ) 
		Stg_Component_Initialise( self->feVariableList[feVar_I], data, False );

	Stg_Component_Initialise( self->materialSwarm, data, False );
}

void _Underworld_SwarmOutput_Execute( void* uwSwarmOutput, void* data ) {
	Underworld_SwarmOutput*	self = (Underworld_SwarmOutput*)uwSwarmOutput;
	PICelleratorContext*    context = (PICelleratorContext*)data;
	MaterialPointsSwarm*    materialSwarm;
	FeVariable*             feVar;
	char*                   filename;
	unsigned int            feVar_I, feVarNum, timeStep;

	FILE*              outputFile;
	MPI_Comm	   comm;
	int                myRank;
	int                nProcs;
	MPI_Status         status;
	const int          FINISHED_WRITING_TAG = 100;
	int                confirmation = 0;

	MPI_Comm_size( comm, (int*)&nProcs );
	MPI_Comm_rank( comm, (int*)&myRank );
	
	timeStep = context->timeStep;
	feVarNum = self->sizeList;
	materialSwarm = self->materialSwarm;
	
	/* for each field create file and then write */
	for( feVar_I = 0 ; feVar_I < feVarNum; feVar_I++ ) {
		feVar = self->feVariableList[feVar_I];
		comm = Comm_GetMPIComm( Mesh_GetCommTopology( feVar->feMesh, MT_VERTEX ) );

		/*                                                FieldName             SwarmName                . time  .  dat \0 */
		filename = Memory_Alloc_Array_Unnamed( char, strlen(feVar->name) + strlen(materialSwarm->name) + 1 + 5 + 1 + 3 + 1 );
		sprintf( filename, "%s.%.5u.dat", filename, timeStep );
		/* wait for go-ahead from process ranked lower than me, to avoid competition writing to file */
		if ( myRank != 0 ) {
			MPI_Recv( &confirmation, 1, MPI_INT, myRank - 1, FINISHED_WRITING_TAG, comm, &status );
		}	

		/* if myRank is 0, create file, otherwise append to file */
		(myRank == 0) ?
		       	( outputFile = fopen( filename, "w" ) ) :
		       	( outputFile = fopen( filename, "a" ) ) ;

		self->_getFeValuesFunc( self, feVar, materialSwarm, outputFile );
		fclose( outputFile );
	}
	Memory_Free( filename );
}

void _Underworld_SwarmOutput_Destroy( void* uwSwarmOutput, void* data ) {
	Underworld_SwarmOutput*	self = (Underworld_SwarmOutput*)uwSwarmOutput;
	
	Memory_Free( self->feVariableList );
}

void _Underworld_SwarmOutput_PrintStandardFormat( MaterialPoint* particle, unsigned particleID, double* result, unsigned fieldComponentCount, FILE* outputFile ) {
		unsigned dof_I;
		fprintf( outputFile, "%u ", particleID );
                fprintf( outputFile, "%.15g %.15g %.15g ", particle->coord[0],particle->coord[1], particle->coord[2] ); 
		for ( dof_I = 0; dof_I < fieldComponentCount; dof_I++ ) 
			fprintf( outputFile, "%.15g ", result[dof_I] );
			
		fprintf( outputFile, "\n" );
}
	
void _Underworld_SwarmOutput_GetFeVariableValues(Underworld_SwarmOutput* uwSwarmOutput, FeVariable* feVariable, MaterialPointsSwarm* swarm, FILE* outputFile) {
	Underworld_SwarmOutput*	self = (Underworld_SwarmOutput*)uwSwarmOutput;
	MaterialPoint*          particle;
	unsigned int            localElementCount, lElement_I; 
	unsigned int            cellID, cellParticleCount, cParticle_I, particleID;
	unsigned int            fieldComponentCount = feVariable->fieldComponentCount;
	FeMesh*                 feMesh = feVariable->feMesh;
	double*                 result = Memory_Alloc_Array( double, fieldComponentCount, "Answer to Life?" );

	localElementCount = FeMesh_GetElementLocalSize( feMesh );

	/* Loop over local elements in problem */
	for( lElement_I = 0 ; lElement_I < localElementCount ; lElement_I++ ) {
		cellID            = CellLayout_MapElementIdToCellId( swarm->cellLayout, lElement_I );
		cellParticleCount = swarm->cellParticleCountTbl[ cellID ];

		/* Loop over materialPoints in element */
		for ( cParticle_I = 0 ; cParticle_I < cellParticleCount ; cParticle_I++ ) {
			particle = (MaterialPoint*)Swarm_ParticleInCellAt( swarm, cellID, cParticle_I );
			particleID = swarm->cellParticleTbl[cellID][cParticle_I];
			/* Interpolate field to particle location (is Global because it uses the
			 * material Points Swarm, could be a bit slow ???*/
			_FeVariable_InterpolateValueAt( feVariable, particle->coord, result );

			self->_printFunc( particle, particleID, result, fieldComponentCount, outputFile );

		}
	}
	Memory_Free( result );
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

