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
#include <StgDomain/StgDomain.h>
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

Underworld_SwarmOutput* _Underworld_SwarmOutput_New(  UNDERWORLD_SWARMOUTPUT_DEFARGS  )
{
	Underworld_SwarmOutput* self;
	
	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( _sizeOfSelf >= sizeof(Underworld_SwarmOutput) );
	/* The following terms are parameters that have been passed into this function but are being set before being passed onto the parent */
	/* This means that any values of these parameters that are passed into this function are not passed onto the parent function
	   and so should be set to ZERO in any children of this class. */
	nameAllocationType = NON_GLOBAL;

	self = (Underworld_SwarmOutput*)_Stg_Component_New(  STG_COMPONENT_PASSARGS  );

	
	/* General info */

	/* Virtual Info */
	
	return self;
}


/*------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _Underworld_SwarmOutput_Delete( void* uwSwarmOutput ) {
	Underworld_SwarmOutput* self = (Underworld_SwarmOutput*)uwSwarmOutput;

	Memory_Free( self->feVariableList );
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
	/* Variables set in this function */
	SizeT                                              _sizeOfSelf = sizeof(Underworld_SwarmOutput);
	Type                                                      type = Underworld_SwarmOutput_Type;
	Stg_Class_DeleteFunction*                              _delete = _Underworld_SwarmOutput_Delete;
	Stg_Class_PrintFunction*                                _print = _Underworld_SwarmOutput_Print;
	Stg_Class_CopyFunction*                                  _copy = _Underworld_SwarmOutput_Copy;
	Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _Underworld_SwarmOutput_DefaultNew;
	Stg_Component_ConstructFunction*                    _construct = _Underworld_SwarmOutput_AssignFromXML;
	Stg_Component_BuildFunction*                            _build = _Underworld_SwarmOutput_Build;
	Stg_Component_InitialiseFunction*                  _initialise = _Underworld_SwarmOutput_Initialise;
	Stg_Component_ExecuteFunction*                        _execute = _Underworld_SwarmOutput_Execute;
	Stg_Component_DestroyFunction*                        _destroy = _Underworld_SwarmOutput_Destroy;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

	return (void*) _Underworld_SwarmOutput_New(  UNDERWORLD_SWARMOUTPUT_PASSARGS  );
}

void _Underworld_SwarmOutput_Init( Underworld_SwarmOutput* self,
				PICelleratorContext*  context,
				MaterialPointsSwarm*  materialSwarm,
				unsigned int          listCount,
				FeVariable**          feVariableList) {
	self->materialSwarm = materialSwarm;
	self->sizeList = listCount;
	self->_getFeValuesFunc = _Underworld_SwarmOutput_GetFeVariableValues;
	self->_printFunc = _Underworld_SwarmOutput_PrintStandardFormat;
	self->feVariableList = feVariableList;

	/* my Swarm output will run on the SaveClass EP - the same time as standard checkpointing */
	EP_AppendClassHook( Context_GetEntryPoint( context, AbstractContext_EP_SaveClass ), _Underworld_SwarmOutput_Execute, self );
}

void _Underworld_SwarmOutput_AssignFromXML( void* uwSwarmOutput, Stg_ComponentFactory* cf, void* data ) {
	Underworld_SwarmOutput*  self          = (Underworld_SwarmOutput*) uwSwarmOutput;

	PICelleratorContext*    context;
	MaterialPointsSwarm*    materialSwarm;
	Dictionary*             dictionary = Dictionary_GetDictionary( cf->componentDict, self->name );
	Dictionary_Entry_Value* list;
	unsigned int            listCount, feVar_I;
	char*                   varName;
	Stream                  *errorStream = Journal_Register( Error_Type, "_Underworld_SwarmOutput_Construct" );
   FeVariable**            feVariableList;
	
	context      =  Stg_ComponentFactory_ConstructByName(  cf,  "context", PICelleratorContext,  True, data ) ;
	materialSwarm = (MaterialPointsSwarm*)Stg_ComponentFactory_ConstructByKey( cf, self->name, "Swarm", MaterialPointsSwarm, True, data );
	
	/* Get all Swarms specified in input file, under swarms */
	list = Dictionary_Get( dictionary, "FeVariables" );
	
	Journal_Firewall(
			list != NULL,
			errorStream,
			"Error in %s:\n"
			"You must specify a list of fevariables to interpolate onto the swarm in your xml. Example:\n"
			"<list name=\"FeVariables\">\n"
				"\t<param>PressureField</param>\n"
			"</list>\n", __func__ );
	
	listCount = Dictionary_Entry_Value_GetCount( list );

   Journal_Firewall(
			listCount != 0,
			errorStream,
			"Error in %s:\n"
			"You have no FeVariables defined in the FeVariable list. At least one FeVariable must be included. Example:\n"
			"<list name=\"FeVariables\">\n"
				"\t<param>PressureField</param>\n"
			"</list>\n", __func__ );
	
	/* Allocate the memory to store pointers to them */
	feVariableList = Memory_Alloc_Array( FeVariable*, listCount, "List FeVariables" );

	for( feVar_I = 0 ; feVar_I < listCount ; feVar_I++ ) {
		varName = Dictionary_Entry_Value_AsString( Dictionary_Entry_Value_GetElement( list, feVar_I ) );
		feVariableList[ feVar_I ] = Stg_ComponentFactory_ConstructByName( cf, varName, FeVariable, True, data );
	}

	_Underworld_SwarmOutput_Init( self,
				 context,
				 materialSwarm,
				 listCount,
				 feVariableList);

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
	char*                   outputPath = context->outputPath;
	unsigned int            feVar_I, feVarNum, timeStep;

	FILE*              outputFile;
	MPI_Comm	   comm;
	int                myRank;
	int                nProcs;
	MPI_Status         status;
	const int          FINISHED_WRITING_TAG = 100;
	int                canExecute = 0;

	
	timeStep = context->timeStep;
	feVarNum = self->sizeList;
	materialSwarm = self->materialSwarm;
	
	/* for each field create file and then write */
	for( feVar_I = 0 ; feVar_I < feVarNum; feVar_I++ ) {
		feVar = self->feVariableList[feVar_I];

		comm = Comm_GetMPIComm( Mesh_GetCommTopology( feVar->feMesh, MT_VERTEX ) );

		MPI_Comm_size( comm, (int*)&nProcs );
		MPI_Comm_rank( comm, (int*)&myRank );

		filename = Memory_Alloc_Array_Unnamed( char,
		        /*           OutputPath     .      FieldName        .      SwarmName                . time . dat\0 */
			       	strlen(outputPath) +1+ strlen(feVar->name) +1+ strlen(materialSwarm->name) +1+ 5 + 1+ 3 +1 );
		sprintf( filename, "%s/%s.%s.%.5u.dat", outputPath, feVar->name, materialSwarm->name, timeStep );
		/* wait for go-ahead from process ranked lower than me, to avoid competition writing to file */
		if ( myRank != 0  && canExecute == 0 ) {
			MPI_Recv( &canExecute, 1, MPI_INT, myRank - 1, FINISHED_WRITING_TAG, comm, &status );
		}	

		/* if myRank is 0, create or append file, otherwise append to file */
		if (myRank == 0) {
			/* append to file if restarting from checkpoint */
			if( context->loadFromCheckPoint )
				outputFile = fopen( filename, "a" );
			else
				outputFile = fopen( filename, "w" );
			fprintf( outputFile, "# FORMAT is:\n# MaterialID, Xpos, Ypos, Zpos, (values of field)\n" );
		} else {
		       	outputFile = fopen( filename, "a" );
		}

		self->_getFeValuesFunc( self, feVar, materialSwarm, outputFile );
		fclose( outputFile );
	}

	/* confirms this processor is finshed */
	canExecute = 1;
	/* send go-ahead from process ranked lower than me, to avoid competition writing to file */
	if ( myRank != nProcs - 1 ) {
		MPI_Ssend( &canExecute, 1, MPI_INT, myRank + 1, FINISHED_WRITING_TAG, comm );
	}

	Memory_Free( filename );
}

void _Underworld_SwarmOutput_Destroy( void* uwSwarmOutput, void* data ) {
	Underworld_SwarmOutput*	self = (Underworld_SwarmOutput*)uwSwarmOutput;
	unsigned int feVar_I, feVarNum;

	/* Destroy all StGermain based data structure this component uses */
	feVarNum = self->sizeList;	
	for( feVar_I = 0 ; feVar_I < feVarNum; feVar_I++ ) 
		Stg_Component_Destroy( self->feVariableList[feVar_I], data, False );

	Stg_Component_Destroy( self->materialSwarm, data, False );
	
}

void _Underworld_SwarmOutput_PrintStandardFormat( MaterialPoint* particle, double* result, unsigned fieldComponentCount, FILE* outputFile ) {
		unsigned dof_I;
		fprintf( outputFile, "%u ", particle->materialIndex );
                fprintf( outputFile, "%.15g %.15g %.15g ", particle->coord[0], particle->coord[1], particle->coord[2] ); 
		for ( dof_I = 0; dof_I < fieldComponentCount; dof_I++ ) 
			fprintf( outputFile, "%.15g ", result[dof_I] );
			
		fprintf( outputFile, "\n" );
}
	
void _Underworld_SwarmOutput_GetFeVariableValues(Underworld_SwarmOutput* uwSwarmOutput, FeVariable* feVariable, MaterialPointsSwarm* swarm, FILE* outputFile) {
	Underworld_SwarmOutput*	self = (Underworld_SwarmOutput*)uwSwarmOutput;
	MaterialPoint*          particle;
	unsigned int            localElementCount, lElement_I; 
	unsigned int            cellID, cellParticleCount, cParticle_I;
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
			/* Interpolate field to particle location (is Global because it uses the
			 * material Points Swarm, could be a bit slow ???*/
			_FeVariable_InterpolateValueAt( feVariable, particle->coord, result );

			self->_printFunc( particle, result, fieldComponentCount, outputFile );

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



