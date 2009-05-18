/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003-2006, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street,
**	Melbourne, 3053, Australia.
**
** Primary Contributing Organisations:
**	Victorian Partnership for Advanced Computing Ltd, Computational Software Development - http://csd.vpac.org
**	Australian Computational Earth Systems Simulator - http://www.access.edu.au
**	Monash Cluster Computing - http://www.mcc.monash.edu.au
**	Computational Infrastructure for Geodynamics - http://www.geodynamics.org
**
** Contributors:
**	Patrick D. Sunter, Software Engineer, VPAC. (pds@vpac.org)
**	Robert Turnbull, Research Assistant, Monash University. (robert.turnbull@sci.monash.edu.au)
**	Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**	David May, PhD Student, Monash University (david.may@sci.monash.edu.au)
**	Louis Moresi, Associate Professor, Monash University. (louis.moresi@sci.monash.edu.au)
**	Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**	Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
**	Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
**	Julian Giordani, Research Assistant, Monash University. (julian.giordani@sci.monash.edu.au)
**	Vincent Lemiale, Postdoctoral Fellow, Monash University. (vincent.lemiale@sci.monash.edu.au)
**
**  This library is free software; you can redistribute it and/or
**  modify it under the terms of the GNU Lesser General Public
**  License as published by the Free Software Foundation; either
**  version 2.1 of the License, or (at your option) any later version.
**
**  This library is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
**  Lesser General Public License for more details.
**
**  You should have received a copy of the GNU Lesser General Public
**  License along with this library; if not, write to the Free Software
**  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
**
** $Id: SolutionVector.c 1071 2008-03-12 02:23:49Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include "StgFEM/Discretisation/Discretisation.h"
#include "units.h"
#include "types.h"
#include "shortcuts.h"
#include "SolutionVector.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

/* Textual name of this class */
const Type SolutionVector_Type = "SolutionVector";

static const int VALUE_TAG = 1;
static const int VALUE_REQUEST_TAG = 2;

typedef struct RequestInfo {
	Node_LocalIndex		lNode_I;
	Dof_Index		nodeLocalDof_I;
} RequestInfo;

void _SolutionVector_ShareValuesNotStoredLocally(
	SolutionVector*		self,
	Index*			reqFromOthersCounts,
	RequestInfo**		reqFromOthersInfos,
	Dof_EquationNumber**	reqFromOthers,
	double*			localSolnVecValues );


void* SolutionVector_DefaultNew( Name name ) {
	return _SolutionVector_New( 
		sizeof(SolutionVector), 
		SolutionVector_Type, 
		_SolutionVector_Delete,
		_SolutionVector_Print, 
		_SolutionVector_Copy,
		SolutionVector_DefaultNew,
		_SolutionVector_Construct,
		_SolutionVector_Build, 
		_SolutionVector_Initialise, 
		_SolutionVector_Execute, 
		_SolutionVector_Destroy,
		name,
		False,
		MPI_COMM_WORLD, 
		NULL );
}

SolutionVector* SolutionVector_New(
		Name 						name,
		MPI_Comm					comm,
		FeVariable*					feVariable )
{
	return _SolutionVector_New( 
		sizeof(SolutionVector), 
		SolutionVector_Type, 
		_SolutionVector_Delete,
		_SolutionVector_Print, 
		_SolutionVector_Copy,
		SolutionVector_DefaultNew,
		_SolutionVector_Construct,
		_SolutionVector_Build, 
		_SolutionVector_Initialise, 
		_SolutionVector_Execute, 
		_SolutionVector_Destroy, 
		name,
		True,
		comm, 
		feVariable );
}


void SolutionVector_Init(
		SolutionVector*					self,
		Name 						name,
		MPI_Comm					comm,
		FeVariable*					feVariable )
{
	/* General info */
	self->type = SolutionVector_Type;
	self->_sizeOfSelf = sizeof(SolutionVector);
	self->_deleteSelf = False;
	
	/* Virtual info */
	self->_delete = _SolutionVector_Delete;
	self->_print = _SolutionVector_Print;
	self->_copy = _SolutionVector_Copy;
	self->_build = _SolutionVector_Build;

	_Stg_Class_Init( (Stg_Class*)self );
	_Stg_Object_Init( (Stg_Object*)self, name, NON_GLOBAL );
	_Stg_Component_Init( (Stg_Component*)self );

	_SolutionVector_Init( self, comm, feVariable );
}


SolutionVector* _SolutionVector_New( 
		SizeT						_sizeOfSelf,
		Type						type,
		Stg_Class_DeleteFunction*			_delete,
		Stg_Class_PrintFunction*			_print,
		Stg_Class_CopyFunction*				_copy, 
		Stg_Component_DefaultConstructorFunction*	_defaultConstructor,
		Stg_Component_ConstructFunction*		_construct,
		Stg_Component_BuildFunction*			_build,
		Stg_Component_InitialiseFunction*		_initialise,
		Stg_Component_ExecuteFunction*			_execute,
		Stg_Component_DestroyFunction*			_destroy,
		Name 						name, 
		Bool						initFlag,
		MPI_Comm					comm,
		FeVariable*					feVariable )
{
	SolutionVector*		self;
	
	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(SolutionVector) );
	self = (SolutionVector*)_Stg_Component_New( _sizeOfSelf, type, _delete, _print, _copy, _defaultConstructor, _construct,
			_build, _initialise, _execute, _destroy, name, NON_GLOBAL );
	
	/* General info */
	
	/* Virtual functions */
	
	if( initFlag ){
		_SolutionVector_Init( self, comm, feVariable );
	}
	
	return self;
}

void _SolutionVector_Init(
		SolutionVector*					self,
		MPI_Comm					comm,
		FeVariable*					feVariable )
{
	/* General and Virtual info should already be set */
	
	/* SolutionVector info */
	self->isConstructed = True;
	self->debug = Stream_RegisterChild( StgFEM_SLE_SystemSetup_Debug, self->type );
	self->comm = comm;
	self->feVariable = feVariable;
}


void _SolutionVector_Delete( void* solutionVector ) {
	SolutionVector* self = (SolutionVector*)solutionVector;
	
	Journal_DPrintf( self->debug, "In %s - for soln. vector %s\n", __func__, self->name );
	Stream_IndentBranch( StgFEM_Debug );
	//FreeObject( self->vector );
	if( self->vector != PETSC_NULL )
		VecDestroy( self->vector );
	
	/* Stg_Class_Delete parent*/
	_Stg_Component_Delete( self );
	Stream_UnIndentBranch( StgFEM_Debug );
}


void _SolutionVector_Print( void* solutionVector, Stream* stream ) {
	SolutionVector* self = (SolutionVector*)solutionVector;
	
	/* Set the Journal for printing informations */
	Stream* solutionVectorStream = stream;
	
	/* General info */
	Journal_Printf( solutionVectorStream, "SolutionVector (ptr): %p\n", self );
	
	/* Print parent */
	_Stg_Component_Print( self, solutionVectorStream );
	
	/* Virtual info */
	
	/* SolutionVector info */

	Stg_Class_Print( self->feVariable, solutionVectorStream );
	Journal_Printf( solutionVectorStream, "\tComm: %u\n", self->comm );
}


void* _SolutionVector_Copy( void* solutionVector, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
	SolutionVector*	self = (SolutionVector*)solutionVector;
	SolutionVector*	newSolutionVector;
	PtrMap*		map = ptrMap;
	Bool		ownMap = False;
	
	if( !map ) {
		map = PtrMap_New( 10 );
		ownMap = True;
	}
	
	newSolutionVector = _Stg_Component_Copy( self, dest, deep, nameExt, map );
	
	/* TODO: copy vector? */
	newSolutionVector->vector = self->vector;
	newSolutionVector->comm = self->comm;
	
	if( deep ) {
		newSolutionVector->debug = (Stream*)Stg_Class_Copy( self->debug, NULL, deep, nameExt, map );
		newSolutionVector->feVariable = (FeVariable*)Stg_Class_Copy( self->feVariable, NULL, deep, nameExt, map );
	}
	else {
		newSolutionVector->debug = self->debug;
		newSolutionVector->feVariable = self->feVariable;
	}
	
	if( ownMap ) {
		Stg_Class_Delete( map );
	}
	
	return (void*)newSolutionVector;
}


void _SolutionVector_Construct( void* solutionVector, Stg_ComponentFactory* cf, void* data ) 
{
	SolutionVector*    self = (SolutionVector*)solutionVector;
	FeVariable*        feVariable = NULL;
	
	feVariable = Stg_ComponentFactory_ConstructByKey( cf, self->name, "FeVariable", FeVariable, True, data ) ;
	_SolutionVector_Init( self, MPI_COMM_WORLD, (FeVariable*)feVariable );
}


void _SolutionVector_Build( void* solutionVector, void* data ) {
	SolutionVector* self = (SolutionVector*)solutionVector;
	
	Journal_DPrintf( self->debug, "In %s - for \"%s\"\n", __func__, self->name );
	Stream_IndentBranch( StgFEM_Debug );

	/* ensure variables are built */
	if( self->feVariable )
		Stg_Component_Build( self->feVariable, 0, False );

	/* Allocate the vector */
//	self->vector = PETScVector_New( "" );
//	Vector_SetComm( self->vector, self->comm );
//	Vector_SetLocalSize( self->vector, self->feVariable->eqNum->localEqNumsOwnedCount );
	VecCreate( self->comm, &self->vector );
	VecSetSizes( self->vector, self->feVariable->eqNum->localEqNumsOwnedCount, PETSC_DECIDE );
	VecSetFromOptions( self->vector );
#if( PETSC_VERSION_MAJOR <= 2 && PETSC_VERSION_MINOR >= 3 && PETSC_VERSION_SUBMINOR >= 3 )
	VecSetOption( self->vector, VEC_IGNORE_NEGATIVE_INDICES );
#elif( PETSC_VERSION_MAJOR >= 3 )
	VecSetOption( self->vector, VEC_IGNORE_NEGATIVE_INDICES, PETSC_TRUE );
#endif

	Stream_UnIndentBranch( StgFEM_Debug );
}


void _SolutionVector_Initialise( void* solutionVector, void* data ) {
	SolutionVector*          self = (SolutionVector *)solutionVector;
	
	Journal_DPrintf( self->debug, "In %s - for \"%s\"\n", __func__, self->name );
	Stream_IndentBranch( StgFEM_Debug );
	/* ensure variables are initialised */
	if( self->feVariable ) {
		Stg_Component_Initialise( self->feVariable, 0, False );
	}
		
	Stream_UnIndentBranch( StgFEM_Debug );
}


void _SolutionVector_Execute( void* solutionVector, void* data ) {
}

void _SolutionVector_Destroy( void* solutionVector, void* data ) {
}


void SolutionVector_ApplyBCsToVariables( void* solutionVector, void* data ) {
	SolutionVector*		self = (SolutionVector *)solutionVector;

	FeVariable_ApplyBCs( self->feVariable, data );
}

/* from the depreciated Vector class */
_SolutionVector_VectorView( Vec v, Stream* stream ) {
	unsigned	entry_i;
	PetscInt	size;
	PetscScalar*	array;

	VecGetSize( v, &size );
	VecGetArray( v, &array );

	Journal_Printf( stream, "%p = [", v );
	for( entry_i = 0; entry_i < size; entry_i++ ) 
		Journal_Printf( stream, "\t%u: \t %.12g\n", entry_i, array[entry_i] );
	Journal_Printf( stream, "];\n" );

	VecRestoreArray( v, &array );
}

void SolutionVector_UpdateSolutionOntoNodes( void* solutionVector ) {
	SolutionVector*		self = (SolutionVector *)solutionVector;
	double*			localSolnVecValues;
	Node_LocalIndex 	lNode_I = 0;
	Dof_Index		currNodeNumDofs;
	Dof_Index		nodeLocalDof_I;
	Partition_Index		ownerProc;
	FeVariable*		feVar = self->feVariable;
	FeMesh*			feMesh = feVar->feMesh;
	MPI_Comm		mpiComm;
	FeEquationNumber*	eqNum = feVar->eqNum;
	Dof_EquationNumber	currEqNum;
	Index			indexIntoLocalSolnVecValues;
	Index*			reqFromOthersCounts;
	Index*			reqFromOthersSizes;
	RequestInfo**		reqFromOthersInfos;
	Dof_EquationNumber**	reqFromOthers;
	Comm*			comm;
	Partition_Index		nProc;
	Partition_Index		myRank;
	Partition_Index		proc_I;
	double			initialGuessAtNonLocalEqNumsRatio = 0.1;
	double			ratioToIncreaseRequestArraySize = 1.5;
	Index			newReqFromOthersSize;
	
	Journal_DPrintf( self->debug, "In %s - for \"%s\"\n", __func__, self->name );
	Stream_IndentBranch( StgFEM_Debug );

	#if DEBUG
	if ( Stream_IsPrintableLevel( self->debug, 3 ) ) {
		Journal_DPrintf( self->debug, "Vector data:\n" );
		_SolutionVector_VectorView( self->vector, self->debug );
	}
	#endif

	comm = Mesh_GetCommTopology( feMesh, MT_VERTEX );
	mpiComm = Comm_GetMPIComm( comm );
	MPI_Comm_size( mpiComm, (int*)&nProc );
	MPI_Comm_rank( mpiComm, (int*)&myRank );

	/* allocate arrays for nodes that I want on each processor */
	reqFromOthersCounts = Memory_Alloc_Array( Index, nProc, "reqFromOthersCounts" );
	reqFromOthersSizes = Memory_Alloc_Array( Index, nProc, "reqFromOthersSizes" );
	reqFromOthersInfos = Memory_Alloc_Array( RequestInfo*, nProc, "reqFromOthersInfos" );
	reqFromOthers = Memory_Alloc_Array( Dof_EquationNumber*, nProc, "reqFromOthers" );
	/* Allocate the arrays of req. values from others independently, as we don't know how large they'll be */
	for ( proc_I=0; proc_I < nProc; proc_I++ ) {
		reqFromOthersCounts[proc_I] = 0;

		if (proc_I == myRank) continue;

		/* Our initial guess at number of non-local eqNums is a small ratio of the number of local dofs */
		reqFromOthersSizes[proc_I] = eqNum->localEqNumsOwnedCount * initialGuessAtNonLocalEqNumsRatio;
		/* Special case for really small meshes: make sure it's at least 1 */
		if (0 == reqFromOthersSizes[proc_I] ) {
			reqFromOthersSizes[proc_I]++;
		}
		reqFromOthersInfos[proc_I] = Memory_Alloc_Array( RequestInfo, reqFromOthersSizes[proc_I], 
			"reqFromOthersInfos[proc_I]" );
		reqFromOthers[proc_I] = Memory_Alloc_Array( Dof_EquationNumber, reqFromOthersSizes[proc_I],
			"reqFromOthers[proc_I]" );
	}
	
	/* Get the locally held part of the vector */
	//Vector_GetArray( self->vector, &localSolnVecValues );
	VecGetArray( self->vector, &localSolnVecValues );
	
	for( lNode_I=0; lNode_I < Mesh_GetLocalSize( feMesh, MT_VERTEX ); lNode_I++ ) {
		currNodeNumDofs = feVar->dofLayout->dofCounts[ lNode_I ];
		Journal_DPrintfL( self->debug, 3, "getting solutions for local node %d, has %d dofs.\n", lNode_I, currNodeNumDofs );
		
		/* process each dof */
		for ( nodeLocalDof_I = 0; nodeLocalDof_I < currNodeNumDofs; nodeLocalDof_I++ ) {
			Journal_DPrintfL( self->debug, 3, "\tdof %d: ", nodeLocalDof_I );
			
			currEqNum = eqNum->destinationArray[lNode_I][nodeLocalDof_I];
			if( currEqNum != -1 ) {
				Journal_DPrintfL( self->debug, 3, "is unconstrained, eqNum %d:", currEqNum );

				if( STreeMap_HasKey( eqNum->ownedMap, &currEqNum ) ) {
					indexIntoLocalSolnVecValues = *(int*)STreeMap_Map( eqNum->ownedMap, &currEqNum );
					Journal_DPrintfL( self->debug, 3, "local -> just copying value %f\n",
						localSolnVecValues[indexIntoLocalSolnVecValues] );
					DofLayout_SetValueDouble( feVar->dofLayout, lNode_I, nodeLocalDof_I,
						localSolnVecValues[indexIntoLocalSolnVecValues] );
				}
				else {
					RequestInfo*	requestInfo;
					
					Journal_DPrintfL( self->debug, 3, "nonlocal -> add to req list " );
					ownerProc = FeEquationNumber_CalculateOwningProcessorOfEqNum( eqNum, currEqNum ); 
					Journal_DPrintfL( self->debug, 3, "from proc %d\n", ownerProc );
					/* first check count & realloc if necessary */
					if (reqFromOthersCounts[ownerProc] == reqFromOthersSizes[ownerProc] ) {
						newReqFromOthersSize = reqFromOthersSizes[ownerProc] * ratioToIncreaseRequestArraySize;
						if ( newReqFromOthersSize == reqFromOthersSizes[ownerProc] ) {
							/* Special case: always increase by at least 1 */
							newReqFromOthersSize++;
						}
						reqFromOthersSizes[ownerProc] = newReqFromOthersSize;
						
						Journal_DPrintfL( self->debug, 3, "req list from proc %d count %d now "
							"equal to size, so reallocing to size %d\n", 
							ownerProc, reqFromOthersCounts[ownerProc],
							reqFromOthersSizes[ownerProc] );
							
						reqFromOthersInfos[ownerProc] = Memory_Realloc_Array(
							reqFromOthersInfos[ownerProc], RequestInfo, reqFromOthersSizes[ownerProc] );
						reqFromOthers[ownerProc] = Memory_Realloc_Array(
							reqFromOthers[ownerProc], Dof_EquationNumber, reqFromOthersSizes[ownerProc] );
					}
					requestInfo = &reqFromOthersInfos[ownerProc][ reqFromOthersCounts[ownerProc] ];
					requestInfo->lNode_I = lNode_I;
					requestInfo->nodeLocalDof_I = nodeLocalDof_I;
					reqFromOthers[ownerProc][reqFromOthersCounts[ownerProc]] = currEqNum;
					(reqFromOthersCounts[ownerProc])++;
				}	
			}
			else {
				Journal_DPrintfL( self->debug, 3, "is a BC, so skipping...\n" );
			}
		}
	}				
	
	if ( nProc > 1 ) {
		_SolutionVector_ShareValuesNotStoredLocally( self, reqFromOthersCounts, reqFromOthersInfos, reqFromOthers,
			localSolnVecValues );
	}	

	for ( proc_I=0; proc_I < nProc; proc_I++ ) {
		if (proc_I == myRank) continue;
		Memory_Free( reqFromOthers[proc_I] );
		Memory_Free( reqFromOthersInfos[proc_I] );
	}
	Memory_Free( reqFromOthers );
	Memory_Free( reqFromOthersInfos );
	Memory_Free( reqFromOthersCounts );
	Memory_Free( reqFromOthersSizes );

	//Vector_RestoreArray( self->vector, &localSolnVecValues );
	VecRestoreArray( self->vector, &localSolnVecValues );

	/*
	** Syncronise the FEVariable in question.
	*/

	FeVariable_SyncShadowValues( feVar );

	Stream_UnIndentBranch( StgFEM_Debug );
}


void _SolutionVector_ShareValuesNotStoredLocally(
	SolutionVector*		self,
	Index*			reqFromOthersCounts,
	RequestInfo**		reqFromOthersInfos,
	Dof_EquationNumber**	reqFromOthers,
	double*			localSolnVecValues )
{

	FeVariable*		feVar = self->feVariable;
	FeMesh*			feMesh = feVar->feMesh;
	FeEquationNumber*	eqNum = feVar->eqNum;
	Comm*			comm;
	MPI_Comm		mpiComm;
	Partition_Index		nProc;
	Partition_Index		myRank;
	Partition_Index		proc_I;
	Index			req_I;
	Index			indexIntoLocalSolnVecValues;
	MPI_Status		status;
	Index*			reqFromMeCounts;
	Dof_EquationNumber**	reqFromMe;
	double**		reqValuesFromMe;
	MPI_Request**		reqValuesFromMeHandles;
	MPI_Request**		reqFromOthersHandles;
	double**		reqValuesFromOthers;
	MPI_Request**		reqValuesFromOthersHandles;
	Bool*			reqValuesFromOthersReceived;
	Partition_Index	     reqValueSetsFromOthersNotYetReceivedCount;
	Dof_EquationNumber   totalRequestedFromOthers = 0;
	Dof_EquationNumber   totalRequestedFromMe = 0;

	Journal_DPrintf( self->debug, "In %s - for \"%s\"\n", __func__, self->name );
	Stream_IndentBranch( StgFEM_Debug );

	comm = Mesh_GetCommTopology( feMesh, MT_VERTEX );
	mpiComm = Comm_GetMPIComm( comm );
	MPI_Comm_size( mpiComm, (int*)&nProc );
	MPI_Comm_rank( mpiComm, (int*)&myRank );

	reqFromMeCounts = Memory_Alloc_Array( Index, nProc, "reqFromMeCounts" );
	reqFromOthersHandles = Memory_Alloc_Array_Unnamed( MPI_Request*, nProc );
	reqValuesFromOthersHandles = Memory_Alloc_Array_Unnamed( MPI_Request*, nProc );
	reqValuesFromMeHandles = Memory_Alloc_Array_Unnamed( MPI_Request*, nProc );
	reqValuesFromOthers = Memory_Alloc_2DComplex( double, nProc, reqFromOthersCounts, "reqValuesFromOthers" );
	reqValuesFromOthersReceived = Memory_Alloc_Array_Unnamed( Bool, nProc );

	#if DEBUG
	if ( Stream_IsPrintableLevel( self->debug, 2 ) ) {
		Journal_DPrintf( self->debug, "Final list of vec values I need from other procs:\n" );
		for ( proc_I=0; proc_I < nProc; proc_I++ ) {
			if ( proc_I == myRank ) continue;
			Journal_DPrintf( self->debug, "\t%d[0-%d]: ", proc_I, reqFromOthersCounts[proc_I] );
			for ( req_I=0; req_I < reqFromOthersCounts[proc_I]; req_I++ ) {
				RequestInfo* reqInfo = &reqFromOthersInfos[proc_I][req_I];
				Journal_DPrintf( self->debug, "(lnode %d, dof %d -> %d ), ",
					reqInfo->lNode_I, reqInfo->nodeLocalDof_I,
					reqFromOthers[proc_I][req_I] );
			}
			Journal_DPrintf( self->debug, "\n" );
		}
	}	
	#endif
	
	/* send out my request counts, receive the req. counts others want from me */
	MPI_Alltoall( reqFromOthersCounts, 1, MPI_UNSIGNED,
		      reqFromMeCounts, 1, MPI_UNSIGNED, mpiComm );

	Journal_DPrintf( self->debug, "After MPI_Alltoall- counts are:\n" );
	totalRequestedFromOthers = 0;
	totalRequestedFromMe = 0;
	Stream_Indent( self->debug );
	Journal_DPrintf( self->debug, "reqFromOthersCounts: " );
	for ( proc_I=0; proc_I < nProc; proc_I++ ) {
		if ( proc_I == myRank ) continue;
		Journal_DPrintf( self->debug, "\tp%d:%d, ", proc_I, reqFromOthersCounts[proc_I] );
		totalRequestedFromOthers += reqFromOthersCounts[proc_I];
	}
	Journal_DPrintf( self->debug, "\n" );
	Journal_DPrintf( self->debug, "reqFromMeCounts: " );
	for ( proc_I=0; proc_I < nProc; proc_I++ ) {
		if ( proc_I == myRank ) continue;
		Journal_DPrintf( self->debug, "\tp%d:%d, ", proc_I, reqFromMeCounts[proc_I] );
		totalRequestedFromMe += reqFromMeCounts[proc_I];
	}
	Journal_DPrintf( self->debug, "\n" );
	Stream_UnIndent( self->debug );

	if ( ( totalRequestedFromOthers == 0) && (totalRequestedFromMe == 0) ) 
	{	
		Journal_DPrintf( self->debug, "No vector values either required from others or "
			"required by others from me, therefore cleaning up memory and returning.\n" );
		Memory_Free( reqFromMeCounts );
		Memory_Free( reqFromOthersHandles );
		Memory_Free( reqValuesFromOthersHandles );
		Memory_Free( reqValuesFromMeHandles );
		Memory_Free( reqValuesFromOthers ); 
		Memory_Free( reqValuesFromOthersReceived );
		Stream_UnIndentBranch( StgFEM_Debug );
		return;
	}

	Journal_DPrintfL( self->debug, 2, "Starting non-blocking sends of my lists of vector entry indices I want from others:\n" );
	Stream_Indent( self->debug );
	for( proc_I=0; proc_I < nProc; proc_I++) {
		if ( proc_I == myRank ) continue; 
/* Journal_Printf( Journal_Register( Info_Type, "mpi" ),  "!!! line %d, proc_I %d: count = %u\n", __LINE__, proc_I, reqFromOthersCounts[proc_I] ); */
		if ( reqFromOthersCounts[proc_I] > 0 ) {
			Journal_DPrintfL( self->debug, 2, "Sending to proc %d the list of %d vector entry indices I want from it:\n"
				"\t(tracking via reqFromOthersHandles[%d], tag %d)\n", proc_I,
				reqFromOthersCounts[proc_I], proc_I, VALUE_REQUEST_TAG );

			reqFromOthersHandles[proc_I] = Memory_Alloc_Unnamed( MPI_Request );
			MPI_Isend( reqFromOthers[proc_I], reqFromOthersCounts[proc_I], MPI_UNSIGNED,
				proc_I, VALUE_REQUEST_TAG, mpiComm, reqFromOthersHandles[proc_I] );
		}	
	}
	Stream_UnIndent( self->debug );


	Journal_DPrintfL( self->debug, 2, "Starting non-blocking receive of the vector entries I want from others:\n" );
	Stream_Indent( self->debug );
	for( proc_I=0; proc_I < nProc; proc_I++) {
		if ( proc_I == myRank ) continue; 
		if ( reqFromOthersCounts[proc_I] > 0 ) {
			Journal_DPrintfL( self->debug, 2, "Posting recv reqst from proc %d for the %d vector entries I want from it:\n"
				"\t(tracking via reqValuesFromOthersHandles[%d], tag %d)\n", proc_I,
				reqFromOthersCounts[proc_I], proc_I, VALUE_TAG );
			reqValuesFromOthersHandles[proc_I] = Memory_Alloc_Unnamed( MPI_Request );
			MPI_Irecv( reqValuesFromOthers[proc_I], reqFromOthersCounts[proc_I], MPI_DOUBLE,
				proc_I, VALUE_TAG, mpiComm, reqValuesFromOthersHandles[proc_I] );
		}	
	}
	Stream_UnIndent( self->debug );

	Journal_DPrintfL( self->debug, 2, "Starting blocking receive of the lists of vector entry indices "
		"others want from me:\n" );
	Stream_Indent( self->debug );
	reqFromMe = Memory_Alloc_2DComplex( Dof_EquationNumber, nProc, reqFromMeCounts, "reqFromMe" );
	reqValuesFromMe = Memory_Alloc_2DComplex( double, nProc, reqFromMeCounts, "reqValuesFromMe" );
	for( proc_I=0; proc_I < nProc; proc_I++) {
		if ( proc_I == myRank ) continue; 
/* /Journal_Printf( Journal_Register( Info_Type, "mpi" ),  "!!! line %d, proc_I %d: count = %u\n", __LINE__, proc_I, reqFromMeCounts[proc_I] ); */
		if ( reqFromMeCounts[proc_I] > 0 ) {
			MPI_Recv( reqFromMe[proc_I], reqFromMeCounts[proc_I], MPI_UNSIGNED,
				proc_I, VALUE_REQUEST_TAG, mpiComm, &status );
			Journal_DPrintfL( self->debug, 3, "Received a list of %u requested vector entry indices from proc %u, "
				"with tag %d\n", reqFromMeCounts[proc_I], proc_I, status.MPI_TAG );
		}	
	}
	Stream_UnIndent( self->debug );
	
	#if DEBUG
	if ( Stream_IsPrintableLevel( self->debug, 2 ) ) {
		Journal_DPrintf( self->debug, "Final lists of vector entry indices other procs want from me are:\n" );
		Stream_Indent( self->debug );
		for ( proc_I=0; proc_I < nProc; proc_I++ ) {
			if ( proc_I == myRank ) continue;
			if ( reqFromMeCounts[proc_I] > 0 ) {
				Journal_DPrintf( self->debug, "%d[0-%d]: ", proc_I, reqFromMeCounts[proc_I] );
				for ( req_I=0; req_I < reqFromMeCounts[proc_I]; req_I++ ) {
					Journal_DPrintf( self->debug, "(eqNum %d), ", reqFromMe[proc_I][req_I] );
				}
				Journal_DPrintf( self->debug, "\n" );
			}
		}
		Stream_UnIndent( self->debug );
	}	
	#endif
	
	/* for all those requested from me, non-blocking send out values */
	Journal_DPrintfL( self->debug, 2, "Beginning non-blocking send out of vector entry lists requested by others:\n" );
	Stream_Indent( self->debug );
	for( proc_I=0; proc_I < nProc; proc_I++) {
		if ( proc_I == myRank ) continue; 
		if ( reqFromMeCounts[proc_I] > 0 ) {
			Journal_DPrintfL( self->debug, 3, "list to proc %d is: ", proc_I );
			for ( req_I=0; req_I < reqFromMeCounts[proc_I]; req_I++ ) {
				/* look up and fill in correct value in array */
				indexIntoLocalSolnVecValues = *(int*)STreeMap_Map( eqNum->ownedMap,
										   reqFromMe[proc_I] + req_I );
				reqValuesFromMe[proc_I][req_I] = localSolnVecValues[indexIntoLocalSolnVecValues];
				Journal_DPrintfL( self->debug, 3, "%d=%f, ", reqFromMe[proc_I][req_I],
					reqValuesFromMe[proc_I][req_I] );
			}
			Journal_DPrintfL( self->debug, 3, "\n" );
			/* Non-blocking send out the now-complete list to this processor */
			reqValuesFromMeHandles[proc_I] = Memory_Alloc_Unnamed( MPI_Request );
			Journal_DPrintfL( self->debug, 2, "Sending to proc %d the list of %d vector entries they want:\n"
				"\t(tracking via reqValuesFromMe[%d], tag %d)\n", proc_I,
				reqFromMeCounts[proc_I], proc_I, VALUE_TAG );
			MPI_Isend( reqValuesFromMe[proc_I], reqFromMeCounts[proc_I], MPI_DOUBLE,
				proc_I, VALUE_TAG, mpiComm, reqValuesFromMeHandles[proc_I] );
		}	
	}
	Stream_UnIndent( self->debug );
	
	Journal_DPrintfL( self->debug, 1, "Starting iterative-test receive of the vector entries I "
		"requested from others:\n" );
	/* Set up an array for keeping track of who we've received things from
	 * already */
	reqValueSetsFromOthersNotYetReceivedCount = nProc-1;
	for( proc_I=0; proc_I < nProc; proc_I++) {
		if ( proc_I == myRank ) continue; 
		reqValuesFromOthersReceived[proc_I] = False;
		if ( reqFromOthersCounts[proc_I] == 0 ) {
			reqValueSetsFromOthersNotYetReceivedCount--;
		}	
	}	

	#if DEBUG
	Journal_DPrintfL( self->debug, 2, "(Expecting %d receives from procs: ",
		reqValueSetsFromOthersNotYetReceivedCount );
	for( proc_I=0; proc_I < nProc; proc_I++) {
		if ( proc_I == myRank ) continue; 
		if ( reqFromOthersCounts[proc_I] > 0 ) {
			Journal_DPrintfL( self->debug, 2, "%d, ", proc_I );
		}	
	}	
	Journal_DPrintfL( self->debug, 2, ")\n" );
	#endif

	Stream_Indent( self->debug );
	/* now update the values at nodes that I requested from others, as they come in */
	while ( reqValueSetsFromOthersNotYetReceivedCount ) {	
		int flag = 0;

		Journal_DPrintfL( self->debug, 3, "%d sets still to go...\n", reqValueSetsFromOthersNotYetReceivedCount );
		for( proc_I=0; proc_I < nProc; proc_I++) {
			if ( proc_I == myRank ) continue; 

			if ( (reqFromOthersCounts[proc_I] > 0) && (False == reqValuesFromOthersReceived[proc_I]) ) {
				MPI_Test( reqValuesFromOthersHandles[proc_I], &flag, &status );
				if ( !flag ) {
					/* No results yet from this proc -> continue to next. */
					continue;
				}
				else {
					RequestInfo* reqInfo;
					Journal_DPrintfL( self->debug, 2, "received some requested "
						"values (using reqValuesFromOthersHandles) from proc %d "
						"(with tag %d, exp %d):", proc_I, status.MPI_TAG, VALUE_TAG );
					/* go through each value received from that proc & update onto node */
					for ( req_I=0; req_I < reqFromOthersCounts[proc_I]; req_I++ ) {
						reqInfo = &reqFromOthersInfos[proc_I][req_I];
						Journal_DPrintfL( self->debug, 3, "(lnode %d, dof %d -> %d )=%f, ",
							reqInfo->lNode_I, reqInfo->nodeLocalDof_I,
							reqFromOthers[proc_I][req_I], reqValuesFromOthers[proc_I][req_I] );
						DofLayout_SetValueDouble( feVar->dofLayout, reqInfo->lNode_I, reqInfo->nodeLocalDof_I,
							reqValuesFromOthers[proc_I][req_I] );
					}
					Journal_DPrintfL( self->debug, 2, "\n" );
					reqValuesFromOthersReceived[proc_I] = True;
					reqValueSetsFromOthersNotYetReceivedCount--;
					Memory_Free( reqValuesFromOthersHandles[proc_I] );
				}	
			}	
		}
	}	
	Stream_UnIndent( self->debug );

	/* MPI_Wait to be sure all sends to others have completed */
	Journal_DPrintfL( self->debug, 2, "Making sure all comms of this function finished:...\n" );
	Stream_Indent( self->debug );

	Journal_DPrintfL( self->debug, 2, "Confirming completion of my sends of "
		"vector entry index lists I wanted from others were received:\n" );
	Stream_Indent( self->debug );
	for( proc_I=0; proc_I < nProc; proc_I++) {
		if ( proc_I == myRank ) continue; 
		if ( reqFromOthersCounts[proc_I] > 0 ) {
			MPI_Wait( reqFromOthersHandles[proc_I], MPI_STATUS_IGNORE ); 
			Journal_DPrintfL( self->debug, 2, "Confirmed wait on reqFromOthersHandles[%u]"
				"\n", proc_I );
			Memory_Free( reqFromOthersHandles[proc_I] );
		}	
	}	
	Stream_UnIndent( self->debug );
	Journal_DPrintfL( self->debug, 2, "done.\n" );

	Journal_DPrintfL( self->debug, 2, "Confirming completion of my sends of "
		"vector entry values requested by others were received:\n" );
	Stream_Indent( self->debug );
	for( proc_I=0; proc_I < nProc; proc_I++) {
		if ( proc_I == myRank ) continue; 
		if ( reqFromMeCounts[proc_I] > 0 ) {
			MPI_Wait( reqValuesFromMeHandles[proc_I], MPI_STATUS_IGNORE );
			Journal_DPrintfL( self->debug, 2, "Confirmed wait on reqValuesFromMeHandles[%u]"
				"\n", proc_I );
			Memory_Free( reqValuesFromMeHandles[proc_I] );
		}
	}	
	Stream_UnIndent( self->debug );
	Journal_DPrintfL( self->debug, 2, "done.\n" );

	Stream_UnIndent( self->debug );
	Journal_DPrintfL( self->debug, 2, "done.\n" );

	Memory_Free( reqFromMeCounts );
	Memory_Free( reqFromMe );
	Memory_Free( reqValuesFromMe );
	Memory_Free( reqValuesFromOthers );
	Memory_Free( reqValuesFromOthersReceived );
	Memory_Free( reqFromOthersHandles );
	Memory_Free( reqValuesFromOthersHandles );
	Memory_Free( reqValuesFromMeHandles );

	Stream_UnIndentBranch( StgFEM_Debug );
	return;
}


void SolutionVector_LoadCurrentFeVariableValuesOntoVector( void* solutionVector ) {
	SolutionVector*		self = (SolutionVector*)solutionVector;
	FeVariable*		feVar = self->feVariable;
	FeMesh*			feMesh = feVar->feMesh;
	Node_LocalIndex		node_lI = 0;
	Dof_Index		dof_I = 0;
	double			value = 0;
	Index			insertionIndex = 0;
	
	for ( node_lI = 0; node_lI < FeMesh_GetNodeLocalSize( feMesh ); node_lI++ ) {
		for ( dof_I = 0; dof_I < feVar->dofLayout->dofCounts[node_lI]; dof_I++ ) {
			value = DofLayout_GetValueDouble( feVar->dofLayout, node_lI, dof_I );
			insertionIndex = feVar->eqNum->destinationArray[node_lI][dof_I];
			//Vector_InsertEntries( self->vector, 1, &insertionIndex, &value );
			VecSetValues( self->vector, 1, &insertionIndex, &value, INSERT_VALUES );
		}	
	}

	//Vector_AssemblyBegin( self->vector );
	//Vector_AssemblyEnd( self->vector ); 
	VecAssemblyBegin( self->vector );
	VecAssemblyEnd( self->vector ); 
}
