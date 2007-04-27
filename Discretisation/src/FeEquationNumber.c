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
** $Id: FeEquationNumber.c 822 2007-04-27 06:20:35Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include "units.h"
#include "types.h"
#include "shortcuts.h"
#include "ElementType.h"
#include "ElementType_Register.h"
#include "Element.h"
#include "FeMesh.h"
#include "FeEquationNumber.h"
#include "LinkedDofInfo.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <mpi.h>

/*###### Typedefs and Structs ######*/

/** Textual name of this class */
const Type FeEquationNumber_Type = "FeEquationNumber";

/** struct to store sub-totals: what is the equation number up to at a given
    node? These can then be exchanged between processors */
typedef struct CritPointInfo {
	Node_GlobalIndex	index;
	Dof_EquationNumber	eqNum;
} CritPointInfo;


/** An element of linked list of critical point info. Several of the functions
    use this to keep track of key points */
typedef struct AddListEntry {
	CritPointInfo*		critPointInfo;
	struct AddListEntry*	next;
} AddListEntry;

/** Enum to say whetehr values at crit. nodes should be printed */
typedef enum PrintValuesFlag {
	DONT_PRINT_VALUES,
	PRINT_VALUES
} PrintValuesFlag;

/**	MPI datatype handle for efficiently exchanging CritPointInfo structures.
   \see FeEquationNumber_Create_CritPointInfo_MPI_Datatype() for where this
   handle is defined. */
MPI_Datatype MPI_critPointInfoType;

/*###### Private Function Declarations ######*/

#if 0
static void _FeEquationNumber_BuildRemappedNodeInfoTable( void* feEquationNumber );

static void _FeEquationNumber_CalculateDomainKnownCriticalPoints( FeEquationNumber* self, Node_DomainIndex nodeDomainCount,
								  CritPointInfo* mySetEnds, Index* const mySetEndsTotal,
								  Node_GlobalIndex* myWantedCriticalPoints, Index* const myWantedCriticalPointsTotal );
	
static void _FeEquationNumber_HandleNode( FeEquationNumber* self, const Node_DomainIndex dNode_I,
					  Dof_EquationNumber* const currEqNum );
	
static void _FeEquationNumber_PostProcessLinkedDofs( FeEquationNumber* self );

static void _FeEquationNumber_CalculateCritPointsIHave( FeEquationNumber* self,
							Node_GlobalIndex** const myWantedCriticalPointsPtr, Node_GlobalIndex myWantedCriticalPointsTotal,
							CritPointInfo** const critPointsIHave, Index* const critPointsIHaveTotal,
							CritPointInfo* const critPointsToSend, Index* const critPointsToSendTotal );
	
static void _FeEquationNumber_ShareCriticalPoints(
	FeEquationNumber* self,
	Node_GlobalIndex** const myCriticalPoints,
	Node_GlobalIndex myCriticalPointsTotal,
	Node_GlobalIndex** allCriticalPoints,
	Index** procCritPointsTotals, 
	Node_GlobalIndex* const maxCritPointsPerProc );

static void _FeEquationNumber_ShareCritPointInfo(
	FeEquationNumber* self,
	CritPointInfo** const myCritPointInfo,
	Index myCritPointInfoTotal,
	CritPointInfo** allCritPointInfo,
	Index** procCritPointInfoTotals,
	Index* const maxCritPointInfoPerProc,
	PrintValuesFlag printValuesFlag );

static void _FeEquationNumber_DoPartialTotals( FeEquationNumber* self,
					       CritPointInfo* const critPointsIHave, Index critPointsIHaveTotal,
					       CritPointInfo* const critPointsToSend, Index critPointsToSendTotal );

static void _FeEquationNumber_AddAllPartialTotals( FeEquationNumber* self, CritPointInfo* mySubTotals,
						   Index myCritPointInfoTotal, CritPointInfo* allSubTotals, Index* procCritPointInfoTotals,
						   Index maxSubTotalsPerProc );

Node_RemappedGlobalIndex _FeEquationNumber_RemapNode( 
	Mesh* mesh, 
	Index newDimOrder[3],
	Node_GlobalIndex gNode_I );

/** Tests if the critical point from another processor is held by ours.
    Is complicated by the possibility of remapping. */
static Bool _FeEquationNumber_IHaveCritPoint(
	FeEquationNumber* self,
	Node_RemappedGlobalIndex critPoint );
#endif

/*###### Function Definitions ######*/

/** Public constructor */
void* FeEquationNumber_DefaultNew( Name name )
{
	return _FeEquationNumber_New( sizeof(FeEquationNumber), FeEquationNumber_Type, _FeEquationNumber_Delete,
				      _FeEquationNumber_Print, _FeEquationNumber_Copy,
				      (Stg_Component_DefaultConstructorFunction*)FeEquationNumber_DefaultNew,
				      _FeEquationNumber_Construct, (Stg_Component_BuildFunction*)_FeEquationNumber_Build, 
				      (Stg_Component_InitialiseFunction*)_FeEquationNumber_Initialise,
				      _FeEquationNumber_Execute, _FeEquationNumber_Destroy, name, False, NULL, NULL, NULL, NULL );
}

FeEquationNumber* FeEquationNumber_New(
	Name						name,
	void* 						mesh,
	DofLayout*					dofLayout,
	VariableCondition*				bcs,
	LinkedDofInfo*					linkedDofInfo )
{
	return _FeEquationNumber_New( sizeof(FeEquationNumber), FeEquationNumber_Type, _FeEquationNumber_Delete,
				      _FeEquationNumber_Print, _FeEquationNumber_Copy, 
				      (Stg_Component_DefaultConstructorFunction*)FeEquationNumber_DefaultNew,
				      _FeEquationNumber_Construct, (Stg_Component_BuildFunction*)_FeEquationNumber_Build, 
				      (Stg_Component_InitialiseFunction*)_FeEquationNumber_Initialise,
				      _FeEquationNumber_Execute, _FeEquationNumber_Destroy, name, True, mesh, dofLayout, bcs, linkedDofInfo );
}


/** Constructor for when the FeEquationNumber is already allocated. */
void FeEquationNumber_Init(
	FeEquationNumber*				self, 
	Name						name,
	void* 						mesh,
	DofLayout*					dofLayout,
	VariableCondition*				bcs,
	LinkedDofInfo*					linkedDofInfo )
{
	/* General info */
	self->type = FeEquationNumber_Type;
	self->_sizeOfSelf = sizeof(FeEquationNumber);
	self->_deleteSelf = False;
	
	/* Virtual info */
	self->_delete = _FeEquationNumber_Delete;
	self->_print = _FeEquationNumber_Print;
	self->_copy = _FeEquationNumber_Copy;
	self->_defaultConstructor = (Stg_Component_DefaultConstructorFunction*)FeEquationNumber_DefaultNew;
	self->_construct = _FeEquationNumber_Construct;
	self->_build = (Stg_Component_BuildFunction*)_FeEquationNumber_Build;
	self->_initialise = (Stg_Component_InitialiseFunction*)_FeEquationNumber_Initialise;
	self->_execute = _FeEquationNumber_Execute;
	self->_destroy = _FeEquationNumber_Destroy;
	
	_Stg_Class_Init( (Stg_Class*)self );
	_Stg_Object_Init( (Stg_Object*)self, name, NON_GLOBAL );
	_Stg_Component_Init( (Stg_Component*)self );
	
	/* FeEquationNumber info */
	_FeEquationNumber_Init( self, mesh, dofLayout, bcs, linkedDofInfo );
}


/** Constructor implementation. */
FeEquationNumber* _FeEquationNumber_New(
	SizeT						_sizeOfSelf, 
	Type						type,
	Stg_Class_DeleteFunction*				_delete,
	Stg_Class_PrintFunction*				_print, 
	Stg_Class_CopyFunction*				_copy, 
	Stg_Component_DefaultConstructorFunction*	_defaultConstructor,
	Stg_Component_ConstructFunction*			_construct,
	Stg_Component_BuildFunction*		_build,
	Stg_Component_InitialiseFunction*		_initialise,
	Stg_Component_ExecuteFunction*		_execute,
	Stg_Component_DestroyFunction*		_destroy,
	Name							name,
	Bool						initFlag,
	void*						mesh,
	DofLayout*					dofLayout,
	VariableCondition*				bcs,
	LinkedDofInfo*					linkedDofInfo )
{
	FeEquationNumber* self;
	
	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(FeEquationNumber) );
	self = (FeEquationNumber*)_Stg_Component_New( _sizeOfSelf, type, _delete, _print, _copy, _defaultConstructor, _construct, _build, 
						      _initialise, _execute, _destroy, name, NON_GLOBAL );
	
	/* General info */
	
	/* Virtual info */
	self->_build = _build;
	self->_initialise = _initialise;
	
	/* Mesh info */
	if( initFlag ){
		_FeEquationNumber_Init( self, mesh, dofLayout, bcs, linkedDofInfo );
	}
	
	return self;
}


/** Constructor variables initialisation. Doesn't allocate any
    memory, just saves arguments and sets counters to 0. */
void _FeEquationNumber_Init(
	FeEquationNumber*				self, 
	void*						feMesh,
	DofLayout*					dofLayout,
	VariableCondition*				bcs,
	LinkedDofInfo*					linkedDofInfo )
{
	
	/* General and Virtual info should already be set */
	
	/* FinteElementMesh info */
	self->isConstructed = True;
	self->feMesh = (FeMesh*)feMesh;
	self->globalSumUnconstrainedDofs = 0;
	self->isBuilt = False;
	self->locationMatrixBuilt = False;
	self->_lowestLocalEqNum = -1;
	self->_lowestGlobalEqNums = NULL;
	self->_highestLocalEqNum = -1;
	self->dofLayout = dofLayout;
	self->bcs = bcs;
	self->linkedDofInfo = linkedDofInfo;
	self->remappingActivated = False;
	/* register streams */
	self->debug = Stream_RegisterChild( StgFEM_Discretisation_Debug, FeEquationNumber_Type );
	self->debugLM = Stream_RegisterChild( self->debug, "LM" );
	self->warning = Stream_RegisterChild( StgFEM_Warning, FeEquationNumber_Type );
	self->removeBCs = True;
	self->bcEqNums = NULL;
}

void _FeEquationNumber_Construct( void* feEquationNumber, Stg_ComponentFactory *cf, void* data ){
	
}
	
void _FeEquationNumber_Execute( void* feEquationNumber, void *data ){
	
}
	
void _FeEquationNumber_Destroy( void* feEquationNumber, void *data ){
	
}
	
/* Copy */

/** Stg_Class_Delete implementation. */
void _FeEquationNumber_Delete( void* feEquationNumber ) {
	FeEquationNumber* self = (FeEquationNumber*) feEquationNumber;

	Journal_DPrintfL( self->debug, 1, "In %s\n",  __func__ );
	Stream_IndentBranch( StgFEM_Debug );

	FreeArray( self->remappedNodeInfos );
	/* free destination array memory */
	Journal_DPrintfL( self->debug, 2, "Freeing I.D. Array\n" );
	FreeArray( self->destinationArray );
	
	if (self->locationMatrix) {
		Journal_DPrintfL( self->debug, 2, "Freeing Full L.M. Array\n" );
		FreeArray( self->locationMatrix );
	}

	if( self->bcEqNums ) {
		Stg_Class_Delete( self->bcEqNums );
	}
	
	/* Stg_Class_Delete parent */
	_Stg_Class_Delete( self );
	Stream_UnIndentBranch( StgFEM_Debug );
}


/** Print implementation */
void _FeEquationNumber_Print( void* mesh, Stream* stream ) {
	FeEquationNumber* self = (FeEquationNumber*)mesh;
	
	/* General info */
	Journal_Printf( stream, "FeEquationNumber (ptr): %p\n", self );
	
	/* Print parent */
	_Stg_Class_Print( self, stream );
	
	/* Virtual info */
	Journal_Printf( stream,  "\t_build (func ptr): %p\n", self->_build );
	Journal_Printf( stream,  "\t_intialise (func ptr): %p\n", self->_initialise );
	
	/* FeEquationNumber info */
	/* Don't print dofs or bcs as these will be printed by FeVariable */
	
	if ( self->destinationArray ) {
		FeEquationNumber_PrintDestinationArray( self, stream );
		FeEquationNumber_PrintLocationMatrix( self, stream );
	}
	else {
		Journal_Printf( stream, "\tdestinationArray: (null)... not built yet\n" );
		Journal_Printf( stream, "\tlocationMatrix: (null)... not built yet\n" );
	}
}


void* _FeEquationNumber_Copy( void* feEquationNumber, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
#if 0
	FeEquationNumber*	self = (FeEquationNumber*)feEquationNumber;
	Mesh*			mesh;
	FeEquationNumber*	newFeEquationNumber;
	PtrMap*			map = ptrMap;
	Bool			ownMap = False;

	mesh = (Mesh*)self->feMesh;
	
	if( !map ) {
		map = PtrMap_New( 10 );
		ownMap = True;
	}
	
	newFeEquationNumber = _Stg_Class_Copy( self, dest, deep, nameExt, map );
	
	/* Virtual methods */
	newFeEquationNumber->_build = self->_build;
	newFeEquationNumber->_initialise = self->_initialise;
	
	newFeEquationNumber->_highestLocalEqNum = self->_highestLocalEqNum;
	newFeEquationNumber->_lowestLocalEqNum = self->_lowestLocalEqNum;
	newFeEquationNumber->_eqNumsPerProcDivisor = self->_eqNumsPerProcDivisor;
	newFeEquationNumber->_eqNumsRemainder = self->_eqNumsRemainder;
	newFeEquationNumber->_remNotAddedChangeover = self->_remNotAddedChangeover;
	newFeEquationNumber->firstOwnedEqNum = self->firstOwnedEqNum;
	newFeEquationNumber->lastOwnedEqNum = self->lastOwnedEqNum;
	newFeEquationNumber->localEqNumsOwnedCount = self->localEqNumsOwnedCount;
	newFeEquationNumber->globalSumUnconstrainedDofs = self->globalSumUnconstrainedDofs;
	newFeEquationNumber->locationMatrixBuilt = self->locationMatrixBuilt;
	newFeEquationNumber->remappingActivated = self->remappingActivated;
	
	if( deep ) {
		newFeEquationNumber->debug = (Stream*)Stg_Class_Copy( self->debug, NULL, deep, nameExt, map );
		newFeEquationNumber->debugLM = (Stream*)Stg_Class_Copy( self->debugLM, NULL, deep, nameExt, map );
		newFeEquationNumber->warning = (Stream*)Stg_Class_Copy( self->warning, NULL, deep, nameExt, map );
		newFeEquationNumber->feMesh = (FeMesh*)Stg_Class_Copy( self->feMesh, NULL, deep, nameExt, map );
		newFeEquationNumber->dofLayout = (DofLayout*)Stg_Class_Copy( self->dofLayout, NULL, deep, nameExt, map );
		newFeEquationNumber->bcs = (VariableCondition*)Stg_Class_Copy( self->bcs, NULL, deep, nameExt, map );
		newFeEquationNumber->linkedDofInfo = (LinkedDofInfo*)Stg_Class_Copy( self->linkedDofInfo, NULL, deep, nameExt, map );
		
		if( (newFeEquationNumber->remappedNodeInfos = PtrMap_Find( map, self->remappedNodeInfos )) == NULL && self->remappedNodeInfos ) {
			Node_DomainIndex	nodeDomainCount;

			nodeDomainCount = Mesh_GetDomainSize( mesh, MT_VERTEX );
			
			newFeEquationNumber->remappedNodeInfos = Memory_Alloc_Array( RemappedNodeInfo, nodeDomainCount, "FeEquationNumber->remappedNodeInfos" );
			memcpy( newFeEquationNumber->remappedNodeInfos, self->remappedNodeInfos, sizeof(RemappedNodeInfo) * nodeDomainCount );
			PtrMap_Append( map, self->remappedNodeInfos, newFeEquationNumber->remappedNodeInfos );
		}
		
		if( (newFeEquationNumber->destinationArray = PtrMap_Find( map, self->destinationArray )) == NULL && self->destinationArray ) {
			Node_DomainIndex	nodeDomainCount;
			Node_DomainIndex	node_dI;

			nodeDomainCount = Mesh_GetDomainSize( mesh, MT_VERTEX );
			
			newFeEquationNumber->destinationArray = Memory_Alloc_2DComplex( Dof_EquationNumber, nodeDomainCount, self->dofLayout->dofCounts, "FeEquationNumber->destinationArray" );
			for( node_dI = 0; node_dI < nodeDomainCount; node_dI++ ) {
				memcpy( newFeEquationNumber->destinationArray[node_dI], self->destinationArray[node_dI], sizeof(Dof_EquationNumber) * self->dofLayout->dofCounts[node_dI] );
			}
			PtrMap_Append( map, self->destinationArray, newFeEquationNumber->destinationArray );
		}
		
		if( (newFeEquationNumber->_lowestGlobalEqNums = PtrMap_Find( map, self->_lowestGlobalEqNums )) == NULL && self->_lowestGlobalEqNums ) {
			Partition_Index		nProc;

			nProc = Mesh_GetCommTopology( mesh, MT_VERTEX )->nProcs;
			
			newFeEquationNumber->_lowestGlobalEqNums = Memory_Alloc_Array( Dof_EquationNumber, nProc, "FeEquationNumber->_lowestGlobalEqNums" );
			memcpy( newFeEquationNumber->_lowestGlobalEqNums, self->_lowestGlobalEqNums, sizeof(Dof_EquationNumber) * nProc );
			PtrMap_Append( map, self->_lowestGlobalEqNums, newFeEquationNumber->_lowestGlobalEqNums );
		}
		
		if( (newFeEquationNumber->locationMatrix = PtrMap_Find( map, self->locationMatrix )) == NULL && self->locationMatrix ) {
			FeMesh*	mesh = self->feMesh;
			Element_LocalIndex	lElement_I;
			Node_LocalIndex		numNodesThisElement;
			Node_LocalIndex		elLocalNode_I;
			Dof_Index		numDofsThisNode;
			Element_LocalIndex	elementLocalCount = mesh->elementLocalCount;
			Dof_Index**		dofCountsAtElementNodesArray;
			
			dofCountsAtElementNodesArray = Memory_Alloc_3DSetup( elementLocalCount, mesh->elementNodeCountTbl );
			for ( lElement_I = 0; lElement_I < elementLocalCount; lElement_I++ ) {
				numNodesThisElement = mesh->elementNodeCountTbl[lElement_I];
				
				for( elLocalNode_I = 0; elLocalNode_I < numNodesThisElement; elLocalNode_I++) {
					Node_LocalIndex		dNode_I = mesh->elementNodeTbl[lElement_I][elLocalNode_I];
					
					numDofsThisNode = self->dofLayout->dofCounts[dNode_I];
					dofCountsAtElementNodesArray[lElement_I][elLocalNode_I] = numDofsThisNode;
				}
			}
			
			newFeEquationNumber->locationMatrix = Memory_Alloc_3DComplex( Dof_EquationNumber, elementLocalCount, mesh->elementNodeCountTbl, dofCountsAtElementNodesArray, "FeEquationNumber->locationMatrix" );
			for( lElement_I = 0; lElement_I < elementLocalCount; lElement_I++ ) {
				for( elLocalNode_I = 0; elLocalNode_I < mesh->elementNodeCountTbl[lElement_I]; elLocalNode_I++ ) {
					memcpy( newFeEquationNumber->locationMatrix[lElement_I][elLocalNode_I], self->locationMatrix[lElement_I][elLocalNode_I], sizeof(Dof_EquationNumber) * dofCountsAtElementNodesArray[lElement_I][elLocalNode_I] );
				}
			}
			
			Memory_Free( dofCountsAtElementNodesArray );
			PtrMap_Append( map, self->locationMatrix, newFeEquationNumber->locationMatrix );
		}
	}
	else {
		newFeEquationNumber->debug = self->debug;
		newFeEquationNumber->debugLM = self->debugLM;
		newFeEquationNumber->warning = self->warning;
		newFeEquationNumber->feMesh = self->feMesh;
		newFeEquationNumber->dofLayout = self->dofLayout;
		newFeEquationNumber->bcs = self->bcs;
		newFeEquationNumber->linkedDofInfo = self->linkedDofInfo;
		newFeEquationNumber->remappedNodeInfos = self->remappedNodeInfos;
		newFeEquationNumber->destinationArray = self->destinationArray;
		newFeEquationNumber->_lowestGlobalEqNums = self->_lowestGlobalEqNums;
		newFeEquationNumber->locationMatrix = self->locationMatrix;
	}
	
	if( ownMap ) {
		Stg_Class_Delete( map );
	}
	
	return (void*)newFeEquationNumber;
#endif
	abort();
}


void FeEquationNumber_Build( void* feEquationNumber ) {
	FeEquationNumber* self = (FeEquationNumber*)feEquationNumber;
	
	if ( False == self->isBuilt ) {
		self->_build( self, NULL );
		self->isBuilt = True;
	}
}


void _FeEquationNumber_Build( void* feEquationNumber ) {
	FeEquationNumber* self = (FeEquationNumber*) feEquationNumber;

	assert(self);
	
	Journal_DPrintf( self->debug, "In %s:\n",  __func__ );
	Stream_IndentBranch( StgFEM_Debug );

#if 0
	/* If we have new mesh topology information, do this differently. */
	if( self->feMesh->topo->domains && self->feMesh->topo->domains[MT_VERTEX] ) {
#endif
		FeEquationNumber_BuildWithTopology( self );
#if 0
	}
	else {
		_FeEquationNumber_BuildRemappedNodeInfoTable( self );
		_FeEquationNumber_BuildDestinationArray( self );
		_FeEquationNumber_CalculateGlobalUnconstrainedDofTotal( self );
		_FeEquationNumber_CalculateEqNumsDecomposition( self );
	}
#endif

	if ( Stream_IsPrintableLevel( self->debug, 3 ) ) {
		FeEquationNumber_PrintDestinationArray( self, self->debug );
	}

	Stream_UnIndentBranch( StgFEM_Debug );
}

void FeEquationNumber_Initialise( void* feEquationNumber ) {
	FeEquationNumber* self = (FeEquationNumber*)feEquationNumber;
	
	self->_initialise( self, NULL );
	self->isInitialised = True;
}


/** Initialise implementation. Currently does nothing. */
void _FeEquationNumber_Initialise( void* feEquationNumber ) {
#if DEBUG
	FeEquationNumber* self = (FeEquationNumber*) feEquationNumber;
	assert(self);
#endif
}


#if 0
Node_RemappedGlobalIndex _FeEquationNumber_RemapNode( 
	HexaMD* hexaMD,
	Index newDimOrder[3],
	Node_GlobalIndex gNode_I )
{
	Node_RemappedGlobalIndex remappedGlobalIndex = 0;
	Index dim_I = 0;
	Index currNewDim = 0;
	IJK defaultIJK;
	IJK newOrderIJK;
	IJK	newOrderNodeCounts;
	
	Dimension_1DTo3D( gNode_I, hexaMD->nodeGlobal3DCounts, defaultIJK );

	for ( dim_I = I_AXIS; dim_I < 3; dim_I++ ) {
		currNewDim = newDimOrder[dim_I];
		newOrderIJK[dim_I] = defaultIJK[currNewDim]; 
		newOrderNodeCounts[dim_I] = hexaMD->nodeGlobal3DCounts[currNewDim];
	}

	Dimension_3DTo1D( newOrderIJK, newOrderNodeCounts, &remappedGlobalIndex );

	return remappedGlobalIndex;
}


static void _FeEquationNumber_BuildRemappedNodeInfoTable( void* feEquationNumber ) {
	FeEquationNumber* self = (FeEquationNumber*) feEquationNumber;
	FeMesh* mesh = self->feMesh;
	Node_DomainIndex nodeDomainCount = Mesh_GetDomainSize( mesh, MT_VERTEX );
	CommTopology*	commTopo = Mesh_GetCommTopology( mesh, MT_VERTEX );
	Node_GlobalIndex gNode_I = 0;
	Node_DomainIndex dNode_I = 0;
/*#if DEBUG*/
	RemappedNodeInfo_Index rNodeInfo_I;
/*#endif*/
	
	Journal_DPrintfL( self->debug, 1, "In %s:\n",  __func__ );
	Stream_IndentBranch( StgFEM_Debug );

	self->remappedNodeInfos = Memory_Alloc_Array( RemappedNodeInfo, nodeDomainCount, 
						      "FeEquationNumber->remappedNodeInfos" );
	
	if ( 1 == commTopo->nProcs ) {
		Journal_DPrintfL( self->debug, 1, "Serial code: No remapping required...filling remapping table with normal values.\n" );
		for (dNode_I = 0; dNode_I < nodeDomainCount; dNode_I++ ) {
			gNode_I = mesh->nodeD2G[dNode_I];
			self->remappedNodeInfos[dNode_I].remappedGlobal = gNode_I;
			self->remappedNodeInfos[dNode_I].domain = dNode_I;
		}		
	}
	else if ( meshDecomp->type != HexaMD_Type ) {
		Node_GlobalIndex	prevGNode_I = (unsigned int)-1;
		RemappedNodeInfo_Index	insertIndex;
	
		Journal_DPrintfL( self->debug, 1, "Parallel, but non-hexa decomp: building the remapping table by "
				  "insertion-sorting based on global node number order.\n" );
		for (dNode_I = 0; dNode_I < nodeDomainCount; dNode_I++ ) {
			gNode_I = mesh->nodeD2G[dNode_I];

			if ( gNode_I > prevGNode_I ) {
				self->remappedNodeInfos[dNode_I].remappedGlobal = gNode_I;
				self->remappedNodeInfos[dNode_I].domain = dNode_I;
				prevGNode_I = gNode_I;
				continue;
			}
			else {
				insertIndex = 0;
				while ( gNode_I > self->remappedNodeInfos[insertIndex].remappedGlobal ) insertIndex++;
				memmove( &self->remappedNodeInfos[insertIndex+1], &self->remappedNodeInfos[insertIndex],
					 sizeof(RemappedNodeInfo) * (dNode_I - insertIndex + 1) );
				self->remappedNodeInfos[insertIndex].remappedGlobal = gNode_I;
				self->remappedNodeInfos[insertIndex].domain = dNode_I;
			}
		}
	}
	else if ( ((HexaMD*)meshDecomp)->numPartitionedDims > 1 )
	{
		Node_GlobalIndex	nextShadowGlobalMap;
		Node_LocalIndex		lNode_I;
		Node_ShadowIndex	sNode_I = mesh->nodeLocalCount;
		RemappedNodeInfo_Index	insertIndex = 0;
	
		Journal_DPrintfL( self->debug, 1, "Parallel, hexa decomp but multiple partition dims:\n"
				  "filling remapping table with insertion sort\n"
				  "of local & shadow elements.\n" );
		
		if ( meshDecomp->shadowDepth > 0 ) {
			nextShadowGlobalMap = mesh->nodeD2G[sNode_I];
		}
		else {
			nextShadowGlobalMap = mesh->nodeGlobalCount;
		}

		insertIndex = 0;
		for (lNode_I = 0; lNode_I < mesh->nodeLocalCount; lNode_I++ ) {
			gNode_I = mesh->nodeL2G[lNode_I];

			while ( nextShadowGlobalMap < gNode_I ) {
				self->remappedNodeInfos[insertIndex].remappedGlobal = nextShadowGlobalMap;
				self->remappedNodeInfos[insertIndex++].domain = sNode_I++;
				if ( sNode_I < mesh->nodeDomainCount ) {
					nextShadowGlobalMap = mesh->nodeD2G[sNode_I];
				}
				else {
					nextShadowGlobalMap = mesh->nodeGlobalCount;
				}
			}
			self->remappedNodeInfos[insertIndex].remappedGlobal = gNode_I;
			self->remappedNodeInfos[insertIndex++].domain = lNode_I;
		}		

		/* now know all local elements inserted, so ensure all shadow elements have also been inserted */
		for ( ; sNode_I < mesh->nodeDomainCount; sNode_I++ ) {
			gNode_I = mesh->nodeD2G[sNode_I];
			self->remappedNodeInfos[insertIndex].remappedGlobal = gNode_I;
			self->remappedNodeInfos[insertIndex++].domain = sNode_I;
		}
		Journal_Firewall( (insertIndex == mesh->nodeDomainCount), self->debug,
				  "Stuffup: should have inserted exactly the right number of values by here.\n" );
	}
	else {
		/* Otherwise, we remap each of the global node values, so the eqNums
		   will cause minimal communication during mesh assembly */
		Node_RemappedGlobalIndex currRemappedGNode_I = 0;
		Node_RemappedGlobalIndex firstRemappedGNode_I = 0;
		Index newDimOrder[3] = {0,0,0};
		Index newOrder_I = 0;
		Index partitionedAxis = I_AXIS;
		Index dim_I;
		RemappedNodeInfo_Index insertIndex;
		Bool lowestIsShadow = False;

		self->remappingActivated = True;

		Journal_DPrintfL( self->debug, 1, "Parallel, hexa decomp with 1 partition dim:\n"
				  "remapping eqNum traversal order to be aligned with mesh decomposition:\n" );

		/* Work out the new dimension order for the EqNums. The partitioned index
		   goes last */
		for ( dim_I = I_AXIS; dim_I < 3; dim_I++ ) {
			if (True == ((HexaMD*)meshDecomp)->partitionedAxis[dim_I]) {
				partitionedAxis = dim_I;
				newDimOrder[2] = partitionedAxis;
			} else {
				newDimOrder[newOrder_I++] = dim_I;
			}
		}
		Journal_DPrintfL( self->debug, 1, "Given partition axis is %c, calc. remapping dimension order as %c,%c,%c\n",
				  IJKTopology_DimNumToDimLetter[partitionedAxis],
				  IJKTopology_DimNumToDimLetter[newDimOrder[0]],
				  IJKTopology_DimNumToDimLetter[newDimOrder[1]],
				  IJKTopology_DimNumToDimLetter[newDimOrder[2]] );

		Journal_DPrintfL( self->debug, 4, "Calculating over %d domain nodes, %d local, %d shadow.\n",
				  mesh->nodeDomainCount, mesh->nodeLocalCount, mesh->nodeShadowCount );

		/* Work out the lowest remapped global index (so we know where to insert into the table) */
		gNode_I = mesh->nodeD2G[0];
		firstRemappedGNode_I = _FeEquationNumber_RemapNode( (HexaMD*)meshDecomp, newDimOrder, gNode_I );
		/* if shadow elements are enabled, check if the first shadow element has lower
		   global index than first local element */
		if ( meshDecomp->shadowDepth > 0 ) {
			Node_RemappedGlobalIndex firstRemappedShadowGNode_I;
		
			gNode_I = mesh->nodeD2G[mesh->nodeLocalCount];
			firstRemappedShadowGNode_I = _FeEquationNumber_RemapNode( (HexaMD*)meshDecomp,
										  newDimOrder, gNode_I );
			
			if ( firstRemappedShadowGNode_I < firstRemappedGNode_I ) {
				firstRemappedGNode_I = firstRemappedShadowGNode_I;
				lowestIsShadow = True;
			}
		}
		Journal_DPrintfL( self->debug, 4, "Calculated lowest remapped global index is %d",
				  firstRemappedGNode_I );
		if ( lowestIsShadow ) {
			Journal_DPrintfL( self->debug, 4, " (a shadow node.)\n" );
		}
		else {
			Journal_DPrintfL( self->debug, 4, " (a local node.)\n" );
		}

		Stream_IndentBranch( StgFEM_Debug );
		/* Now add the rest, relative to the first */
		for (dNode_I = 0; dNode_I < nodeDomainCount; dNode_I++ ) {
			gNode_I = mesh->nodeD2G[dNode_I];
			currRemappedGNode_I = _FeEquationNumber_RemapNode( (HexaMD*)meshDecomp, newDimOrder, gNode_I );
			insertIndex = currRemappedGNode_I - firstRemappedGNode_I;
			self->remappedNodeInfos[insertIndex].remappedGlobal = currRemappedGNode_I;
			self->remappedNodeInfos[insertIndex].domain = dNode_I;
			Journal_DPrintfL( self->debug, 4, "Processing domain node %d: original global = %d, remaps to %d.\n",
					  dNode_I, gNode_I, currRemappedGNode_I );
		}	
		Stream_UnIndentBranch( StgFEM_Debug );
	}

/*#if DEBUG*/
	if ( Stream_IsPrintableLevel( self->debug, 3 ) ) {
		Journal_DPrintf( self->debug, "Calculated remapped node info table as:\n{\n" );
		for (rNodeInfo_I = 0; rNodeInfo_I < nodeDomainCount; rNodeInfo_I++ ) {
			Journal_DPrintf( self->debug, " %d:(remap=%d,domain=%d),\n", rNodeInfo_I,
					 self->remappedNodeInfos[rNodeInfo_I].remappedGlobal,
					 self->remappedNodeInfos[rNodeInfo_I].domain );
		}
		Journal_DPrintf( self->debug, "}\n");
	}
/*#endif*/
	Stream_UnIndentBranch( StgFEM_Debug );
}	


/* Build the processor's ID (destination) array */
void _FeEquationNumber_BuildDestinationArray( FeEquationNumber* self ) {
	MeshDecomp* meshDecomp = self->feMesh->layout->decomp;
	Node_DomainIndex nodeDomainCount = meshDecomp->nodeDomainCount;
	CritPointInfo* critPointsIHave = NULL;
	Index critPointsIHaveTotal = 0;
	CritPointInfo* critPointsToSend = NULL;
	Index critPointsToSendTotal = 0;
	Stream* errorStream = Journal_Register( Error_Type, self->type );

	Journal_DPrintfL( self->debug, 1, "In %s:\n",  __func__ );
	Stream_IndentBranch( StgFEM_Debug );

	Journal_Firewall( ( nodeDomainCount == self->dofLayout->_numItemsInLayout ),
			  errorStream, "Error: In %s: DofLayout's size is %d, which isn't equal to "
			  "Node domain count of %d. Did you create the dofLayout of size node local count? "
			  "It should be aset up to size node domain count.\n",
			  __func__, self->dofLayout->_numItemsInLayout, nodeDomainCount );
		

	self->destinationArray = Memory_Alloc_2DComplex( Dof_EquationNumber, nodeDomainCount,
							 self->dofLayout->dofCounts, "FeEquationNumber->destinationArray" );


	if (meshDecomp->procsInUse == 1)
	{
		/* for serial jobs, don't worry about the critical points, just calculate
		   the totals normally */
		critPointsIHave = Memory_Alloc( CritPointInfo, "critPointsIHave (BuildDestinationArray)" );
		critPointsToSend = Memory_Alloc( CritPointInfo, "critPointsToSend (BuildDestinationArray)" );
		critPointsIHave[0].index = meshDecomp->nodeGlobalCount;
		critPointsToSend[0].index = meshDecomp->nodeGlobalCount;
		_FeEquationNumber_DoPartialTotals( self, critPointsIHave, 0, critPointsToSend, 0 );
		self->_lowestLocalEqNum = 0;
	}
	else {
		Node_GlobalIndex* myWantedCriticalPoints = NULL;
		Index myWantedCriticalPointsTotal = 0;
		CritPointInfo* allSendCriticalPoints = NULL;
		Index* procSendCritPointsTotals = NULL;
		Index maxSendCritPointsPerProc = 0;

		critPointsIHave = Memory_Alloc_Array( CritPointInfo, nodeDomainCount,
						      "critPointsIHave (BuildDestinationArray)" );
		critPointsToSend = Memory_Alloc_Array( CritPointInfo, nodeDomainCount,
						       "critPointsToSend (BuildDestinationArray)" );
		myWantedCriticalPoints = Memory_Alloc_Array( Node_GlobalIndex, nodeDomainCount,
							     "myWantedCritialPoints (BuildDestinationArray)" );

		/* Go through all domain nodes, work out end of runs, and interfaces needed */
		_FeEquationNumber_CalculateDomainKnownCriticalPoints( self, nodeDomainCount, critPointsIHave, &critPointsIHaveTotal,
								      myWantedCriticalPoints, &myWantedCriticalPointsTotal );
	
		/* initialise the below to the set end critical nodes I already worked out I have */
		critPointsToSend = (CritPointInfo*) memcpy( critPointsToSend, critPointsIHave, (critPointsIHaveTotal * sizeof(CritPointInfo)) );
		critPointsToSendTotal = critPointsIHaveTotal;
		
		/* work out critical points (inc. those needed by others) I hold, and those to send.
		   The list I have could be greater than my end of runs, in an irregular mesh case */
		_FeEquationNumber_CalculateCritPointsIHave( self, 
							    &myWantedCriticalPoints, myWantedCriticalPointsTotal, &critPointsIHave, &critPointsIHaveTotal,
							    critPointsToSend, &critPointsToSendTotal );
	
		Memory_Free( myWantedCriticalPoints );

		/* OK: We now know all the critical points on this processor, so go ahead and work out the partial ID number
		   at all of them. */

		Journal_DPrintfL( self->debug, 2, "Looping through nodes again, calculating eq nums (re-setting to 0 at critical points)\n,"
				  "and storing subtotal of each each critical node:\n" );
		_FeEquationNumber_DoPartialTotals( self, critPointsIHave, critPointsIHaveTotal,
						   critPointsToSend, critPointsToSendTotal );

		Journal_DPrintfL( self->debug, 2, "Now sub-totals have been calculated, share all of my critical point sub-totals "
				  "and get any that I need:\n" );
		_FeEquationNumber_ShareCritPointInfo( self, &critPointsToSend, critPointsToSendTotal,
						      &allSendCriticalPoints, &procSendCritPointsTotals, &maxSendCritPointsPerProc, PRINT_VALUES );

		Journal_DPrintfL( self->debug, 2, "Final loop: add the sub-totals at all the critical points to the existing values:\n" );
		_FeEquationNumber_AddAllPartialTotals( self, critPointsIHave, critPointsIHaveTotal,
						       allSendCriticalPoints, procSendCritPointsTotals, maxSendCritPointsPerProc );

		/* Post process the linked dofs */
		if ( self->linkedDofInfo ) {
			_FeEquationNumber_PostProcessLinkedDofs( self );
		}	

		Memory_Free( allSendCriticalPoints );
		Memory_Free( procSendCritPointsTotals );
	}

	/* If not removing BCs, construct a table of which equation numbers are actually BCs. */
	if( !self->removeBCs ) {
		FeMesh*	feMesh = self->feMesh;
		DofLayout*		dofLayout = self->dofLayout;
		VariableCondition*	bcs = self->bcs;
		unsigned		lNode_i;

		/* New index set. */
		self->bcEqNums = IndexSet_New( self->_highestLocalEqNum - self->_lowestLocalEqNum );

		/* Fill it up. */
		for( lNode_i = 0; lNode_i < feMesh->nodeLocalCount; lNode_i++ ) {
			unsigned	nDofs = dofLayout->dofCounts[lNode_i];
			unsigned	dof_i;

			for( dof_i = 0; dof_i < nDofs; dof_i++ ) {
				unsigned	varInd = dofLayout->varIndices[lNode_i][dof_i];

				if( bcs && VariableCondition_IsCondition( bcs, lNode_i, varInd ) ) {
					IndexSet_Add( self->bcEqNums, self->destinationArray[lNode_i][dof_i] );
				}
			}
		}
	}

	Memory_Free( critPointsIHave );
	Memory_Free( critPointsToSend );
	Stream_UnIndentBranch( StgFEM_Debug );
}




/** Works out 'critical nodes' I know from domain information (those I want are the start
    of my runs -1, those I know others want are the end of my runs. */
static void _FeEquationNumber_CalculateDomainKnownCriticalPoints(
	FeEquationNumber* self,
	Node_DomainIndex nodeDomainCount,
	CritPointInfo* mySetEnds,
	Index* const mySetEndsTotal,
	Node_GlobalIndex* myWantedCriticalPoints,
	Index* const myWantedCriticalPointsTotal )
{
	RemappedNodeInfo_Index rNodeInfo_I = 0;
	Node_RemappedGlobalIndex remappedGlobal_I = 0;
	MeshDecomp* meshDecomp = self->feMesh->layout->decomp;

	Journal_DPrintf( self->debug, "In %s\n", __func__ );
	Stream_IndentBranch( StgFEM_Debug );

	if ( self->remappingActivated ) {
		/* if remapping has been done, the only possible crit points are the first and last
		   points */
		remappedGlobal_I = self->remappedNodeInfos[0].remappedGlobal;
		if (remappedGlobal_I != 0) {
			Journal_DPrintfL( self->debug, 3, "Adding interface c.p. to myWantedCriticalPoints[%d] = %d\n", 
					  (*myWantedCriticalPointsTotal), remappedGlobal_I-1 );
			myWantedCriticalPoints[(*myWantedCriticalPointsTotal)++] = remappedGlobal_I-1;
		}	

		remappedGlobal_I = self->remappedNodeInfos[nodeDomainCount-1].remappedGlobal;
		if ( remappedGlobal_I != (meshDecomp->nodeGlobalCount-1) ) {
			Journal_DPrintfL( self->debug, 3, "Adding set end c.p. mySetEnds[%d] = %d\n", 
					  (*mySetEndsTotal), remappedGlobal_I );
			mySetEnds[(*mySetEndsTotal)++].index = remappedGlobal_I;
		}
		Stream_UnIndentBranch( StgFEM_Debug );
		return;
	}
	
	while ( rNodeInfo_I < nodeDomainCount ) {
		/* save start of set */
		remappedGlobal_I = self->remappedNodeInfos[rNodeInfo_I].remappedGlobal;

		if (remappedGlobal_I != 0) {
			Journal_DPrintfL( self->debug, 3, "Adding interface c.p. to myWantedCriticalPoints[%d] = %d\n", 
					  (*myWantedCriticalPointsTotal), remappedGlobal_I-1 );
			myWantedCriticalPoints[(*myWantedCriticalPointsTotal)++] = remappedGlobal_I-1;
		}	
		rNodeInfo_I++;

		/* skip to end of contiguous set */
		while ( (rNodeInfo_I < nodeDomainCount)
			&& (self->remappedNodeInfos[rNodeInfo_I].remappedGlobal ==
			    (self->remappedNodeInfos[rNodeInfo_I-1].remappedGlobal + 1) ) )
		{
			Journal_DPrintfL( self->debug, 4, "skipping remapped gNode[%d] domain = %d\n", 
					  self->remappedNodeInfos[rNodeInfo_I].remappedGlobal,
					  self->remappedNodeInfos[rNodeInfo_I].domain );
			rNodeInfo_I++;
		}	

		/* and add the critical point, if not the last node */
		remappedGlobal_I = self->remappedNodeInfos[rNodeInfo_I-1].remappedGlobal;
		if ( remappedGlobal_I != (meshDecomp->nodeGlobalCount-1) ) {
			Journal_DPrintfL( self->debug, 3, "Adding set end c.p. mySetEnds[%d] = %d\n", 
					  (*mySetEndsTotal), remappedGlobal_I );
			mySetEnds[(*mySetEndsTotal)++].index = remappedGlobal_I;
		}
	}

/*#if DEBUG*/
	if ( Stream_IsPrintableLevel( self->debug, 2 ) ) {
		Index cPoint_I;

		Journal_DPrintf( self->debug, "My set end Crit Points:[0-%d][ ",  (*mySetEndsTotal) );
		for (cPoint_I = 0 ; cPoint_I < (*mySetEndsTotal) ; cPoint_I++) {
			Journal_DPrintf( self->debug, "%2d, ", mySetEnds[cPoint_I].index );
		}	
		Journal_DPrintf( self->debug, "]\n" );

		Journal_DPrintf( self->debug, "interface Crit Points I want:[0-%d][ ",  (*myWantedCriticalPointsTotal) );
		for (cPoint_I = 0 ; cPoint_I < (*myWantedCriticalPointsTotal) ; cPoint_I++) {
			Journal_DPrintf( self->debug, "%2d, ", myWantedCriticalPoints[cPoint_I] );
		}
		Journal_DPrintf( self->debug, "]\n" );
	}
/*#endif*/
	Stream_UnIndentBranch( StgFEM_Debug );
}


/* work out critical points (inc. those needed by others) I hold, and those to send.
   The list I have could be greater than my end of runs, in an irregular mesh case */
static void _FeEquationNumber_CalculateCritPointsIHave( FeEquationNumber* self,
							Node_GlobalIndex** const myWantedCriticalPointsPtr, Node_GlobalIndex myWantedCriticalPointsTotal,
							CritPointInfo** const critPointsIHavePtr, Index* const critPointsIHaveTotal,
							CritPointInfo* const critPointsToSend, Index* const critPointsToSendTotal )
{
	MeshDecomp* meshDecomp = self->feMesh->layout->decomp;
	Partition_Index proc_I; 
	Node_GlobalIndex* allWantedCriticalPoints = NULL;
	CritPointInfo* allHaveCriticalPoints = NULL;
	Node_GlobalIndex* procWantedCritPointsTotals = NULL;
	Index* procHaveCritPointsTotals = NULL;
	Node_GlobalIndex maxWantedCritPointsPerProc = 0;
	Node_GlobalIndex maxHaveCritPointsPerProc = 0;
	Index myCritPoint_I = 0;
	PartitionIndex myRank = meshDecomp->rank;

	Journal_DPrintf( self->debug, "In %s\n",  __func__ );
	Stream_IndentBranch( StgFEM_Debug );

	Journal_DPrintfL( self->debug, 2, "Globally sharing locally known wanted C.N. lists:\n",  __func__ );
	_FeEquationNumber_ShareCriticalPoints( self, myWantedCriticalPointsPtr, myWantedCriticalPointsTotal,
					       &allWantedCriticalPoints, &procWantedCritPointsTotals, &maxWantedCritPointsPerProc );

	Journal_DPrintfL( self->debug, 2, "Globally sharing locally known 'have' C.N. lists:\n",  __func__ );
	_FeEquationNumber_ShareCritPointInfo( self, critPointsIHavePtr, (*critPointsIHaveTotal),
					      &allHaveCriticalPoints, &procHaveCritPointsTotals, &maxHaveCritPointsPerProc, DONT_PRINT_VALUES );

	/* now we build our lists of crit points on my processor, and crit points on my processor to send to others */
	for (proc_I = 0; proc_I < meshDecomp->nproc; proc_I++) {
		Index otherCritPoint_I = proc_I * maxWantedCritPointsPerProc;
		Index otherCritPoint_End_I = (proc_I * maxWantedCritPointsPerProc) + procWantedCritPointsTotals[proc_I];
		myCritPoint_I = 0;

		/* don't try and add my own */
		if ( myRank == proc_I ) continue;

		/* for each crit point they want, if I have it add it to my have and send lists */
		for (; otherCritPoint_I < otherCritPoint_End_I; otherCritPoint_I++) {
			Node_GlobalIndex otherCritPoint = allWantedCriticalPoints[otherCritPoint_I];
			
			if ( _FeEquationNumber_IHaveCritPoint( self, otherCritPoint ) ) {

				while ( (myCritPoint_I < (*critPointsIHaveTotal) ) &&
					(otherCritPoint > (*critPointsIHavePtr)[myCritPoint_I].index) ) {
					/* skip ahead to right spot */
					myCritPoint_I++;
				}
				if ( myCritPoint_I == (*critPointsIHaveTotal) ) {
					(*critPointsIHavePtr)[(*critPointsIHaveTotal)++].index = otherCritPoint; 
					critPointsToSend[(*critPointsToSendTotal)++].index = otherCritPoint;
				}
				else if ( otherCritPoint != (*critPointsIHavePtr)[myCritPoint_I].index ) {
					memmove( &((*critPointsIHavePtr)[myCritPoint_I+1]),
						 &((*critPointsIHavePtr)[myCritPoint_I]),
						 ((*critPointsIHaveTotal) - myCritPoint_I) * sizeof(CritPointInfo) );
					(*critPointsIHavePtr)[myCritPoint_I].index = otherCritPoint;
					(*critPointsIHaveTotal)++;
					memmove( &critPointsToSend[myCritPoint_I+1], &critPointsToSend[myCritPoint_I],
						 ((*critPointsToSendTotal) - myCritPoint_I) * sizeof(CritPointInfo) );
					critPointsToSend[myCritPoint_I].index = otherCritPoint;
					(*critPointsToSendTotal)++;
				}
				/* move to the next point in any case */
				myCritPoint_I++;
			}	
		}		
	}			

	/* For each crit point the others have, if I have it add it to my have list */
	for (proc_I = 0; proc_I < meshDecomp->nproc; proc_I++) {
		Index otherCritPointInfo_I = proc_I * maxHaveCritPointsPerProc;
		Index otherCritPointInfo_End_I = (proc_I * maxHaveCritPointsPerProc) + procHaveCritPointsTotals[proc_I];
		myCritPoint_I = 0;

		/* don't try and add my own */
		if ( myRank == proc_I ) continue;

		for (; otherCritPointInfo_I < otherCritPointInfo_End_I; otherCritPointInfo_I++) {
			Node_GlobalIndex otherCritPoint = allHaveCriticalPoints[otherCritPointInfo_I].index;
			
			if ( _FeEquationNumber_IHaveCritPoint( self, otherCritPoint ) ) {
			
				while ( (myCritPoint_I < (*critPointsIHaveTotal) ) &&
					(otherCritPoint > (*critPointsIHavePtr)[myCritPoint_I].index) ) {
					myCritPoint_I++;
					/* skip ahead to right spot */
				}
				if ( myCritPoint_I == (*critPointsIHaveTotal) ) {
					(*critPointsIHavePtr)[(*critPointsIHaveTotal)++].index = otherCritPoint; 
				}
				else if ( otherCritPoint != (*critPointsIHavePtr)[myCritPoint_I].index ) {
					memmove( &((*critPointsIHavePtr)[myCritPoint_I+1]),
						 &((*critPointsIHavePtr)[myCritPoint_I]),
						 ((*critPointsIHaveTotal) - myCritPoint_I) * sizeof(CritPointInfo) );
					(*critPointsIHavePtr)[myCritPoint_I].index = otherCritPoint;
					(*critPointsIHaveTotal)++;
				}
				/* move to the next point in any case */
				myCritPoint_I++;
			}	
		}		
	}			
/*#if DEBUG*/
	if ( Stream_IsPrintableLevel( self->debug, 2 ) ) {
		Journal_DPrintf( self->debug, "Calculated crit points I Have:\n[" );
		for ( myCritPoint_I = 0; myCritPoint_I < (*critPointsIHaveTotal); myCritPoint_I++ ) {
			Journal_DPrintf( self->debug, "%3d, ", (*critPointsIHavePtr)[myCritPoint_I].index );
		}
		Journal_DPrintf( self->debug, "]\nCalculated crit points To Send:\n[" );
		for ( myCritPoint_I = 0; myCritPoint_I < (*critPointsToSendTotal); myCritPoint_I++ ) {
			Journal_DPrintf( self->debug, "%3d, ", critPointsToSend[myCritPoint_I].index );
		}
		Journal_DPrintf( self->debug, "]\n");
	}	
/*#endif*/
	
	Memory_Free( allWantedCriticalPoints );
	Memory_Free( procWantedCritPointsTotals );
	Memory_Free( allHaveCriticalPoints );
	Memory_Free( procHaveCritPointsTotals );
	Stream_UnIndentBranch( StgFEM_Debug );
}


static Bool _FeEquationNumber_IHaveCritPoint(
	FeEquationNumber* self,
	Node_RemappedGlobalIndex critPoint )
{
	MeshDecomp* meshDecomp = self->feMesh->layout->decomp;
	Node_DomainIndex nodeDomainCount = meshDecomp->nodeDomainCount;
	PartitionIndex myRank = meshDecomp->rank;

	/* Case 1: for remapped hexa mesh, just check if its between max and min numbers */	
	if ( self->remappingActivated ) {
		if ( ( self->remappedNodeInfos[0].remappedGlobal <= critPoint ) &&
		     ( critPoint <= self->remappedNodeInfos[nodeDomainCount-1].remappedGlobal ) )
		{
			return True;
		}
		else {
			return False;
		}	
	}
	/* Case 2: 2D+ hexa decomps or other decomp strategy, use local (and shadow if enabled) node IndexSet */
	else {
		if ( IndexSet_IsMember( meshDecomp->localNodeSets[myRank], critPoint ) ) {
			return True;
		}
		else if ( (meshDecomp->shadowDepth > 0) 
			  && IndexSet_IsMember( meshDecomp->shadowNodeSets[myRank], critPoint ) ) {
			return True;
		}
		else {
			return False;
		}
	}
}


/** Performs an AllGather on a set of critical points: Each processer will afterwards have all the critical points
    of all the processors, in one large array. The user must then index into the array by processor carefully. */
static void _FeEquationNumber_ShareCriticalPoints(
	FeEquationNumber* self,
	Node_GlobalIndex** const myCriticalPoints,
	Node_GlobalIndex myCriticalPointsTotal,
	Node_GlobalIndex** allCriticalPoints,
	Index** procCritPointsTotals, 
	Node_GlobalIndex* const maxCritPointsPerProc)
{
	MeshDecomp* meshDecomp = self->feMesh->layout->decomp;
	Partition_Index proc_I; 
/*#if DEBUG*/
	Index point_I = 0;
/*#endif*/

	Journal_DPrintf( self->debug, "In %s\n", __func__ );
	Stream_IndentBranch( StgFEM_Debug );
	
	(*maxCritPointsPerProc) = 0;
	
/*#if DEBUG*/
	if ( Stream_IsPrintableLevel( self->debug, 2 ) ) {
		Journal_DPrintf( self->debug, "myCriticalPointsTotal=%u\n",  myCriticalPointsTotal);
		Journal_DPrintf( self->debug, "myCritPoints:(" );
		for ( point_I = 0; point_I < myCriticalPointsTotal; point_I++)
		{
			Journal_DPrintf( self->debug, "%u, ", (*myCriticalPoints)[point_I]);
		}	
		Journal_DPrintf( self->debug, ")\n");
	}	
/*#endif*/

	/* First, we allocate and calculate how many critical points are on each node (this way, we don't have to guess
	   how much memory to allocate in the main array.) */
	(*procCritPointsTotals) = Memory_Alloc_Array( Node_GlobalIndex, meshDecomp->nproc,
						      "*procCritPointsTotals (ShareCriticalPoints)" ); 
	
	MPI_Allgather( &myCriticalPointsTotal, 1, MPI_UNSIGNED,
		       (*procCritPointsTotals), 1, MPI_UNSIGNED,
		       meshDecomp->communicator );
	
	for (proc_I = 0; proc_I < meshDecomp->nproc; proc_I++) {
		if ( (*procCritPointsTotals)[proc_I] > (*maxCritPointsPerProc) )
			(*maxCritPointsPerProc) = (*procCritPointsTotals)[proc_I];
	}

/*#if DEBUG*/
	Journal_DPrintfL( self->debug, 2, "procCritPoints totals:(");
	for (proc_I = 0; proc_I < meshDecomp->nproc; proc_I++) {
		Journal_DPrintfL( self->debug, 2, "%d, ", (*procCritPointsTotals)[proc_I]);
	}
	Journal_DPrintfL( self->debug, 2, ")\n");
	Journal_DPrintfL( self->debug, 2, "MaxCritPointsPerProc = %d\n", (*maxCritPointsPerProc));
/*#endif*/

	/* Now share the actual point values */
	(*myCriticalPoints) = Memory_Realloc_Array( (*myCriticalPoints), Node_GlobalIndex, *maxCritPointsPerProc );
	(*allCriticalPoints) = Memory_Alloc_Array( Node_GlobalIndex, meshDecomp->nproc * (*maxCritPointsPerProc),
						   "*allCriticalPoints (ShareCriticalPoints)" );

	MPI_Allgather( (*myCriticalPoints), (*maxCritPointsPerProc), MPI_UNSIGNED,
		       (*allCriticalPoints), (*maxCritPointsPerProc), MPI_UNSIGNED,
		       meshDecomp->communicator );

/*#if DEBUG*/
	Journal_DPrintfL( self->debug, 2, "procCritPoints:(" );
	for (proc_I = 0; proc_I < meshDecomp->nproc; proc_I++) {
		for ( point_I = 0; point_I < (*procCritPointsTotals)[proc_I]; point_I++)
		{
			Journal_DPrintfL( self->debug, 2, "%u, ", 
					  (*allCriticalPoints)[(*maxCritPointsPerProc)*proc_I + point_I]);
		}	
	}
	Journal_DPrintfL( self->debug, 2, ")\n");
/*#endif*/
	Stream_UnIndentBranch( StgFEM_Debug );
}



/** Performs an AllGather on a set of CritPointInfo structs: Each processer will afterwards have all the CritPointInfo
    structs of all the processors, in one large array. The user must then index into the array by processor
    carefully. */
static void _FeEquationNumber_ShareCritPointInfo(
	FeEquationNumber* self,
	CritPointInfo** const myCritPointInfo,
	Index myCritPointInfoTotal,
	CritPointInfo** allCritPointInfo,
	Index** procCritPointInfoTotals,
	Index* const maxCritPointInfoPerProc,
	PrintValuesFlag printValuesFlag )
{
	MeshDecomp* meshDecomp = self->feMesh->layout->decomp;
	Partition_Index proc_I; 
/*#if DEBUG*/
	Index point_I = 0;
/*#endif*/

	Journal_DPrintfL( self->debug, 1, "In %s\n",  __func__ );
	Stream_IndentBranch( StgFEM_Debug );
	(*maxCritPointInfoPerProc) = 0;

	/* allgather the totals (to save allocating unnecessary memory) */
	(*procCritPointInfoTotals) = Memory_Alloc_Array( Node_GlobalIndex, meshDecomp->nproc, "*procCritPointInfoTotals" ); 

	MPI_Allgather( &myCritPointInfoTotal, 1, MPI_INT, (*procCritPointInfoTotals), 1, MPI_INT,
		       meshDecomp->communicator );

	for (proc_I = 0; proc_I < meshDecomp->nproc; proc_I++) {
		if ( (*procCritPointInfoTotals)[proc_I] > (*maxCritPointInfoPerProc) )
			(*maxCritPointInfoPerProc) = (*procCritPointInfoTotals)[proc_I];
	}

/*#if DEBUG*/
	if ( Stream_IsPrintableLevel( self->debug, 2 ) ) {
		Journal_DPrintf( self->debug, "procCritPointInfo totals:(" );
		for (proc_I = 0; proc_I < meshDecomp->nproc; proc_I++) {
			Journal_DPrintf( self->debug, "%d, ", (*procCritPointInfoTotals)[proc_I]);
		}
		Journal_DPrintf( self->debug, ")\n");
		Journal_DPrintf( self->debug, "MaxCritPointInfoPerProc = %d\n", (*maxCritPointInfoPerProc));
	}	
/*#endif*/

	/* Now share the actual critPointInfo arrays */
	(*allCritPointInfo) = Memory_Alloc_Array( CritPointInfo, meshDecomp->nproc * (*maxCritPointInfoPerProc), 
						  "allCritPointInfo (ShareCritPointInfo)" );

	/* below changed by PatrickSunter 27/1/2004 to work on grendel */
	MPI_Allgather( (*myCritPointInfo), (*maxCritPointInfoPerProc), MPI_critPointInfoType,
		       (*allCritPointInfo), (*maxCritPointInfoPerProc), MPI_critPointInfoType,
		       meshDecomp->communicator );
	
/*#if DEBUG*/
	if ( Stream_IsPrintableLevel( self->debug, 2 ) ) {	
		Journal_DPrintf( self->debug, "procCritPointInfos" );
		if ( DONT_PRINT_VALUES == printValuesFlag ) {
			Journal_DPrintf( self->debug, "(indices only)" );
		}
		Journal_DPrintf( self->debug, ": (" );
		for (proc_I = 0; proc_I < meshDecomp->nproc; proc_I++) {
			for ( point_I = 0; point_I < (*procCritPointInfoTotals)[proc_I]; point_I++)
			{
				Index current = (*maxCritPointInfoPerProc)*proc_I + point_I;
				Journal_DPrintf( self->debug, "%u", (*allCritPointInfo)[current].index );
				if ( PRINT_VALUES == printValuesFlag ) {
					Journal_DPrintf( self->debug, "= %d", (*allCritPointInfo)[current].eqNum );
				}	
				Journal_DPrintf( self->debug, ", " );
			}	
		}
		Journal_DPrintf( self->debug, ")\n");
	}	
/*#endif*/
	Stream_UnIndentBranch( StgFEM_Debug );
}


/** For each domain node, calculate the Equation number, restarting at 0 whenever a critical node is reached and
    saving the value I was up to. */
static void _FeEquationNumber_DoPartialTotals( FeEquationNumber* self,
					       CritPointInfo* const critPointsIHave, Index critPointsIHaveTotal,
					       CritPointInfo* const critPointsToSend, Index critPointsToSendTotal )
{
	MeshDecomp* meshDecomp = self->feMesh->layout->decomp;
	Node_DomainIndex nodeDomainCount = meshDecomp->nodeDomainCount;
	RemappedNodeInfo_Index rNodeInfo_I = 0;
	Node_RemappedGlobalIndex remappedGlobal_I = 0;
	Node_DomainIndex dNode_I = 0;
	Index haveCritPoint_I = 0;
	Index sendCritPoint_I = 0;
	Dof_EquationNumber currSubTotalEqNum = 0;

	Journal_DPrintf( self->debug, "In %s:\n", __func__ );
	Stream_IndentBranch( StgFEM_Debug );
	
	while ( (rNodeInfo_I < nodeDomainCount) ) {
		remappedGlobal_I = self->remappedNodeInfos[rNodeInfo_I].remappedGlobal;
		dNode_I = self->remappedNodeInfos[rNodeInfo_I].domain;
	
		_FeEquationNumber_HandleNode( self, dNode_I, &currSubTotalEqNum );

		/* if the node is next critical point */
		if ( ( haveCritPoint_I < critPointsIHaveTotal ) 
		     && ( remappedGlobal_I == critPointsIHave[haveCritPoint_I].index ) ) 
		{
			Journal_DPrintfL( self->debug, 3, "storing c.p. subtotal at remapped global index %u = %d\n", 
					  remappedGlobal_I, currSubTotalEqNum );
			/* store current value */
			critPointsIHave[haveCritPoint_I++].eqNum = currSubTotalEqNum;

			/* check if its also a to send point */
			if ( ( sendCritPoint_I < critPointsToSendTotal )
			     && ( remappedGlobal_I == critPointsToSend[sendCritPoint_I].index ) )
			{
                          /* store current value */
				critPointsToSend[sendCritPoint_I++].eqNum = currSubTotalEqNum;
			}
			/* reset counter and move to next point */
			currSubTotalEqNum = 0;
		}	
		rNodeInfo_I++;
	}

/*#if DEBUG*/
	if ( Stream_IsPrintableLevel( self->debug, 2 ) ) {
		Journal_DPrintf( self->debug, "Totals I have:[ " );
		haveCritPoint_I = 0;
		for ( haveCritPoint_I = 0; haveCritPoint_I < critPointsIHaveTotal; haveCritPoint_I++ )
		{	
			Journal_DPrintf( self->debug, "%2u=%d, ", critPointsIHave[haveCritPoint_I].index, critPointsIHave[haveCritPoint_I].eqNum  );
		}	
		Journal_DPrintf( self->debug, "]\n" );
		Journal_DPrintf( self->debug, "Totals to send:[ " );
		sendCritPoint_I = 0;
		for ( sendCritPoint_I = 0; sendCritPoint_I < critPointsToSendTotal; sendCritPoint_I++ )
		{	
			Journal_DPrintf( self->debug, "%2u=%d, ", critPointsToSend[sendCritPoint_I].index, critPointsToSend[sendCritPoint_I].eqNum  );
		}	
		Journal_Printf( self->debug, "]\n" );
	}
/*#endif*/
	Stream_UnIndentBranch( StgFEM_Debug );
}


/** Called right at the end of calculating ID array, to add the locally calculated partial totals with the final
    global information, and arrive at the complete figures. */
static void _FeEquationNumber_AddAllPartialTotals( FeEquationNumber* self, CritPointInfo* mySubTotals,
						   Index myCritPointInfoTotal, CritPointInfo* allSubTotals, Index* procCritPointInfoTotals,
						   Index maxSubTotalsPerProc )
{	
	Partition_Index proc_I;
	FeMesh* mesh = self->feMesh;
	MeshDecomp* meshDecomp = self->feMesh->layout->decomp;
	Index myRank = meshDecomp->rank;
	Node_GlobalIndex nodeDomainCount = meshDecomp->nodeDomainCount;
	RemappedNodeInfo_Index rNodeInfo_I = 0;
	Node_RemappedGlobalIndex remappedGlobal_I = 0;
	Node_DomainIndex dNode_I = 0;
	Node_GlobalIndex gNode_I = 0;
	Index myCritPointInfo_I = 0;
	Index* procCritPointInfo_Is = NULL;
	Dof_EquationNumber adjustSum = 0;
	Dof_EquationNumber currEqNum;
	Dof_Index currNodeNumDofs;
	Dof_Index nodeLocalDof_I;
	AddListEntry* addListStart;
	AddListEntry* addListPtr;
	AddListEntry* addListPrev;
	Bool firstOfSet;
	Bool* haveUpdatedLinkedDofSetTbl = NULL;
	Index linkedDof_I;
	Index linkedDofSet;
#define NEXT_PROCESSOR_SECTION_START ( maxSubTotalsPerProc * (proc_I) + procCritPointInfoTotals[proc_I] )

	Journal_DPrintfL( self->debug, 1, "In %s:\n",  __func__ );
	Stream_IndentBranch( StgFEM_Debug );
	
	gNode_I += 0;
	
	self->_lowestLocalEqNum = -1;
	self->_highestLocalEqNum = -1;
	procCritPointInfo_Is = Memory_Alloc_Array( Index, meshDecomp->nproc, "procCritPointInfo_Is (AddAllPartialTotals)" );
	
	/* set iterators into processors */
	for ( proc_I = 0; proc_I < meshDecomp->nproc; proc_I++ ) {
		procCritPointInfo_Is[proc_I] = maxSubTotalsPerProc * proc_I;
	}
							
	if ( self->linkedDofInfo ) {
		haveUpdatedLinkedDofSetTbl = Memory_Alloc_Array( Bool, self->linkedDofInfo->linkedDofSetsCount,
								 "haveUpdatedLinkedDofSetTbl" );
		for ( linkedDof_I=0; linkedDof_I < self->linkedDofInfo->linkedDofSetsCount; linkedDof_I++ ) {
			haveUpdatedLinkedDofSetTbl[linkedDof_I] = False;
		}	
	}
							

	Journal_DPrintfL( self->debug, 3, "beginning sequential run through nodes...\n" );

	while (rNodeInfo_I < nodeDomainCount) {
		Node_GlobalIndex newCritInfoIndex;
		/* reset the add list information */
		addListPtr = addListStart = NULL;
		remappedGlobal_I = self->remappedNodeInfos[rNodeInfo_I].remappedGlobal;
		dNode_I = self->remappedNodeInfos[rNodeInfo_I].domain;
	
		Journal_DPrintfL( self->debug, 3, "At remapped gNode[%d]: ",  remappedGlobal_I );
		Journal_DPrintfL( self->debug, 3, "...adding critical node totals from the block we just skipped...\n" );

		/* Ok, now add critical node totals from other processors in the block we just skipped. */
		for ( proc_I = 0; proc_I < meshDecomp->nproc; proc_I++ ) {
			if ( proc_I == myRank ) continue;

			/* reset the add list back to start */
			addListPtr = addListStart;
			addListPrev = NULL;
			

			/* Only add critical node totals from other processors lower than current index */
			while (( procCritPointInfo_Is[proc_I] < NEXT_PROCESSOR_SECTION_START ) &&
			       (newCritInfoIndex = allSubTotals[(procCritPointInfo_Is[proc_I])].index) < remappedGlobal_I )
			{
				/* See if we are scheduled to add this critical point yet, and if not where it should go */
				while ( addListPtr && ( newCritInfoIndex > addListPtr->critPointInfo->index ) ) {
					addListPrev = addListPtr; 
					addListPtr = addListPtr->next;
				}
				/* Note: The 2nd condition guards against trying to add a duplicate. In that case we
				 * just progress on. */
				if ( (addListPtr == NULL) || ( addListPtr->critPointInfo->index != newCritInfoIndex )) {
					/*
					  if inside this condition we are either:
					  1. Adding the first crit node to be added
					  2. Have found a non-duplicate crit node to be added, either at the end or
					  partway through.
					  So go ahead and add to the addList, reconnecting both the prev and next.
					*/
					AddListEntry* newAddListEntry = Memory_Alloc( AddListEntry,
										      "newAddListEntry (AddAllPartialTotals)" );
					adjustSum += allSubTotals[procCritPointInfo_Is[proc_I]].eqNum;
					Journal_DPrintfL( self->debug, 3, "at other c.p. (%d):incrementing adjustSum by %d to = %d\n", 
							  allSubTotals[procCritPointInfo_Is[proc_I]].index,
							  allSubTotals[procCritPointInfo_Is[proc_I]].eqNum,
							  adjustSum );
					newAddListEntry->critPointInfo = &allSubTotals[procCritPointInfo_Is[proc_I]];
					newAddListEntry->next = addListPtr;
					if ( addListPrev ) {
						addListPrev->next = newAddListEntry;
					}
					else {
						addListStart = newAddListEntry;
					}
					addListPtr = newAddListEntry;
				}	
				/* Move to the next potential crit node from that processor, regardless of whether we
				 * added the last one or not. */
				(procCritPointInfo_Is[proc_I])++;
			}	
		}

		/* free the add list */
		while ( addListStart ) {
			addListPrev = addListStart;
			addListStart = addListStart->next;
			Memory_Free( addListPrev );
		}
		
		/* Now go through a run of nodes this processor owns, and add subtotals locally calculated */
		firstOfSet = True;

		while ( (firstOfSet) ||
			( (rNodeInfo_I < nodeDomainCount)
			  && (self->remappedNodeInfos[rNodeInfo_I].remappedGlobal ==
			      (self->remappedNodeInfos[rNodeInfo_I-1].remappedGlobal + 1) ) ) )
		{
			remappedGlobal_I = self->remappedNodeInfos[rNodeInfo_I].remappedGlobal;
			dNode_I = self->remappedNodeInfos[rNodeInfo_I].domain;
			gNode_I = mesh->nodeD2G[dNode_I];
			firstOfSet = False;
			
			Journal_DPrintfL( self->debug, 3, "At remapped gNode[%d]: I own it, so adjusting EqNums\n", 
					  remappedGlobal_I );

			/* increment the value at each non-bc dof by adjustSum */
			currNodeNumDofs = self->dofLayout->dofCounts[ dNode_I ];

			for ( nodeLocalDof_I = 0; nodeLocalDof_I < currNodeNumDofs; nodeLocalDof_I++ ) {
				if( !self->bcs || !VariableCondition_IsCondition( self->bcs, dNode_I,
										  self->dofLayout->varIndices[dNode_I][nodeLocalDof_I] ) )
				{
                                  if ( (NULL == self->linkedDofInfo) ||	/* TAG: bcs */
					     ( self->linkedDofInfo->linkedDofTbl[dNode_I][nodeLocalDof_I] == -1 ) )
					{
						/* A normal point, so increment */
						self->destinationArray[dNode_I][nodeLocalDof_I] += adjustSum;
					}
					else {
						/* Handling of linked dofs: the first time we encounter one of a set,
						   add the adjust sum. Afterwards, just set to the updated value. */
						linkedDofSet = self->linkedDofInfo->linkedDofTbl[dNode_I][nodeLocalDof_I];
						if ( False == haveUpdatedLinkedDofSetTbl[linkedDofSet] ) {
							self->destinationArray[dNode_I][nodeLocalDof_I] += adjustSum;
							self->linkedDofInfo->eqNumsOfLinkedDofs[linkedDofSet] += adjustSum;
							haveUpdatedLinkedDofSetTbl[linkedDofSet] = True;
						}
						else {
							self->destinationArray[dNode_I][nodeLocalDof_I] =
								self->linkedDofInfo->eqNumsOfLinkedDofs[linkedDofSet];
						}
					}
						
					currEqNum = self->destinationArray[dNode_I][nodeLocalDof_I];
					/* If current node not a shadow, then check totals */
					if ( dNode_I < mesh->nodeLocalCount ) {
                                          if ( currEqNum != -1 ) {	/* TAG: bcs */
							if ( self->_lowestLocalEqNum == -1 ) {
								/* We need this if statement to initialise this */
								self->_lowestLocalEqNum = currEqNum;
							}	
							else if ( currEqNum < self->_lowestLocalEqNum ) {
								self->_lowestLocalEqNum = currEqNum;
							}

							if ( currEqNum > self->_highestLocalEqNum ) {
								self->_highestLocalEqNum = currEqNum;
							}	
						}
					}
				}	
			}		

			/* if the node is my next critical point, add the sub-total */
			if ( ( myCritPointInfo_I < myCritPointInfoTotal ) &&
			     ( remappedGlobal_I == mySubTotals[myCritPointInfo_I].index)  ) 
			{
				adjustSum += mySubTotals[myCritPointInfo_I].eqNum;
				Journal_DPrintfL( self->debug, 3, "at   my    c.p. %d(%d):incrementing adjustSum by %d to = %d\n", 
						  myCritPointInfo_I, remappedGlobal_I, mySubTotals[myCritPointInfo_I].eqNum, adjustSum );
				myCritPointInfo_I++;
			}

			rNodeInfo_I++;
		}	

		/* If any of the other processors had sub-totals I've just added, jump past them */
		if (rNodeInfo_I < nodeDomainCount) {
			for ( proc_I = 0; proc_I < meshDecomp->nproc; proc_I++ ) {

				while (( procCritPointInfo_Is[proc_I] < NEXT_PROCESSOR_SECTION_START ) &&
				       allSubTotals[(procCritPointInfo_Is[proc_I])].index <= remappedGlobal_I )
				{
					(procCritPointInfo_Is[proc_I])++;
				}

			}
		}
	}

	if ( self->linkedDofInfo ) {
		Memory_Free( haveUpdatedLinkedDofSetTbl );
	}	
	Memory_Free( procCritPointInfo_Is );
	Stream_UnIndentBranch( StgFEM_Debug );
}


/** Updates the currEqNum argument for the number of dofs at the current node, and increments the ID Eq. Num for each
    dof where appropriate (setting it to -1 if its a B.C.) */
static void _FeEquationNumber_HandleNode( FeEquationNumber* self, const Node_DomainIndex dNode_I,
					  Dof_EquationNumber* const currEqNum )
{
	Dof_Index	currNodeNumDofs;
	Dof_Index	nodeLocalDof_I;
/*#if DEBUG*/
	Stream*		error;
/*#endif*/
	
	Journal_DPrintfL( self->debug, 3, "In %s:",  __func__ );
/*#if DEBUG*/
	error = Journal_Register( Error_Type, self->type );
	Journal_Firewall( (dNode_I < self->feMesh->nodeDomainCount), error, "Stuffup: %s() asked to operate on node %d, bigger "
			  "than domain count %d. Exiting.\n", __func__, dNode_I, self->feMesh->nodeDomainCount );
/*#endif*/
	
	/* get the number of dofs */
	currNodeNumDofs = self->dofLayout->dofCounts[ dNode_I ];
	Journal_DPrintfL( self->debug, 3, " domain node %d, has %d dofs.\n", dNode_I, currNodeNumDofs );
	/* process each dof, at the same time updating the updatePtr for the next node. */
	for ( nodeLocalDof_I = 0; nodeLocalDof_I < currNodeNumDofs; nodeLocalDof_I++ ) {
		Journal_DPrintfL( self->debug, 3, "dof %d: ", nodeLocalDof_I );

		/* if dof is a boundary condition, set eq num=-1. Otherwise increment counter. */
		/* Also, only set to -1 if we want the BCs removed from the matrix. - Luke */
		if( self->removeBCs && 
		    self->bcs && VariableCondition_IsCondition( self->bcs, dNode_I, self->dofLayout->varIndices[dNode_I][nodeLocalDof_I] ) ) 
		{
			Journal_DPrintfL( self->debug, 3, "is a BC, setting ID[%d][%d] = -1.\n", dNode_I, nodeLocalDof_I );
			self->destinationArray[dNode_I][nodeLocalDof_I] = -1;
		}	
		/* Check if node,dof is a linked dof */
		else if ( self->linkedDofInfo && ( self->linkedDofInfo->linkedDofTbl[dNode_I][nodeLocalDof_I] != -1 ) ) {
			Index linkedDofSet = self->linkedDofInfo->linkedDofTbl[dNode_I][nodeLocalDof_I];
			
			if ( self->linkedDofInfo->eqNumsOfLinkedDofs[linkedDofSet] != -1 ) {
				Journal_DPrintfL( self->debug, 3, "is a linked Dof, so setting ID[%d][%d] to the "
						  "previous value = %d.\n", dNode_I, nodeLocalDof_I, 
						  self->linkedDofInfo->eqNumsOfLinkedDofs[linkedDofSet] );
				/* value is same as the first part of the linked node set */
				self->destinationArray[dNode_I][nodeLocalDof_I] =
					self->linkedDofInfo->eqNumsOfLinkedDofs[linkedDofSet];
			}
			else /* This is the first time we've hit this linked eq num, so treat as normal */
			{
				Journal_DPrintfL( self->debug, 3, "is a linked Dof hit the first time, so setting "
						  "ID[%d][%d] = %d.\n", dNode_I, nodeLocalDof_I, *currEqNum );
			
				self->destinationArray[dNode_I][nodeLocalDof_I] = (*currEqNum)++;
				if ( dNode_I < self->feMesh->nodeLocalCount ) {
					self->_highestLocalEqNum = self->destinationArray[dNode_I][nodeLocalDof_I];
				}	
				/* And save the value at the linked eq num for later */
				self->linkedDofInfo->eqNumsOfLinkedDofs[linkedDofSet] = 
					self->destinationArray[dNode_I][nodeLocalDof_I];
			}
		}
		else {	
			Journal_DPrintfL( self->debug, 3, "is a normal node, setting ID[%d][%d] = %d.\n", dNode_I, nodeLocalDof_I, *currEqNum );
			self->destinationArray[dNode_I][nodeLocalDof_I] = (*currEqNum)++;
			if ( dNode_I < self->feMesh->nodeLocalCount ) {
				self->_highestLocalEqNum = self->destinationArray[dNode_I][nodeLocalDof_I];
			}	
		}
	}
}


void _FeEquationNumber_PostProcessLinkedDofs( FeEquationNumber* self ) {
	Index			linkedDof_I;
	int			valueToReduce;
	int			minimumValue;
	Node_Index		rNodeInfo_I;
	Node_DomainIndex	dNode_I;
	MeshDecomp*		meshDecomp = self->feMesh->layout->decomp;
	Dof_Index		nodeLocalDof_I;
	Dof_Index		currNodeNumDofs;
	Index			linkedDofSet;
	Bool*			adjustDueToNonLocalLinkedDofsTbl = NULL;
	unsigned int		subtractionAdjustmentTotal;
	Dof_EquationNumber      currEqNum;
	
	Journal_DPrintfL( self->debug, 1, "In Func %s():\n", __func__ );
	adjustDueToNonLocalLinkedDofsTbl = Memory_Alloc_Array( Bool, self->linkedDofInfo->linkedDofSetsCount,
							       "adjustDueToNonLocalLinkedDofsTbl" );
		
	/* We need to re-calculate these, since the post processing may alter both the lowest and higest eqNum values */
	self->_lowestLocalEqNum = -1;
	self->_highestLocalEqNum = -1;
		
	Journal_DPrintfL( self->debug, 2, "Note: in reductions to follow, value of %d, the node global count, "
			  "denotes that this processor doesn't hold any of that linked dof.\n", self->feMesh->nodeGlobalCount );

	Stream_Indent( self->debug );
	for ( linkedDof_I=0; linkedDof_I < self->linkedDofInfo->linkedDofSetsCount; linkedDof_I++ ) {
		if ( -1 != self->linkedDofInfo->eqNumsOfLinkedDofs[linkedDof_I] ) {
		 	valueToReduce = self->linkedDofInfo->eqNumsOfLinkedDofs[linkedDof_I];
		}
		else {
			valueToReduce = self->feMesh->nodeGlobalCount;
		}	
		
		Journal_DPrintfL( self->debug, 2, "Reducing linked Dof %d: this proc has value %d\n", linkedDof_I, valueToReduce );
		MPI_Allreduce( &valueToReduce, &minimumValue, 1, MPI_INT, MPI_MIN, meshDecomp->communicator );
		
		Journal_DPrintfL( self->debug, 2, "\tMinimum val = %d\n", minimumValue );
		self->linkedDofInfo->eqNumsOfLinkedDofs[linkedDof_I] = minimumValue;

		if ( valueToReduce != minimumValue ) {
			adjustDueToNonLocalLinkedDofsTbl[linkedDof_I] = True;
		}
		else {
			adjustDueToNonLocalLinkedDofsTbl[linkedDof_I] = False;
		}
	}	
	Stream_UnIndent( self->debug );
	
	Journal_DPrintfL( self->debug, 2, "Post-processing totals to correct based on non-local linked dofs\n" );
	subtractionAdjustmentTotal = 0;
	
	Stream_Indent( self->debug );
	for ( rNodeInfo_I = 0; rNodeInfo_I < self->feMesh->nodeDomainCount; rNodeInfo_I++ ) {
		dNode_I = self->remappedNodeInfos[rNodeInfo_I].domain;
		Journal_DPrintfL( self->debug, 3, "At remapped node %d (dNode %d):\n", rNodeInfo_I, dNode_I );
		currNodeNumDofs = self->dofLayout->dofCounts[ dNode_I ];
		
		for ( nodeLocalDof_I = 0; nodeLocalDof_I < currNodeNumDofs; nodeLocalDof_I++ ) {

			Stream_Indent( self->debug );
			Journal_DPrintfL( self->debug, 3, "At dof %d: ", nodeLocalDof_I );
		
			if( self->bcs && VariableCondition_IsCondition( self->bcs, dNode_I,
									self->dofLayout->varIndices[dNode_I][nodeLocalDof_I] ) )
			{
				Journal_DPrintfL( self->debug, 3, "is a BC: ignoring.\n" );
			}
			else if ( self->linkedDofInfo->linkedDofTbl[dNode_I][nodeLocalDof_I] != -1 ) {
				linkedDofSet = self->linkedDofInfo->linkedDofTbl[dNode_I][nodeLocalDof_I];

				Journal_DPrintfL( self->debug, 3, "is a linked Dof, so setting value to %d\n",
						  self->linkedDofInfo->eqNumsOfLinkedDofs[linkedDofSet] );

				self->destinationArray[dNode_I][nodeLocalDof_I] -= subtractionAdjustmentTotal;

				if ( True == adjustDueToNonLocalLinkedDofsTbl[linkedDofSet] ) {
					/* We now need to test, after the subtraction adjustment has been made, if the 
					   value of this linked dof now matches the minimum. If so, don't subtract further. */
					if ( self->destinationArray[dNode_I][nodeLocalDof_I] != 
					     self->linkedDofInfo->eqNumsOfLinkedDofs[linkedDofSet] )
					{
						subtractionAdjustmentTotal++;
						Journal_DPrintfL( self->debug, 3, "And since first was non-local, increasing "
								  "subtraction adj. to %d\n", subtractionAdjustmentTotal );
					}	
						
					/* Make sure we only adjust once */
					adjustDueToNonLocalLinkedDofsTbl[linkedDofSet] = False;
				}

				self->destinationArray[dNode_I][nodeLocalDof_I] =
					self->linkedDofInfo->eqNumsOfLinkedDofs[linkedDofSet];
				
			}
			else {
				Journal_DPrintfL( self->debug, 3, "subtracting %d\n", subtractionAdjustmentTotal );
				self->destinationArray[dNode_I][nodeLocalDof_I] -= subtractionAdjustmentTotal;
			}

			currEqNum = self->destinationArray[dNode_I][nodeLocalDof_I];
			if ( dNode_I < self->feMesh->nodeLocalCount ) {
                          if ( currEqNum != -1 ) {	/* TAG: bcs */
					if ( self->_lowestLocalEqNum == -1 ) {
						/* We need this if statement to initialise this */
						self->_lowestLocalEqNum = currEqNum;
					}	
					else if ( currEqNum < self->_lowestLocalEqNum ) {
						self->_lowestLocalEqNum = currEqNum;
					}

					if ( currEqNum > self->_highestLocalEqNum ) {
						self->_highestLocalEqNum = currEqNum;
					}	
				}
			}	
			Stream_UnIndent( self->debug );
		}	
	}		
	Stream_UnIndent( self->debug );

	Memory_Free( adjustDueToNonLocalLinkedDofsTbl );
}
#endif


Index FeEquationNumber_CalculateActiveEqCountAtNode(
	void*			feEquationNumber,
	Node_DomainIndex	dNode_I,
	Dof_EquationNumber*	lowestActiveEqNumAtNodePtr )
{
	FeEquationNumber*	self = (FeEquationNumber*) feEquationNumber;
	Dof_Index		nodalDof_I = 0;
	Index			activeEqsAtCurrRowNode = 0;		
	Dof_EquationNumber	currEqNum;
	Bool			foundLowest = False;

	for ( nodalDof_I = 0; nodalDof_I < self->dofLayout->dofCounts[dNode_I]; nodalDof_I++ ) {
		currEqNum = self->destinationArray[dNode_I][nodalDof_I];
		if ( currEqNum != -1 ) {
			activeEqsAtCurrRowNode++;
			if ( False == foundLowest ) {
				(*lowestActiveEqNumAtNodePtr) = currEqNum;
				foundLowest = True;
			}
		}
	}

	return activeEqsAtCurrRowNode;
}


void FeEquationNumber_BuildLocationMatrix( void* feEquationNumber ) {
	FeEquationNumber*	self = (FeEquationNumber*)feEquationNumber;
	FeMesh*			feMesh;
	unsigned		nDims;
	unsigned		nDomainEls;
	unsigned		nLocalNodes;
	unsigned*		nNodalDofs;
	unsigned		nElNodes;
	unsigned*		elNodes;
	int**			dstArray;
	int***			locMat;
	unsigned		e_i, n_i, dof_i;

	assert( self );

	/* Don't build if already done. */
	if( self->locationMatrixBuilt ) {
		Journal_DPrintf( self->debugLM, "In %s: LM already built, so just returning.\n",  __func__ );
		Stream_UnIndentBranch( StgFEM_Debug );
		return;
	}

	/* Shortcuts. */
	feMesh = self->feMesh;
	nDims = Mesh_GetDimSize( feMesh );
	nDomainEls = FeMesh_GetElementDomainSize( feMesh );
	nLocalNodes = FeMesh_GetNodeLocalSize( feMesh );
	nNodalDofs = self->dofLayout->dofCounts;
	dstArray = self->destinationArray;

	/* Allocate for the location matrix. */
	locMat = AllocArray( int**, nDomainEls );
	for( e_i = 0; e_i < nDomainEls; e_i++ ) {
		FeMesh_GetElementNodes( feMesh, e_i, &nElNodes, &elNodes );
		locMat[e_i] = AllocArray( int*, nElNodes );
		for( n_i = 0; n_i < nElNodes; n_i++ )
			locMat[e_i][n_i] = AllocArray( int, nNodalDofs[elNodes[n_i]] );
	}

	/* Build location matrix. */
	for( e_i = 0; e_i < nDomainEls; e_i++ ) {
		FeMesh_GetElementNodes( feMesh, e_i, &nElNodes, &elNodes );
		for( n_i = 0; n_i < nElNodes; n_i++ ) {
			for( dof_i = 0; dof_i < nNodalDofs[elNodes[n_i]]; dof_i++ )
				locMat[e_i][n_i][dof_i] = dstArray[elNodes[n_i]][dof_i];
		}
	}

	/* Store result. */
	self->locationMatrix = locMat;
}


#if 0
/** build the element location matrix mapping elements, element node, dof -> eq num */
void FeEquationNumber_BuildLocationMatrix( FeEquationNumber* self ) {
	FeMesh* mesh = self->feMesh;
	Element_LocalIndex lElement_I;
	Node_LocalIndex numNodesThisElement = 0;
	Node_LocalIndex elLocalNode_I = 0;
	Dof_Index numDofsThisNode = 0;
	Dof_Index** dofCountsAtElementNodesArray = NULL;
	Element_LocalIndex elementLocalCount = Mesh_GetLocalSize( mesh, Mesh_GetDimSize( mesh ) );

	Journal_DPrintf( self->debug, "In %s():\n", __func__ );
	Stream_IndentBranch( StgFEM_Debug );
	
	if (self->locationMatrixBuilt) {
		Journal_DPrintf( self->debugLM, "In %s: LM already built, so just returning.\n",  __func__ );
		Stream_UnIndentBranch( StgFEM_Debug );
		return;
	}

	Journal_DPrintf( self->debugLM, "In %s: building over %d elements.\n",  __func__, elementLocalCount );

	/* Allocate the LM 3D array using the Memory module, in 2 stage process */
	dofCountsAtElementNodesArray = Memory_Alloc_3DSetup( elementLocalCount, 
							     mesh->topo->nIncEls[nDims][MT_VERTEX] );

	for ( lElement_I = 0; lElement_I < elementLocalCount; lElement_I++ ) {
		numNodesThisElement = mesh->elementNodeCountTbl[lElement_I];

		for( elLocalNode_I = 0; elLocalNode_I < numNodesThisElement; elLocalNode_I++) {
			Node_LocalIndex dNode_I = mesh->elementNodeTbl[lElement_I][elLocalNode_I];
			numDofsThisNode = self->dofLayout->dofCounts[dNode_I];
			dofCountsAtElementNodesArray[lElement_I][elLocalNode_I] = numDofsThisNode;
		}
	}

	self->locationMatrix = Memory_Alloc_3DComplex( Dof_EquationNumber, elementLocalCount, mesh->elementNodeCountTbl,
						       dofCountsAtElementNodesArray, "FeEquationNumber->locationMatrix" );
	/* Free the dof counts array:- we have to look up domain node numbers anyway later, might as
	   well just use the dof counts array then. */
	Memory_Free( dofCountsAtElementNodesArray );

	for ( lElement_I = 0; lElement_I < elementLocalCount; lElement_I++ ) {
		FeEquationNumber_BuildOneElementLocationMatrix( self, lElement_I );
	}	

	self->locationMatrixBuilt = True;
	Stream_UnIndentBranch( StgFEM_Debug );
}
#endif


/** Build an element's local location matrix */
Dof_EquationNumber** FeEquationNumber_BuildOneElementLocationMatrix( void* feEquationNumber, Element_LocalIndex lElement_I ) {
	FeEquationNumber* self = (FeEquationNumber*)feEquationNumber;
	Node_DomainIndex elLocalNode_I;
	Node_DomainIndex numNodesThisElement, *elInc;
	Dof_EquationNumber** localLocationMatrix = NULL;
	FeMesh* feMesh = self->feMesh;
	Dof_Index numDofsThisNode = 0;

	FeMesh_GetElementNodes( feMesh, lElement_I, &numNodesThisElement, &elInc );

	/* HACK: Make sure global element location matrix is built. */
	if( !self->locationMatrixBuilt )
		abort();

	/* if ( big LM allocated ) set pointer into it correctly */
	if ( self->locationMatrix ) {
		/* set ptr to correct set of local nodes ptrs */
		localLocationMatrix = self->locationMatrix[lElement_I];
		Journal_DPrintfL( self->debugLM, 3, "set localLocationMatrix to pt. to big LM[%d] = %p\n",  lElement_I, self->locationMatrix[lElement_I] ) ;
	}

#if 0
	else {
		Dof_Index* numDofsEachNode = NULL;

		/* allocate memory for local LM to return */
		numDofsEachNode = Memory_Alloc_Array_Unnamed( Dof_Index, numNodesThisElement );

		for( elLocalNode_I = 0; elLocalNode_I < numNodesThisElement; elLocalNode_I++) {
			Node_LocalIndex localNode = mesh->elementNodeTbl[lElement_I][elLocalNode_I];
			numDofsEachNode[elLocalNode_I] = self->dofLayout->dofCounts[localNode];
		}

		localLocationMatrix = Memory_Alloc_2DComplex( Dof_EquationNumber, numNodesThisElement, numDofsEachNode,
							      "localLocationMatrix (set of ptrs to dof lists, indexed by element-local node)" );
		Memory_Free( numDofsEachNode );
	}
#endif

#if 0
	/* If we haven't yet built full LM, copy ID values across */
	if ( False == self->locationMatrixBuilt ) {
		/* for (each el-local node) */
		for ( elLocalNode_I = 0; elLocalNode_I < numNodesThisElement; elLocalNode_I++ ) {
			/* look up processor local node number. */
			Node_LocalIndex procDomainNode = mesh->elementNodeTbl[lElement_I][elLocalNode_I];
			numDofsThisNode = self->dofLayout->dofCounts[procDomainNode];

			/* copy pointers to dof eq nums from ID array relevant to that node */
			Journal_DPrintfL( self->debugLM, 3, "copying %d dof eq. numbers from ID[%d] to LM[%d][%d]\n", 
					  numDofsThisNode, procDomainNode, lElement_I, elLocalNode_I );
			memcpy( localLocationMatrix[elLocalNode_I],
				self->destinationArray[procDomainNode], numDofsThisNode * sizeof(Dof_EquationNumber) );
		}
	}
#endif

	return localLocationMatrix;
}


void FeEquationNumber_PrintDestinationArray( void* feFeEquationNumber, Stream* stream ) {
	FeEquationNumber* self = (FeEquationNumber*) feFeEquationNumber;
	FeMesh*		feMesh = self->feMesh;
	MPI_Comm comm = CommTopology_GetComm( Mesh_GetCommTopology( feMesh, MT_VERTEX ) );
	unsigned rank;
	Node_GlobalIndex gNode_I;
	Node_GlobalIndex nodeGlobalCount = FeMesh_GetNodeGlobalSize( feMesh );

	MPI_Comm_rank( comm, (int*)&rank );
	Journal_Printf( stream, "%d: *** Printing destination array ***\n", rank );

	for (gNode_I =0; gNode_I < nodeGlobalCount; gNode_I++) {
		Node_DomainIndex dNode_I;

		if ( !Mesh_GlobalToDomain( feMesh, MT_VERTEX, gNode_I, &dNode_I ) ) {
			Journal_Printf( stream, "\tdestinationArray[(gnode)%2d]: on another proc\n", gNode_I);
		}
		else {
			Dof_Index currNodeNumDofs = self->dofLayout->dofCounts[ dNode_I ];
			Dof_Index nodeLocalDof_I;

			Journal_Printf( stream, "\tdestinationArray[(gnode)%2d][(dof)0-%d]:",gNode_I, currNodeNumDofs );
			for( nodeLocalDof_I = 0; nodeLocalDof_I < currNodeNumDofs; nodeLocalDof_I++ ) {
				Journal_Printf( stream, "%3d, ", self->destinationArray[dNode_I][nodeLocalDof_I] );
			}
			Journal_Printf( stream, "\n" );
		}
	}		
}


void FeEquationNumber_PrintLocationMatrix( void* feFeEquationNumber, Stream* stream ) {
	FeEquationNumber* self = (FeEquationNumber*) feFeEquationNumber;
	FeMesh*		feMesh = self->feMesh;
	MPI_Comm comm = CommTopology_GetComm( Mesh_GetCommTopology( feMesh, MT_VERTEX ) );
	unsigned rank;
	Element_GlobalIndex gEl_I;
	unsigned nDims = Mesh_GetDimSize( feMesh );
	Element_GlobalIndex elementGlobalCount = FeMesh_GetElementGlobalSize( feMesh );
	unsigned nLocalEls = FeMesh_GetElementLocalSize( feMesh );

	Journal_Printf( stream, "%d: *** Printing location matrix ***\n", rank  );

	MPI_Comm_rank( comm, (int*)&rank );
	
	for (gEl_I =0; gEl_I < elementGlobalCount; gEl_I++ ) {
		Element_LocalIndex lEl_I;

		if ( !Mesh_GlobalToDomain( feMesh, nDims, gEl_I, &lEl_I ) || lEl_I >= nLocalEls ) {
			Journal_Printf( stream, "\tLM[(g/l el)%2d/XXX]: on another proc\n", gEl_I);
		}
		else {
			Node_LocalIndex numNodesAtElement;
			Node_LocalIndex elLocalNode_I;
			unsigned*	incNodes;

			FeMesh_GetElementNodes( self->feMesh, lEl_I, &numNodesAtElement, &incNodes );

			Journal_Printf( stream, "\tLM[(g/l el)%2d/%2d][(enodes)0-%d]", gEl_I, lEl_I, numNodesAtElement);	
			/* print the nodes and dofs */
			for ( elLocalNode_I = 0; elLocalNode_I < numNodesAtElement; elLocalNode_I++ ) {
				/* look up processor local node number. */
				Element_LocalIndex currNode = incNodes[elLocalNode_I == 2 ? 3 : 
								       elLocalNode_I == 3 ? 2 : 
								       elLocalNode_I == 6 ? 7 : 
								       elLocalNode_I == 7 ? 6 : 
								       elLocalNode_I];
				/* get the number of dofs at current node */
				Dof_Index currNodeNumDofs = self->dofLayout->dofCounts[ currNode ];
				Dof_Index nodeLocalDof_I;

				Journal_Printf( stream, "({%2d}", currNode );
				for( nodeLocalDof_I = 0; nodeLocalDof_I < currNodeNumDofs; nodeLocalDof_I++ ) {
					Journal_Printf( stream, "%3d,", self->destinationArray[currNode][nodeLocalDof_I] );
				}
				Journal_Printf( stream, "), " );
			}	

			Journal_Printf( stream, "\n" );
		}

	}
}


#if 0
void FeEquationNumber_PrintElementLocationMatrix(
	void*			feEquationNumber,
	Dof_EquationNumber**	elementLM,
	Element_LocalIndex	element_lI,
	Stream*			stream )
{
	FeEquationNumber*	self = (FeEquationNumber*)feEquationNumber;
	Dof_Index		dof_elLocalI;
	Node_ElementLocalIndex	nodeCountThisEl = self->feMesh->elementNodeCountTbl[element_lI];
	Node_LocalIndex		node_lI = self->feMesh->elementNodeTbl[element_lI][nodeCountThisEl-1];
	Dof_Index		dofCountLastNode = self->dofLayout->dofCounts[node_lI];
	Dof_Index		totalDofsThisElement = 0;
	
	totalDofsThisElement = &elementLM[nodeCountThisEl-1][dofCountLastNode-1] - &elementLM[0][0] + 1;

	Journal_DPrintf( stream, "LM[ el %d ], dofs[0-%d] = {", element_lI, totalDofsThisElement );

	for( dof_elLocalI=0; dof_elLocalI < totalDofsThisElement; dof_elLocalI++ ) {
		Journal_DPrintf( stream, "%d, ", elementLM[0][dof_elLocalI] ); 
	}
	Journal_DPrintf( stream, "}\n" );
}


/* Calculates global sum unconstrained dofs */
void _FeEquationNumber_CalculateGlobalUnconstrainedDofTotal( FeEquationNumber* self ) {
	int	globalSumUnconstrainedDofs;
	MeshDecomp* meshDecomp = self->feMesh->layout->decomp;
	
	Journal_DPrintfL( self->debug, 1, "In %s:\n",  __func__ );
	Stream_IndentBranch( StgFEM_Debug );
	MPI_Allreduce( &self->_highestLocalEqNum, &globalSumUnconstrainedDofs, 1, MPI_INT, MPI_MAX, meshDecomp->communicator );
	self->globalSumUnconstrainedDofs = (unsigned)(globalSumUnconstrainedDofs+1);
	
	Journal_DPrintf( self->debug, "Calculated total (across all processors) unconstrained dofs as:%d\n", self->globalSumUnconstrainedDofs  );
	Stream_UnIndentBranch( StgFEM_Debug );
}


/* calculate the minimum and maximum parts that my processor is responsible for holding */
void _FeEquationNumber_CalculateEqNumsDecomposition( FeEquationNumber* self ) {
	MeshDecomp*		meshDecomp = self->feMesh->layout->decomp;
	Partition_Index		myRank = meshDecomp->rank;
	Partition_Index		nProc = meshDecomp->nproc;

	Journal_DPrintfL( self->debug, 1, "In %s:\n",  __func__ );
	Stream_IndentBranch( StgFEM_Debug );

	if ( (self->remappingActivated) && ( (self->linkedDofInfo == NULL) || (self->linkedDofInfo->linkedDofSetsCount == 0 ) ) ) {
		/* If the remapping is activated, and things aren't complicated by periodic BCs,
		   then each processor should hold the Matrix/Vector
		   component corresponding to the lowest local eqNum, to the lowest eqNum of the next
		   processor. This means that _only shared boundary nodes_ will need to be communicated,
		   and the last processor will have no communication. */
		Journal_DPrintfL( self->debug, 2, "Remapping active and no periodic bcs: using lowest local eqNums as boundaries.\n");

		self->_lowestGlobalEqNums = Memory_Alloc_Array( Dof_EquationNumber, nProc,
								"FeEquationNumber->_lowestGlobalEqNums" );
		MPI_Allgather(
			&self->_lowestLocalEqNum, 1, MPI_INT,
			self->_lowestGlobalEqNums, 1, MPI_INT,
			meshDecomp->communicator );

		self->firstOwnedEqNum = self->_lowestLocalEqNum;
		if (myRank == nProc-1) {
			self->lastOwnedEqNum = self->_highestLocalEqNum;
		}
		else {
			Node_LocalIndex nextProcLowestEqNum = self->_lowestGlobalEqNums[myRank+1] - 1;
			if ( (unsigned int)-1 == nextProcLowestEqNum ) {
				/* Pathological case of next proc having all B.C.s */
				self->lastOwnedEqNum = self->_highestLocalEqNum;
			}	
			else {
				self->lastOwnedEqNum = nextProcLowestEqNum;
			}
		}

		self->localEqNumsOwnedCount = self->lastOwnedEqNum - self->firstOwnedEqNum + 1;	
	}	
	else {
		/* If the remapping isn't activated and the eqNum ordering isn't aligned with the mesh
		   decomposition, or there are periodic BCs, we can't get a clear idea of where the processor boundaries are:
		   therefore just split up the EqNums equally between processors: still should be 
		   fairly good alignment. */
		Journal_DPrintfL( self->debug, 2, "Remapping inactive and/or periodic bcs used: just dividing eqNums as evenly as possible.\n");

		self->_eqNumsPerProcDivisor = self->globalSumUnconstrainedDofs / nProc;
		self->_eqNumsRemainder = self->globalSumUnconstrainedDofs % nProc;
		Journal_DPrintfL( self->debug, 2, "Calculated %d eqNums per proc, with %d remainder\n",
				  self->_eqNumsPerProcDivisor, self->_eqNumsRemainder );
		self->_remNotAddedChangeover = (self->_eqNumsPerProcDivisor+1) * self->_eqNumsRemainder;

		self->localEqNumsOwnedCount = self->_eqNumsPerProcDivisor;
		if ( myRank < self->_eqNumsRemainder ) {
			self->localEqNumsOwnedCount++;
		}

		if ( myRank < self->_eqNumsRemainder ) {
			self->firstOwnedEqNum = myRank * (self->_eqNumsPerProcDivisor+1);
		}
		else {
			self->firstOwnedEqNum = self->_remNotAddedChangeover
				+ (myRank - self->_eqNumsRemainder) * self->_eqNumsPerProcDivisor;
		}
		self->lastOwnedEqNum = self->firstOwnedEqNum + self->localEqNumsOwnedCount - 1;
	}	

	Journal_DPrintfL( self->debug, 1, "Calculated I own %d eqNums, between indices %d to %d\n",
			  self->localEqNumsOwnedCount, self->firstOwnedEqNum, self->lastOwnedEqNum );
	Journal_DPrintfL( self->debug, 1, "(Range of eqNums on local mesh segment is %d to %d)\n",
			  self->_lowestLocalEqNum, self->_highestLocalEqNum );

	Stream_UnIndentBranch( StgFEM_Debug );
}
#endif


Partition_Index FeEquationNumber_CalculateOwningProcessorOfEqNum( void* feEquationNumber, Dof_EquationNumber eqNum ) {
	FeEquationNumber* self = (FeEquationNumber*)feEquationNumber;
	Partition_Index ownerProc = (unsigned int)-1;
	CommTopology*	commTopo = Mesh_GetCommTopology( self->feMesh, MT_VERTEX );
	MPI_Comm	comm = CommTopology_GetComm( commTopo );
	unsigned	nProcs;
	unsigned	p_i;

	MPI_Comm_size( comm, (int*)&nProcs );
	for( p_i = 1; p_i < nProcs; p_i++ ) {
		if( eqNum < self->_lowestGlobalEqNums[p_i] )
			break;
	}

	return p_i - 1;

#if 0
		if ( (self->remappingActivated) && ( (self->linkedDofInfo == NULL) || (self->linkedDofInfo->linkedDofSetsCount == 0 ) ) ) {
			MeshDecomp*		meshDecomp = self->feMesh->layout->decomp;
			Partition_Index		myRank = meshDecomp->rank;
			Partition_Index		nProc = meshDecomp->nproc;

			/* Expect it to be on the next processor, so try there first */
			if ( eqNum > self->lastOwnedEqNum ) {
				ownerProc = myRank + 1;
				while ( (ownerProc+1) < nProc ) {
					if ( eqNum >= self->_lowestGlobalEqNums[ownerProc+1] ) {
						ownerProc++;
					}	
					else {
						break;
					}	
				}
			}
			/* otherwise count back from current */
			else {
				ownerProc = myRank;
				while ( ownerProc > 0 ) {
					if ( eqNum < self->_lowestGlobalEqNums[ownerProc] ) {
						ownerProc--;
					}
					else {
						break;
					}	
				}
			}
		}
		else {
			if ( eqNum < self->_remNotAddedChangeover ) {
				ownerProc = eqNum / (self->_eqNumsPerProcDivisor+1);
			}
			else {
				ownerProc = self->_eqNumsRemainder + (eqNum - self->_remNotAddedChangeover) / self->_eqNumsPerProcDivisor;
			}
		}
	}

	return ownerProc;
#endif
}	


#if 0
void FeEquationNumber_Create_CritPointInfo_MPI_Datatype( void ) {
#define CRIT_POINT_INFO_NBLOCKS 2
	MPI_Aint indexExtent = 0;
	MPI_Datatype critPointInfoTypes[CRIT_POINT_INFO_NBLOCKS] = {MPI_UNSIGNED, MPI_INT };
	MPI_Aint critPointInfoBlockDisplacements[CRIT_POINT_INFO_NBLOCKS];
	int critPointInfoBlockLengths[CRIT_POINT_INFO_NBLOCKS] = { 1, 1 };

	MPI_Type_extent(MPI_UNSIGNED, &indexExtent);
	critPointInfoBlockDisplacements[0] = 0;
	critPointInfoBlockDisplacements[1] = indexExtent;

	MPI_Type_struct( CRIT_POINT_INFO_NBLOCKS, critPointInfoBlockLengths, critPointInfoBlockDisplacements,
			 critPointInfoTypes, &MPI_critPointInfoType );
	
	MPI_Type_commit( &MPI_critPointInfoType );
}
#endif


void FeEquationNumber_BuildWithTopology( FeEquationNumber* self ) {
	Stream*			stream;
	double			startTime, endTime;
	FeMesh*			feMesh;
	Decomp_Sync*		sync;
	CommTopology*		commTopo;
	MPI_Comm		comm;
	unsigned		rank, nProcs;
	unsigned		nDims;
	unsigned		nDomainNodes, nDomainEls;
	unsigned		nLocalNodes;
	unsigned*		nNodalDofs;
	unsigned		nElNodes, *elNodes;
	int**			dstArray;
	int			*nLocMatDofs, ***locMat;
	unsigned		varInd;
	unsigned		curEqNum;
	unsigned		base;
	unsigned		subTotal;
	MPI_Status		status;
	unsigned		maxDofs;
	unsigned*		tuples;
	Decomp_Sync_Array*	array;
	LinkedDofInfo*		links;
	unsigned		e_i, n_i, dof_i;

	assert( self );

	stream = Journal_Register( Info_Type, self->type );
	Journal_Printf( stream, "FeEquationNumber: '%s'\n", self->name );
	Stream_Indent( stream );
	Journal_Printf( stream, "Generating equation numbers...\n" );
	Stream_Indent( stream );
	if( self->removeBCs )
		Journal_Printf( stream, "BCs set to be removed.\n" );
	else
		Journal_Printf( stream, "BCs will not be removed.\n" );

	startTime = MPI_Wtime();

	/* Shortcuts. */
	feMesh = self->feMesh;
	commTopo = Mesh_GetCommTopology( feMesh, MT_VERTEX );
	comm = CommTopology_GetComm( commTopo );
	MPI_Comm_size( comm, (int*)&nProcs );
	MPI_Comm_rank( comm, (int*)&rank );
	nDims = Mesh_GetDimSize( feMesh );
	nDomainNodes = FeMesh_GetNodeDomainSize( feMesh );
	nDomainEls = FeMesh_GetElementDomainSize( feMesh );
	nLocalNodes = FeMesh_GetNodeLocalSize( feMesh );
	nNodalDofs = self->dofLayout->dofCounts;
	links = self->linkedDofInfo;

	/* Allocate for destination array. */
	dstArray = Memory_Alloc_2DComplex( int, nDomainNodes, nNodalDofs, 
					   "FeEquationNumber::destinationArray" );

	/* If needed, allocate for linked equation numbers. */
	if( links ) {
		unsigned	s_i;

		links->eqNumsOfLinkedDofs = ReallocArray( links->eqNumsOfLinkedDofs, int, links->linkedDofSetsCount );
		for( s_i = 0; s_i < links->linkedDofSetsCount; s_i++ )
			links->eqNumsOfLinkedDofs[s_i] = -1;
	}

	/* Allocate for the location matrix. */
	nLocMatDofs = NULL;
	locMat = AllocArray( int**, nDomainEls );
	for( e_i = 0; e_i < nDomainEls; e_i++ ) {
		FeMesh_GetElementNodes( feMesh, e_i, &nElNodes, &elNodes );
		nLocMatDofs = ReallocArray( nLocMatDofs, int, nElNodes );
		for( n_i = 0; n_i < nElNodes; n_i++ )
			nLocMatDofs[n_i] = nNodalDofs[elNodes[n_i]];
		locMat[e_i] = AllocComplex2D( int, nElNodes, nLocMatDofs );
	}
	FreeArray( nLocMatDofs );

	/* Build initial destination array and store max dofs. */
	curEqNum = 0;
	maxDofs = 0;
	for( n_i = 0; n_i < nLocalNodes; n_i++ ) {
		if( nNodalDofs[n_i] > maxDofs )
			maxDofs = nNodalDofs[n_i];

		for( dof_i = 0; dof_i < nNodalDofs[n_i]; dof_i++ ) {
			varInd = self->dofLayout->varIndices[n_i][dof_i];
			if( !self->bcs || !VariableCondition_IsCondition( self->bcs, n_i, varInd ) || 
			    !self->removeBCs )
			{
				if( links && links->linkedDofTbl[n_i][dof_i] != -1 ) {
					if( links->eqNumsOfLinkedDofs[links->linkedDofTbl[n_i][dof_i]] == -1 )
						links->eqNumsOfLinkedDofs[links->linkedDofTbl[n_i][dof_i]] = curEqNum++;
					dstArray[n_i][dof_i] = links->eqNumsOfLinkedDofs[links->linkedDofTbl[n_i][dof_i]];
				}
				else
					dstArray[n_i][dof_i] = curEqNum++;
			}
			else
				dstArray[n_i][dof_i] = -1;
		}
	}

	/* Order the equation numbers based on processor rank; cascade counts forward. */
	base = 0;
	subTotal = curEqNum;
	if( rank > 0 ) {
		MPI_Recv( &base, 1, MPI_UNSIGNED, rank - 1, 6669, comm, &status );
		subTotal = base + curEqNum;
	}
	if( rank < nProcs - 1 )
		MPI_Send( &subTotal, 1, MPI_UNSIGNED, rank + 1, 6669, comm );

	/* Reduce to find lowest linked DOFs. */
	if( links ) {
		unsigned	lowest, highest;
		unsigned	s_i;

		for( s_i = 0; s_i < links->linkedDofSetsCount; s_i++ ) {
			if( links->eqNumsOfLinkedDofs[s_i] != -1 )
				links->eqNumsOfLinkedDofs[s_i] += base;
			MPI_Allreduce( links->eqNumsOfLinkedDofs + s_i, &lowest, 1, MPI_UNSIGNED, MPI_MIN, comm );
			MPI_Allreduce( links->eqNumsOfLinkedDofs + s_i, &highest, 1, MPI_UNSIGNED, MPI_MIN, comm );
			assert( (lowest == (unsigned)-1) ? lowest == highest : 1 );
			links->eqNumsOfLinkedDofs[s_i] = lowest;
		}
	}

	/* Modify existing destination array and dump to a tuple array. */
	tuples = AllocArray( unsigned, nDomainNodes * maxDofs );
	for( n_i = 0; n_i < nLocalNodes; n_i++ ) {
		for( dof_i = 0; dof_i < nNodalDofs[n_i]; dof_i++ ) {
			varInd = self->dofLayout->varIndices[n_i][dof_i];
			if( !self->bcs || !VariableCondition_IsCondition( self->bcs, n_i, varInd ) || 
			    !self->removeBCs )
			{
				if( links && links->linkedDofTbl[n_i][dof_i] != -1 )
					dstArray[n_i][dof_i] = links->eqNumsOfLinkedDofs[links->linkedDofTbl[n_i][dof_i]];
				else
					dstArray[n_i][dof_i] += base;
			}
			tuples[n_i * maxDofs + dof_i] = dstArray[n_i][dof_i];
		}
	}

	/* Update all other procs. */
	sync = Mesh_GetSync( feMesh, MT_VERTEX );
	array = Decomp_Sync_Array_New();
	Decomp_Sync_Array_SetSync( array, sync );
	Decomp_Sync_Array_SetMemory( array, tuples, tuples + nLocalNodes * maxDofs, 
				     maxDofs * sizeof(unsigned), maxDofs * sizeof(unsigned), 
				     maxDofs * sizeof(unsigned) );
	Decomp_Sync_Array_Sync( array );
	FreeObject( array );

	/* Update destination array's domain indices. */
	for( n_i = nLocalNodes; n_i < nDomainNodes; n_i++ ) {
		for( dof_i = 0; dof_i < nNodalDofs[n_i]; dof_i++ ) {
			varInd = self->dofLayout->varIndices[n_i][dof_i];
			if( !self->bcs || !VariableCondition_IsCondition( self->bcs, n_i, varInd ) || 
			    !self->removeBCs )
			{
				dstArray[n_i][dof_i] = tuples[n_i * maxDofs + dof_i];
			}
			else
				dstArray[n_i][dof_i] = -1;
		}
	}

	/* Destroy tuple array. */
	FreeArray( tuples );

	/* Build location matrix. */
	for( e_i = 0; e_i < nDomainEls; e_i++ ) {
		FeMesh_GetElementNodes( feMesh, e_i, &nElNodes, &elNodes );
		for( n_i = 0; n_i < nElNodes; n_i++ ) {
			for( dof_i = 0; dof_i < nNodalDofs[elNodes[n_i]]; dof_i++ )
				locMat[e_i][n_i][dof_i] = dstArray[elNodes[n_i]][dof_i];
		}
	}

	/* Store stuff on class. */
	self->destinationArray = dstArray;
	self->locationMatrix = locMat;
	self->locationMatrixBuilt = True;
	self->remappingActivated = False;
	self->localEqNumsOwnedCount = curEqNum;
	self->firstOwnedEqNum = base;
	self->lastOwnedEqNum = subTotal - 1;
	self->_lowestLocalEqNum = self->firstOwnedEqNum;

	/* Bcast global sum from highest rank. */
	if( rank == nProcs - 1 )
		self->globalSumUnconstrainedDofs = self->lastOwnedEqNum + 1;
	MPI_Bcast( &self->globalSumUnconstrainedDofs, 1, MPI_UNSIGNED, nProcs - 1, comm );

	/* Construct lowest global equation number list. */
	self->_lowestGlobalEqNums = AllocArray( int, nProcs );
	MPI_Allgather( &self->firstOwnedEqNum, 1, MPI_UNSIGNED, self->_lowestGlobalEqNums, 1, MPI_UNSIGNED, comm );

	endTime = MPI_Wtime();

	Journal_Printf( stream, "Assigned %d global equation numbers.\n", self->globalSumUnconstrainedDofs );
	Journal_Printf( stream, "Assigned %d local equation numbers.\n", self->lastOwnedEqNum - self->firstOwnedEqNum + 1 );
	Journal_Printf( stream, "Local equation numbers in the range %d to %d.\n", 
			self->firstOwnedEqNum, self->lastOwnedEqNum + 1 );
	Stream_UnIndent( stream );
	Journal_Printf( stream, "... Completed in %g seconds.\n", endTime - startTime );
	Stream_UnIndent( stream );
}

#if 0
void FeEquationNumber_Invert( void* feEqNum, int equation, unsigned* node, unsigned* dof ) {
	FeEquationNumber*	self = (FeEquationNumber*)feEqNum;

	assert( self && Stg_CheckType( self, FeEquationNumber ) );
	assert( equation - self->firstOwnedEqNum < self->localEqNumsOwnedCount );
	assert( node );
	assert( dof );

	eq = equation - self->firstOwnedEqNum;
	*node = self->eqToNode[eq];
	*dof = self->eqToDof[eq];
}

Bool FeEquationNumber_IsKnown( void* feEqNum, int equation ) {
	FeEquationNumber*	self = (FeEquationNumber*)feEqNum;
	unsigned		node, dof;
	unsigned		varInd;

	assert( self && Stg_CheckType( self, FeEquationNumber ) );

	if( !self->bcs ) return False;
	FeEquationNumber_Invert( self, equation, &node, &dof );

	assert( self->dofLayout );
	assert( self->dofLayout->varIndices );
	assert( self->dofLayout->varIndices[node] );
	varInd = self->dofLayout->varIndices[node][dof];
	return VariableCondition_IsCondition( self->bcs, n_i, varInd );
}
#endif
