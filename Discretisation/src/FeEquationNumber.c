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
** $Id: FeEquationNumber.c 1191 2008-07-25 03:06:19Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
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
#include <petsc.h>
#include <petscvec.h>

int stgCmpInt( const void *l, const void *r ) {
   return *(int*)l - *(int*)r;
}

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

int GenerateEquationNumbering(
		int NX, int NY, int NZ,
		int nlocal, int g_node_id[],
		int dof, int nglobal,
		PetscTruth periodic_x, PetscTruth periodic_y, PetscTruth periodic_z,
		int npx, int npy, int npz,
		int periodic_x_gnode_id[], int periodic_y_gnode_id[], int periodic_z_gnode_id[],
		int eqnums[], int *neqnums );

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
                                 _FeEquationNumber_AssignFromXML, (Stg_Component_BuildFunction*)_FeEquationNumber_Build, 
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
                                 _FeEquationNumber_AssignFromXML, (Stg_Component_BuildFunction*)_FeEquationNumber_Build, 
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
   self->_construct = _FeEquationNumber_AssignFromXML;
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
   self->bcEqNums = STree_New();
   STree_SetIntCallbacks( self->bcEqNums );
   STree_SetItemSize( self->bcEqNums, sizeof(int) );
   self->ownedMap = STreeMap_New();
   STreeMap_SetItemSize( self->ownedMap, sizeof(int), sizeof(int) );
   STree_SetIntCallbacks( self->ownedMap );

   Stream_SetPrintingRank( self->debug, 0 );
 

}

void _FeEquationNumber_AssignFromXML( void* feEquationNumber, Stg_ComponentFactory *cf, void* data ){
	FeEquationNumber* self = (FeEquationNumber*) feEquationNumber;
	
	self->context = Stg_ComponentFactory_ConstructByKey( cf, self->name, "Context", DomainContext, False, data );
	if( !self->context ) 
		self->context = Stg_ComponentFactory_ConstructByName( cf, "context", DomainContext, True, data );
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
      NewClass_Delete( self->bcEqNums );
   }
   if( self->ownedMap ) {
      NewClass_Delete( self->ownedMap );
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

   /* If we have new mesh topology information, do this differently. */
/*
   if( self->feMesh->topo->domains && self->feMesh->topo->domains[MT_VERTEX] ) {
*/
   if( Mesh_HasExtension( self->feMesh, "vertexGrid" ) )
      FeEquationNumber_BuildWithDave( self );
   else
      FeEquationNumber_BuildWithTopology( self );

   /* If not removing BCs, construct a table of which equation numbers are actually BCs. */
   if( !self->removeBCs ) {
      FeMesh* mesh = self->feMesh;
      DofLayout* dofLayout = self->dofLayout;
      VariableCondition* bcs = self->bcs;
      int nDofs, varInd;
      int ii, jj;

      for( ii = 0; ii < FeMesh_GetNodeLocalSize( mesh ); ii++ ) {
         nDofs = dofLayout->dofCounts[ii];
         for( jj = 0; jj < nDofs; jj++ ) {
            varInd = dofLayout->varIndices[ii][jj];
            if( bcs && VariableCondition_IsCondition( bcs, ii, varInd ) ) {
               if( !STree_Has( self->bcEqNums, self->destinationArray[ii] + jj ) )
                  STree_Insert( self->bcEqNums, self->destinationArray[ii] + jj );
            }
         }
      }
   }

/*
   }
   else {
      _FeEquationNumber_BuildRemappedNodeInfoTable( self );
      _FeEquationNumber_BuildDestinationArray( self );
      _FeEquationNumber_CalculateGlobalUnconstrainedDofTotal( self );
      _FeEquationNumber_CalculateEqNumsDecomposition( self );
   }
*/

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
   IArray*			inc;
   unsigned		e_i, n_i, dof_i;

   assert( self );

   /* Don't build if already done. */
   if( self->locationMatrixBuilt ) {
      Journal_DPrintf( self->debugLM, "In %s: LM already built, so just returning.\n",  __func__ );
      Stream_UnIndentBranch( StgFEM_Debug );
      return;
   }

   inc = IArray_New();

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
      FeMesh_GetElementNodes( feMesh, e_i, inc );
      nElNodes = IArray_GetSize( inc );
      elNodes = IArray_GetPtr( inc );
      locMat[e_i] = AllocArray( int*, nElNodes );
      for( n_i = 0; n_i < nElNodes; n_i++ )
         locMat[e_i][n_i] = AllocArray( int, nNodalDofs[elNodes[n_i]] );
   }

   /* Build location matrix. */
   for( e_i = 0; e_i < nDomainEls; e_i++ ) {
      FeMesh_GetElementNodes( feMesh, e_i, inc );
      nElNodes = IArray_GetSize( inc );
      elNodes = IArray_GetPtr( inc );
      for( n_i = 0; n_i < nElNodes; n_i++ ) {
         for( dof_i = 0; dof_i < nNodalDofs[elNodes[n_i]]; dof_i++ )
            locMat[e_i][n_i][dof_i] = dstArray[elNodes[n_i]][dof_i];
      }
   }

   NewClass_Delete( inc );

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
   IArray* inc;

   inc = IArray_New();
   FeMesh_GetElementNodes( feMesh, lElement_I, inc );
   numNodesThisElement = IArray_GetSize( inc );
   elInc = IArray_GetPtr( inc );

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

   NewClass_Delete( inc );

   return localLocationMatrix;
}


void FeEquationNumber_PrintDestinationArray( void* feFeEquationNumber, Stream* stream ) {
   FeEquationNumber* self = (FeEquationNumber*) feFeEquationNumber;
   FeMesh*		feMesh = self->feMesh;
   MPI_Comm comm = Comm_GetMPIComm( Mesh_GetCommTopology( feMesh, MT_VERTEX ) );
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

void FeEquationNumber_PrintDestinationArrayBox( void* feFeEquationNumber, Stream* stream ) {
	FeEquationNumber*   self = (FeEquationNumber*) feFeEquationNumber;
	FeMesh*             feMesh = self->feMesh;
	MPI_Comm            comm = Comm_GetMPIComm( Mesh_GetCommTopology( feMesh, MT_VERTEX ) );
	unsigned            rank;
	Grid*               vGrid;
	int                 ijk[3];
	Element_LocalIndex  lNode_I;
	Node_GlobalIndex    gNode_I;
	Node_Index          sizeI, sizeJ, sizeK;
	Dof_Index           nDofs;
	Dof_Index           dof_I;

	MPI_Comm_rank( comm, (int*)&rank );
	Journal_Printf( stream, "%d: *** Printing destination array ***\n", rank );
	
	vGrid = *Mesh_GetExtension( feMesh, Grid**, "vertexGrid" );
	nDofs = self->dofLayout->dofCounts[0];
	sizeI = vGrid->sizes[ I_AXIS ];
	sizeJ = vGrid->sizes[ J_AXIS ];
	sizeK = vGrid->sizes[ K_AXIS ] ? vGrid->sizes[ K_AXIS ] : 1;

	for ( ijk[2] = 0 ; ijk[2] < sizeK ; ijk[2]++ ) {
		if ( sizeK != 1 )
			Journal_Printf( stream, "\nk = %d\n", ijk[2] );
		for ( ijk[1] = sizeJ - 1 ; ijk[1] >= 0 ; ijk[1]-- ) {
			Journal_Printf( stream, "%2d - ", ijk[1] );
			for ( ijk[0] = 0 ; ijk[0] < sizeI ; ijk[0]++ ) {
				gNode_I = Grid_Project( vGrid, ijk );
				Journal_Printf( stream, "{ " );
				if ( Mesh_GlobalToDomain( feMesh, MT_VERTEX, gNode_I, &lNode_I ) ) {
					for ( dof_I = 0 ; dof_I < nDofs ; dof_I++ )
						Journal_Printf( stream, "%3d ", self->destinationArray[lNode_I][dof_I] );
				}
				else {
					for ( dof_I = 0 ; dof_I < nDofs ; dof_I++ )
						Journal_Printf( stream, " XX " );
				}
				Journal_Printf( stream, "}" );
			}
			Journal_Printf( stream, "\n" );
		}
		/* Bottom row */
		Journal_Printf( stream, "    " );
		for ( ijk[0] = 0 ; ijk[0] < sizeI ; ijk[0]++ ) {
			Journal_Printf( stream, "    %3d    ", ijk[0] );
			if ( nDofs == 3 )
				Journal_Printf( stream, "    " );
		}
		Journal_Printf( stream, "\n" );
	}	
}

void FeEquationNumber_PrintLocationMatrix( void* feFeEquationNumber, Stream* stream ) {
   FeEquationNumber* self = (FeEquationNumber*) feFeEquationNumber;
   FeMesh*		feMesh = self->feMesh;
   MPI_Comm comm = Comm_GetMPIComm( Mesh_GetCommTopology( feMesh, MT_VERTEX ) );
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
         IArray*		inc;

         inc = IArray_New();
         FeMesh_GetElementNodes( self->feMesh, lEl_I, inc );
         numNodesAtElement = IArray_GetSize( inc );
         incNodes = IArray_GetPtr( inc );

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

         NewClass_Delete( inc );
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
   Comm*	comm = Mesh_GetCommTopology( self->feMesh, MT_VERTEX );
   MPI_Comm	mpiComm = Comm_GetMPIComm( comm );
   unsigned	nProcs;
   unsigned	p_i;

   MPI_Comm_size( mpiComm, (int*)&nProcs );
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
   Stream*		stream;
   double		startTime, endTime, time, tmin, tmax;
   FeMesh*		feMesh;
   Sync*		sync;
   Comm*		comm;
   MPI_Comm		mpiComm;
   unsigned		rank, nProcs;
   unsigned		nDims;
   unsigned		nDomainNodes, nDomainEls;
   unsigned		nLocalNodes;
   unsigned*		nNodalDofs;
   unsigned		nElNodes, *elNodes;
   int**		dstArray;
   int			*nLocMatDofs, ***locMat;
   unsigned		varInd;
   unsigned		curEqNum;
   unsigned		base;
   unsigned		subTotal;
   MPI_Status		status;
   unsigned		maxDofs;
   unsigned*		tuples;
   LinkedDofInfo*	links;
   unsigned		highest;
   IArray*		inc;
   unsigned             e_i, n_i, dof_i, s_i;
   int			ii;

   assert( self );

   inc = IArray_New();

   stream = Journal_Register( Info_Type, self->type );
   Stream_SetPrintingRank( stream, 0 );

   Journal_RPrintf( stream, "FeEquationNumber: '%s'\n", self->name );
   Stream_Indent( stream );
   Journal_RPrintf( stream, "Generating equation numbers...\n" );
   Stream_Indent( stream );
   if( self->removeBCs )
      Journal_RPrintf( stream, "BCs set to be removed.\n" );
   else
      Journal_RPrintf( stream, "BCs will not be removed.\n" );

   startTime = MPI_Wtime();

   /* Shortcuts. */
   feMesh = self->feMesh;
   comm = Mesh_GetCommTopology( feMesh, MT_VERTEX );
   mpiComm = Comm_GetMPIComm( comm );
   MPI_Comm_size( mpiComm, (int*)&nProcs );
   MPI_Comm_rank( mpiComm, (int*)&rank );
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
      FeMesh_GetElementNodes( feMesh, e_i, inc );
      nElNodes = IArray_GetSize( inc );
      elNodes = IArray_GetPtr( inc );
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
               if( rank > 0 ) {
                  dstArray[n_i][dof_i] = -2;
                  continue;
               }
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
   if( rank > 0 )
      MPI_Recv( &base, 1, MPI_UNSIGNED, rank - 1, 6669, mpiComm, &status );
   subTotal = base + curEqNum;
   if( rank < nProcs - 1 )
      MPI_Send( &subTotal, 1, MPI_UNSIGNED, rank + 1, 6669, mpiComm );

   if( links ) {
      /* Reduce to find lowest linked DOFs. */
      for( s_i = 0; s_i < links->linkedDofSetsCount; s_i++ ) {
         if( links->eqNumsOfLinkedDofs[s_i] != -1 )
            links->eqNumsOfLinkedDofs[s_i] += base;
/*
  MPI_Allreduce( links->eqNumsOfLinkedDofs + s_i, &lowest, 1, MPI_UNSIGNED, MPI_MAX, mpiComm );
*/
         MPI_Allreduce( links->eqNumsOfLinkedDofs + s_i, &highest, 1, MPI_INT, MPI_MAX, mpiComm );
/*
  assert( (lowest == (unsigned)-1) ? lowest == highest : 1 );
*/
         links->eqNumsOfLinkedDofs[s_i] = highest;
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
            if( links && links->linkedDofTbl[n_i][dof_i] != -1 ) {
               highest = links->eqNumsOfLinkedDofs[links->linkedDofTbl[n_i][dof_i]];
               dstArray[n_i][dof_i] = highest;
            }
            else
               dstArray[n_i][dof_i] += base;
         }
         tuples[n_i * maxDofs + dof_i] = dstArray[n_i][dof_i];
      }
   }

   /* Update all other procs. */
   sync = Mesh_GetSync( feMesh, MT_VERTEX );
   Sync_SyncArray( sync, tuples, maxDofs * sizeof(unsigned), 
                   tuples + nLocalNodes * maxDofs, maxDofs * sizeof(unsigned), 
                   maxDofs * sizeof(unsigned) );

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
      FeMesh_GetElementNodes( feMesh, e_i, inc );
      nElNodes = IArray_GetSize( inc );
      elNodes = IArray_GetPtr( inc );
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

   /* Setup owned mapping. */
   STree_Clear( self->ownedMap );
   for( ii = self->firstOwnedEqNum; ii <= self->lastOwnedEqNum; ii++ ) {
      int val = ii - self->firstOwnedEqNum;
      STreeMap_Insert( self->ownedMap, &ii, &val );
   }

   /* Bcast global sum from highest rank. */
   if( rank == nProcs - 1 )
      self->globalSumUnconstrainedDofs = self->lastOwnedEqNum + 1;
   MPI_Bcast( &self->globalSumUnconstrainedDofs, 1, MPI_UNSIGNED, nProcs - 1, mpiComm );

   /* Construct lowest global equation number list. */
   self->_lowestGlobalEqNums = AllocArray( int, nProcs );
   MPI_Allgather( &self->firstOwnedEqNum, 1, MPI_UNSIGNED, self->_lowestGlobalEqNums, 1, MPI_UNSIGNED, mpiComm );

   endTime = MPI_Wtime();

   Journal_RPrintf( stream, "Assigned %d global equation numbers.\n", self->globalSumUnconstrainedDofs );
   Journal_Printf( stream, "[%u] Assigned %d local equation numbers, within range %d to %d.\n", 
                   rank, self->lastOwnedEqNum - self->firstOwnedEqNum + 1, self->firstOwnedEqNum, self->lastOwnedEqNum + 1 );
   Stream_UnIndent( stream );

   time = endTime - startTime;
   MPI_Reduce( &time, &tmin, 1, MPI_DOUBLE, MPI_MIN, 0, mpiComm );
   MPI_Reduce( &time, &tmax, 1, MPI_DOUBLE, MPI_MAX, 0, mpiComm );
   Journal_RPrintf( stream, "... Completed in %g [min] / %g [max] seconds.\n", tmin, tmax );
   Stream_UnIndent( stream );

   NewClass_Delete( inc );
}

#if 0
void FeEquationNumber_BuildVariableIndices( FeEquationNumber* self, int* nInds, int** inds, int* maxDofs ) {
   int maxDofs;
   ISet indSetObj, *indSet = &indSetObj;

   *maxDofs = 0;
   for( n_i = 0; n_i < nLocalNodes; n_i++ ) {
      if( nNodalDofs[n_i] > *maxDofs )
         *maxDofs = nNodalDofs[n_i];
   }

   ISet_Construct( indSet );
   ISet_SetMaxSize( indSet, *maxDofs );
   for( n_i = 0; n_i < nLocalNodes; n_i++ ) {
      for( dof_i = 0; dof_i < nNodalDofs[n_i]; dof_i++ ) {
         varInd = self->dofLayout->varIndices[n_i][dof_i];
         ISet_TryInsert( indSet, varInd );
      }
   }

   *nInds = ISet_GetSize( indSet );
   *inds = MemArray( int, *nInds, FeEquationNumber_Type );
   ISet_GetArray( indSet, NULL, inds );
   ISet_Destruct( indSet );
}
#endif


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



void FeEquationNumber_BuildWithDave( FeEquationNumber* self ) {
   int nLocals, *locals;
   Grid *vGrid;
   int varInd;
   int nEqNums, **dstArray;
   IArray *inc;
   int nDofs;
   int *periodic;
   int ***locMat;
   int nDims;
   int *elNodes;
   Comm *comm;
   MPI_Comm mpiComm;
   int nRanks, rank;
   Sync *sync;
   Bool isCond;
   int nPeriodicInds[3];
   int *periodicInds[3];
   int inds[3];
   Bool usePeriodic;
   int *tmpArray, nLocalEqNums;
   int lastOwnedEqNum, ind;
   STree *doneSet;
   int ii, jj, kk;

   comm = Mesh_GetCommTopology( self->feMesh, 0 );
   mpiComm = Comm_GetMPIComm( comm );
   MPI_Comm_size( mpiComm, &nRanks );
   MPI_Comm_rank( mpiComm, &rank );

   /* Setup an array containing global indices of all locally owned nodes. */
   nLocals = Mesh_GetLocalSize( self->feMesh, 0 );
   locals = AllocArray( int, nLocals );
   for( ii = 0; ii < nLocals; ii++ )
      locals[ii] = Mesh_DomainToGlobal( self->feMesh, 0, ii );

   /* Allocate for destination array. */
   nDofs = self->dofLayout->dofCounts[0];
   dstArray = AllocArray2D( int, Mesh_GetDomainSize( self->feMesh, 0 ), nDofs );

   /* Get the vertex grid extension and any periodicity. */
   nDims = Mesh_GetDimSize( self->feMesh );
   vGrid = *Mesh_GetExtension( self->feMesh, Grid**, "vertexGrid" );
   periodic = Mesh_GetExtension( self->feMesh, int*, "periodic" );

   /* Fill destination array with initial values, setting dirichlet BCs as we go. */
   for( ii = 0; ii < nLocals; ii++ ) {
      Grid_Lift( vGrid, locals[ii], inds );
      usePeriodic = False;
      for( jj = 0; jj < nDims; jj++ ) {
	 if( periodic[jj] && (inds[jj] == 0 || inds[jj] == vGrid->sizes[jj] - 1) ) {
            usePeriodic = True;
	    break;
	 }
      }
      for( jj = 0; jj < nDofs; jj++ ) {
         varInd = self->dofLayout->varIndices[ii][jj];
	 if( self->bcs )
	    isCond = VariableCondition_IsCondition( self->bcs, ii, varInd );
	 else
	    isCond = False;
         if( isCond && self->removeBCs )
            dstArray[ii][jj] = -1;
         else
            dstArray[ii][jj] = 0;
      }
   }

   /* Generate opposing indices for periodicity. */
   for( ii = 0; ii < nDims; ii++ ) {
      nPeriodicInds[ii] = 0;
      periodicInds[ii] = NULL;
      if( periodic[ii] ) {
	 periodicInds[ii] = AllocArray( int, nLocals );
	 for( jj = 0; jj < nLocals; jj++ ) {
	    Grid_Lift( vGrid, locals[jj], inds );
	    if( inds[ii] != vGrid->sizes[ii] - 1 ) continue;
/*
            for( kk = 0; kk < nDofs; kk++ )
               if( dstArray[jj][kk] == -1 ) break;
            if( kk < nDofs ) continue;
*/
            periodicInds[ii][nPeriodicInds[ii]++] = locals[jj];
	 }
      }
   }

   /* Call Dave's equation number generation routine. */
   if( nDims == 2 ) {
      GenerateEquationNumbering( vGrid->sizes[0], vGrid->sizes[1], 1,
				 nLocals, locals,
				 nDofs, Mesh_GetGlobalSize( self->feMesh, 0 ),
				 periodic[0], periodic[1], False,
				 nPeriodicInds[0], nPeriodicInds[1], 0,
				 periodicInds[0], periodicInds[1], NULL,
				 dstArray[0], &nEqNums );
   }
   else {
      GenerateEquationNumbering( vGrid->sizes[0], vGrid->sizes[1], vGrid->sizes[2],
				 nLocals, locals,
				 nDofs, Mesh_GetGlobalSize( self->feMesh, 0 ),
				 periodic[0], periodic[1], periodic[2],
				 nPeriodicInds[0], nPeriodicInds[1], nPeriodicInds[2],
				 periodicInds[0], periodicInds[1], periodicInds[2],
				 dstArray[0], &nEqNums );
   }

   /* Free periodic arrays. */
   for( ii = 0; ii < nDims; ii++ ) {
      if( periodicInds[ii] )
	 FreeArray( periodicInds[ii] );
   }

   /* Setup owned mapping part 1. */
   STree_Clear( self->ownedMap );
   for( ii = 0; ii < nLocals; ii++ ) {
      Grid_Lift( vGrid, locals[ii], inds );
      for( jj = 0; jj < nDims; jj++ ) {
	 if( periodic[jj] && inds[jj] == vGrid->sizes[jj] - 1 ) {
            inds[jj] = 0;
            ind = Grid_Project( vGrid, inds );
            if( !FeMesh_NodeGlobalToDomain( self->feMesh, ind, &ind ) )
               break;
         }
      }
      if( jj < nDims ) continue;
      for( jj = 0; jj < nDofs; jj++ ) {
	 if( dstArray[ii][jj] == -1 || STreeMap_HasKey( self->ownedMap, dstArray[ii] + jj ) )
            continue;
	 STreeMap_Insert( self->ownedMap, dstArray[ii] + jj, &ii );
      }
   }

   /* Setup owned mapping. */
   tmpArray = AllocArray( int, nLocals * nDofs );
   memcpy( tmpArray, dstArray[0], nLocals * nDofs * sizeof(int) );
   qsort( tmpArray, nLocals * nDofs, sizeof(int), stgCmpInt );
   doneSet = STree_New();
   STree_SetItemSize( doneSet, sizeof(int) );
   STree_SetIntCallbacks( doneSet );
   for( nLocalEqNums = 0, ii = 0; ii < nLocals * nDofs; ii++ ) {
      if( tmpArray[ii] != -1 &&
          STreeMap_HasKey( self->ownedMap, tmpArray + ii ) &&
          !STree_Has( doneSet, tmpArray + ii ) )
      {
	 if( !nLocalEqNums )
	    self->_lowestLocalEqNum = tmpArray[ii];
	 *(int*)STreeMap_Map( self->ownedMap, tmpArray + ii ) = nLocalEqNums;
         STree_Insert( doneSet, tmpArray + ii );
	 nLocalEqNums++;
      }
   }
   lastOwnedEqNum = -1; /* Don't need this anymore. */
   FreeArray( tmpArray );
   NewClass_Delete( doneSet );

   /* Transfer remote equation numbers. */
   sync = Mesh_GetSync( self->feMesh, 0 );
   Sync_SyncArray( sync, dstArray[0], nDofs * sizeof(int),
                   dstArray[0] + nLocals * nDofs, nDofs * sizeof(int),
                   nDofs * sizeof(int) );

   /* Allocate for location matrix. */
   locMat = AllocArray( int**, Mesh_GetDomainSize( self->feMesh, nDims ) );
   for( ii = 0; ii < Mesh_GetDomainSize( self->feMesh, nDims ); ii++ )
      locMat[ii] = AllocArray2D( int, FeMesh_GetElementNodeSize( self->feMesh, 0 ), nDofs );

   /* Fill in location matrix. */
   inc = IArray_New();
   for( ii = 0; ii < Mesh_GetDomainSize( self->feMesh, nDims ); ii++ ) {
      FeMesh_GetElementNodes( self->feMesh, ii, inc );
      elNodes = IArray_GetPtr( inc );
      for( jj = 0; jj < FeMesh_GetElementNodeSize( self->feMesh, 0 ); jj++ ) {
         for( kk = 0; kk < nDofs; kk++ )
            locMat[ii][jj][kk] = dstArray[elNodes[jj]][kk];
      }
   }
   NewClass_Delete( inc );

   /* Fill in our other weird values. */
   self->destinationArray = dstArray;
   self->locationMatrix = locMat;
   self->locationMatrixBuilt = True;
   self->remappingActivated = False;
   self->localEqNumsOwnedCount = nLocalEqNums;

   /* Bcast global sum from highest rank. */
   self->globalSumUnconstrainedDofs = nEqNums;

   /* Construct lowest global equation number list. */
   self->_lowestGlobalEqNums = AllocArray( int, nRanks );
   MPI_Allgather( &self->_lowestLocalEqNum, 1, MPI_UNSIGNED, 
                  self->_lowestGlobalEqNums, 1, MPI_UNSIGNED,
                  mpiComm );

   FreeArray( locals );

   /* 
   printf( "%d: localEqNumsOwned = %d\n", rank, self->localEqNumsOwnedCount );
   printf( "%d: globalSumUnconstrainedDofs = %d\n", rank, self->globalSumUnconstrainedDofs );
   */
}




/*

Input:
  nlocal - number of locally ownded nodes
  g_node_id - global indices of nodes owned locally. Size nlocal
  dof - degrees of freedom per node
  nglobal - number of global nodes
  npx - number of consider to be periodic in x (local to this proc)
  npy - number of consider to be periodic in y (local to this proc)
  periodic_x_gnode_id - global indices of nodes (on this proc) which are on right hand side boundary
  periodic_y_gnode_id - global indices of nodes (on this proc) which are on top boundary
  eqnums - contains any dirichlet boundary conditions. Size nlocal*dof

Output;
  eqnums - contains full list of eqnums

Assumptions: 
- Ordering eqnums[] = { (node_0,[dof_0,dof_1,..,dof_x]), (node_1,[dof_0,dof_1,..,dof_x]), ... }
- Any dirichlet set along a boundary deemed to be periodic will be clobbered.
- Processors may have duplicate nodes in the g_node_id[] list.
- A number in the corner is considered part of both boundaries (horiz and vert)
- If npx is not 0, then periodicity is assumed in x
- If npy is not 0, then periodicity is assumed in y
- Dofs constrained to be dirichlet must be marked with a negative number.
- Dofs are NOT split across processors.
- We can define a logical i,j,k ordering to uniquely identify nodes.

*/

PetscErrorCode _VecScatterBeginEnd( VecScatter vscat, Vec FROM, Vec TO, InsertMode addv,ScatterMode mode )
{
#if( (PETSC_VERSION_MAJOR==2) && (PETSC_VERSION_MINOR==3) && (PETSC_VERSION_SUBMINOR==2) )
	// 2.3.2 ordering of args
	VecScatterBegin( FROM, TO, addv, mode, vscat );
	VecScatterEnd( FROM, TO, addv, mode, vscat );
#else
	// 2.3.3 or 3.0.0
	VecScatterBegin( vscat, FROM, TO, addv, mode );
	VecScatterEnd( vscat, FROM, TO, addv, mode );
#endif
	
	PetscFunctionReturn(0);
}

int GenerateEquationNumbering(
		int NX, int NY, int NZ,
		int nlocal, int g_node_id[],
		int dof, int nglobal,
		PetscTruth periodic_x, PetscTruth periodic_y, PetscTruth periodic_z,
		int npx, int npy, int npz,
		int periodic_x_gnode_id[], int periodic_y_gnode_id[], int periodic_z_gnode_id[],
		int eqnums[], int *neqnums )
{
	PetscErrorCode ierr;
	PetscInt periodic_mask;
	Vec global_eqnum, g_ownership;
	PetscInt i;
	PetscMPIInt rank;        /* processor rank */
	PetscMPIInt size;        /* size of communicator */
	Vec local_ownership, local_eqnum;
	PetscInt *_g_node_id;
	IS is_gnode, is_eqnum;
	VecScatter vscat_ownership, vscat_eqnum;
	PetscScalar *_local_ownership, *_local_eqnum;
	PetscInt local_eqnum_count,global_eqnum_count;
	PetscScalar val[10];
	PetscInt d,idx[10];
	PetscInt *to_fetch,cnt,number_to_fetch;
	PetscInt eq_cnt;
	Vec offset_list;
	VecScatter vscat_offset;
	Vec seq_offset_list;
	PetscInt offset, inc;
	
	PetscInt spanx,spany,spanz,total;
	PetscInt loc;
	PetscReal max;
	PetscInt n_inserts;
	
	
	ierr = MPI_Comm_rank(PETSC_COMM_WORLD,&rank);CHKERRQ(ierr);
	ierr = MPI_Comm_size(PETSC_COMM_WORLD,&size);CHKERRQ(ierr);
	
	if( dof>=10 ) {
		SETERRQ(PETSC_ERR_SUP, "Max allowable degrees of freedom per node = 10. Change static size" );
	}
	
	
	/* 
	Claim locally owned nodes. Duplicate nodes on the interior will be resolved by the processor
	which inserts last.
	*/
	VecCreate( PETSC_COMM_WORLD, &g_ownership );
	VecSetSizes( g_ownership, PETSC_DECIDE, nglobal );
	VecSetFromOptions( g_ownership );
	
	for( i=0; i<nlocal; i++ ) {
		VecSetValue( g_ownership, g_node_id[i], rank, INSERT_VALUES );
	}
	VecAssemblyBegin(g_ownership);
	VecAssemblyEnd(g_ownership);
	
	
	/* Mask out the periodic boundaries. */
	periodic_mask = -6699.0;
	if (periodic_x_gnode_id!=NULL) {
		for( i=0; i<npx; i++ ) {
			VecSetValue( g_ownership, periodic_x_gnode_id[i], periodic_mask, INSERT_VALUES );
		}
	}
	if (periodic_y_gnode_id!=NULL) {
		for( i=0; i<npy; i++ ) {
			VecSetValue( g_ownership, periodic_y_gnode_id[i], periodic_mask, INSERT_VALUES );
		}
	}
	if (periodic_z_gnode_id!=NULL) {
		for( i=0; i<npz; i++ ) {
			VecSetValue( g_ownership, periodic_z_gnode_id[i], periodic_mask, INSERT_VALUES );
		}
	}
	VecAssemblyBegin(g_ownership);
	VecAssemblyEnd(g_ownership);
	
	/*
	PetscPrintf(PETSC_COMM_WORLD, "g_ownership \n");
	VecView( g_ownership, PETSC_VIEWER_STDOUT_WORLD );
	*/
	
	/* Get all locally owned nodes */
	VecCreate( PETSC_COMM_SELF, &local_ownership );
	VecSetSizes( local_ownership, PETSC_DECIDE, nlocal );
	VecSetFromOptions( local_ownership );
	
	PetscMalloc( sizeof(PetscInt)*nlocal, &_g_node_id);
	for( i=0; i<nlocal; i++ ) {
		_g_node_id[i] = g_node_id[i]; 
	}
	ISCreateGeneralWithArray( PETSC_COMM_WORLD, nlocal, _g_node_id, &is_gnode );
	VecScatterCreate( g_ownership, is_gnode, local_ownership, PETSC_NULL, &vscat_ownership );
	
	
	/* assign unique equation numbers */
	VecSet( local_ownership, -6699 );
	_VecScatterBeginEnd( vscat_ownership, g_ownership, local_ownership, INSERT_VALUES, SCATTER_FORWARD );
	
	
	/* Count instances of rank in the local_ownership vector */
	VecGetArray( local_ownership, &_local_ownership );
	local_eqnum_count = 0;
	for( i=0; i<nlocal; i++ ) {
		if( ((PetscInt)_local_ownership[i]) == rank ) {
			local_eqnum_count++;
		}
	}
	
	MPI_Allreduce( &local_eqnum_count, &global_eqnum_count, 1, MPI_INT, MPI_SUM, PETSC_COMM_WORLD );
	/* PetscPrintf( PETSC_COMM_SELF, 
	    "[%d] number of local,global equations (without dofs) %d,%d \n", rank, local_eqnum_count, global_eqnum_count ); */
	/* check */
	spanx = NX;
	spany = NY;
	spanz = NZ;
	if( periodic_x==PETSC_TRUE ) {
		spanx--;
	}
	if( periodic_y==PETSC_TRUE ) {
		spany--;
	}
	if( periodic_z==PETSC_TRUE ) {
		spanz--;
	}
	total = spanx*spany*spanz;
	if( total!=global_eqnum_count ) {
		SETERRQ(PETSC_ERR_SUP, "Something stinks. Computed global size for nodes does not match expected" );
	}
	
	
	
	VecCreate( PETSC_COMM_WORLD, &global_eqnum );
	VecSetSizes( global_eqnum, PETSC_DECIDE, nglobal*dof );
	VecSetFromOptions( global_eqnum );
	VecSet( global_eqnum, 0.0 );
	
	/* Load existing eqnums in */
	for( i=0; i<nlocal; i++ ) {
		n_inserts = 0;
		for( d=0; d<dof; d++ ) {
/*
			idx[d] = -(g_node_id[i]*dof + d);
			val[d] = 0.0;
*/
			if( eqnums[ i*dof + d ] < 0.0 ) {
				idx[n_inserts] = g_node_id[i]*dof + d;
				val[n_inserts] = eqnums[ i*dof + d ];
				n_inserts++;
			}
		}
		/* only insert dirichlet bc's */
		VecSetValues( global_eqnum, n_inserts, idx, val, INSERT_VALUES );
	}
	VecAssemblyBegin(global_eqnum);
	VecAssemblyEnd(global_eqnum);
	
	
	
	
	
	/* Generate list of eqnums to get */
	PetscMalloc( sizeof(PetscInt)*nlocal*dof, &to_fetch );
	cnt = 0;
	for( i=0; i<nlocal; i++ ) {
		if( _local_ownership[i]==rank ) {
			for( d=0; d<dof; d++ ) {
				to_fetch[cnt] = g_node_id[i]*dof + d;
				cnt++;
			}
		}
	}
	number_to_fetch = cnt;
	
	
	VecCreate( PETSC_COMM_SELF, &local_eqnum );
	VecSetSizes( local_eqnum, PETSC_DECIDE, number_to_fetch);
	VecSetFromOptions( local_eqnum );
	
	ISCreateGeneralWithArray( PETSC_COMM_SELF, number_to_fetch, to_fetch, &is_eqnum );
	VecScatterCreate( global_eqnum, is_eqnum, local_eqnum, PETSC_NULL, &vscat_eqnum );
	
	VecSet( local_eqnum, -6699 );
	_VecScatterBeginEnd( vscat_eqnum, global_eqnum, local_eqnum, INSERT_VALUES, SCATTER_FORWARD );
	
	
	/* compute offset */
	/* count how many entries there are */
	VecGetArray( local_eqnum, &_local_eqnum );
	eq_cnt = 0;
	for( i=0; i<number_to_fetch; i++ ) {
		if( (PetscInt)_local_eqnum[i]==0 ) {
			eq_cnt++;
		}
	}
	VecRestoreArray( local_eqnum, &_local_eqnum );
	
	VecCreate( PETSC_COMM_WORLD, &offset_list );
	VecSetSizes( offset_list, PETSC_DECIDE, (size+1) );
	VecSetFromOptions( offset_list );
	for( i=rank; i<size; i++ ) {
		VecSetValue( offset_list, i+1, eq_cnt, ADD_VALUES );
	}
	VecAssemblyBegin(offset_list);
	VecAssemblyEnd(offset_list);
	/*
	PetscPrintf(PETSC_COMM_WORLD, "offset_list \n");
	VecView( offset_list, PETSC_VIEWER_STDOUT_WORLD );
	*/
	
	VecScatterCreateToAll(offset_list,&vscat_offset,&seq_offset_list);
	_VecScatterBeginEnd( vscat_offset, offset_list, seq_offset_list, INSERT_VALUES, SCATTER_FORWARD );
	
	{
		PetscScalar *_seq_offset_list;
		
		VecGetArray( seq_offset_list, &_seq_offset_list );
		offset = _seq_offset_list[ rank ];
		VecRestoreArray( seq_offset_list, &_seq_offset_list );
	}
	VecScatterDestroy(vscat_offset);
	VecDestroy(offset_list);
	VecDestroy(seq_offset_list);
	
	/* PetscPrintf( PETSC_COMM_SELF, "[%d]: offset = %d \n", rank, offset ); */
	
	
	VecGetArray( local_eqnum, &_local_eqnum );
	inc = 0;
	for( i=0; i<number_to_fetch; i++ ) {
		if( (PetscInt)_local_eqnum[i]==0 ) {
			_local_eqnum[i] = offset+inc;
			inc++;
		}
	}
	VecRestoreArray( local_eqnum, &_local_eqnum );
	
	_VecScatterBeginEnd( vscat_eqnum, local_eqnum, global_eqnum, INSERT_VALUES, SCATTER_REVERSE );
	/*
	PetscPrintf(PETSC_COMM_WORLD, "global_eqnum \n");
	VecView( global_eqnum, PETSC_VIEWER_STDOUT_WORLD );
	*/
	
	/* For each periodic boundary, get the mapped nodes */
	
	if( periodic_x==PETSC_TRUE ) {
		VecScatter vscat_p;
		IS is_from;
		PetscInt *from, *to;
		Vec mapped;
		PetscScalar *_mapped;
		PetscInt c;
		
		PetscMalloc( sizeof(PetscInt)*npx*dof, &from );
		PetscMalloc( sizeof(PetscInt)*npx*dof, &to );
		
		c = 0;
		for( i=0; i<npx; i++ ) {
			PetscInt I,J,K,gid,from_gid;
			
			gid = periodic_x_gnode_id[i];
			K = gid/(NX*NY);
			J = (gid - K*(NX*NY))/NX;
			I = gid - K*(NX*NY) - J*NX;
			from_gid = (I-(NX-1)) + J*NX + K*(NX*NY);
			
			for( d=0; d<dof; d++ ) {
				to[c] = gid * dof + d;
				from[c] = from_gid * dof + d;
				c++;
			}
		}
		
		
		VecCreate( PETSC_COMM_SELF, &mapped );
		VecSetSizes( mapped, PETSC_DECIDE, npx*dof );
		VecSetFromOptions( mapped );
		
		ISCreateGeneralWithArray( PETSC_COMM_SELF, npx*dof, from, &is_from );
		VecScatterCreate( global_eqnum, is_from, mapped, PETSC_NULL, &vscat_p );
		
		_VecScatterBeginEnd( vscat_p, global_eqnum, mapped, INSERT_VALUES, SCATTER_FORWARD );
		if( npx>0 ) {
			VecGetArray( mapped, &_mapped );
			VecSetValues( global_eqnum, npx*dof, to, _mapped,  INSERT_VALUES );
			VecRestoreArray( mapped, &_mapped );
		}
		
		VecAssemblyBegin(global_eqnum);
		VecAssemblyEnd(global_eqnum);
		
		VecScatterDestroy( vscat_p );
		ISDestroy( is_from );
		VecDestroy( mapped );
		PetscFree( from );
		PetscFree( to );
	}
	
	
	if( periodic_y==PETSC_TRUE ) {
		VecScatter vscat_p;
		IS is_from;
		PetscInt *from, *to;
		Vec mapped;
		PetscScalar *_mapped;
		PetscInt c;
		
		PetscMalloc( sizeof(PetscInt)*npy*dof, &from );
		PetscMalloc( sizeof(PetscInt)*npy*dof, &to );
		
		c = 0;
		for( i=0; i<npy; i++ ) {
			PetscInt I,J,K,gid,from_gid;
			
			gid = periodic_y_gnode_id[i];
			K = gid/(NX*NY);
			J = (gid - K*(NX*NY))/NX;
			I = gid - K*(NX*NY) - J*NX;
			from_gid = I + (J - (NY - 1))*NX + K*(NX*NY);
			
			for( d=0; d<dof; d++ ) {
				to[c] = gid * dof + d;
				from[c] = from_gid * dof + d;
				c++;
			}
		}
		
		
		VecCreate( PETSC_COMM_SELF, &mapped );
		VecSetSizes( mapped, PETSC_DECIDE, npy*dof );
		VecSetFromOptions( mapped );
		
		ISCreateGeneralWithArray( PETSC_COMM_SELF, npy*dof, from, &is_from );
		VecScatterCreate( global_eqnum, is_from, mapped, PETSC_NULL, &vscat_p );
		
		_VecScatterBeginEnd( vscat_p, global_eqnum, mapped, INSERT_VALUES, SCATTER_FORWARD );
		if( npy>0 ) {
			VecGetArray( mapped, &_mapped );
			VecSetValues( global_eqnum, npy*dof, to, _mapped,  INSERT_VALUES );
			VecRestoreArray( mapped, &_mapped );
		}
		
		VecAssemblyBegin(global_eqnum);
		VecAssemblyEnd(global_eqnum);
		
		VecScatterDestroy( vscat_p );
		ISDestroy( is_from );
		VecDestroy( mapped );
		PetscFree( from );
		PetscFree( to );
	}

	if( periodic_z==PETSC_TRUE ) {
		VecScatter vscat_p;
		IS is_from;
		PetscInt *from, *to;
		Vec mapped;
		PetscScalar *_mapped;
		PetscInt c;
		
		PetscMalloc( sizeof(PetscInt)*npz*dof, &from );
		PetscMalloc( sizeof(PetscInt)*npz*dof, &to );
		
		c = 0;
		for( i=0; i<npz; i++ ) {
			PetscInt I,J,K,gid,from_gid;
			
			gid = periodic_z_gnode_id[i];
			K = gid/(NX*NY);
			J = (gid - K*(NX*NY))/NX;
			I = gid - K*(NX*NY) - J*NX;
			from_gid = I + J*NX + (K - (NZ-1))*(NX*NY);
			
			for( d=0; d<dof; d++ ) {
				to[c] = gid * dof + d;
				from[c] = from_gid * dof + d;
				c++;
			}
		}
		
		
		VecCreate( PETSC_COMM_SELF, &mapped );
		VecSetSizes( mapped, PETSC_DECIDE, npz*dof );
		VecSetFromOptions( mapped );
		
		ISCreateGeneralWithArray( PETSC_COMM_SELF, npz*dof, from, &is_from );
		VecScatterCreate( global_eqnum, is_from, mapped, PETSC_NULL, &vscat_p );
		
		_VecScatterBeginEnd( vscat_p, global_eqnum, mapped, INSERT_VALUES, SCATTER_FORWARD );
		if( npz>0 ) {
			VecGetArray( mapped, &_mapped );
			VecSetValues( global_eqnum, npz*dof, to, _mapped,  INSERT_VALUES );
			VecRestoreArray( mapped, &_mapped );
		}
		
		VecAssemblyBegin(global_eqnum);
		VecAssemblyEnd(global_eqnum);
		
		VecScatterDestroy( vscat_p );
		ISDestroy( is_from );
		VecDestroy( mapped );
		PetscFree( from );
		PetscFree( to );
	}

	
	/*
	PetscPrintf(PETSC_COMM_WORLD, "global_eqnum following periodic \n");
	VecView( global_eqnum, PETSC_VIEWER_STDOUT_WORLD );
	*/
	
	
	/* Finally, scatter stuff from global_eqnums into MY array */
	{
		IS is_all_my_eqnums;
		PetscInt *all_my_eqnum_index;
		PetscInt _I,_D,CNT;
		Vec all_my_eqnums;
		PetscScalar *_all_my_eqnums;
		VecScatter vscat_p;
		
		PetscMalloc( sizeof(PetscInt)*dof*nlocal, &all_my_eqnum_index );
		
		CNT = 0;
		for( _I=0; _I<nlocal; _I++ ) {
			for( _D=0; _D<dof; _D++ ) {
				all_my_eqnum_index[CNT] = g_node_id[_I]*dof + _D;
				CNT++;
			}
		}
		
		ISCreateGeneralWithArray( PETSC_COMM_SELF, nlocal*dof, all_my_eqnum_index, &is_all_my_eqnums );
		VecCreate( PETSC_COMM_SELF, &all_my_eqnums );
		VecSetSizes( all_my_eqnums, PETSC_DECIDE, nlocal*dof );
		VecSetFromOptions( all_my_eqnums );
		VecScatterCreate( global_eqnum, is_all_my_eqnums, all_my_eqnums, PETSC_NULL, &vscat_p );
		_VecScatterBeginEnd( vscat_p, global_eqnum, all_my_eqnums, INSERT_VALUES, SCATTER_FORWARD );
		VecGetArray( all_my_eqnums, &_all_my_eqnums );
		for( i=0; i<nlocal*dof; i++ ) {
			eqnums[i] = (int)_all_my_eqnums[i];
		}
		VecRestoreArray( all_my_eqnums, &_all_my_eqnums );
		
		VecScatterDestroy( vscat_p );
		VecDestroy( all_my_eqnums );
		ISDestroy( is_all_my_eqnums );
		PetscFree( all_my_eqnum_index );
	}
	
	VecMax( global_eqnum, &loc, &max );
	*neqnums = (int)max;
	(*neqnums)++;
	
	
	
	/* tidy up */
	VecRestoreArray( local_ownership, &_local_ownership );
	
	
	VecDestroy( g_ownership );
	VecDestroy( local_ownership );
	PetscFree( _g_node_id );
	ISDestroy( is_gnode );
	VecScatterDestroy( vscat_ownership );
	
	VecDestroy( global_eqnum );
	VecDestroy( local_eqnum );
	PetscFree( to_fetch );
	ISDestroy( is_eqnum );
	VecScatterDestroy( vscat_eqnum );
	
	
	
	return 0;
}
