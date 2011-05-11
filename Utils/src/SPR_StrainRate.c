#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include <Underworld/Rheology/Rheology.h>

#include "types.h"
#include "BaseRecoveryFeVar.h"
#include "SPR_StrainRate.h"
#include "REP_Algorithm.h"

#include <math.h>
#include <assert.h>
#include <string.h>

const Type SPR_StrainRate_Type = "SPR_StrainRate";
char*  SPR_NameOfPatch;


SPR_StrainRate* SPR_StrainRate_New(
	Name							name,
	DomainContext*				context,
	void*                   feMesh,
   void*                   geometryMesh,
	DofLayout*              dofLayout,                                                                                 
   void*                   bcs,
   void*                   ics,
   void*                   linkedDofInfo,
   void*                   templateFeVariable,    
	Index							fieldComponentCount,
	Dimension_Index			dim,
	Bool                    isCheckpointedAndReloaded,                                                                 
   Bool                    isReferenceSolution,                                                                       
   Bool                    loadReferenceEachTimestep,    
 	MPI_Comm						communicator,
	FieldVariable_Register*	fV_Register,
	Variable_Register*		vr, 
	FeVariable*					rawField, 
	int							rawOrderOfInterpolation,
	Bool							coeffInterpolation )
{
  SPR_StrainRate* self = (SPR_StrainRate*)_SPR_StrainRate_DefaultNew( name );

	self->isConstructed = True;
   _FieldVariable_Init( (FieldVariable*)self, context, fieldComponentCount, dim, isCheckpointedAndReloaded, communicator, fV_Register );
	_FeVariable_Init( (FeVariable*) self, feMesh, geometryMesh, dofLayout, bcs, ics, linkedDofInfo, templateFeVariable, isReferenceSolution, loadReferenceEachTimestep  );       
   _BaseRecoveryFeVar_Init( (BaseRecoveryFeVar*)self, vr, rawField, rawOrderOfInterpolation, coeffInterpolation );
   _SPR_StrainRate_Init( self );

	return self;
}

void* _SPR_StrainRate_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                              _sizeOfSelf = sizeof(SPR_StrainRate);
	Type                                                      type = SPR_StrainRate_Type;
	Stg_Class_DeleteFunction*                              _delete = _SPR_StrainRate_Delete;
	Stg_Class_PrintFunction*                                _print = _SPR_StrainRate_Print;
	Stg_Class_CopyFunction*                                  _copy = _SPR_StrainRate_Copy;
	Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _SPR_StrainRate_DefaultNew;
	Stg_Component_ConstructFunction*                    _construct = _SPR_StrainRate_AssignFromXML;
	Stg_Component_BuildFunction*                            _build = _SPR_StrainRate_Build;
	Stg_Component_InitialiseFunction*                  _initialise = _SPR_StrainRate_Initialise;
	Stg_Component_ExecuteFunction*                        _execute = _SPR_StrainRate_Execute;
	Stg_Component_DestroyFunction*                        _destroy = _SPR_StrainRate_Destroy;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType                                         nameAllocationType = (AllocationType)ZERO;
	FieldVariable_InterpolateValueAtFunction*             _interpolateValueAt = ZERO;
	FieldVariable_GetValueFunction*               _getMinGlobalFieldMagnitude = ZERO;
	FieldVariable_GetValueFunction*               _getMaxGlobalFieldMagnitude = ZERO;
	FieldVariable_GetCoordFunction*                  _getMinAndMaxLocalCoords = ZERO;
	FieldVariable_GetCoordFunction*                 _getMinAndMaxGlobalCoords = ZERO;
	FeVariable_InterpolateWithinElementFunction*    _interpolateWithinElement = ZERO;
	FeVariable_GetValueAtNodeFunction*                        _getValueAtNode = ZERO;
	FeVariable_SyncShadowValuesFunc*                        _syncShadowValues = ZERO;

	return (void*)_SPR_StrainRate_New(  SPR_STRAINRATE_PASSARGS  );
}

SPR_StrainRate* _SPR_StrainRate_New(  SPR_STRAINRATE_DEFARGS  )
{
	SPR_StrainRate* self;
	
	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(SPR_StrainRate) );
	self = (SPR_StrainRate*) _BaseRecoveryFeVar_New(  BASERECOVERYFEVAR_PASSARGS  );

	return self;
}

void _SPR_StrainRate_Delete(  void* sprVar  ) {
   _BaseRecoveryFeVar_Delete( sprVar );
}

/* --- Virtual Function Implementations --- */
void _SPR_StrainRate_Destroy( void* sprVar, void *data ) {
	_BaseRecoveryFeVar_Destroy( sprVar, data );
}

void _SPR_StrainRate_Print( void* sprVar, Stream* stream ) {
	SPR_StrainRate* self = (SPR_StrainRate*) sprVar;
	
	/* General info */
	Journal_Printf( stream, "SPR_StrainRate (ptr): %p\n", self );
	
	/* Print parent */
	_FeVariable_Print( self, stream );
}


void* _SPR_StrainRate_Copy( const void* sprVar, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
   return NULL;
}

void _SPR_StrainRate_Init( SPR_StrainRate* self ) {
   /* initialise the fieldComponent count here */
   self->fieldComponentCount = StGermain_nSymmetricTensorVectorComponents(self->dim);

	/* Setup basic pointers and functionPtrs that are specific to problem spec  */
   ContextEP_Append_AlwaysLast( self->context, "stokesEqn-execute", _SPR_StrainRate_Execute );

	SPR_NameOfPatch = self->name;
}

void _SPR_StrainRate_AssignFromXML( void* sprVar, Stg_ComponentFactory* cf, void* data ) {
	SPR_StrainRate* self = (SPR_StrainRate*) sprVar;

  /* Construct the parent */
  _BaseRecoveryFeVar_AssignFromXML( self, cf, data );

  _SPR_StrainRate_Init( self );
}

void _SPR_StrainRate_Build( void* sprVar, void* data ) {
	_BaseRecoveryFeVar_Build( sprVar, data );	
}

void _SPR_StrainRate_Initialise( void* sprVar, void* data ) {
	_BaseRecoveryFeVar_Initialise( sprVar, data );	
}

void _SPR_StrainRate_Execute( void* patch, void* data ) {
  /* Fuction Description: is the driver of the main algorithm
   * All repFields associated with this component will be
   * are executed in the below functions 
   */
	PICelleratorContext* context = (PICelleratorContext*)data;
	SPR_StrainRate* self = (SPR_StrainRate*) LiveComponentRegister_Get( context->CF->LCRegister, (Name)SPR_NameOfPatch );
	double startTime;
	assert( self );

	FeVariable_SyncShadowValues( self->rawField );
	MPI_Barrier(MPI_COMM_WORLD);

	startTime = MPI_Wtime( );
	Journal_RPrintf( Journal_Register( Info_Type, (Name)"REP" ), "Start Recovery Method\n" );
		
 	/* Phases of the SPR_StrainRate in this code */

	/* 1: Assemble A and b see eq. 14.31, REFERENCE and solve for [a] */
 	_SPR_StrainRate_AssembleSolveLocalPatchs( self );	
 	/* 2: Communicate [a]  */
	FeVariable_SyncShadowValues( self );
 	/* 3: Calculate the Boundary values */
	_BaseRecoveryFeVar_Boundaries( self );
 	/* 4: Perform a final sync for new proc boundary domain edge nodes */
	FeVariable_SyncShadowValues( self );

	Journal_RPrintf( Journal_Register( Info_Type, (Name)"REP" ), "Time Taken for Recovery Method %f\n", MPI_Wtime()-startTime);
}


void _SPR_StrainRate_AssembleSolveLocalPatchs( void* sprVar ) {
	SPR_StrainRate*		self = (SPR_StrainRate*)sprVar;
	FeVariable*				rawField = self->rawField;
 	BoundaryNodesInfo*	bninfo = self->bninfo;
	Index						dofThatExist = rawField->fieldComponentCount;
	Index						orderOfInterpolation = self->orderOfInterpolation;
	Index						dof_I, nLocalNodes, node_I;
	double					**AMatrix;   /* holds the AMatrix */     
	double					**bVector;   /* holds bVec for each dof */
	double					*patchCoeff;
	

	/* Allocate memory to AMatrix and bVectors */
	AMatrix = Memory_Alloc_2DArray( double, orderOfInterpolation, orderOfInterpolation, (Name)"A matrix"  );
  /* multiple bVectors are needed for each dof */
	bVector = Memory_Alloc_2DArray( double, dofThatExist, orderOfInterpolation, (Name)"b Vector * dof"  );
	patchCoeff = Memory_Alloc_Array( double, orderOfInterpolation * dofThatExist, "tmp coefficient array" );
	
  nLocalNodes = FeMesh_GetNodeLocalSize( self->rawField->feMesh );

  /* loop through patchable nodes assembling AMatrix and bVectors and solving for patch coefficients */
  for( node_I = 0 ; node_I < nLocalNodes ; node_I++ ) {
    if( bninfo[node_I].onMeshBoundary )
      continue;
    ZeroMatrix( AMatrix, orderOfInterpolation, orderOfInterpolation );
    ZeroMatrix( bVector, dofThatExist, orderOfInterpolation );

    _SPR_StrainRate_AssemblePatch( self, node_I, AMatrix, bVector);

    for( dof_I = 0 ; dof_I < dofThatExist ; dof_I++ ) {
      _REP_Algorithm_Solver( AMatrix, bVector[dof_I], orderOfInterpolation );
      /* copy bVector into the corresponding patchCoeff array */
      memcpy( &patchCoeff[dof_I*orderOfInterpolation], bVector[dof_I], orderOfInterpolation*sizeof(double) );
    }
    /* write all dof patch coefficients on the node */
    FeVariable_SetValueAtNode( self, node_I , patchCoeff );
  }

	Memory_Free( AMatrix );
	Memory_Free( bVector );
	Memory_Free( patchCoeff );
	
}

void _SPR_StrainRate_AssemblePatch( SPR_StrainRate* self, int node_I, double** AMatrix, double** bVector) {
/* Functiona Details:
 * This assembles the AMatrix and the bVectors, see eq. 14.31, REFERENCE.
 * First the elements around the patch node are found, then for each element the
 * cooordinate and the strain rate at the super convergent locations are stored.
 */

	FeVariable*			rawField = self->rawField;
	FeMesh*				feMesh = (FeMesh*)rawField->feMesh;
	Index					dofThatExist = rawField->fieldComponentCount;
	Index					nbrEl_I, nbrElementID, nbrElCount;
	Index					count_i, count_j, dof_I;
	double				center[3] = {0.0,0.0,0.0};
	int					orderOfInterpolation = self->orderOfInterpolation;
	IArray*				inc = self->inc; 

	SymmetricTensor	*scp_eps; 
	double				**globalCoord;
	double				**pVec;
	int					*nbrElList;    

	/* 1) find the elements around the node and point to them via the nbrElList*/
	FeMesh_GetNodeElements( feMesh, node_I, inc );
	nbrElCount = IArray_GetSize( inc );
	nbrElList = IArray_GetPtr( inc );

	/* 2) Memory Allocations + Initialisations */
	scp_eps  = Memory_Alloc_Array( SymmetricTensor, nbrElCount, "StrainRate at superconvergent points" );

	globalCoord = Memory_Alloc_2DArray( double, nbrElCount, self->dim, (Name)"Global Coords of superconvergent points"  );

	pVec = Memory_Alloc_2DArray( double, nbrElCount, orderOfInterpolation, (Name)"Holds transformed global coord polynomials" );

	/* 3 ) Now collect information, to go find each elements contribution to the patch 
   * So find the p vectors and strain-rate pseudo vectors for the Ax=b equation
   * */
	for( nbrEl_I = 0 ; nbrEl_I < nbrElCount ; nbrEl_I++ ) {
		nbrElementID = nbrElList[ nbrEl_I ];

		FeMesh_CoordLocalToGlobal( feMesh, nbrElementID, center, globalCoord[nbrEl_I] );
		self->_makePoly( globalCoord[nbrEl_I], pVec[nbrEl_I] );

		FeVariable_InterpolateWithinElement( rawField , nbrElList[nbrEl_I], center, scp_eps[nbrEl_I] );
	}

 	/* Construct A Matrix (Geometric based) and b Vectors (tensor based) */
	for( nbrEl_I = 0 ; nbrEl_I < nbrElCount ; nbrEl_I++ ) {
		for( count_i = 0 ; count_i < orderOfInterpolation ; count_i++ ) {
			for( count_j = 0; count_j < orderOfInterpolation ; count_j++ ) {
				AMatrix[count_i][count_j] += ( pVec[nbrEl_I][count_i] * pVec[nbrEl_I][count_j] );
			}
		}
		for(dof_I = 0 ; dof_I < dofThatExist ; dof_I++ ) {
			for( count_i = 0 ; count_i < orderOfInterpolation ; count_i++ ) {
				bVector[dof_I][ count_i ] += (scp_eps[nbrEl_I][dof_I] * pVec[nbrEl_I][count_i]);
			}
 		}
	}
	Memory_Free( scp_eps );
	Memory_Free( globalCoord );
	Memory_Free( pVec );
}


