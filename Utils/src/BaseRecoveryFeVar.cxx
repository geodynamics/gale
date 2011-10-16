#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include "Underworld/Rheology/Rheology.h"

#include "types.h"
#include "BaseRecoveryFeVar.h"

#include <math.h>
#include <assert.h>
#include <string.h>

const Type BaseRecoveryFeVar_Type = "BaseRecoveryFeVar";

void* _BaseRecoveryFeVar_DefaultNew( Name name ) {
  /* this is an abstract call so this function should never be called */
	assert(0);
}
BaseRecoveryFeVar* _BaseRecoveryFeVar_New(  BASERECOVERYFEVAR_DEFARGS  ) 
{
	BaseRecoveryFeVar* self;
	
	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(BaseRecoveryFeVar) );
	/* The following terms are parameters that have been passed into this function but are being set before being passed onto the parent */
	/* This means that any values of these parameters that are passed into this function are not passed onto the parent function
	   and so should be set to ZERO in any children of this class. */
	nameAllocationType          = NON_GLOBAL;
	_interpolateValueAt         = _FeVariable_InterpolateValueAt;
	_getMinGlobalFieldMagnitude = _FeVariable_GetMinGlobalFieldMagnitude;
	_getMaxGlobalFieldMagnitude = _FeVariable_GetMaxGlobalFieldMagnitude;
	_getMinAndMaxLocalCoords    = _FeVariable_GetMinAndMaxLocalCoords;
	_getMinAndMaxGlobalCoords   = _FeVariable_GetMinAndMaxGlobalCoords;
	_interpolateWithinElement   = _FeVariable_InterpolateNodeValuesToElLocalCoord;
	_getValueAtNode             = _FeVariable_GetValueAtNode;
	_syncShadowValues           = _FeVariable_SyncShadowValues;

	self = (BaseRecoveryFeVar*)
		_FeVariable_New(  FEVARIABLE_PASSARGS  );

	self->_getCoeffAtNode = _BaseRecoveryFeVar_GetCoeffAtNode;
	self->_interpolateWithinElement = NULL;

	return self;
}

void _BaseRecoveryFeVar_Delete( void* _self ) {
	BaseRecoveryFeVar* self = (BaseRecoveryFeVar*)_self;

	_FeVariable_Delete( self );
}
  
void _BaseRecoveryFeVar_Destroy( void* _self, void* data ) { 
	BaseRecoveryFeVar* self = (BaseRecoveryFeVar*)_self;
	Memory_Free( self->pVec );
	Memory_Free( self->bninfo );

	Stg_Component_Destroy( self->dataVariable, data, False );
	Stg_Component_Destroy( self->dofLayout, data, False );
	Stg_Component_Destroy( self->rawField, data, False );

	_FeVariable_Destroy( _self, data ); 
} 

void _BaseRecoveryFeVar_Execute( void* _self, void* data ) { _FeVariable_Execute( _self, data ); }

void _BaseRecoveryFeVar_Init(
	BaseRecoveryFeVar*	self, 
	Variable_Register*	vr, 
	FeVariable*				rawField, 
	int						rawOrderOfInterpolation,
	Bool						coeffInterpolation )
{

   self->inc = IArray_New(); /* Must call here because this component never calls _FeVariable_Init */
   self->vr = vr;
   self->rawField = rawField;

   /* record the rawField that is required on the repRequired_RawField register */
   if( (unsigned)-1 == Stg_ObjectList_GetIndex( repRequiredRawFields_Reg, rawField->name ) )
      Stg_ObjectList_Append( repRequiredRawFields_Reg, rawField ); 

  /* we save the rawOrderOfInterpolation here and allow the instances of this class
   * to change their orderOfInterpolation values individually in there Init */
  self->orderOfInterpolation = rawOrderOfInterpolation;

  /* set this components mesh to equal the rawField's mesh */
  self->feMesh = rawField->feMesh;

	/* set polynomial length: linear only available because we only use
	 * linear elements and patches containing immdeiately adjacent elements only */
	if( self->dim == 2 ) {
		self->orderOfInterpolation = 3;
		self->_makePoly = _BaseRecoveryFeVar_pVec_2Dorder1;
	} else {
		self->orderOfInterpolation = 4;
		self->_makePoly = _BaseRecoveryFeVar_pVec_3Dorder1;
	}

	/* override FeVariable functions here */
	self->_getValueAtNode = _BaseRecoveryFeVar_GetValueAtNode;

	/* select interpolation technique */
	if( coeffInterpolation ) 
		self->_interpolateWithinElement = _BaseRecoveryFeVar_GetValueInElementWithCoeffInterpolation;
	else
		self->_interpolateWithinElement = _BaseRecoveryFeVar_GetValueInElementWithStdInterpolation;
}

void _BaseRecoveryFeVar_AssignFromXML( void* _self, Stg_ComponentFactory* cf, void* data ) {
  BaseRecoveryFeVar*	self = (BaseRecoveryFeVar*) _self;
  FeVariable*			rawField = NULL;
  int						rawOrderOfInterpolation;
  Bool					coeffInterpolation;
  Variable_Register*	vr;

  /* Construct FieldVariable, NOT FeVariable parent */
  _FieldVariable_AssignFromXML( self, cf, data );

  vr = self->context->variable_Register;

  /* get the initial field */
  rawField = Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"RawField", FeVariable, True, data  ); 

  /* get the order of interpolation */
  rawOrderOfInterpolation = Stg_ComponentFactory_GetInt( cf, self->name , (Dictionary_Entry_Key)"orderOfInterpolation", 1  );

  coeffInterpolation = Stg_ComponentFactory_GetBool( cf, self->name , (Dictionary_Entry_Key)"coeffInterpolation", True );

  Journal_Firewall( (rawOrderOfInterpolation < 2  ) ,
      Journal_Register( Error_Type, (Name)"BadNews"  ),
      "Error in %s, the orderOfInterpolation must be: 1, higher order coefficients haven't been implemented. Currently it's %d\n"
      "The larger the number the orderOfInterpolation:\n"
      " * the more accurate the recoveredField\n"
      " * the longer an iteration will take\n\n", __func__, rawOrderOfInterpolation );

  _BaseRecoveryFeVar_Init( self, vr, rawField, rawOrderOfInterpolation, coeffInterpolation );
}
void _BaseRecoveryFeVar_Build( void* _self, void* data ) {
  BaseRecoveryFeVar* self = (BaseRecoveryFeVar*) _self;
	Sync*	               sync;
	uint                  componentsCount = 0;
	uint                  variable_I;
	int                  *nodeDomainCountPtr = NULL;
	char                 **variableName, *tmpName;
	Variable*            dataVariable;
	Variable_Register*	 variable_Register = self->vr;

	Stg_Component_Build( self->feMesh, data, False );
	/**** Create FeVariable to store the recovered field ****/

	/* 1st Make sure rawField is built */
	Stg_Component_Build( self->rawField, data, False );

	/* 2nd create vector to store the repField */ 
	componentsCount = self->orderOfInterpolation * self->fieldComponentCount;
  /* create an variable names - create an extra variableName for the dofLayout */
	variableName = Memory_Alloc_Array( char*, componentsCount * sizeof(char*) + 1, "names");
	tmpName = Stg_Object_AppendSuffix( self, (Name)self->name  );

  /* setup the sync and the nodeDomainCount */
	sync = Mesh_GetSync( self->feMesh, MT_VERTEX );
	assert( Class_IsSuper( self->feMesh->topo, IGraph ) );
	nodeDomainCountPtr = &((IGraph*)self->feMesh->topo)->remotes[MT_VERTEX]->nDomains;
	assert( componentsCount <= 100 );

	for ( variable_I = 0 ; variable_I < componentsCount ; variable_I++ ) {
		Stg_asprintf( &variableName[ variable_I ], "%s_%s-%d", self->type, self->name, variable_I );
	}
	Stg_asprintf( &variableName[ componentsCount ], "%s-dofLayout", self->name );
	dataVariable = Variable_NewVector2( tmpName,
												(AbstractContext*)self->context,
												Variable_DataType_Double, 
												componentsCount,
												(unsigned*)&sync->nDomains, 
												NULL, 
												(void**)NULL,
												variable_Register,
                                            (const char**)variableName );

	for( variable_I = 0; variable_I < dataVariable->dataTypeCounts[0]; variable_I++ )
		dataVariable->components[variable_I]->allocateSelf = True;

	dataVariable->allocateSelf = True;

	/* 3rd create the dofLayout */
	self->dofLayout = DofLayout_New( variableName[componentsCount], self->context, variable_Register, *nodeDomainCountPtr, self->feMesh );
	for ( variable_I = 0 ; variable_I < componentsCount ; variable_I++ ) {

			/* We have to set the array ptr ptr for these guys manually - this should be fixed */
			Variable* variable = Variable_Register_GetByName( variable_Register, variableName[ variable_I ] );
			variable->arrayPtrPtr = &dataVariable->arrayPtr;

			/* Assign variable to each node */
			for( int node_I = 0; node_I < *nodeDomainCountPtr; node_I++ ) {
				DofLayout_AddDof_ByVarName( self->dofLayout, variableName[variable_I], node_I );
			}
			/* Free Name */
			Memory_Free( variableName[ variable_I ] );
		}
	
	Stg_Component_Build( self->dofLayout, NULL, False );
	Stg_Component_Initialise( self->dofLayout, NULL, False );

	LiveComponentRegister_Add( LiveComponentRegister_GetLiveComponentRegister(), (Stg_Component*) dataVariable );
	LiveComponentRegister_Add( LiveComponentRegister_GetLiveComponentRegister(), (Stg_Component*) self->dofLayout );

   self->dataVariable = dataVariable;
	/**** END building repField ****/

  self->bninfo = Memory_Alloc_Array( BoundaryNodesInfo, 
						      FeMesh_GetNodeDomainSize( self->feMesh ),
						      "BoundaryNodesInfo" );
}
void _BaseRecoveryFeVar_Initialise( void* _self, void* data ) {
	BaseRecoveryFeVar*      self = (BaseRecoveryFeVar*) _self;

	/* Initialise parent */
	_FieldVariable_Initialise( self, data );

	/* Initialise class specific stuff */
	Stg_Component_Initialise( self->rawField, data, False );
	Stg_Component_Initialise( self->dofLayout, data, False );

	self->pVec = Memory_Alloc_Array( double, self->orderOfInterpolation, "position polynomial" );
  
	/* BaseUtils_PopulateBoundaryNodesInfo should be called each timestep if AMR is used.
	 * Here it's assumed the domain is not REMESHED after the it's created */
	BaseUtils_PopulateBoundaryNodesInfo( self->rawField->feMesh, self->bninfo );
}

void _BaseRecoveryFeVar_Boundaries( void* _self ) {
	BaseRecoveryFeVar* self         = (BaseRecoveryFeVar*)_self;
	FeMesh*             feMesh   = self->feMesh;
  BoundaryNodesInfo*  bninfo = self->bninfo;
	Index               coeffCount, nPatches, nLocalNodes, coeff_I, node_I, patch_I, nbrNodeID;
  double *coord, *nbrCoord, distance, weight; 
  double numerator[50], denominator[50], coeff[50], nbrCoeff[50];

  coeffCount = self->orderOfInterpolation * self->fieldComponentCount;
  nLocalNodes = FeMesh_GetNodeLocalSize( self->rawField->feMesh );

  for( node_I = 0 ; node_I < nLocalNodes ; node_I++ ) {
    if( !bninfo[node_I].onMeshBoundary )
      continue;
    memset( numerator, 0, sizeof(double)*50 );
    memset( denominator, 0, sizeof(double)*50 );
    coord =  Mesh_GetVertex( feMesh , node_I );
    nPatches = bninfo[node_I].numOfPatches2use;
    /* go through nodes with patches that will contribute to this boundary node */
    for( patch_I = 0 ; patch_I < nPatches; patch_I++ ) {
      nbrNodeID = bninfo[node_I].patchNodes[patch_I];
      nbrCoord = Mesh_GetVertex( feMesh, nbrNodeID );

      /* weight the coefficients to be interpolated */
      distance = StGermain_DistanceBetweenPoints( nbrCoord, coord , (self->dim) );
      weight = 1/(distance*distance);
      _BaseRecoveryFeVar_GetCoeffAtNode( self, nbrNodeID, nbrCoeff );

      for( coeff_I = 0 ; coeff_I < coeffCount ; coeff_I++ ) {
        numerator[coeff_I] += weight*nbrCoeff[coeff_I];
        denominator[coeff_I] += weight;
      }
    }

    for( coeff_I = 0 ; coeff_I < coeffCount ; coeff_I++ ) 
      coeff[coeff_I] = numerator[coeff_I]/denominator[coeff_I];

    /* set value at boundary node */
    FeVariable_SetValueAtNode( self, node_I , coeff );
  }
}

void _BaseRecoveryFeVar_GetCoeffAtNode( void* feVariable, Node_DomainIndex dNode_I, double* coeff ) {
	/**** Get Coefficients ****/
	BaseRecoveryFeVar* self = (BaseRecoveryFeVar*) feVariable;
	Variable*	currVariable = NULL;
	Dof_Index	dofCountThisNode = 0;
	Dof_Index	nodeLocalDof_I;

	dofCountThisNode = self->dofLayout->dofCounts[dNode_I];
	
	for ( nodeLocalDof_I=0; nodeLocalDof_I < dofCountThisNode; nodeLocalDof_I++ ) {
		currVariable = DofLayout_GetVariable( self->dofLayout, dNode_I, nodeLocalDof_I );
		coeff[ nodeLocalDof_I ] = Variable_GetValueDouble( currVariable, dNode_I );
	}
}

void _BaseRecoveryFeVar_GetValueAtNode( void* feVariable, Node_DomainIndex dNode_I, double* value ) {
	BaseRecoveryFeVar* self = (BaseRecoveryFeVar*) feVariable;
	int dof_I, order, dofThatExist;
	double coeff[50];
	double *ptr = coeff;
	double *coord;

	/**** Get Coefficients ****/
	Variable*	currVariable = NULL;
	Dof_Index	dofCountThisNode = 0;
	Dof_Index	nodeLocalDof_I;

	dofCountThisNode = self->dofLayout->dofCounts[dNode_I];
	
	for ( nodeLocalDof_I=0; nodeLocalDof_I < dofCountThisNode; nodeLocalDof_I++ ) {
		currVariable = DofLayout_GetVariable( self->dofLayout, dNode_I, nodeLocalDof_I );
		ptr[ nodeLocalDof_I ] = Variable_GetValueDouble( currVariable, dNode_I );
	}

	/**** Apply Coefficients ****/
	order = self->orderOfInterpolation;
	dofThatExist = self->fieldComponentCount;
	coord = Mesh_GetVertex( self->feMesh , dNode_I );

	for( dof_I = 0 ; dof_I < dofThatExist ; dof_I++ ) 
			value[dof_I] = _BaseRecoveryFeVar_ApplyCoeff( self, &(coeff[dof_I*order]), coord, order, self->pVec );
}

void _BaseRecoveryFeVar_GetValueInElementWithStdInterpolation( void* feVariable, Element_Index lEl_I, double* xi, double* value ) {
	BaseRecoveryFeVar  *self = (BaseRecoveryFeVar*) feVariable;
	FeMesh             *mesh = self->feMesh;
	ElementType        *elementType = NULL;
	IArray*            inc  = self->inc;
	double             valueAtNode[27][9];  /* TODO: make non static */
	double             Ni[27]; 
	int nNodes, *nodes, node_I, dof_I, dofThatExist;

	dofThatExist = self->fieldComponentCount;

	/* Get nodes in element lEl_I */
	FeMesh_GetElementNodes( mesh, lEl_I, self->inc );
	nNodes = IArray_GetSize( inc );
	nodes = IArray_GetPtr( inc );

	/* Calculate shapefunction for interpolation */
	elementType = FeMesh_GetElementType( mesh, lEl_I );
	ElementType_EvaluateShapeFunctionsAt( elementType, xi, Ni );

	/* zero the value vector */
	memset( value, 0, sizeof(double)*dofThatExist );

	/* Get value of self an the nodes and store in valueAtNode */ 
	for( node_I = 0 ; node_I < nNodes ; node_I++ ) {
		_BaseRecoveryFeVar_GetValueAtNode( self, nodes[ node_I ], valueAtNode[ node_I ] );
		/* Interpolate valueAtNode to xi manually */
		for( dof_I = 0 ; dof_I < dofThatExist ; dof_I++ ) {
				value[dof_I] += Ni[node_I]*valueAtNode[node_I][dof_I];
		}
	}
}

void _BaseRecoveryFeVar_GetValueInElementWithCoeffInterpolation( void* feVariable, Element_Index lEl_I, double* xi, double* value ) {
	BaseRecoveryFeVar* self = (BaseRecoveryFeVar*) feVariable;
	FeMesh*              mesh = self->feMesh;
	double               globalCoord[3];
	double               coeff[50]; 
	int order, dof_I, dofThatExist;

	order = self->orderOfInterpolation;
	dofThatExist = self->fieldComponentCount;
	FeMesh_CoordLocalToGlobal( mesh, lEl_I, xi, globalCoord );

	/* Get nodes in element lEl_I */
	FeMesh_GetElementNodes( mesh, lEl_I, self->inc );

	/** Get the coeffieients of each node **/
	_FeVariable_InterpolateNodeValuesToElLocalCoord( self, lEl_I, xi, coeff );

	for( dof_I = 0 ; dof_I < dofThatExist ; dof_I++ ) 
			value[dof_I] = _BaseRecoveryFeVar_ApplyCoeff( self, &(coeff[dof_I*order]), globalCoord, order, self->pVec );
}

double _BaseRecoveryFeVar_ApplyCoeff( BaseRecoveryFeVar* self, double* coeff, double* coord, int order, double* pVec ) {
	double value = 0;
	int order_I;

	self->_makePoly( coord, pVec );
	for( order_I = 0 ; order_I < order ; order_I++ ) 
		value += coeff[order_I] * pVec[order_I];

	return value;
}

void BaseUtils_Add2LmStruct( LmStruct* lmStruct, int nodeID ) {
	int count_I;
	for( count_I = 0 ; count_I < lmStruct->numberOfNodes ; count_I++ ) {
		/* To prevent duplication in the list*/
		if( lmStruct->nodeIDList[count_I] == nodeID )
			return;
	}
	lmStruct->nodeIDList[ lmStruct->numberOfNodes ] = nodeID;
	lmStruct->numberOfNodes++;
}

void BaseUtils_PopulateBoundaryNodesInfo( FeMesh* mesh, BoundaryNodesInfo* bninfo ) {
/* Function Description:
 * Populates a mesh-related data structure, called self->boundaryNodeInfo.
 * An array of size n, where n is the number of DOMAIN nodes. 
 * each element contain:
 * Bool onMeshBoundary   ... is node on the Mesh boundary
 * int  numOfPatches2use ... number of Patches which will contribute to node
 *  This will aid the patch algorithms, which needs to know if a node
 *  is on the domain boundary or not (i.e. does not have a full set of valid adjacent elements).
 *  Algorithm to assign onBoundary is
 *  To go through each node
 *  	set onBoundary = false
 *  	check if a single neighbour element is invalid (thus outside the domain)
 *  		if so, set onBoundary = true 
 * 
 */
	Node_Index              tmpNodeID, tmpNode_I;
	Element_Index           nbrElementID;
        int                     *els;
	uint                    nNodes;
	LmStruct                list;
	Sync*                   sync;
	uint                    nLocalNodes;
	int                     nVerts;
        const int               *verts;
	IArray*			            inc[2];

	sync = Mesh_GetSync( mesh, MT_VERTEX );
	IGraph_GetBoundaryElements( mesh->topo, MT_VERTEX, &nVerts, &verts );
	nLocalNodes = FeMesh_GetNodeLocalSize( mesh );
	int domainNodes = FeMesh_GetNodeDomainSize( mesh );

	assert( nVerts <= domainNodes );

  /* initialise boundaryNodesInfo data structure:
   * for all nodes
   * onMeshBoundary = false,
   * numOfPatches2use = 0
   * patchNodes list = [-1,-1,-1,....,-1]
   * */
	for( int dNodes_I = 0 ; dNodes_I < domainNodes ; dNodes_I++ ) {
		bninfo[dNodes_I].onMeshBoundary = False;
		bninfo[dNodes_I].numOfPatches2use = 0;
		memset( bninfo[dNodes_I].patchNodes, -1, sizeof(int)*REP_MAXNODESPERPATCH );
	}

	/* First flag boundary nodes */
	for( int dNodes_I = 0 ; dNodes_I < nVerts ; dNodes_I++ )
	       bninfo[verts[dNodes_I]].onMeshBoundary = True;

	/* Update all other procs. */
    /* Argument Describtion:
     * 1) start of local address
     * 2) size of chunk in array
     * 3) location and size of local information in array
     * 4) size of chunk to send???
     * 5) size of chunk to receive??
     */
	Sync_SyncArray( sync, bninfo, sizeof(BoundaryNodesInfo), 
			bninfo + nLocalNodes , sizeof(BoundaryNodesInfo), 
			sizeof(BoundaryNodesInfo) );	

	/* Now calculate the numOfPatches2use for each onMeshBoundary node
	 * If node is on boundary, find how many patches will contribute to recovery at that location
	 * 	 ONLY LOOP OVER LOCAL NODES HERE
	 */
	inc[0] = IArray_New();
	inc[1] = IArray_New();
	for( uint dNodes_I = 0 ; dNodes_I < FeMesh_GetNodeDomainSize( mesh ) ; dNodes_I++ ) {
		if( !bninfo[dNodes_I].onMeshBoundary ) 
			continue;

		list.numberOfNodes = 0;
				
		/* Go through all neighbour elements */
		FeMesh_GetNodeElements( mesh, dNodes_I, inc[0] );
		int nEls = IArray_GetSize( inc[0] );
		els = IArray_GetPtr( inc[0] );

		for(int nbrElement_I = 0; nbrElement_I < nEls; nbrElement_I++) {
			nbrElementID = els[nbrElement_I];
			/* Go through nodes on elements and see if they're valid patchs */
			FeMesh_GetElementNodes( mesh, nbrElementID, inc[1] );
			nNodes = IArray_GetSize( inc[1] );
			int *nodes = IArray_GetPtr( inc[1] );

			for( tmpNode_I = 0 ; tmpNode_I < nNodes ; tmpNode_I++ ) {
				tmpNodeID = nodes[tmpNode_I];
				if( !bninfo[tmpNodeID].onMeshBoundary ) 
					BaseUtils_Add2LmStruct( &list, tmpNodeID );
			}
		}
		bninfo[dNodes_I].numOfPatches2use = list.numberOfNodes;
		assert( list.numberOfNodes < 20 ); /* TODO: not pretty */
		memcpy( bninfo[dNodes_I].patchNodes, list.nodeIDList, list.numberOfNodes*sizeof(int) );
	}
	NewClass_Delete( inc[0] );
	NewClass_Delete( inc[1] );
}


void _BaseRecoveryFeVar_pVec_2Dorder1( double *global, double *pVec ) {
	pVec[0] = 1;
	pVec[1] = global[0];
	pVec[2] = global[1];
}

void _BaseRecoveryFeVar_pVec_3Dorder1( double *global, double *pVec ) {
	pVec[0] = 1;
	pVec[1] = global[0];
	pVec[2] = global[1];
	pVec[3] = global[2];
}


