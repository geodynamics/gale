/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
** Copyright (c) 2007, Monash Cluster Computing 
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
**
** Author:
**		Julian Giordani - Julian.Giordani@sci.monash.edu.au
**
** Contributors:
*+		Robert Turnbull
*+		Vincent Lemiale
*+		Louis Moresi
*+		David May
*+		David Stegman
*+		Mirko Velic
*+		Patrick Sunter
*+
** $Id:  $
** 
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/*
 * The REP algorithm is desciped in the paper
 * 	B.Boroomand & O.C.Zienkiewicz, "An Improved REP Recovery and the Effectivity Robustness Test",
 * 	Int. J. for Numerical Methods in Engineering, vol. 40, pages 3247-3277, 1997.
 * 
 *
 * Phases of the REP_Algorithm in this code
 * 	1: Assemble element H and F
 * 	2: Communicate H and F elemental objects
 * 	3: Assemble patch matrix, A and rightside vector b
 * 	4: Solve Ax=b, where x is the patch Coefficients
 * 	5: Communicate proc boundary coefficient for completeness with respect to domain recoveries
 *  6: Perform a final sync of the node values at proc-boundary domain-edge nodes <br/>
 */ 

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include <Underworld/Underworld.h>

#include "types.h"
#include "BaseRecoveryFeVar.h"
#include "RecoveredFeVariable.h"
#include "REP_Algorithm.h"

#include <assert.h>
#include <string.h>

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type REP_Algorithm_Type = "REP_Algorithm";
char*  NameOfPatch;
void* _REP_Algorithm_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                              _sizeOfSelf = sizeof(REP_Algorithm);
	Type                                                      type = REP_Algorithm_Type;
	Stg_Class_DeleteFunction*                              _delete = _REP_Algorithm_Delete;
	Stg_Class_PrintFunction*                                _print = _REP_Algorithm_Print;
	Stg_Class_CopyFunction*                                  _copy = NULL;
	Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _REP_Algorithm_DefaultNew;
	Stg_Component_ConstructFunction*                    _construct = _REP_Algorithm_AssignFromXML;
	Stg_Component_BuildFunction*                            _build = _REP_Algorithm_Build;
	Stg_Component_InitialiseFunction*                  _initialise = _REP_Algorithm_Initialise;
	Stg_Component_ExecuteFunction*                        _execute = _REP_Algorithm_Execute;
	Stg_Component_DestroyFunction*                        _destroy = _REP_Algorithm_Destroy;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

	return (void*)_REP_Algorithm_New(  REP_ALGORITHM_PASSARGS  );
}

/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
REP_Algorithm* _REP_Algorithm_New(  REP_ALGORITHM_DEFARGS  )
{
	REP_Algorithm* self;
	
	/* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
	assert( _sizeOfSelf >= sizeof(REP_Algorithm) );
	/* The following terms are parameters that have been passed into this function but are being set before being passed onto the parent */
	/* This means that any values of these parameters that are passed into this function are not passed onto the parent function
	   and so should be set to ZERO in any children of this class. */
	nameAllocationType = NON_GLOBAL;

	self = (REP_Algorithm*) _Stg_Component_New(  STG_COMPONENT_PASSARGS  );

	return self;
}


void _REP_Algorithm_Delete( void* patchRecoveryFeVariable ) {
}

void _REP_Algorithm_Print( void* patchRecoveryFeVariable, Stream* stream ) {

}

void _REP_Algorithm_Init( void* rep, 
			  UnderworldContext* context )
{
   REP_Algorithm* self = (REP_Algorithm*)rep;

#if PATCH_DEBUG
   char* filename;
   self->myStream = Journal_Register( Info_Type, (Name)self->type  );
   Stg_asprintf( &filename, "FieldRecovery-Output.%dof%d.dat", context->rank, context->nproc );
   Stream_RedirectFile_WithPrependedPath( self->myStream, context->outputPath, filename );
   Memory_Free( filename );
#endif

   self->context = context;
   ContextEP_Append_AlwaysLast( context, "stokesEqn-execute", _REP_Algorithm_Execute );
   NameOfPatch = self->name;

   self->incArray = IArray_New();
}


void _REP_Algorithm_AssignFromXML( void* patchRecoveryFeVariable, Stg_ComponentFactory* cf, void* data ) {
   REP_Algorithm*          self     = (REP_Algorithm*)patchRecoveryFeVariable;
   UnderworldContext*      context = NULL;
   Dictionary*             dictionary = Dictionary_GetDictionary( cf->componentDict, self->name );
   Dictionary_Entry_Value* list;
   char*                   varName;
   Stream*                 errorStream = Journal_Register( Error_Type, (Name)"_REP_Algorithm_Construct"  );
   int field_I, listCount;

	self->IPswarm         = NULL;

	context = Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"Context", UnderworldContext, False, data );
	if( !context  )
		context = Stg_ComponentFactory_ConstructByName( cf, (Name)"context", UnderworldContext, True, data  );

  list = Dictionary_Get( dictionary, (Dictionary_Entry_Key)"RepFieldList"  );
  Journal_Firewall(
			list != NULL,
			errorStream,
			"Error in %s:\n"
			"You must specify a list of RecoveredFeVariable in your xml.\nExample:\n"
			"<list name=\"RepFieldList\">\n"
				"\t<param>recoveredSigmaField</param>\n"
			"</list>\n", __func__ );
	
	listCount = Dictionary_Entry_Value_GetCount( list );
  Journal_Firewall(
			listCount != 0,
			errorStream,
			"Error in %s:\n"
			"You have no RecoveredFeVariable defined in the RepFieldList list. At least one RecoveredFeVariable must be included.\nExample:\n"
			"<list name=\"RepFieldList\">\n"
				"\t<param>recoveredSigmaField</param>\n"
			"</list>\n", __func__ );
	
	/* Allocate the memory to store pointers */
	self->repFieldList = Memory_Alloc_Array( RecoveredFeVariable*, listCount, "List RecoveredFeVariable" );

	/* get all the rep fields from the dictionary */
	for( field_I = 0 ; field_I < listCount ; field_I++ ) {
		varName = Dictionary_Entry_Value_AsString( Dictionary_Entry_Value_GetElement( list, field_I ) );
		self->repFieldList[ field_I ] = Stg_ComponentFactory_ConstructByName( cf, (Name)varName, RecoveredFeVariable, True, data  );
	}

	self->constitutiveMatrix = Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"ConstitutiveMatrix", ConstitutiveMatrix, True, data );

	if(self->constitutiveMatrix ) { 
    /* implicitly get from constitutiveMatrix dictionary */
		self->IPswarm = Stg_ComponentFactory_ConstructByKey( cf, self->constitutiveMatrix->name, (Dictionary_Entry_Key)"Swarm", IntegrationPointsSwarm, True, data  );

		/* check if any of the repFields require the storeConstitutiveMatrix, if so turnOn */
		for( field_I = 0 ; field_I < listCount ; field_I++ ) {
			if( self->repFieldList[field_I]->nonLinearProblem == True ) { 
				self->constitutiveMatrix->storeConstitutiveMatrix = True;
				ConstitutiveMatrixCartesian_SetupParticleStorage( (ConstitutiveMatrixCartesian*)self->constitutiveMatrix );
				/* only need to do the above once if we need it */
				break;
			}
		}
	}
	else {
    /* explicitly get from component dictionary */
		self->IPswarm = Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"IntegrationPoints", IntegrationPointsSwarm, True, data  );
	}

  self->repFieldCount = listCount;

	_REP_Algorithm_Init( self, context );
}

void _REP_Algorithm_Build( void* patchRecoveryFeVariable, void* data ) {
   REP_Algorithm* self = (REP_Algorithm*)patchRecoveryFeVariable;
   int xx;
	
   Stg_Component_Build( self->IPswarm, data, False );
   Stg_Component_Build( self->constitutiveMatrix, data, False );
   for( xx = 0 ; xx < self->repFieldCount ; xx++ )
      Stg_Component_Build( self->repFieldList[xx], data, False );
}

void _REP_Algorithm_Initialise( void* patchRecoveryFeVariable, void* data ) {
	REP_Algorithm* self = (REP_Algorithm*)patchRecoveryFeVariable;
	int xx;
	
   Stg_Component_Initialise( self->IPswarm, data, False );
   Stg_Component_Initialise( self->constitutiveMatrix, data, False );
	for( xx = 0 ; xx < self->repFieldCount ; xx++ )
		Stg_Component_Initialise( self->repFieldList[xx], data, False );
}

void _REP_Algorithm_Execute( void* patch, void* data ) {
  /* Fuction Description: is the driver of the main algorithm
   * All repFields associated with this component will be
   * are executed in the below functions 
   */
	UnderworldContext* context = (UnderworldContext*)data;
	REP_Algorithm* self = (REP_Algorithm*)LiveComponentRegister_Get( context->CF->LCRegister, (Name)NameOfPatch );
	double startTime;
	int field_I;
	assert( self );

	/* TODO: Optimisation question. Is this a test and do function or a straigh do function */
   for( field_I = 0 ; field_I < Stg_ObjectList_Count( repRequiredRawFields_Reg ); field_I++ )
      FeVariable_SyncShadowValues( (FeVariable* )Stg_ObjectList_At( repRequiredRawFields_Reg, field_I ) );

	MPI_Barrier(MPI_COMM_WORLD);

	startTime = MPI_Wtime();
	Journal_RPrintf( Journal_Register( Info_Type, (Name)"REP" ),
			"Start Recovery Method\n" );
		
  /* Phases of the REP_Algorithm in this code */

  /* 1: Assemble element H and F */
    _REP_Algorithm_AssembleAllLocalElement( self );	
  /* 2: Communicate H and F elemental objects */
    _REP_Algorithm_Communicate( self );
  /* 3: Assemble patch matrix, A and rightside vector b and solve for x*/
    _REP_Algorithm_SolvePatches( self );
  /* 4: Communicate proc boundary coefficient for completeness with respect to domain recoveries */
    _REP_Algorithm_CommunicateBoundaries( self );
  /* 5: Calculate the Boundary values */
    _REP_Algorithm_Boundaries( self );
  /* 6: Perform a final sync for new proc boundary domain edge nodes */
    _REP_Algorithm_CommunicateBoundaries( self );

		Journal_RPrintf( Journal_Register( Info_Type, (Name)"REP" ),
				"Time Taken for Recovery Method %f\n", MPI_Wtime()-startTime);
}

void _REP_Algorithm_Destroy( void* rep, void* data ) {
  REP_Algorithm* self = (REP_Algorithm*)rep;
  int xx;

   Stg_Component_Destroy( self->IPswarm, data, False );
   Stg_Component_Destroy( self->constitutiveMatrix, data, False );
	NewClass_Delete( self->incArray );
	
	for( xx = 0 ; xx < self->repFieldCount ; xx++ )
		Stg_Component_Destroy( self->repFieldList[xx], data, False );

	Memory_Free( self->repFieldList );
}

void _REP_Algorithm_CommunicateBoundaries( REP_Algorithm* self ) {
	int field_I;
	for( field_I = 0 ; field_I < self->repFieldCount ; field_I++ ) 
		_FeVariable_SyncShadowValues( self->repFieldList[field_I] );
}

void _REP_Algorithm_Boundaries( REP_Algorithm* self ) {
	int field_I;
	for( field_I = 0 ; field_I < self->repFieldCount ; field_I++ ) 
		_BaseRecoveryFeVar_Boundaries( self->repFieldList[field_I] );

}

void _REP_Algorithm_SolvePatches( REP_Algorithm* self ) {
	int elList[20]/*TODO is maxNumberOfElementAroundApatch*/, nEl, node_I, nLocalNodes, field_I;
	LmStruct lmStruct;

	nLocalNodes = FeMesh_GetNodeLocalSize( self->repFieldList[0]->feMesh );

	for( field_I = 0 ; field_I < self->repFieldCount ; field_I++ ) 
		RecoveredFeVariable_SetupWorkSpace( self->repFieldList[field_I] );

	for( node_I = 0 ; node_I < nLocalNodes ; node_I++ ) {
		if( self->repFieldList[0]->bninfo[node_I].onMeshBoundary ) /*TODO: only uses 1st mesh*/
			continue;
		REP_Algorithm_MakeLM( self->repFieldList[0]->feMesh, node_I, elList, &nEl, &lmStruct );
		/* for each recovered field */
		for( field_I = 0 ; field_I < self->repFieldCount ; field_I++ )
			RecoveredFeVariable_SolvePatch( self->repFieldList[field_I], node_I, elList, nEl, &lmStruct );
	}
	for( field_I = 0 ; field_I < self->repFieldCount ; field_I++ )
		RecoveredFeVariable_RemoveWorkSpace( self->repFieldList[field_I] );
}
void _REP_Algorithm_Communicate( REP_Algorithm* self ) {
	int field_I;
	for( field_I = 0 ; field_I < self->repFieldCount ; field_I++ )
		RecoveredFeVariable_CommunicateHF( self->repFieldList[field_I] );
}

void REP_Algorithm_MakeLM( FeMesh* mesh, int lPatchID, int *elList, int* nEl, LmStruct* lmStruct ) {
	/* Fuction Description: populates the 
	 * elList = the elements surrounding patchNode
	 * nEl = the number of elements surrounding 
	 * lm  = location matrix (ordering) needed of Sys Linear Eq. REP uses
	 */
	int   el_I, node_I, nodesPerEl, *nodeList;
	Index dElID, dNodeID;
	IArray *inc_E, *inc_N;

	/* 1) Initialisation of Location Matrix datastructure */
	lmStruct->numberOfNodes = 0;
	memset( lmStruct->nodeIDList, -1, sizeof(int)*REP_MAXNODESPERPATCH );

	inc_E = IArray_New();
	inc_N = IArray_New();
  /* 2) find the elements around the node and put in elList */
	FeMesh_GetNodeElements( mesh, lPatchID, inc_E );
	*nEl = IArray_GetSize( inc_E );
	memcpy(elList, IArray_GetPtr( inc_E ), *nEl*sizeof(unsigned) );

	/* 3) foreach element find each node which will form patch
	 * and save in lmStruct */
	for( el_I = 0 ; el_I < *nEl ; el_I++ ) {
		dElID = elList[ el_I ];
		FeMesh_GetElementNodes( mesh, dElID, inc_N );
		nodesPerEl = IArray_GetSize( inc_N );
		nodeList = IArray_GetPtr( inc_N );
		for( node_I = 0 ; node_I < nodesPerEl ; node_I++ ) {
			dNodeID = nodeList[node_I];
			if( el_I == 0 ) {
				/* This is just a speed-up condition for large meshes
				 * it can be ignored altogether and you can just run the else condition,
				 * but it's a poor search algorithm */
				lmStruct->nodeIDList[ lmStruct->numberOfNodes ] = dNodeID;
				lmStruct->numberOfNodes++;
			} else {
				/* if nodeID NOT in table
					   then add nodeID to table
				*/
				if( _REP_Algorithm_locateInPatchListStruct( lmStruct, dNodeID ) == -1 ) {
					lmStruct->nodeIDList[ lmStruct->numberOfNodes ] = dNodeID;
					lmStruct->numberOfNodes++;
				}
			}
		}	 
	}
	NewClass_Delete( inc_N );
	NewClass_Delete( inc_E );
}

int _REP_Algorithm_locateInPatchListStruct( LmStruct* lmStruct, int nodeID ) {
	/*TODO: This is a linear search, could take a while for 3-D meshes, especially if non-square elements are used 
	 * should change it to quicker search */
	int node_I;
	for( node_I = 0 ; node_I < lmStruct->numberOfNodes ; node_I++ ){
		if( lmStruct->nodeIDList[ node_I ] == nodeID ) 
			return node_I;
	}
	return -1;
}

void _REP_Algorithm_ZeroHiFi( REP_Algorithm* self, double*** elHi_Mat, double*** elFi_Mat ) {
  /* Function Describption: zeros the Hi and Fi memory spaces */
	/* TODO: Can only do one repFieldList element now */
	int components           = self->repFieldList[0]->fieldComponentCount;
	int dim                  = self->repFieldList[0]->dim;
	int order = self->repFieldList[0]->orderOfInterpolation;
	int nodesPerEl           = self->repFieldList[0]->nodesPerEl;
	int comp_I, dim_I, node_I;
	
	for( comp_I = 0 ; comp_I < components ; comp_I++ ) {
		for( dim_I = 0 ; dim_I < dim; dim_I++ ) 
			memset( elHi_Mat[comp_I][dim_I], 0, order*sizeof(double) );
		for( node_I = 0 ; node_I < nodesPerEl; node_I++ ) 
			memset( elFi_Mat[comp_I][node_I], 0, dim*sizeof(double) );
	}
}
void _REP_Algorithm_AssembleAllLocalElement( REP_Algorithm* self ) {
  /* Function Describption: assembles all local element H and F objects for each
   * repField associated with this REP_Algorithm
   */
  double**** elHi_Mat;
  double*** elFi_Mat;
  double** p_pVec, **GNx;
  int field_I, lElement_I, repFieldCount, numberLocalElements;

  repFieldCount = self->repFieldCount;
	elHi_Mat = Memory_Alloc_Array( double***, repFieldCount, "localElement_HiPerField" );
	elFi_Mat = Memory_Alloc_Array( double**, repFieldCount, "localElement_FiPerField" );
	p_pVec = Memory_Alloc_Array( double*, repFieldCount, "position polynomialPerField" );
	/* Assume GNx this the same for all fields */
	GNx = Memory_Alloc_2DArray( double, self->repFieldList[0]->dim, REP_MAXNODESPERPATCH, (Name)"GNx" );
	numberLocalElements = FeMesh_GetElementLocalSize( self->repFieldList[0]->feMesh );

	for( field_I = 0 ; field_I < self->repFieldCount ; field_I++  ) {
		/* TODO assumes constant nodesPerEl across domain */
		elHi_Mat[field_I] = Memory_Alloc_3DArray( double, self->repFieldList[field_I]->fieldComponentCount, self->repFieldList[field_I]->dim * self->repFieldList[field_I]->nodesPerEl, self->repFieldList[field_I]->orderOfInterpolation, (Name)"localElement_Hi"  );
		elFi_Mat[field_I] = Memory_Alloc_2DArray( double, self->repFieldList[field_I]->fieldComponentCount, self->repFieldList[field_I]->dim * self->repFieldList[field_I]->nodesPerEl, (Name)"localElement_Fi" );
		p_pVec[field_I] = Memory_Alloc_Array( double, self->repFieldList[field_I]->orderOfInterpolation, "position polynomial" );
	}

	/* go throught all local elements */
	for( lElement_I = 0 ; lElement_I < numberLocalElements ; lElement_I++ ) {
		/* Zero temp data structures */
		for( field_I = 0 ; field_I < self->repFieldCount ; field_I++ )
			_RecoveredFeVariable_ZeroHiFi( self->repFieldList[field_I], elHi_Mat[field_I], elFi_Mat[field_I] );

		_REP_Algorithm_AssembleElement( self, lElement_I, elHi_Mat, elFi_Mat, p_pVec, GNx );

		/* incorporate Hi_Mat and Fi_Mat into proc size arrays */
		for( field_I = 0 ; field_I < self->repFieldCount ; field_I++ ) 
			self->repFieldList[field_I]->_putElIntoProc( self->repFieldList[field_I], lElement_I, elHi_Mat[field_I], elFi_Mat[field_I] );
	}

	for( field_I = 0 ; field_I < self->repFieldCount ; field_I++ ) {
		Memory_Free( elHi_Mat[field_I] );
		Memory_Free( elFi_Mat[field_I] ); 
		Memory_Free( p_pVec[field_I] );	
	}
	Memory_Free( elHi_Mat );
	Memory_Free( elFi_Mat ); 
	Memory_Free( p_pVec );
	Memory_Free( GNx );
}

void _REP_Algorithm_AssembleElement( REP_Algorithm* self, int lElement_I, double**** Hi_Mat, double*** Fi_Mat, double** p_pVec, double** GNx ) {
	IntegrationPointsSwarm* swarm = self->IPswarm;
	FeMesh*                 mesh = self->repFieldList[0]->rawField->feMesh;
	ElementType*            elementType;
	IntegrationPoint*       particle;
 	double globalCoord[3], detJac;
	int cell_I, cellParticleCount, nodesPerEl, cParticle_I;
	/* Only need one */
	int dim = self->repFieldList[0]->dim;
	int field_I;

	/* Get the element type */
	elementType = FeMesh_GetElementType( mesh, lElement_I );

	/* Get the number of nodes per element */
	FeMesh_GetElementNodes( mesh, lElement_I, self->incArray );
	nodesPerEl = IArray_GetSize( self->incArray );

	/* Get number of particles per element */
	cell_I = CellLayout_MapElementIdToCellId( swarm->cellLayout, lElement_I );
	cellParticleCount = swarm->cellParticleCountTbl[ cell_I ];

	for( cParticle_I = 0 ; cParticle_I < cellParticleCount ; cParticle_I++ ) {
		particle = (IntegrationPoint*) Swarm_ParticleInCellAt( swarm, cell_I, cParticle_I );	

		/* calculate derivatives on particles */
		ElementType_ShapeFunctionsGlobalDerivs( elementType, 
							mesh, 
							lElement_I,
							particle->xi,
							dim,
							&detJac, GNx );
		/* put particle co-ords into polynomial */
		FeMesh_CoordLocalToGlobal( mesh, lElement_I, particle->xi, globalCoord );

		for( field_I = 0 ; field_I < self->repFieldCount ; field_I++ )  {
			_RecoveredFeVariable_AssembleAtParticle( self->repFieldList[field_I], self->constitutiveMatrix, 
					self->constitutiveMatrix->currentParticleIndex, particle, lElement_I, 
					globalCoord, (double**)GNx, detJac, Hi_Mat[field_I], Fi_Mat[field_I]);

			/* could be a useful functionPtr later on, JG 12May09
			self->repFieldList[field_I]->_assembleOnParticle( self->repFieldList[field_I], self->constitutiveMatrix, 
					self->constitutiveMatrix->currentParticleIndex, particle, lElement_I, 
					globalCoord, (double**)GNx, detJac, Hi_Mat[field_I], Fi_Mat[field_I]);
					*/
		}

	}
}

void _REP_Algorithm_Solver(double **array, double* bVec, int n) {
	/* Function Description: Calls the actual lapack function to solve the small dense linear system Ax=b */
	int i,j, info, LDA, LDB, N, NRHS;
	
	int* IPIV = Memory_Alloc_Array( int, n, "IPIV" ); 
	double* AT = Memory_Alloc_Array( double, n * n , "AT" ); 
	assert( IPIV );
	assert( AT );
	/* must transform, because LAPACK uses Fortran 
	 * transform AMatrix in A Column */
	for( i = 0 ; i < n ; i++ ) {
		for( j = 0 ; j < n ; j++ ) {
			AT[ (i * n) + j ] = array[j][i];
		}
	}
		
	NRHS = 1;
	N = LDA = LDB = n;
	dgesv_( &N, &NRHS, AT, &LDA, IPIV, bVec, &LDB, &info);

	Journal_Firewall( info == 0, Journal_Register( Error_Type, (Name)"error_REP"  ), "Error: In %s looks like the lapack solver (DGESV) died with the error code %d. Could be due to ill-conditioned matrix ... I advise that you manually print the results of the matrices that lapack uses or contact a developer.\n", __func__, info);
	Memory_Free(IPIV);
	Memory_Free(AT);
}



