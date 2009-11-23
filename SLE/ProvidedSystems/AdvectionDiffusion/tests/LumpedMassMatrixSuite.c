/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street, Melbourne, 3053, Australia.
**
** Authors:
**   Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**   Patrick D. Sunter, Software Engineer, VPAC. (pds@vpac.org)
**   Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**   Siew-Ching Tan, Software Engineer, VPAC. (siew@vpac.org)
**   Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
**   Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
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
** Role:
**   Tests the LumpedMassMatrixSuite
**
** $Id: testTemplate.c 3462 2006-02-19 06:53:24Z WalterLandry $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>

#include "pcu/pcu.h"
#include <StGermain/StGermain.h> 
#include <StgDomain/StgDomain.h>
#include "StgFEM/Discretisation/Discretisation.h"
#include "StgFEM/SLE/SystemSetup/SystemSetup.h"
#include "StgFEM/SLE/ProvidedSystems/AdvectionDiffusion/AdvectionDiffusion.h"

#include "LumpedMassMatrixSuite.h"

#define CURR_MODULE_NAME "StgFEMContext.c"
	
typedef struct Node Node;
typedef struct Element Element;

struct Node {
	double phi;
};

struct Element {
	__FiniteElement_Element
};

struct _Particle {
	Coord coord;
};


typedef struct {
	MPI_Comm	comm;
	unsigned rank;
	unsigned nProcs;
} LumpedMassMatrixSuiteData;

void LumpedMassMatrixSuite_quadratic(Index index, Variable_Index var_I, void* context, void* result) {
	*(double *)result = 20.0;
}

FeMesh* LumpedMassMatrixSuite_buildFeMesh( unsigned nDims, unsigned* size, double* minCrds, double* maxCrds, ExtensionManager_Register* emReg, ElementType_Register* etReg ) {
	CartesianGenerator*	gen;
	FeMesh*					feMesh;
	unsigned					maxDecomp[3] = {0, 1, 1};

	gen = CartesianGenerator_New( "", NULL );
	gen->shadowDepth = 0;
	CartesianGenerator_SetDimSize( gen, nDims );
	CartesianGenerator_SetTopologyParams( gen, size, 0, NULL, maxDecomp );
	CartesianGenerator_SetGeometryParams( gen, minCrds, maxCrds );

	feMesh = FeMesh_New( "", NULL );
	Mesh_SetExtensionManagerRegister( feMesh, emReg );
	Mesh_SetGenerator( feMesh, gen );
	FeMesh_SetElementFamily( feMesh, "linear" );

	Mesh_SetTopologyDataSize( feMesh, MT_VERTEX, sizeof(Node) );
	Mesh_SetTopologyDataSize( feMesh, nDims, sizeof(Element) );

	Stg_Component_Build( feMesh, NULL, False );
	Stg_Component_Initialise( feMesh, NULL, False );

	return feMesh;
}

void LumpedMassMatrixSuite_Setup( LumpedMassMatrixSuiteData* data ) {
	/* MPI Initializations */
	data->comm = MPI_COMM_WORLD;
	MPI_Comm_rank( data->comm, &data->rank );
	MPI_Comm_size( data->comm, &data->nProcs );
}

void LumpedMassMatrixSuite_Teardown( LumpedMassMatrixSuiteData* data ) {
}

void LumpedMassMatrixSuite_TestLumpedMassMatrix( LumpedMassMatrixSuiteData* data ) {
	Dictionary*						dictionary;
	Dictionary_Entry_Value*		currBC;
	Dictionary_Entry_Value*		bcList;
	unsigned							nDims = 2;
	unsigned							meshSize[3] = {2, 2, 0};
	unsigned							nDomainVerts;
	double							minCrds[3] = {0.0, 0.0, 0.0};
	double							maxCrds[3] = {1.2, 1.2, 1.2};
	FeMesh*							feMesh;
	Node*								nodes;
	DofLayout*						dofs;
	ElementType_Register*		elementType_Register;
	ExtensionManager_Register*	extensionMgr_Register;
	WallVC*							wallVC;
	FieldVariable_Register*		fV_Register;
	FeVariable*						feVariable;
	Variable_Register*			variableRegister;
	FiniteElementContext*		context;
	Index								i;
	Dimension_Index				dim;
	Stream*							stream;
	/* Swarm Stuff */
	CellLayout*						singleCellLayout;
	ParticleLayout*				gaussParticleLayout;
	Swarm*							swarm;
	unsigned							dimExists[] = { True, True, False };
	/* Mass Matrix Stuff */
	ForceVector*					massMatrix;
	Vec								copyMatrix;
	Vec								expectedMatrix;
	PetscViewer						viewer;
	PetscTruth						flg;
	LumpedMassMatrixForceTerm*	massMatrixForceTerm;
	Particle_InCellIndex			particlesPerDim[] = {2,2,2};
	char								expected_file[PCU_PATH_MAX];

	stream = Journal_Register (Info_Type, "LumpedMassMatrixStream");

	/* Read input */
	dictionary = Dictionary_New();
	Dictionary_Add( dictionary, "outputPath", Dictionary_Entry_Value_FromString( "./output" ) );
	Dictionary_Add( dictionary, "rank", Dictionary_Entry_Value_FromUnsignedInt( data->rank ) );
	Dictionary_Add( dictionary, "numProcessors", Dictionary_Entry_Value_FromUnsignedInt( data->nProcs ) );
	Dictionary_Add( dictionary, "gaussParticlesX", Dictionary_Entry_Value_FromUnsignedInt( 2 ) );
	Dictionary_Add( dictionary, "gaussParticlesY", Dictionary_Entry_Value_FromUnsignedInt( 2 ) );

	bcList = Dictionary_Entry_Value_NewList();
	currBC = Dictionary_Entry_Value_NewStruct();
	Dictionary_Entry_Value_AddMember( currBC, "name", Dictionary_Entry_Value_FromString( "phi" ) );
	Dictionary_Entry_Value_AddMember( currBC, "type", Dictionary_Entry_Value_FromString( "double" ) );
	Dictionary_Entry_Value_AddMember( currBC, "value", Dictionary_Entry_Value_FromDouble( -1.0f ) );
	Dictionary_Entry_Value_AddElement( bcList, currBC );
	currBC = Dictionary_Entry_Value_NewStruct();
	Dictionary_Entry_Value_AddMember( currBC, "wall", Dictionary_Entry_Value_FromString( "left" ) );
	Dictionary_Entry_Value_AddMember( currBC, "variables", bcList );
	Dictionary_Add( dictionary, "boundaryCondition", currBC );

	/* Create Context */
	context = FiniteElementContext_New( "context", 0,0, data->comm, dictionary );

	dim = context->dim = Dictionary_GetUnsignedInt_WithDefault( dictionary, "dim", 2 );
	Journal_Enable_TypedStream( DebugStream_Type, True );
	Stream_SetLevelBranch( StgFEM_Debug, 3 );
	Stream_EnableBranch( StgFEM_Debug, True );
	
	/* create the layout, dof and mesh to use */
	extensionMgr_Register = ExtensionManager_Register_New();
	elementType_Register = ElementType_Register_New("elementTypeRegister");
	ElementType_Register_Add( elementType_Register, (ElementType*)ConstantElementType_New("constant") );
	ElementType_Register_Add( elementType_Register, (ElementType*)BilinearElementType_New("bilinear") );
	ElementType_Register_Add( elementType_Register, (ElementType*)TrilinearElementType_New("trilinear") );
	
	feMesh = (FeMesh*) LumpedMassMatrixSuite_buildFeMesh( nDims, meshSize, minCrds, maxCrds, extensionMgr_Register, elementType_Register );
	nDomainVerts = Mesh_GetDomainSize( feMesh, MT_VERTEX );
	nodes = Mesh_GetTopologyData( feMesh, MT_VERTEX );
	
	/* Create variable register */
	variableRegister = Variable_Register_New();
	
	/* Create variables */
	Variable_NewScalar( "phi", Variable_DataType_Double, &nDomainVerts, NULL, (void**)&nodes, variableRegister );

	dofs = DofLayout_New( "dofLayout", context, variableRegister, Mesh_GetDomainSize( feMesh, MT_VERTEX ), NULL );
	for (i = 0; i < Mesh_GetDomainSize( feMesh, MT_VERTEX ); i++)
		DofLayout_AddDof_ByVarName(dofs, "phi", i);
	
	wallVC = WallVC_New( "WallVC", "boundaryCondition", variableRegister, NULL, dictionary, feMesh );

	/* Create the finite element field variable*/
	fV_Register = FieldVariable_Register_New();
	feVariable = FeVariable_New( "phi", NULL, feMesh, NULL, dofs, wallVC, NULL, NULL, context->dim, False, False, False, fV_Register );
	
	/* Create Swarm */
	if ( 3 == dim ) 
		dimExists[K_AXIS] = True;
	singleCellLayout= (CellLayout*)SingleCellLayout_New( "SingleCellLayout", (AbstractContext*)context, dimExists, NULL, NULL );
	gaussParticleLayout = (ParticleLayout*)GaussParticleLayout_New( "GaussParticleLayout", NULL, LocalCoordSystem, True, dim, particlesPerDim );
	swarm = Swarm_New( "gaussSwarm", (AbstractContext*)context, singleCellLayout, gaussParticleLayout,
		dim, sizeof(IntegrationPoint), extensionMgr_Register, context->variable_Register, MPI_COMM_WORLD, NULL );
	
	/* Lumping of Mass Matrix */
	massMatrix = ForceVector_New( "MassMatrix", context, feVariable, dim, context->entryPoint_Register, MPI_COMM_WORLD );
	massMatrixForceTerm = LumpedMassMatrixForceTerm_New( "forceTerm", context, massMatrix, swarm );
	EP_ReplaceAll( massMatrix->assembleForceVector, ForceVector_GlobalAssembly_General );
	ForceTerm_SetAssembleElementFunction( massMatrixForceTerm, _LumpedMassMatrixForceTerm_AssembleElement_General );

	/* Build */
	Stg_Component_Build( feMesh, context, False );
	Variable_Register_BuildAll( variableRegister );
	Stg_Component_Build( wallVC, context, False );
	Stg_Component_Build( dofs, context, False);
	Stg_Component_Build( feVariable, context, False );
	FeEquationNumber_BuildLocationMatrix( feVariable->eqNum );
	Stg_Component_Build( singleCellLayout,    context, False );
	Stg_Component_Build( gaussParticleLayout, context, False );
	Stg_Component_Build( swarm,               context, False );
	Stg_Component_Build( massMatrix,          context, False );

	/* Initialise */
	Stg_Component_Initialise( feMesh,     context, False );
	Stg_Component_Initialise( feVariable, context, False );
	Stg_Component_Initialise( swarm,      context, False );
	Stg_Component_Initialise( massMatrix, context, False );

	/* Assemble */
	if( data->rank == 0 ) {
		ForceVector_Assemble( massMatrix );

		PetscViewerCreate( MPI_COMM_WORLD, &viewer );
		PetscViewerSetType( viewer, PETSC_VIEWER_BINARY );

		pcu_filename_expected( "testLumpedMassMatrix.expected", expected_file );
		PetscViewerBinaryOpen( MPI_COMM_WORLD, expected_file, FILE_MODE_READ, &viewer );
	
		VecLoad( viewer, VECMPI, &expectedMatrix );	

		/* Try out optimised one */
		VecSet( massMatrix->vector, 0.0 );
	
		/* Assemble */
		ForceTerm_SetAssembleElementFunction( massMatrixForceTerm, _LumpedMassMatrixForceTerm_AssembleElement_Box );
		ForceVector_Assemble( massMatrix );

		/* Check vector */
		VecEqual( (Vec)massMatrix->vector, expectedMatrix, &flg );
		pcu_check_true( flg );
	}

	/* Destroy stuff */
	Stg_Class_Delete( massMatrix );
	Stg_Class_Delete( massMatrixForceTerm );
	Stg_Class_Delete( fV_Register );
	Stg_Class_Delete( wallVC );
	Stg_Class_Delete( feMesh );
	Stg_Class_Delete( elementType_Register );
	Stg_Class_Delete( extensionMgr_Register );
	Stg_Class_Delete( dofs );
	Stg_Class_Delete( dictionary );
}

void LumpedMassMatrixSuite( pcu_suite_t* suite ) {
	pcu_suite_setData( suite, LumpedMassMatrixSuiteData );
   pcu_suite_setFixtures( suite, LumpedMassMatrixSuite_Setup, LumpedMassMatrixSuite_Teardown );
   pcu_suite_addTest( suite, LumpedMassMatrixSuite_TestLumpedMassMatrix );
}
