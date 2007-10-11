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
** $Id: testLumpedMassMatrix.c 964 2007-10-11 08:03:06Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include "StgFEM/Discretisation/Discretisation.h"
#include "StgFEM/SLE/LinearAlgebra/LinearAlgebra.h"
#include "StgFEM/SLE/SystemSetup/SystemSetup.h"
#include "StgFEM/SLE/ProvidedSystems/AdvectionDiffusion/AdvectionDiffusion.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "petsc.h"


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


FeMesh* buildFeMesh( unsigned nDims, unsigned* size, 
		     double* minCrds, double* maxCrds, 
		     ExtensionManager_Register* emReg, 
		     ElementType_Register* etReg )
{
	CartesianGenerator*	gen;
	FeMesh*			feMesh;
	unsigned		maxDecomp[3] = {0, 1, 1};

	gen = CartesianGenerator_New( "" );
	gen->shadowDepth = 0;
	CartesianGenerator_SetDimSize( gen, nDims );
	CartesianGenerator_SetTopologyParams( gen, size, 0, NULL, maxDecomp );
	CartesianGenerator_SetGeometryParams( gen, minCrds, maxCrds );

	feMesh = FeMesh_New( "" );
	Mesh_SetExtensionManagerRegister( feMesh, emReg );
	Mesh_SetGenerator( feMesh, gen );
	FeMesh_SetElementFamily( feMesh, "linear" );

	Mesh_SetTopologyDataSize( feMesh, MT_VERTEX, sizeof(Node) );
	Mesh_SetTopologyDataSize( feMesh, nDims, sizeof(Element) );

	Stg_Component_Build( feMesh, NULL, False );
	Stg_Component_Initialise( feMesh, NULL, False );

	return feMesh;
}


int main( int argc, char* argv[] ) {
	MPI_Comm                   CommWorld;
	int                        rank;
	int                        numProcessors;
	int                        procToWatch;
	Dictionary*                dictionary;
	Dictionary_Entry_Value*    currBC;
	Dictionary_Entry_Value*	   bcList;

	unsigned	nDims = 2;
	unsigned	meshSize[3] = {2, 2, 0};
	double		minCrds[3] = {0.0, 0.0, 0.0};
	double		maxCrds[3] = {1.2, 1.2, 1.2};
	FeMesh*		feMesh;
	unsigned	nDomainVerts;
	Node*		nodes;

	DofLayout*                 dofs;
	ElementType_Register*      elementType_Register;
	ExtensionManager_Register* extensionMgr_Register;
	WallVC*                    wallVC;
	FieldVariable_Register*    fV_Register;
	FeVariable*                feVariable;
	Variable_Register*         variableRegister;
	FiniteElementContext*      context;
	Index                      i;
	Dimension_Index            dim;
	Stream*	                   stream;
	/* Swarm Stuff */
	CellLayout*                singleCellLayout;
	ParticleLayout*            gaussParticleLayout;
	Swarm*                     swarm;
	unsigned int               dimExists[]     = { True, True, False };
	/* Mass Matrix Stuff */
	ForceVector*               massMatrix;
	LumpedMassMatrixForceTerm* massMatrixForceTerm;
	/* Stream Stuff */
	Stream*                    outputStream;
	Particle_InCellIndex       particlesPerDim[] = {2,2,2};
	
	/* Initialise MPI, get world info */
	MPI_Init( &argc, &argv );
	MPI_Comm_dup( MPI_COMM_WORLD, &CommWorld );
	MPI_Comm_size( CommWorld, &numProcessors );
	MPI_Comm_rank( CommWorld, &rank );
	
	StGermain_Init( &argc, &argv );
	StgDomain_Init( &argc, &argv );
	StgFEM_Discretisation_Init( &argc, &argv );
	StgFEM_SLE_LinearAlgebra_Init( &argc, &argv );
	StgFEM_SLE_SystemSetup_Init( &argc, &argv );
	StgFEM_SLE_ProvidedSystems_AdvectionDiffusion_Init( &argc, &argv );
	MPI_Barrier( CommWorld ); /* Ensures copyright info always come first in output */
	
	stream = Journal_Register (Info_Type, "myStream");

	if( argc >= 2 ) {
		procToWatch = atoi( argv[1] );
	}
	else {
		procToWatch = 0;
	}
	if( rank == procToWatch ) printf( "Watching rank: %i\n", rank );
	
	/* Read input */
	dictionary = Dictionary_New();
	Dictionary_Add( dictionary, "outputPath", Dictionary_Entry_Value_FromString( "./output" ) );
	Dictionary_Add( dictionary, "rank", Dictionary_Entry_Value_FromUnsignedInt( rank ) );
	Dictionary_Add( dictionary, "numProcessors", Dictionary_Entry_Value_FromUnsignedInt( numProcessors ) );
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
	context = FiniteElementContext_New( "context", 0,0, MPI_COMM_WORLD, dictionary );
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

	feMesh = buildFeMesh( nDims, meshSize, minCrds, maxCrds, 
			      extensionMgr_Register, elementType_Register );
	nDomainVerts = Mesh_GetDomainSize( feMesh, MT_VERTEX );
	nodes = Mesh_GetTopologyData( feMesh, MT_VERTEX );
	
	/* Create variable register */
	variableRegister = Variable_Register_New();
	
	/* Create variables */
	Variable_NewScalar( 
		"phi", 
		Variable_DataType_Double, 
		&nDomainVerts, 
		NULL,
		(void**)&nodes, 
		variableRegister );

	dofs = DofLayout_New( "dofLayout", variableRegister, Mesh_GetDomainSize( feMesh, MT_VERTEX ), NULL );
	for (i = 0; i < Mesh_GetDomainSize( feMesh, MT_VERTEX ); i++)
		DofLayout_AddDof_ByVarName(dofs, "phi", i);
	
	wallVC = WallVC_New( "WallVC", "boundaryCondition", variableRegister, NULL, dictionary, feMesh );

	/* Create the finite element field variable*/
	fV_Register = FieldVariable_Register_New();
	feVariable = FeVariable_New( "phi", feMesh, NULL, dofs, wallVC, NULL, NULL, context->dim,
		False, StgFEM_Native_ImportExportType, StgFEM_Native_ImportExportType, NULL, NULL,
		False, False, fV_Register );
	
	/* Create Stream */
	outputStream = Journal_Register( InfoStream_Type, CURR_MODULE_NAME );
	Stream_RedirectFile( outputStream, "output/output.dat" );

	/* Create Swarm */
	if ( 3 == dim ) 
		dimExists[K_AXIS] = True;
	singleCellLayout    = (CellLayout*)SingleCellLayout_New( "SingleCellLayout", dimExists, NULL, NULL );
	gaussParticleLayout = (ParticleLayout*)GaussParticleLayout_New( "GaussParticleLayout", dim, particlesPerDim );
	swarm = Swarm_New( "gaussSwarm", singleCellLayout, gaussParticleLayout, dim,
		sizeof(IntegrationPoint), context->extensionMgr_Register, context->variable_Register, MPI_COMM_WORLD, NULL );
	
	/* Lumping of Mass Matrix */
	massMatrix = ForceVector_New( "MassMatrix", feVariable, dim, context->entryPoint_Register, MPI_COMM_WORLD );
	massMatrixForceTerm = LumpedMassMatrixForceTerm_New( "forceTerm", massMatrix, swarm );
	EP_ReplaceAll( massMatrix->assembleForceVector, ForceVector_GlobalAssembly_General );
	ForceTerm_SetAssembleElementFunction( massMatrixForceTerm, _LumpedMassMatrixForceTerm_AssembleElement_General );

	/* Build */
	Stg_Component_Build( feMesh, context, False );
	Variable_Register_BuildAll(variableRegister);
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
	ForceVector_Assemble( massMatrix );

	/* Print out vector */
	Vector_View( massMatrix->vector, outputStream );

	/* Try out optimised one */
	Vector_Zero( massMatrix->vector );

	/* Assemble */
	ForceTerm_SetAssembleElementFunction( massMatrixForceTerm, _LumpedMassMatrixForceTerm_AssembleElement_Box );
	ForceVector_Assemble( massMatrix );

	/* Print out vector */
	Vector_View( massMatrix->vector, outputStream );	

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
	
	StgFEM_Discretisation_Finalise();
	StgDomain_Finalise();
	StGermain_Finalise();
	
	/* Close off MPI */
	MPI_Finalize();
	
	return 0; /* success */
}
