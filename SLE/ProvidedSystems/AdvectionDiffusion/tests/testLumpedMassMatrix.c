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
** $Id: testLumpedMassMatrix.c 732 2007-02-07 00:38:34Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include "StgFEM/Discretisation/Discretisation.h"
#include "StgFEM/SLE/LinearAlgebra/LinearAlgebra.h"
#include "StgFEM/SLE/SystemSetup/SystemSetup.h"
#include "StgFEM/SLE/ProvidedSystems/AdvectionDiffusion/AdvectionDiffusion.h"

#include <stdio.h>
#include <stdlib.h>
#include "petsc.h"

struct _Node {
	double phi;
};

struct _Element {
	__FiniteElement_Element
};

struct _Particle {
	Coord coord;
};

int main( int argc, char* argv[] ) {
	MPI_Comm                   CommWorld;
	int                        rank;
	int                        numProcessors;
	int                        procToWatch;
	Dictionary*                dictionary;
	Dictionary_Entry_Value*    currBC;
	Dictionary_Entry_Value*	   bcList;
	Topology*                  nTopology;
	ElementLayout*             eLayout;
	NodeLayout*                nLayout;
	MeshDecomp*                decomp;
	MeshLayout*	               meshLayout;
	DofLayout*                 dofs;
	ElementType_Register*      elementType_Register;
	ExtensionManager_Register* extensionMgr_Register;
	FiniteElement_Mesh*        feMesh;
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
	Dictionary_Add( dictionary, "meshSizeI", Dictionary_Entry_Value_FromUnsignedInt( 3 ) );
	Dictionary_Add( dictionary, "meshSizeJ", Dictionary_Entry_Value_FromUnsignedInt( 3 ) );
	Dictionary_Add( dictionary, "meshSizeK", Dictionary_Entry_Value_FromUnsignedInt( 1 ) );
	Dictionary_Add( dictionary, "minX", Dictionary_Entry_Value_FromDouble( 0.0f ) );
	Dictionary_Add( dictionary, "minY", Dictionary_Entry_Value_FromDouble( 0.0f ) );
	Dictionary_Add( dictionary, "minZ", Dictionary_Entry_Value_FromDouble( 0.0f ) );
	Dictionary_Add( dictionary, "maxX", Dictionary_Entry_Value_FromDouble( 1.2f ) );
	Dictionary_Add( dictionary, "maxY", Dictionary_Entry_Value_FromDouble( 1.2f ) );
	Dictionary_Add( dictionary, "maxZ", Dictionary_Entry_Value_FromDouble( 1.2f ) );
	Dictionary_Add( dictionary, "allowPartitionOnElement", Dictionary_Entry_Value_FromBool( False ) );
	Dictionary_Add( dictionary, "buildElementNodeTbl", Dictionary_Entry_Value_FromBool( True ) );
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
	nTopology = (Topology*)IJK6Topology_New( "IJK6Topology", dictionary );
	eLayout = (ElementLayout*)ParallelPipedHexaEL_New( "ElementLayout", dim, dictionary );
	nLayout = (NodeLayout*)CornerNL_New( "CornerNL", dictionary, eLayout, nTopology );
	decomp = (MeshDecomp*)HexaMD_New( "HexaMD", dictionary, MPI_COMM_WORLD, eLayout, nLayout );
	meshLayout = MeshLayout_New( "MeshLayout", eLayout, nLayout, decomp );
	
	extensionMgr_Register = ExtensionManager_Register_New();
	elementType_Register = ElementType_Register_New("elementTypeRegister");
	ElementType_Register_Add( elementType_Register, (ElementType*)ConstantElementType_New("constant") );
	ElementType_Register_Add( elementType_Register, (ElementType*)BilinearElementType_New("bilinear") );
	ElementType_Register_Add( elementType_Register, (ElementType*)TrilinearElementType_New("trilinear") );
	feMesh = FiniteElement_Mesh_New( "testMesh", meshLayout, sizeof(Node), sizeof(Element), extensionMgr_Register,
		elementType_Register, dictionary );
	
	/* Create variable register */
	variableRegister = Variable_Register_New();
	
	/* Create variables */
	Variable_NewScalar( 
		"phi", 
		Variable_DataType_Double, 
		&feMesh->nodeDomainCount, 
		(void**)&feMesh->node, 
		variableRegister );

	dofs = DofLayout_New( "dofLayout", variableRegister, decomp->nodeLocalCount );
	for (i = 0; i < decomp->nodeLocalCount; i++)
		DofLayout_AddDof_ByVarName(dofs, "phi", i);
	
	wallVC = WallVC_New( "WallVC", "boundaryCondition", variableRegister, NULL, dictionary, feMesh );

	/* Create the finite element field variable*/
	fV_Register = FieldVariable_Register_New();
	feVariable = FeVariable_New( "phi", feMesh, NULL, dofs, wallVC, NULL, NULL, context->dim,
		False, StgFEM_Native_ImportExportType, StgFEM_Native_ImportExportType, fV_Register );
	
	/* Create Stream */
	outputStream = Journal_Register( InfoStream_Type, CURR_MODULE_NAME );
	Stream_RedirectFile( outputStream, "output/output.dat" );

	/* Create Swarm */
	if ( 3 == dim ) 
		dimExists[K_AXIS] = True;
	singleCellLayout    = (CellLayout*)    SingleCellLayout_New( "SingleCellLayout", dimExists, NULL, NULL );
	gaussParticleLayout = (ParticleLayout*)GaussParticleLayout_New( "GaussParticleLayout", dim, particlesPerDim );
	swarm = Swarm_New( "gaussSwarm", singleCellLayout, gaussParticleLayout, dim,
		sizeof(IntegrationPoint), context->extensionMgr_Register, context->variable_Register, MPI_COMM_WORLD );
	
	/* Lumping of Mass Matrix */
	massMatrix = ForceVector_New( "MassMatrix", feVariable, dim, context->entryPoint_Register, MPI_COMM_WORLD );
	massMatrixForceTerm = LumpedMassMatrixForceTerm_New( "forceTerm", massMatrix, swarm );
	EP_ReplaceAll( massMatrix->assembleForceVector, ForceVector_GlobalAssembly_General );
	ForceTerm_SetAssembleElementFunction( massMatrixForceTerm, _LumpedMassMatrixForceTerm_AssembleElement_General );

	/* Build */
	Build( feMesh, context, False );
	Variable_Register_BuildAll(variableRegister);
	Build( wallVC, context, False );
	Build( dofs, context, False);
	Build( feVariable, context, False );
	FeEquationNumber_BuildLocationMatrix( feVariable->eqNum );
	Build( singleCellLayout,    context, False );
	Build( gaussParticleLayout, context, False );
	Build( swarm,               context, False );
	Build( massMatrix,          context, False );

	/* Initialise */
	Initialise( feMesh,     context, False );
	Initialise( feVariable, context, False );
	Initialise( swarm,      context, False );
	Initialise( massMatrix, context, False );

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
	Stg_Class_Delete( meshLayout );
	Stg_Class_Delete( decomp );
	Stg_Class_Delete( nLayout );
	Stg_Class_Delete( eLayout );
	Stg_Class_Delete( nTopology );
	Stg_Class_Delete( dictionary );
	
	StgFEM_Discretisation_Finalise();
	StGermain_Finalise();
	
	/* Close off MPI */
	MPI_Finalize();
	
	return 0; /* success */
}
