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
** $Id: testStiffnessMatrix-nonZeroCalculation.c 964 2007-10-11 08:03:06Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include "StgFEM/Discretisation/Discretisation.h"
//#include "StgFEM/SLE/LinearAlgebra/LinearAlgebra.h"
#include "StgFEM/SLE/SystemSetup/SystemSetup.h"

#include <stdio.h>
#include <stdlib.h>

struct _Node {
	Coord coord;
	double vel[3];
};

struct _Element {
	__FiniteElement_Element
};

struct _Particle {
	__IntegrationPoint
};


int main( int argc, char* argv[] ) {
#if 0
	MPI_Comm			CommWorld;
	int				rank;
	int				numProcessors;
	int				procToWatch;
	Dictionary*			dictionary;
	Dictionary_Entry_Value*		bcVarList;
	Dictionary_Entry_Value*		currBC_Var;
	Dictionary_Entry_Value*		currBC;
	Dictionary_Entry_Value*		bcList;
	Dictionary_Entry_Value*		compositeVC_Entry;
	Topology*			nTopology;
	ElementLayout*			eLayout;
	NodeLayout*			nLayout;
	MeshDecomp*			decomp;
	MeshLayout*			meshLayout;
	DofLayout*			dofLayout;
	ElementType_Register*		elementType_Register;
	ExtensionManager_Register*	extensionMgr_Register;
	EntryPoint_Register*		entryPoint_Register;
	FiniteElement_Mesh*		feMesh;
	Dictionary*			compositeDict;
	CompositeVC*			compositeVC;
	FieldVariable_Register*		fV_Register;
	FeVariable*			feVariable;
	Variable_Register*		variableRegister;
	ConditionFunction_Register*	condFunc_Register;
	Index				i;
	Name				velocityName = "velocity";
	Name				velCompNames[3] = { "vx", "vy", "vz" };
	Stream*				stream;
	Dimension_Index			numDims = 2;
	Dimension_Index			dim_I;
	ForceVector*			fVector;
	StiffnessMatrix*		kMatrix;
	void*				swarm;
	
	/* Initialise MPI, get world info */
	MPI_Init( &argc, &argv );
	MPI_Comm_dup( MPI_COMM_WORLD, &CommWorld );
	MPI_Comm_size( CommWorld, &numProcessors );
	MPI_Comm_rank( CommWorld, &rank );
	
	StGermain_Init( &argc, &argv );
	StgFEM_Discretisation_Init( &argc, &argv );
	//StgFEM_SLE_LinearAlgebra_Init( &argc, &argv );
	StgFEM_SLE_SystemSetup_Init( &argc, &argv );
	MPI_Barrier( CommWorld ); /* Ensures copyright info always come first in output */
	
	stream = Journal_Register (Info_Type, "myStream");
	/* set line below to true to True to activate outputting */
	Journal_Enable_TypedStream( DebugStream_Type, False );
	Stream_SetLevelBranch( StgFEM_Debug, 2 );
	Stream_EnableBranch( StgFEM_Debug, True );
	Stream_SetLevel( Journal_Register( DebugStream_Type, "StgFEM.SLE.SystemSetup.StiffnessMatrix" ), 3 );
	Stream_SetLevel( Journal_Register( DebugStream_Type, "StgFEM.Discretisation.FeEquationNumber" ), 1 );

	if( argc >= 2 ) {
		procToWatch = atoi( argv[1] );
	}
	else {
		procToWatch = 0;
	}
	if( rank == procToWatch ) printf( "Watching rank: %i\n", rank );
	
	/* Read input */
	dictionary = Dictionary_New();
	Dictionary_Add( dictionary, "rank", Dictionary_Entry_Value_FromUnsignedInt( rank ) );
	Dictionary_Add( dictionary, "numProcessors", Dictionary_Entry_Value_FromUnsignedInt( numProcessors ) );
	Dictionary_Add( dictionary, "meshSizeI", Dictionary_Entry_Value_FromUnsignedInt( 5 ) );
	Dictionary_Add( dictionary, "meshSizeJ", Dictionary_Entry_Value_FromUnsignedInt( 5 ) );
	Dictionary_Add( dictionary, "minX", Dictionary_Entry_Value_FromDouble( 0.0f ) );
	Dictionary_Add( dictionary, "minY", Dictionary_Entry_Value_FromDouble( 0.0f ) );
	Dictionary_Add( dictionary, "minZ", Dictionary_Entry_Value_FromDouble( 0.0f ) );
	Dictionary_Add( dictionary, "maxX", Dictionary_Entry_Value_FromDouble( 300.0f ) );
	Dictionary_Add( dictionary, "maxY", Dictionary_Entry_Value_FromDouble( 12.0f ) );
	Dictionary_Add( dictionary, "maxZ", Dictionary_Entry_Value_FromDouble( 300.0f ) );
	Dictionary_Add( dictionary, "allowPartitionOnElement", Dictionary_Entry_Value_FromBool( False ) );
	Dictionary_Add( dictionary, "buildElementNodeTbl", Dictionary_Entry_Value_FromBool( True ) );
	Dictionary_Add( dictionary, "shadowDepth", Dictionary_Entry_Value_FromUnsignedInt( 1 ) );
	
	bcList = Dictionary_Entry_Value_NewList();
	
	bcVarList = Dictionary_Entry_Value_NewList();
	currBC_Var = Dictionary_Entry_Value_NewStruct();
	Dictionary_Entry_Value_AddMember( currBC_Var, "name", Dictionary_Entry_Value_FromString( "vx" ) );
	Dictionary_Entry_Value_AddMember( currBC_Var, "type", Dictionary_Entry_Value_FromString( "double" ) );
	Dictionary_Entry_Value_AddMember( currBC_Var, "value", Dictionary_Entry_Value_FromDouble( 0.0f ) );
	Dictionary_Entry_Value_AddElement( bcVarList, currBC_Var );
	currBC = Dictionary_Entry_Value_NewStruct();
	Dictionary_Entry_Value_AddMember( currBC, "type", Dictionary_Entry_Value_FromString( "WallVC" ) );
	Dictionary_Entry_Value_AddMember( currBC, "name", Dictionary_Entry_Value_FromString( "vxWallVC" ) );
	Dictionary_Entry_Value_AddMember( currBC, "wall", Dictionary_Entry_Value_FromString( "left" ) );
	Dictionary_Entry_Value_AddMember( currBC, "variables", bcVarList );
	Dictionary_Entry_Value_AddElement( bcList, currBC );

	bcVarList = Dictionary_Entry_Value_NewList();
	currBC_Var = Dictionary_Entry_Value_NewStruct();
	Dictionary_Entry_Value_AddMember( currBC_Var, "name", Dictionary_Entry_Value_FromString( "vy" ) );
	Dictionary_Entry_Value_AddMember( currBC_Var, "type", Dictionary_Entry_Value_FromString( "double" ) );
	Dictionary_Entry_Value_AddMember( currBC_Var, "value", Dictionary_Entry_Value_FromDouble( 0.0f ) );
	Dictionary_Entry_Value_AddElement( bcVarList, currBC_Var );
	currBC = Dictionary_Entry_Value_NewStruct();
	Dictionary_Entry_Value_AddMember( currBC, "type", Dictionary_Entry_Value_FromString( "WallVC" ) );
	Dictionary_Entry_Value_AddMember( currBC, "name", Dictionary_Entry_Value_FromString( "vyWallVC" ) );
	Dictionary_Entry_Value_AddMember( currBC, "wall", Dictionary_Entry_Value_FromString( "bottom" ) );
	Dictionary_Entry_Value_AddMember( currBC, "variables", bcVarList );
	Dictionary_Entry_Value_AddElement( bcList, currBC );

	compositeVC_Entry = Dictionary_Entry_Value_NewStruct();
	Dictionary_Entry_Value_AddMember( compositeVC_Entry, "type", Dictionary_Entry_Value_FromString( "CompositeVC" ) );
	Dictionary_Entry_Value_AddMember( compositeVC_Entry, "vcList", bcList );
	Dictionary_Add( dictionary, "bcs", compositeVC_Entry );
	
	/* create the layout, dof and mesh to use */
	nTopology = (Topology*)IJK6Topology_New( "IJK6Topology", dictionary );
	eLayout = (ElementLayout*)ParallelPipedHexaEL_New( "ElementLaoyut", 2, dictionary );
	nLayout = (NodeLayout*)CornerNL_New( "CornerNL", dictionary, eLayout, nTopology );
	decomp = (MeshDecomp*)HexaMD_New( "HexaMD", dictionary, MPI_COMM_WORLD, eLayout, nLayout );
	meshLayout = MeshLayout_New( "MeshLayout", eLayout, nLayout, decomp );
	
	extensionMgr_Register = ExtensionManager_Register_New();
	elementType_Register = ElementType_Register_New("elementTypeRegister");
	ElementType_Register_Add( elementType_Register, (ElementType*)ConstantElementType_New("constant") );
	ElementType_Register_Add( elementType_Register, (ElementType*)BilinearElementType_New("bilinear") );
	ElementType_Register_Add( elementType_Register, (ElementType*)TrilinearElementType_New("trilinear") );
	entryPoint_Register = EntryPoint_Register_New();
	condFunc_Register = ConditionFunction_Register_New();
	feMesh = FiniteElement_Mesh_New( "testMesh", meshLayout, sizeof(Node), sizeof(Element), extensionMgr_Register,
		elementType_Register, dictionary );
	Build( feMesh, 0, False );
	
	/* Create variable register */
	variableRegister = Variable_Register_New();
	
	/* Create variables */
	Variable_NewVector( 
		"velocity", 
		Variable_DataType_Double, 
		3, 
		&feMesh->nodeDomainCount, 
		(void**)&feMesh->nodeCoord, 
		variableRegister, 
		"vx", 
		"vy", 
		"vz" );

	dofLayout = DofLayout_New( "dofLayout", variableRegister, decomp->nodeDomainCount );
	for (i = 0; i < decomp->nodeDomainCount; i++)
	{
		for ( dim_I = 0; dim_I < numDims; dim_I++ ) {
			DofLayout_AddDof_ByVarName(dofLayout, velCompNames[dim_I], i);
		}
	}
	
	compositeDict = Dictionary_Entry_Value_AsDictionary( Dictionary_Get( dictionary, "bcs" ) );
	compositeVC = CompositeVC_New( "CompositeVC", variableRegister, condFunc_Register, compositeDict, feMesh );

	/* Create the fe variable */
	fV_Register = FieldVariable_Register_New();
	feVariable = FeVariable_New( velocityName, feMesh, NULL, dofLayout, compositeVC, NULL, NULL, numDims,False,
		StgFEM_Native_ImportExportType, StgFEM_Native_ImportExportType, fV_Register );
	
	/* Create a hack swarm ptr - we won't actually use it in this test */
	swarm = feVariable;
	
	/* Create a Stiffness Matrix to use for testing */
	fVector = ForceVector_New( "F", feVariable, numDims, entryPoint_Register, CommWorld );	
	kMatrix = StiffnessMatrix_New( "K", feVariable, feVariable, fVector, NULL, numDims, False, True,
		entryPoint_Register, CommWorld );
	
	/* Build and initialise system */
	Build( feVariable, 0, False );
	Variable_Register_BuildAll(variableRegister);

	/* Ok, so now test the stiffness matrix non-zero count functions. */
	/* The left wall has had BCs applied to it, so this will affect the counts */
	kMatrix->rowLocalSize = feVariable->eqNum->localEqNumsOwnedCount;
	_StiffnessMatrix_CalculateNonZeroEntries( kMatrix );

	if ( rank == procToWatch ) {
		Index	matRow_I = 0;
		Node_LocalIndex 	lNode_I = 0;
		Node_GlobalIndex 	gNode_I = 0;
		Dof_Index		currNodeNumDofs;
		Dof_Index		nodeLocalDof_I;
		Dof_EquationNumber	currEqNum;
		Variable*		currVariable;
		
		Journal_Printf( stream, "Printing non-zero entries calculation:\n" );
		for( lNode_I=0; lNode_I < feMesh->nodeLocalCount; lNode_I++ ) {
			gNode_I = feMesh->nodeL2G[lNode_I];
			currNodeNumDofs = dofLayout->dofCounts[ lNode_I ];
			Journal_Printf( stream, "node %d (global index %d):\n", lNode_I, gNode_I, currNodeNumDofs );
			
			/* Print each dof */
			for ( nodeLocalDof_I = 0; nodeLocalDof_I < currNodeNumDofs; nodeLocalDof_I++ ) {
				currVariable = DofLayout_GetVariable( dofLayout, lNode_I, nodeLocalDof_I );
				Journal_Printf( stream, "\tdof %d \"%s\": ", nodeLocalDof_I, currVariable->name );
				currEqNum = feVariable->eqNum->destinationArray[lNode_I][nodeLocalDof_I];
				if ( currEqNum == -1 ) {
					Journal_Printf( stream, "(BC)", currEqNum );
				}
				else {
					Journal_Printf( stream, "(eq num %d)", currEqNum );
					matRow_I = currEqNum - feVariable->eqNum->firstOwnedEqNum;
					if ( matRow_I >= feVariable->eqNum->localEqNumsOwnedCount ) {
						Journal_Printf( stream, "\tnot owned by this proc." );
					}
					else {
						Journal_Printf( stream, "\tlocal matRow %d: diag entries = %d, off diag = %d",
							matRow_I,
							kMatrix->diagonalNonZeroIndices[matRow_I],
							kMatrix->offDiagonalNonZeroIndices[matRow_I] );
					}		
				}
				Journal_Printf( stream, "\n", currEqNum );
			}
		}
	}

	/* Destroy stuff */
	Stg_Class_Delete( kMatrix );
	Stg_Class_Delete( fVector );
	Stg_Class_Delete( fV_Register );
	Stg_Class_Delete( compositeVC );
	Stg_Class_Delete( feMesh );
	Stg_Class_Delete( entryPoint_Register );
	Stg_Class_Delete( elementType_Register );
	Stg_Class_Delete( extensionMgr_Register );
	Stg_Class_Delete( dofLayout );
	Stg_Class_Delete( meshLayout );
	Stg_Class_Delete( decomp );
	Stg_Class_Delete( nLayout );
	Stg_Class_Delete( eLayout );
	Stg_Class_Delete( nTopology );
	Stg_Class_Delete( dictionary );
	
	StgFEM_SLE_SystemSetup_Finalise();
	//StgFEM_SLE_LinearAlgebra_Finalise();
	StgFEM_Discretisation_Finalise();
	StGermain_Finalise();
	
	/* Close off MPI */
	MPI_Finalize();
	
	return 0; /* success */
#endif
}
