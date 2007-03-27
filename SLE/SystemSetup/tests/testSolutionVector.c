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
** $Id: testSolutionVector.c 656 2006-10-18 06:45:50Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include "StgFEM/Discretisation/Discretisation.h"
#include "StgFEM/SLE/LinearAlgebra/LinearAlgebra.h"
#include "StgFEM/SLE/SystemSetup/SystemSetup.h"

#include <stdio.h>
#include <stdlib.h>

struct _Node {
	Coord coord;
	double temp;
};

struct _Element {
	__FiniteElement_Element
};

struct _Particle {
	__IntegrationPoint
};

int main( int argc, char* argv[] ) {
	MPI_Comm			CommWorld;
	int				rank;
	int				numProcessors;
	int				procToWatch;
	Dictionary*			dictionary;
	Dictionary_Entry_Value*		currBC;
	Dictionary_Entry_Value*		bcList;
	Topology*			nTopology;
	ElementLayout*			eLayout;
	NodeLayout*			nLayout;
	MeshDecomp*			decomp;
	MeshLayout*			meshLayout;
	DofLayout*			dofLayout;
	ElementType_Register*		elementType_Register;
	ExtensionManager_Register*		extensionMgr_Register;
	FiniteElement_Mesh*		feMesh;
	WallVC*				wallVC;
	FieldVariable_Register*		fV_Register;
	FeVariable*			feVariable;
	Variable_Register*		variableRegister;
	SolutionVector*			solnVec;
	Index				i;
	Name				velocityName = "velocity";
	Stream*				stream;
	Index				addCount;
	Index*				addIndices;
	double*				addValues;
	
	/* Initialise MPI, get world info */
	MPI_Init( &argc, &argv );
	MPI_Comm_dup( MPI_COMM_WORLD, &CommWorld );
	MPI_Comm_size( CommWorld, &numProcessors );
	MPI_Comm_rank( CommWorld, &rank );
	
	StGermain_Init( &argc, &argv );
	StgFEM_Discretisation_Init( &argc, &argv );
	StgFEM_SLE_LinearAlgebra_Init( &argc, &argv );
	StgFEM_SLE_SystemSetup_Init( &argc, &argv );
	MPI_Barrier( CommWorld ); /* Ensures copyright info always come first in output */
	
	stream = Journal_Register (Info_Type, "myStream");
	/* set line below to true to True to activate outputting */
	Journal_Enable_TypedStream( DebugStream_Type, False );
	Stream_SetLevelBranch( StgFEM_Debug, 2 );
	Stream_EnableBranch( StgFEM_Debug, True );
	Stream_SetLevel( Journal_Register( DebugStream_Type, "FE.FE_Vectors.SolutionVector" ), 3 );
	Stream_SetLevel( Journal_Register( DebugStream_Type, "FE.FE_Matrices.FeEquationNumber.LM" ), 1 );


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
	Dictionary_Add( dictionary, "meshSizeI", Dictionary_Entry_Value_FromUnsignedInt( 4 ) );
	Dictionary_Add( dictionary, "meshSizeJ", Dictionary_Entry_Value_FromUnsignedInt( 2 ) );
	Dictionary_Add( dictionary, "meshSizeK", Dictionary_Entry_Value_FromUnsignedInt( 11 ) );
	Dictionary_Add( dictionary, "minX", Dictionary_Entry_Value_FromDouble( 0.0f ) );
	Dictionary_Add( dictionary, "minY", Dictionary_Entry_Value_FromDouble( 0.0f ) );
	Dictionary_Add( dictionary, "minZ", Dictionary_Entry_Value_FromDouble( 0.0f ) );
	Dictionary_Add( dictionary, "maxX", Dictionary_Entry_Value_FromDouble( 300.0f ) );
	Dictionary_Add( dictionary, "maxY", Dictionary_Entry_Value_FromDouble( 12.0f ) );
	Dictionary_Add( dictionary, "maxZ", Dictionary_Entry_Value_FromDouble( 300.0f ) );
	Dictionary_Add( dictionary, "allowPartitionOnElement", Dictionary_Entry_Value_FromBool( False ) );
	Dictionary_Add( dictionary, "buildElementNodeTbl", Dictionary_Entry_Value_FromBool( True ) );
	
	bcList = Dictionary_Entry_Value_NewList();
	currBC = Dictionary_Entry_Value_NewStruct();
	Dictionary_Entry_Value_AddMember( currBC, "name", Dictionary_Entry_Value_FromString( "vx" ) );
	Dictionary_Entry_Value_AddMember( currBC, "type", Dictionary_Entry_Value_FromString( "double" ) );
	Dictionary_Entry_Value_AddMember( currBC, "value", Dictionary_Entry_Value_FromDouble( -1.0f ) );
	Dictionary_Entry_Value_AddElement( bcList, currBC );
	currBC = Dictionary_Entry_Value_NewStruct();
	Dictionary_Entry_Value_AddMember( currBC, "wall", Dictionary_Entry_Value_FromString( "left" ) );
	Dictionary_Entry_Value_AddMember( currBC, "variables", bcList );
	Dictionary_Add( dictionary, "boundaryCondition", currBC );
	
	/* create the layout, dof and mesh to use */
	nTopology = (Topology*)IJK6Topology_New( "IJK6Topology", dictionary );
	eLayout = (ElementLayout*)ParallelPipedHexaEL_New( "PPHexaEL", 3, dictionary );
	nLayout = (NodeLayout*)CornerNL_New( "CornerNL", dictionary, eLayout, nTopology );
	decomp = (MeshDecomp*)HexaMD_New( "HexaMD", dictionary, MPI_COMM_WORLD, eLayout, nLayout );
	meshLayout = MeshLayout_New( "MeshLayout", eLayout, nLayout, decomp );
	
	extensionMgr_Register = ExtensionManager_Register_New();
	elementType_Register = ElementType_Register_New("testElementTypeRegister");
	ElementType_Register_Add( elementType_Register, (ElementType*)ConstantElementType_New("constant") );
	ElementType_Register_Add( elementType_Register, (ElementType*)BilinearElementType_New("bilinear") );
	ElementType_Register_Add( elementType_Register, (ElementType*)TrilinearElementType_New("trilinear") );
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
		&feMesh->nodeLocalCount, 
		(void**)&feMesh->nodeCoord, 
		variableRegister, 
		"vx", 
		"vy", 
		"vz" );
	Variable_Register_BuildAll(variableRegister);

	dofLayout = DofLayout_New( "dofLayout", variableRegister, decomp->nodeLocalCount );
	for (i = 0; i < decomp->nodeLocalCount; i++)
	{
		DofLayout_AddDof_ByVarName(dofLayout, "vx", i);
	}
	Build(dofLayout, 0, False);
	
	wallVC = WallVC_New( "WallVC", "boundaryCondition", variableRegister, NULL, dictionary, feMesh );

	/* Create the fe variable and solution vector */
	fV_Register = FieldVariable_Register_New();
	feVariable = FeVariable_New( velocityName, feMesh, NULL, dofLayout, wallVC, NULL, NULL, 3,
		False, StgFEM_Native_ImportExportType, StgFEM_Native_ImportExportType, fV_Register );
	solnVec = SolutionVector_New( velocityName, CommWorld, feVariable);
	
	/* Build and initialise system */
	Build( wallVC, 0, False );
	Build( solnVec, 0, False );
	Initialise( solnVec, 0, False );
	
	/* Set some test values into the Vector */
	addIndices = Memory_Alloc_Array_Unnamed( Index, feMesh->nodeLocalCount );
	addValues = Memory_Alloc_Array_Unnamed( double, feMesh->nodeLocalCount );
	for ( addCount=0; addCount < feMesh->nodeLocalCount; addCount++) {
		addIndices[addCount] = feVariable->eqNum->destinationArray[addCount][0];
		addValues[addCount] = (rank+1)*100 + addCount;
	}
	Vector_AddTo( solnVec->vector, addCount, addIndices, addValues );
	Vector_AssemblyBegin( solnVec->vector );
	Vector_AssemblyEnd( solnVec->vector );
	Memory_Free( addIndices );
	Memory_Free( addValues );
	
	/* Interpolate those values back to the Mesh */
	SolutionVector_UpdateSolutionOntoNodes( solnVec );
	
	if( rank == procToWatch ) {
		FeEquationNumber_PrintDestinationArray( feVariable->eqNum, stream );
		FeVariable_PrintLocalDiscreteValues( feVariable, stream );
	}
	
	/* Ok, now set some new values on the mesh, then try upload them onto the vector */
	{
		Node_LocalIndex		node_lI;
		Node_GlobalIndex	node_gI;
		Dof_Index		dof_I;
		double			value;
		
		for ( node_lI=0; node_lI < feMesh->nodeLocalCount; node_lI++ ) {
			node_gI = feMesh->nodeL2G[node_lI];
			for ( dof_I=0; dof_I < dofLayout->dofCounts[node_lI]; dof_I++ ) {
				value = 10 * node_gI + dof_I;
				DofLayout_SetValueDouble( dofLayout, node_lI, dof_I, value );
			}
		}	
	}

	if( rank == procToWatch ) {
		printf("2nd time round: testing new values on mesh loaded correctly onto vec:\n\n" );
		FeVariable_PrintLocalDiscreteValues( feVariable, stream );
	}	
	SolutionVector_LoadCurrentFeVariableValuesOntoVector( solnVec );

	if( rank == procToWatch ) {
		Dof_Index		dof_I;
		double*			vecValues = NULL;

		
		printf( "Local vector values are:\n" );
		Vector_Get( solnVec->vector, &vecValues );
		for ( dof_I = 0; dof_I < feVariable->eqNum->localEqNumsOwnedCount; dof_I++ ) {
			printf( "Vec unknown %d = %g\n", dof_I, vecValues[dof_I] );
		}

		Vector_Restore( solnVec->vector, &vecValues );
	}	
	
	/* Destroy stuff */
	Stg_Class_Delete( fV_Register );
	Stg_Class_Delete( wallVC );
	Stg_Class_Delete( feMesh );
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
	StgFEM_SLE_LinearAlgebra_Finalise();
	StgFEM_Discretisation_Finalise();
	StGermain_Finalise();
	
	/* Close off MPI */
	MPI_Finalize();
	
	return 0; /* success */
}
