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
** $Id: testFeEquationNumber-LinkedDofs2.c 656 2006-10-18 06:45:50Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include "StgFEM/Discretisation/Discretisation.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <petsc.h>

struct _Node {
	Coord coord;
	double temp;
};

struct _Element {
	__FiniteElement_Element
};

struct _Particle {
	Coord coord;
};


void Test_FeEquationNumberRun_Regular( Dictionary* dictionary, const char* nLayoutType, IJK elSizes, Partition_Index rank, Partition_Index procToWatch );

int main( int argc, char* argv[] ) {
	MPI_Comm			CommWorld;
	int				rank;
	int				numProcessors;
	int				procToWatch;
	Dictionary*			dictionary;
	Dictionary_Entry_Value*		wallBC;
	Dictionary_Entry_Value*		variablesList;
	Dictionary_Entry_Value*		xVariable;
	XML_IO_Handler*			geomIO_Handler;
	Index					elSizes[3];
	Stream* 				feDebugStream;
	#if 0
	JournalFile*			file;
	#endif
	
	/* Initialise MPI, get world info */
	MPI_Init( &argc, &argv );
	MPI_Comm_dup( MPI_COMM_WORLD, &CommWorld );
	MPI_Comm_size( CommWorld, &numProcessors );
	MPI_Comm_rank( CommWorld, &rank );
	
	StGermain_Init( &argc, &argv );
	StgFEM_Discretisation_Init( &argc, &argv );
	MPI_Barrier( CommWorld ); /* Ensures copyright info always come first in output */

	Journal_Enable_TypedStream( DebugStream_Type, False );
	feDebugStream = Journal_Register( DebugStream_Type, "StgFEM" );
	Stream_EnableBranch( feDebugStream, True );
	Stream_SetLevelBranch( feDebugStream, 3 );
	Stream_EnableBranch( Journal_Register( DebugStream_Type, "StgFEM.StgFEM_Discretisation" ), True );
	#if 0
	file = CFile_New();
	if ( JournalFile_Open( file, "debugInfo" ) )
	{
		Journal_RegisterFile( file );
	}
	Stream_SetFile( Journal_GetTypedStream( DebugStream_Type ), file );
	Stream_SetFile( feDebugStream, file );
	Stream_SetFile( Journal_Register( DebugStream_Type, "StgFEM.StgFEM_Discretisation.FeEquationNumber" ), file );
	Stream_SetFile( Journal_Register( DebugStream_Type, "StgFEM.StgFEM_Discretisation.FeEquationNumber.LM" ), file );
	#endif
	
	if( argc >= 2 ) {
		procToWatch = atoi( argv[1] );
	}
	else {
		procToWatch = 0;
	}
	if( rank == procToWatch ) printf( "Watching rank: %i\n", rank );
	
	/* Read input */
	dictionary = Dictionary_New();
	geomIO_Handler = XML_IO_Handler_New();
	IO_Handler_ReadAllFromFile( geomIO_Handler, "data/triMeshAndBCs.xml", dictionary );
	Dictionary_Add( dictionary, "rank", Dictionary_Entry_Value_FromUnsignedInt( rank ) );
	Dictionary_Add( dictionary, "numProcessors", Dictionary_Entry_Value_FromUnsignedInt( numProcessors ) );
	Dictionary_Add( dictionary, "minX", Dictionary_Entry_Value_FromDouble( 0.0f ) );
	Dictionary_Add( dictionary, "minY", Dictionary_Entry_Value_FromDouble( 0.0f ) );
	Dictionary_Add( dictionary, "minZ", Dictionary_Entry_Value_FromDouble( 0.0f ) );
	Dictionary_Add( dictionary, "maxX", Dictionary_Entry_Value_FromDouble( 300.0f ) );
	Dictionary_Add( dictionary, "maxY", Dictionary_Entry_Value_FromDouble( 12.0f ) );
	Dictionary_Add( dictionary, "maxZ", Dictionary_Entry_Value_FromDouble( 300.0f ) );
	Dictionary_Add( dictionary, "allowPartitionOnElement", Dictionary_Entry_Value_FromBool( False ) );
	Dictionary_Add( dictionary, "allowUnbalancing", Dictionary_Entry_Value_FromBool( True ) );
	Dictionary_Add( dictionary, "buildElementNodeTbl", Dictionary_Entry_Value_FromBool( True ) );
	
	wallBC = Dictionary_Entry_Value_NewStruct();
	Dictionary_Add( dictionary, "wallBC", wallBC );
	Dictionary_Entry_Value_AddMember( wallBC, "type", Dictionary_Entry_Value_FromString( "WallVC" ) );
	Dictionary_Entry_Value_AddMember( wallBC, "wall", Dictionary_Entry_Value_FromString( "left" ) );
	variablesList = Dictionary_Entry_Value_NewList();
	Dictionary_Entry_Value_AddMember( wallBC, "variables", variablesList );
	xVariable = Dictionary_Entry_Value_NewStruct();
	Dictionary_Entry_Value_AddElement( variablesList, xVariable );
	Dictionary_Entry_Value_AddMember( xVariable, "name", Dictionary_Entry_Value_FromString( "x" ) );
	Dictionary_Entry_Value_AddMember( xVariable, "type", Dictionary_Entry_Value_FromString( "double" ) );
	Dictionary_Entry_Value_AddMember( xVariable, "value", Dictionary_Entry_Value_FromDouble( -1.0f ) );

	Dictionary_Add( dictionary, "meshSizeI", Dictionary_Entry_Value_FromUnsignedInt( 0 ) );
	Dictionary_Add( dictionary, "meshSizeJ", Dictionary_Entry_Value_FromUnsignedInt( 0 ) );
	Dictionary_Add( dictionary, "meshSizeK", Dictionary_Entry_Value_FromUnsignedInt( 0 ) );

	if( rank == procToWatch ) printf( "\n***  REGULAR node/element layout tests ***\n" );
	if( rank == procToWatch ) printf( "\n***  Balanced: 6*3*1 elements ***\n" );

	MPI_Barrier( CommWorld );

	elSizes[I_AXIS] = 6; elSizes[J_AXIS] = 3; elSizes[K_AXIS] = 1;
	Test_FeEquationNumberRun_Regular( dictionary, "corner", elSizes, rank, procToWatch );

	Stg_Class_Delete( dictionary );
	Stg_Class_Delete( geomIO_Handler );
	//if( rank == procToWatch ) Memory_Print();

	StgFEM_Discretisation_Finalise();
	StGermain_Finalise();

	/* Close off MPI */
	MPI_Finalize();
	
	return 0; /* success */
}


void Test_FeEquationNumberRun_Regular( Dictionary* dictionary, const char* nLayoutType, IJK elSizes,
	Partition_Index rank, Partition_Index procToWatch )
{
	Topology*			nTopology;
	ElementLayout*			eLayout;
	NodeLayout*			nLayout;
	MeshDecomp*			decomp;
	MeshLayout*			meshLayout;
	ElementType_Register*		elementType_Register;
	ExtensionManager_Register*		extensionMgr_Register;
	FiniteElement_Mesh*		feMesh;
	VariableCondition*		vc;
	FeEquationNumber*		feEquationNumber;
	char*				varNames[] = {"x", "y", "z"};
	Index				i, j;
	Node_GlobalIndex		nodeCount;
	Index				numVars = 3;
	Variable_Register*		variableRegister;
	DofLayout*			dofs;
	Stream*				stream;
	LinkedDofInfo*			linkedDofInfo = NULL;
	Dictionary_Entry_Value*		linkedDofSpecsList;
	Dictionary_Entry_Value*		linkSpecEntry;

	if( rank == procToWatch ) printf("Creating Geometry, Decomps and Layouts:\n");
	Dictionary_Set( dictionary, "meshSizeI", Dictionary_Entry_Value_FromUnsignedInt( elSizes[0]+1 ) );
	Dictionary_Set( dictionary, "meshSizeJ", Dictionary_Entry_Value_FromUnsignedInt( elSizes[1]+1 ) );
	Dictionary_Set( dictionary, "meshSizeK", Dictionary_Entry_Value_FromUnsignedInt( elSizes[2]+1 ) );
	/* create the layout, dof and mesh to use */
	eLayout = (ElementLayout*)ParallelPipedHexaEL_New( "PPHexaEL", 3, dictionary );
	if ("corner" == nLayoutType) {
		nTopology = (Topology*)IJK6Topology_New( "IJK6Topology", dictionary );
		nLayout = (NodeLayout*)CornerNL_New( "CornerNL", dictionary, eLayout, nTopology );
	}
	decomp = (MeshDecomp*)HexaMD_New( "HexaMD", dictionary, MPI_COMM_WORLD, eLayout, nLayout );
	meshLayout = MeshLayout_New( "MeshLayout", eLayout, nLayout, decomp );

	stream = Journal_Register (Info_Type, "myStream");

	extensionMgr_Register = ExtensionManager_Register_New();
	elementType_Register = ElementType_Register_New("elementTypeRegister");
	ElementType_Register_Add( elementType_Register, (ElementType*)ConstantElementType_New("constant") );
	ElementType_Register_Add( elementType_Register, (ElementType*)BilinearElementType_New("bilinear") );
	ElementType_Register_Add( elementType_Register, (ElementType*)TrilinearElementType_New("trilinear") );
	ElementType_Register_Add( elementType_Register, (ElementType*)LinearTriangleElementType_New("linear") );
	if( rank == procToWatch ) printf("Creating F.E. Mesh:\n");
	feMesh = FiniteElement_Mesh_New( "testMesh", meshLayout, sizeof(Node), sizeof(Element), extensionMgr_Register,
		elementType_Register, dictionary );
	
	/* Create variable register */
	variableRegister = Variable_Register_New();
	
	/* Create variables */
	if( rank == procToWatch ) printf("Creating Vars:\n");
	Variable_NewVector( 
		"coords", 
		Variable_DataType_Double, 
		3, 
		&feMesh->nodeLocalCount, 
		(void**)&feMesh->nodeCoord, 
		variableRegister, 
		varNames[0], 
		varNames[1], 
		varNames[2] );


	nodeCount = feMesh->layout->decomp->nodeLocalCount;
	dofs = DofLayout_New( "dofLayout", variableRegister, nodeCount );

	for (i = 0; i < nodeCount; i++)
	{
		for (j = 0; j < numVars; j++) {
			DofLayout_AddDof_ByVarName(dofs, varNames[j], i);
		}
	}
	Build(dofs, 0, False);
	
	if( rank == procToWatch ) printf("Creating VC:\n");
	if ( IrregEL_Type == meshLayout->elementLayout->type) {
		vc = (VariableCondition*) SetVC_New( "SetVC", "setBC", variableRegister, NULL, dictionary );
	}
	else {
		vc = (VariableCondition*) WallVC_New( "WallVC", "wallBC", variableRegister, NULL, dictionary, feMesh );
	}
	
	if( rank == procToWatch ) printf("Building:\n");
	/* Build and initialise system */
	if( rank == procToWatch ) printf("Building mesh:\n");
	Build( feMesh, 0, False );
	if( rank == procToWatch ) printf("Building Variable Conditions:\n");
	Build( vc, 0, False );
	
	linkedDofSpecsList = Dictionary_Entry_Value_NewList();
	Dictionary_Add( dictionary, "linkSpecifications", linkedDofSpecsList );
	#if 0
	linkSpecEntry = Dictionary_Entry_Value_NewStruct();
	Dictionary_Entry_Value_AddElement( linkedDofSpecsList, linkSpecEntry );
	Dictionary_Entry_Value_AddMember( linkSpecEntry, "dofSet", Dictionary_Entry_Value_FromUnsignedInt( 0 ) );
	Dictionary_Entry_Value_AddMember( linkSpecEntry, "wall", Dictionary_Entry_Value_FromString( "bottom" ) );
	Dictionary_Entry_Value_AddMember( linkSpecEntry, "dof", Dictionary_Entry_Value_FromUnsignedInt( 2 ) );
	#endif
	linkSpecEntry = Dictionary_Entry_Value_NewStruct();
	Dictionary_Entry_Value_AddElement( linkedDofSpecsList, linkSpecEntry );
	Dictionary_Entry_Value_AddMember( linkSpecEntry, "wallPair", Dictionary_Entry_Value_FromString( "left-right" ) );
	Dictionary_Entry_Value_AddMember( linkSpecEntry, "dof", Dictionary_Entry_Value_FromUnsignedInt( 1 ) );
	
	linkedDofInfo = LinkedDofInfo_New( "linkedDofInfo", feMesh, dofs, dictionary );
	Build( linkedDofInfo, 0, False );
	/* Demonstration of setting up linked dof info through code */

	if( rank == procToWatch ) printf("Creating EQ num:\n");
	/* Create the finite element equation number utility */
	feEquationNumber = FeEquationNumber_New( "feEquationNumber", feMesh, dofs, vc, linkedDofInfo );

	
	if( rank == procToWatch ) printf("Building FE Eq num:\n");
	FeEquationNumber_Build( feEquationNumber );
	
	if( rank == procToWatch ) printf("Initialising:\n");
	Initialise( feMesh, 0, False );
	FeEquationNumber_Initialise( feEquationNumber );
	
	if( rank == procToWatch ) printf("Building LM:\n");
	FeEquationNumber_BuildLocationMatrix( feEquationNumber );
	
	if( rank == procToWatch ) {
	
		Journal_Printf( stream, "V.C. applied: " );
		VariableCondition_PrintConcise( vc, stream );
		FeEquationNumber_PrintDestinationArray( feEquationNumber, stream );
		//FeEquationNumber_PrintLocationMatrix( feEquationNumber, stream );
		Print( linkedDofInfo, stream );
	}
	
	
	/* Destroy stuff */
	if( rank == procToWatch ) printf("Cleaning Up:\n");
	Stg_Class_Delete( feEquationNumber );
	Stg_Class_Delete( vc );
	Stg_Class_Delete( feMesh );
	Stg_Class_Delete( extensionMgr_Register );
	Stg_Class_Delete( elementType_Register );
	Stg_Class_Delete( variableRegister );
	Stg_Class_Delete( dofs );


	Stg_Class_Delete( meshLayout );
	Stg_Class_Delete( decomp );
	Stg_Class_Delete( nLayout );
	Stg_Class_Delete( eLayout );
	Stg_Class_Delete( nTopology );
}	
