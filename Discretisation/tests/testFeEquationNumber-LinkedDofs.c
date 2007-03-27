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
** $Id: testFeEquationNumber-LinkedDofs.c 732 2007-02-07 00:38:34Z PatrickSunter $
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

void _SetDt( void* context, double dt ) {
}

void Test_FeEquationNumberRun_Regular( Dictionary* dictionary, void* context, const char* nLayoutType, IJK elSizes, Partition_Index rank, Partition_Index procToWatch );

int main( int argc, char* argv[] ) {
	MPI_Comm			CommWorld;
	int				rank;
	int				numProcessors;
	int				procToWatch;
	Dictionary*			dictionary;
	Dictionary*                     componentDict;
	XML_IO_Handler*			io_Handler;
	Index				elSizes[3];
	Stream* 			feDebugStream;
	Stream* 			stream;
	DiscretisationContext*          context;
	
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

	if( argc >= 2 ) {
		procToWatch = atoi( argv[1] );
	}
	else {
		procToWatch = 0;
	}
	stream = Journal_Register( Info_Type, __FILE__ );
	Stream_SetPrintingRank( stream, procToWatch );
	Stream_SetAutoFlush( stream, True );
	Journal_Printf( stream, "Watching rank: %i\n", rank );
	
	/* Make sure print statements come in order */
	MPI_Barrier( CommWorld );

	/* Read input */
	dictionary = Dictionary_New();
	Dictionary_Add( dictionary, "outputPath", Dictionary_Entry_Value_FromString( "./output" ) );
	io_Handler = XML_IO_Handler_New();
	IO_Handler_ReadAllFromFile( io_Handler, "./data/testFeEquationNumber-LinkedDofs.xml", dictionary );

	Dictionary_Add( dictionary, "meshSizeI", Dictionary_Entry_Value_FromUnsignedInt( 0 ) );
	Dictionary_Add( dictionary, "meshSizeJ", Dictionary_Entry_Value_FromUnsignedInt( 0 ) );
	Dictionary_Add( dictionary, "meshSizeK", Dictionary_Entry_Value_FromUnsignedInt( 0 ) );

	/* Create context */
	context = _DiscretisationContext_New(
		sizeof(DiscretisationContext),
		DiscretisationContext_Type,
		_DiscretisationContext_Delete,
		_DiscretisationContext_Print,
		NULL,
		NULL,
		NULL,
		_AbstractContext_Build,
		_AbstractContext_Initialise,
		_AbstractContext_Execute,
		_AbstractContext_Destroy,
		"context",
		True,
		_SetDt,
		0,
		0,
		CommWorld,
		dictionary );
	
	componentDict = Dictionary_GetDictionary( dictionary, "components" );
	context->CF = Stg_ComponentFactory_New( dictionary, componentDict, context->register_Register );
	Stream_SetPrintingRank( context->CF->infoStream, procToWatch );

	LiveComponentRegister_Add( context->CF->LCRegister, (Stg_Component*) context );
	
	Stg_ComponentFactory_CreateComponents( context->CF );
	Stg_ComponentFactory_ConstructComponents( context->CF, 0 /* dummy */ );

	/* Make sure print statements come in order */
	MPI_Barrier( CommWorld );

	Journal_Printf( stream, "\n***  REGULAR node/element layout tests ***\n" );
	Journal_Printf( stream, "\n***  Balanced: 6*1*1 elements ***\n" );
	MPI_Barrier( CommWorld );

	elSizes[I_AXIS] = 3; elSizes[J_AXIS] = 3; elSizes[K_AXIS] = 3;
	Test_FeEquationNumberRun_Regular( dictionary, context, "corner", elSizes, rank, procToWatch );

	Stg_Class_Delete( context );
	Stg_Class_Delete( dictionary );
	Stg_Class_Delete( io_Handler );

	StgFEM_Discretisation_Finalise();
	StGermain_Finalise();

	/* Close off MPI */
	MPI_Finalize();
	
	return 0; /* success */
}


void Test_FeEquationNumberRun_Regular( Dictionary* dictionary, void* context, const char* nLayoutType, IJK elSizes,
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
	int                             dim = Dictionary_GetUnsignedInt_WithDefault( dictionary, "dim", 3 );
	#if 0
	IndexSet*			bottomSet;
	IndexSet*			leftSet;
	IndexSet*			rightSet;
	#endif

	if( rank == procToWatch ) printf("Creating Geometry, Decomps and Layouts:\n");
	Dictionary_Set( dictionary, "meshSizeI", Dictionary_Entry_Value_FromUnsignedInt( elSizes[0]+1 ) );
	Dictionary_Set( dictionary, "meshSizeJ", Dictionary_Entry_Value_FromUnsignedInt( elSizes[1]+1 ) );
	Dictionary_Set( dictionary, "meshSizeK", Dictionary_Entry_Value_FromUnsignedInt( elSizes[2]+1 ) );
	/* create the layout, dof and mesh to use */
	eLayout = (ElementLayout*)ParallelPipedHexaEL_New( "PPHexaEL", dim, dictionary );
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
	Build( feMesh, context, False );
	if( rank == procToWatch ) printf("Building Variable Conditions:\n");
	Build( vc, context, False );
	
	linkedDofInfo = LinkedDofInfo_New( "linkedDofInfo", feMesh, dofs, dictionary );
	Build( linkedDofInfo, context, False );
	/* Demonstration of setting up linked dof info through code */
	#if 0
	LinkedDofInfo_AddDofSet( linkedDofInfo );
	LinkedDofInfo_AddDofSet( linkedDofInfo );
	bottomSet = RegularMeshUtils_CreateGlobalBottomSet( feMesh );
	leftSet = RegularMeshUtils_CreateGlobalLeftSet( feMesh );
	rightSet = RegularMeshUtils_CreateGlobalRightSet( feMesh );
	
	LinkedDofInfo_AddDofsToSet_FromIndexSet( linkedDofInfo, 0, bottomSet, 1 );
	LinkedDofInfo_AddDofsToSet_FromIndexSet( linkedDofInfo, 1, leftSet, 2 );
	LinkedDofInfo_AddDofsToSet_FromIndexSet( linkedDofInfo, 1, rightSet, 2 );
	#endif

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
		FeEquationNumber_PrintLocationMatrix( feEquationNumber, stream );
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
