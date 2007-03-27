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
** $Id: testTrilinearShapeFuncLocalDerivs.c 656 2006-10-18 06:45:50Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include "StgFEM/Discretisation/Discretisation.h"

#include <stdio.h>
#include <stdlib.h>
#include "petsc.h"

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
	DofLayout*			dofs;
	ElementType_Register*		elementType_Register;
	ExtensionManager_Register*		extensionMgr_Register;
	FiniteElement_Mesh*		feMesh;
	WallVC*	wallVC;
	Variable_Register*		variableRegister;
	
	int test=0;
	unsigned int i,d,nElT;
	double xi[3];
	double** GNi;
	double xGNi[3][8];
	Stream* feDebugStream;
	
	/* Initialise MPI, get world info */
	MPI_Init( &argc, &argv );
	MPI_Comm_dup( MPI_COMM_WORLD, &CommWorld );
	MPI_Comm_size( CommWorld, &numProcessors );
	MPI_Comm_rank( CommWorld, &rank );
	
	StGermain_Init( &argc, &argv );
	StgFEM_Discretisation_Init( &argc, &argv );
	MPI_Barrier( CommWorld ); /* Ensures copyright info always come first in output */
	
	/* debugging stuff */
	Journal_Enable_TypedStream( DebugStream_Type, False );
	feDebugStream = Journal_Register( DebugStream_Type, "FE" );
	Stream_SetLevelBranch( feDebugStream, 2 );
	Stream_EnableBranch( feDebugStream, True );

	
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
	Dictionary_Add( dictionary, "meshSizeI", Dictionary_Entry_Value_FromUnsignedInt( 3 ) );
	Dictionary_Add( dictionary, "meshSizeJ", Dictionary_Entry_Value_FromUnsignedInt( 3 ) );
	Dictionary_Add( dictionary, "meshSizeK", Dictionary_Entry_Value_FromUnsignedInt( 3 ) );
	Dictionary_Add( dictionary, "minX", Dictionary_Entry_Value_FromDouble( 0.0f ) );
	Dictionary_Add( dictionary, "minY", Dictionary_Entry_Value_FromDouble( 0.0f ) );
	Dictionary_Add( dictionary, "minZ", Dictionary_Entry_Value_FromDouble( 0.0f ) );
	Dictionary_Add( dictionary, "maxX", Dictionary_Entry_Value_FromDouble( 100.0f ) );
	Dictionary_Add( dictionary, "maxY", Dictionary_Entry_Value_FromDouble( 200.0f ) );
	Dictionary_Add( dictionary, "maxZ", Dictionary_Entry_Value_FromDouble( 400.0f ) );
	Dictionary_Add( dictionary, "buildElementNodeTbl", Dictionary_Entry_Value_FromBool( True ) );
	
	Dictionary_Add( dictionary, "nElementTypes", Dictionary_Entry_Value_FromUnsignedInt( 3 ) );
	nElT = Dictionary_Entry_Value_AsUnsignedInt( Dictionary_Get( dictionary,"nElementTypes") );
	printf("my val = %u \n",nElT);
			
	bcList = Dictionary_Entry_Value_NewList();
	currBC = Dictionary_Entry_Value_NewStruct();
	Dictionary_Entry_Value_AddMember( currBC, "name", Dictionary_Entry_Value_FromString( "x" ) );
	Dictionary_Entry_Value_AddMember( currBC, "type", Dictionary_Entry_Value_FromString( "double" ) );
	Dictionary_Entry_Value_AddMember( currBC, "value", Dictionary_Entry_Value_FromDouble( -1.0f ) );
	Dictionary_Entry_Value_AddElement( bcList, currBC );
	currBC = Dictionary_Entry_Value_NewStruct();
	Dictionary_Entry_Value_AddMember( currBC, "wall", Dictionary_Entry_Value_FromString( "left" ) );
	Dictionary_Entry_Value_AddMember( currBC, "variables", bcList );
	Dictionary_Add( dictionary, "boundaryCondition", currBC );
	
	/* create the layout, dof and mesh to use */
	nTopology = (Topology*)IJK6Topology_New( "IJK6Topology", dictionary );
	eLayout = (ElementLayout*)ParallelPipedHexaEL_New( "PPHexa_El", 3, dictionary );
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
	Variable_NewVector( 
		"coords", 
		Variable_DataType_Double, 
		3, 
		&feMesh->nodeLocalCount, 
		(void**)&feMesh->nodeCoord, 
		variableRegister, 
		"x", 
		"y", 
		"z" );

	dofs = DofLayout_New( "dofLayout", variableRegister, decomp->nodeLocalCount );
	for (i = 0; i < decomp->nodeLocalCount; i++)
	{
		DofLayout_AddDof_ByVarName(dofs, "x", i);
		DofLayout_AddDof_ByVarName(dofs, "y", i);
		DofLayout_AddDof_ByVarName(dofs, "z", i);
	}
	Build(dofs, 0, False);
	
	wallVC = WallVC_New( "WallVC", "boundaryCondition", variableRegister, NULL, dictionary, feMesh );

	
	/* Build and initialise system */
	Build( feMesh, 0, False );
	Build( wallVC, 0, False );
	
	Initialise( feMesh, 0, False );
	
	xi[0] = xi[1] = xi[2] = 0.25;
		
	/* Alloc */
	GNi = Memory_Alloc_2DArray( double, 3, 8, "GNi" );
	
	/* Test the local shape function derivatives */
	ElementType_EvaluateShapeFunctionLocalDerivsAt( FeMesh_ElementTypeAt( feMesh, 0 ), xi, GNi );
	
	if( test == 0 ) {
		for( d=0; d<3; d++ ) {
			for( i=0; i<8; i++ ) {
				printf("GNi[%d][%d] = %f \n",d,i,GNi[d][i]);
			}
			printf("\n");
		}	
	}
	else {
		/* Hand calc */
		printf("Into hand calc \n");
		xGNi[0][0]=-0.125*(0.75)*(0.75);
		xGNi[0][3]=-0.125*(1.25)*(0.75);
		xGNi[0][2]= 0.125*(1.25)*(0.75);
		xGNi[0][1]= 0.125*(0.75)*(0.75);
		xGNi[0][4]=-0.125*(0.75)*(1.25);
		xGNi[0][7]=-0.125*(1.25)*(1.25);
		xGNi[0][6]= 0.125*(1.25)*(1.25);
		xGNi[0][5]= 0.125*(0.75)*(1.25);
		
		xGNi[1][0]=-0.125*(0.75)*(0.75);
		xGNi[1][3]= 0.125*(0.75)*(0.75);
		xGNi[1][2]= 0.125*(1.25)*(0.75);
		xGNi[1][1]=-0.125*(1.25)*(0.75);
		xGNi[1][4]=-0.125*(0.75)*(1.25);
		xGNi[1][7]= 0.125*(0.75)*(1.25);
		xGNi[1][6]= 0.125*(1.25)*(1.25);
		xGNi[1][5]=-0.125*(1.25)*(1.25);
		
		xGNi[2][0]=-0.125*(0.75)*(0.75);
		xGNi[2][3]=-0.125*(0.75)*(1.25);
		xGNi[2][2]=-0.125*(1.25)*(1.25);
		xGNi[2][1]=-0.125*(1.25)*(0.75);
		xGNi[2][4]= 0.125*(0.75)*(0.75);
		xGNi[2][7]= 0.125*(0.75)*(1.25);
		xGNi[2][6]= 0.125*(1.25)*(1.25);
		xGNi[2][5]= 0.125*(1.25)*(0.75);
		
		
		for( d=0; d<3; d++ ) {			
				for( i=0; i<8; i++ ) {
				printf("GNi[%d][%d] = %f \n",d,i,xGNi[d][i]);
			}
			printf("\n");	
		}
	}		
	
	Memory_Free(GNi); 
	
	
	
	/* Destroy stuff */
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
