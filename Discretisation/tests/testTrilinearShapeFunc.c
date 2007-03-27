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
** $Id: testTrilinearShapeFunc.c 656 2006-10-18 06:45:50Z SteveQuenette $
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
	double temp;
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
	ElementCellLayout*		elementCellLayout;
	ElementType_Register*		elementType_Register;
	ExtensionManager_Register*	extensionMgr_Register;
	FiniteElement_Mesh*		feMesh;
	WallVC*				wallVC;
	FeEquationNumber*		feEquationNumber;
	Variable_Register*		variableRegister;
	
	int				test=0;
	unsigned int			i,nElT;
	double				xi[3];
	double				Ni[8];
	double				xNi[8];	
	
	const Coord**			globalNodeCoordPtrs = NULL;
	Coord				globalCoord = {0,0,0};
	Coord				elLocalCoord = {0,0,0};
	ElementType*			elementType=NULL;
	Element_DomainIndex		element_dI=0;
	
	/* Initialise MPI, get world info */
	MPI_Init( &argc, &argv );
	MPI_Comm_dup( MPI_COMM_WORLD, &CommWorld );
	MPI_Comm_size( CommWorld, &numProcessors );
	MPI_Comm_rank( CommWorld, &rank );
	
	StGermain_Init( &argc, &argv );
	StgFEM_Discretisation_Init( &argc, &argv );
	MPI_Barrier( CommWorld ); /* Ensures copyright info always come first in output */
	
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
	Dictionary_Add( dictionary, "feMeshSizeI", Dictionary_Entry_Value_FromUnsignedInt( 3 ) );
	Dictionary_Add( dictionary, "feMeshSizeJ", Dictionary_Entry_Value_FromUnsignedInt( 3 ) );
	Dictionary_Add( dictionary, "feMeshSizeK", Dictionary_Entry_Value_FromUnsignedInt( 3 ) );
	Dictionary_Add( dictionary, "meshSizeI", Dictionary_Entry_Value_FromUnsignedInt( 3 ) );
	Dictionary_Add( dictionary, "meshSizeJ", Dictionary_Entry_Value_FromUnsignedInt( 3 ) );
	Dictionary_Add( dictionary, "meshSizeK", Dictionary_Entry_Value_FromUnsignedInt( 3 ) );
	Dictionary_Add( dictionary, "minX", Dictionary_Entry_Value_FromDouble( 0.0f ) );
	Dictionary_Add( dictionary, "minY", Dictionary_Entry_Value_FromDouble( 0.0f ) );
	Dictionary_Add( dictionary, "minZ", Dictionary_Entry_Value_FromDouble( 0.0f ) );
	Dictionary_Add( dictionary, "maxX", Dictionary_Entry_Value_FromDouble( 2.0f ) );
	Dictionary_Add( dictionary, "maxY", Dictionary_Entry_Value_FromDouble( 2.0f ) );
	Dictionary_Add( dictionary, "maxZ", Dictionary_Entry_Value_FromDouble( 2.0f ) );
	Dictionary_Add( dictionary, "buildElementNodeTbl", Dictionary_Entry_Value_FromBool( True ) );
	Dictionary_Add( dictionary, "allowPartitionOnElement", Dictionary_Entry_Value_FromBool( False ) );
	
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
	eLayout = (ElementLayout*)ParallelPipedHexaEL_New( "PPHexaEL", 3, dictionary );
	nLayout = (NodeLayout*)CornerNL_New( "CornerNL", dictionary, eLayout, nTopology );
	decomp = (MeshDecomp*)HexaMD_New( "HexaMD", dictionary, MPI_COMM_WORLD, eLayout, nLayout );
	meshLayout = MeshLayout_New( "MeshLayout", eLayout, nLayout, decomp );
	
	extensionMgr_Register = ExtensionManager_Register_New();
	elementType_Register = ElementType_Register_New("elementTypeRegister");
	ElementType_Register_Add( elementType_Register, (ElementType*)ConstantElementType_New("constElementType") );
	ElementType_Register_Add( elementType_Register, (ElementType*)BilinearElementType_New("BilinearElementType") );
	ElementType_Register_Add( elementType_Register, (ElementType*)TrilinearElementType_New("TrilinearElementType") );
	feMesh = FiniteElement_Mesh_New( "testMesh", meshLayout, sizeof(Node), sizeof(Element),
		extensionMgr_Register, elementType_Register, dictionary );
	Build( feMesh, 0, False );
	
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
	Variable_Register_BuildAll(variableRegister);


	dofs = DofLayout_New( "dofLayout", variableRegister, decomp->nodeLocalCount );
	for (i = 0; i < decomp->nodeLocalCount; i++)
	{
		DofLayout_AddDof_ByVarName(dofs, "x", i);
		DofLayout_AddDof_ByVarName(dofs, "y", i);
		DofLayout_AddDof_ByVarName(dofs, "z", i);
	}
	Build(dofs, 0, False);
	
	wallVC = WallVC_New( "WallVC", "boundaryCondition", variableRegister, NULL, dictionary, feMesh );

	/* create the element type register and feMesh to use */
	elementCellLayout = ElementCellLayout_New( "ElementCellLayout", feMesh );
	
	/* Create the finite element equation number utility */
	feEquationNumber = FeEquationNumber_New( "feEquationNumber", feMesh, dofs, (VariableCondition*)wallVC, NULL );
	
	/* Build and initialise system */
	Build( wallVC, 0, False );
	FeEquationNumber_Build( feEquationNumber );

	Initialise( feMesh, 0, False );
	FeEquationNumber_Initialise( feEquationNumber );
	
	xi[0] = xi[1] = xi[2] = 0.25;

	/* Test the shape functions */
	if( test == 0 ) {
		ElementType_EvaluateShapeFunctionsAt( FeMesh_ElementTypeAt( feMesh, 0 ), xi, Ni );
		for( i=0; i<8; i++ ) {
			printf("Ni[%d] = %f \n",i,Ni[i]);
		}
	}
	else {
		printf("Into hand calc. \n ");
		/* 
		Expected result in feEquationNumber ordering 
		Hand calculation
		*/
		xNi[0]=0.125*(0.75)*(0.75)*(0.75);
		xNi[3]=0.125*(0.75)*(1.25)*(0.75);
		xNi[2]=0.125*(1.25)*(1.25)*(0.75);
		xNi[1]=0.125*(1.25)*(0.75)*(0.75);
		xNi[4]=0.125*(0.75)*(0.75)*(1.25);
		xNi[7]=0.125*(0.75)*(1.25)*(1.25);
		xNi[6]=0.125*(1.25)*(1.25)*(1.25);
		xNi[5]=0.125*(1.25)*(0.75)*(1.25);
		for( i=0; i<8; i++ ) {
			printf("Ni[%d] = %f \n",i,xNi[i]);
		}
	}	

	printf( "Testing Global Coord -> ElLocal coord conversion:\n" );

	/* Ok, test several co-ordinates within the local element */
	for ( element_dI=0; element_dI < 2; element_dI++ ) {
		double ii, jj, kk;
		Coord refCoord;

		printf( "Element %d:\n", element_dI );
		elementType = FeMesh_ElementTypeAt( feMesh, element_dI );
		globalNodeCoordPtrs = Memory_Alloc_Array( const Coord*, feMesh->elementNodeCountTbl[element_dI], "globalNodeCoordPtrs" );
		Mesh_GetNodeCoordPtrsOfElement( feMesh, element_dI, (Coord**) globalNodeCoordPtrs );

		refCoord[0] = (*(globalNodeCoordPtrs[0]))[0];
		refCoord[1] = (*(globalNodeCoordPtrs[0]))[1];
		refCoord[2] = (*(globalNodeCoordPtrs[0]))[2];
		
		for ( ii = 0; ii < 1.01; ii += 0.2 ) {
			globalCoord[I_AXIS] = refCoord[I_AXIS] + ii;
			for ( jj = 0; jj < 1.01; jj += 0.2 ) {
				globalCoord[J_AXIS] = refCoord[J_AXIS] + jj;
				for ( kk = 0; kk < 1.01; kk += 0.2 ) {
					globalCoord[K_AXIS] = refCoord[K_AXIS] + kk;

					ElementType_ConvertGlobalCoordToElLocal(
						elementType,
						eLayout,
						globalNodeCoordPtrs,
						globalCoord,
						elLocalCoord );

					printf( "Global coord (%.1f,%.1f,%.1f) maps to elLocal (%.1f,%.1f,%.1f)\n",
						globalCoord[0], globalCoord[1], globalCoord[2],
						elLocalCoord[0], elLocalCoord[1], elLocalCoord[2] );	
				}		
			}
		}	
	}

	
	Memory_Free( globalNodeCoordPtrs );

	/* Destroy stuff */
	Stg_Class_Delete( feEquationNumber );
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
