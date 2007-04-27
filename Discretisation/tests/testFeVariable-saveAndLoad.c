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
** $Id: testFeVariable-saveAndLoad.c 822 2007-04-27 06:20:35Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include "StgFEM/Discretisation/Discretisation.h"

#include <stdio.h>
#include <stdlib.h>
#include "petsc.h"

struct _Node {
	double velocity[3];
};

struct _Element {
	__FiniteElement_Element
};

struct _Particle {
	Coord coord;
};

int main( int argc, char* argv[] ) {
#if 0
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
	WallVC*				wallVC;
	FieldVariable_Register*		fV_Register;
	FeVariable*			feVariable;
	Variable_Register*		variableRegister;
	DiscretisationContext* context;
	Index				i;
	Stream*				stream;
	
	/* Initialise MPI, get world info */
	MPI_Init( &argc, &argv );
	MPI_Comm_dup( MPI_COMM_WORLD, &CommWorld );
	MPI_Comm_size( CommWorld, &numProcessors );
	MPI_Comm_rank( CommWorld, &rank );
	
	StGermain_Init( &argc, &argv );
	StgFEM_Discretisation_Init( &argc, &argv );
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
	Dictionary_Add( dictionary, "rank", Dictionary_Entry_Value_FromUnsignedInt( rank ) );
	Dictionary_Add( dictionary, "numProcessors", Dictionary_Entry_Value_FromUnsignedInt( numProcessors ) );
	Dictionary_Add( dictionary, "meshSizeI", Dictionary_Entry_Value_FromUnsignedInt( 7 ) );
	Dictionary_Add( dictionary, "meshSizeJ", Dictionary_Entry_Value_FromUnsignedInt( 7 ) );
	Dictionary_Add( dictionary, "meshSizeK", Dictionary_Entry_Value_FromUnsignedInt( 1 ) );
	Dictionary_Add( dictionary, "minX", Dictionary_Entry_Value_FromDouble( 0.0f ) );
	Dictionary_Add( dictionary, "minY", Dictionary_Entry_Value_FromDouble( 0.0f ) );
	Dictionary_Add( dictionary, "minZ", Dictionary_Entry_Value_FromDouble( 0.0f ) );
	Dictionary_Add( dictionary, "maxX", Dictionary_Entry_Value_FromDouble( 1.2f ) );
	Dictionary_Add( dictionary, "maxY", Dictionary_Entry_Value_FromDouble( 1.2f ) );
	Dictionary_Add( dictionary, "maxZ", Dictionary_Entry_Value_FromDouble( 1.2f ) );
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

	/* Create Context */
	context = _DiscretisationContext_New( 
			sizeof(DiscretisationContext), 
			DiscretisationContext_Type, 
			_DiscretisationContext_Delete, 
			_DiscretisationContext_Print,
			NULL,
			NULL,
			NULL,
			NULL,
			NULL,
			NULL,
			NULL,
			"discretisationContext",
			True,
			NULL,
			0,0,
			CommWorld, dictionary );
	context->dim = Dictionary_GetUnsignedInt_WithDefault( dictionary, "dim", 2 );
	
	/* create the layout, dof and mesh to use */
	nTopology = (Topology*)IJK6Topology_New( "IJK6Topology", dictionary );
	eLayout = (ElementLayout*)ParallelPipedHexaEL_New( "PPHexaEL", context->dim, dictionary );
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
	Build( feMesh, 0, False );
	Initialise( feMesh, 0, False );
	
	/* Create variable register */
	variableRegister = Variable_Register_New();
	
	/* Create variables */
	Variable_NewVector( 
		"velocity", 
		Variable_DataType_Double, 
		3, 
		&feMesh->nodeDomainCount, 
		(void**)&feMesh->node, 
		variableRegister, 
		"vx", 
		"vy", 
		"vz" );
	Variable_Register_BuildAll(variableRegister);

	dofs = DofLayout_New( "dofLayout", variableRegister, decomp->nodeLocalCount );
	for (i = 0; i < decomp->nodeLocalCount; i++)
	{
		DofLayout_AddDof_ByVarName(dofs, "vx", i);
		DofLayout_AddDof_ByVarName(dofs, "vy", i);
	}
	Build(dofs, 0, False);
	
	wallVC = WallVC_New( "WallVC", "boundaryCondition", variableRegister, NULL, dictionary, feMesh );

	/* Create the finite element field variable*/
	fV_Register = FieldVariable_Register_New();
	feVariable = FeVariable_New( "velocity", feMesh, NULL, dofs, wallVC, NULL, NULL, context->dim,
		True,StgFEM_Native_ImportExportType, StgFEM_Native_ImportExportType, fV_Register );
	
	/* Apply some arbitrary initial conditions */
	for ( i = 0; i < decomp->nodeDomainCount; i++ ) {
		const double pi=acos(-1.0);
		double* pos = Mesh_CoordAt( feMesh, i );
		feMesh->node[i].velocity[0] = ( 1.0 - pos[1] )+ 0.1 * cos( pi * pos[0] ) * sin(  pi * pos[1]  );
		feMesh->node[i].velocity[1] = ( 1.0 - pos[1] )+ 0.1 * cos( pi * pos[0] ) * sin(  pi * pos[1]  );
		feMesh->node[i].velocity[2] = 0;
	}
	
	/* Build and initialise system */
	Build( wallVC, 0, False );
	Build( feVariable, 0, False );
	/* Initialise( feVariable, 0, False ); */
	
	if( rank == procToWatch ) {
		FeVariable_PrintLocalDiscreteValues_2dBox( feVariable, stream );
		Journal_Printf( stream, "\n" );
	}	

	FeVariable_SaveToFile( feVariable, "output/", 0 );
	
	for ( i = 0; i < decomp->nodeDomainCount; i++ ) {
		feMesh->node[i].velocity[0] = 0;
		feMesh->node[i].velocity[1] = 0;
		feMesh->node[i].velocity[2] = 0;
	}
	
	FeVariable_ReadFromFile( feVariable, "output/", 0 );
	
	if( rank == procToWatch ) {
		Journal_Printf( stream, "**After save, zero then reload**\n\n" );
		FeVariable_PrintLocalDiscreteValues_2dBox( feVariable, stream );
	}
	
		
	/* Destroy stuff */
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
#endif
}
