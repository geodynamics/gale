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
** $Id: testFeEquationNumber-LinkedDofs.c 860 2007-06-07 05:47:20Z LukeHodkinson $
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

void Test_FeEquationNumberRun_Regular( Dictionary* dictionary, void* context, IJK elSizes, 
				       unsigned rank, unsigned nProcs, unsigned procToWatch );

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
	Test_FeEquationNumberRun_Regular( dictionary, context, elSizes, rank, numProcessors, procToWatch );

	Stg_Class_Delete( context );
	Stg_Class_Delete( dictionary );
	Stg_Class_Delete( io_Handler );

	StgFEM_Discretisation_Finalise();
	StGermain_Finalise();

	/* Close off MPI */
	MPI_Finalize();
	
	return 0; /* success */
}


FeEquationNumber* buildEqNum( unsigned nProcs, unsigned* sizes, Dictionary* dict, Context* context, 
			      VariableCondition** bcs )
{
	CartesianGenerator*		gen;
	FeMesh*				feMesh;
	DofLayout*			dofs;
	FeEquationNumber*		eqNum;
	Variable_Register*		varReg;
	unsigned			maxDecomp[3] = {0, 1, 1};
	double				minCrd[3];
	double				maxCrd[3];
	char*				varNames[] = {"x", "y", "z"};
	LinkedDofInfo*			linkedDofInfo;

	minCrd[0] = minCrd[1] = minCrd[2] = 0.0;
	maxCrd[0] = maxCrd[1] = maxCrd[2] = 1.0;

	gen = CartesianGenerator_New( "" );
	CartesianGenerator_SetDimSize( gen, 3 );
	CartesianGenerator_SetTopologyParams( gen, sizes, 0, NULL, maxDecomp );
	CartesianGenerator_SetGeometryParams( gen, minCrd, maxCrd );
	CartesianGenerator_SetShadowDepth( gen, 1 );

	feMesh = FeMesh_New( "" );
	Mesh_SetGenerator( feMesh, gen );
	FeMesh_SetElementFamily( feMesh, "linear" );
	Stg_Component_Build( feMesh, NULL, False );

	varReg = Variable_Register_New();

	Variable_NewVector( "coords", Variable_DataType_Double, 3, 
			    &feMesh->topo->remotes[0]->nDomains, (void**)&feMesh->verts, 
			    varReg, 
			    varNames[0], 
			    varNames[1], 
			    varNames[2] );

	dofs = DofLayout_New( "", varReg, 0, feMesh );
	dofs->nBaseVariables = 3;
	dofs->baseVariables = Memory_Alloc_Array_Unnamed( Variable*, 3 );
	dofs->baseVariables[0] = Variable_Register_GetByName( varReg, varNames[0] );
	dofs->baseVariables[1] = Variable_Register_GetByName( varReg, varNames[1] );
	dofs->baseVariables[2] = Variable_Register_GetByName( varReg, varNames[2] );
	Stg_Component_Build( dofs, NULL, False );
	Stg_Component_Initialise( dofs, NULL, False );

	*bcs = (VariableCondition*)WallVC_New( "WallVC", "wallBC", varReg, NULL, dict, feMesh );
	Stg_Component_Build( *bcs, NULL, False );
	Stg_Component_Initialise( *bcs, NULL, False );

	linkedDofInfo = LinkedDofInfo_New( "linkedDofInfo", feMesh, dofs, dict );
	Stg_Component_Build( linkedDofInfo, context, False );

	eqNum = FeEquationNumber_New( "feEquationNumber", feMesh, dofs, *bcs, linkedDofInfo );
	Stg_Component_Build( eqNum, NULL, False );
	Stg_Component_Initialise( eqNum, NULL, False );

	return eqNum;
}


void Test_FeEquationNumberRun_Regular( Dictionary* dictionary, void* context, IJK elSizes, 
				       unsigned rank, unsigned nProcs, unsigned procToWatch )
{
	VariableCondition*		vc;
	FeEquationNumber*		feEquationNumber;
	Stream*				stream;

	if( rank == procToWatch ) printf("Creating Geometry, Decomps and Layouts:\n");
	Dictionary_Set( dictionary, "meshSizeI", Dictionary_Entry_Value_FromUnsignedInt( elSizes[0]+1 ) );
	Dictionary_Set( dictionary, "meshSizeJ", Dictionary_Entry_Value_FromUnsignedInt( elSizes[1]+1 ) );
	Dictionary_Set( dictionary, "meshSizeK", Dictionary_Entry_Value_FromUnsignedInt( elSizes[2]+1 ) );

	stream = Journal_Register (Info_Type, "myStream");

	if( rank == procToWatch ) printf("Creating F.E. Mesh:\n");
	if( rank == procToWatch ) printf("Creating Vars:\n");
	if( rank == procToWatch ) printf("Creating VC:\n");
	if( rank == procToWatch ) printf("Building:\n");
	if( rank == procToWatch ) printf("Building mesh:\n");
	if( rank == procToWatch ) printf("Building Variable Conditions:\n");

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
	if( rank == procToWatch ) printf("Building FE Eq num:\n");
	if( rank == procToWatch ) printf("Initialising:\n");
	if( rank == procToWatch ) printf("Building LM:\n");

	feEquationNumber = buildEqNum( nProcs, elSizes, dictionary, context, &vc );
	
	if( rank == procToWatch ) {
		Journal_Printf( stream, "V.C. applied: " );
		VariableCondition_PrintConcise( vc, stream );
		FeEquationNumber_PrintDestinationArray( feEquationNumber, stream );
		FeEquationNumber_PrintLocationMatrix( feEquationNumber, stream );
		Stg_Class_Print( feEquationNumber->linkedDofInfo, stream );
	}
	
	/* Destroy stuff */
	if( rank == procToWatch ) printf("Cleaning Up:\n");
	Stg_Class_Delete( feEquationNumber );
}
