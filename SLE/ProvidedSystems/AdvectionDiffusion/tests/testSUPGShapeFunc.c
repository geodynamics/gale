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
** $Id: testSUPGShapeFunc.c 656 2006-10-18 06:45:50Z SteveQuenette $
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

struct _Node {
	double velocity[3];
};

struct _Element {
	__FiniteElement_Element
};

struct _Particle {
	Coord coord;
};

void PutParticlesOnBoundary( Swarm* swarm, double value, Dimension_Index boundaryAxis, Dimension_Index valueAxis ) {
	Particle_InCellIndex       cParticle_I;
	for ( cParticle_I = 0 ; cParticle_I < swarm->cellParticleCountTbl[ 0 ] ; cParticle_I++ ) {
		IntegrationPoint*          particle = (IntegrationPoint*)Swarm_ParticleInCellAt( swarm, 0, cParticle_I );

		particle->xi[boundaryAxis] = 2.0 * (double)cParticle_I/(double)(swarm->cellParticleCountTbl[ 0 ] - 1) - 1.0;
		particle->xi[valueAxis]    = value;
	}
}
void CheckShapeFunc( AdvectionDiffusionSLE* advDiffSLE, AdvDiffResidualForceTerm* residualForceTerm, Swarm* gaussSwarm, Element_Index lElement_I, Dimension_Index dim, FiniteElement_Mesh* feMesh, Stream* stream ) {
	double**                   shapeFunc;
	Node_Index                 nodeCount;
	ElementType*               elementType;
	Particle_InCellIndex       cParticle_I;
	Cell_Index                 cell_I;
	IntegrationPoint*          particle;
	double                     Ni[8];
	Node_DomainIndex           node_I;

	shapeFunc = AdvDiffResidualForceTerm_BuildSUPGShapeFunctions( residualForceTerm, advDiffSLE, gaussSwarm, lElement_I, dim );

	/* Print Shape Funcions */
	cell_I = CellLayout_MapElementIdToCellId( gaussSwarm->cellLayout, lElement_I );
	elementType = FeMesh_ElementTypeAt( feMesh, lElement_I );
	nodeCount = elementType->nodeCount;
	for ( node_I = 0 ; node_I < nodeCount ; node_I++ ) {
		for ( cParticle_I = 0 ; cParticle_I < gaussSwarm->cellParticleCountTbl[ cell_I ] ; cParticle_I++ ) {
			particle = (IntegrationPoint*)Swarm_ParticleInCellAt( gaussSwarm, cell_I, cParticle_I );

			ElementType_EvaluateShapeFunctionsAt( elementType, particle->xi, Ni );

			Journal_Printf( stream, "%12.6g %12.6g %12.6g %12.6g %12.6g\n", (double)node_I, particle->xi[0], particle->xi[1], shapeFunc[ cParticle_I ][ node_I ], Ni[ node_I ] );
		}
	}
}

int main( int argc, char* argv[] ) {
	MPI_Comm                   CommWorld;
	int                        rank;
	int                        numProcessors;
	int                        procToWatch;
	Dictionary*                dictionary;
	Topology*                  nTopology;
	ElementLayout*             eLayout;
	NodeLayout*                nLayout;
	MeshDecomp*                decomp;
	MeshLayout*                meshLayout;
	DofLayout*                 dofs;
	ElementType_Register*      elementType_Register;
	ExtensionManager_Register* extensionMgr_Register;
	FiniteElement_Mesh*        feMesh;
	FeVariable*                feVariable;
	Variable_Register*         variableRegister;
	DiscretisationContext*     context;
	Node_DomainIndex           node_I;
	Stream*	                   stream;
	/* Swarm Stuff */
	Swarm*                     gaussSwarm;
	CellLayout*                singleCellLayout;
	ParticleLayout*            gaussParticleLayout;
	Dimension_Index            dimExists[] = {True,True,False};
	/* Advection Diffusion Stuff */
	AdvectionDiffusionSLE*     advDiffSLE;
	Element_LocalIndex         lElement_I;
	Particle_InCellIndex       particlesPerDim[] = {2,2,2};
	ForceVector*               residual;
	AdvDiffResidualForceTerm*  residualForceTerm;


	
	/* Initialise MPI, get world info */
	MPI_Init( &argc, &argv );
	MPI_Comm_dup( MPI_COMM_WORLD, &CommWorld );
	MPI_Comm_size( CommWorld, &numProcessors );
	MPI_Comm_rank( CommWorld, &rank );
	
	StGermain_Init( &argc, &argv );
	StgFEM_Discretisation_Init( &argc, &argv );
	StgFEM_SLE_SystemSetup_Init( &argc, &argv );
	StgFEM_SLE_ProvidedSystems_AdvectionDiffusion_Init( &argc, &argv );
	MPI_Barrier( CommWorld ); /* Ensures copyright info always come first in output */
	
	stream = Journal_Register (Info_Type, "myStream");
	Stream_RedirectFile( stream, "output/output.dat" );

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
	Dictionary_Add( dictionary, "dim", Dictionary_Entry_Value_FromUnsignedInt( 2 ) );
	Dictionary_Add( dictionary, "meshSizeI", Dictionary_Entry_Value_FromUnsignedInt( 7 ) );
	Dictionary_Add( dictionary, "meshSizeJ", Dictionary_Entry_Value_FromUnsignedInt( 7 ) );
	Dictionary_Add( dictionary, "meshSizeK", Dictionary_Entry_Value_FromUnsignedInt( 1 ) );
	Dictionary_Add( dictionary, "minX", Dictionary_Entry_Value_FromDouble( 0.0f ) );
	Dictionary_Add( dictionary, "minY", Dictionary_Entry_Value_FromDouble( 0.0f ) );
	Dictionary_Add( dictionary, "minZ", Dictionary_Entry_Value_FromDouble( 0.0f ) );
	Dictionary_Add( dictionary, "maxX", Dictionary_Entry_Value_FromDouble( 1.2f ) );
	Dictionary_Add( dictionary, "maxY", Dictionary_Entry_Value_FromDouble( 2.2f ) );
	Dictionary_Add( dictionary, "maxZ", Dictionary_Entry_Value_FromDouble( 1.2f ) );
	Dictionary_Add( dictionary, "allowPartitionOnElement", Dictionary_Entry_Value_FromBool( False ) );
	Dictionary_Add( dictionary, "buildElementNodeTbl", Dictionary_Entry_Value_FromBool( True ) );
	Dictionary_Add( dictionary, "gaussParticlesX", Dictionary_Entry_Value_FromUnsignedInt( 2 ) );
	Dictionary_Add( dictionary, "gaussParticlesY", Dictionary_Entry_Value_FromUnsignedInt( 2 ) );
	
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
	
	/* create the layout, dof and mesh to use */
	nTopology = (Topology*)IJK6Topology_New( "IJK6Topology", dictionary );
	eLayout = (ElementLayout*)ParallelPipedHexaEL_New( "ElementLayout", 2, dictionary );
	nLayout = (NodeLayout*)CornerNL_New( "CornerNL", dictionary, eLayout, nTopology );
	decomp = (MeshDecomp*)HexaMD_New( "HexaMD", dictionary, MPI_COMM_WORLD, eLayout, nLayout );
	meshLayout = MeshLayout_New( "MeshLayout", eLayout, nLayout, decomp );
	
	extensionMgr_Register = ExtensionManager_Register_New();
	elementType_Register = ElementType_Register_New("elementTypeRegister");
	ElementType_Register_Add( elementType_Register, (ElementType*)ConstantElementType_New("constantElementType") );
	ElementType_Register_Add( elementType_Register, (ElementType*)BilinearElementType_New("BilinearElementType") );
	ElementType_Register_Add( elementType_Register, (ElementType*)TrilinearElementType_New("TrilinearElementType") );
	feMesh = FiniteElement_Mesh_New( "testMesh", meshLayout, sizeof(Node), sizeof(Element), extensionMgr_Register,
		elementType_Register, dictionary );
	
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


	dofs = DofLayout_New( "dofLayout", variableRegister, decomp->nodeLocalCount );
	for (node_I = 0; node_I < decomp->nodeLocalCount; node_I++)	{
		DofLayout_AddDof_ByVarName(dofs, "vx", node_I);
		DofLayout_AddDof_ByVarName(dofs, "vy", node_I);
	}
	
	/* Create the finite element field variable*/
	feVariable    = FeVariable_New( "VelocityField", feMesh, NULL, dofs, NULL, NULL, NULL, context->dim, False, 
		StgFEM_Native_ImportExportType, StgFEM_Native_ImportExportType, context->fieldVariable_Register );

	/* Swarm stuff */
	if ( context->dim == 3 ) 
		dimExists[K_AXIS] = True;

	/* setup layout for swarm */
	singleCellLayout    = (CellLayout*)     SingleCellLayout_New( "SingleCellLayout", dimExists, NULL, NULL );
	gaussParticleLayout = (ParticleLayout*) GaussParticleLayout_New( "GaussParticleLayout", context->dim, particlesPerDim );

	/* construct swarm */
	gaussSwarm = Swarm_New(  
			"gaussSwarm",
			singleCellLayout, 
			gaussParticleLayout,
			context->dim,
			sizeof(IntegrationPoint), 
			context->extensionMgr_Register, 
			context->variable_Register, 
			context->communicator );

	residual = ForceVector_New( 
			"residual",
			feVariable,
			context->dim,
			context->entryPoint_Register,
			MPI_COMM_WORLD );

	residualForceTerm = AdvDiffResidualForceTerm_New( 
			"residualForceTerm",
			residual,
			gaussSwarm,
			NULL,
			feVariable,
			NULL /* Diffusivity Variable */,
			1.0,
			Exact );
	
	/* Create Advection Diffusion Solver */
	advDiffSLE = AdvectionDiffusionSLE_New( 
			"TestSLE",
			NULL,
			NULL,
			False,
			0,
			0,
			False,
			context->entryPoint_Register,
			context->communicator,
			feVariable, /* phi field */
			residual /* Residual Solution Vector */,
			NULL /* Mass Matrix */,
			context->dim,
			0.5,
			variableRegister,
			context->fieldVariable_Register );
			
	/* Build and initialise system */
	Build(dofs, 0, False);
	Build( feVariable, 0, False );
	/* Build Swarm */
	Build( singleCellLayout, 0, False );
	Build( gaussParticleLayout, 0, False );
	Build( gaussSwarm, 0, False );
	Build( feMesh, 0, False );
	Variable_Register_BuildAll(variableRegister);
	FeEquationNumber_BuildLocationMatrix( feVariable->eqNum );

	Initialise( feMesh, 0, False );
	Initialise( singleCellLayout, 0, False );
	Initialise( gaussParticleLayout, 0, False );
	Initialise( gaussSwarm, 0, False );
	Initialise(dofs, 0, False);
	Initialise( feVariable, 0, False );
	Variable_Update( Variable_Register_GetByName( variableRegister, "vx" ) );
	Variable_Update( Variable_Register_GetByName( variableRegister, "vy" ) );
	
	lElement_I = 0;

	Journal_Printf(stream, "#Checking pure diffusion\n");
	/* Apply some arbitrary initial conditions */
	for ( node_I = 0; node_I < decomp->nodeDomainCount; node_I++ ) {
		feMesh->node[node_I].velocity[0] = 0.0;
		feMesh->node[node_I].velocity[1] = 0.0;
		feMesh->node[node_I].velocity[2] = 0.0;
	}
	residualForceTerm->defaultDiffusivity = 1.0;
	CheckShapeFunc( advDiffSLE, residualForceTerm, gaussSwarm, lElement_I, context->dim, feMesh, stream );
	
	Journal_Printf(stream, "#Checking pure advection vx = 1.0\n");
	for ( node_I = 0; node_I < decomp->nodeDomainCount; node_I++ ) {
		feMesh->node[node_I].velocity[0] = 1.0;
		feMesh->node[node_I].velocity[1] = 0.0;
		feMesh->node[node_I].velocity[2] = 0.0;
	}
	residualForceTerm->defaultDiffusivity = 0.0;
	CheckShapeFunc( advDiffSLE, residualForceTerm, gaussSwarm, lElement_I, context->dim, feMesh, stream );
	
	Journal_Printf(stream, "#Checking pure advection vy = 1.0\n");
	for ( node_I = 0; node_I < decomp->nodeDomainCount; node_I++ ) {
		feMesh->node[node_I].velocity[0] = 0.0;
		feMesh->node[node_I].velocity[1] = 1.0;
		feMesh->node[node_I].velocity[2] = 0.0;
	}
	residualForceTerm->defaultDiffusivity = 0.0;
	CheckShapeFunc( advDiffSLE, residualForceTerm, gaussSwarm, lElement_I, context->dim, feMesh, stream );

	/* Destroy stuff */
	Stg_Class_Delete( gaussSwarm );
	Stg_Class_Delete( residual );
	Stg_Class_Delete( residualForceTerm );
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
	
	StgFEM_SLE_SystemSetup_Finalise();
	StgFEM_Discretisation_Finalise();
	StGermain_Finalise();
	
	/* Close off MPI */
	MPI_Finalize();
	
	return 0; /* success */
}
