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
** $Id: testIntegration.c 656 2006-10-18 06:45:50Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include "StgFEM/Discretisation/Discretisation.h"

#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include "petsc.h"

struct _Node {
	double value;
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
	Dictionary_Entry_Value*    bcList;
	Topology*                  nTopology;
	ElementLayout*             eLayout;
	NodeLayout*                nLayout;
	MeshDecomp*                decomp;
	MeshLayout*                meshLayout;
	BlockGeometry*             geometry;
	DofLayout*                 dofs;
	ElementType_Register*      elementType_Register;
	ExtensionManager_Register* extensionMgr_Register;
	FiniteElement_Mesh*        feMesh;
	WallVC*                    wallVC;
	FieldVariable_Register*    fV_Register;
	FeVariable*                feVariable;
	Variable_Register*         variableRegister;
	DiscretisationContext*     context;
	Node_DomainIndex           node_I;
	Stream*	                   stream;
	/* Swarm Stuff */
	Swarm*                     gaussSwarm;
	CellLayout*                singleCellLayout;
	ParticleLayout*            gaussParticleLayout;
	Bool                       dimExists[] = {True,True,False};
	double                     integral;
	double                     error;
	double                     planeValue;
	double                     analytic;
	Index                      plane_I;
	Index                      planeCount = 100;
	Particle_InCellIndex       particlesPerDim[] = {2,2,2};

	
	/* Initialise MPI, get world info */
	MPI_Init( &argc, &argv );
	MPI_Comm_dup( MPI_COMM_WORLD, &CommWorld );
	MPI_Comm_size( CommWorld, &numProcessors );
	MPI_Comm_rank( CommWorld, &rank );
	
	StGermain_Init( &argc, &argv );
	StgFEM_Discretisation_Init( &argc, &argv );
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
	Dictionary_Add( dictionary, "meshSizeI", Dictionary_Entry_Value_FromUnsignedInt( 17 ) );
	Dictionary_Add( dictionary, "meshSizeJ", Dictionary_Entry_Value_FromUnsignedInt( 17 ) );
	Dictionary_Add( dictionary, "meshSizeK", Dictionary_Entry_Value_FromUnsignedInt( 1 ) );
	Dictionary_Add( dictionary, "minX", Dictionary_Entry_Value_FromDouble( 0.0f ) );
	Dictionary_Add( dictionary, "minY", Dictionary_Entry_Value_FromDouble( 0.0f ) );
	Dictionary_Add( dictionary, "minZ", Dictionary_Entry_Value_FromDouble( 0.0f ) );
	Dictionary_Add( dictionary, "maxX", Dictionary_Entry_Value_FromDouble( 1.2f ) );
	Dictionary_Add( dictionary, "maxY", Dictionary_Entry_Value_FromDouble( 2.2f ) );
	Dictionary_Add( dictionary, "maxZ", Dictionary_Entry_Value_FromDouble( 1.2f ) );
	Dictionary_Add( dictionary, "allowPartitionOnElement", Dictionary_Entry_Value_FromBool( False ) );
	Dictionary_Add( dictionary, "buildElementNodeTbl", Dictionary_Entry_Value_FromBool( True ) );
	
	bcList = Dictionary_Entry_Value_NewList();
	currBC = Dictionary_Entry_Value_NewStruct();
	Dictionary_Entry_Value_AddMember( currBC, "name", Dictionary_Entry_Value_FromString( "testVariable" ) );
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
	
	/* create the layout, dof and mesh to use */
	nTopology = (Topology*)IJK6Topology_New( "IJK6Topology", dictionary );
	eLayout = (ElementLayout*)ParallelPipedHexaEL_New( "elementLayout", 2, dictionary );
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
	Variable_NewScalar( 
		"testVariable", 
		Variable_DataType_Double, 
		&feMesh->nodeDomainCount, 
		(void**)&feMesh->node, 
		variableRegister );
	Variable_Register_BuildAll(variableRegister);

	dofs = DofLayout_New( "dofLayout", variableRegister, decomp->nodeLocalCount );
	for (node_I = 0; node_I < decomp->nodeLocalCount; node_I++) 	{
		DofLayout_AddDof_ByVarName(dofs, "testVariable", node_I);
	}
	Build(dofs, 0, False);
	
	wallVC = WallVC_New( "WallVC", "boundaryCondition", variableRegister, NULL, dictionary, feMesh );

	/* Create the finite element field variable*/
	fV_Register = FieldVariable_Register_New();
	feVariable = FeVariable_New( "testField", feMesh, NULL, dofs, wallVC, NULL, NULL, context->dim,
		False, StgFEM_Native_ImportExportType, StgFEM_Native_ImportExportType, fV_Register );
	
	/* simulate solid body rotation */
	geometry = (BlockGeometry*)eLayout->geometry; 
	for ( node_I = 0; node_I < decomp->nodeDomainCount; node_I++ ) {
		double* pos = Mesh_CoordAt( feMesh, node_I );
		feMesh->node[node_I].value = 
			(1.0 - cos( 2.0 * M_PI /(geometry->max[ I_AXIS ] - geometry->min[ I_AXIS ]) * (pos[I_AXIS] - geometry->min[ I_AXIS ]) )) * 
			(1.0 - cos( 2.0 * M_PI /(geometry->max[ J_AXIS ] - geometry->min[ J_AXIS ]) * (pos[J_AXIS] - geometry->min[ J_AXIS ]) )) ; 
	}
	
	/* Build and initialise system */
	Build( wallVC, 0, False );
	Build( feVariable, 0, False );
	//Initialise( feVariable, 0, False );
	
	FeEquationNumber_BuildLocationMatrix( feVariable->eqNum );

	/* Swarm stuff */
	if ( context->dim == 3 ) 
		dimExists[K_AXIS] = True;

	/* setup layout for swarm */
	singleCellLayout    = (CellLayout*)     SingleCellLayout_New( "SingleCellLayout", dimExists, NULL, NULL );
	gaussParticleLayout = (ParticleLayout*) GaussParticleLayout_New( "GaussParticleCellLayout", context->dim, particlesPerDim );

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

	/* Build Swarm */
	Build( singleCellLayout, 0, False );
	Build( gaussParticleLayout, 0, False );
	Build( gaussSwarm, 0, False );
	Initialise( singleCellLayout, 0, False );
	Initialise( gaussParticleLayout, 0, False );
	Initialise( gaussSwarm, 0, False );
	
	/* Do Integration */
	integral = FeVariable_Integrate( feVariable, gaussSwarm );

	error = fabs( integral - (geometry->max[ I_AXIS ] - geometry->min[ I_AXIS ]) * 
		(geometry->max[ J_AXIS ] - geometry->min[ J_AXIS ])) ;

	if ( error > 1.0e-6 ) {
		Journal_Printf( stream, "Test Failed - Integral = %g Error = %g\n", integral, error );
	}	
	else {
		Journal_Printf( stream, "Test Passed\n");
	}	
	
	for ( plane_I = 0 ; plane_I < planeCount ; plane_I++ ) {
		planeValue = geometry->min[ I_AXIS ] + 
			(double)plane_I/(double)(planeCount - 1) * (geometry->max[ I_AXIS ] - geometry->min[ I_AXIS ]);

		integral = FeVariable_IntegratePlane( feVariable, I_AXIS, planeValue );
		analytic = 
			(1.0 - cos( 2.0*M_PI/(geometry->max[I_AXIS]-geometry->min[I_AXIS])*(planeValue - geometry->min[ I_AXIS ]) )) 
			* (geometry->max[ J_AXIS ] - geometry->min[ J_AXIS ]);

		error = fabs( integral - analytic );
		/*
		Journal_Printf( 
			stream, 
			"plane_I = %d, planeCount = %d, integral %f, an %f , planeValue %f\n", 
			plane_I, 
			planeCount, 
			integral, 
			analytic, 
			planeValue );
		*/

		if ( error > 0.045 ) {
			Journal_Printf( stream, "Test Failed FeVariable_IntegratePlane at plane = %g - error = %.g\n", 
					planeValue, error );
		}			
	}
	
	for ( plane_I = 0 ; plane_I < planeCount ; plane_I++ ) {
		planeValue = geometry->min[ J_AXIS ] + 
			(double)plane_I/(double)(planeCount - 1) * (geometry->max[ J_AXIS ] - geometry->min[ J_AXIS ]);

		integral = FeVariable_IntegratePlane( feVariable, J_AXIS, planeValue );
		analytic = 
			(1.0 - cos( 2.0*M_PI/(geometry->max[J_AXIS]-geometry->min[J_AXIS])*(planeValue - geometry->min[ J_AXIS ]) )) 
			* (geometry->max[ I_AXIS ] - geometry->min[ I_AXIS ]);

		error = fabs( integral - analytic );
		if ( error > 0.025 ) {
			Journal_Printf( stream, "Test Failed FeVariable_IntegratePlane at plane = %g - error = %.g\n", 
					planeValue, error );
		}			
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
}
