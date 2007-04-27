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
** $Id: testElementType.c 822 2007-04-27 06:20:35Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include "StgFEM/Discretisation/Discretisation.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
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

FeMesh* buildFeMesh( unsigned nDims ) {
	CartesianGenerator*	gen;
	FeMesh*			feMesh;
	unsigned		maxDecomp[3] = {0, 1, 1};
	unsigned		sizes[3];
	double			minCrd[3];
	double			maxCrd[3];

	sizes[0] = sizes[1] = 6;
	sizes[2] = 1;
	minCrd[0] = minCrd[1] = minCrd[2] = 0.0;
	maxCrd[0] = maxCrd[1] = maxCrd[2] = 1.2;

	gen = CartesianGenerator_New( "" );
	CartesianGenerator_SetDimSize( gen, nDims );
	CartesianGenerator_SetTopologyParams( gen, sizes, 0, NULL, maxDecomp );
	CartesianGenerator_SetGeometryParams( gen, minCrd, maxCrd );
	CartesianGenerator_SetShadowDepth( gen, 0 );

	feMesh = FeMesh_New( "" );
	Mesh_SetGenerator( feMesh, gen );
	FeMesh_SetElementFamily( feMesh, "linear" );
	Build( feMesh, NULL, False );

	return feMesh;
}

int main( int argc, char* argv[] ) {
	MPI_Comm			CommWorld;
	int				rank;
	int				numProcessors;
	int				procToWatch;
	Dictionary*			dictionary;
	ExtensionManager_Register*		extensionMgr_Register;
	FeMesh*		feMesh;
	DiscretisationContext* context;
	Index				test_I;
	Stream*				stream;
	Coord               globalCoord;
	Coord               elLocalCoord;
	Coord               elLocalCoordGeneral;
	Element_Index       elementCoordIn;
	Node_Index          currElementNodeCount;
	const Index         maxTests            = 40;
	ElementType*        elementType;
	double              error               = 0.0;
	const double        tolerance           = 1.0e-8;
	XYZ                 errorVector;
	
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
	Dictionary_ReadAllParamFromCommandLine( dictionary, argc, argv );
	
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
	extensionMgr_Register = ExtensionManager_Register_New();
	feMesh = buildFeMesh( context->dim );
	Build( feMesh, 0, False );
	Initialise( feMesh, 0, False );

	srand48(0);
	for ( test_I = 0 ; test_I < maxTests ; test_I++ ) {
		globalCoord[ I_AXIS ] = drand48();
		globalCoord[ J_AXIS ] = drand48();
		globalCoord[ K_AXIS ] = drand48();
		Mesh_Algorithms_SearchElements( feMesh->algorithms, globalCoord, &elementCoordIn );
		currElementNodeCount = FeMesh_GetElementNodeSize( feMesh, elementCoordIn );
		elementType = FeMesh_GetElementType( feMesh, elementCoordIn );

		_ElementType_ConvertGlobalCoordToElLocal( elementType, feMesh, elementCoordIn, 
							  globalCoord, elLocalCoordGeneral );

		if ( context->dim == 2 ) {
			_BilinearElementType_ConvertGlobalCoordToElLocal( elementType, feMesh, elementCoordIn, 
									  globalCoord, elLocalCoord );
		}
		else {
			_TrilinearElementType_ConvertGlobalCoordToElLocal( elementType, feMesh, elementCoordIn, 
									   globalCoord, elLocalCoord );
		}

		StGermain_VectorSubtraction( errorVector, elLocalCoordGeneral, elLocalCoord, context->dim );

		error += StGermain_VectorMagnitude( errorVector, context->dim );
	}

	Journal_PrintBool( stream, error < tolerance );
	
	/* Destroy stuff */
	Stg_Class_Delete( feMesh );
	Stg_Class_Delete( extensionMgr_Register );
	Stg_Class_Delete( dictionary );
	
	StgFEM_Discretisation_Finalise();
	StGermain_Finalise();
	
	/* Close off MPI */
	MPI_Finalize();
	
	return 0; /* success */
}
