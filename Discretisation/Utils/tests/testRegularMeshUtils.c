/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street, Melbourne, 3053, Australia.
**
** Authors:
**	Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**	Patrick D. Sunter, Software Engineer, VPAC. (pds@vpac.org)
**	Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**	Siew-Ching Tan, Software Engineer, VPAC. (siew@vpac.org)
**	Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
**	Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
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
** $Id: testRegularMeshUtils.c 4081 2007-04-27 06:20:07Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include "Base/Base.h"

#include "Discretisation/Geometry/Geometry.h"
#include "Discretisation/Shape/Shape.h"
#include "Discretisation/Mesh/Mesh.h"
#include "Discretisation/Utils/Utils.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct _Node
{
	double temp;
};

struct _Element
{
	double temp;
};


Mesh* buildMesh( unsigned nDims, unsigned* size, 
		     double* minCrds, double* maxCrds, 
		     ExtensionManager_Register* emReg )
{
	CartesianGenerator*	gen;
	Mesh*			mesh;
	unsigned		maxDecomp[3] = {0, 1, 1};

	gen = CartesianGenerator_New( "" );
	CartesianGenerator_SetDimSize( gen, nDims );
	CartesianGenerator_SetTopologyParams( gen, size, 0, NULL, maxDecomp );
	CartesianGenerator_SetGeometryParams( gen, minCrds, maxCrds );

	mesh = Mesh_New( "" );
	Mesh_SetExtensionManagerRegister( mesh, emReg );
	Mesh_SetGenerator( mesh, gen );

	Build( mesh, NULL, False );
	Initialise( mesh, NULL, False );

	KillObject( mesh->generator );

	return mesh;
}


int main(int argc, char *argv[])
{
	MPI_Comm			CommWorld;
	int				rank;
	int				procCount;
	int				procToWatch;
	ExtensionManager_Register*	extensionMgr_Register;
	Stream*				stream;

	unsigned	nDims = 3;
	unsigned	meshSize[3] = {6, 6, 6};
	double		minCrds[3] = {0.0, 0.0, 0.0};
	double		maxCrds[3] = {1.0, 1.0, 1.0};
	Mesh*		mesh;
	
	/* Initialise MPI, get world info */
	MPI_Init(&argc, &argv);
	MPI_Comm_dup( MPI_COMM_WORLD, &CommWorld );
	MPI_Comm_size(CommWorld, &procCount);
	MPI_Comm_rank(CommWorld, &rank);

	Base_Init( &argc, &argv );
	
	DiscretisationGeometry_Init( &argc, &argv );
	DiscretisationShape_Init( &argc, &argv );
	DiscretisationMesh_Init( &argc, &argv );

	stream = Journal_Register (Info_Type, "myStream");
	procToWatch = argc >= 2 ? atoi(argv[1]) : 0;
	
	extensionMgr_Register = ExtensionManager_Register_New();
	mesh = buildMesh( nDims, meshSize, minCrds, maxCrds, extensionMgr_Register );
	
	if (rank == procToWatch)
	{
		unsigned		currElementNodesCount=0;	
		unsigned*         	currElementNodes = NULL;
		unsigned          	element_dI = 0;
		unsigned             	refNode_eI = 0;
		unsigned		node_Diagonal = 0;
		unsigned		node_Diagonal_gI = 0;
		
		/* only use this while setting up the test
		   Print(mesh, stream);
				
		   Some tests involving RegularMeshUtils_GetDiagOppositeAcrossElementNodeIndex() */
		
		for (element_dI=0; element_dI < Mesh_GetDomainSize( mesh, nDims ); element_dI++) {
			Mesh_GetIncidence( mesh, nDims, element_dI, MT_VERTEX, 
					   &currElementNodesCount, &currElementNodes );
			
			for (refNode_eI = 0; refNode_eI < currElementNodesCount; refNode_eI++ ) {
				
				node_Diagonal = RegularMeshUtils_GetDiagOppositeAcrossElementNodeIndex(mesh, element_dI, 
					currElementNodes[refNode_eI]) ;
				node_Diagonal_gI = Mesh_DomainToGlobal( mesh, MT_VERTEX, node_Diagonal );
				/*print message stating: Element #, curr node #, diag opp node #*/
				printf("Element #: %d, Current Node #: %d, Diagonal Node #: %d, (%d) \n",
					element_dI, currElementNodes[refNode_eI], node_Diagonal, node_Diagonal_gI);
				
			}
		}
	}
	
	Stg_Class_Delete(mesh);
	
	DiscretisationMesh_Finalise();
	DiscretisationShape_Finalise();
	DiscretisationGeometry_Finalise();
	
	Base_Finalise();
	
	/* Close off MPI */
	MPI_Finalize();
	
	return 0; /* success */
}
