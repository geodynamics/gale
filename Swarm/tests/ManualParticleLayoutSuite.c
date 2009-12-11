/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street, Melbourne, 3053, Australia.
**
** Authors:
**   Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**   Patrick D. Sunter, Software Engineer, VPAC. (pds@vpac.org)
**   Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**   Siew-Ching Tan, Software Engineer, VPAC. (siew@vpac.org)
**   Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
**   Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
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
** Role:
**   Tests the ManualParticleLayoutSuite
**
** $Id: testTemplate.c 3462 2006-02-19 06:53:24Z WalterLandry $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>

#include "pcu/pcu.h"
#include <StGermain/StGermain.h>
#include "StgDomain/Geometry/Geometry.h"
#include "StgDomain/Shape/Shape.h"
#include "StgDomain/Mesh/Mesh.h"
#include "StgDomain/Utils/Utils.h"
#include "StgDomain/Swarm/Swarm.h"

#include "ManualParticleLayoutSuite.h"

struct _Particle {
	__IntegrationPoint
};

typedef struct {
	MPI_Comm comm;
	unsigned rank;
	unsigned nProcs;
} ManualParticleLayoutSuiteData;

Mesh* ManualParticleLayoutSuite_BuildMesh( unsigned nDims, unsigned* size, double* minCrds, double* maxCrds, ExtensionManager_Register* emReg ) {
	CartesianGenerator*	gen;
	Mesh*						mesh;

	gen = CartesianGenerator_New( "", NULL );
	CartesianGenerator_SetDimSize( gen, nDims );
	CartesianGenerator_SetTopologyParams( gen, size, 0, NULL, NULL );
	CartesianGenerator_SetGeometryParams( gen, minCrds, maxCrds );

	mesh = Mesh_New( "", NULL );
	Mesh_SetExtensionManagerRegister( mesh, emReg );
	Mesh_SetGenerator( mesh, gen );

	Stg_Component_Build( mesh, NULL, False );
	Stg_Component_Initialise( mesh, NULL, False );

	FreeObject( mesh->generator );

	return mesh;
}

void ManualParticleLayoutSuite_Setup( ManualParticleLayoutSuiteData* data ) {
	Journal_Enable_AllTypedStream( False );

	/* MPI Initializations */
	data->comm = MPI_COMM_WORLD;  
	MPI_Comm_rank( data->comm, &data->rank );
	MPI_Comm_size( data->comm, &data->nProcs );
}

void ManualParticleLayoutSuite_Teardown( ManualParticleLayoutSuiteData* data ) {
	Journal_Enable_AllTypedStream( True );
}

void ManualParticleLayoutSuite_TestManualParticle( ManualParticleLayoutSuiteData* data ) {
	unsigned							nDims = 3;
	unsigned							meshSize[3] = {4, 2, 1};
	double							minCrds[3] = {0.0, 0.0, 0.0};
	double							maxCrds[3] = {1.0, 1.0, 1.0};
	ExtensionManager_Register*	extensionMgr_Register;
	Mesh*								mesh;
	GaussParticleLayout*			gaussParticleLayout;
	ElementCellLayout*			elementCellLayout;
	Swarm*							swarm;
	Dictionary*						dictionary;
	Dictionary_Entry_Value*		particlePositionsList;
	Dictionary_Entry_Value*		particlePositionEntry;
	Dimension_Index				dim;
	ManualParticleLayout*		particleLayout;
	int								procToWatch = data->nProcs > 1 ? 1 : 0;
	Stream*							stream;
	char								expected_file[PCU_PATH_MAX];
	
	if( data->rank == procToWatch ) {
		/* Dictionary Initialization */
		dictionary = Dictionary_New();
		particlePositionsList = Dictionary_Entry_Value_NewList();
		Dictionary_Add( dictionary, "manualParticlePositions", particlePositionsList );
	
		particlePositionEntry = Dictionary_Entry_Value_NewStruct();
		Dictionary_Entry_Value_AddElement( particlePositionsList, particlePositionEntry );
		Dictionary_Entry_Value_AddMember( particlePositionEntry, "x", Dictionary_Entry_Value_FromDouble( 0.4 ) );
		Dictionary_Entry_Value_AddMember( particlePositionEntry, "y", Dictionary_Entry_Value_FromDouble( 0.3 ) );
		Dictionary_Entry_Value_AddMember( particlePositionEntry, "z", Dictionary_Entry_Value_FromDouble( 0.2 ) );
	
		particlePositionEntry = Dictionary_Entry_Value_NewStruct();
		Dictionary_Entry_Value_AddElement( particlePositionsList, particlePositionEntry );
		Dictionary_Entry_Value_AddMember( particlePositionEntry, "x", Dictionary_Entry_Value_FromDouble( 0.7 ) );
		Dictionary_Entry_Value_AddMember( particlePositionEntry, "y", Dictionary_Entry_Value_FromDouble( 0.6 ) );
		Dictionary_Entry_Value_AddMember( particlePositionEntry, "z", Dictionary_Entry_Value_FromDouble( 0.5 ) );
	
		particlePositionEntry = Dictionary_Entry_Value_NewStruct();
		Dictionary_Entry_Value_AddElement( particlePositionsList, particlePositionEntry );
		Dictionary_Entry_Value_AddMember( particlePositionEntry, "x", Dictionary_Entry_Value_FromDouble( 0.8 ) );
		Dictionary_Entry_Value_AddMember( particlePositionEntry, "y", Dictionary_Entry_Value_FromDouble( 0.1 ) );
		Dictionary_Entry_Value_AddMember( particlePositionEntry, "z", Dictionary_Entry_Value_FromDouble( 0.3 ) );
	
		particlePositionEntry = Dictionary_Entry_Value_NewStruct();
		Dictionary_Entry_Value_AddElement( particlePositionsList, particlePositionEntry );
		Dictionary_Entry_Value_AddMember( particlePositionEntry, "x", Dictionary_Entry_Value_FromDouble( 0.9 ) );
		Dictionary_Entry_Value_AddMember( particlePositionEntry, "y", Dictionary_Entry_Value_FromDouble( 0.4 ) );
		Dictionary_Entry_Value_AddMember( particlePositionEntry, "z", Dictionary_Entry_Value_FromDouble( 0.1 ) );
	
		/* Init mesh */
		extensionMgr_Register = ExtensionManager_Register_New();
		mesh = ManualParticleLayoutSuite_BuildMesh( nDims, meshSize, minCrds, maxCrds, extensionMgr_Register );
	
		/* Configure the element-cell-layout */
		elementCellLayout = ElementCellLayout_New( "elementCellLayout", NULL, mesh );
	
		/* Build the mesh */
		Stg_Component_Build( mesh, 0, False );
		Stg_Component_Initialise( mesh, 0, False );
	
		/* Configure the gauss-particle-layout */
		particleLayout = ManualParticleLayout_New( "manualParticleLayout", NULL, GlobalCoordSystem, False,
         0, 0.0, dictionary );
	
		swarm = Swarm_New( "manualParticleSwarm", NULL, elementCellLayout, particleLayout, dim, sizeof(Particle),
			extensionMgr_Register, NULL, data->comm, NULL );
	
		/* Build the swarm */
		Stg_Component_Build( swarm, 0, False );
		Stg_Component_Initialise( swarm, 0, False );

		Journal_Enable_AllTypedStream( True );
		stream = Journal_Register( Info_Type, "ManualParticle" );

		/* Print out the particles on all cells */
		Stream_RedirectFile( stream, "testManualParticle.dat" );
		Swarm_PrintParticleCoords_ByCell( swarm, stream );
		Journal_Enable_AllTypedStream( False );

		pcu_filename_expected( "testManualParticleLayoutOutput.expected", expected_file );
		pcu_check_fileEq( "testManualParticle.dat", expected_file );

		Stg_Class_Delete( extensionMgr_Register );
		_Stg_Component_Delete( particleLayout );
		_Stg_Component_Delete( elementCellLayout );
		_Stg_Component_Delete( swarm );
		Stg_Class_Delete( dictionary );	
		remove( "testManualParticle.dat" );
	}
}

void ManualParticleLayoutSuite( pcu_suite_t* suite ) {
	pcu_suite_setData( suite, ManualParticleLayoutSuiteData );
	pcu_suite_setFixtures( suite, ManualParticleLayoutSuite_Setup, ManualParticleLayoutSuite_Teardown );
	pcu_suite_addTest( suite, ManualParticleLayoutSuite_TestManualParticle );
}


