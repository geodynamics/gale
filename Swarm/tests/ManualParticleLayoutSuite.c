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
	unsigned							nDims;
	unsigned							meshSize[3];
	double							minCrds[3];
	double							maxCrds[3];
	ExtensionManager_Register*	extensionMgr_Register;
	Mesh*								mesh;
	GaussParticleLayout*			gaussParticleLayout;
	ElementCellLayout*			elementCellLayout;
	Swarm*							swarm;
	Dictionary*                dictionary;
	Dictionary_Entry_Value*    particlePositionsList;
	Dictionary_Entry_Value*    particlePositionEntry;
	ManualParticleLayout*		particleLayout;
	DomainContext*   				context;
	MPI_Comm       				comm;
   unsigned int   				rank;
   unsigned int   				nProcs;
} ManualParticleLayoutSuiteData;

Mesh* ManualParticleLayoutSuite_BuildMesh( unsigned nDims, unsigned* size, double* minCrds, double* maxCrds, ExtensionManager_Register* emReg ) {
	CartesianGenerator*	gen;
	Mesh*			mesh;

	gen = CartesianGenerator_New( "" );
	CartesianGenerator_SetDimSize( gen, nDims );
	CartesianGenerator_SetTopologyParams( gen, size, 0, NULL, NULL );
	CartesianGenerator_SetGeometryParams( gen, minCrds, maxCrds );

	mesh = Mesh_New( "" );
	Mesh_SetExtensionManagerRegister( mesh, emReg );
	Mesh_SetGenerator( mesh, gen );

	Stg_Component_Build( mesh, NULL, False );
	Stg_Component_Initialise( mesh, NULL, False );

	KillObject( mesh->generator );

	return mesh;
}

void ManualParticleLayoutSuite_Setup( ManualParticleLayoutSuiteData* data ) {
	Dimension_Index	dim;
	char					input_file[PCU_PATH_MAX];
	
	/* MPI Initializations */
	data->comm = MPI_COMM_WORLD;  
   MPI_Comm_rank( data->comm, &data->rank );
   MPI_Comm_size( data->comm, &data->nProcs );
   
   /* Read in the input xml file. */
	/*pcu_filename_input( "testManualParticleLayoutInput.xml", input_file );
	data->context = (DomainContext*)stgMainInitFromXML( input_file, data->comm );
	dictionary = data->context->dictionary;*/
   
	data->nDims = 3;
	data->meshSize[0] = 2;	data->meshSize[1] = 3;	data->meshSize[2] = 2;
	data->minCrds[0] = 0.0; data->minCrds[1] = 0.0; data->minCrds[2] = 0.0;
	data->maxCrds[0] = 300.0; data->maxCrds[1] = 12.0; data->maxCrds[2] = 300.0;
	
	/* Dictionary Initialization */
	data->dictionary = Dictionary_New();
	data->particlePositionsList = Dictionary_Entry_Value_NewList();
	Dictionary_Add( data->dictionary, "manualParticlePositions", data->particlePositionsList );
	
	data->particlePositionEntry = Dictionary_Entry_Value_NewStruct();
	Dictionary_Entry_Value_AddElement( data->particlePositionsList, data->particlePositionEntry );
	Dictionary_Entry_Value_AddMember( data->particlePositionEntry, "x", Dictionary_Entry_Value_FromDouble( 0.4 ) );
	Dictionary_Entry_Value_AddMember( data->particlePositionEntry, "y", Dictionary_Entry_Value_FromDouble( 0.3 ) );
	Dictionary_Entry_Value_AddMember( data->particlePositionEntry, "z", Dictionary_Entry_Value_FromDouble( 0.2 ) );
	
	data->particlePositionEntry = Dictionary_Entry_Value_NewStruct();
	Dictionary_Entry_Value_AddElement( data->particlePositionsList, data->particlePositionEntry );
	Dictionary_Entry_Value_AddMember( data->particlePositionEntry, "x", Dictionary_Entry_Value_FromDouble( 0.7 ) );
	Dictionary_Entry_Value_AddMember( data->particlePositionEntry, "y", Dictionary_Entry_Value_FromDouble( 0.6 ) );
	Dictionary_Entry_Value_AddMember( data->particlePositionEntry, "z", Dictionary_Entry_Value_FromDouble( 0.5 ) );
	
	data->particlePositionEntry = Dictionary_Entry_Value_NewStruct();
	Dictionary_Entry_Value_AddElement( data->particlePositionsList, data->particlePositionEntry );
	Dictionary_Entry_Value_AddMember( data->particlePositionEntry, "x", Dictionary_Entry_Value_FromDouble( 0.8 ) );
	Dictionary_Entry_Value_AddMember( data->particlePositionEntry, "y", Dictionary_Entry_Value_FromDouble( 0.1 ) );
	Dictionary_Entry_Value_AddMember( data->particlePositionEntry, "z", Dictionary_Entry_Value_FromDouble( 0.3 ) );
	
	data->particlePositionEntry = Dictionary_Entry_Value_NewStruct();
	Dictionary_Entry_Value_AddElement( data->particlePositionsList, data->particlePositionEntry );
	Dictionary_Entry_Value_AddMember( data->particlePositionEntry, "x", Dictionary_Entry_Value_FromDouble( 0.9 ) );
	Dictionary_Entry_Value_AddMember( data->particlePositionEntry, "y", Dictionary_Entry_Value_FromDouble( 0.4 ) );
	Dictionary_Entry_Value_AddMember( data->particlePositionEntry, "z", Dictionary_Entry_Value_FromDouble( 0.1 ) );
	
	/* Init mesh */
	data->extensionMgr_Register = ExtensionManager_Register_New();
	data->mesh = ManualParticleLayoutSuite_BuildMesh( data->nDims, data->meshSize, data->minCrds, data->maxCrds, data->extensionMgr_Register );
	
	/* Configure the element-cell-layout */
	data->elementCellLayout = ElementCellLayout_New( "elementCellLayout", data->mesh );
	
	/* Build the mesh */
	Stg_Component_Build( data->mesh, 0, False );
	Stg_Component_Initialise( data->mesh, 0, False );
	
	/* Configure the gauss-particle-layout */
	data->particleLayout = ManualParticleLayout_New( "manualParticleLayout", data->dictionary );
	
	data->swarm = Swarm_New( "testSwarm", data->elementCellLayout, data->particleLayout, dim, sizeof(Particle),
		data->extensionMgr_Register, NULL, data->comm, NULL );
	
	/* Build the swarm */
	Stg_Component_Build( data->swarm, 0, False );
	Stg_Component_Initialise( data->swarm, 0, False );
}

void ManualParticleLayoutSuite_Teardown( ManualParticleLayoutSuiteData* data ) {
	/* Destroy stuff */
	Stg_Class_Delete( data->particleLayout );
	Stg_Class_Delete( data->elementCellLayout );
	Stg_Class_Delete( data->swarm );
	Stg_Class_Delete( data->mesh );
	Stg_Class_Delete( data->extensionMgr_Register );
	Stg_Class_Delete( data->dictionary );	
}

void ManualParticleLayoutSuite_TestManualParticle( ManualParticleLayoutSuiteData* data ) {
	double 						x,y,z;
	unsigned int 				p, i, len;
	int							procToWatch;
	Stream*						stream;
	char 							expected_file[PCU_PATH_MAX];
	
	if( data->nProcs >= 2 ) {
		procToWatch = 1;
	}
	else {
		procToWatch = 0;
	}
	if( data->rank == procToWatch ) printf( "Watching rank: %i\n", data->rank );
	
	stream = Journal_Register( Info_Type, "ManualParticle" );
	
	if( data->rank == procToWatch ) {
		Stg_Class_Print( data->particleLayout, stream );
		/* Print out the particles on all cells */
		Stream_RedirectFile( stream, "manualParticle.dat" );
		Swarm_PrintParticleCoords_ByCell( data->swarm, stream );
	}
}

void ManualParticleLayoutSuite( pcu_suite_t* suite ) {
   pcu_suite_setData( suite, ManualParticleLayoutSuiteData );
   pcu_suite_setFixtures( suite, ManualParticleLayoutSuite_Setup, ManualParticleLayoutSuite_Teardown );
   pcu_suite_addTest( suite, ManualParticleLayoutSuite_TestManualParticle );
}
