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
**   Tests the CornerVCSuite
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

#include "CornerVCSuite.h"

typedef struct {
	MPI_Comm	comm;
	unsigned rank;
	unsigned nProcs;
} CornerVCSuiteData;

void CornerVCSuite_quadratic(Index index, Variable_Index var_I, void* context, void* result) {
	*(double *)result = 20.0;
}

void CornerVCSuite_exponential(Index index, Variable_Index var_I, void* context, void* result) {
	*(double *)result = 30.0;
}


Mesh* CornerVCSuite_buildMesh( unsigned nDims, unsigned* size, double* minCrds, double* maxCrds, ExtensionManager_Register* emReg ) {
	CartesianGenerator*	gen;
	Mesh*						mesh;

	gen = CartesianGenerator_New( "", NULL );
	gen->shadowDepth = 0;
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

void CornerVCSuite_Setup( CornerVCSuiteData* data ) {
	/* MPI Initializations */
	data->comm = MPI_COMM_WORLD;
	MPI_Comm_rank( data->comm, &data->rank );
	MPI_Comm_size( data->comm, &data->nProcs );
}

void CornerVCSuite_Teardown( CornerVCSuiteData* data ) {
}

void CornerVCSuite_TestCornerVC( CornerVCSuiteData* data ) {
	unsigned								nDomains;
	unsigned								nDims = 3;
	unsigned								meshSize[3] = {3, 3, 3};
	int									procToWatch;
	double								minCrds[3] = {0.0, 0.0, 0.0};
	double								maxCrds[3] = {1.0, 1.0, 1.0};
	double*								array[7];
	char*									vcKey[] = { "CornerVC_BottomLeftFront", "CornerVC_BottomRightFront",
															"CornerVC_TopLeftFront", "CornerVC_TopRightFront",
															"CornerVC_BottomLeftBack", "CornerVC_BottomRightBack",
															"CornerVC_TopLeftBack", "CornerVC_TopRightBack" };
	char*									vcKeyName[] = { "CornerVC_BottomLeftFrontName", "CornerVC_BottomRightFrontName",
															"CornerVC_TopLeftFrontName", "CornerVC_TopRightFrontName",
															"CornerVC_BottomLeftBackName", "CornerVC_BottomRightBackName",
															"CornerVC_TopLeftBackName", "CornerVC_TopRightBackName" };
	char*									varName[] = {"x", "y", "z", "vx", "vy", "vz", "temp"};
	char									input_file[PCU_PATH_MAX];
	char									expected_file[PCU_PATH_MAX];
	Mesh*									mesh;
	Variable_Register*				variable_Register;
	ConditionFunction*				quadCF;
	ConditionFunction_Register*	conFunc_Register;
	ExtensionManager_Register*		extensionMgr_Register;
	DomainContext*						context;
	Dictionary*							dictionary;
	Stream*								stream;
	XML_IO_Handler*					io_handler;
	Variable*							var[9];
	VariableCondition*				vc; 
	Index									i;

	procToWatch = data->nProcs >=2 ? 1 : 0;

   io_handler = XML_IO_Handler_New();

 	stream = Journal_Register( Info_Type, "CornerVCStream" );
   Stream_RedirectFile( stream, "testCornerVC.dat" );

   dictionary = Dictionary_New();
   Dictionary_Add(dictionary, "outputPath", Dictionary_Entry_Value_FromString("./output"));

	/* Input file */
	pcu_filename_input( "cornerVC.xml", input_file );
   IO_Handler_ReadAllFromFile(io_handler, input_file, dictionary);
   fflush(stdout);

   extensionMgr_Register = ExtensionManager_Register_New(); 

 	/* Create a mesh. */
   mesh = (Mesh*) CornerVCSuite_buildMesh( nDims, meshSize, minCrds, maxCrds, extensionMgr_Register );
   nDomains = Mesh_GetDomainSize( mesh, MT_VERTEX );

   /* Create CF stuff */
   conFunc_Register = ConditionFunction_Register_New();

	/* Create variable register */
   variable_Register = Variable_Register_New();

	/* Create variables */
	for (i = 0; i < 7; i++) {
		array[i] = Memory_Alloc_Array( double, nDomains, "array[i]" );
      var[i] = Variable_NewScalar( varName[i], Variable_DataType_Double, &nDomains, NULL, (void**)&array[i], 0 );
		Variable_Register_Add(variable_Register, var[i]);
	}

	Variable_Register_BuildAll(variable_Register);
   /* Create CornerVC */
	for (i = 0; i < 8; i++) {
 		Index j, k;

		vc = (VariableCondition*) CornerVC_New( vcKeyName[i], vcKey[i], variable_Register, conFunc_Register, dictionary, mesh );
		_CornerVC_ReadDictionary(vc, dictionary);
		Stg_Component_Build( vc, 0, False );
		for (j = 0; j < 7; j++) { 
			memset(array[j], 0, sizeof(double)* nDomains );
		}
		VariableCondition_Apply(vc, NULL);

		if (data->rank == procToWatch) {
			Journal_Printf( stream, "Testing for %s\n", vcKey[i] );
			for (j = 0; j < 7; j++) {
				Journal_Printf( stream,"\nvar[%u]: %.2lf", j, array[j][0]);
				for (k = 1; k < nDomains; k++)
					Journal_Printf( stream,", %.2lf", array[j][k]);
			} Journal_Printf( stream,"\n\n");

			for (j = 0; j < 6; j++) {
				for (k = 0; k < nDomains; k++)
					Journal_Printf( stream,"%s ", VariableCondition_IsCondition(vc, k, j) ? "True " : "False");
				Journal_Printf( stream,"\n");
         } Journal_Printf( stream,"\n");

			for (j = 0; j < 6; j++) {
				for (k = 0; k < nDomains; k++) {
					VariableCondition_ValueIndex  valIndex;
					valIndex = VariableCondition_GetValueIndex(vc, k, j);
					if (valIndex != (unsigned)-1)
						Journal_Printf( stream,"%03u ", valIndex);
					else
						Journal_Printf( stream,"XXX ");
				} Journal_Printf( stream,"\n");
			} Journal_Printf( stream,"\n");
	
      }
		Stg_Class_Delete(vc);
	}
	if (data->rank == procToWatch) {
		pcu_filename_expected( "testCornerVC.expected", expected_file );
		pcu_check_fileEq( "testCornerVC.dat", expected_file );
		remove( "testCornerVC.dat" );
	}

	Stg_Class_Delete(variable_Register);
	for (i = 0; i < 7; i++) {
		Stg_Class_Delete(var[i]);
		if (array[i]) Memory_Free(array[i]);
	}
	Stg_Class_Delete(conFunc_Register);
	Stg_Class_Delete(dictionary);
	FreeObject( mesh );
}

void CornerVCSuite( pcu_suite_t* suite ) {
	pcu_suite_setData( suite, CornerVCSuiteData );
   pcu_suite_setFixtures( suite, CornerVCSuite_Setup, CornerVCSuite_Teardown );
   pcu_suite_addTest( suite, CornerVCSuite_TestCornerVC );
}


