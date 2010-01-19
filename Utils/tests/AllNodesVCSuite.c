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
**   Tests the AllNodesVCSuite
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

#include "AllNodesVCSuite.h"

typedef struct {
	MPI_Comm	comm;
	int		rank;
	int		nProcs;
} AllNodesVCSuiteData;

void AllNodesVCSuite_quadratic(Index index, Variable_Index var_I, void* context, void* result) {
	*(double *)result = 20.0;
}

Mesh* AllNodesVCSuite_buildMesh( unsigned nDims, unsigned* size, double* minCrds, double* maxCrds, ExtensionManager_Register* emReg ) {
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

void AllNodesVCSuite_Setup( AllNodesVCSuiteData* data ) {
	/* MPI Initializations */
	data->comm = MPI_COMM_WORLD;
	MPI_Comm_rank( data->comm, &data->rank );
	MPI_Comm_size( data->comm, &data->nProcs );
}

void AllNodesVCSuite_Teardown( AllNodesVCSuiteData* data ) {
}

/* Prototype because I can not find the prototype in any headers. */
Mesh* buildMesh( unsigned nDims, unsigned* size, double* minCrds, double* maxCrds, ExtensionManager_Register* emReg );

void AllNodesVCSuite_TestAllNodesVC( AllNodesVCSuiteData* data ) {
	unsigned								nDomains;
	unsigned								nDims = 3;
	unsigned								meshSize[3] = {3, 3, 3};
	int									procToWatch;
	double								minCrds[3] = {0.0, 0.0, 0.0};
	double								maxCrds[3] = {1.0, 1.0, 1.0};
	double*								array[7];
	char*									vcKey = "AllNodesVC";
	char*									varName[] = {"x", "y", "z", "vx", "vy", "vz", "temp"};
	char									input_file[PCU_PATH_MAX];
	char									expected_file[PCU_PATH_MAX];
	Mesh*									mesh;
	Variable_Register*				variable_Register;
	ConditionFunction*				quadCF;
	ConditionFunction_Register*	conFunc_Register;
	ExtensionManager_Register*		extensionMgr_Register;
	Dictionary*							dictionary;
	Stream*								stream;
	XML_IO_Handler*					io_handler;
	Variable*							var[7];
	VariableCondition*				vc; 
	Index									i, j, k;

	procToWatch = data->nProcs >=2 ? 1 : 0;

   io_handler = XML_IO_Handler_New();

 	stream = Journal_Register( Info_Type, (Name)"AllNodesVCStream"  );
   Stream_RedirectFile( stream, "testAllNodesVC.dat" );

   dictionary = Dictionary_New();
   Dictionary_Add( dictionary, (Dictionary_Entry_Key)"outputPath", Dictionary_Entry_Value_FromString("./output") );

	/* Input file */
	pcu_filename_input( "allVC.xml", input_file );
   IO_Handler_ReadAllFromFile(io_handler, input_file, dictionary);
   fflush(stdout);

   extensionMgr_Register = ExtensionManager_Register_New(); 

 	/* Create a mesh. */
   mesh = (Mesh*) AllNodesVCSuite_buildMesh( nDims, meshSize, minCrds, maxCrds, extensionMgr_Register );
   nDomains = Mesh_GetDomainSize( mesh, MT_VERTEX );

   /* Create CF stuff */
   quadCF = ConditionFunction_New( AllNodesVCSuite_quadratic, (Name)"quadratic");
   conFunc_Register = ConditionFunction_Register_New( );
   ConditionFunction_Register_Add(conFunc_Register, quadCF);

	/* Create variable register */
   variable_Register = Variable_Register_New();

   /* Create variables */
   for (i = 0; i < 6; i++) {
		array[i] = Memory_Alloc_Array( double, nDomains, "array[i]" );
      var[i] = Variable_NewScalar( varName[i], NULL, Variable_DataType_Double, (Index*)&nDomains, NULL, (void**)&array[i], 0  );
      Variable_Register_Add(variable_Register, var[i]);
   }
   array[6] = Memory_Alloc_Array( double, nDomains*5, "array[6]" );
   var[6] = Variable_NewVector( varName[6], NULL, Variable_DataType_Double, 5, &nDomains, NULL, (void**)&array[6], 0 );
   Variable_Register_Add(variable_Register, var[6]);
   Variable_Register_BuildAll(variable_Register);

   /* Create AllVC */
   vc = (VariableCondition*)AllNodesVC_New( "AllNodesVC", NULL, vcKey, variable_Register, conFunc_Register, dictionary, mesh );
   Stg_Component_Build( vc, 0, False );

   for (j = 0; j < 6; j++)
      memset(array[j], 0, sizeof(double)*nDomains);
   memset(array[6], 0, sizeof(double)*nDomains*5);
   VariableCondition_Apply(vc, NULL);

	if (data->rank == procToWatch) {
		Journal_Printf( stream,"Testing for %s\n", vcKey);
		for (j = 0; j < 6; j++) {
			Journal_Printf( stream, "\nvar[%u]: %.2lf", j, array[j][0]) ;
			for (k = 1; k < nDomains; k++)
				Journal_Printf( stream, ", %.2lf", array[j][k] );
      }
		Journal_Printf( stream, "\nvar[6]: %.2lf", array[6][0] );
		for (j = 1; j < nDomains*5; j++)
			Journal_Printf( stream, ", %.2lf", array[6][j] );
		Journal_Printf( stream, "\n\n" );

		for (j = 0; j < 7; j++) {
			for (k = 0; k < nDomains; k++)
				Journal_Printf( stream, "%s ", VariableCondition_IsCondition(vc, k, j) ? "True " : "False" );
			Journal_Printf( stream, "\n" );
		} Journal_Printf( stream, "\n" );

		for (j = 0; j < 7; j++) {
			for (k = 0; k < nDomains; k++) {
				VariableCondition_ValueIndex  valIndex;

            valIndex = VariableCondition_GetValueIndex(vc, k, j);
            if (valIndex != (unsigned)-1)
					Journal_Printf( stream, "%03u ", valIndex);
            else
					Journal_Printf( stream, "XXX ");
			} Journal_Printf( stream, "\n" );
		} Journal_Printf( stream, "\n" );

		pcu_filename_expected( "testAllNodesVC.expected", expected_file );
		pcu_check_fileEq( "testAllNodesVC.dat", expected_file );
		remove( "testAllNodesVC.dat" );
	}

   Stg_Class_Delete(vc);
   Stg_Class_Delete(variable_Register);

   for (i = 0; i < 7; i++) {
      Stg_Class_Delete(var[i]);
      if (array[i]) Memory_Free(array[i]);
   }

	Stg_Class_Delete(extensionMgr_Register);
	Stg_Class_Delete(io_handler);
   Stg_Class_Delete(conFunc_Register);
   Stg_Class_Delete(dictionary);
   FreeObject( mesh );
}

void AllNodesVCSuite( pcu_suite_t* suite ) {
	pcu_suite_setData( suite, AllNodesVCSuiteData );
   pcu_suite_setFixtures( suite, AllNodesVCSuite_Setup, AllNodesVCSuite_Teardown );
   pcu_suite_addTest( suite, AllNodesVCSuite_TestAllNodesVC );
}


