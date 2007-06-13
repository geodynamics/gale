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
** $Id: testCompositeVC.c 4144 2007-06-13 04:33:59Z LukeHodkinson $
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


void quadratic(Index index, Variable_Index var_I, void* context, void* result)
{
	*(double *)result = 20.0;
}


void exponential(Index index, Variable_Index var_I, void* context, void* result)
{
	*(double *)result = 30.0;
}


Mesh* buildMesh( unsigned nDims, unsigned* size, 
		     double* minCrds, double* maxCrds, 
		     ExtensionManager_Register* emReg )
{
	CartesianGenerator*	gen;
	Mesh*			mesh;
	unsigned		maxDecomp[3] = {0, 1, 1};

	gen = CartesianGenerator_New( "" );
	gen->shadowDepth = 0;
	CartesianGenerator_SetDimSize( gen, nDims );
	CartesianGenerator_SetTopologyParams( gen, size, 0, NULL, maxDecomp );
	CartesianGenerator_SetGeometryParams( gen, minCrds, maxCrds );

	mesh = Mesh_New( "" );
	Mesh_SetExtensionManagerRegister( mesh, emReg );
	Mesh_SetGenerator( mesh, gen );

	Stg_Component_Build( mesh, NULL, False );
	Stg_Component_Initialise( mesh, NULL, False );

	KillObject( mesh->generator );

	return mesh;
}


int main(int argc, char *argv[])
{
	MPI_Comm		CommWorld;
	int		rank;
	int		procCount;
	int		procToWatch;
	Stream*		stream;
	
	Dictionary*		dictionary;
	XML_IO_Handler*		io_handler;
	
	unsigned	nDims = 3;
	unsigned	meshSize[3] = {3, 3, 3};
	double		minCrds[3] = {0.0, 0.0, 0.0};
	double		maxCrds[3] = {1.0, 1.0, 1.0};
	Mesh*		mesh;
	
	Variable*			var[7];
	Variable_Register*		variable_Register;
	WallVC*				vc;
	CompositeVC*			cvc;
	ConditionFunction*		quadCF;
	ConditionFunction*		expCF;
	ConditionFunction_Register*	conFunc_Register;
	ExtensionManager_Register*		extensionMgr_Register;
	
	double*		array[7];
	char*		vcKey[] = {"WallVC_Front", "WallVC_Back", "WallVC_Left", "WallVC_Right",
				"WallVC_Top", "WallVC_Bottom"};
	char*		vcKeyName[] = {"WallVC_FrontName", "WallVC_BackName", "WallVC_LeftName", "WallVC_RightName",
				"WallVC_TopName", "WallVC_BottomName"};
	char*		varName[] = {"x", "y", "z", "vx", "vy", "vz", "temp"};

	unsigned	nDomains;
	
	Index	i, j, k;
	
	/* Initialise MPI, get world info */
	MPI_Init(&argc, &argv);
	MPI_Comm_dup( MPI_COMM_WORLD, &CommWorld );
	MPI_Comm_size(CommWorld, &procCount);
	MPI_Comm_rank(CommWorld, &rank);
	
	Base_Init( &argc, &argv );
	
	DiscretisationGeometry_Init( &argc, &argv );
	DiscretisationShape_Init( &argc, &argv );
	DiscretisationMesh_Init( &argc, &argv );
	DiscretisationUtils_Init( &argc, &argv );
	MPI_Barrier( CommWorld ); /* Ensures copyright info always come first in output */

	Journal_Enable_NamedStream( Info_Type, CartesianGenerator_Type, False );

	io_handler = XML_IO_Handler_New();
	
	stream = Journal_Register (Info_Type, "myStream");

	procToWatch = argc >= 2 ? atoi(argv[1]) : 0;
	
	dictionary = Dictionary_New();
	Dictionary_Add(dictionary, "outputPath", Dictionary_Entry_Value_FromString("./output"));
	IO_Handler_ReadAllFromFile(io_handler, "data/wallVC.xml", dictionary);

	extensionMgr_Register = ExtensionManager_Register_New();
	
	/* Create a mesh. */
	mesh = buildMesh( nDims, meshSize, minCrds, maxCrds, extensionMgr_Register );
	nDomains = Mesh_GetDomainSize( mesh, MT_VERTEX );
	
	/* Create CF stuff */
	quadCF = ConditionFunction_New(quadratic, "quadratic");
	expCF = ConditionFunction_New(exponential, "exponential");
	conFunc_Register = ConditionFunction_Register_New();
	ConditionFunction_Register_Add(conFunc_Register, quadCF);
	ConditionFunction_Register_Add(conFunc_Register, expCF);
	
	/* Create variable register */
	variable_Register = Variable_Register_New();
	
	/* Create variables */
	for (i = 0; i < 6; i++) {
		array[i] = Memory_Alloc_Array( double, nDomains, "array[i]" );
		var[i] = Variable_NewScalar( varName[i], Variable_DataType_Double, &nDomains, (void**)&array[i], 0 ); 
		Variable_Register_Add(variable_Register, var[i]);
	}
	array[6] = Memory_Alloc_Array( double, nDomains*5, "array[6]" );
	var[6] = Variable_NewVector( varName[6], Variable_DataType_Double, 5, &nDomains, (void**)&array[6], 0 );
	Variable_Register_Add(variable_Register, var[6]);
	Variable_Register_BuildAll(variable_Register);
	
	/* Create CompositeVC */
	cvc = CompositeVC_New( "CompositeVC", variable_Register, conFunc_Register, dictionary, mesh );
	
	for (i = 0; i < 6; i++)
	{
		vc = WallVC_New( vcKeyName[i], vcKey[i], variable_Register, conFunc_Register, dictionary, mesh );
		Stg_Component_Build( vc, 0, False );
		CompositeVC_Add(cvc, vc, True);
	}
	
	Stg_Component_Build( cvc, 0, False );
	
	for (j = 0; j < 6; j++)
		memset(array[j], 0, sizeof(double)*nDomains);
	memset(array[6], 0, sizeof(double)*nDomains*5);
	VariableCondition_Apply(cvc, NULL);
	
	if (rank == procToWatch)
	{
		Stg_Class_Print(cvc, stream);
		printf("\n");
		for (j = 0; j < 6; j++)
		{
			printf("\nvar[%u]: %.2lf", j, array[j][0]);
			for (k = 1; k < nDomains; k++)
				printf(", %.2lf", array[j][k]);
		}
		printf("\nvar[7]: %.2lf", array[6][0]);
		for (j = 1; j < nDomains*5; j++)
			printf(", %.2lf", array[6][j]);
		printf("\n\n");
			
		for (j = 0; j < 7; j++)
		{
			for (k = 0; k < nDomains; k++)
				printf("%s ", VariableCondition_IsCondition(cvc, k, j) ? "True " : "False");
			printf("\n");
		}
		printf("\n");
			
		for (j = 0; j < 7; j++)
		{
			for (k = 0; k < nDomains; k++)
			{
				VariableCondition_ValueIndex	valIndex;
				
				valIndex = VariableCondition_GetValueIndex(cvc, k, j);
				if (valIndex != (unsigned)-1)
					printf("%03u ", valIndex);
				else
					printf("XXX ");
			}
			printf("\n");
		}
		printf("\n");
	}

	Stg_Class_Delete(cvc);
	Stg_Class_Delete(variable_Register);
	for (i = 0; i < 7; i++)
	{
		Stg_Class_Delete(var[i]);
		if (array[i]) Memory_Free(array[i]);
	}
	Stg_Class_Delete(conFunc_Register);
	Stg_Class_Delete(quadCF);
	Stg_Class_Delete(expCF);
	Stg_Class_Delete(dictionary);
	FreeObject( mesh );
	
	DiscretisationUtils_Finalise();
	DiscretisationMesh_Finalise();
	DiscretisationShape_Finalise();
	DiscretisationGeometry_Finalise();
	
	Base_Finalise();
	
	/* Close off MPI */
	MPI_Finalize();
	
	return 0; /* success */
}
