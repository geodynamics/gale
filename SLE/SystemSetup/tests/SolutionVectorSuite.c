#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pcu/pcu.h"
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include "StgFEM/Discretisation/Discretisation.h"
#include "StgFEM/SLE/SystemSetup/SystemSetup.h"

typedef struct {
} SolutionVectorSuiteData;

void SolutionVectorSuite_Setup( SolutionVectorSuiteData* data ) { }

void SolutionVectorSuite_Teardown( SolutionVectorSuiteData* data ) {}


FeVariable* SolutionVectorSuite_buildFeVar() {
   CartesianGenerator*				gen;
   FeMesh*								feMesh;
   DofLayout*							dofs;
   FeEquationNumber*					eqNum;
   Variable_Register*				varReg;
   int									maxDecomp[3] = {0, 1, 1};
   int									sizes[3];
   double								minCrd[3];
   double								maxCrd[3];
   static int							arraySize;
   static double*						arrayPtrs[2];
   int									nRanks;
   Variable*							var;
   VariableCondition*				bcs;
   ConditionFunction_Register*	cfReg;
   Dictionary*							dict;
   XML_IO_Handler*					ioHandler;
   FieldVariable_Register*			fieldReg;
   FeVariable*							feVar;
   int									n_i;
	char									xml_input[PCU_PATH_MAX];

	pcu_filename_input( "testSolutionVector.xml", xml_input );

   MPI_Comm_size( MPI_COMM_WORLD, &nRanks );
   sizes[0] = nRanks * 2;
   sizes[1] = sizes[2] = 2;
   minCrd[0] = minCrd[1] = minCrd[2] = 0.0;
   maxCrd[0] = maxCrd[1] = maxCrd[2] = (double)nRanks;

   gen = CartesianGenerator_New( "", NULL );
   CartesianGenerator_SetDimSize( gen, 2 );
   CartesianGenerator_SetTopologyParams( gen, (unsigned*)sizes, 0, NULL, (unsigned*)maxDecomp );
   CartesianGenerator_SetGeometryParams( gen, minCrd, maxCrd );
   CartesianGenerator_SetShadowDepth( gen, 0 );

   feMesh = FeMesh_New( "", NULL );
   Mesh_SetGenerator( feMesh, gen );
   FeMesh_SetElementFamily( feMesh, "linear" );
   Stg_Component_Build( feMesh, NULL, False );

   varReg = Variable_Register_New();
   cfReg = ConditionFunction_Register_New();

   arraySize = Mesh_GetDomainSize( feMesh, MT_VERTEX );
   arrayPtrs[0] = Memory_Alloc_Array_Unnamed( double, arraySize * 2 );
/*
   arrayPtrs[1] = Memory_Alloc_Array_Unnamed( double, arraySize );
*/
   var = Variable_NewVector( "velocity", Variable_DataType_Double, 2, (unsigned*)&arraySize, NULL, 
		(void**)arrayPtrs, varReg, 
		"vx", "vy" );
   Variable_Register_BuildAll( varReg );

   dofs = DofLayout_New( "", NULL, varReg, 0, feMesh );
   dofs->nBaseVariables = 2;
   dofs->baseVariables = Memory_Alloc_Array_Unnamed( Variable*, 2 );
   dofs->baseVariables[0] = var->components[0];
   dofs->baseVariables[1] = var->components[1];
   Stg_Component_Build( dofs, NULL, False );
   Stg_Component_Initialise( dofs, NULL, False );

   ioHandler = XML_IO_Handler_New();
   dict = Dictionary_New();
   IO_Handler_ReadAllFromFile( ioHandler, xml_input, dict );
   bcs = (VariableCondition*)WallVC_New( "", "wallVC", varReg, cfReg, dict, feMesh );
   Stg_Component_Build( bcs, NULL, False );
   Stg_Component_Initialise( bcs, NULL, False );

   eqNum = FeEquationNumber_New( "", NULL, feMesh, dofs, bcs, NULL );
   Stg_Component_Build( eqNum, NULL, False );
   Stg_Component_Initialise( eqNum, NULL, False );

   fieldReg = FieldVariable_Register_New();
   feVar = FeVariable_New( "velocity", NULL, feMesh, NULL, dofs, bcs, NULL, NULL, 2, True, 
		False, False, fieldReg );

   for( n_i = 0; n_i < Mesh_GetLocalSize( feMesh, 0 ); n_i++ ) {
      /*const double pi=acos(-1.0);*/
      double* pos = Mesh_GetVertex( feMesh, n_i );

      Variable_SetValue( var, n_i, pos );
		/*
      arrayPtrs[0][n_i] = (1.0 - pos[1]) + 0.1 * cos( pi * pos[0] ) * sin(  pi * pos[1]  );
      arrayPtrs[1][n_i] = (1.0 - pos[1]) + 0.1 * cos( pi * pos[0] ) * sin(  pi * pos[1]  );
		*/
   }

   /* Build and initialise system */
   Stg_Component_Build( bcs, 0, False );
   Stg_Component_Build( feVar, 0, False );
   Stg_Component_Initialise( feVar, 0, False );

   return feVar;
}


void SolutionVectorSuite_TestSolutionVector( SolutionVectorSuiteData* data ) {
	FeVariable*				feVar;
	FeMesh*					mesh;
	int						nEls, nVerts, nDims;
	const int				*verts;
	double*					vert;
	double					val[3];
	InterpolationResult	ret;
	SolutionVector*		sol;
	int						lSize;
	double*					array;
	IArray*					incArray;
	int						e_i, v_i, d_i, a_i;

	feVar = SolutionVectorSuite_buildFeVar();
	FeVariable_SyncShadowValues( feVar );
	sol = SolutionVector_New( "velocity", NULL, MPI_COMM_WORLD, feVar );
	/* Check solution vector created */
	pcu_check_true(sol);
	Stg_Component_Build( sol, NULL, False );

	SolutionVector_LoadCurrentFeVariableValuesOntoVector( sol );
	VecGetLocalSize( sol->vector, &lSize );
	VecGetArray( sol->vector, &array );

	for( a_i = 0; a_i < lSize; a_i++ )
		array[a_i] += 1.0;
	VecRestoreArray( sol->vector, &array );
	SolutionVector_UpdateSolutionOntoNodes( sol );

	mesh = feVar->feMesh;
	nDims = Mesh_GetDimSize( mesh );
	nEls = Mesh_GetDomainSize( mesh, nDims );
	incArray = IArray_New();

	for( e_i = 0; e_i < nEls; e_i++ ) {
		Mesh_GetIncidence( mesh, nDims, e_i, 0, incArray );
		nVerts = IArray_GetSize( incArray );
		verts = IArray_GetPtr( incArray );
		for( v_i = 0; v_i < nVerts; v_i++ ) {
			vert = Mesh_GetVertex( mesh, verts[v_i] );
			ret = FieldVariable_InterpolateValueAt( feVar, vert, val );
			if( ret != LOCAL && ret != SHADOW )
				continue;
			for( d_i = 0; d_i < nDims; d_i++ ) {
				if( !Num_Approx( vert[d_i] + 1.0, val[d_i] ) )
					break;
			}
			if( d_i < nDims )
				break;
		}
		if( v_i < nVerts )
			break;
	}

	/* Check all elements processed */
	pcu_check_true(e_i == nEls);

	NewClass_Delete( incArray );

	FreeObject( feVar );
	FreeObject( sol );
}

void SolutionVectorSuite( pcu_suite_t* suite ) {
   pcu_suite_setData( suite, SolutionVectorSuiteData );
   pcu_suite_setFixtures( suite, SolutionVectorSuite_Setup, SolutionVectorSuite_Teardown );
   pcu_suite_addTest( suite, SolutionVectorSuite_TestSolutionVector );
}


