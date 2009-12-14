#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "pcu/pcu.h"
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include "StgFEM/Discretisation/Discretisation.h"
#include "StgFEM/SLE/SystemSetup/SystemSetup.h"

typedef struct {
} StiffnessMatrixSuiteData;

void StiffnessMatrixSuite_Setup( StiffnessMatrixSuiteData* data ) { 
	Journal_Enable_AllTypedStream( False );
}

void StiffnessMatrixSuite_Teardown( StiffnessMatrixSuiteData* data ) {
	Journal_Enable_AllTypedStream( True );
}

FeVariable* buildFeVar() {
   CartesianGenerator*				gen;
   FeMesh*								feMesh;
   DofLayout*							dofs;
   FeEquationNumber*					eqNum;
   Variable_Register*				varReg;
   int									maxDecomp[3] = {0, 1, 1};
   int									sizes[3];
   double								minCrd[3];
   double								maxCrd[3];
   char									xml_input[PCU_PATH_MAX];
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
   var = Variable_NewVector( "velocity", NULL, Variable_DataType_Double, 2, (unsigned*)&arraySize, NULL, 
		(void**)arrayPtrs, varReg, "vx", "vy" );
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
   pcu_filename_input( "velWallVC.xml", xml_input );
   IO_Handler_ReadAllFromFile( ioHandler, xml_input, dict );
   bcs = (VariableCondition*)WallVC_New( "", NULL, "wallVC", varReg, cfReg, dict, feMesh );
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

      /* assign the position vector as the nodes vector */
      Variable_SetValue( var, n_i, pos );
   }

   /* Build and initialise system */
   Stg_Component_Build( bcs, 0, False );
   Stg_Component_Build( feVar, 0, False );
   Stg_Component_Initialise( feVar, 0, False );

   return feVar;
}


void StiffnessMatrixSuite_TestStiffnessMatrix( StiffnessMatrixSuiteData* data ) {
  FeVariable*				feVar;
  FeMesh*					mesh;
  int							nEls, nVerts, nDims;
  const int					*verts;
  double*					vert;
  double						val[3];
  InterpolationResult	ret;
  EntryPoint_Register*	ep_reg;
  StiffnessMatrix*		mat;
  ForceVector*				vec;
  IArray*					incArray = NULL;
  int							e_i, v_i, d_i;
  MPI_Comm					comm;

  pcu_docstring( "This test just creates a Stiffness matrix data structure, builds, itialises, refreshes and then destroys it.\n" );

  /* we need an EP register here for the construction of the ForceVector and StiffnessMatrix */
  ep_reg = EntryPoint_Register_New();

  /* here we build a feVar with the position vector stored on the nodes */
  feVar = buildFeVar();
  mesh = feVar->feMesh;
  FeVariable_SyncShadowValues( feVar );

  /* create SiffnessMatrix, it requires a ForceVector */
  comm = Comm_GetMPIComm( Mesh_GetCommTopology( mesh, MT_VERTEX ) );
  vec = ForceVector_New( "testVector", NULL, feVar, 2, ep_reg, comm );
  mat = StiffnessMatrix_New( "testMatrix", feVar, feVar, vec, NULL, 2, False, False, ep_reg, comm);

  /* build & initialise the mat, this should build and initialise the vec */
  Stg_Component_Build( mat, NULL, False );
  Stg_Component_Initialise( mat, NULL, False );

  StiffnessMatrix_RefreshMatrix( mat );

  pcu_check_true( mat->matrix != NULL );

  Stg_Component_Destroy( feVar, NULL, True );
  Stg_Component_Destroy( vec, NULL, True );
  Stg_Component_Destroy( mat, NULL, True );
  _Stg_Component_Delete( feVar );
  _Stg_Component_Delete( vec );
  /*_Stg_Component_Delete( mat );*/
  Stg_Class_Delete( ep_reg );
}

void StiffnessMatrixSuite( pcu_suite_t* suite ) {
   pcu_suite_setData( suite, StiffnessMatrixSuiteData );
   pcu_suite_setFixtures( suite, StiffnessMatrixSuite_Setup, StiffnessMatrixSuite_Teardown );
   pcu_suite_addTest( suite, StiffnessMatrixSuite_TestStiffnessMatrix );
}


