/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003, Victorian Partnership for Advanced Computing (VPAC) Ltd, 
** 110 Victoria Street, Melbourne, 3053, Australia.
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
** $Id: testSolutionVector.c 964 2007-10-11 08:03:06Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include "StgFEM/Discretisation/Discretisation.h"
#include "StgFEM/SLE/SystemSetup/SystemSetup.h"

#include "StGermain/Base/Foundation/TestBegin.h"


FeVariable* buildFeVar() {
   CartesianGenerator* gen;
   FeMesh* feMesh;
   DofLayout* dofs;
   FeEquationNumber* eqNum;
   Variable_Register* varReg;
   int maxDecomp[3] = {0, 1, 1};
   int sizes[3];
   double minCrd[3];
   double maxCrd[3];
   static int arraySize;
   static double* arrayPtrs[2];
   int nRanks;
   Variable* var;
   VariableCondition* bcs;
   ConditionFunction_Register* cfReg;
   Dictionary* dict;
   XML_IO_Handler* ioHandler;
   FieldVariable_Register* fieldReg;
   FeVariable* feVar;
   int n_i;

   MPI_Comm_size( MPI_COMM_WORLD, &nRanks );
   sizes[0] = nRanks * 2;
   sizes[1] = sizes[2] = 2;
   minCrd[0] = minCrd[1] = minCrd[2] = 0.0;
   maxCrd[0] = maxCrd[1] = maxCrd[2] = (double)nRanks;

   gen = CartesianGenerator_New( "" );
   CartesianGenerator_SetDimSize( gen, 2 );
   CartesianGenerator_SetTopologyParams( gen, (unsigned*)sizes, 0, 
					 NULL, (unsigned*)maxDecomp );
   CartesianGenerator_SetGeometryParams( gen, minCrd, maxCrd );
   CartesianGenerator_SetShadowDepth( gen, 0 );

   feMesh = FeMesh_New( "" );
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

   dofs = DofLayout_New( "", varReg, 0, feMesh );
   dofs->nBaseVariables = 2;
   dofs->baseVariables = Memory_Alloc_Array_Unnamed( Variable*, 2 );
   dofs->baseVariables[0] = var->components[0];
   dofs->baseVariables[1] = var->components[1];
   Stg_Component_Build( dofs, NULL, False );
   Stg_Component_Initialise( dofs, NULL, False );

   ioHandler = XML_IO_Handler_New();
   dict = Dictionary_New();
   IO_Handler_ReadAllFromFile( ioHandler, "data/velWallVC.xml", dict );
   bcs = (VariableCondition*)WallVC_New( "", "wallVC", varReg, cfReg, 
					 dict, feMesh );
   Stg_Component_Build( bcs, NULL, False );
   Stg_Component_Initialise( bcs, NULL, False );

   eqNum = FeEquationNumber_New( "", feMesh, dofs, bcs, NULL );
   Stg_Component_Build( eqNum, NULL, False );
   Stg_Component_Initialise( eqNum, NULL, False );

   fieldReg = FieldVariable_Register_New();
   feVar = FeVariable_New( "velocity", feMesh, NULL, dofs, bcs, NULL, NULL, 2, True, 
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


void testSetup( int* argc, char** argv[] ) {
   StGermain_Init( argc, argv );
   StgDomain_Init( argc, argv );
   StgFEM_Discretisation_Init( argc, argv );
   //StgFEM_SLE_LinearAlgebra_Init( argc, argv );
   StgFEM_SLE_SystemSetup_Init( argc, argv );
   Stream_Enable( Journal_GetTypedStream( Debug_Type ), False );
   Stream_Enable( Journal_GetTypedStream( Info_Type ), False );
   Stream_Enable( Journal_GetTypedStream( Dump_Type ), False );
   Stream_Enable( Journal_GetTypedStream( Error_Type ), False );
}

void testTeardown() {
   StgFEM_SLE_SystemSetup_Finalise();
   //StgFEM_SLE_LinearAlgebra_Finalise();
   StgFEM_Discretisation_Finalise();
   StgDomain_Finalise();
   StGermain_Finalise();
}

TestBegin( All ) {
   FeVariable* feVar;
   FeMesh* mesh;
   int nEls, nVerts, nDims;
   const int *verts;
   double* vert;
   double val[3];
   InterpolationResult ret;
   SolutionVector* sol;
   int lSize;
   double* array;
   IArray* incArray;
   int e_i, v_i, d_i, a_i;

   feVar = buildFeVar();
   FeVariable_SyncShadowValues( feVar );
   sol = SolutionVector_New( "velocity", MPI_COMM_WORLD, feVar );
   TestTrue( sol );
   Stg_Component_Build( sol, NULL, False );

   SolutionVector_LoadCurrentFeVariableValuesOntoVector( sol );
   lSize = Vector_GetLocalSize( sol->vector );
   Vector_GetArray( sol->vector, &array );
   for( a_i = 0; a_i < lSize; a_i++ )
      array[a_i] += 1.0;
   Vector_RestoreArray( sol->vector, &array );
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
   TestTrue( e_i == nEls );

   NewClass_Delete( incArray );

  done:
   FreeObject( feVar );
   FreeObject( sol );
}
TestEnd


#define nTests 1
TestSuite_Test tests[nTests] = {{"all", testAll}};


#include "StGermain/Base/Foundation/TestEnd.h"
