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
**   Tests the FeEquationNumberSuite
**
** $Id: testFeEquationNumber.c 3462 2006-02-19 06:53:24Z WalterLandry $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>

#include "pcu/pcu.h"
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/Discretisation/Discretisation.h>
#include "FeEquationNumberSuite.h"

typedef struct {
} FeEquationNumberSuiteData;

FeEquationNumber* buildEqNum() {
   CartesianGenerator*	gen;
   FeMesh*					feMesh;
   DofLayout*				dofs;
   FeEquationNumber*		eqNum;
   Variable_Register*	varReg;
   Variable*				vars[2];
   int						maxDecomp[3] = {0, 1, 1};
   int						sizes[3];
   double					minCrd[3];
   double					maxCrd[3];
   SizeT						dataOffs = 1;
   Variable_DataType		dataType = Variable_DataType_Double;
   int						nDataTypes = 1;
   char*						dataNames = "nothing";
   static SizeT			structSize = sizeof(double);
   static int				arraySize;
   static void*			arrayPtrs[2];
   int						nRanks;

   MPI_Comm_size( MPI_COMM_WORLD, &nRanks );
   sizes[0] = nRanks * 2;
   sizes[1] = sizes[2] = 2;
   minCrd[0] = minCrd[1] = minCrd[2] = 0.0;
   maxCrd[0] = maxCrd[1] = maxCrd[2] = (double)nRanks;

   gen = CartesianGenerator_New( "", NULL );
   CartesianGenerator_SetDimSize( gen, 3 );
   CartesianGenerator_SetTopologyParams( gen, (unsigned*)sizes, 0, NULL, (unsigned*)maxDecomp );
   CartesianGenerator_SetGeometryParams( gen, minCrd, maxCrd );
   CartesianGenerator_SetShadowDepth( gen, 0 );

   feMesh = FeMesh_New( "", NULL );
   Mesh_SetGenerator( feMesh, gen );
   FeMesh_SetElementFamily( feMesh, "linear" );
   Stg_Component_Build( feMesh, NULL, False );

   varReg = Variable_Register_New();

   arraySize = Mesh_GetDomainSize( feMesh, MT_VERTEX );
   arrayPtrs[0] = Memory_Alloc_Array_Unnamed( double, arraySize );
   arrayPtrs[1] = Memory_Alloc_Array_Unnamed( double, arraySize );
   vars[0] = Variable_New( "one", 1, &dataOffs, &dataType, (unsigned*)&nDataTypes, 
		&dataNames, &structSize, (unsigned*)&arraySize, NULL, arrayPtrs, varReg );
   vars[1] = Variable_New( "two", 1, &dataOffs, &dataType, (unsigned*)&nDataTypes, 
		&dataNames, &structSize, (unsigned*)&arraySize, NULL, arrayPtrs + 1, varReg );

   dofs = DofLayout_New( "", NULL, varReg, 0, feMesh );
   dofs->nBaseVariables = 2;
   dofs->baseVariables = Memory_Alloc_Array_Unnamed( Variable*, 2 );
   dofs->baseVariables[0] = vars[0];
   dofs->baseVariables[1] = vars[1];
   Stg_Component_Build( dofs, NULL, False );
   Stg_Component_Initialise( dofs, NULL, False );

   eqNum = FeEquationNumber_New( "", feMesh, dofs, NULL, NULL );
   Stg_Component_Build( eqNum, NULL, False );
   Stg_Component_Initialise( eqNum, NULL, False );

   return eqNum;
}

FeEquationNumber* buildEqNumBCs() {
   CartesianGenerator* gen;
   FeMesh* feMesh;
   DofLayout* dofs;
   FeEquationNumber* eqNum;
   Variable_Register* varReg;
   Variable* vars[2];
   int maxDecomp[3] = {0, 1, 1};
   int sizes[3];
   double minCrd[3];
   double maxCrd[3];
   SizeT dataOffs = 1;
   Variable_DataType dataType = Variable_DataType_Double;
   int nDataTypes = 1;
   char* dataNames = "nothing";
   static SizeT structSize = sizeof(double);
   static int arraySize;
   static void* arrayPtrs[2];
   int nRanks;
   VariableCondition* bcs;
   ConditionFunction_Register* cfReg;
   Dictionary* dict;
   XML_IO_Handler* ioHandler;
   char filename[PCU_PATH_MAX];

   MPI_Comm_size( MPI_COMM_WORLD, &nRanks );
   sizes[0] = nRanks * 2;
   sizes[1] = sizes[2] = 2;
   minCrd[0] = minCrd[1] = minCrd[2] = 0.0;
   maxCrd[0] = maxCrd[1] = maxCrd[2] = (double)nRanks;

   gen = CartesianGenerator_New( "", NULL );
   CartesianGenerator_SetDimSize( gen, 3 );
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
   arrayPtrs[0] = Memory_Alloc_Array_Unnamed( double, arraySize );
   arrayPtrs[1] = Memory_Alloc_Array_Unnamed( double, arraySize );
   vars[0] = Variable_New( "one", 1, &dataOffs, &dataType, (unsigned*)&nDataTypes, 
		&dataNames, &structSize, (unsigned*)&arraySize, NULL, arrayPtrs, varReg );
	vars[1] = Variable_New( "two", 1, &dataOffs, &dataType, (unsigned*)&nDataTypes, 
		&dataNames, &structSize, (unsigned*)&arraySize, NULL, arrayPtrs + 1, varReg );

   dofs = DofLayout_New( "", NULL, varReg, 0, feMesh );
   dofs->nBaseVariables = 2;
   dofs->baseVariables = Memory_Alloc_Array_Unnamed( Variable*, 2 );
   dofs->baseVariables[0] = vars[0];
   dofs->baseVariables[1] = vars[1];
   Stg_Component_Build( dofs, NULL, False );
   Stg_Component_Initialise( dofs, NULL, False );

   ioHandler = XML_IO_Handler_New();
   dict = Dictionary_New();
   pcu_filename_input( "wallVC.xml", filename );
   IO_Handler_ReadAllFromFile( ioHandler, filename, dict );
   bcs = (VariableCondition*)WallVC_New( "", "wallVC", varReg, cfReg, dict, feMesh );
   Stg_Component_Build( bcs, NULL, False );
   Stg_Component_Initialise( bcs, NULL, False );

   eqNum = FeEquationNumber_New( "", feMesh, dofs, bcs, NULL );
   Stg_Component_Build( eqNum, NULL, False );
   Stg_Component_Initialise( eqNum, NULL, False );

   return eqNum;
}

void FeEquationNumberSuite_Setup( FeEquationNumberSuiteData* data ) {
}

void FeEquationNumberSuite_Teardown( FeEquationNumberSuiteData* data ) {
}


void FeEquationNumberSuite_TestLocal( FeEquationNumberSuiteData* data ) {
   FeEquationNumber* eqNum;
   FeMesh* feMesh;
   int eqNumsPerProc;
   int curEqNum;
   int nDofs;
   int rank;
   int n_i, dof_i;

   eqNum = buildEqNum();
   pcu_check_true( eqNum );

   MPI_Comm_rank( MPI_COMM_WORLD, &rank );
   feMesh = eqNum->feMesh;
   eqNumsPerProc = (Mesh_GetDimSize( feMesh ) == 3) ? 27 : (Mesh_GetDimSize( feMesh ) == 2) ? 9 : 3;
   curEqNum = eqNumsPerProc * rank;
   if( rank >= 1 )
      curEqNum -= (eqNumsPerProc / 3) * (rank - 1);
   curEqNum *= 2;

   for( n_i = 0; n_i < Mesh_GetLocalSize( feMesh, MT_VERTEX ); n_i++ ) {
      nDofs = eqNum->dofLayout->dofCounts[n_i];
      for( dof_i = 0; dof_i < nDofs; dof_i++ ) {
	 if( eqNum->destinationArray[n_i][dof_i] != curEqNum++ )
	    break;
      }
      if( dof_i < nDofs )
	 break;
   }
   pcu_check_true( n_i == Mesh_GetLocalSize( feMesh, 0 ) );

   FreeObject( eqNum );
}

void FeEquationNumberSuite_TestShadow( FeEquationNumberSuiteData* data ) {
   FeEquationNumber* eqNum;
   FeMesh* feMesh;
   int eqNumsPerProc;
   int curEqNum;
   int nDofs;
   int rank;
   int nLocalNodes, nDomainNodes;
   int n_i, dof_i;

   eqNum = buildEqNum();
   pcu_check_true( eqNum );

   MPI_Comm_rank( MPI_COMM_WORLD, &rank );
   feMesh = eqNum->feMesh;
   eqNumsPerProc = (Mesh_GetDimSize( feMesh ) == 3) ? 27 : (Mesh_GetDimSize( feMesh ) == 2) ? 9 : 3;
   if( rank >= 2 ) {
      curEqNum = eqNumsPerProc;
      curEqNum += 2 * (eqNumsPerProc / 3) * (rank - 2);
      curEqNum += 1;
   }
   else if( rank == 1 )
      curEqNum = 2;
   curEqNum *= 2;

   nLocalNodes = Mesh_GetLocalSize( feMesh, MT_VERTEX );
   nDomainNodes = Mesh_GetDomainSize( feMesh, MT_VERTEX );
   for( n_i = nLocalNodes; n_i < nDomainNodes; n_i++ ) {
      nDofs = eqNum->dofLayout->dofCounts[n_i];
      for( dof_i = 0; dof_i < nDofs; dof_i++ ) {
	 if( eqNum->destinationArray[n_i][dof_i] != curEqNum++ )
	    break;
      }
      if( dof_i < nDofs )
	 break;

      if( rank == 1 )
	 curEqNum += 4;
      else
	 curEqNum += 2;
   }
   pcu_check_true( n_i == Mesh_GetDomainSize( feMesh, 0 ) );

   FreeObject( eqNum );
}

void FeEquationNumberSuite_TestBCs( FeEquationNumberSuiteData* data ) {
   FeEquationNumber* eqNum;
   FeMesh* feMesh;
   unsigned nLocalNodes, nDomainNodes;
   unsigned eqNumsPerProc;
   unsigned curEqNum;
   unsigned nDofs;
   unsigned dof_i;
   unsigned eq;
   int rank;
   unsigned n_i;

   eqNum = buildEqNumBCs();
   pcu_check_true( eqNum );

   MPI_Comm_rank( MPI_COMM_WORLD, &rank );
   feMesh = eqNum->feMesh;
   eqNumsPerProc = (Mesh_GetDimSize( feMesh ) == 3) ? 18 : (Mesh_GetDimSize( feMesh ) == 2) ? 6 : 0;
   if( rank >= 2 ) {
      curEqNum = eqNumsPerProc;
      curEqNum += 2 * (eqNumsPerProc / 3) * (rank - 2);
      curEqNum += 1;
   }
   else if( rank == 1 )
      curEqNum = 2;
   curEqNum *= 2;

   nLocalNodes = Mesh_GetLocalSize( feMesh, MT_VERTEX );
   nDomainNodes = Mesh_GetDomainSize( feMesh, MT_VERTEX );
   for( n_i = nLocalNodes; n_i < nDomainNodes; n_i++ ) {
      nDofs = eqNum->dofLayout->dofCounts[n_i];
      for( dof_i = 0; dof_i < nDofs; dof_i++ ) {
	 eq = eqNum->destinationArray[n_i][dof_i];

	 if( n_i % 3 == 0 ) {
	    if( eq != (unsigned)-1 )
	       break;
	 }
	 else if( eq != curEqNum++ )
	    break;
      }
      if( dof_i < nDofs )
	 break;

      if( n_i % 3 != 0 ) {
	 if( rank == 1 )
	    curEqNum += 4;
	 else
	    curEqNum += 2;
      }
   }
   pcu_check_true( n_i == nDomainNodes );

   FreeObject( eqNum );
}

void FeEquationNumberSuite( pcu_suite_t* suite ) {
   pcu_suite_setData( suite, FeEquationNumberSuiteData );
   pcu_suite_setFixtures( suite, FeEquationNumberSuite_Setup, FeEquationNumberSuite_Teardown );
   pcu_suite_addTest( suite, FeEquationNumberSuite_TestLocal );
   pcu_suite_addTest( suite, FeEquationNumberSuite_TestShadow );
   pcu_suite_addTest( suite, FeEquationNumberSuite_TestBCs );
}
