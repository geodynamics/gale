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
** $Id: testSROpGenerator.c 678 2006-12-22 01:45:53Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <mpi.h>

#include "StGermain/StGermain.h"
#include "Discretisation/Discretisation.h"
#include "SLE/LinearAlgebra/LinearAlgebra.h"


FeEquationNumber* buildEqNumWithBCs( unsigned nProcs ) {
	CartesianGenerator*		gen;
	FeMesh*				feMesh;
	DofLayout*			dofs;
	FeEquationNumber*		eqNum;
	VariableCondition*		bcs;
	Variable_Register*		varReg;
	ConditionFunction_Register*	cfReg;
	Variable*			vars[2];
	unsigned			maxDecomp[3] = {0, 1, 1};
	unsigned			sizes[3];
	double				minCrd[3];
	double				maxCrd[3];
	SizeT				dataOffs = 1;
	Variable_DataType		dataType = Variable_DataType_Double;
	unsigned			nDataTypes = 1;
	char*				dataNames = "nothing";
	static SizeT			structSize = sizeof(double);
	static unsigned			arraySize;
	static void*			arrayPtrs[2];
	Dictionary*			dict;
	XML_IO_Handler*			ioHandler;

	sizes[1] = sizes[2] = sizes[0] = nProcs * 8;
	minCrd[0] = minCrd[1] = minCrd[2] = 0.0;
	maxCrd[0] = maxCrd[1] = maxCrd[2] = (double)nProcs;

	gen = CartesianGenerator_New( "" );
	CartesianGenerator_SetDimSize( gen, 3 );
	CartesianGenerator_SetTopologyParams( gen, sizes, 0, NULL, maxDecomp );
	CartesianGenerator_SetGeometryParams( gen, minCrd, maxCrd );
	CartesianGenerator_SetShadowDepth( gen, 0 );

	feMesh = FeMesh_New( "" );
	Mesh_SetGenerator( feMesh, gen );
	FeMesh_SetElementFamily( feMesh, "linear" );
	Build( feMesh, NULL, False );

	varReg = Variable_Register_New();
	cfReg = ConditionFunction_Register_New();

	arraySize = Mesh_GetDomainSize( feMesh, MT_VERTEX );
	arrayPtrs[0] = Memory_Alloc_Array_Unnamed( double, arraySize );
	arrayPtrs[1] = Memory_Alloc_Array_Unnamed( double, arraySize );
	vars[0] = Variable_New( "one", 1, &dataOffs, &dataType, &nDataTypes, &dataNames, 
				&structSize, &arraySize, arrayPtrs, varReg );
	vars[1] = Variable_New( "two", 1, &dataOffs, &dataType, &nDataTypes, &dataNames, 
				&structSize, &arraySize, arrayPtrs + 1, varReg );

	dofs = DofLayout_New( "", varReg, 0, feMesh );
	dofs->nBaseVariables = 2;
	dofs->baseVariables = Memory_Alloc_Array_Unnamed( Variable*, 2 );
	dofs->baseVariables[0] = vars[0];
	dofs->baseVariables[1] = vars[1];
	Build( dofs, NULL, False );
	Initialise( dofs, NULL, False );

	ioHandler = XML_IO_Handler_New();
	dict = Dictionary_New();
	IO_Handler_ReadAllFromFile( ioHandler, "data/wallVC.xml", dict );
	bcs = (VariableCondition*)WallVC_New( "", "wallVC", varReg, cfReg, dict, feMesh );
	Build( bcs, NULL, False );
	Initialise( bcs, NULL, False );

	eqNum = FeEquationNumber_New( "", feMesh, dofs, bcs, NULL );
	Build( eqNum, NULL, False );
	Initialise( eqNum, NULL, False );

	return eqNum;
}


Bool testOps( unsigned rank, unsigned nProcs, unsigned watch ) {
	Bool			result = True;
	FeEquationNumber*	eqNum;
	MultigridSolver*	solver;
	SROpGenerator*		opGen;
	Matrix			**pOps, **rOps;
	Matrix*			mat;

	mat = (Matrix*)PETScMatrix_New( "" );
	Matrix_SetLocalSize( mat, nProcs, nProcs );
	Matrix_AssemblyBegin( mat );
	Matrix_AssemblyEnd( mat );

	eqNum = buildEqNumWithBCs( nProcs );
	solver = MultigridSolver_New( "" );
	MultigridSolver_SetLevels( solver, 4 );

	MatrixSolver_SetMatrix( solver, mat );
	opGen = SROpGenerator_New( "" );
	MGOpGenerator_SetMatrixSolver( opGen, solver );
	/*SROpGenerator_SetVariable( opGen, eqNum ); TODO */
	Build( opGen, NULL, False );
	SROpGenerator_Generate( opGen, &pOps, &rOps );

	if( rank == watch ) {
	}

done:
	FreeObject( opGen );
	FreeObject( mat );
	FreeObject( solver );
	FreeObject( eqNum );

	return result;
}


#define nTests	1

TestSuite_Test	tests[nTests] = {{"test operator generation", testOps, 1}};


int main( int argc, char* argv[] ) {
	TestSuite*	suite;

	/* Initialise MPI, get world info. */
	MPI_Init( &argc, &argv );

	/* Initialise StGermain. */
	StGermain_Init( &argc, &argv );
	StgFEM_Discretisation_Init( &argc, &argv );
	StgFEM_SLE_LinearAlgebra_Init( &argc, &argv );

	/* Disable all journaling. */
	Journal_Enable_TypedStream( Debug_Type, False );
	Journal_Enable_TypedStream( Info_Type, False );

	/* Create the test suite. */
	suite = TestSuite_New();
	TestSuite_SetProcToWatch( suite, (argc >= 2) ? atoi( argv[1] ) : 0 );
	TestSuite_SetTests( suite, nTests, tests );

	/* Run the tests. */
	TestSuite_Run( suite );

	/* Destroy test suites. */
	FreeObject( suite );

	/* Finalise StGermain. */
	BaseContainer_Finalise();
	BaseIO_Finalise();
	BaseFoundation_Finalise();

	/* Close off MPI */
	MPI_Finalize();

	return MPI_SUCCESS;
}
