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
** $Id: testFeEquationNumber.c 822 2007-04-27 06:20:35Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <mpi.h>

#include "StGermain/StGermain.h"
#include "Discretisation/Discretisation.h"


FeEquationNumber* buildEqNum( unsigned nProcs ) {
	CartesianGenerator*	gen;
	FeMesh*			feMesh;
	DofLayout*		dofs;
	FeEquationNumber*	eqNum;
	Variable_Register*	varReg;
	Variable*		vars[2];
	unsigned		maxDecomp[3] = {0, 1, 1};
	unsigned		sizes[3];
	double			minCrd[3];
	double			maxCrd[3];
	SizeT			dataOffs = 1;
	Variable_DataType	dataType = Variable_DataType_Double;
	unsigned		nDataTypes = 1;
	char*			dataNames = "nothing";
	static SizeT		structSize = sizeof(double);
	static unsigned		arraySize;
	static void*		arrayPtrs[2];

	sizes[0] = nProcs * 2;
	sizes[1] = sizes[2] = 2;
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

	eqNum = FeEquationNumber_New( "", feMesh, dofs, NULL, NULL );
	Build( eqNum, NULL, False );
	Initialise( eqNum, NULL, False );

	return eqNum;
}

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

	sizes[0] = nProcs * 2;
	sizes[1] = sizes[2] = 2;
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


Bool testLocalDstArray( unsigned rank, unsigned nProcs, unsigned watch ) {
	Bool			result = True;
	FeEquationNumber*	eqNum;

	eqNum = buildEqNum( nProcs );

	if( rank == watch ) {
		FeMesh*		feMesh;
		unsigned	eqNumsPerProc;
		unsigned	curEqNum;
		unsigned	n_i;

		feMesh = eqNum->feMesh;
		eqNumsPerProc = (Mesh_GetDimSize( feMesh ) == 3) ? 27 :
			(Mesh_GetDimSize( feMesh ) == 2) ? 9 :
			3;
		curEqNum = eqNumsPerProc * rank;
		if( rank >= 1 )
			curEqNum -= (eqNumsPerProc / 3) * (rank - 1);
		curEqNum *= 2;

		for( n_i = 0; n_i < Mesh_GetLocalSize( feMesh, MT_VERTEX ); n_i++ ) {
			unsigned	nDofs;
			unsigned	dof_i;

			nDofs = eqNum->dofLayout->dofCounts[n_i];
			for( dof_i = 0; dof_i < nDofs; dof_i++ ) {
				if( eqNum->destinationArray[n_i][dof_i] != curEqNum++ ) {
					result = False;
					goto done;
				}
			}
		}
	}

done:
	FreeObject( eqNum );

	return result;
}

Bool testShadowDstArray( unsigned rank, unsigned nProcs, unsigned watch ) {
	Bool			result = True;
	FeEquationNumber*	eqNum;

	eqNum = buildEqNum( nProcs );

	if( rank == watch ) {
		FeMesh*		feMesh;
		unsigned	nLocalNodes, nDomainNodes;
		unsigned	eqNumsPerProc;
		unsigned	curEqNum;
		unsigned	n_i;

		feMesh = eqNum->feMesh;
		eqNumsPerProc = (Mesh_GetDimSize( feMesh ) == 3) ? 27 :
			(Mesh_GetDimSize( feMesh ) == 2) ? 9 :
			3;
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
			unsigned	nDofs;
			unsigned	dof_i;

			nDofs = eqNum->dofLayout->dofCounts[n_i];
			for( dof_i = 0; dof_i < nDofs; dof_i++ ) {
				if( eqNum->destinationArray[n_i][dof_i] != curEqNum++ ) {
					result = False;
					goto done;
				}
			}

			if( rank == 1 )
				curEqNum += 4;
			else
				curEqNum += 2;
		}
	}

done:
	FreeObject( eqNum );

	return result;
}

Bool testBCs( unsigned rank, unsigned nProcs, unsigned watch ) {
	Bool			result = True;
	FeEquationNumber*	eqNum;

	eqNum = buildEqNumWithBCs( nProcs );

	if( rank == watch ) {
		FeMesh*		feMesh;
		unsigned	nLocalNodes, nDomainNodes;
		unsigned	eqNumsPerProc;
		unsigned	curEqNum;
		unsigned	n_i;

		feMesh = eqNum->feMesh;
		eqNumsPerProc = (Mesh_GetDimSize( feMesh ) == 3) ? 18 :
			(Mesh_GetDimSize( feMesh ) == 2) ? 6 :
			0;
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
			unsigned	nDofs;
			unsigned	dof_i;

			nDofs = eqNum->dofLayout->dofCounts[n_i];
			for( dof_i = 0; dof_i < nDofs; dof_i++ ) {
				unsigned	eq;

				eq = eqNum->destinationArray[n_i][dof_i];

				if( n_i % 3 == 0 ) {
					if( eq != (unsigned)-1 ) {
						result = False;
						goto done;
					}
				}
				else if( eq != curEqNum++ ) {
					result = False;
					goto done;
				}
			}

			if( n_i % 3 != 0 ) {
				if( rank == 1 )
					curEqNum += 4;
				else
					curEqNum += 2;
			}
		}
	}

done:
	FreeObject( eqNum );

	return result;
}


#define nTests	3

TestSuite_Test	tests[nTests] = {{"test local destination array", testLocalDstArray, 1}, 
				 {"test shadow destination array", testShadowDstArray, 1}, 
				 {"test destination array with BCs", testBCs, 1}};


int main( int argc, char* argv[] ) {
	TestSuite*	suite;

	/* Initialise MPI, get world info. */
	MPI_Init( &argc, &argv );

	/* Initialise StGermain. */
	StGermain_Init( &argc, &argv );
	StgFEM_Discretisation_Init( &argc, &argv );

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
