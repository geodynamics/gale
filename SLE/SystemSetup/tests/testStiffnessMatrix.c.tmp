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
** $Id: testStiffnessMatrix.c 3664 2006-07-04 04:26:57Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <mpi.h>

#include "StGermain/StGermain.h"
#include "StgFEM/Discretisation/Discretisation.h"
#include "StgFEM/SLE/LinearAlgebra/LinearAlgebra.h"
#include "StgFEM/SLE/SystemSetup/SystemSetup.h"


struct _Particle {
	__IntegrationPoint
};


void identityAssembly( void* term, StiffnessMatrix* stiffMat, 
		       unsigned element, 
		       SystemLinearEquations* sle, 
		       double** elStiffMat )
{
	FeMesh*			feMesh = stiffMat->rowVariable->feMesh;
	FeEquationNumber*	eqNum = stiffMat->rowVariable->eqNum;
	DofLayout*		dofs = stiffMat->rowVariable->dofLayout;
	unsigned		nIncNodes, *incNodes;
	unsigned		nIncEls;
	unsigned		row, col;
	unsigned		eqI, eqJ;
	unsigned		n_i, n_j;
	unsigned		dof_i, dof_j;

	FeMesh_GetElementNodes( feMesh, element, &nIncNodes, &incNodes );

	row = 0;
	for( n_i = 0; n_i < nIncNodes; n_i++ ) {
		for( dof_i = 0; dof_i < dofs->dofCounts[n_i]; dof_i++ ) {
			eqI = eqNum->destinationArray[n_i][dof_i];
			col = 0;

			for( n_j = 0; n_j < nIncNodes; n_j++ ) {
				for( dof_j = 0; dof_j < dofs->dofCounts[n_j]; dof_j++ ) {
					eqJ = eqNum->destinationArray[n_j][dof_j];

					if( eqI == eqJ ) {
						nIncEls = FeMesh_GetNodeElementSize( feMesh, n_i );
						elStiffMat[row][col] = 1.0 / (double)nIncEls;
					}
					else
						elStiffMat[row][col] = 0.0;

					col++;
				}
			}

			row++;
		}
	}
}


StiffnessMatrix* buildStiffMat( unsigned nProcs ) {
	CartesianGenerator*		gen;
	FeMesh*				feMesh;
	DofLayout*			dofs;
	VariableCondition*		bcs;
	FeEquationNumber*		eqNum;
	FeVariable*			feVar;
	ForceVector*			rhs;
	StiffnessMatrix*		stiffMat;
	CellLayout*			cellLayout;
	ParticleLayout*			particleLayout;
	Swarm*				swarm;
	StiffnessMatrixTerm*		term;
	ExtensionManager_Register*	emReg;
	Variable_Register*		varReg;
	ConditionFunction_Register*	cfReg;
	FieldVariable_Register*		fvReg;
	EntryPoint_Register*		epReg;
	Variable*			vars[2];
	unsigned			sizes[3];
	double				minCrd[3];
	double				maxCrd[3];
	Bool				dimExists[3];
	unsigned			particlesPerDim[3];
	SizeT				dataOffs = 1;
	Variable_DataType		dataType = Variable_DataType_Double;
	unsigned			nDataTypes = 1;
	char*				dataNames = "nothing";
	static SizeT			structSize = sizeof(double);
	static unsigned			arraySize;
	static void*			arrayPtrs[2];
	Dictionary*			dict;
	XML_IO_Handler*			ioHandler;

	sizes[0] = sizes[1] = sizes[2] = nProcs * 2;
	minCrd[0] = minCrd[1] = minCrd[2] = 0.0;
	maxCrd[0] = maxCrd[1] = maxCrd[2] = (double)nProcs;
	dimExists[0] = dimExists[1] = dimExists[2] = True;
	particlesPerDim[0] = particlesPerDim[1] = particlesPerDim[2] = 2;

	gen = CartesianGenerator_New( "" );
	CartesianGenerator_SetDimSize( gen, 3 );
	CartesianGenerator_SetTopologyParams( gen, sizes, 0, NULL, NULL );
	CartesianGenerator_SetGeometryParams( gen, minCrd, maxCrd );
	CartesianGenerator_SetShadowDepth( gen, 0 );

	feMesh = FeMesh_New( "" );
	Mesh_SetGenerator( feMesh, gen );
	FeMesh_SetElementFamily( feMesh, "linear" );
	Build( feMesh, NULL, False );

	emReg = ExtensionManager_Register_New();
	varReg = Variable_Register_New();
	fvReg = FieldVariable_Register_New();
	epReg = EntryPoint_Register_New();

	arraySize = FeMesh_GetNodeDomainSize( feMesh );
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

	feVar = FeVariable_New( "", 
				feMesh, NULL, 
				dofs, bcs, NULL, 
				NULL, Mesh_GetDimSize( feMesh ), False, 
				StgFEM_Native_ImportExportType, StgFEM_Native_ImportExportType, fvReg );
	Build( feVar, NULL, False );
	Initialise( feVar, NULL, False );

	rhs = ForceVector_New( "", 
			       feVar, Mesh_GetDimSize( feMesh ), 
			       epReg, MPI_COMM_WORLD );
	Build( rhs, NULL, False );
	Initialise( rhs, NULL, False );

	stiffMat = StiffnessMatrix_New( "", 
					feVar, feVar, 
					rhs, 
					NULL, Mesh_GetDimSize( feMesh ), False, False, epReg, MPI_COMM_WORLD );

	cellLayout = (CellLayout*)SingleCellLayout_New( "", dimExists, minCrd, maxCrd );
	particleLayout = (ParticleLayout*)GaussParticleLayout_New( "", Mesh_GetDimSize( feMesh ), particlesPerDim );
	swarm = Swarm_New( "", cellLayout, particleLayout, Mesh_GetDimSize( feMesh ), sizeof(Particle), 
			   emReg, varReg, MPI_COMM_WORLD );

	term = StiffnessMatrixTerm_New( "", stiffMat, swarm, NULL );
	term->_assembleElement = identityAssembly;
	StiffnessMatrix_AddStiffnessMatrixTerm( stiffMat, term );

	Build( stiffMat, NULL, False );
	Initialise( stiffMat, NULL, False );

	return stiffMat;
}


Bool testStiffMat( unsigned rank, unsigned nProcs, unsigned watch ) {
	Bool			result = True;
	StiffnessMatrix*	stiffMat;

	stiffMat = buildStiffMat( nProcs );
	StiffnessMatrix_Assemble( stiffMat, True, NULL, NULL );

	if( rank == watch ) {
	}

done:
	FreeObject( stiffMat );

	return result;
}


#define nTests	1

TestSuite_Test	tests[nTests] = {{"test stiffness matrix", testStiffMat, 1}};


int main( int argc, char* argv[] ) {
	TestSuite*	suite;

	/* Initialise MPI, get world info. */
	MPI_Init( &argc, &argv );

	/* Initialise StGermain. */
	StGermain_Init( &argc, &argv );
	StgFEM_Discretisation_Init( &argc, &argv );
	StgFEM_SLE_LinearAlgebra_Init( &argc, &argv );
	StgFEM_SLE_SystemSetup_Init( &argc, &argv );

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
