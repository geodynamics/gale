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
**   Tests the FeVariableSuite
**
** $Id: testFeVariable.c 3462 2006-02-19 06:53:24Z WalterLandry $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>

#include "pcu/pcu.h"
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/Discretisation/Discretisation.h>
#include "FeVariableSuite.h"

#define EPSILON 1.0E-6

typedef struct {
} FeVariableSuiteData;

struct _Particle {
	__IntegrationPoint;
};

FeVariable* BuildFeVariable_AsPosition( unsigned dim ) {
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
	static double* arrayPtrs[3];
	int nRanks;
	Variable* var;
	FieldVariable_Register* fieldReg;
	FeVariable* feVar;
	int n_i;

	MPI_Comm_size( MPI_COMM_WORLD, &nRanks );
	sizes[0] = nRanks * 3;
	sizes[1] = sizes[2] = 3;
	minCrd[0] = minCrd[1] = minCrd[2] = 0.0;
	maxCrd[0] = maxCrd[1] = maxCrd[2] = (double)nRanks;

	gen = CartesianGenerator_New( "", NULL );
	CartesianGenerator_SetDimSize( gen, dim );
	CartesianGenerator_SetTopologyParams( gen, (unsigned*)sizes, 0, NULL, (unsigned*)maxDecomp );
	CartesianGenerator_SetGeometryParams( gen, minCrd, maxCrd );
	CartesianGenerator_SetShadowDepth( gen, 0 );

	feMesh = FeMesh_New( "", NULL );
	Mesh_SetGenerator( feMesh, gen );
	FeMesh_SetElementFamily( feMesh, "linear" );
	Stg_Component_Build( feMesh, NULL, False );

	varReg = Variable_Register_New();

	arraySize = Mesh_GetDomainSize( feMesh, MT_VERTEX );
	arrayPtrs[0] = Memory_Alloc_Array_Unnamed( double, arraySize * dim );

	var = Variable_NewVector( "velocity", Variable_DataType_Double, dim, (unsigned*)&arraySize, NULL,
	                          (void**)arrayPtrs, varReg, "vx", "vy", "vz" );
	Variable_Register_BuildAll( varReg );

	dofs = DofLayout_New( "", varReg, 0, feMesh );
	dofs->nBaseVariables = dim;
	dofs->baseVariables = Memory_Alloc_Array_Unnamed( Variable*, dim );
	dofs->baseVariables[0] = var->components[0];
	dofs->baseVariables[1] = var->components[1];
	dofs->baseVariables[2] = var->components[2];
	Stg_Component_Build( dofs, NULL, False );
	Stg_Component_Initialise( dofs, NULL, False );
	
	eqNum = FeEquationNumber_New( "", feMesh, dofs, NULL, NULL );
	Stg_Component_Build( eqNum, NULL, False );
	Stg_Component_Initialise( eqNum, NULL, False );

	fieldReg = FieldVariable_Register_New();
	feVar = FeVariable_New( "velocity", feMesh, NULL, dofs, NULL, NULL, NULL, dim, True, False, False, fieldReg );

	for( n_i = 0; n_i < Mesh_GetLocalSize( feMesh, 0 ); n_i++ ) {
		double* pos = Mesh_GetVertex( feMesh, n_i );
		Variable_SetValue( var, n_i, pos );
	}

	/* Build and initialise system */
	Stg_Component_Build( feVar, 0, False );
	Stg_Component_Initialise( feVar, 0, False );

	return feVar;
}

FeVariable* BuildFeVariable_AsConstant( unsigned dim ) {
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
	static double* arrayPtr;
	int nRanks;
	Variable* var;
	FieldVariable_Register* fieldReg;
	FeVariable* feVar;
	int n_i;
	double constant[3] = {1.0, 1.0, 1.0};

	MPI_Comm_size( MPI_COMM_WORLD, &nRanks );
	sizes[0] = nRanks * 3;
	sizes[1] = sizes[2] = 3;
	minCrd[0] = minCrd[1] = minCrd[2] = 0.0;
	maxCrd[0] = maxCrd[1] = maxCrd[2] = (double)nRanks;

	gen = CartesianGenerator_New( "", NULL );
	CartesianGenerator_SetDimSize( gen, dim );
	CartesianGenerator_SetTopologyParams( gen, (unsigned*)sizes, 0, NULL, (unsigned*)maxDecomp );
	CartesianGenerator_SetGeometryParams( gen, minCrd, maxCrd );
	CartesianGenerator_SetShadowDepth( gen, 0 );

	feMesh = FeMesh_New( "", NULL );
	Mesh_SetGenerator( feMesh, gen );
	FeMesh_SetElementFamily( feMesh, "linear" );
	Stg_Component_Build( feMesh, NULL, False );

	varReg = Variable_Register_New();

	arraySize = Mesh_GetDomainSize( feMesh, MT_VERTEX );
	arrayPtr = Memory_Alloc_Array_Unnamed( double, arraySize );

	var = Variable_NewScalar( "pressure", Variable_DataType_Double, (unsigned*)&arraySize, NULL,
	                          (void**)&arrayPtr, varReg );
	Variable_Register_BuildAll( varReg );

	dofs = DofLayout_New( "", varReg, 0, feMesh );
	dofs->nBaseVariables = 1;
	dofs->baseVariables = Memory_Alloc_Array_Unnamed( Variable*, 1 );
	dofs->baseVariables[0] = var;
	Stg_Component_Build( dofs, NULL, False );
	Stg_Component_Initialise( dofs, NULL, False );
	
	eqNum = FeEquationNumber_New( "", feMesh, dofs, NULL, NULL );
	Stg_Component_Build( eqNum, NULL, False );
	Stg_Component_Initialise( eqNum, NULL, False );

	fieldReg = FieldVariable_Register_New();
	feVar = FeVariable_New( "pressure", feMesh, NULL, dofs, NULL, NULL, NULL, dim, True, False, False, fieldReg );

	for( n_i = 0; n_i < Mesh_GetLocalSize( feMesh, 0 ); n_i++ ) {
		Variable_SetValue( var, n_i, constant );
	}

	/* Build and initialise system */
	Stg_Component_Build( feVar, 0, False );
	Stg_Component_Initialise( feVar, 0, False );

	return feVar;
}

Swarm* BuildSwarm( FeMesh* mesh ) {
	ElementCellLayout*		elCellLayout;
	GaussParticleLayout* 		gaussLayout;
	ExtensionManager_Register*	extMgr_Reg;
	Swarm*				swarm;
	unsigned			dim		= Mesh_GetDimSize( mesh );
	unsigned			partPerDim[3]	= { 2, 2, 2 };
	
	extMgr_Reg = ExtensionManager_Register_New();
	elCellLayout = ElementCellLayout_New( "elementCellLayout", mesh );
	gaussLayout = GaussParticleLayout_New( "gaussParticleLayout", NULL, LocalCoordSystem, True, dim, partPerDim );
	swarm = Swarm_New( "gaussSwarm", NULL, elCellLayout, gaussLayout, dim, sizeof(Particle), extMgr_Reg, NULL, MPI_COMM_WORLD, NULL );

	Stg_Component_Build( swarm, NULL, True );
	Stg_Component_Initialise( swarm, NULL, True );

	return swarm;
}

void FeVariableSuite_Setup( FeVariableSuiteData* data ) {
}

void FeVariableSuite_Teardown( FeVariableSuiteData* data ) {
}

void FeVariableSuite_Interpolate( FeVariableSuiteData* data ) {
	FeVariable*		feVar		= BuildFeVariable_AsPosition( 3 );
	IArray*			inc		= IArray_New();
	unsigned		el_i;
	unsigned		vert_i;
	unsigned		nVerts;
	unsigned*		verts;
	double*			vert;
	double			value[3];
	InterpolationResult	interpRes;
	unsigned		dim_i;

	for( el_i = 0; el_i < Mesh_GetDomainSize( feVar->feMesh, 3 ); el_i++ ) {
		Mesh_GetIncidence( feVar->feMesh, 3, el_i, MT_VERTEX, inc );
		nVerts = IArray_GetSize( inc );
		verts = IArray_GetPtr( inc );

		for( vert_i = 0; vert_i < nVerts; vert_i++ ) {
			vert = Mesh_GetVertex( feVar->feMesh, verts[vert_i] );
			interpRes = FieldVariable_InterpolateValueAt( feVar, vert, value );

			if( interpRes != LOCAL && interpRes != SHADOW )
				continue;

			for( dim_i = 0; dim_i < 3; dim_i++ ) {
				pcu_check_true( fabs( vert[dim_i] - value[dim_i] ) < EPSILON );
			}
		}
	}

	NewClass_Delete( inc );
	Stg_Component_Destroy( feVar, NULL, True );
}

void FeVariableSuite_Integrate( FeVariableSuiteData* data ) {
	FeVariable*			feVar;
	Swarm*				swarm;
	double				globalMin[3];
	double				globalMax[3];
	double				volFromSize	= 1.0;
	double				volFromInt;
	unsigned			dim_i;

	feVar = BuildFeVariable_AsConstant( 3 );
	swarm = BuildSwarm( feVar->feMesh );

	Mesh_GetGlobalCoordRange( feVar->feMesh, globalMin, globalMax );
	for( dim_i = 0; dim_i < 3; dim_i++ )
		volFromSize *= globalMax[dim_i] - globalMin[dim_i];

	volFromInt = FeVariable_Integrate( feVar, swarm );

	pcu_check_true( fabs( volFromSize - volFromInt ) < EPSILON );
	
	Stg_Component_Destroy( feVar, NULL, True );
	Stg_Component_Destroy( swarm, NULL, True );
}

void FeVariableSuite( pcu_suite_t* suite ) {
   pcu_suite_setData( suite, FeVariableSuiteData );
   pcu_suite_setFixtures( suite, FeVariableSuite_Setup, FeVariableSuite_Teardown );
   pcu_suite_addTest( suite, FeVariableSuite_Interpolate );
   pcu_suite_addTest( suite, FeVariableSuite_Integrate );
}
