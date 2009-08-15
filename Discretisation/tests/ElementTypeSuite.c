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
**   Tests the ElementTypeSuite
**
** $Id: testElementType.c 3462 2006-02-19 06:53:24Z WalterLandry $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>

#include "pcu/pcu.h"
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/Discretisation/Discretisation.h>
#include "ElementTypeSuite.h"

typedef struct {
	Dictionary*	dict;
} ElementTypeSuiteData;

FeMesh* BuildMesh( unsigned nDims ) {
	CartesianGenerator*	gen;
	FeMesh*			feMesh;
	unsigned		maxDecomp[3] = {0, 1, 1};
	unsigned		sizes[3];
	double			minCrd[3];
	double			maxCrd[3];

	sizes[0] = sizes[1] = 6;
	sizes[2] = 1;
	minCrd[0] = minCrd[1] = minCrd[2] = 0.0;
	maxCrd[0] = maxCrd[1] = maxCrd[2] = 1.2;

	gen = CartesianGenerator_New( "" );
	CartesianGenerator_SetDimSize( gen, nDims );
	CartesianGenerator_SetTopologyParams( gen, sizes, 0, NULL, maxDecomp );
	CartesianGenerator_SetGeometryParams( gen, minCrd, maxCrd );
	CartesianGenerator_SetShadowDepth( gen, 0 );

	feMesh = FeMesh_New( "" );
	Mesh_SetGenerator( feMesh, gen );
	FeMesh_SetElementFamily( feMesh, "linear" );
	Stg_Component_Build( feMesh, NULL, False );
	Stg_Component_Initialise( feMesh, NULL, False );

	return feMesh;
}

void ElementTypeSuite_Setup( ElementTypeSuiteData* data ) {
	data->dict = Dictionary_New();
}

void ElementTypeSuite_Teardown( ElementTypeSuiteData* data ) {
	Stg_Class_Delete( data->dict );
}

#define TOLERANCE 1.0e-6
void ElementTypeSuite_Test2D( ElementTypeSuiteData* data ) {
	FeMesh*		feMesh;
	unsigned	maxTests	= 40;
	unsigned	test_i;
	Coord		gCoord, lCoord, gCoord_fromLocal;
	unsigned	el, elNodeCount, elNode_i;
	ElementType*	elType;
	IArray*		inc		= IArray_New();
	unsigned*	elNodes;
	double		Ni[9];
	double		vecNorm;

	feMesh = BuildMesh( 2 );

	srand48( 0 );
	for( test_i = 0; test_i < maxTests; test_i++ ) {
		gCoord[I_AXIS] = drand48();
		gCoord[J_AXIS] = drand48();
	
		Mesh_Algorithms_SearchElements( feMesh->algorithms, gCoord, &el );
		elNodeCount = FeMesh_GetElementNodeSize( feMesh, el );
		elType = FeMesh_GetElementType( feMesh, el );
		_ElementType_ConvertGlobalCoordToElLocal( elType, feMesh, el, gCoord, lCoord );
		ElementType_EvaluateShapeFunctionsAt( elType, lCoord, Ni );

		Mesh_GetIncidence( feMesh, 2, el, MT_VERTEX, inc );
		elNodes = IArray_GetPtr( inc );
		memset( gCoord_fromLocal, 0, sizeof( double ) * 2 );
		for( elNode_i = 0; elNode_i < elNodeCount; elNode_i++ ) {
			gCoord_fromLocal[I_AXIS] += Ni[elNode_i] * feMesh->verts[elNodes[elNode_i]][I_AXIS];
			gCoord_fromLocal[J_AXIS] += Ni[elNode_i] * feMesh->verts[elNodes[elNode_i]][J_AXIS];
		}

		vecNorm = sqrt( (gCoord[I_AXIS] - gCoord_fromLocal[I_AXIS])*(gCoord[I_AXIS] - gCoord_fromLocal[I_AXIS]) + 
				(gCoord[J_AXIS] - gCoord_fromLocal[J_AXIS])*(gCoord[J_AXIS] - gCoord_fromLocal[J_AXIS]) );

		pcu_check_true( vecNorm < TOLERANCE );
	}

	NewClass_Delete( inc );
	Stg_Component_Destroy( feMesh, NULL, True );	
}


void ElementTypeSuite( pcu_suite_t* suite ) {
   pcu_suite_setData( suite, ElementTypeSuiteData );
   pcu_suite_setFixtures( suite, ElementTypeSuite_Setup, ElementTypeSuite_Teardown );
   pcu_suite_addTest( suite, ElementTypeSuite_Test2D );
}
