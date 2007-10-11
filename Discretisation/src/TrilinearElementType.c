/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003-2006, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street,
**	Melbourne, 3053, Australia.
**
** Primary Contributing Organisations:
**	Victorian Partnership for Advanced Computing Ltd, Computational Software Development - http://csd.vpac.org
**	Australian Computational Earth Systems Simulator - http://www.access.edu.au
**	Monash Cluster Computing - http://www.mcc.monash.edu.au
**	Computational Infrastructure for Geodynamics - http://www.geodynamics.org
**
** Contributors:
**	Patrick D. Sunter, Software Engineer, VPAC. (pds@vpac.org)
**	Robert Turnbull, Research Assistant, Monash University. (robert.turnbull@sci.monash.edu.au)
**	Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**	David May, PhD Student, Monash University (david.may@sci.monash.edu.au)
**	Louis Moresi, Associate Professor, Monash University. (louis.moresi@sci.monash.edu.au)
**	Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**	Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
**	Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
**	Julian Giordani, Research Assistant, Monash University. (julian.giordani@sci.monash.edu.au)
**	Vincent Lemiale, Postdoctoral Fellow, Monash University. (vincent.lemiale@sci.monash.edu.au)
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
** $Id: TrilinearElementType.c 964 2007-10-11 08:03:06Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include "units.h"
#include "types.h"
#include "shortcuts.h"
#include "ElementType.h"
#include "TrilinearElementType.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

const Type TrilinearElementType_Type = "TrilinearElementType";
#define _TrilinearElementType_NodeCount 8

void* TrilinearElementType_DefaultNew( Name name ) {
	return _TrilinearElementType_New( sizeof(TrilinearElementType), TrilinearElementType_Type, _TrilinearElementType_Delete,
		_TrilinearElementType_Print, NULL, TrilinearElementType_DefaultNew, _TrilinearElementType_Construct,
		_TrilinearElementType_Build, _TrilinearElementType_Initialise, _TrilinearElementType_Execute, _TrilinearElementType_Destroy,
		name, False, _TrilinearElementType_SF_allNodes, 
		_TrilinearElementType_SF_allLocalDerivs_allNodes, _TrilinearElementType_ConvertGlobalCoordToElLocal,
		_TrilinearElementType_NodeCount );
}

TrilinearElementType* TrilinearElementType_New( Name name ) {
	return _TrilinearElementType_New( sizeof(TrilinearElementType), TrilinearElementType_Type, _TrilinearElementType_Delete,
		_TrilinearElementType_Print, NULL, TrilinearElementType_DefaultNew, _TrilinearElementType_Construct,
		_TrilinearElementType_Build, _TrilinearElementType_Initialise, _TrilinearElementType_Execute, _TrilinearElementType_Destroy,
		name, True, _TrilinearElementType_SF_allNodes, 
		_TrilinearElementType_SF_allLocalDerivs_allNodes, _TrilinearElementType_ConvertGlobalCoordToElLocal,
		_TrilinearElementType_NodeCount );
}


TrilinearElementType* _TrilinearElementType_New( 
		SizeT								_sizeOfSelf,
		Type								type,
		Stg_Class_DeleteFunction*						_delete,
		Stg_Class_PrintFunction*						_print,
		Stg_Class_CopyFunction*						_copy, 
		Stg_Component_DefaultConstructorFunction*	_defaultConstructor,
		Stg_Component_ConstructFunction*			_construct,
		Stg_Component_BuildFunction*		_build,
		Stg_Component_InitialiseFunction*		_initialise,
		Stg_Component_ExecuteFunction*		_execute,
		Stg_Component_DestroyFunction*		_destroy,
		Name							name,
		Bool							initFlag,
		ElementType_EvaluateShapeFunctionsAtFunction*			_evaluateShapeFunctionsAt,
		ElementType_EvaluateShapeFunctionLocalDerivsAtFunction*		_evaluateShapeFunctionLocalDerivsAt,
		ElementType_ConvertGlobalCoordToElLocalFunction*		_convertGlobalCoordToElLocal,
		Index								nodeCount )
{
	TrilinearElementType*		self;
	
	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(TrilinearElementType) );
	self = (TrilinearElementType*)_ElementType_New( _sizeOfSelf, type, _delete, _print, _copy, _defaultConstructor,
			_construct, _build, _initialise, _execute, _destroy, name, initFlag, _evaluateShapeFunctionsAt,
		_evaluateShapeFunctionLocalDerivsAt, _convertGlobalCoordToElLocal, nodeCount );
	
	/* General info */
	
	/* Virtual functions */
	
	/* TrilinearElementType info */
	if( initFlag ){
		_TrilinearElementType_Init( self );
	}
	
	return self;
}



void _TrilinearElementType_Init( TrilinearElementType* self ) {
	Dimension_Index dim_I=0;

	/* General and Virtual info should already be set */
	
	/* TrilinearElementType info */
	self->isConstructed = True;
	for ( dim_I = 0; dim_I < 3; dim_I++ ) {
		self->minElLocalCoord[dim_I] = -1;
		self->maxElLocalCoord[dim_I] = 1;
		self->elLocalLength[dim_I] = self->maxElLocalCoord[dim_I] - self->minElLocalCoord[dim_I];
	}

	/* Set up the tetrahedral indices. */
	self->tetInds = Memory_Alloc_2DArray( unsigned, 10, 4, "Mesh_HexType::tetInds" );
	self->tetInds[0][0] = 0; self->tetInds[0][1] = 1; self->tetInds[0][2] = 2; self->tetInds[0][3] = 4;
	self->tetInds[1][0] = 1; self->tetInds[1][1] = 2; self->tetInds[1][2] = 3; self->tetInds[1][3] = 7;
	self->tetInds[2][0] = 1; self->tetInds[2][1] = 4; self->tetInds[2][2] = 5; self->tetInds[2][3] = 7;
	self->tetInds[3][0] = 2; self->tetInds[3][1] = 4; self->tetInds[3][2] = 6; self->tetInds[3][3] = 7;
	self->tetInds[4][0] = 1; self->tetInds[4][1] = 2; self->tetInds[4][2] = 4; self->tetInds[4][3] = 7;
	self->tetInds[5][0] = 0; self->tetInds[5][1] = 1; self->tetInds[5][2] = 3; self->tetInds[5][3] = 5;
	self->tetInds[6][0] = 0; self->tetInds[6][1] = 4; self->tetInds[6][2] = 5; self->tetInds[6][3] = 6;
	self->tetInds[7][0] = 0; self->tetInds[7][1] = 2; self->tetInds[7][2] = 3; self->tetInds[7][3] = 6;
	self->tetInds[8][0] = 3; self->tetInds[8][1] = 5; self->tetInds[8][2] = 6; self->tetInds[8][3] = 7;
	self->tetInds[9][0] = 0; self->tetInds[9][1] = 3; self->tetInds[9][2] = 5; self->tetInds[9][3] = 6;
}

void _TrilinearElementType_Delete( void* elementType ) {
	TrilinearElementType* self = (TrilinearElementType*)elementType;
	Journal_DPrintf( self->debug, "In %s\n", __func__ );

	FreeArray( self->tetInds );
	
	/* Stg_Class_Delete parent*/
	_ElementType_Delete( self );
}

void _TrilinearElementType_Print( void* elementType, Stream* stream ) {
	TrilinearElementType* self = (TrilinearElementType*)elementType;
	
	/* Set the Journal for printing informations */
	Stream* trilinearElementTypeStream = stream;
	
	/* General info */
	Journal_Printf( trilinearElementTypeStream, "TrilinearElementType (ptr): %p\n", self );
	
	/* Print parent */
	_ElementType_Print( self, trilinearElementTypeStream );
	
	/* Virtual info */
	
	/* TrilinearElementType info */
}

void _TrilinearElementType_Construct( void* elementType, Stg_ComponentFactory *cf, void* data ){
	
}
	
void _TrilinearElementType_Initialise( void* elementType, void *data ){
	
}
	
void _TrilinearElementType_Execute( void* elementType, void *data ){
	
}
	
void _TrilinearElementType_Destroy( void* elementType, void *data ){
	
}

void _TrilinearElementType_Build( void* elementType, void *data ) {
	
}

#if 0
void _TrilinearElementType_ConvertGlobalCoordToElementLocal( void* elementType, Element_DomainIndex element,const Coord globalCoord, Coord localCoord ) 
{
	TrilinearElementType*	self = (TrilinearElementType*)elementType;
	Dimension_Index		dim_I;

	for ( dim_I=0; dim_I < 3; dim_I++ ) {
	}
}
#endif


/*

 - Shape function definitions
 - Local node numbering convention for billinear, trilinear element (xi, eta, zeta)
 - Local coordinate domain spans  -1 <= xi,eta,zeta <= 1

    eta
     |
     |____ xi
    /
   /
 zeta


  eta
   |
3-----2
|  |__|___xi
|     |
0-----1
(zeta = -1 plane)


  eta
   |
7-----6
|  |__|___xi
|     |
4-----5
(zeta = +1 plane)


*/
void _TrilinearElementType_SF_allNodes( void* elementType, const double localCoord[], double* const evaluatedValues ) {
	double xi, eta, zeta;
	
	xi   = localCoord[0];
	eta  = localCoord[1];
	zeta = localCoord[2];	
	
	evaluatedValues[0] = 0.125*( 1.0-xi )*( 1.0-eta )*( 1.0-zeta );
	evaluatedValues[2] = 0.125*( 1.0-xi )*( 1.0+eta )*( 1.0-zeta );
	evaluatedValues[3] = 0.125*( 1.0+xi )*( 1.0+eta )*( 1.0-zeta );
	evaluatedValues[1] = 0.125*( 1.0+xi )*( 1.0-eta )*( 1.0-zeta );
	
	evaluatedValues[4] = 0.125*( 1.0-xi )*( 1.0-eta )*( 1.0+zeta );
	evaluatedValues[6] = 0.125*( 1.0-xi )*( 1.0+eta )*( 1.0+zeta );
	evaluatedValues[7] = 0.125*( 1.0+xi )*( 1.0+eta )*( 1.0+zeta );
	evaluatedValues[5] = 0.125*( 1.0+xi )*( 1.0-eta )*( 1.0+zeta );
}


void _TrilinearElementType_SF_allLocalDerivs_allNodes( void* elementType, const double localCoord[],
		double** const evaluatedDerivatives )
{		
	double xi, eta, zeta;
	
	xi   = localCoord[0];
	eta  = localCoord[1];
	zeta = localCoord[2];	
	
	/* derivatives wrt xi */
	evaluatedDerivatives[0][0] = - 0.125*( 1.0-eta )*( 1.0-zeta );
	evaluatedDerivatives[0][2] = - 0.125*( 1.0+eta )*( 1.0-zeta );
	evaluatedDerivatives[0][3] =   0.125*( 1.0+eta )*( 1.0-zeta );
	evaluatedDerivatives[0][1] =   0.125*( 1.0-eta )*( 1.0-zeta );
	evaluatedDerivatives[0][4] = - 0.125*( 1.0-eta )*( 1.0+zeta );
	evaluatedDerivatives[0][6] = - 0.125*( 1.0+eta )*( 1.0+zeta );
	evaluatedDerivatives[0][7] =   0.125*( 1.0+eta )*( 1.0+zeta );
	evaluatedDerivatives[0][5] =   0.125*( 1.0-eta )*( 1.0+zeta );
	
	/* derivatives wrt eta */	
	evaluatedDerivatives[1][0] = - 0.125*( 1.0-xi )*( 1.0-zeta );
	evaluatedDerivatives[1][2] =   0.125*( 1.0-xi )*( 1.0-zeta );
	evaluatedDerivatives[1][3] =   0.125*( 1.0+xi )*( 1.0-zeta );
	evaluatedDerivatives[1][1] = - 0.125*( 1.0+xi )*( 1.0-zeta );
	evaluatedDerivatives[1][4] = - 0.125*( 1.0-xi )*( 1.0+zeta );
	evaluatedDerivatives[1][6] =   0.125*( 1.0-xi )*( 1.0+zeta );
	evaluatedDerivatives[1][7] =   0.125*( 1.0+xi )*( 1.0+zeta );
	evaluatedDerivatives[1][5] = - 0.125*( 1.0+xi )*( 1.0+zeta );
	
	/* derivatives wrt zeta */		
	evaluatedDerivatives[2][0] = -0.125*( 1.0-xi )*( 1.0-eta );
	evaluatedDerivatives[2][2] = -0.125*( 1.0-xi )*( 1.0+eta );
	evaluatedDerivatives[2][3] = -0.125*( 1.0+xi )*( 1.0+eta );
	evaluatedDerivatives[2][1] = -0.125*( 1.0+xi )*( 1.0-eta );
	evaluatedDerivatives[2][4] =  0.125*( 1.0-xi )*( 1.0-eta );
	evaluatedDerivatives[2][6] =  0.125*( 1.0-xi )*( 1.0+eta );
	evaluatedDerivatives[2][7] =  0.125*( 1.0+xi )*( 1.0+eta );
	evaluatedDerivatives[2][5] =  0.125*( 1.0+xi )*( 1.0-eta );
}


void _TrilinearElementType_ConvertGlobalCoordToElLocal(
		void*		elementType,
		void*		_mesh, 
		unsigned	element, 
		const double*	globalCoord,
		double*		elLocalCoord )
{
	TrilinearElementType*	self = (TrilinearElementType*)elementType;
	Mesh*			mesh = (Mesh*)_mesh;
	unsigned		inside;
	double			bc[4];
	static double		lCrds[8][3] = {{-1.0, -1.0, -1.0}, {1.0, -1.0, -1.0}, 
					       {-1.0, 1.0, -1.0}, {1.0, 1.0, -1.0}, 
					       {-1.0, -1.0, 1.0}, {1.0, -1.0, 1.0}, 
					       {-1.0, 1.0, 1.0}, {1.0, 1.0, 1.0}};
	unsigned		nInc, *inc;
	unsigned		bc_i;

	Mesh_GetIncidence( mesh, MT_VOLUME, element, MT_VERTEX, self->inc );
	nInc = IArray_GetSize( self->inc );
	inc = IArray_GetPtr( self->inc );
	assert( nInc == 8 );

	insist( Simplex_Search3D( mesh->verts, inc, 10, self->tetInds, (double*)globalCoord, bc, &inside ), == True );

	elLocalCoord[0] = bc[0] * lCrds[self->tetInds[inside][0]][0];
	elLocalCoord[1] = bc[0] * lCrds[self->tetInds[inside][0]][1];
	elLocalCoord[2] = bc[0] * lCrds[self->tetInds[inside][0]][2];
	for( bc_i = 1; bc_i < 4; bc_i++ ) {
		elLocalCoord[0] += bc[bc_i] * lCrds[self->tetInds[inside][bc_i]][0];
		elLocalCoord[1] += bc[bc_i] * lCrds[self->tetInds[inside][bc_i]][1];
		elLocalCoord[2] += bc[bc_i] * lCrds[self->tetInds[inside][bc_i]][2];
	}
}
