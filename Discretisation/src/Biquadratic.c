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
** $Id: Biquadratic.c 3584 2006-05-16 11:11:07Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <mpi.h>

#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include "Discretisation.h"


/* Textual name of this class */
const Type Biquadratic_Type = "Biquadratic";


/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/

#define BIQUADRATICNODECOUNT 9

Biquadratic* Biquadratic_New( Name name ) {
	return _Biquadratic_New( sizeof(Biquadratic), 
				 Biquadratic_Type, 
				 _Biquadratic_Delete, 
				 _Biquadratic_Print, 
				 NULL, 
				 (void* (*)(Name))_Biquadratic_New, 
				 _Biquadratic_AssignFromXML, 
				 _Biquadratic_Build, 
				 _Biquadratic_Initialise, 
				 _Biquadratic_Execute, 
				 _Biquadratic_Destroy, 
				 name, 
				 True, 
				 Biquadratic_EvalBasis, 
				 Biquadratic_EvalLocalDerivs, 
				 _ElementType_ConvertGlobalCoordToElLocal, 
				 Biquadratic_JacobianDeterminantSurface,
				 _ElementType_SurfaceNormal,
				 BIQUADRATICNODECOUNT );
}

Biquadratic* _Biquadratic_New( BIQUADRATIC_DEFARGS ) {
	Biquadratic*	self;

	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(Biquadratic) );
	self = (Biquadratic*)_ElementType_New( ELEMENTTYPE_PASSARGS );

	/* Virtual info */

	/* Biquadratic info */
	_Biquadratic_Init( self );

	return self;
}

void _Biquadratic_Init( Biquadratic* self ) {
	assert( self && Stg_CheckType( self, Biquadratic ) );

	self->dim = 2;
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _Biquadratic_Delete( void* elementType ) {
	Biquadratic*	self = (Biquadratic*)elementType;

	/* Delete the parent. */
	_ElementType_Delete( self );
}

void _Biquadratic_Print( void* elementType, Stream* stream ) {
	Biquadratic*	self = (Biquadratic*)elementType;
	
	/* Set the Journal for printing informations */
	Stream* elementTypeStream;
	elementTypeStream = Journal_Register( InfoStream_Type, "BiquadraticStream" );

	/* Print parent */
	Journal_Printf( stream, "Biquadratic (ptr): (%p)\n", self );
	_ElementType_Print( self, stream );
}

void _Biquadratic_AssignFromXML( void* elementType, Stg_ComponentFactory* cf, void* data ) {
}

void _Biquadratic_Build( void* elementType, void* data ) {
}

void _Biquadratic_Initialise( void* elementType, void* data ) {
	Biquadratic*	self = (Biquadratic*)elementType;

	self->faceNodes = Memory_Alloc_2DArray( unsigned, 4, 3, "node indices for element faces" );

	self->faceNodes[0][0] = 0; self->faceNodes[0][1] = 1; self->faceNodes[0][2] = 2;
	self->faceNodes[1][0] = 6; self->faceNodes[1][1] = 7; self->faceNodes[1][2] = 8;
	self->faceNodes[2][0] = 0; self->faceNodes[2][1] = 3; self->faceNodes[2][2] = 6;
	self->faceNodes[3][0] = 2; self->faceNodes[3][1] = 5; self->faceNodes[3][2] = 8;

	self->evaluatedShapeFunc = Memory_Alloc_Array( double, self->nodeCount, "evaluatedShapeFuncs" );
	self->GNi = Memory_Alloc_2DArray( double, self->dim, self->nodeCount, "localShapeFuncDerivatives" );
}

void _Biquadratic_Execute( void* elementType, void* data ) {
}

void _Biquadratic_Destroy( void* elementType, void* data ) {
	Biquadratic*	self = (Biquadratic*)elementType;

	Memory_Free( self->faceNodes );
	Memory_Free( self->evaluatedShapeFunc );
	Memory_Free( self->GNi );

	_ElementType_Destroy( elementType, data );
}


/*--------------------------------------------------------------------------------------------------------------------------
** Public Functions
*/

void Biquadratic_EvalBasis( void* elementType, const double* localCoord, double* basis ) {
	double	xi = localCoord[0], eta = localCoord[1];
	double	a0 = xi - 1.0, b0 = eta - 1.0;
	double	a1 = 1.0 - xi * xi, b1 = 1.0 - eta * eta;
	double	a2 = xi + 1.0, b2 = eta + 1.0;
	double	m0 = 0.5 * xi;
	double	m1 = 0.5 * eta;
	double	m2 = 0.25 * xi * eta;

	basis[0] = m2 * a0 * b0;
	basis[1] = m1 * a1 * b0;
	basis[2] = m2 * a2 * b0;

	basis[3] = m0 * a0 * b1;
	basis[4] = a1 * b1;
	basis[5] = m0 * a2 * b1;

	basis[6] = m2 * a0 * b2;
	basis[7] = m1 * a1 * b2;
	basis[8] = m2 * a2 * b2;
}

void Biquadratic_EvalLocalDerivs( void* elementType, const double* localCoord, double** derivs ) {
	double	xi = localCoord[0], eta = localCoord[1];
	double	a0 = xi - 1.0, b0 = eta - 1.0;
	double	a1 = xi + 1.0, b1 = eta + 1.0;
	double	a2 = 2.0 * xi - 1.0, b2 = 2.0 * eta - 1.0;
	double	a3 = 2.0 * xi + 1.0, b3 = 2.0 * eta + 1.0;
	double	a4 = 1.0 - xi * xi, b4 = 1.0 - eta * eta;
	double	m0 = 0.25 * xi;
	double	m1 = 0.25 * eta;
	double	m2 = -xi * eta;

	/* Corner nodes. */
	derivs[0][0] = m1 * a2 * b0;
	derivs[0][2] = m1 * a3 * b0;
	derivs[0][6] = m1 * a2 * b1;
	derivs[0][8] = m1 * a3 * b1;
	derivs[1][0] = m0 * a0 * b2;
	derivs[1][2] = m0 * a1 * b2;
	derivs[1][6] = m0 * a0 * b3;
	derivs[1][8] = m0 * a1 * b3;

	/* Side nodes. */
	derivs[0][1] = m2 * b0;
	derivs[0][7] = m2 * b1;
	derivs[0][3] = 0.5 * a2 * b4;
	derivs[0][5] = 0.5 * a3 * b4;
	derivs[1][1] = 0.5 * a4 * b2;
	derivs[1][7] = 0.5 * a4 * b3;
	derivs[1][3] = m2 * a0;
	derivs[1][5] = m2 * a1;

	/* Center node. */
	derivs[0][4] = -2.0 * xi * b4;
	derivs[1][4] = -2.0 * eta * a4;
}

double Biquadratic_JacobianDeterminantSurface( void* elementType, void* _mesh, unsigned element_I, const double* localCoord, 
						unsigned face_I, unsigned norm )
{
	Biquadratic*	self		= (Biquadratic*) elementType;
	Mesh*		mesh		= (Mesh*)_mesh;
	unsigned	surfaceDim	= ( norm + 1 ) % 2;
	double		x[3];
	double		detJac;
	unsigned	nodes[3];

	ElementType_GetFaceNodes( elementType, mesh, element_I, face_I, 3, nodes );

	x[0] = Mesh_GetVertex( mesh, nodes[0] )[surfaceDim];
	x[1] = Mesh_GetVertex( mesh, nodes[1] )[surfaceDim];
	x[2] = Mesh_GetVertex( mesh, nodes[2] )[surfaceDim];

	detJac = ( localCoord[surfaceDim] - 0.5 ) * x[0] - 2.0 * localCoord[surfaceDim] * x[1] + 
		 ( localCoord[surfaceDim] + 0.5 ) * x[2];

	return fabs( detJac );
}

/*----------------------------------------------------------------------------------------------------------------------------------
** Private Functions
*/
