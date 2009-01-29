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
** $Id: P1.c 3584 2006-05-16 11:11:07Z PatrickSunter $
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
const Type P1_Type = "P1";


/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/

P1* P1_New( Name name ) {
	return _P1_New( sizeof(P1), 
			P1_Type, 
			_P1_Delete, 
			_P1_Print, 
			NULL, 
			(void* (*)(Name))_P1_New, 
			_P1_Construct, 
			_P1_Build, 
			_P1_Initialise, 
			_P1_Execute, 
			_P1_Destroy, 
			name, 
			True, 
			P1_EvalBasis, 
			P1_EvalLocalDerivs, 
			NULL/*P1_CoordGlobalToLocal*/, 
			_ElementType_JacobianDeterminantSurface,
			_P1_SurfaceNormal,
			3 );
}

P1* _P1_New( P1_DEFARGS ) {
	P1*	self;

	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(P1) );
	self = (P1*)_ElementType_New( ELEMENTTYPE_PASSARGS );

	/* Virtual info */

	/* P1 info */
	_P1_Init( self );

	return self;
}

void _P1_Init( P1* self ) {
	assert( self && Stg_CheckType( self, P1 ) );
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _P1_Delete( void* elementType ) {
	P1*	self = (P1*)elementType;

	/* Delete the parent. */
	_ElementType_Delete( self );
}

void _P1_Print( void* elementType, Stream* stream ) {
	P1*	self = (P1*)elementType;
	
	/* Set the Journal for printing informations */
	Stream* elementTypeStream;
	elementTypeStream = Journal_Register( InfoStream_Type, "P1Stream" );

	/* Print parent */
	Journal_Printf( stream, "P1 (ptr): (%p)\n", self );
	_ElementType_Print( self, stream );
}

void _P1_Construct( void* elementType, Stg_ComponentFactory* cf, void* data ) {
}

void _P1_Build( void* elementType, void* data ) {
}

void _P1_Initialise( void* elementType, void* data ) {
}

void _P1_Execute( void* elementType, void* data ) {
}

void _P1_Destroy( void* elementType, void* data ) {
}


/*--------------------------------------------------------------------------------------------------------------------------
** Public Functions
*/

void P1_EvalBasis( void* elementType, const double* localCoord, double* basis ) {
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

void P1_EvalLocalDerivs( void* elementType, const double* localCoord, double** derivs ) {
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

int _P1_SurfaceNormal( void* elementType, unsigned element_I, unsigned dim, double* xi, double* normal ) {
	Stream* errStream = Journal_Register( ErrorStream_Type, ElementType_Type );

	Journal_Printf( errStream, "surface normal function not implemented for this element type.\n" );
	assert( 0 );

	normal = NULL;

	return -1;
}

/*----------------------------------------------------------------------------------------------------------------------------------
** Private Functions
*/
