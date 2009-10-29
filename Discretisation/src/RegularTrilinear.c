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
** $Id: RegularTrilinear.c 3584 2006-05-16 11:11:07Z PatrickSunter $
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
const Type RegularTrilinear_Type = "RegularTrilinear";


/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/

RegularTrilinear* RegularTrilinear_New( Name name ) {
	return _RegularTrilinear_New( sizeof(RegularTrilinear), 
				      RegularTrilinear_Type, 
				      _RegularTrilinear_Delete, 
				      _RegularTrilinear_Print, 
				      NULL, 
				      (void* (*)(Name))_RegularTrilinear_New, 
				      _RegularTrilinear_AssignFromXML, 
				      _RegularTrilinear_Build, 
				      _RegularTrilinear_Initialise, 
				      _RegularTrilinear_Execute, 
				      _RegularTrilinear_Destroy, 
				      name, 
				      True, 
				      _TrilinearElementType_SF_allNodes, 
				      _TrilinearElementType_SF_allLocalDerivs_allNodes, 
				      _ElementType_ConvertGlobalCoordToElLocal, 
				      _TrilinearElementType_JacobianDeterminantSurface,
				      _ElementType_SurfaceNormal,
				      8 );
}

RegularTrilinear* _RegularTrilinear_New( REGULARTRILINEAR_DEFARGS ) {
	RegularTrilinear*	self;

	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(RegularTrilinear) );
	self = (RegularTrilinear*)_TrilinearElementType_New( ELEMENTTYPE_PASSARGS );

	/* Virtual info */

	/* RegularTrilinear info */
	_RegularTrilinear_Init( self );

	return self;
}

void _RegularTrilinear_Init( RegularTrilinear* self ) {
	assert( self && Stg_CheckType( self, RegularTrilinear ) );
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _RegularTrilinear_Delete( void* elementType ) {
	RegularTrilinear*	self = (RegularTrilinear*)elementType;

	/* Delete the parent. */
	_TrilinearElementType_Delete( self );
}

void _RegularTrilinear_Print( void* elementType, Stream* stream ) {
	RegularTrilinear*	self = (RegularTrilinear*)elementType;
	
	/* Set the Journal for printing informations */
	Stream* elementTypeStream;
	elementTypeStream = Journal_Register( InfoStream_Type, "RegularTrilinearStream" );

	/* Print parent */
	Journal_Printf( stream, "RegularTrilinear (ptr): (%p)\n", self );
	_TrilinearElementType_Print( self, stream );
}

void _RegularTrilinear_AssignFromXML( void* elementType, Stg_ComponentFactory* cf, void* data ) {
}

void _RegularTrilinear_Build( void* elementType, void* data ) {
	_TrilinearElementType_Build( elementType, data );
}

void _RegularTrilinear_Initialise( void* elementType, void* data ) {
	_TrilinearElementType_Initialise( elementType, data );
}

void _RegularTrilinear_Execute( void* elementType, void* data ) {
}

void _RegularTrilinear_Destroy( void* elementType, void* data ) {
}


/*--------------------------------------------------------------------------------------------------------------------------
** Public Functions
*/

#if 0
void RegularTrilinear_ConvertGlobalCoordToElLocal( void* elementType, void* mesh, unsigned element, 
						   const double* globalCoord, double* localCoord )
{
	RegularTrilinear*	self = (RegularTrilinear*)elementType;
	unsigned		nInc, *inc;
	double*			vert[2];
	double			w;

	assert( self && Stg_CheckType( self, RegularTrilinear ) );

	Mesh_GetIncidence( mesh, MT_VOLUME, element, MT_VERTEX, self->inc );
	nInc = IArray_GetSize( self->inc );
	inc = IArray_GetSize( self->inc );

	vert[0] = Mesh_GetVertex( mesh, inc[0] );
	vert[1] = Mesh_GetVertex( mesh, inc[7] );
	w = vert[1][0] - vert[0][0];
	localCoord[0] = 2.0 * (globalCoord[0] - vert[0][0]) / w - 1.0;
	w = vert[1][1] - vert[0][1];
	localCoord[1] = 2.0 * (globalCoord[1] - vert[0][1]) / w - 1.0;
	w = vert[1][2] - vert[0][2];
	localCoord[2] = 2.0 * (globalCoord[2] - vert[0][2]) / w - 1.0;

	assert( Num_InRange( localCoord[0], -1.0, 1.0 ) );
	assert( Num_InRange( localCoord[1], -1.0, 1.0 ) );
	assert( Num_InRange( localCoord[2], -1.0, 1.0 ) );
}
#endif


/*----------------------------------------------------------------------------------------------------------------------------------
** Private Functions
*/
