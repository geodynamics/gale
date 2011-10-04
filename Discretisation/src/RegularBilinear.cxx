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
** $Id: RegularBilinear.c 3584 2006-05-16 11:11:07Z PatrickSunter $
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
const Type RegularBilinear_Type = "RegularBilinear";

/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/
#define REGULARBILINEAR_NODECOUNT 4

RegularBilinear* RegularBilinear_New( Name name ) {
	/* Variables set in this function */
	SizeT                                                                            _sizeOfSelf = sizeof(RegularBilinear);
	Type                                                                                    type = RegularBilinear_Type;
	Stg_Class_DeleteFunction*                                                            _delete = _RegularBilinear_Delete;
	Stg_Class_PrintFunction*                                                              _print = _RegularBilinear_Print;
	Stg_Class_CopyFunction*                                                                _copy = NULL;
	Stg_Component_DefaultConstructorFunction*                                _defaultConstructor = (void* (*)(Name))_RegularBilinear_New;
	Stg_Component_ConstructFunction*                                                  _construct = _RegularBilinear_AssignFromXML;
	Stg_Component_BuildFunction*                                                          _build = _RegularBilinear_Build;
	Stg_Component_InitialiseFunction*                                                _initialise = _RegularBilinear_Initialise;
	Stg_Component_ExecuteFunction*                                                      _execute = _RegularBilinear_Execute;
	Stg_Component_DestroyFunction*                                                      _destroy = _RegularBilinear_Destroy;
	AllocationType                                                            nameAllocationType = NON_GLOBAL;
	ElementType_EvaluateShapeFunctionsAtFunction*                      _evaluateShapeFunctionsAt = _BilinearElementType_SF_allNodes;
	ElementType_EvaluateShapeFunctionLocalDerivsAtFunction*  _evaluateShapeFunctionLocalDerivsAt = _BilinearElementType_SF_allLocalDerivs_allNodes;
	ElementType_ConvertGlobalCoordToElLocalFunction*                _convertGlobalCoordToElLocal = _ElementType_ConvertGlobalCoordToElLocal;
	ElementType_JacobianDeterminantSurfaceFunction*                  _jacobianDeterminantSurface = _BilinearElementType_JacobianDeterminantSurface;
	ElementType_SurfaceNormalFunction*                                            _surfaceNormal = _ElementType_SurfaceNormal;

	return _RegularBilinear_New(  REGULARBILINEAR_PASSARGS  );
}

RegularBilinear* _RegularBilinear_New(  REGULARBILINEAR_DEFARGS  ) {
	RegularBilinear*	self;

	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(RegularBilinear) );
	self = (RegularBilinear*)_BilinearElementType_New(  BILINEARELEMENTTYPE_PASSARGS  );

	/* Virtual info */

	/* RegularBilinear info */
	self->isConstructed = True;
	_ElementType_Init( (ElementType*)self, REGULARBILINEAR_NODECOUNT );
   _BilinearElementType_Init( (BilinearElementType*)self );
	_RegularBilinear_Init( self );

	return self;
}

void _RegularBilinear_Init( RegularBilinear* self ) {
	assert( self && Stg_CheckType( self, RegularBilinear ) );
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _RegularBilinear_Delete( void* elementType ) {
	RegularBilinear* self = (RegularBilinear*)elementType;

	/* Delete the parent. */
	_BilinearElementType_Delete( self );
}

void _RegularBilinear_Print( void* elementType, Stream* stream ) {
	RegularBilinear*	self = (RegularBilinear*)elementType;
	
	/* Print parent */
	Journal_Printf( stream, "RegularBilinear (ptr): (%p)\n", self );
	_BilinearElementType_Print( self, stream );
}

void _RegularBilinear_AssignFromXML( void* elementType, Stg_ComponentFactory* cf, void* data ) {
}

void _RegularBilinear_Build( void* elementType, void* data ) {
	_BilinearElementType_Build( elementType, data );
}

void _RegularBilinear_Initialise( void* elementType, void* data ) {
	_BilinearElementType_Initialise( elementType, data );
}

void _RegularBilinear_Execute( void* elementType, void* data ) {
}

void _RegularBilinear_Destroy( void* elementType, void* data ) {
	RegularBilinear* self = (RegularBilinear*)elementType;

	_BilinearElementType_Destroy( self, data );
}


/*--------------------------------------------------------------------------------------------------------------------------
** Public Functions
*/

#if 0
void RegularBilinear_ConvertGlobalCoordToElLocal( void* elementType, void* mesh, unsigned element, 
						  const double* globalCoord, double* localCoord )
{
	RegularBilinear*	self = (RegularBilinear*)elementType;
	unsigned		nInc, *inc;
	double*			vert[2];
	double			w;

	assert( self && Stg_CheckType( self, RegularBilinear ) );

	Mesh_GetIncidence( mesh, MT_VOLUME, element, MT_VERTEX, self->inc );
	nInc = IArray_GetSize( self->inc );
	inc = IArray_GetSize( self->inc );

	vert[0] = Mesh_GetVertex( mesh, inc[0] );
	vert[1] = Mesh_GetVertex( mesh, inc[3] );
	w = vert[1][0] - vert[0][0];
	localCoord[0] = 2.0 * (globalCoord[0] - vert[0][0]) / w - 1.0;
	w = vert[1][1] - vert[0][1];
	localCoord[1] = 2.0 * (globalCoord[1] - vert[0][1]) / w - 1.0;

	assert( localCoord[0] >= -1.0 && localCoord[0] <= 1.0 );
	assert( localCoord[1] >= -1.0 && localCoord[1] <= 1.0 );
}
#endif


/*----------------------------------------------------------------------------------------------------------------------------------
** Private Functions
*/




