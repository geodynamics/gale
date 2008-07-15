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
** $Id: ConstantElementType.c 1177 2008-07-15 01:29:58Z DavidLee $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include "units.h"
#include "types.h"
#include "shortcuts.h"
#include "ElementType.h"
#include "ConstantElementType.h"
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

const Type ConstantElementType_Type = "ConstantElementType";
#define _ConstantElementType_NodeCount 1

void* ConstantElementType_DefaultNew( Name name ) {
	return _ConstantElementType_New( 
			sizeof(ConstantElementType),
			ConstantElementType_Type,
			_ConstantElementType_Delete,
			_ConstantElementType_Print,
			NULL, 
			ConstantElementType_DefaultNew,
			_ConstantElementType_Construct,
			_ConstantElementType_Build,
			_ConstantElementType_Initialise,
			_ConstantElementType_Execute,
			_ConstantElementType_Destroy, 
			name,
			False,
			_ConstantElementType_SF_allNodes,
			_ConstantElementType_SF_allLocalDerivs_allNodes,
			_ConstantElementType_ConvertGlobalCoordToElLocal,
			_ElementType_JacobianDeterminantSurface,
			_ConstantElementType_NodeCount );
}

ConstantElementType* ConstantElementType_New( Name name ) {
	return _ConstantElementType_New( 
			sizeof(ConstantElementType),
			ConstantElementType_Type,
			_ConstantElementType_Delete,
			_ConstantElementType_Print, 
			NULL, 
			ConstantElementType_DefaultNew,
			_ConstantElementType_Construct,
			_ConstantElementType_Build,
			_ConstantElementType_Initialise,
			_ConstantElementType_Execute,
			_ConstantElementType_Destroy,
			name, 
			True,
			_ConstantElementType_SF_allNodes,
			_ConstantElementType_SF_allLocalDerivs_allNodes,
			_ConstantElementType_ConvertGlobalCoordToElLocal,
			_ElementType_JacobianDeterminantSurface,
			_ConstantElementType_NodeCount );
}


ConstantElementType* _ConstantElementType_New( 
		SizeT								_sizeOfSelf,
		Type								type,
		Stg_Class_DeleteFunction*					_delete,
		Stg_Class_PrintFunction*					_print,
		Stg_Class_CopyFunction*						_copy, 
		Stg_Component_DefaultConstructorFunction*			_defaultConstructor,
		Stg_Component_ConstructFunction*				_construct,
		Stg_Component_BuildFunction*					_build,
		Stg_Component_InitialiseFunction*				_initialise,
		Stg_Component_ExecuteFunction*					_execute,
		Stg_Component_DestroyFunction*					_destroy,
		Name								name,
		Bool								initFlag,
		ElementType_EvaluateShapeFunctionsAtFunction*			_evaluateShapeFunctionsAt,
		ElementType_EvaluateShapeFunctionLocalDerivsAtFunction*		_evaluateShapeFunctionLocalDerivsAt,
		ElementType_ConvertGlobalCoordToElLocalFunction*		_convertGlobalCoordToElLocal,
		ElementType_JacobianDeterminantSurfaceFunction*			_jacobianDeterminantSurface,
		Index								nodeCount )
{
	ConstantElementType*		self;
	
	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(ConstantElementType) );
	self = (ConstantElementType*)_ElementType_New( _sizeOfSelf, type, _delete, _print, _copy,
		_defaultConstructor, _construct, _build, _initialise, _execute, _destroy, name, initFlag,
		_evaluateShapeFunctionsAt, _evaluateShapeFunctionLocalDerivsAt, _convertGlobalCoordToElLocal,
		_jacobianDeterminantSurface, nodeCount );
	
	/* General info */
	
	/* Virtual functions */
	
	/* ConstantElementType info */
	if( initFlag ){
		_ConstantElementType_Init( self );
	}
	
	return self;
}


void _ConstantElementType_Init( ConstantElementType* self ) {
	/* General and Virtual info should already be set */
	self->isConstructed = True;
	/* ConstantElementType info */
}

void _ConstantElementType_Delete( void* elementType ) {
	ConstantElementType* self = (ConstantElementType*)elementType;
	Journal_DPrintf( self->debug, "In %s\n", __func__ );

	/* Stg_Class_Delete parent*/
	_ElementType_Delete( self );
}

void _ConstantElementType_Print( void* elementType, Stream* stream ) {
	ConstantElementType* self = (ConstantElementType*)elementType;
	
	/* Set the Journal for printing informations */
	Stream* constantElementTypeStream = stream;
	
	/* General info */
	Journal_Printf( constantElementTypeStream, "ConstantElementType (ptr): %p\n", self );
	
	/* Print parent */
	_ElementType_Print( self, constantElementTypeStream );
	
	/* Virtual info */
	
	/* ConstantElementType info */
}

void _ConstantElementType_Construct( void* elementType, Stg_ComponentFactory *cf, void* data ){
	
}
	
void _ConstantElementType_Initialise( void* elementType, void *data ){
	
}
	
void _ConstantElementType_Execute( void* elementType, void *data ){
	
}
	
void _ConstantElementType_Destroy( void* elementType, void *data ){
	
}

void _ConstantElementType_Build( void* elementType, void *data ) {
	/* NOTHING */

}

/*

 - Shape function definition
 - Constant shape function in 2d/3d has only one node at centroid of element.
 - Node id is 0 in both cases.
 - Local coordinate domain spans  -1 <= xi,eta <= 1 in 2d
 - Local coordinate domain spans  -1 <= xi,eta,zeta <= 1 in 3d

*/
void _ConstantElementType_SF_allNodes( void* elementType, const double localCoord[], double* const evaluatedValues ) {
	evaluatedValues[0] = 1.0;
}


/*
Since we use only have one constant shape func for 2d and 3d quads, then
if we want to return the zero derivatives we need to pass in "dim" so we know
whether to fill in GNi[0][0] = GNi[1][0] = 0.0 + GNi[2][0] = 0.0 if dim == 3
Should just return error if we try to take deriv of constant. No one would want
to do this!
*/
void _ConstantElementType_SF_allLocalDerivs_allNodes( void* elementType, const double localCoord[],
		double** const evaluatedDerivatives )
{		
	Stream* error = Journal_Register( ErrorStream_Type, ConstantElementType_Type );
	/* all derivatives are zero*/
	Journal_Printf( error, "Error: Trying take derivative of a constant shape function \n" );
	assert( 0 );
}


void _ConstantElementType_ConvertGlobalCoordToElLocal(
		void*		elementType,
		Mesh*		mesh, 
		unsigned	element, 
		const double*	globalCoord,
		double*		elLocalCoord )
{
	/* See header file function introduction for explanation... */
	elLocalCoord[0] = elLocalCoord[1] = elLocalCoord[2] = 0;
}
