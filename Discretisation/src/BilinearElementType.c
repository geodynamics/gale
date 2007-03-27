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
** $Id: BilinearElementType.c 659 2006-10-26 02:03:34Z KathleenHumble $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>

#include "units.h"
#include "types.h"
#include "shortcuts.h"
#include "ElementType.h"
#include "BilinearElementType.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

const Type BilinearElementType_Type = "BilinearElementType";

#define _BilinearElementType_NodeCount 4

void* BilinearElementType_DefaultNew( Name name ) {
	return _BilinearElementType_New(
			sizeof(BilinearElementType), 
			BilinearElementType_Type,
			_BilinearElementType_Delete,
			_BilinearElementType_Print,
			NULL, 
			BilinearElementType_DefaultNew,
			_BilinearElementType_Construct,
			_BilinearElementType_Build,
			_BilinearElementType_Initialise,
			_BilinearElementType_Execute,
			_BilinearElementType_Destroy, 
			name, 
			False,
			_BilinearElementType_SF_allNodes,
			_BilinearElementType_SF_allLocalDerivs_allNodes,
			_BilinearElementType_ConvertGlobalCoordToElLocal,
			_BilinearElementType_NodeCount );
}

BilinearElementType* BilinearElementType_New( Name name ) {
	return _BilinearElementType_New( sizeof(BilinearElementType), BilinearElementType_Type, _BilinearElementType_Delete,
		_BilinearElementType_Print, NULL, BilinearElementType_DefaultNew, _BilinearElementType_Construct, _BilinearElementType_Build,
		_BilinearElementType_Initialise, _BilinearElementType_Execute, _BilinearElementType_Destroy, name, True, _BilinearElementType_SF_allNodes, 
		_BilinearElementType_SF_allLocalDerivs_allNodes, _BilinearElementType_ConvertGlobalCoordToElLocal,
		_BilinearElementType_NodeCount );
}


BilinearElementType* _BilinearElementType_New( 
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
	BilinearElementType*		self;
	
	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(BilinearElementType) );
	self = (BilinearElementType*)_ElementType_New( _sizeOfSelf, type, _delete, _print, _copy, _defaultConstructor,
			_construct, _build, _initialise, _execute, _destroy, name, initFlag,
		_evaluateShapeFunctionsAt, _evaluateShapeFunctionLocalDerivsAt, _convertGlobalCoordToElLocal,
		nodeCount );
	
	/* General info */
	
	/* Virtual functions */
	
	/* BilinearElementType info */
	if( initFlag ){
		_BilinearElementType_Init( self );
	}
	
	return self;
}



void _BilinearElementType_Init( BilinearElementType* self ) {
	Dimension_Index dim_I=0;
	/* General and Virtual info should already be set */
	
	/* BilinearElementType info */
	self->isConstructed = True;
	for ( dim_I = 0; dim_I < 2; dim_I++ ) {
		self->minElLocalCoord[dim_I] = -1;
		self->maxElLocalCoord[dim_I] = 1;
		self->elLocalLength[dim_I] = self->maxElLocalCoord[dim_I] - self->minElLocalCoord[dim_I];
	}
}

void _BilinearElementType_Delete( void* elementType ) {
	BilinearElementType* self = (BilinearElementType*)elementType;
	
	Journal_DPrintf( self->debug, "In %s\n", __func__ );
	/* Stg_Class_Delete parent*/
	_ElementType_Delete( self );
}


void _BilinearElementType_Print( void* elementType, Stream* stream ) {
	BilinearElementType* self = (BilinearElementType*)elementType;
	Dimension_Index dim_I=0;

	/* General info */
	Journal_Printf( stream, "BilinearElementType (ptr): %p\n", self );
	
	/* Print parent */
	_ElementType_Print( self, stream );
	
	/* Virtual info */
	
	/* BilinearElementType info */
	Journal_Printf( stream, "self->minElLocalCoord: (", self );
	for ( dim_I = 0; dim_I < 2; dim_I++ ) {
		Journal_Printf( stream, "%0.5f,", self->minElLocalCoord[dim_I] );
	}
	Journal_Printf( stream, ")\n", self );
	Journal_Printf( stream, "self->maxElLocalCoord: (", self );
	for ( dim_I = 0; dim_I < 2; dim_I++ ) {
		Journal_Printf( stream, "%0.5f,", self->maxElLocalCoord[dim_I] );
	}	
	Journal_Printf( stream, ")\n", self );
	Journal_Printf( stream, "self->elLocalLength: (", self );
	for ( dim_I = 0; dim_I < 2; dim_I++ ) {
		Journal_Printf( stream, "%0.5f,", self->elLocalLength[dim_I] );
	}
	Journal_Printf( stream, ")\n", self );
}

void _BilinearElementType_Construct( void* elementType, Stg_ComponentFactory *cf, void* data ){
	
}
	
void _BilinearElementType_Initialise( void* elementType, void *data ){
	
}
	
void _BilinearElementType_Execute( void* elementType, void *data ){
	
}
	
void _BilinearElementType_Destroy( void* elementType, void *data ){
	
}

void _BilinearElementType_Build( void* elementType, void *data ) {

}

/*

 - Shape function definitions
 - Local node numbering convention for billinear element (xi, eta)
 - Local coordinate domain spans  -1 <= xi,eta <= 1

  eta
   |
3-----2
|  |__|___xi
|     |
0-----1

*/
void _BilinearElementType_SF_allNodes( void* elementType, const double localCoord[], double* const evaluatedValues ) {
	double xi, eta;
	
	xi  = localCoord[0];
	eta = localCoord[1];
	
	evaluatedValues[0] = 0.25*( 1.0-xi )*( 1.0-eta );
	evaluatedValues[1] = 0.25*( 1.0+xi )*( 1.0-eta );
	evaluatedValues[2] = 0.25*( 1.0+xi )*( 1.0+eta );
	evaluatedValues[3] = 0.25*( 1.0-xi )*( 1.0+eta );
}


void _BilinearElementType_SF_allLocalDerivs_allNodes( void* elementType, const double localCoord[],
		double** const evaluatedDerivatives ) 
{		
	double xi, eta;
	
	xi  = localCoord[0];
	eta = localCoord[1];
	
	/* derivatives wrt xi */
	evaluatedDerivatives[0][0] = - 0.25*( 1.0 - eta );
	evaluatedDerivatives[0][1] =   0.25*( 1.0 - eta );
	evaluatedDerivatives[0][2] =   0.25*( 1.0 + eta );
	evaluatedDerivatives[0][3] = - 0.25*( 1.0 + eta );
	
	/* derivatives wrt eta */
	evaluatedDerivatives[1][0] = - 0.25*( 1.0 - xi );
	evaluatedDerivatives[1][1] = - 0.25*( 1.0 + xi );
	evaluatedDerivatives[1][2] =   0.25*( 1.0 + xi );
	evaluatedDerivatives[1][3] =   0.25*( 1.0 - xi );
}


/*
** Calculates the barycenter of a triangle with respect to some point.
*/

#if 0
void _BilinearElementType_TriBarycenter( Coord tri[3], const Coord pnt, Coord dst ) {
	double	a = tri[1][0] - tri[0][0];
	double	b = tri[2][0] - tri[0][0];
	double	c = tri[1][1] - tri[0][1];
	double	d = tri[2][1] - tri[0][1];
	double	e = pnt[0] - tri[0][0];
	double	f = pnt[1] - tri[0][1];
	double	aInv = 1.0 / a;

	dst[2] = (f - (c * e * aInv)) / (d - (c * b * aInv));
	dst[1] = (e - b * dst[1]) * aInv;
	dst[0] = 1.0 - dst[1] - dst[2];
}
#endif

void _BilinearElementType_ConvertGlobalCoordToElLocal(
		void*		elementType,
		ElementLayout*	elementLayout,
		const Coord**	globalNodeCoordPtrsInElement,
		const Coord	globalCoord,
		Coord		elLocalCoord )
{
	BilinearElementType*	self = (BilinearElementType*)elementType;
	Dimension_Index		dim_I = 0;
	double			globalToElLocalScaling[2] = {0.0,0.0};
	Coord			relToBottomLeftGlobalCoord = {0.0,0.0,0.0};

	if ( elementLayout->type == ParallelPipedHexaEL_Type ) {
		double	elLen[2];
		elLen[0] = (*(globalNodeCoordPtrsInElement[1]))[0] - (*(globalNodeCoordPtrsInElement[0]))[0];
		elLen[1] = (*(globalNodeCoordPtrsInElement[3]))[1] - (*(globalNodeCoordPtrsInElement[0]))[1];
		
		/* Initially set elLocalCoord to (0,0,0) */
		memset( elLocalCoord, 0, sizeof( Coord ) );
	
		for( dim_I=0; dim_I < 2; dim_I++ ) {
			globalToElLocalScaling[dim_I] = self->elLocalLength[dim_I] / elLen[dim_I];

			/* The bottom left node is always at index zero */
			relToBottomLeftGlobalCoord[dim_I] = globalCoord[dim_I] - (*(globalNodeCoordPtrsInElement[0]))[dim_I];
			
			elLocalCoord[dim_I] =  
				self->minElLocalCoord[dim_I] + relToBottomLeftGlobalCoord[dim_I] * globalToElLocalScaling[dim_I];
			assert( elLocalCoord[dim_I] >= -1.0 || elLocalCoord[dim_I] <= 1.0 );
		}
	}
	else if( elementLayout->type == HexaEL_Type ) {
		Coord		crds[4];
		Coord		bc;
		unsigned	inds[3];
		Coord		lCrds[4] = {{-1.0, -1.0, 0.0}, {1.0, -1.0, 0.0}, 
					    {-1.0, 1.0, 0.0}, {1.0, 1.0, 0.0}};
		/*unsigned	gElInd;*/
		unsigned	bc_i;

		/* Check first triangle. */
		memcpy( crds[0], *(globalNodeCoordPtrsInElement[0]), sizeof(Coord) );
		memcpy( crds[1], *(globalNodeCoordPtrsInElement[1]), sizeof(Coord) );
		memcpy( crds[2], *(globalNodeCoordPtrsInElement[3]), sizeof(Coord) );
		memcpy( crds[3], *(globalNodeCoordPtrsInElement[2]), sizeof(Coord) );
#ifndef NDEBUG
		if( !_HexaEL_FindTriBarycenter( (const Coord*)crds, globalCoord, bc, inds, INCLUSIVE_UPPER_BOUNDARY, NULL, 0 ) )
			assert( 0 );
#else
		_HexaEL_FindTriBarycenter( crds, globalCoord, bc, inds, INCLUSIVE_UPPER_BOUNDARY, NULL, 0 );
#endif

		/* Interpolate. */
		memset( elLocalCoord, 0, sizeof(Coord) );
		for( bc_i = 0; bc_i < 3; bc_i++ ) {
			elLocalCoord[0] += bc[bc_i] * lCrds[inds[bc_i]][0];
			elLocalCoord[1] += bc[bc_i] * lCrds[inds[bc_i]][1];
		}
	}
	else {
		/* Not a rectangular element -> Just use the general version */
		_ElementType_ConvertGlobalCoordToElLocal( self, elementLayout, globalNodeCoordPtrsInElement,
			globalCoord, elLocalCoord );
	}
}
