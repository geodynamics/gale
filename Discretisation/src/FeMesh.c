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
** $Id: FeMesh.c 822 2007-04-27 06:20:35Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <mpi.h>

#include "StGermain/StGermain.h"
#include "Discretisation.h"


/* Textual name of this class */
const Type FeMesh_Type = "FeMesh";


/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/

FeMesh* FeMesh_New( Name name ) {
	return _FeMesh_New( sizeof(FeMesh), 
			  FeMesh_Type, 
			  _FeMesh_Delete, 
			  _FeMesh_Print, 
			  NULL, 
			  (void* (*)(Name))_FeMesh_New, 
			  _FeMesh_Construct, 
			  _FeMesh_Build, 
			  _FeMesh_Initialise, 
			  _FeMesh_Execute, 
			  _FeMesh_Destroy, 
			  name, 
			  NON_GLOBAL );
}

FeMesh* _FeMesh_New( FEMESH_DEFARGS ) {
	FeMesh*	self;

	/* Allocate memory */
	assert( sizeOfSelf >= sizeof(FeMesh) );
	self = (FeMesh*)_Mesh_New( MESH_PASSARGS );

	/* Virtual info */

	/* FeMesh info */
	_FeMesh_Init( self );

	return self;
}

void _FeMesh_Init( FeMesh* self ) {
	Stream*	stream;

	assert( self && Stg_CheckType( self, FeMesh ) );

	stream = Journal_Register( Info_Type, self->type );
	Stream_SetPrintingRank( stream, 0 );

	self->feElType = NULL;
	self->feElFamily = NULL;
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _FeMesh_Delete( void* feMesh ) {
	FeMesh*	self = (FeMesh*)feMesh;

	FeMesh_Destruct( self );

	/* Delete the parent. */
	_Mesh_Delete( self );
}

void _FeMesh_Print( void* feMesh, Stream* stream ) {
	FeMesh*	self = (FeMesh*)feMesh;
	
	/* Set the Journal for printing informations */
	Stream* feMeshStream;
	feMeshStream = Journal_Register( InfoStream_Type, "FeMeshStream" );

	/* Print parent */
	Journal_Printf( stream, "FeMesh (ptr): (%p)\n", self );
	_Mesh_Print( self, stream );
}

void _FeMesh_Construct( void* feMesh, Stg_ComponentFactory* cf, void* data ) {
	FeMesh*		self = (FeMesh*)feMesh;
	char*		family;

	assert( self );

	_Mesh_Construct( self, cf, data );

	family = Stg_ComponentFactory_GetString( cf, self->name, "elementType", "linear" );
	FeMesh_SetElementFamily( self, family );
}

void _FeMesh_Build( void* feMesh, void* data ) {
	FeMesh*		self = (FeMesh*)feMesh;
	Stream*		stream;
	ElementType*	elType;

	assert( self );

	stream = Journal_Register( Info_Type, self->type );

	_Mesh_Build( self, data );

	Stream_Indent( stream );
	Journal_Printf( stream, "Assigning FeMesh element types...\n" );
	Stream_Indent( stream );

	if( !strcmp( self->feElFamily, "quadratic" ) ) {
		unsigned	nDims;

		nDims = Mesh_GetDimSize( self );
		if( nDims == 2 )
			elType = (ElementType*)Biquadratic_New( "" );
		else
			abort();
	}
	else if( !strcmp( self->feElFamily, "linear" ) ) {
		unsigned	nDims;

		nDims = Mesh_GetDimSize( self );
		if( nDims == 3 )
			elType = (ElementType*)TrilinearElementType_New( "" );
		else if( nDims == 2 )
			elType = (ElementType*)BilinearElementType_New( "" );
		else
			abort();
	}
	else if( !strcmp( self->feElFamily, "linear-regular" ) ) {
		unsigned	nDims;

		nDims = Mesh_GetDimSize( self );
		if( nDims == 3 )
			elType = (ElementType*)RegularTrilinear_New( "" );
		else if( nDims == 2 )
			elType = (ElementType*)RegularBilinear_New( "" );
		else
			abort();
	}
	else if( !strcmp( self->feElFamily, "constant" ) ) {
		elType = (ElementType*)ConstantElementType_New( "" );
	}
	else if( !strcmp( self->feElFamily, "p1" ) ) {
		elType = (ElementType*)P1_New( "" );
	}
	else
		abort();
	FeMesh_SetElementType( self, elType );
	if( self->feElType )
		Stg_Component_Build( self->feElType, data, False );

	Journal_Printf( stream, "... FE element types are '%s',\n", elType->type );
	Journal_Printf( stream, "... done.\n" );
	Stream_UnIndent( stream );
	Stream_UnIndent( stream );
}

void _FeMesh_Initialise( void* feMesh, void* data ) {
	FeMesh*	self = (FeMesh*)feMesh;

	assert( self );

	_Mesh_Initialise( self, data );

	if( self->feElType )
		Stg_Component_Initialise( self->feElType, data, False );
}

void _FeMesh_Execute( void* feMesh, void* data ) {
}

void _FeMesh_Destroy( void* feMesh, void* data ) {
}


/*--------------------------------------------------------------------------------------------------------------------------
** Public Functions
*/

void FeMesh_SetElementType( void* feMesh, ElementType* elType ) {
	FeMesh*	self = (FeMesh*)feMesh;

	assert( self );

	FreeObject( self->feElType );
	self->feElType = elType;
}

void FeMesh_SetElementFamily( void* feMesh, const char* family ) {
	FeMesh*	self = (FeMesh*)feMesh;

	assert( self );

	self->feElFamily = (char*)family;
}

ElementType* FeMesh_GetElementType( void* feMesh, unsigned element ) {
	FeMesh*	self = (FeMesh*)feMesh;

	assert( self );

	return self->feElType;
}

unsigned FeMesh_GetNodeLocalSize( void* feMesh ) {
	return Mesh_GetLocalSize( feMesh, MT_VERTEX );
}

unsigned FeMesh_GetNodeRemoteSize( void* feMesh ) {
	return Mesh_GetRemoteSize( feMesh, MT_VERTEX );
}

unsigned FeMesh_GetNodeDomainSize( void* feMesh ) {
	return Mesh_GetDomainSize( feMesh, MT_VERTEX );
}

unsigned FeMesh_GetNodeGlobalSize( void* feMesh ) {
	return Mesh_GetGlobalSize( feMesh, MT_VERTEX );
}

unsigned FeMesh_GetElementLocalSize( void* feMesh ) {
	return Mesh_GetLocalSize( feMesh, Mesh_GetDimSize( feMesh ) );
}

unsigned FeMesh_GetElementRemoteSize( void* feMesh ) {
	return Mesh_GetRemoteSize( feMesh, Mesh_GetDimSize( feMesh ) );
}

unsigned FeMesh_GetElementDomainSize( void* feMesh ) {
	return Mesh_GetDomainSize( feMesh, Mesh_GetDimSize( feMesh ) );
}

unsigned FeMesh_GetElementGlobalSize( void* feMesh ) {
	return Mesh_GetGlobalSize( feMesh, Mesh_GetDimSize( feMesh ) );
}

unsigned FeMesh_GetElementNodeSize( void* feMesh, unsigned element ) {
	return Mesh_GetIncidenceSize( feMesh, Mesh_GetDimSize( feMesh ), element, MT_VERTEX );
}

unsigned FeMesh_GetNodeElementSize( void* feMesh, unsigned node ) {
	return Mesh_GetIncidenceSize( feMesh, MT_VERTEX, node, Mesh_GetDimSize( feMesh ) );
}

void FeMesh_GetElementNodes( void* feMesh, unsigned element, unsigned* nNodes, unsigned** nodes ) {
	Mesh_GetIncidence( feMesh, Mesh_GetDimSize( feMesh ), element, MT_VERTEX, nNodes, nodes );
}

void FeMesh_GetNodeElements( void* feMesh, unsigned node, unsigned* nElements, unsigned** elements ) {
	Mesh_GetIncidence( feMesh, MT_VERTEX, node, Mesh_GetDimSize( feMesh ), nElements, elements );
}

unsigned FeMesh_ElementDomainToGlobal( void* feMesh, unsigned domain ) {
	return Mesh_DomainToGlobal( feMesh, Mesh_GetDimSize( feMesh ), domain );
}

Bool FeMesh_ElementGlobalToDomain( void* feMesh, unsigned global, unsigned* domain ) {
	return Mesh_GlobalToDomain( feMesh, Mesh_GetDimSize( feMesh ), global, domain );
}

unsigned FeMesh_NodeDomainToGlobal( void* feMesh, unsigned domain ) {
	return Mesh_DomainToGlobal( feMesh, MT_VERTEX, domain );
}

Bool FeMesh_NodeGlobalToDomain( void* feMesh, unsigned global, unsigned* domain ) {
	return Mesh_GlobalToDomain( feMesh, MT_VERTEX, global, domain );
}

void FeMesh_CoordGlobalToLocal( void* feMesh, unsigned element, double* global, double* local ) {
	FeMesh*		self = (FeMesh*)feMesh;
	ElementType*	elType;

	assert( self );
	assert( element < FeMesh_GetElementDomainSize( self ) );
	assert( global );
	assert( local );

	elType = FeMesh_GetElementType( self, element );
	ElementType_ConvertGlobalCoordToElLocal( elType, self, element, global, local );
}

void FeMesh_CoordLocalToGlobal( void* feMesh, unsigned element, double* local, double* global ) {
	FeMesh*		self = (FeMesh*)feMesh;
	unsigned	nDims;
	ElementType*	elType;
	double*		basis;
	unsigned	nElNodes, *elNodes;
	double		dimBasis;
	double*		vert;
	unsigned	n_i, d_i;

	assert( self );
	assert( element < FeMesh_GetElementDomainSize( self ) );
	assert( global );
	assert( local );

	nDims = Mesh_GetDimSize( self );
	elType = FeMesh_GetElementType( self, element );
	FeMesh_GetElementNodes( self, element, &nElNodes, &elNodes );
	basis = AllocArray( double, nElNodes );
	ElementType_EvaluateShapeFunctionsAt( elType, local, basis );

	memset( global, 0, nDims * sizeof(double) );
	for( n_i = 0; n_i < nElNodes; n_i++ ) {
		dimBasis = basis[n_i];
		vert = Mesh_GetVertex( self, elNodes[n_i] );
		for( d_i = 0; d_i < nDims; d_i++ )
			global[d_i] += dimBasis * vert[d_i];
	}

	FreeArray( basis );
}

void FeMesh_EvalBasis( void* feMesh, unsigned element, double* localCoord, double* basis ) {
	FeMesh*		self = (FeMesh*)feMesh;
	ElementType*	elType;

	assert( self );
	assert( localCoord );

	elType = FeMesh_GetElementType( self, element );
	ElementType_EvaluateShapeFunctionsAt( elType, localCoord, basis );
}

void FeMesh_EvalLocalDerivs( void* feMesh, unsigned element, double* localCoord, double** derivs ) {
	FeMesh*		self = (FeMesh*)feMesh;
	ElementType*	elType;

	assert( self );
	assert( localCoord );
	assert( derivs );

	elType = FeMesh_GetElementType( self, element );
	ElementType_EvaluateShapeFunctionLocalDerivsAt( elType, localCoord, derivs );
}

void FeMesh_EvalGlobalDerivs( void* feMesh, unsigned element, double* localCoord, double** derivs, double* jacDet ) {
	FeMesh*		self = (FeMesh*)feMesh;
	unsigned	nDims;
	ElementType*	elType;
	double		jd;

	assert( self );
	assert( localCoord );
	assert( derivs );

	nDims = Mesh_GetDimSize( self );
	elType = FeMesh_GetElementType( self, element );
	ElementType_ShapeFunctionsGlobalDerivs( elType, self, element, localCoord, nDims, 
						&jd, derivs );
	if( jacDet )
		*jacDet = jd;
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Private Functions
*/

void FeMesh_Destruct( FeMesh* self ) {
	self->feElFamily = NULL;
	KillObject( self->feElType );
}
