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
** $Id: SemiRegDeform.c 2192 2004-10-15 02:45:38Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <mpi.h>

#include <StGermain/StGermain.h>
#include <StgDomain/Geometry/Geometry.h>
#include <StgDomain/Shape/Shape.h>
#include <StgDomain/Mesh/Mesh.h>

#include "types.h"
#include "SemiRegDeform.h"


/* Textual name of this class */
const Type SemiRegDeform_Type = "SemiRegDeform";


/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/

SemiRegDeform* SemiRegDeform_DefaultNew( Name name ) {
	return _SemiRegDeform_New( 
		sizeof(SemiRegDeform), 
		SemiRegDeform_Type, 
		_SemiRegDeform_Delete, 
		_SemiRegDeform_Print, 
		NULL, 
		(Stg_Component_DefaultConstructorFunction*)SemiRegDeform_DefaultNew, 
		_SemiRegDeform_Construct, 
		_SemiRegDeform_Build, 
		_SemiRegDeform_Initialise, 
		_SemiRegDeform_Execute, 
		_SemiRegDeform_Destroy, 
		name );
}

SemiRegDeform* SemiRegDeform_New( Name name ) {
	return _SemiRegDeform_New( 
		sizeof(SemiRegDeform), 
		SemiRegDeform_Type, 
		_SemiRegDeform_Delete, 
		_SemiRegDeform_Print, 
		NULL, 
		(Stg_Component_DefaultConstructorFunction*)SemiRegDeform_DefaultNew, 
		_SemiRegDeform_Construct, 
		_SemiRegDeform_Build, 
		_SemiRegDeform_Initialise, 
		_SemiRegDeform_Execute, 
		_SemiRegDeform_Destroy, 
		name );
}

SemiRegDeform* _SemiRegDeform_New( SizeT					_sizeOfSelf, 
				   Type						type,
				   Stg_Class_DeleteFunction*			_delete,
				   Stg_Class_PrintFunction*			_print, 
				   Stg_Class_CopyFunction*			_copy, 
				   Stg_Component_DefaultConstructorFunction*	_defaultConstructor,
				   Stg_Component_ConstructFunction*		_construct,
				   Stg_Component_BuildFunction*			_build,
				   Stg_Component_InitialiseFunction*		_initialise,
				   Stg_Component_ExecuteFunction*		_execute,
				   Stg_Component_DestroyFunction*		_destroy, 
				   Name						name )
{
	SemiRegDeform*	self;
	
	/* Allocate memory. */
	self = (SemiRegDeform*)_Stg_Component_New(
		_sizeOfSelf,
		type,
		_delete,
		_print, 
		_copy, 
		_defaultConstructor, 
		_construct, 
		_build, 
		_initialise, 
		_execute, 
		_destroy, 
		name, 
		NON_GLOBAL );
	
	/* General info */
	
	/* Virtual info */
	
	/* SemiRegDeform info */
	_SemiRegDeform_Init( self );
	
	return self;
}

void SemiRegDeform_Init( SemiRegDeform* self, Name name ) {
	/* General info */
	self->type = SemiRegDeform_Type;
	self->_sizeOfSelf = sizeof(SemiRegDeform);
	self->_deleteSelf = False;
	
	/* Virtual info */
	self->_delete = _SemiRegDeform_Delete;
	self->_print = _SemiRegDeform_Print;
	self->_copy = NULL;
	self->_defaultConstructor = (Stg_Component_DefaultConstructorFunction*)SemiRegDeform_DefaultNew;
	self->_construct = _SemiRegDeform_Construct;
	self->_build = _SemiRegDeform_Build;
	self->_execute = _SemiRegDeform_Execute;
	self->_destroy = _SemiRegDeform_Destroy;

	_Stg_Class_Init( (Stg_Class*)self );
	_Stg_Object_Init( (Stg_Object*)self, name, NON_GLOBAL );
	_Stg_Component_Init( (Stg_Component*)self );
	
	/* SemiRegDeform info */
	_SemiRegDeform_Init( self );
}

void _SemiRegDeform_Init( SemiRegDeform* self ) {
	self->nStrips = 0;
	self->beginInds = NULL;
	self->endInds = NULL;
	self->conDims = NULL;
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _SemiRegDeform_Delete( void* srd ) {
	SemiRegDeform*	self = (SemiRegDeform*)srd;
	
	/* Delete the class itself */
	SemiRegDeform_Destruct( self );
	
	/* Delete parent */
	_Stg_Class_Delete( self );
}

void _SemiRegDeform_Print( void* srd, Stream* stream ) {
	SemiRegDeform*	self = (SemiRegDeform*)srd;
	Stream*	myStream;
	
	/* Set the Journal for printing informations */
	myStream = Journal_Register( InfoStream_Type, "SemiRegDeformStream" );
	
	/* Print parent */
	_Stg_Class_Print( self, stream );
	
	/* General info */
	Journal_Printf( myStream, "SemiRegDeform (ptr): (%p)\n", self );
	
	/* Virtual info */
	
	/* SemiRegDeform info */
}

void _SemiRegDeform_Construct( void* srd, Stg_ComponentFactory* cf, void* data ) {
}

void _SemiRegDeform_Build( void* srd, void* data ) {
}

void _SemiRegDeform_Initialise( void* srd, void* data ) {
	SemiRegDeform*	self = (SemiRegDeform*)srd;

	assert( self );


	/*
	** Validate the provided strips.  No two strips can occupy any of the same nodes in the same
	** dimension.
	*/

	/* TODO */


	/*
	** Initialise the synchronisation.
	*/

	SemiRegDeform_InitSync( self );
}

void _SemiRegDeform_Execute( void* srd, void* data ) {
}

void _SemiRegDeform_Destroy( void* srd, void* data ) {
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Public Functions
*/

void SemiRegDeform_SetMesh( void* srd, Mesh* mesh ) {
	SemiRegDeform*	self = (SemiRegDeform*)srd;

	assert( !self->isInitialised );

	SemiRegDeform_Destruct( self );
	self->mesh = mesh;
}

void SemiRegDeform_AddStrip( void* srd, unsigned begin, unsigned end ) {
	SemiRegDeform*	self = (SemiRegDeform*)srd;
	Grid*		vertGrid;
	unsigned short	conDim;
	IJK		inds[2];
	Bool		store;
	unsigned	s_i;

	assert( self->mesh );
	assert( !self->isInitialised );

	/* Get the vertex grid. */
	vertGrid = *(Grid**)ExtensionManager_Get( self->mesh->info, self->mesh, 
						  ExtensionManager_GetHandle( self->mesh->info, "vertexGrid" ) );

	/*
	** Ensure the specified strip has not already been added.
	*/

#ifndef NDEBUG
	for( s_i = 0; s_i < self->nStrips; s_i++ ) {
		if( self->beginInds[s_i] == begin && self->endInds[s_i] == end ) {
			assert( 0 );
		}
	}
#endif


	/*
	** Ensure the specified strip is one dimensionally valid.
	*/

	{
		Bool		found;
		unsigned	d_i;

		Grid_Lift( vertGrid, begin, inds[0] );
		Grid_Lift( vertGrid, end, inds[1] );

		/* Find the one dimension that is not in-line. */
		found = False;
		for( d_i = 0; d_i < Mesh_GetDimSize( self->mesh ); d_i++ ) {
			if( inds[0][d_i] != inds[1][d_i] ) {
				/* Check if we have found multiple connected dimensions. */
				assert( found == False );

				found = True;
				conDim = d_i;
			}
		}
	}


	/*
	** If the strip has no points stored locally then don't store on this processor.
	*/

	{
		unsigned	len = inds[1][conDim] - inds[0][conDim] + 1;
		IJK		cur;
		unsigned	n_i;

		store = False;
		memcpy( cur, inds[0], sizeof(IJK) );
		for( n_i = 0; n_i < len; n_i++ ) {
			unsigned	gInd, dInd;

			gInd = Grid_Project( vertGrid, cur );
			if( Mesh_GlobalToDomain( self->mesh, MT_VERTEX, gInd, &dInd ) && 
			    dInd < Mesh_GetLocalSize( self->mesh, MT_VERTEX ) )
			{
				store = True;
				break;
			}
			cur[conDim]++;
		}
	}

	if( !store )
		return;

	/*
	** Store.
	*/

	self->beginInds = Memory_Realloc_Array( self->beginInds, unsigned, self->nStrips + 1 );
	self->endInds = Memory_Realloc_Array( self->endInds, unsigned, self->nStrips + 1 );
	self->conDims = Memory_Realloc_Array( self->conDims, unsigned, self->nStrips + 1 );
	self->beginInds[self->nStrips] = begin;
	self->endInds[self->nStrips] = end;
	self->conDims[self->nStrips] = conDim;
	self->nStrips++;
}

#define GET_VAL( ind )							\
	(((ind) < Mesh_GetLocalSize( self->mesh, MT_VERTEX )) ? self->mesh->verts[ind] : \
	 self->remVerts[ind - Mesh_GetLocalSize( self->mesh, MT_VERTEX )])

void SemiRegDeform_Deform( void* srd ) {
	SemiRegDeform*	self = (SemiRegDeform*)srd;
	Grid*		vertGrid;

	assert( self );

	/* Get the vertex grid. */
	vertGrid = *(Grid**)ExtensionManager_Get( self->mesh->info, self->mesh, 
						  ExtensionManager_GetHandle( self->mesh->info, "vertexGrid" ) );

	/*
	** Actually deform the specified strips.
	*/

	/* Import remote values. */
	Decomp_Sync_SyncArray( self->sync, self->syncArray );

	/* Interpolate each strip. */
	{
		unsigned	nDims;
		unsigned*	begin;
		unsigned*	end;
		unsigned	strip_i;

		/* Get dimensionality. */
		nDims = Mesh_GetDimSize( self->mesh );

		/* Allocate for the dimensions. */
		begin = Memory_Alloc_Array( unsigned, nDims, "SemiRegDeform" );
		end = Memory_Alloc_Array( unsigned, nDims, "SemiRegDeform" );

		for( strip_i = 0; strip_i < self->nStrips; strip_i++ ) {
			unsigned	len;
			unsigned	conDim;
			double		first, step;
			unsigned	dInd;
			unsigned	node_i;

			/* Extract the basics. */
			Grid_Lift( vertGrid, self->beginInds[strip_i], begin );
			Grid_Lift( vertGrid, self->endInds[strip_i], end );
			conDim = self->conDims[strip_i];
			len = end[conDim] - begin[conDim] + 1;
			assert( len > 1 );

			insist( Sync_GlobalToDomain( self->sync, self->beginInds[strip_i], &dInd ), == True );
			first = GET_VAL( dInd )[conDim];
			insist( Sync_GlobalToDomain( self->sync, self->endInds[strip_i], &dInd ), == True );
			step = GET_VAL( dInd )[conDim];
			step = (step - first) / (len - 1);

			/* Loop and interpolate. */
			for( node_i = 1; node_i < len - 1; node_i++ ) {
				unsigned	ind;

				begin[conDim]++;
				ind = Grid_Project( vertGrid, begin );
				if( Sync_GlobalToDomain( self->sync, ind, &dInd ) ) {
					GET_VAL( dInd )[conDim] = first + (double)node_i * step;
				}
			}
		}
		FreeArray( begin );
		FreeArray( end );
	}
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Private Functions
*/


void SemiRegDeform_InitSync( SemiRegDeform* self ) {
	unsigned	nRequired;
	unsigned*	required;
	unsigned	nDims;
	unsigned	strip_i;

	assert( self );

	/*
	** Setup the synchronisation component.
	*/

	/* Build required indices. */
	nRequired = 0;
	required = Memory_Alloc_Array( unsigned, self->nStrips * 2, "SemiRegDeform" );
	for( strip_i = 0; strip_i < self->nStrips; strip_i++ ) {
		unsigned	dInd;

		if( !Mesh_GlobalToDomain( self->mesh, MT_VERTEX, self->beginInds[strip_i], &dInd ) || 
		    dInd >= Mesh_GetLocalSize( self->mesh, MT_VERTEX ) )
		{
			required[nRequired++] = self->beginInds[strip_i];
		}
		if( !Mesh_GlobalToDomain( self->mesh, MT_VERTEX, self->endInds[strip_i], &dInd ) || 
		    dInd >= Mesh_GetLocalSize( self->mesh, MT_VERTEX ) )
		{
			required[nRequired++] = self->endInds[strip_i];
		}
	}
	required = Memory_Realloc_Array( required, unsigned, nRequired );

	self->sync = Decomp_Sync_New( "" );
	Decomp_Sync_SetDecomp( self->sync, self->mesh->topo->domains[MT_VERTEX]->decomp );
	Decomp_Sync_SetRequired( self->sync, nRequired, required );

	/* Free arrays. */
	FreeArray( required );

	/* Allocate for sources. */
	nDims = Mesh_GetDimSize( self->mesh );
	self->remVerts = AllocNamedArray2D( double, Decomp_Sync_GetRemoteSize( self->sync ), nDims, 
					    "SemiRegDeform::remVerts" );

	/* Initialise array. */
	self->syncArray = Decomp_Sync_Array_New();
	Decomp_Sync_Array_SetSync( self->syncArray, self->sync );
	Decomp_Sync_Array_SetMemory( self->syncArray, 
				     self->mesh->verts[0], self->remVerts ? self->remVerts[0] : NULL, 
				     sizeof(double) * nDims, sizeof(double) * nDims, 
				     sizeof(double) * nDims );
}

void SemiRegDeform_Destruct( SemiRegDeform* self ) {
	assert( self );

	self->mesh = NULL;
	KillObject( self->syncArray );
	KillObject( self->sync );

	KillArray( self->beginInds );
	KillArray( self->endInds );
	KillArray( self->conDims );
}
