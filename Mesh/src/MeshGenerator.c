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
** $Id: MeshGenerator.c 3584 2006-05-16 11:11:07Z PatrickSunter $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <mpi.h>
#include <StGermain/StGermain.h>

#include <StgDomain/Geometry/Geometry.h>
#include <StgDomain/Shape/Shape.h>

#include "types.h"
#include "shortcuts.h"
#include "Sync.h"
#include "MeshTopology.h"
#include "MeshClass.h"
#include "MeshGenerator.h"


/* Textual name of this class */
const Type MeshGenerator_Type = "MeshGenerator";

/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/

MeshGenerator* _MeshGenerator_New( MESHGENERATOR_DEFARGS ) {
	MeshGenerator*	self;

	/* Allocate memory */
	assert( sizeOfSelf >= sizeof(MeshGenerator) );
	self = (MeshGenerator*)_Stg_Component_New( STG_COMPONENT_PASSARGS );

	/* Virtual info */
	self->setDimSizeFunc = setDimSizeFunc;
	self->generateFunc = generateFunc;

	/* MeshGenerator info */
	_MeshGenerator_Init( self );

	return self;
}

void _MeshGenerator_Init( MeshGenerator* self ) {
	self->mpiComm = MPI_COMM_WORLD;
	self->nMeshes = 0;
	self->meshes = NULL;
	self->nDims = 0;
	self->enabledDims = NULL;
	self->enabledInc = NULL;
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _MeshGenerator_Delete( void* meshGenerator ) {
	MeshGenerator*	self = (MeshGenerator*)meshGenerator;

	MeshGenerator_Destruct( self );

	/* Delete the parent. */
	_Stg_Component_Delete( self );
}

void _MeshGenerator_Print( void* meshGenerator, Stream* stream ) {
	MeshGenerator*	self = (MeshGenerator*)meshGenerator;
	
	/* Set the Journal for printing informations */
	Stream* meshGeneratorStream;
	meshGeneratorStream = Journal_Register( InfoStream_Type, "MeshGeneratorStream" );

	/* Print parent */
	Journal_Printf( stream, "MeshGenerator (ptr): (%p)\n", self );
	_Stg_Component_Print( self, stream );
}

void _MeshGenerator_AssignFromXML( void* meshGenerator, Stg_ComponentFactory* cf, void* data ) {
	MeshGenerator*		self = (MeshGenerator*)meshGenerator;
	Dictionary*		dict;
	unsigned		nDims;
	Dictionary_Entry_Value*	meshList;
	Dictionary_Entry_Value	*enabledDimsList, *enabledIncList;
	Mesh*			mesh;

	assert( self );
	assert( cf );

	/* Rip out the components structure as a dictionary. */
	dict = Dictionary_Entry_Value_AsDictionary( Dictionary_Get( cf->componentDict, self->name ) );

	/* Set the communicator to a default. */
	MeshGenerator_SetMPIComm( self, MPI_COMM_WORLD );

	self->context = Stg_ComponentFactory_ConstructByKey( cf, self->name, "Context", AbstractContext, False, data );
	if( !self->context )
		self->context = Stg_ComponentFactory_ConstructByName( cf, "context", AbstractContext, True, data );

	/* Read the individual mesh if specified. */
	mesh = Stg_ComponentFactory_ConstructByKey( cf, self->name, "mesh", Mesh, False, data );
	if( mesh )
		MeshGenerator_AddMesh( self, mesh );

	/* Read the mesh list, if it's there. */
	meshList = Dictionary_Get( dict, "meshes" );
	if( meshList ) {
		unsigned	nMeshes;
		char*		name;
		unsigned	m_i;

		nMeshes = Dictionary_Entry_Value_GetCount( meshList );
		for( m_i = 0; m_i < nMeshes; m_i++ ) {
			Mesh*	mesh;

			name = Dictionary_Entry_Value_AsString( Dictionary_Entry_Value_GetElement( meshList, m_i ) );
			mesh = Stg_ComponentFactory_ConstructByName( cf, name, Mesh, True, data );
			MeshGenerator_AddMesh( self, mesh );
		}
	}

	/* Read dimensions and state. */
	nDims = Stg_ComponentFactory_GetUnsignedInt( cf, self->name, "dims", 2 );
	MeshGenerator_SetDimSize( self, nDims );
	enabledDimsList = Dictionary_Get( dict, "enabledDims" );
	enabledIncList = Dictionary_Get( dict, "enabledIncidence" );

    /* Clear dims/incidence flags */
    unsigned    d_i;
    memset( self->enabledDims, 0, (nDims + 1) * sizeof(Bool) );
    for( d_i = 0; d_i <= nDims; d_i++ )
        memset( self->enabledInc[d_i], 0, (nDims + 1) * sizeof(Bool) );

    if( enabledDimsList ) {
        unsigned    dim;
        unsigned    nEnabledDims;
        nEnabledDims = Dictionary_Entry_Value_GetCount( enabledDimsList );
        for( d_i = 0; d_i < nEnabledDims; d_i++ ) {
            dim = Dictionary_Entry_Value_AsUnsignedInt( Dictionary_Entry_Value_GetElement( enabledDimsList, d_i ) );
            if (dim > nDims + 1)
                Journal_Printf(Mesh_Warning, "Warning - in %s: *** Skipping out of range dimension: %d\n", __func__, dim);
            else
                MeshGenerator_SetDimState( self, dim, True );
        }
    }
    else
    {
        /* Default to all dimensions enabled */
        for( d_i = 0; d_i < nDims + 1; d_i++ ) 
            MeshGenerator_SetDimState( self, d_i, True );
    }

    if( enabledIncList ) {
        unsigned    nEnabledInc;
        unsigned    fromDim, toDim;
        nEnabledInc = Dictionary_Entry_Value_GetCount( enabledIncList );
        assert( nEnabledInc % 2 == 0 );
        for( d_i = 0; d_i < nEnabledInc; d_i += 2 ) {
            fromDim = Dictionary_Entry_Value_AsUnsignedInt( Dictionary_Entry_Value_GetElement( enabledIncList, d_i ) );
            toDim = Dictionary_Entry_Value_AsUnsignedInt( Dictionary_Entry_Value_GetElement( enabledIncList, d_i + 1 ) );
            if (fromDim > nDims || toDim > nDims) 
                Journal_Printf( Mesh_Warning, "Warning - in %s: *** Skipping out of range incidence: %d to %d\n", 
                            __func__ , fromDim, toDim);
            else
                MeshGenerator_SetIncidenceState( self, fromDim, toDim, True );
        }
    }
    else
    {
        /* Default incidence setup 0->1,2,3  1->0,2  2->0,1  3->0 */
        MeshGenerator_SetIncidenceState( self, 0, 0, True );
        for( d_i = 1; d_i <= nDims; d_i ++ ) {
            MeshGenerator_SetIncidenceState( self, 0, d_i, True );
            MeshGenerator_SetIncidenceState( self, d_i, 0, True );
        }
        if (nDims >= 2) {
            MeshGenerator_SetIncidenceState( self, 1, 2, True );
            MeshGenerator_SetIncidenceState( self, 2, 1, True );
        }
    }
}

void _MeshGenerator_Build( void* meshGenerator, void* data ) {
}

void _MeshGenerator_Initialise( void* meshGenerator, void* data ) {
}

void _MeshGenerator_Execute( void* meshGenerator, void* data ) {
}

void _MeshGenerator_Destroy( void* meshGenerator, void* data ) {
}

void MeshGenerator_SetFullIncidence( void* meshGenerator ) {
	MeshGenerator* self = (MeshGenerator*)meshGenerator;
	int d_i, d_j;

	assert( self );

	for( d_i = 0; d_i <= self->nDims; d_i++ ) {
		self->enabledDims[d_i] = True;
		for( d_j = 0; d_j <= self->nDims; d_j++ )
			self->enabledInc[d_i][d_j] = True;
	}
}

void _MeshGenerator_SetDimSize( void* meshGenerator, unsigned nDims ) {
	MeshGenerator*	self = (MeshGenerator*)meshGenerator;
	unsigned	d_i, d_j;

	assert( self && Stg_CheckType( self, MeshGenerator ) );

	self->nDims = nDims;
	self->enabledDims = ReallocArray( self->enabledDims, Bool, nDims + 1 );
	self->enabledInc = ReallocArray2D( self->enabledInc, Bool, nDims + 1, nDims + 1 );
	for( d_i = 0; d_i <= nDims; d_i++ ) {
		if( d_i == 0 || d_i == nDims )
			self->enabledDims[d_i] = True;
		else
			self->enabledDims[d_i] = False;
		for( d_j = 0; d_j <= nDims; d_j++ ) {
			if( (d_i == 0 || d_i == nDims) && (d_j == 0 || d_j == nDims) ) {
				if( d_i == d_j && d_j == nDims )
					self->enabledInc[d_i][d_j] = False;
				else
					self->enabledInc[d_i][d_j] = True;
			}
			else
				self->enabledInc[d_i][d_j] = False;
		}
	}
}

/*--------------------------------------------------------------------------------------------------------------------------
** Public Functions
*/

void MeshGenerator_SetMPIComm( void* meshGenerator, MPI_Comm mpiComm ) {
	MeshGenerator*	self = (MeshGenerator*)meshGenerator;

	/* Sanity check. */
	assert( self );

	/* Should probably kill some stuff when I do this. Oh well. */
	self->mpiComm = mpiComm;
}

void MeshGenerator_AddMesh( void* meshGenerator, void* mesh ) {
	MeshGenerator*	self = (MeshGenerator*)meshGenerator;

	/* Sanity check. */
	assert( self );

	/* If not already allocated, allocate now. */
	if( !self->meshes ) {
		self->nMeshes = 1;
		self->meshes = Memory_Alloc_Array( Mesh*, 1, "MeshGenerator::meshes" );
	}
	else {
		/* Otherwise reallocate. */
		self->nMeshes++;
		self->meshes = Memory_Realloc_Array( self->meshes, Mesh*, self->nMeshes );
	}

	/* Add the new mesh. */
	self->meshes[self->nMeshes - 1] = mesh;
	((Mesh*)mesh)->generator = self;
}

void MeshGenerator_SetDimState( void* meshGenerator, unsigned dim, Bool state ) {
	MeshGenerator*	self = (MeshGenerator*)meshGenerator;

	assert( self && Stg_CheckType( self, MeshGenerator ) );
	assert( dim <= self->nDims );
	assert( self->enabledDims );

	self->enabledDims[dim] = state;
}

void MeshGenerator_ClearIncidenceStates( void* meshGenerator ) {
	MeshGenerator*	self = (MeshGenerator*)meshGenerator;
	unsigned	d_i, d_j;

	assert( self && Stg_CheckType( self, MeshGenerator ) );

	for( d_i = 0; d_i <= self->nDims; d_i++ ) {
		for( d_j = 0; d_j <= self->nDims; d_j++ )
			self->enabledInc[d_i][d_j] = False;
	}
}

void MeshGenerator_SetIncidenceState( void* meshGenerator, unsigned fromDim, unsigned toDim, Bool state ) {
	MeshGenerator*	self = (MeshGenerator*)meshGenerator;

	assert( self && Stg_CheckType( self, MeshGenerator ) );
	assert( fromDim <= self->nDims );
	assert( toDim <= self->nDims );
	assert( self->enabledInc );

	self->enabledInc[fromDim][toDim] = state;
}

/*----------------------------------------------------------------------------------------------------------------------------------
** Private Functions
*/

void MeshGenerator_Destruct( MeshGenerator* self ) {
	assert( self && Stg_CheckType( self, MeshGenerator ) );

	KillArray( self->enabledDims );
	KillArray( self->enabledInc );
	KillArray( self->meshes );
	self->nMeshes = 0;
	self->nDims = 0;
}
