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
** $Id: FeMesh_Algorithms.c 3584 2006-05-16 11:11:07Z PatrickSunter $
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
const Type FeMesh_Algorithms_Type = "FeMesh_Algorithms";

/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/

FeMesh_Algorithms* FeMesh_Algorithms_New( Name name ) {
	return _FeMesh_Algorithms_New( sizeof(FeMesh_Algorithms), 
				     FeMesh_Algorithms_Type, 
				     _FeMesh_Algorithms_Delete, 
				     _FeMesh_Algorithms_Print, 
				     NULL, 
				     (void* (*)(Name))_FeMesh_Algorithms_New, 
				     _FeMesh_Algorithms_AssignFromXML, 
				     _FeMesh_Algorithms_Build, 
				     _FeMesh_Algorithms_Initialise, 
				     _FeMesh_Algorithms_Execute, 
				     _FeMesh_Algorithms_Destroy, 
				     name, 
				     NON_GLOBAL, 
				     _Mesh_Algorithms_SetMesh, 
				     _Mesh_Algorithms_Update, 
				     _Mesh_Algorithms_NearestVertex, 
				     _FeMesh_Algorithms_Search, 
				     _FeMesh_Algorithms_SearchElements, 
				     _Mesh_Algorithms_GetMinimumSeparation, 
				     _Mesh_Algorithms_GetLocalCoordRange, 
				     _Mesh_Algorithms_GetDomainCoordRange, 
				     _Mesh_Algorithms_GetGlobalCoordRange );
}

FeMesh_Algorithms* _FeMesh_Algorithms_New( FEMESH_ALGORITHMS_DEFARGS ) {
	FeMesh_Algorithms* self;
	
	/* Allocate memory */
	assert( sizeOfSelf >= sizeof(FeMesh_Algorithms) );
	self = (FeMesh_Algorithms*)_Mesh_Algorithms_New( MESH_ALGORITHMS_PASSARGS );

	/* Virtual info */

	/* FeMesh_Algorithms info */
	_FeMesh_Algorithms_Init( self );

	return self;
}

void _FeMesh_Algorithms_Init( FeMesh_Algorithms* self ) {
	_Mesh_Algorithms_Init( self );
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _FeMesh_Algorithms_Delete( void* algorithms ) {
	FeMesh_Algorithms*	self = (FeMesh_Algorithms*)algorithms;

	/* Delete the parent. */
	_Mesh_Algorithms_Delete( self );
}

void _FeMesh_Algorithms_Print( void* algorithms, Stream* stream ) {
	FeMesh_Algorithms*	self = (FeMesh_Algorithms*)algorithms;
	
	/* Set the Journal for printing informations */
	Stream* algorithmsStream;
	algorithmsStream = Journal_Register( InfoStream_Type, "FeMesh_AlgorithmsStream" );

	/* Print parent */
	Journal_Printf( stream, "FeMesh_Algorithms (ptr): (%p)\n", self );
	_Stg_Component_Print( self, stream );
}

void _FeMesh_Algorithms_AssignFromXML( void* algorithms, Stg_ComponentFactory* cf, void* data ) {
}

void _FeMesh_Algorithms_Build( void* algorithms, void* data ) {
}

void _FeMesh_Algorithms_Initialise( void* algorithms, void* data ) {
}

void _FeMesh_Algorithms_Execute( void* algorithms, void* data ) {
}

void _FeMesh_Algorithms_Destroy( void* algorithms, void* data ) {
}

Bool _FeMesh_Algorithms_Search( void* algorithms, double* point, 
			      MeshTopology_Dim* dim, unsigned* ind )
{
	FeMesh_Algorithms*	self = (FeMesh_Algorithms*)algorithms;

	return FeMesh_Algorithms_SearchWithTree( self, point, dim, ind );
}

Bool _FeMesh_Algorithms_SearchElements( void* algorithms, double* point, 
				      unsigned* elInd )
{
	FeMesh_Algorithms*	self = (FeMesh_Algorithms*)algorithms;
	unsigned		dim;

	assert( self );

	return FeMesh_Algorithms_SearchWithTree( self, point, &dim, elInd );
}



/*--------------------------------------------------------------------------------------------------------------------------
** Public Functions
*/

Bool FeMesh_Algorithms_SearchWithTree( void* _self, double* pnt, unsigned* dim, unsigned* el ) {
   FeMesh_Algorithms* self = (FeMesh_Algorithms*)_self;
   int nEls, *els;
   int curDim, curRank, curEl;
   int nLocals, owner;
   FeMesh_ElementType* elType;
   int ii;

   *dim = Mesh_GetDimSize( self->mesh );
   MPI_Comm_size( MPI_COMM_WORLD, &curRank );
   nLocals = Mesh_GetLocalSize( self->mesh, *dim );
   if( !SpatialTree_Search( self->tree, pnt, &nEls, &els ) )
      return False;

   *el = nLocals;
   elType = Mesh_GetElementType( self->mesh, 0 );
   for( ii = 0; ii < nEls; ii++ ) {
      if( FeMesh_ElementType_ElementHasPoint( elType, els[ii], pnt, &curDim, &curEl ) ) {
	 if( curEl >= nLocals ) {
	    owner = Mesh_GetOwner( self->mesh, curDim, curEl - nLocals );
	    owner = Comm_RankLocalToGlobal( self->mesh->topo->comm, owner );
	    if( owner <= curRank ) {
	       curRank = owner;
	       *el = curEl;
	    }
	 }
	 else if( self->rank <= curRank && curEl < *el ) {
	    curRank = self->rank;
	    *el = curEl;
	 }
      }
   }

   return True;
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Private Functions
*/
