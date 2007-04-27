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
** $Id:  $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <mpi.h>
#include <StGermain/StGermain.h>
#include "StgFEM/Discretisation/Discretisation.h"
#include "StgFEM/SLE/LinearAlgebra/LinearAlgebra.h"
#include "StgFEM/SLE/SystemSetup/SystemSetup.h"

#include "units.h"
#include "types.h"
#include "shortcuts.h"
#include "MGCoarsener.h"
#include "MGCoarsener_RegCartesian.h"


/* Textual name of this class */
const Type MGCoarsener_RegCartesian_Type = "MGCoarsener_RegCartesian";


/*----------------------------------------------------------------------------------------------------------------------------------
** Constructors
*/

MGCoarsener_RegCartesian* MGCoarsener_RegCartesian_New( void ) {
	return _MGCoarsener_RegCartesian_New( 
		sizeof(MGCoarsener_RegCartesian), 
		MGCoarsener_RegCartesian_Type, 
		_MGCoarsener_RegCartesian_Delete, 
		_MGCoarsener_RegCartesian_Print, 
		_MGCoarsener_RegCartesian_Copy, 
		_MGCoarsener_SetMesh, 
		_MGCoarsener_RegCartesian_Coarsen );
}


MGCoarsener_RegCartesian* _MGCoarsener_RegCartesian_New(
		SizeT					_sizeOfSelf, 
		Type						type,
		Stg_Class_DeleteFunction*	_delete,
		Stg_Class_PrintFunction*		_print, 
		Stg_Class_CopyFunction*		_copy, 
		MGCoarsener_SetMeshFunc*		_setMesh, 
		MGCoarsener_CoarsenFunc*		_coarsen )
{
	MGCoarsener_RegCartesian*	self;
	
	/* Allocate memory. */
	self = (MGCoarsener_RegCartesian*)_MGCoarsener_New(
		_sizeOfSelf,
		type,
		_delete,
		_print, 
		_copy, 
		_setMesh, 
		_coarsen );
	
	/* General info */
	
	/* Virtual info */
	
	/* MGCoarsener_RegCartesian info */
	_MGCoarsener_RegCartesian_Init( self );
	
	return self;
}


void MGCoarsener_RegCartesian_Init( MGCoarsener_RegCartesian* self ) {
	/* General info */
	self->type = MGCoarsener_RegCartesian_Type;
	self->_sizeOfSelf = sizeof(MGCoarsener_RegCartesian);
	self->_deleteSelf = False;
	
	/* Virtual info */
	self->_delete = _MGCoarsener_RegCartesian_Delete;
	self->_print = _MGCoarsener_RegCartesian_Print;
	self->_copy = _MGCoarsener_RegCartesian_Copy;
	_Stg_Class_Init( (Stg_Class*)self );
	_MGCoarsener_Init( (MGCoarsener*)self );
	
	/* MGCoarsener_RegCartesian info */
	_MGCoarsener_RegCartesian_Init( self );
}


void _MGCoarsener_RegCartesian_Init( MGCoarsener_RegCartesian* self ) {
	/* General and Virtual info should already be set */
	
	/* MGCoarsener_RegCartesian info */
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Virtual functions
*/

void _MGCoarsener_RegCartesian_Delete( void* coarsener ) {
	MGCoarsener_RegCartesian*	self = (MGCoarsener_RegCartesian*)coarsener;
	
	/* Delete the class itself */
	
	/* Delete parent */
	_MGCoarsener_Delete( self );
}


void _MGCoarsener_RegCartesian_Print( void* coarsener, Stream* stream ) {
	MGCoarsener_RegCartesian*	self = (MGCoarsener_RegCartesian*)coarsener;
	Stream*	myStream;
	
	/* Set the Journal for printing informations */
	myStream = Journal_Register( InfoStream_Type, "MGCoarsener_RegCartesianStream" );
	
	/* Print parent */
	_MGCoarsener_Print( self, stream );
	
	/* General info */
	Journal_Printf( myStream, "MGCoarsener_RegCartesian (ptr): (%p)\n", self );
	
	/* Virtual info */
	
	/* MGCoarsener_RegCartesian info */
}


void* _MGCoarsener_RegCartesian_Copy( void* coarsener, void* dest, Bool deep, Name nameExt, PtrMap* ptrMap ) {
#if 0
	MGCoarsener_RegCartesian*	self = (MGCoarsener_RegCartesian*)coarsener;
	MGCoarsener_RegCartesian*	newMGCoarsener_RegCartesian;
	PtrMap*	map = ptrMap;
	Bool		ownMap = False;
	
	if( !map ) {
		map = PtrMap_New( 10 );
		ownMap = True;
	}
	
	/* Copy goes here. */
	
	if( ownMap ) {
		Stg_Class_Delete( map );
	}
	
	return (void*)newMGCoarsener_RegCartesian;
#endif
	return NULL;
}


/*----------------------------------------------------------------------------------------------------------------------------------
** Public Functions
*/


/*----------------------------------------------------------------------------------------------------------------------------------
** Private Functions
*/

void _MGCoarsener_RegCartesian_Coarsen( void* coarsener, unsigned level, MGMapping* dstMap ) {
	MGCoarsener_RegCartesian*	self = (MGCoarsener_RegCartesian*)coarsener;
	IJK				newElSizes;
	IJK				prevElSizes;
	IJK				topElSizes;
	unsigned			denom;
	unsigned			nDims;
	Mesh*				mesh;
	Grid*				vertGrid;
	unsigned			*localOrigin, *localRange;

	void _MGCoarsener_RegCartesian_Coarsen3D( MGCoarsener_RegCartesian* self, 
						  IJK newElSizes, IJK prevElSizes, 
						  unsigned denom, 
						  MGMapping* dstMap );
	void _MGCoarsener_RegCartesian_Coarsen2D( MGCoarsener_RegCartesian* self, 
						  IJK newElSizes, IJK prevElSizes, IJK topElSizes, 
						  unsigned denom, 
						  MGMapping* dstMap );
	void _MGCoarsener_RegCartesian_Coarsen1D( MGCoarsener_RegCartesian* self, 
						  IJK newElSizes, IJK prevElSizes, 
						  unsigned denom, 
						  MGMapping* dstMap );

	mesh = self->mesh;
	nDims = Mesh_GetDimSize( mesh );
	vertGrid = *(Grid**)ExtensionManager_Get( mesh->info, mesh, 
						  ExtensionManager_GetHandle( mesh->info, "vertexGrid" ) );
	localOrigin = (unsigned*)ExtensionManager_Get( mesh->info, mesh, 
						       ExtensionManager_GetHandle( mesh->info, "localOrigin" ) );
	localRange = (unsigned*)ExtensionManager_Get( mesh->info, mesh, 
						      ExtensionManager_GetHandle( mesh->info, "localRange" ) );

	/*
	** Calculate the new dimensional element counts as well as a few other things.
	*/
	
	{
		unsigned	el_k, el_j, el_i;
		unsigned	dim_i;
		
		/* Quick hack. */
		denom = 2;
		for( dim_i = 0; dim_i < level; dim_i++ ) {
			denom *= 2;
		}

		for( dim_i = 0; dim_i < nDims; dim_i++ ) {
			div_t	res;

			if( vertGrid->sizes[dim_i] == 1 ) {
				newElSizes[dim_i] = 0;
				prevElSizes[dim_i] = 0;
				continue;
			}
			
			res = div( localRange[dim_i], denom );
			assert( res.rem == 0 );
			
			newElSizes[dim_i] = res.quot;
			prevElSizes[dim_i] = newElSizes[dim_i] * 2;
			topElSizes[dim_i] = localRange[dim_i];
		}

		if( nDims <= 2 ) {
			newElSizes[2] = 0;
			prevElSizes[2] = 0;
			topElSizes[2] = 0;
		}

		if( nDims <= 1 ) {
			newElSizes[1] = 0;
			prevElSizes[1] = 0;
			topElSizes[1] = 0;
		}

		/* Calc the number of coarse elements. */
		dstMap->nEls = 1;
		for( dim_i = 0; dim_i < nDims; dim_i++ ) {
			dstMap->nEls *= newElSizes[dim_i];
		}

		/* ... and the number of coarse nodes. */
		dstMap->nNodes = 1;
		for( dim_i = 0; dim_i < nDims; dim_i++ ) {
			dstMap->nNodes *= (newElSizes[dim_i] + 1);
		}

		/* Map the coarse nodes to top-level nodes. */
		dstMap->nodesTop = Memory_Alloc_Array( unsigned, dstMap->nNodes, "MultiGrid" );
		for( el_k = 0; el_k < (newElSizes[2] + 1); el_k++ ) {
			for( el_j = 0; el_j < (newElSizes[1] + 1); el_j++ ) {
				for( el_i = 0; el_i < (newElSizes[0] + 1); el_i++ ) {
					unsigned	ind;
					unsigned	topInds[3];

					ind = el_k * (newElSizes[1] + 1) * (newElSizes[0] + 1) + el_j * (newElSizes[0] + 1) + el_i;
					topInds[0] = el_i * denom + localOrigin[0];
					if( nDims >= 2 )
						topInds[1] = el_j * denom + localOrigin[1];
					if( nDims >= 3 )
						topInds[2] = el_k * denom + localOrigin[2];
					dstMap->nodesTop[ind] = RegularMeshUtils_Node_3DTo1D( mesh, topInds );
					insist( Mesh_GlobalToDomain( mesh, MT_VERTEX, dstMap->nodesTop[ind], 
								     dstMap->nodesTop + ind ) );
				}
			}
		}
	}


	/*
	** Based on the number of active dimensions, call the appropriate coarsener.
	*/

	if( nDims == 3 ) {
		_MGCoarsener_RegCartesian_Coarsen3D( self, 
						     newElSizes, prevElSizes, 
						     denom, 
						     dstMap );
	}
	else if( nDims == 2 ) {
		_MGCoarsener_RegCartesian_Coarsen2D( self, 
						     newElSizes, prevElSizes, topElSizes, 
						     denom, 
						     dstMap );
	}
	else {
		_MGCoarsener_RegCartesian_Coarsen1D( self, 
						     newElSizes, prevElSizes, 
						     denom, 
						     dstMap );
	}
}


void _MGCoarsener_RegCartesian_Coarsen3D( MGCoarsener_RegCartesian* self, 
					  IJK newElSizes, IJK prevElSizes, 
					  unsigned denom, 
					  MGMapping* dstMap )
{
	unsigned	el_i, el_j, el_k;
	unsigned	node_i, node_j, node_k;


	/*
	** Construct local incident nodes.
	*/

	dstMap->nIncNodesLocal = Memory_Alloc_Array( unsigned, dstMap->nEls, "MultiGrid" );
	dstMap->incNodesLocal = Memory_Alloc_Array( unsigned*, dstMap->nEls, "MultiGrid" );
	dstMap->nIncNodesPrev = Memory_Alloc_Array( unsigned, dstMap->nEls, "MultiGrid" );
	dstMap->incNodesPrev = Memory_Alloc_Array( unsigned*, dstMap->nEls, "MultiGrid" );
	dstMap->incElTop = Memory_Alloc_Array( unsigned, dstMap->nEls, "MultiGrid" );
	for( el_k = 0; el_k < newElSizes[2]; el_k++ ) {
		for( el_j = 0; el_j < newElSizes[1]; el_j++ ) {
			for( el_i = 0; el_i < newElSizes[0]; el_i++ ) {
				unsigned	elInd;
				unsigned	nodeInd;

				elInd = el_k * newElSizes[1] * newElSizes[0] + el_j * newElSizes[0] + el_i;
				nodeInd = el_k * (newElSizes[1] + 1) * (newElSizes[0] + 1) + el_j * (newElSizes[0] + 1) + 
					el_i;

				/* Calculate an incident top-level element. */
				dstMap->incElTop[elInd] = RegularMeshUtils_Element_Local3DTo1D( decomp, 
												el_i * denom, 
												el_j * denom, 
												el_k * denom );

				/* Calculate local incidence. */
				dstMap->nIncNodesLocal[elInd] = 8;
				dstMap->incNodesLocal[elInd] = Memory_Alloc_Array( unsigned, 8, "MultiGrid" );

				dstMap->incNodesLocal[elInd][0] = nodeInd;
				dstMap->incNodesLocal[elInd][1] = nodeInd + 1;
				dstMap->incNodesLocal[elInd][2] = nodeInd + (newElSizes[0] + 1) + 1;
				dstMap->incNodesLocal[elInd][3] = nodeInd + (newElSizes[0] + 1);
				dstMap->incNodesLocal[elInd][4] = (newElSizes[1] + 1) * (newElSizes[0] + 1) + nodeInd;
				dstMap->incNodesLocal[elInd][5] = (newElSizes[1] + 1) * (newElSizes[0] + 1) + nodeInd + 1;
				dstMap->incNodesLocal[elInd][6] = (newElSizes[1] + 1) * (newElSizes[0] + 1) + 
					nodeInd + (newElSizes[0] + 1) + 1;
				dstMap->incNodesLocal[elInd][7] = (newElSizes[1] + 1) * (newElSizes[0] + 1) + 
					nodeInd + (newElSizes[0] + 1);

				/* Calculate previous incidence. */
				nodeInd = el_k * 2 * (prevElSizes[1] + 1) * (prevElSizes[0] + 1) + 
					el_j * 2 * (prevElSizes[0] + 1) +
					el_i * 2;
				dstMap->nIncNodesPrev[elInd] = 27;
				dstMap->incNodesPrev[elInd] = Memory_Alloc_Array( unsigned, 27, "MultiGrid" );

				dstMap->incNodesPrev[elInd][0] = nodeInd;
				dstMap->incNodesPrev[elInd][1] = nodeInd + 1;
				dstMap->incNodesPrev[elInd][2] = nodeInd + 2;
				dstMap->incNodesPrev[elInd][3] = (prevElSizes[0] + 1) + nodeInd;
				dstMap->incNodesPrev[elInd][4] = (prevElSizes[0] + 1) + nodeInd + 1;
				dstMap->incNodesPrev[elInd][5] = (prevElSizes[0] + 1) + nodeInd + 2;
				dstMap->incNodesPrev[elInd][6] = (prevElSizes[0] + 1) * 2 + nodeInd;
				dstMap->incNodesPrev[elInd][7] = (prevElSizes[0] + 1) * 2 + nodeInd + 1;
				dstMap->incNodesPrev[elInd][8] = (prevElSizes[0] + 1) * 2 + nodeInd + 2;

				dstMap->incNodesPrev[elInd][9] = (prevElSizes[0] + 1) * (prevElSizes[1] + 1) + 
					nodeInd;
				dstMap->incNodesPrev[elInd][10] = (prevElSizes[0] + 1) * (prevElSizes[1] + 1) + 
					nodeInd + 1;
				dstMap->incNodesPrev[elInd][11] = (prevElSizes[0] + 1) * (prevElSizes[1] + 1) + 
					nodeInd + 2;
				dstMap->incNodesPrev[elInd][12] = (prevElSizes[0] + 1) * (prevElSizes[1] + 1) + 
					(prevElSizes[0] + 1) + nodeInd;
				dstMap->incNodesPrev[elInd][13] = (prevElSizes[0] + 1) * (prevElSizes[1] + 1) + 
					(prevElSizes[0] + 1) + nodeInd + 1;
				dstMap->incNodesPrev[elInd][14] = (prevElSizes[0] + 1) * (prevElSizes[1] + 1) + 
					(prevElSizes[0] + 1) + nodeInd + 2;
				dstMap->incNodesPrev[elInd][15] = (prevElSizes[0] + 1) * (prevElSizes[1] + 1) + 
					(prevElSizes[0] + 1) * 2 + nodeInd;
				dstMap->incNodesPrev[elInd][16] = (prevElSizes[0] + 1) * (prevElSizes[1] + 1) + 
					(prevElSizes[0] + 1) * 2 + nodeInd + 1;
				dstMap->incNodesPrev[elInd][17] = (prevElSizes[0] + 1) * (prevElSizes[1] + 1) + 
					(prevElSizes[0] + 1) * 2 + nodeInd + 2;

				dstMap->incNodesPrev[elInd][18] = (prevElSizes[0] + 1) * (prevElSizes[1] + 1) * 2 + 
					nodeInd;
				dstMap->incNodesPrev[elInd][19] = (prevElSizes[0] + 1) * (prevElSizes[1] + 1) * 2 + 
					nodeInd + 1;
				dstMap->incNodesPrev[elInd][20] = (prevElSizes[0] + 1) * (prevElSizes[1] + 1) * 2 + 
					nodeInd + 2;
				dstMap->incNodesPrev[elInd][21] = (prevElSizes[0] + 1) * (prevElSizes[1] + 1) * 2 + 
					(prevElSizes[0] + 1) + nodeInd;
				dstMap->incNodesPrev[elInd][22] = (prevElSizes[0] + 1) * (prevElSizes[1] + 1) * 2 + 
					(prevElSizes[0] + 1) + nodeInd + 1;
				dstMap->incNodesPrev[elInd][23] = (prevElSizes[0] + 1) * (prevElSizes[1] + 1) * 2 + 
					(prevElSizes[0] + 1) + nodeInd + 2;
				dstMap->incNodesPrev[elInd][24] = (prevElSizes[0] + 1) * (prevElSizes[1] + 1) * 2 + 
					(prevElSizes[0] + 1) * 2 + nodeInd;
				dstMap->incNodesPrev[elInd][25] = (prevElSizes[0] + 1) * (prevElSizes[1] + 1) * 2 + 
					(prevElSizes[0] + 1) * 2 + nodeInd + 1;
				dstMap->incNodesPrev[elInd][26] = (prevElSizes[0] + 1) * (prevElSizes[1] + 1) * 2 + 
					(prevElSizes[0] + 1) * 2 + nodeInd + 2;
			}
		}
	}


	/*
	** Build the local elements incident on local nodes.
	*/

	dstMap->nIncElsLocal = Memory_Alloc_Array( unsigned, dstMap->nNodes, "MultiGrid" );
	dstMap->incElsLocal = Memory_Alloc_Array( unsigned*, dstMap->nNodes, "MultiGrid" );
	dstMap->elNodeInds = Memory_Alloc_Array( unsigned*, dstMap->nNodes, "MultiGrid" );
	for( node_k = 0; node_k < newElSizes[2] + 1; node_k++ ) {
		for( node_j = 0; node_j < newElSizes[1] + 1; node_j++ ) {
			for( node_i = 0; node_i < newElSizes[0] + 1; node_i++ ) {
				unsigned	nodeInd;
				unsigned	curInd;

				nodeInd = (newElSizes[0] + 1) * (newElSizes[1] + 1) * node_k + 
					(newElSizes[0] + 1) * node_j + 
					+ node_i;

				dstMap->nIncElsLocal[nodeInd] = 1;
				if( node_i > 0 && node_i < newElSizes[0] ) {
					dstMap->nIncElsLocal[nodeInd] *= 2;
				}
				if( node_j > 0 && node_j < newElSizes[1] ) {
					dstMap->nIncElsLocal[nodeInd] *= 2;
				}
				if( node_k > 0 && node_k < newElSizes[2] ) {
					dstMap->nIncElsLocal[nodeInd] *= 2;
				}
				dstMap->incElsLocal[nodeInd] = Memory_Alloc_Array( unsigned, dstMap->nIncElsLocal[nodeInd], 
										   "MultiGrid" );
				dstMap->elNodeInds[nodeInd] = Memory_Alloc_Array( unsigned, dstMap->nIncElsLocal[nodeInd], 
										  "MultiGrid" );

				curInd = 0;
				if( node_i > 0 && node_j > 0 && node_k > 0 ) {
					dstMap->elNodeInds[nodeInd][curInd] = 6;
					dstMap->incElsLocal[nodeInd][curInd++] = newElSizes[0] * newElSizes[1] * (node_k - 1) + 
						newElSizes[0] * (node_j - 1) +
						(node_i - 1);
				}
				if( node_i < newElSizes[0] && node_j > 0 && node_k > 0 ) {
					dstMap->elNodeInds[nodeInd][curInd] = 7;
					dstMap->incElsLocal[nodeInd][curInd++] = newElSizes[0] * newElSizes[1] * (node_k - 1) + 
						newElSizes[0] * (node_j - 1) +
						(node_i);
				}
				if( node_i > 0 && node_j < newElSizes[1] && node_k > 0 ) {
					dstMap->elNodeInds[nodeInd][curInd] = 5;
					dstMap->incElsLocal[nodeInd][curInd++] = newElSizes[0] * newElSizes[1] * (node_k - 1) + 
						newElSizes[0] * (node_j) +
						(node_i - 1);
				}
				if( node_i < newElSizes[0] && node_j < newElSizes[1] && node_k > 0 ) {
					dstMap->elNodeInds[nodeInd][curInd] = 4;
					dstMap->incElsLocal[nodeInd][curInd++] = newElSizes[0] * newElSizes[1] * (node_k - 1) + 
						newElSizes[0] * (node_j) +
						(node_i);
				}

				if( node_i > 0 && node_j > 0 && node_k < newElSizes[2] ) {
					dstMap->elNodeInds[nodeInd][curInd] = 2;
					dstMap->incElsLocal[nodeInd][curInd++] = newElSizes[0] * newElSizes[1] * (node_k) + 
						newElSizes[0] * (node_j - 1) +
						(node_i - 1);
				}
				if( node_i < newElSizes[0] && node_j > 0 && node_k < newElSizes[2] ) {
					dstMap->elNodeInds[nodeInd][curInd] = 3;
					dstMap->incElsLocal[nodeInd][curInd++] = newElSizes[0] * newElSizes[1] * (node_k) + 
						newElSizes[0] * (node_j - 1) +
						(node_i);
				}
				if( node_i > 0 && node_j < newElSizes[1] && node_k < newElSizes[2] ) {
					dstMap->elNodeInds[nodeInd][curInd] = 1;
					dstMap->incElsLocal[nodeInd][curInd++] = newElSizes[0] * newElSizes[1] * (node_k) + 
						newElSizes[0] * (node_j) +
						(node_i - 1);
				}
				if( node_i < newElSizes[0] && node_j < newElSizes[1] && node_k < newElSizes[2] ) {
					dstMap->elNodeInds[nodeInd][curInd] = 0;
					dstMap->incElsLocal[nodeInd][curInd++] = newElSizes[0] * newElSizes[1] * (node_k) + 
						newElSizes[0] * (node_j) +
						(node_i);
				}
			}
		}
	}
}


void _MGCoarsener_RegCartesian_Coarsen2D( MGCoarsener_RegCartesian* self, 
					  IJK newElSizes, IJK prevElSizes, IJK topElSizes, 
					  unsigned denom, 
					  MGMapping* dstMap )
{
	unsigned	el_i, el_j;
	unsigned	node_i, node_j;


	/*
	** Construct local incident nodes.
	*/

	dstMap->nIncNodesLocal = Memory_Alloc_Array( unsigned, dstMap->nEls, "MultiGrid" );
	dstMap->incNodesLocal = Memory_Alloc_Array( unsigned*, dstMap->nEls, "MultiGrid" );
	dstMap->nIncNodesPrev = Memory_Alloc_Array( unsigned, dstMap->nEls, "MultiGrid" );
	dstMap->incNodesPrev = Memory_Alloc_Array( unsigned*, dstMap->nEls, "MultiGrid" );
	dstMap->incElTop = Memory_Alloc_Array( unsigned, dstMap->nEls, "MultiGrid" );
	for( el_j = 0; el_j < newElSizes[1]; el_j++ ) {
		for( el_i = 0; el_i < newElSizes[0]; el_i++ ) {
			unsigned	elInd;
			unsigned	nodeInd;

			elInd = el_j * newElSizes[0] + el_i;
			nodeInd = el_j * (newElSizes[0] + 1) + el_i;

			/* Calculate an incident top-level element. */
			dstMap->incElTop[elInd] = el_j * denom * topElSizes[0] + el_i * denom;

			/* Calculate local incidence. */
			dstMap->nIncNodesLocal[elInd] = 4;
			dstMap->incNodesLocal[elInd] = Memory_Alloc_Array( unsigned, 4, "MultiGrid" );

			dstMap->incNodesLocal[elInd][0] = nodeInd;
			dstMap->incNodesLocal[elInd][1] = nodeInd + 1;
			dstMap->incNodesLocal[elInd][2] = nodeInd + (newElSizes[0] + 1) + 1;
			dstMap->incNodesLocal[elInd][3] = nodeInd + (newElSizes[0] + 1);

			/* Calculate previous incidence. */
			nodeInd = el_j * 2 * (prevElSizes[0] + 1) + el_i * 2;
			dstMap->nIncNodesPrev[elInd] = 9;
			dstMap->incNodesPrev[elInd] = Memory_Alloc_Array( unsigned, 9, "MultiGrid" );

			dstMap->incNodesPrev[elInd][0] = nodeInd;
			dstMap->incNodesPrev[elInd][1] = nodeInd + 1;
			dstMap->incNodesPrev[elInd][2] = nodeInd + 2;
			dstMap->incNodesPrev[elInd][3] = (prevElSizes[0] + 1) + nodeInd;
			dstMap->incNodesPrev[elInd][4] = (prevElSizes[0] + 1) + nodeInd + 1;
			dstMap->incNodesPrev[elInd][5] = (prevElSizes[0] + 1) + nodeInd + 2;
			dstMap->incNodesPrev[elInd][6] = (prevElSizes[0] + 1) * 2 + nodeInd;
			dstMap->incNodesPrev[elInd][7] = (prevElSizes[0] + 1) * 2 + nodeInd + 1;
			dstMap->incNodesPrev[elInd][8] = (prevElSizes[0] + 1) * 2 + nodeInd + 2;
		}
	}


	/*
	** Build the local elements incident on local nodes.
	*/

	dstMap->nIncElsLocal = Memory_Alloc_Array( unsigned, dstMap->nNodes, "MultiGrid" );
	dstMap->incElsLocal = Memory_Alloc_Array( unsigned*, dstMap->nNodes, "MultiGrid" );
	dstMap->elNodeInds = Memory_Alloc_Array( unsigned*, dstMap->nNodes, "MultiGrid" );
	for( node_j = 0; node_j < newElSizes[1] + 1; node_j++ ) {
		for( node_i = 0; node_i < newElSizes[0] + 1; node_i++ ) {
			unsigned	nodeInd;
			unsigned	curInd;

			nodeInd = (newElSizes[0] + 1) * node_j + node_i;

			dstMap->nIncElsLocal[nodeInd] = 1;
			if( node_i > 0 && node_i < newElSizes[0] ) {
				dstMap->nIncElsLocal[nodeInd] *= 2;
			}
			if( node_j > 0 && node_j < newElSizes[1] ) {
				dstMap->nIncElsLocal[nodeInd] *= 2;
			}
			dstMap->incElsLocal[nodeInd] = Memory_Alloc_Array( unsigned, dstMap->nIncElsLocal[nodeInd], 
									   "MultiGrid" );
			dstMap->elNodeInds[nodeInd] = Memory_Alloc_Array( unsigned, dstMap->nIncElsLocal[nodeInd], 
									  "MultiGrid" );

			curInd = 0;
			if( node_i > 0 && node_j > 0 ) {
				dstMap->elNodeInds[nodeInd][curInd] = 2;
				dstMap->incElsLocal[nodeInd][curInd++] = newElSizes[0] * (node_j - 1) +	(node_i - 1);
			}
			if( node_i < newElSizes[0] && node_j > 0 ) {
				dstMap->elNodeInds[nodeInd][curInd] = 3;
				dstMap->incElsLocal[nodeInd][curInd++] = newElSizes[0] * (node_j - 1) +	(node_i);
			}
			if( node_i > 0 && node_j < newElSizes[1] ) {
				dstMap->elNodeInds[nodeInd][curInd] = 1;
				dstMap->incElsLocal[nodeInd][curInd++] = newElSizes[0] * (node_j) + (node_i - 1);
			}
			if( node_i < newElSizes[0] && node_j < newElSizes[1] ) {
				dstMap->elNodeInds[nodeInd][curInd] = 0;
				dstMap->incElsLocal[nodeInd][curInd++] = newElSizes[0] * (node_j) + (node_i);
			}
		}
	}
}


void _MGCoarsener_RegCartesian_Coarsen1D( MGCoarsener_RegCartesian* self, 
					  IJK newElSizes, IJK prevElSizes, 
					  unsigned denom, 
					  MGMapping* dstMap )
{
	unsigned	el_i;
	unsigned	node_i;


	/*
	** Construct local incident nodes.
	*/

	dstMap->nIncNodesLocal = Memory_Alloc_Array( unsigned, dstMap->nEls, "MultiGrid" );
	dstMap->incNodesLocal = Memory_Alloc_Array( unsigned*, dstMap->nEls, "MultiGrid" );
	dstMap->nIncNodesPrev = Memory_Alloc_Array( unsigned, dstMap->nEls, "MultiGrid" );
	dstMap->incNodesPrev = Memory_Alloc_Array( unsigned*, dstMap->nEls, "MultiGrid" );
	dstMap->incElTop = Memory_Alloc_Array( unsigned, dstMap->nEls, "MultiGrid" );
	for( el_i = 0; el_i < newElSizes[0]; el_i++ ) {
		unsigned	elInd;
		unsigned	nodeInd;

		elInd = el_i;
		nodeInd = el_i;

		/* Calculate an incident top-level element. */
		dstMap->incElTop[elInd] = el_i * denom;

		/* Calculate local incidence. */
		dstMap->nIncNodesLocal[elInd] = 2;
		dstMap->incNodesLocal[elInd] = Memory_Alloc_Array( unsigned, 2, "MultiGrid" );

		dstMap->incNodesLocal[elInd][0] = nodeInd;
		dstMap->incNodesLocal[elInd][1] = nodeInd + 1;

		/* Calculate previous incidence. */
		nodeInd = el_i * 2;
		dstMap->nIncNodesPrev[elInd] = 3;
		dstMap->incNodesPrev[elInd] = Memory_Alloc_Array( unsigned, 3, "MultiGrid" );

		dstMap->incNodesPrev[elInd][0] = nodeInd;
		dstMap->incNodesPrev[elInd][1] = nodeInd + 1;
		dstMap->incNodesPrev[elInd][2] = nodeInd + 2;
	}


	/*
	** Build the local elements incident on local nodes.
	*/

	dstMap->nIncElsLocal = Memory_Alloc_Array( unsigned, dstMap->nNodes, "MultiGrid" );
	dstMap->incElsLocal = Memory_Alloc_Array( unsigned*, dstMap->nNodes, "MultiGrid" );
	dstMap->elNodeInds = Memory_Alloc_Array( unsigned*, dstMap->nNodes, "MultiGrid" );
	for( node_i = 0; node_i < newElSizes[0] + 1; node_i++ ) {
		unsigned	nodeInd;
		unsigned	curInd;

		nodeInd = node_i;

		dstMap->nIncElsLocal[nodeInd] = 1;
		if( node_i > 0 && node_i < newElSizes[0] ) {
			dstMap->nIncElsLocal[nodeInd] *= 2;
		}
		dstMap->incElsLocal[nodeInd] = Memory_Alloc_Array( unsigned, dstMap->nIncElsLocal[nodeInd], 
								   "MultiGrid" );
		dstMap->elNodeInds[nodeInd] = Memory_Alloc_Array( unsigned, dstMap->nIncElsLocal[nodeInd], 
								  "MultiGrid" );

		curInd = 0;
		if( node_i > 0 ) {
			dstMap->elNodeInds[nodeInd][curInd] = 1;
			dstMap->incElsLocal[nodeInd][curInd++] = (node_i - 1);
		}
		if( node_i < newElSizes[0] ) {
			dstMap->elNodeInds[nodeInd][curInd] = 0;
			dstMap->incElsLocal[nodeInd][curInd++] = (node_i);
		}
	}
}
