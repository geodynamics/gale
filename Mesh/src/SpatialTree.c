/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003, Victorian Partnership for Advanced Computing (VPAC) Ltd, 
** 110 Victoria Street, Melbourne, 3053, Australia.
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
** $Id: SpatialTree.c 3952 2007-01-09 06:24:06Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <StGermain/StGermain.h>
#include "types.h"
#include "Mesh.h"
#include "SpatialTree.h"
#include "StGermain/Base/Foundation/ClassDef.h"


#define GETNODECHILDREN( tree, node )		\
  ((void**)(node))

#define GETNODEORIGIN( tree, node )		\
   ((double*)((stgByte*)(node) +		\
	      (tree)->nChld * sizeof(void*)))

#define SETNODEVERTS( tree, node, nVerts, verts )		\
   ((int*)((stgByte*)(node) + (tree)->nChld * sizeof(void*) +	\
	   (tree)->nDims * sizeof(double)))[0] = (nVerts);	\
   ((int**)((stgByte*)(node) + (tree)->nChld * sizeof(void*) +	\
	    (tree)->nDims * sizeof(double) +			\
	    sizeof(int)))[0] = (verts)

#define GETNODEVERTARRAY( tree, node )				\
   ((int**)((stgByte*)(node) + (tree)->nChld * sizeof(void*) +	\
	    (tree)->nDims * sizeof(double) +			\
	    sizeof(int)))[0]

#define GETNODENUMVERTS( tree, node )				\
   ((int*)((stgByte*)(node) + (tree)->nChld * sizeof(void*) +	\
	   (tree)->nDims * sizeof(double)))[0]


void SpatialTree_DestroyNode( SpatialTree* self, void* node );
void SpatialTree_SplitSet( SpatialTree* self, double* orig,
			   int nVerts, int* verts, 
			   int* subSizes, int** subSets,
			   Bool store );
void SpatialTree_SplitNode( SpatialTree* self, void* node, void** parent,
			    double* min, double* max,
			    int nVerts, int* verts );
void SpatialTree_BuildElements( SpatialTree* self, int nVerts, int* verts, 
				int* nEls, int** els );
void SpatialTree_SearchNode( SpatialTree* self, void* node, 
			     const double* pnt, int* nEls, int** els );


void _SpatialTree_Init( void* _self ) {
   SpatialTree* self = Class_Cast( _self, SpatialTree );

   _NewClass_Init( self );
   self->mesh = NULL;
   self->nDims = 0;
   self->min = NULL;
   self->max = NULL;
   self->root = NULL;
   self->tol = 10;
   self->nNodes = 0;
}

void _SpatialTree_Destruct( void* _self ) {
   SpatialTree* self = Class_Cast( _self, SpatialTree );

   if( self->root )
      SpatialTree_DestroyNode( self, self->root );
   _NewClass_Destruct( self );
}

void _SpatialTree_Copy( void* _self, const void* _op ) {
   /*SpatialTree* self = Class_Cast( _self, SpatialTree );*/
   /*const SpatialTree* op = Class_ConstCast( _op, SpatialTree );*/

   abort();
}

void SpatialTree_SetMesh( void* _self, void* mesh ) {
   SpatialTree* self = Class_Cast( _self, SpatialTree );

   SpatialTree_Clear( self );
   self->mesh = mesh;
}

void SpatialTree_Rebuild( void* _self ) {
   SpatialTree* self = Class_Cast( _self, SpatialTree );
   int nVerts, *verts;
   int ii;

   if( !self->mesh )
      return;

   SpatialTree_Clear( self );
   self->nDims = Mesh_GetDimSize( self->mesh );
   self->nChld = 2;
   for( ii = 1; ii < self->nDims; ii++ )
      self->nChld *= 2;

   self->min = Class_Array( self, double, self->nDims );
   self->max = Class_Array( self, double, self->nDims );

   Mesh_GetDomainCoordRange( self->mesh, self->min, self->max );
   nVerts = Mesh_GetDomainSize( self->mesh, 0 );
   verts = Class_Array( self, int, nVerts );
   for( ii = 0; ii < nVerts; ii++ )
      verts[ii] = ii;
   self->root = Class_Array( self, stgByte, self->nChld * sizeof(void*) + 
			     self->nDims * sizeof(double) );
   SpatialTree_SplitNode( self, self->root, &self->root,
			  self->min, self->max, nVerts, verts );
}

Bool SpatialTree_Search( void* _self, const double* pnt, int* nEls, int** els ) {
   SpatialTree* self = Class_Cast( _self, SpatialTree );
   int ii;

   for( ii = 0; ii < self->nDims; ii++ ) {
      if( pnt[ii] < self->min[ii] || pnt[ii] > self->max[ii] )
	 return False;
   }

   SpatialTree_SearchNode( self, self->root, pnt, nEls, els );
   return True;
}

void SpatialTree_Clear( void* _self ) {
   SpatialTree* self = Class_Cast( _self, SpatialTree );

   if( self->root ) {
      SpatialTree_DestroyNode( self, self->root );
      self->root = NULL;
      self->nNodes = 0;
   }


   Class_Free( self, self->min ); self->min = NULL;
   Class_Free( self, self->max ); self->max = NULL;
}

void SpatialTree_SplitNode( SpatialTree* self, void* node, void** parent,
			    double* min, double* max,
			    int nVerts, int* verts )
{
   int ii;

   for( ii = 0; ii < self->nDims; ii++ )
      GETNODEORIGIN( self, node )[ii] = min[ii] + (max[ii] - min[ii]) * 0.5;

   if( nVerts <= self->tol ) {
      int nEls, *els;

      node = (void*)Class_Rearray( self, node, stgByte,
				   self->nChld * sizeof(void*) + 
				   self->nDims * sizeof(double) + 
				   sizeof(int) + sizeof(int*) );
      *parent = node;
      memset( node, 0, self->nChld * sizeof(void*) );

      SpatialTree_BuildElements( self, nVerts, verts, 
				 &nEls, &els );
      Class_Free( self, verts );
      SETNODEVERTS( self, node, nEls, els );
   }
   else {
      void *newNode, *newNodePtr;
      int* subSizes;
      int** subSets;
      double *subMin, *subMax;
      int jj;

      for( ii = 0; ii < self->nChld; ii++ ) {
	 newNode = Class_Array( self, stgByte, self->nChld * sizeof(void*) + 
				self->nDims * sizeof(double) );
	 GETNODECHILDREN( self, node )[ii] = newNode;
      }

      subSizes = Class_Array( self, int, self->nChld );
      SpatialTree_SplitSet( self, GETNODEORIGIN( self, node ),
			    nVerts, verts, subSizes, NULL, False );

      subSets = Class_Array( self, int*, self->nChld );
      for( ii = 0; ii < self->nChld; ii++ )
	 subSets[ii] = Class_Array( self, int, subSizes[ii] );
      SpatialTree_SplitSet( self, GETNODEORIGIN( self, node ),
			    nVerts, verts, subSizes, subSets, True );

      Class_Free( self, verts );

      subMin = Class_Array( self, double, self->nDims );
      subMax = Class_Array( self, double, self->nDims );

      for( ii = 0; ii < self->nChld; ii++ ) {
	 for( jj = 0; jj < self->nDims; jj++ ) {
	    if( ii & (1 << jj) ) {
	       subMin[jj] = GETNODEORIGIN( self, node )[jj];
	       subMax[jj] = max[jj];
	    }
	    else {
	       subMin[jj] = min[jj];
	       subMax[jj] = GETNODEORIGIN( self, node )[jj];
	    }
	 }

	 newNode = GETNODECHILDREN( self, node )[ii];
	 newNodePtr = GETNODECHILDREN( self, node ) + ii;
	 SpatialTree_SplitNode( self, newNode, newNodePtr, subMin, subMax,
				subSizes[ii], subSets[ii] );
      }

      Class_Free( self, subMin );
      Class_Free( self, subMax );
   }
}

void SpatialTree_SplitSet( SpatialTree* self, double* orig,
			   int nVerts, int* verts, 
			   int* subSizes, int** subSets,
			   Bool store )
{
   double* crd;
   int code;
   int ii, jj;

   memset( subSizes, 0, self->nChld * sizeof(int) );
   for( ii = 0; ii < nVerts; ii++ ) {
      crd = Mesh_GetVertex( self->mesh, verts[ii] );
      code = 0;
      for( jj = 0; jj < self->nDims; jj++ ) {
	 if( crd[jj] > orig[jj] )
	    code |= 1 << jj;
      }

      if( store )
	 subSets[code][subSizes[code]++] = verts[ii];
      else
	 subSizes[code]++;
   }
}

void SpatialTree_BuildElements( SpatialTree* self, int nVerts, int* verts, 
				int* nEls, int** els )
{
   int maxEls, *curEls;
   IArray* inc;
   int nIncEls, *incEls;
   int ii, jj, kk;

   maxEls = 0;
   for( ii = 0; ii < nVerts; ii++ )
      maxEls += Mesh_GetIncidenceSize( self->mesh, 0, verts[ii], self->nDims );

   curEls = Class_Array( self, int, maxEls );

   inc = IArray_New();
   maxEls = 0;
   for( ii = 0; ii < nVerts; ii++ ) {
      Mesh_GetIncidence( self->mesh, 0, verts[ii], self->nDims, inc );
      nIncEls = IArray_GetSize( inc );
      incEls = IArray_GetPtr( inc );
      for( jj = 0; jj < nIncEls; jj++ ) {
	 for( kk = 0; kk < maxEls; kk++ ) {
	    if( curEls[kk] == incEls[jj] )
	       break;
	 }
	 if( kk == maxEls )
	    curEls[maxEls++] = incEls[jj];
      }
   }
   NewClass_Delete( inc );

   *nEls = maxEls;
   *els = Class_Rearray( self, curEls, int, maxEls );
}

void SpatialTree_SearchNode( SpatialTree* self, void* node, 
			     const double* pnt, int* nEls, int** els )
{
   if( GETNODECHILDREN( self, node )[0] == NULL ) {
      *nEls = GETNODENUMVERTS( self, node );
      *els = GETNODEVERTARRAY( self, node );
   }
   else {
      double* orig;
      int code;
      int ii;

      orig = GETNODEORIGIN( self, node );
      code = 0;
      for( ii = 0; ii < self->nDims; ii++ ) {
	 if( pnt[ii] > orig[ii] )
	    code |= 1 << ii;
      }

      SpatialTree_SearchNode( self, GETNODECHILDREN( self, node )[code], 
			      pnt, nEls, els );
   }
}

void SpatialTree_DestroyNode( SpatialTree* self, void* node ) {
   Bool leaf;
   int ii;

   if( !node )
      return;

   leaf = True;
   for( ii = 0; ii < self->nChld; ii++ ) {
      if( GETNODECHILDREN( self, node )[ii] )
	 leaf = False;
      SpatialTree_DestroyNode( self, GETNODECHILDREN( self, node )[ii] );
   }

   if( leaf ) {
      Class_Free( self, GETNODEVERTARRAY( self, node ) );
   }

   Class_Free( self, node );
}
