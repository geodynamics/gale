/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003, Victorian Partnership for Advanced Computing (VPAC) Ltd, 
110 Victoria Street, Melbourne, 3053, Australia.
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
** $Id: STree.c 3952 2007-01-09 06:24:06Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#include <StGermain/Base/Foundation/Foundation.h>

#include "types.h"
#include "STree.h"
#include "StGermain/Base/Foundation/ClassDef.h"


void STreeNode_Insert( STreeNode *self, STreeNode *node, STree *tree, STreeNode **par );
void STreeNode_Destroy( STreeNode *self, STree *tree );
int STree_IntCmp( const void* left, const void* right );
void STree_IntDel( void* itm );


int log2i( int x ) {
    int e = 0;
    while((x >> e) != 1 )
	e++;
    return e;
}


void _STree_Init( void* _self ) {
   STree* self = Class_Cast( _self, STree );

   _NewClass_Init( self );
   self->root = NULL;
   self->nNodes = 0;
   self->maxNodes = 0;
   self->alpha = 0.65;
   self->invAlpha = 1.0 / 0.65;
   self->curDepth = 0;
   self->curSize = 0;
   self->flip= 0;
}

void _STree_Destruct( void* _self ) {
   STree* self = Class_Cast( _self, STree );

   if( self->root )
      STreeNode_Destroy( self->root, self );
   _NewClass_Destruct( self );
}

void _STree_Copy( void* _self, const void* _op ) {
   /*STree* self = Class_Cast( _self, STree );*/
   /*const STree* op = Class_ConstCast( _op, STree );*/

   abort();
}

void STree_SetCallbacks( void* _self, STree_CompareCB* cmp, STree_DeleteCB* del ) {
   STree* self = Class_Cast( _self, STree );

   STree_Clear( self );
   self->cmp = cmp;
   self->del = del;
}

void STree_SetIntCallbacks( void* _self ) {
   STree* self = Class_Cast( _self, STree );

   STree_Clear( self );
   self->cmp = STree_IntCmp;
   self->del = STree_IntDel;
}

void STree_SetItemSize( void* _self, int itmSize ) {
   STree* self = Class_Cast( _self, STree );

   STree_Clear( self );
   self->itmSize = itmSize;
}

void STree_SetAlpha( void* _self, float alpha ) {
   STree* self = Class_Cast( _self, STree );

   assert( alpha >= 0.5 && alpha <= 1.0 );
   self->alpha = alpha;
   self->invAlpha = 1.0 / alpha;
}

void STree_Insert( void* _self, const void* itm ) {
   STree* self = Class_Cast( _self, STree );
   STreeNode* node;

   assert( itm );
   node = Class_Alloc( self, STreeNode );
   node->left = NULL;
   node->right = NULL;
   node->data = Class_Array( self, stgByte, self->itmSize );
   memcpy( node->data, itm, self->itmSize );
   if ( self->root ) {
      STreeNode_Insert( self->root, node, self, &self->root );
      self->curDepth = 0;
   }
   else
      self->root = node;
   if ( ++self->nNodes > self->maxNodes )
      self->maxNodes = self->nNodes;
}

void STree_Remove( void* _self, const void* itm ) {
   STree* self = Class_Cast( _self, STree );
   STreeNode *cur = self->root, **pre = &self->root;
   int res;

   assert( itm );
   assert( self->cmp );
   while( (res = self->cmp( itm, cur->data )) ) {
      if ( res < 0 ) {
	 pre = &cur->left;
	 cur = cur->left;
      }
      else {
	 pre = &cur->right;
	 cur = cur->right;
      }
   }
   assert( cur );
   if ( !cur->left )
      *pre = cur->right;
   else if ( !cur->right )
      *pre = cur->left;
   else if ( !cur->left->right ) {
      *pre = cur->left;
      cur->left->right = cur->right;
   }
   else if ( !cur->right->left ) {
      *pre = cur->right;
      cur->right->left = cur->left;
   }
   else if ( self->flip ) {
      STreeNode *last = cur->left, *preLast;
      while ( last->right ) {
	 preLast = last;
	 last = last->right;
      }
      preLast->right = last->left;
      last->left = cur->left;
      last->right = cur->right;
      *pre = last;
      self->flip = 0;
   }
   else {
      STreeNode *last = cur->right, *preLast;
      while ( last->left ) {
	 preLast = last;
	 last = last->left;
      }
      preLast->left = last->right;
      last->right = cur->right;
      last->left = cur->left;
      *pre = last;
      self->flip = 1;
   }
   self->del( cur->data );
   Class_Free( self, cur->data );
   Class_Free( self, cur );
   if ( --self->nNodes <= self->maxNodes / 2 ) {
      self->root = STree_Rebalance( self, self->root, self->nNodes );
      self->maxNodes = self->nNodes;
   }
}

void STree_Clear( void* _self ) {
   STree* self = Class_Cast( _self, STree );

   if ( self->root ) {
      STreeNode_Destroy( self->root, self );
      self->root = NULL;
      self->nNodes = 0;
      self->maxNodes = 0;
      self->curDepth = 0;
      self->curSize = 0;
      self->flip = 0;
   }
}

int STree_GetSize( const void* _self ) {
   const STree* self = Class_ConstCast( _self, STree );

   return self->nNodes;
}

const STreeNode* STree_GetRoot( const void* _self ) {
   const STree* self = Class_ConstCast( _self, STree );

   return self->root;
}

Bool STree_Has( const void* _self, const void* itm ) {
   STree* self = Class_Cast( _self, STree );
   STreeNode *cur = self->root;
   int res;

   assert( self->cmp );
   while ( cur && (res = self->cmp( itm, cur->data )) )
      cur = (res < 0) ? cur->left : cur->right;

   return cur ? True : False;
}

int STree_Size( const STreeNode *node ) {
   if ( node )
      return STree_Size( node->left ) + STree_Size( node->right ) + 1;
   else
      return 0;
}

STreeNode* STree_Rebalance( void* _self, STreeNode *root, int nNodes ) {
   STree* self = Class_Cast( _self, STree );
   STreeNode *pseudo, *tmp;

   pseudo = Class_Alloc( self, STreeNode );
   pseudo->left = NULL;
   pseudo->right = root;
   pseudo->data = NULL;
   STree_Flatten( pseudo );
   STree_Grow( pseudo, nNodes );
   tmp = pseudo->right;
   Class_Free( self, pseudo );

   return tmp;
}

void STree_Flatten( STreeNode *pseudo ) {
   STreeNode *tail, *rem;

   assert( pseudo );
   tail = pseudo;
   rem = tail->right;
   while ( rem ) {
      if ( rem->left ) {
	 STreeNode *tmp = rem->left;
	 rem->left = tmp->right;
	 tmp->right = rem;
	 rem = tmp;
	 tail->right = tmp;
      }
      else {
	 tail = rem;
	 rem = rem->right;
      }
   }
}

void STree_Grow( STreeNode *pseudo, int nNodes ) {
   int nLeaves, nSpineNodes;

   nLeaves = nNodes + 1 - ( 1 << log2i( nNodes + 1 ) );
   STree_Compression( pseudo, nLeaves );
   nSpineNodes = nNodes - nLeaves;
   while ( nSpineNodes > 1 ) {
      nSpineNodes /= 2;
      STree_Compression( pseudo, nSpineNodes );
   }
}

void STree_Compression( STreeNode *pseudo, int nSpineNodes ) {
   STreeNode *scan = pseudo;
   int n_i;

   for ( n_i = 0; n_i < nSpineNodes; n_i++ ) {
      STreeNode *child = scan->right;
      scan->right = child->right;
      scan = scan->right;
      child->right = scan->left;
      scan->left = child;
   }
}

void STreeNode_Insert( STreeNode *self, STreeNode *node, STree *tree, STreeNode **par ) {
   STreeNode **child;
   int res;

   assert( tree->cmp );
   res = tree->cmp( node->data, self->data );
   tree->curDepth++;
   child = (res < 0 || !res) ? &self->left : &self->right;
   if ( *child )
      STreeNode_Insert( *child, node, tree, child );
   else {
      int height;
      *child = node;
      height = ( int )( log2(( float )tree->nNodes + 1 ) /
			log2( tree->invAlpha ) );
      if ( tree->curDepth > height )
	 tree->curSize = 1;
   }
   if ( tree->curSize ) {
      STreeNode **bro = (res < 0) ? &self->right : &self->left;
      int broSize = STree_Size( *bro );
      int nodeSize = tree->curSize + broSize + 1;
      float weight = tree->alpha * ( float )nodeSize;

      if (( float )tree->curSize > weight ||
	  ( float )broSize > weight )
      {
	 *par = STree_Rebalance( tree, self, nodeSize );
	 tree->curSize = 0;
      }
      else
	 tree->curSize = nodeSize;
   }
}

void STreeNode_Destroy( STreeNode *self, STree *tree ) {
   if ( self->left )
      STreeNode_Destroy( self->left, tree );
   if ( self->right )
      STreeNode_Destroy( self->right, tree );
   tree->del( self->data );
   Class_Free( tree, self->data );
   Class_Free( tree, self );
}

int STree_IntCmp( const void* left, const void* right ) {
   return (*((int*)left) < *((int*)right)) ? -1 : 
      (*((int*)left) > *((int*)right)) ? 1 : 0;
}

void STree_IntDel( void* itm ) {
}
