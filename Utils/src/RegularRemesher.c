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
** $Id: RegularRemesher.c 3952 2007-01-09 06:24:06Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifdef HAVE_PETSC

#include <stdlib.h>
#include <string.h>
#include <petsc.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscksp.h>
#include <StGermain/StGermain.h>
#include <StgDomain/Geometry/Geometry.h>
#include <StgDomain/Shape/Shape.h>
#include <StgDomain/Mesh/Mesh.h>

#include "types.h"
#include "NewRemesher.h"
#include "RegularRemesher.h"
#include "StGermain/Base/Foundation/ClassDef.h"


void _RegularRemesher_Init( void* _self ) {
   RegularRemesher* self = Class_Cast( _self, RegularRemesher );

   _NewRemesher_Init( self );
   self->remeshDims = ISet_New();
   self->staticWalls = Class_Array( self, ISet*, 2 );
   self->staticWalls[0] = ISet_New();
   self->staticWalls[1] = ISet_New();
   ISet_SetMaxSize( self->remeshDims, 10 );
   ISet_SetMaxSize( self->staticWalls[0], 10 );
   ISet_SetMaxSize( self->staticWalls[1], 10 );
   self->syncs = NULL;
   self->crds = NULL;
   self->nWallVerts = NULL;
   self->wallVerts = NULL;
   self->wallCrds = NULL;
}

void _RegularRemesher_Copy( void* _self, const void* _op ) {
   RegularRemesher* self = Class_Cast( _self, RegularRemesher );
   const RegularRemesher* op = Class_Cast( _op, RegularRemesher );

   _NewRemesher_Copy( self, op );
   /* TODO */
}

void _RegularRemesher_Destruct( void* _self ) {
   RegularRemesher* self = Class_Cast( _self, RegularRemesher );
   int nDims;
   int d_i;

   if( self->syncs ) {
      nDims = self->nDims;
      for( d_i = 0; d_i < nDims; d_i++ ) {
   	 NewClass_RemoveRef( self->syncs[d_i] );
	 Class_Free( self, self->crds[d_i] );
	 Class_Free( self, self->nWallVerts[d_i] );
	 Class_Free( self, self->wallVerts[d_i][0] );
	 Class_Free( self, self->wallVerts[d_i][1] );
	 Class_Free( self, self->wallVerts[d_i] );
	 Class_Free( self, self->wallCrds[d_i][0] );
	 Class_Free( self, self->wallCrds[d_i][1] );
	 Class_Free( self, self->wallCrds[d_i] );
      }
      Class_Free( self, self->syncs );
      Class_Free( self, self->crds );
      Class_Free( self, self->nWallVerts );
      Class_Free( self, self->wallVerts );
      Class_Free( self, self->wallCrds );
   }
   NewClass_Delete( self->remeshDims );
   NewClass_Delete( self->staticWalls[0] );
   NewClass_Delete( self->staticWalls[1] );
   Class_Free( self, self->staticWalls );
   _NewClass_Destruct( self );
}

void _RegularRemesher_Print( const void* _self, Stream* stream ) {
   RegularRemesher* self = Class_Cast( _self, RegularRemesher );

   _NewRemesher_Print( self, stream );
}

void _RegularRemesher_Remesh( void* _self ) {
   RegularRemesher* self = Class_Cast( _self, RegularRemesher );
   Mesh* mesh;
   Grid *vGrid;
   Sync* meshSync;
   int nDims, nVerts;
   int center, ind, *inds;
   double leftCrd, rightCrd;
   CartesianGenerator* gen;
   int d_i, v_i, w_i;

   assert( self->mesh );

   mesh = self->mesh;
   meshSync = Mesh_GetSync( mesh, 0 );
   nDims = Mesh_GetDimSize( mesh );
   nVerts = Mesh_GetLocalSize( mesh, 0 );
   vGrid = *Mesh_GetExtension( mesh, Grid**, "vertexGrid" );
   inds = Class_Array( self, int, nDims );

   for( d_i = 0; d_i < nDims; d_i++ ) {
      for( w_i = 0; w_i < 2; w_i++ ) {
	 if( !ISet_Has( self->staticWalls[w_i], d_i ) )
	    continue;
	 for( v_i = 0; v_i < self->nWallVerts[d_i][w_i]; v_i++ ) {
	    mesh->verts[self->wallVerts[d_i][w_i][v_i]][d_i] = 
	       self->wallCrds[d_i][w_i][v_i];
	 }
      }
   }

   /* If we have a contact depth set we'll need to manipulate the boundaries
      a little. */
   gen = (CartesianGenerator*)self->mesh->generator;
   if( strcmp( gen->type, CartesianGenerator_Type ) )
      gen = NULL;
   if( gen && (gen->contactDepth[0][0] || gen->contactDepth[0][1] ||
               gen->contactDepth[1][0] || gen->contactDepth[1][1] ||
               gen->contactDepth[2][0] || gen->contactDepth[2][1]) )
   {
      /*int curInd;*/
      int /*ii,*/ d_j;

#if 0
      /* Reset static depths. */
      curInd = 0;
      for( ii = 0; ii < nVerts; ii++ ) {
         Grid_Lift( vGrid, ii, inds );
         if( inds[1] != self->contactDepth ) continue;
         Mesh_GetVertex( mesh, ii )[1] = self->contactVerts[curInd++];
      }
#endif

      /* Also handle contact element boundaries. */
      for( d_i = 0; d_i < nDims; d_i++ ) {
         for( v_i = 0; v_i < nVerts; v_i++ ) {
            ind = Sync_DomainToGlobal( meshSync, v_i );
            Grid_Lift( vGrid, ind, (unsigned*)inds );
            center = inds[d_i];
            if( center == 0 || center == vGrid->sizes[d_i] - 1 ) {

               /* If we're inside the contact depth range, we need to make
                  sure the side coordinates are aligned. */
	       if( d_i == 0 ) {
		  d_j = 1;
	       }
               else if( d_i == 1 ) {
		  d_j = 0;
	       }
               else if( d_i == 2 ) {
		  d_j = 1;
		  // TODO
		  abort();
	       }
               if( inds[d_j] < gen->contactDepth[d_j][0] )
                  inds[d_j] = gen->contactDepth[d_j][0];
               else if( inds[d_j] > vGrid->sizes[d_j] - gen->contactDepth[d_j][1] - 1 )
                  inds[d_j] = vGrid->sizes[d_j] - gen->contactDepth[d_j][1] - 1;
               Mesh_GetVertex( mesh, v_i )[d_i] =
                  Mesh_GetVertex( mesh, Grid_Project( vGrid, inds ) )[d_i];
            }
         }
      }
   }

   for( d_i = 0; d_i < nDims; d_i++ ) {
      if( !ISet_Has( self->remeshDims, d_i ) )
	 continue;

      Sync_SyncArray( self->syncs[d_i], 
		      mesh->verts[0] + d_i, nDims * sizeof(double), 
		      self->crds[d_i], sizeof(double), 
		      sizeof(double) );

      for( v_i = 0; v_i < nVerts; v_i++ ) {
	 ind = Sync_DomainToGlobal( meshSync, v_i );
	 Grid_Lift( vGrid, ind, (unsigned*)inds );
	 center = inds[d_i];
	 if( center == 0 || center == vGrid->sizes[d_i] - 1 )
	    continue;

	 inds[d_i] = 0;
	 ind = Grid_Project( vGrid, (unsigned*)inds );
	 ind = Sync_GlobalToDomain( self->syncs[d_i], ind );
	 if( ind >= nVerts )
	    leftCrd = self->crds[d_i][ind - nVerts];
	 else
	    leftCrd = Mesh_GetVertex( mesh, ind)[d_i];

	 inds[d_i] = vGrid->sizes[d_i] - 1;
	 ind = Grid_Project( vGrid, (unsigned*)inds );
	 ind = Sync_GlobalToDomain( self->syncs[d_i], ind );
	 if( ind >= nVerts )
	    rightCrd = self->crds[d_i][ind - nVerts];
	 else
	    rightCrd = Mesh_GetVertex( mesh, ind)[d_i];

         /* Do interpolation. */
         if( gen ) {
            if( center <= gen->contactDepth[d_i][0] ) {
               mesh->verts[v_i][d_i] = leftCrd;
               if( gen->contactDepth[d_i][0] ) {
                  mesh->verts[v_i][d_i] +=
                     ((double)center / (double)gen->contactDepth[d_i][0]) *
                     gen->contactGeom[d_i];
               }
            }
            else if( center >= vGrid->sizes[d_i] - gen->contactDepth[d_i][1] - 1 ) {
               mesh->verts[v_i][d_i] = rightCrd;
               if( gen->contactDepth[d_i][1] ) {
                  mesh->verts[v_i][d_i] -=
                     ((double)(vGrid->sizes[d_i] - 1 - center) /
                      (double)gen->contactDepth[d_i][1]) *
                     gen->contactGeom[d_i];
               }
            }
            else {
               mesh->verts[v_i][d_i] = leftCrd + (gen->contactDepth[d_i][0] ? gen->contactGeom[d_i] : 0.0) +
                  ((double)(center - gen->contactDepth[d_i][0]) / 
                   (double)(vGrid->sizes[d_i] - (gen->contactDepth[d_i][0] + gen->contactDepth[d_i][1]) - 1)) *
                  ((rightCrd - leftCrd) - ((gen->contactDepth[d_i][1] ? 1.0 : 0.0) + (gen->contactDepth[d_i][0] ? 1.0 : 0.0)) * gen->contactGeom[d_i]);
            }
         }
         else {

               /* Blend coordinate. */
               mesh->verts[v_i][d_i] = leftCrd + 
                  (double)center * (rightCrd - leftCrd) / 
                  (double)(vGrid->sizes[d_i] - 1);

         }

#if 0
         /* Account for contact depth. */
         if( d_i == 1 ) {
            if( center > self->contactDepth ) {

               /* If we're past contact elements, adjust center to be
                  properly smoothed. */
               center -= self->contactDepth;
               inds[1] = self->contactDepth;
               leftCrd = Mesh_GetVertex( mesh, Grid_Project( vGrid, inds ) )[1];

               /* Blend coordinate. */
               mesh->verts[v_i][d_i] = leftCrd + 
                  (double)center * (rightCrd - leftCrd) / 
                  (double)(vGrid->sizes[d_i] - self->contactDepth - 1);

            }
            else if( center < self->contactDepth ) {

               /* If we're inside the contact depth smooth within the
                  contact range. */
               inds[1] = self->contactDepth;
               rightCrd = Mesh_GetVertex( mesh, Grid_Project( vGrid, inds ) )[1];
               mesh->verts[v_i][d_i] = leftCrd + 
                  (double)center * (rightCrd - leftCrd) / 
                  (double)self->contactDepth;

            }
         }
         else {

               /* Blend coordinate. */
               mesh->verts[v_i][d_i] = leftCrd + 
                  (double)center * (rightCrd - leftCrd) / 
                  (double)(vGrid->sizes[d_i] - 1);

         }
#endif
      }
   }

   /* Meddle with the top right corner. */
   if( self->diffuseCorners ) {
      int nodeIndex, innerIndex, inner2Index;
      double /*grad,*/ *nodeVert, *innerVert, *inner2Vert;

      inds[0] = vGrid->sizes[0] - 1;
      inds[1] = vGrid->sizes[1] - 1;
      nodeIndex = Grid_Project( vGrid, inds );
      nodeVert = Mesh_GetVertex( mesh, nodeIndex );
      inds[0]--;
      innerIndex = Grid_Project( vGrid, inds );
      innerVert = Mesh_GetVertex( mesh, innerIndex );
      inds[0]--;
      inner2Index = Grid_Project( vGrid, inds );
      inner2Vert = Mesh_GetVertex( mesh, inner2Index );

/*
      grad = (innerVert[1] - inner2Vert[1]) / (innerVert[0] - inner2Vert[0]);
      nodeVert[1] = innerVert[1] + grad * (innerVert[1] - nodeVert[1]);
*/
      nodeVert[1] = innerVert[1] - (inner2Vert[1] - innerVert[1]);
   }
   else if( self->diffuseSurface ) {
#if 0
      Mat A;
      Vec b, x;
      KSP ksp;
      int nodeIndex, nDofs;
      double *delta;
      double dt, sep;
      double elMat[2][2], elVec[2], Ni[2], Nx[2], GNx[2], xi;
      double jacDet, weightJacDet;
      int nodeInds[2], rows[2];
      double leftCrd, rightCrd;
      double weight = 1.0;
      int ii, jj, kk, ll;

      nDofs = vGrid->sizes[0];
      VecCreateSeq( PETSC_COMM_SELF, nDofs, &x );
      VecDuplicate( x, &b );
      MatCreateSeqAIJ( PETSC_COMM_SELF, nDofs, nDofs, 3, PETSC_NULL, &A );

      inds[1] = vGrid->sizes[1] - 1;

      for( ii = 0; ii < nDofs - 1; ii++ ) {

         inds[0] = ii;
         nodeInds[0] = Grid_Project( vGrid, inds );
         inds[0]++;
         nodeInds[1] = Grid_Project( vGrid, inds );

         mesh->verts[nodeInds[0]][1] = (double)(ii*ii);
         mesh->verts[nodeInds[1]][1] = (double)((ii+1)*(ii+1));

         memset( elMat[0], 0, 4 * sizeof(double) );
         memset( elVec, 0, 2 * sizeof(double) );

         for( jj = 0; jj < 2; jj++ ) {

            xi = (jj == 0) ? -0.57735026919 : 0.57735026919;
            leftCrd = Mesh_GetVertex( self->mesh, nodeInds[0] )[0];
            rightCrd = Mesh_GetVertex( self->mesh, nodeInds[1] )[0];
            jacDet = (rightCrd - leftCrd) * 0.5;
            weightJacDet = weight * jacDet;
            Ni[0] = 0.5 * (1.0 - xi);
            Ni[1] = 0.5 * (1.0 + xi);
            Nx[0] = -0.5;
            Nx[1] = 0.5;
            GNx[0] = Nx[0] / jacDet;
            GNx[1] = Nx[1] / jacDet;

            for( kk = 0; kk < 2; kk++ ) {
               for ( ll = 0 ; ll < 2; ll++ ) {
                  elMat[kk][ll] += weightJacDet * Ni[kk] * Ni[ll];
                  elVec[kk] += -weightJacDet * GNx[kk] * GNx[ll] * (double)nodeInds[ll] * 2.0; /** GNx[kk] * GNx[ll] *
                                                       Mesh_GetVertex( mesh, nodeInds[ll] )[1];*/
               }

               if( kk == 0 ) {
                  elVec[kk] += -1.0 * (GNx[0] * (double)nodeInds[0] * 2.0 + GNx[1] * (double)nodeInds[1] * 2.0);
/*
                  elVec[kk] += -1.0 * (GNx[0] * Mesh_GetVertex( mesh, nodeInds[0] )[1] +
                                       GNx[1] * Mesh_GetVertex( mesh, nodeInds[1] )[1]);
*/
               }
               else {
                  elVec[kk] += 1.0 * (GNx[0] * (double)nodeInds[0] * 2.0 + GNx[1] * (double)nodeInds[1] * 2.0);
/*
                  elVec[kk] += 1.0 * (GNx[0] * Mesh_GetVertex( mesh, nodeInds[0] )[1] +
                                      GNx[1] * Mesh_GetVertex( mesh, nodeInds[1] )[1]);
*/
               }
            }
         }

         rows[0] = ii; rows[1] = ii + 1;
         MatSetValues( A, 2, rows, 2, rows, elMat[0], ADD_VALUES );
         VecSetValues( b, 2, rows, elVec, ADD_VALUES );

      }

      MatAssemblyBegin( A, MAT_FINAL_ASSEMBLY );
      VecAssemblyBegin( b );
      MatAssemblyEnd( A, MAT_FINAL_ASSEMBLY );
      VecAssemblyEnd( b );

      KSPCreate( PETSC_COMM_SELF, &ksp );
      KSPSetOperators( ksp, A, A, DIFFERENT_NONZERO_PATTERN );
      KSPSolve( ksp, b, x );

      VecGetArray( x, &delta );
      for( inds[0] = 0; inds[0] < vGrid->sizes[0]; inds[0]++ ) {
         nodeIndex = Grid_Project( vGrid, inds );
         mesh->verts[nodeIndex][1] -= self->diffusionCoef * self->ctx->dt * delta[inds[0]];
      }
      VecRestoreArray( x, &delta );

      KSPDestroy( ksp );
      MatDestroy( A );
      VecDestroy( b );
      VecDestroy( x );
#endif


      Mat A;
      Vec b, x;
      KSP ksp;
      int cols[3], nDofs, nodeIndex;
      double rhs, coef, *soln;
      double matVals[3], vecVals[3];
      double d, dt, sep;
      int ii;

      nDofs = vGrid->sizes[0];
      VecCreateSeq( PETSC_COMM_SELF, nDofs, &x );
      VecDuplicate( x, &b );
      MatCreateSeqAIJ( PETSC_COMM_SELF, nDofs, nDofs, 3, PETSC_NULL, &A );

      d = self->diffusionCoef;
      dt = self->ctx->dt;
      sep = Mesh_GetVertex( mesh, 1 )[0] - Mesh_GetVertex( mesh, 0 )[0];
      coef = d * dt / (2.0 * (sep * sep));
      matVals[0] = -coef; matVals[1] = 1.0 + 2.0 * coef; matVals[2] = -coef;
      vecVals[0] = coef;  vecVals[1] = 1.0 - 2.0 * coef; vecVals[2] = coef;
      inds[1] = vGrid->sizes[1] - 1;

      for( ii = 0; ii < nDofs; ii++ ) {
         cols[0] = ii - 1; cols[1] = ii; cols[2] = ii + 1;
         if( ii == 0 ) {
            MatSetValues( A, 1, &ii, 2, cols + 1, matVals + 1, ADD_VALUES );
            MatSetValues( A, 1, &ii, 1, &ii, matVals, ADD_VALUES );
         }
         else if( ii == nDofs - 1 ) {
            MatSetValues( A, 1, &ii, 2, cols, matVals, ADD_VALUES );
            MatSetValues( A, 1, &ii, 1, &ii, matVals + 2, ADD_VALUES );
         }
         else
            MatSetValues( A, 1, &ii, 3, cols, matVals, ADD_VALUES );

         inds[0] = ii;
         nodeIndex = Grid_Project( vGrid, inds );
         if( ii == 0 ) {
            rhs = vecVals[0] * Mesh_GetVertex( mesh, nodeIndex + 1 )[1] +
               vecVals[1] * Mesh_GetVertex( mesh, nodeIndex )[1] + 
               vecVals[2] * Mesh_GetVertex( mesh, nodeIndex )[1];
         }
         else if( ii == nDofs - 1 ) {
            rhs = vecVals[0] * Mesh_GetVertex( mesh, nodeIndex )[1] +
               vecVals[1] * Mesh_GetVertex( mesh, nodeIndex )[1] + 
               vecVals[2] * Mesh_GetVertex( mesh, nodeIndex - 1 )[1];
         }
         else {
            rhs = vecVals[0] * Mesh_GetVertex( mesh, nodeIndex + 1 )[1] +
               vecVals[1] * Mesh_GetVertex( mesh, nodeIndex )[1] + 
               vecVals[2] * Mesh_GetVertex( mesh, nodeIndex - 1 )[1];
         }
         VecSetValues( b, 1, &ii, &rhs, ADD_VALUES );
      }

      MatAssemblyBegin( A, MAT_FINAL_ASSEMBLY );
      VecAssemblyBegin( b );
      MatAssemblyEnd( A, MAT_FINAL_ASSEMBLY );
      VecAssemblyEnd( b );

      KSPCreate( PETSC_COMM_SELF, &ksp );
      KSPSetOperators( ksp, A, A, DIFFERENT_NONZERO_PATTERN );
      KSPSolve( ksp, b, x );

      VecGetArray( x, &soln );
      for( inds[0] = 0; inds[0] < vGrid->sizes[0]; inds[0]++ ) {
         nodeIndex = Grid_Project( vGrid, inds );
         mesh->verts[nodeIndex][1] = soln[inds[0]];
      }
      VecRestoreArray( x, &soln );

      KSPDestroy( ksp );
      MatDestroy( A );
      VecDestroy( b );
      VecDestroy( x );


#if 0
      int nodeIndex, inds[2];
      double *newHeights, sep, d, dt;
      int ii;

      newHeights = AllocArray( double, vGrid->sizes[0] );
      d = self->diffusionCoef;
      dt = self->ctx->dt;

      inds[1] = vGrid->sizes[1] - 1;
      inds[0] = 0;
      nodeIndex = Grid_Project( vGrid, inds );
      sep = Mesh_GetVertex( mesh, nodeIndex + 1 )[0] - Mesh_GetVertex( mesh, nodeIndex )[0];
      newHeights[0] = Mesh_GetVertex( mesh, nodeIndex )[1] -
         2.0 * Mesh_GetVertex( mesh, nodeIndex + 1 )[1] + 
         Mesh_GetVertex( mesh, nodeIndex + 2 )[1];
      newHeights[0] /= sep * sep;
      newHeights[0] *= -d;

      for( inds[0] = 1; inds[0] < vGrid->sizes[0] - 1; inds[0]++ ) {
         nodeIndex = Grid_Project( vGrid, inds );
         sep = (Mesh_GetVertex( mesh, nodeIndex + 1 )[0] - Mesh_GetVertex( mesh, nodeIndex - 1)[0]) * 0.5;
         newHeights[inds[0]] = Mesh_GetVertex( mesh, nodeIndex - 1 )[1] -
            2.0 * Mesh_GetVertex( mesh, nodeIndex )[1] + 
            Mesh_GetVertex( mesh, nodeIndex + 1 )[1];
         newHeights[inds[0]] /= sep * sep;
         newHeights[inds[0]] *= -d;
      }

      nodeIndex = Grid_Project( vGrid, inds );
      sep = -1.0 * (Mesh_GetVertex( mesh, nodeIndex - 1 )[0] - Mesh_GetVertex( mesh, nodeIndex )[0]);
      newHeights[inds[0]] = Mesh_GetVertex( mesh, nodeIndex - 2 )[1] -
         2.0 * Mesh_GetVertex( mesh, nodeIndex - 1 )[1] + 
         Mesh_GetVertex( mesh, nodeIndex )[1];
      newHeights[inds[0]] /= sep * sep;
      newHeights[inds[0]] *= -d;

      for( inds[0] = 0; inds[0] < vGrid->sizes[0]; inds[0]++ ) {
         nodeIndex = Grid_Project( vGrid, inds );
         mesh->verts[nodeIndex][1] += newHeights[inds[0]] * dt;
      }

      MemFree( newHeights );
#endif
   }

   Class_Free( self, inds );
   Mesh_Sync( mesh );
   Mesh_DeformationUpdate( mesh );
}

void RegularRemesher_Build( void* _self ) {
   RegularRemesher* self = Class_Cast( _self, RegularRemesher );
   Mesh* mesh;
   Sync* meshSync;
   const Decomp* meshDecomp;
   ISet* wallSet;
   Grid* vGrid;
   int *inds, ind;
   int nDims, nVerts;
   int nRems, *rems;
   int d_i, w_i, v_i;

   mesh = self->mesh;
   Stg_Component_Build( mesh, NULL, False );
   meshSync = Mesh_GetSync( mesh, 0 );
   meshDecomp = Sync_GetDecomp( meshSync );
   vGrid = *Mesh_GetExtension( mesh, Grid**, "vertexGrid" );
   nDims = Mesh_GetDimSize( mesh );
   self->nDims = nDims;
   nVerts = Mesh_GetLocalSize( mesh, 0 );
   inds = Class_Array( self, int, nDims );

   if( !self->syncs ) {
      self->syncs = Class_Array( self, Sync*, nDims );
      self->crds = Class_Array( self, double*, nDims );
      self->nWallVerts = Class_Array( self, int*, nDims );
      self->wallVerts = Class_Array( self, int**, nDims );
      self->wallCrds = Class_Array( self, double**, nDims );
      for( d_i = 0; d_i < nDims; d_i++ ) {
	 self->nWallVerts[d_i] = Class_Array( self, int, 2 );
	 self->wallVerts[d_i] = Class_Array( self, int*, 2 );
	 self->wallCrds[d_i] = Class_Array( self, double*, 2 );
	 memset( self->nWallVerts[d_i], 0, 2 * sizeof(int) );
	 memset( self->wallVerts[d_i], 0, 2 * sizeof(int*) );
	 memset( self->wallCrds[d_i], 0, 2 * sizeof(double*) );
      }
      memset( self->syncs, 0, nDims * sizeof(Sync*) );
      memset( self->crds, 0, nDims * sizeof(double*) );
   }

   wallSet = ISet_New();
   ISet_SetMaxSize( wallSet, nVerts );
   for( d_i = 0; d_i < nDims; d_i++ ) {
      for( w_i = 0; w_i < 2; w_i++ ) {
	 Class_Free( self, self->wallVerts[d_i][w_i] );
	 Class_Free( self, self->wallCrds[d_i][w_i] );
	 if( !ISet_Has( self->staticWalls[w_i], d_i ) ) {
	    self->nWallVerts[d_i][w_i] = 0;
	    self->wallVerts[d_i][w_i] = NULL;
	    self->wallCrds[d_i][w_i] = NULL;
	    continue;
	 }

	 for( v_i = 0; v_i < nVerts; v_i++ ) {
	    ind = Sync_DomainToGlobal( meshSync, v_i );
	    Grid_Lift( vGrid, ind, (unsigned*)inds );
	    if( (w_i == 0 && inds[d_i] == 0) || 
		(w_i == 1 && inds[d_i] == vGrid->sizes[d_i] - 1) )
	    {
	       ISet_Insert( wallSet, v_i );
	    }
	 }

	 self->nWallVerts[d_i][w_i] = ISet_GetSize( wallSet );
	 self->wallVerts[d_i][w_i] = Class_Array( 
	    self, int, self->nWallVerts[d_i][w_i] );
	 self->wallCrds[d_i][w_i] = Class_Array( 
	    self, double, self->nWallVerts[d_i][w_i] );
	 ISet_GetArray( wallSet, self->wallVerts[d_i][w_i] );
	 ISet_Clear( wallSet );
	 for( v_i = 0; v_i < self->nWallVerts[d_i][w_i]; v_i++ ) {
	    self->wallCrds[d_i][w_i][v_i] = 
	       mesh->verts[self->wallVerts[d_i][w_i][v_i]][d_i];
	 }
      }
   }

   for( d_i = 0; d_i < nDims; d_i++ ) {
      NewClass_RemoveRef( self->syncs[d_i] );
      Class_Free( self, self->crds[d_i] );
      if( !ISet_Has( self->remeshDims, d_i ) ) {
	 self->syncs[d_i] = NULL;
	 self->crds[d_i] = NULL;
	 continue;
      }

      for( v_i = 0; v_i < nVerts; v_i++ ) {
	 ind = Sync_DomainToGlobal( meshSync, v_i );
	 Grid_Lift( vGrid, ind, (unsigned*)inds );
	 inds[d_i] = 0;
	 ind = Grid_Project( vGrid, (unsigned*)inds );
	 if( !Sync_TryGlobalToDomain( meshSync, ind, &ind ) )
	    ISet_TryInsert( wallSet, ind );

	 inds[d_i] = vGrid->sizes[d_i] - 1;
	 ind = Grid_Project( vGrid, (unsigned*)inds );
	 if( !Sync_TryGlobalToDomain( meshSync, ind, &ind ) )
	    ISet_TryInsert( wallSet, ind );
      }

      nRems = ISet_GetSize( wallSet );
      rems = Class_Array( self, int, nRems );
      ISet_GetArray( wallSet, rems );
      ISet_Clear( wallSet );
      self->syncs[d_i] = Sync_New();
      NewClass_AddRef( self->syncs[d_i] );
      Sync_SetDecomp( self->syncs[d_i], meshDecomp );
      Sync_FindRemotes( self->syncs[d_i], nRems, rems );
      Class_Free( self, rems );
      self->crds[d_i] = Class_Array( self, double, nRems );
   }

   NewClass_Delete( wallSet );

#if 0
   /* If we have some contact depth, copy the relevant vertex offsets. */
   if( self->contactDepth > 0 ) {
      int curInd;
      Grid* grid;
      int ii;

      /* Get the vertex grid from the mesh. */
      grid = *Mesh_GetExtension( mesh, Grid**, "vertexGrid" );
      assert( grid );

      /* Allocate for all the contact vertices. */
      nVerts = grid->sizes[0];
      if( nDims == 3 ) nVerts *= grid->sizes[1];
      self->contactVerts = MemArray( double, nVerts, "" );

      /* Copy upper strip. */
      nVerts = Mesh_GetLocalSize( mesh, 0 );
      curInd = 0;
      for( ii = 0; ii < nVerts; ii++ ) {
         Grid_Lift( grid, ii, inds );
         if( inds[1] != self->contactDepth ) continue;

	 /* If we were given a contact size, insert that instead. */
	 if( self->contactSize > 0.0 )
	   self->contactVerts[curInd++] = self->contactSize;
	 else
	   self->contactVerts[curInd++] = Mesh_GetVertex( mesh, ii )[1];
      }

   }
#endif

   Class_Free( self, inds );
}

void RegularRemesher_SetRemeshState( void* _self, int dim, Bool state ) {
   RegularRemesher* self = Class_Cast( _self, RegularRemesher );

   if( state )
      ISet_TryInsert( self->remeshDims, dim );
   else if( ISet_Has( self->remeshDims, dim ) )
      ISet_Remove( self->remeshDims, dim );
}

void RegularRemesher_SetStaticWall( void* _self, int dim, int wall, Bool state ) {
   RegularRemesher* self = Class_Cast( _self, RegularRemesher );

   if( state )
      ISet_TryInsert( self->staticWalls[wall], dim );
   else if( ISet_Has( self->staticWalls[wall], dim ) )
      ISet_Remove( self->staticWalls[wall], dim );
}

#endif


