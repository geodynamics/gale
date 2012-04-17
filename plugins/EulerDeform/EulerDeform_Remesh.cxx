/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street, Melbourne, 3053, Australia.
**
** Authors:
**	Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**	Patrick D. Sunter, Software Engineer, VPAC. (pds@vpac.org)
**	David May, PhD Student Monash University, VPAC. (david.may@sci.maths.monash.edu.au)
**	Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
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
** $Id: LevelSetPlg.c 200 2005-07-08 08:24:41Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include <Underworld/Underworld.h>

#include "types.h"
#include "Context.h"
#include "EulerDeform.h"

void EulerDeform_Remesh(TimeIntegrand* crdAdvector, EulerDeform_Context* edCtx)
{
  /* Do not remesh if this is the first step */
  if(crdAdvector->context->timeStep==crdAdvector->context->restartTimestep)
    return;

  Mesh_Algorithms	*tmpAlgs, *oldAlgs;
  unsigned	sys_i;

  assert( edCtx );

  /* We do the second system first, because that is the vertex
     centered system.  We want the cell centered and vertex
     centered systems to be compatible, so the cell centered
     system is based on the vertex centered one. */

  for( sys_i = 1; sys_i < edCtx->nSystems+1; sys_i++ ) {
    EulerDeform_System*	sys = edCtx->systems
      + sys_i%edCtx->nSystems;
    double**		oldCrds;
    double**		newCrds;
    unsigned		nDomainNodes;
    unsigned		nDims;
    unsigned		n_i;
    Grid *grid =
      *(Grid**)ExtensionManager_Get(sys->mesh->info, sys->mesh, 
                                    ExtensionManager_GetHandle( sys->mesh->info,
                                                                "vertexGrid" ));
    /* Set the displacement field. */
    if(sys->dispField) {
      double disp[3];
      int num_verts, num_dims;
      int ii, jj;

      num_dims = Mesh_GetDimSize(sys->mesh);
      num_verts = Mesh_GetLocalSize(sys->mesh, MT_VERTEX);
      for(ii = 0; ii < num_verts; ii++) {
        for(jj = 0; jj < num_dims; jj++)
          disp[jj] = sys->verts[ii*num_dims + jj] - sys->mesh->verts[ii][jj];
        FeVariable_SetValueAtNode(sys->dispField, ii, disp);
      }
    }

    nDims = Mesh_GetDimSize( sys->mesh );
    
    /* Update all local coordinates. */
    for( n_i = 0; n_i < Mesh_GetLocalSize( sys->mesh, MT_VERTEX ); n_i++ )
      memcpy( sys->mesh->verts[n_i], sys->verts + n_i * nDims,
              nDims * sizeof(double));
    
    /* Revert side coordinates if required. */
    if( sys->staticSides ) {
      IndexSet	*tmpIndSet;
      unsigned	nInds, *inds;
      unsigned	ind_i;
      
      /* Collect indices of all the sides. */
      
      tmpIndSet = EulerDeform_CreateStaticSet(sys);
      IndexSet_GetMembers( tmpIndSet, &nInds, &inds );
      
      /* Copy back coords. */
      for( ind_i = 0; ind_i < nInds; ind_i++ )
        {
          std::cout << "side "
                    << ind_i << " "
                    << inds[ind_i] << " "
                    << sys->sideCoords[ind_i][0] << " "
                    << sys->sideCoords[ind_i][1] << " "
                    << sys->sideCoords[ind_i][2] << " "
                    << "\n";
        memcpy( sys->mesh->verts[inds[ind_i]], sys->sideCoords[ind_i],
                nDims * sizeof(double));
        }
      FreeObject( tmpIndSet );
      FreeArray( sys->sideCoords );
      
      if(sys->wrapTop)
        {
          if(sys->staticLeft && !sys->staticLeftTop && !sys->floatLeftTop)
            {
              EulerDeform_Remesh_Corner(sys->mesh,0,1,sys->static_left_coord,
                                        0,1,2);
            }
          if(sys->staticRight && !sys->staticRightTop && !sys->floatRightTop)
            {
              EulerDeform_Remesh_Corner(sys->mesh,grid->sizes[0]-1,
                                        grid->sizes[0]-2,
                                        sys->static_right_coord,0,1,2);
            }
          if(sys->staticFront && !sys->staticTopFront)
            {
              EulerDeform_Remesh_Corner(sys->mesh,0,1,sys->static_front_coord,2,1,0);
            }
          if(sys->staticBack && !sys->staticTopBack)
            {
              EulerDeform_Remesh_Corner(sys->mesh,grid->sizes[2]-1,
                                        grid->sizes[2]-2,sys->static_back_coord,2,1,0);
            }
        }

    }

    /* If we have regular mesh algorithms specified, set the
       algorithms temporarily to an irregular method. */
    if( !strcmp( sys->mesh->algorithms->type, "Mesh_RegularAlgorithms")
        && sys->remesher)
      {
        tmpAlgs = Mesh_Algorithms_New( "", edCtx->ctx );
        oldAlgs = sys->mesh->algorithms;
        sys->mesh->algorithms = NULL;
        Mesh_SetAlgorithms( sys->mesh, tmpAlgs );
      }
    else
      tmpAlgs = NULL;

    /* Every system should synchronise the mesh coordinates. */
    Mesh_Sync( sys->mesh );
    Mesh_DeformationUpdate( sys->mesh );

    /* Only if remesher specified. */
    if( !sys->remesher ) {
      continue;
    }

    /* If a remesh interval is requested, check now. */
    if(sys->interval>0 && edCtx->ctx->timeStep % sys->interval>0)
      {
        Journal_Printf( Underworld_Info,
                        "*** EulerDeform: Not remeshing this timestep.\n" );
        continue;
      }
    Journal_Printf( Underworld_Info, "*** EulerDeform: Remeshing.\n" );

    /* Shrink wrap the top/bottom surface. */
    if( sys->wrapTop )
      EulerDeform_WrapSurface( sys, sys->mesh->verts, true );
    if( sys->wrapBottom )
      EulerDeform_WrapSurface( sys, sys->mesh->verts, false );

    /* Float the top left and right corners if needed. */
    if(sys->floatLeftTop)
      EulerDeform_FloatLeftTop(sys,grid,sys->mesh->verts);
    if(sys->floatRightTop)
      EulerDeform_FloatRightTop(sys,grid,sys->mesh->verts);

    /* Store old coordinates. */
    nDomainNodes = FeMesh_GetNodeDomainSize( sys->mesh );
    oldCrds = AllocArray2D( double, nDomainNodes, nDims );
    for( n_i = 0; n_i < nDomainNodes; n_i++ )
      memcpy( oldCrds[n_i], sys->mesh->verts[n_i], nDims * sizeof(double) );

    for( n_i = 0; n_i < nDomainNodes; n_i++ )
      std::cout << "coord "
                << n_i << " "
                << sys->mesh->verts[n_i][0] << " "
                << sys->mesh->verts[n_i][1] << " "
                << sys->mesh->verts[n_i][2] << " "
                << "\n";

    /* Remesh the system. */
    Stg_Component_Execute( sys->remesher, NULL, True );
    Mesh_Sync( sys->mesh );

    for( n_i = 0; n_i < nDomainNodes; n_i++ )
      std::cout << "Before "
                << n_i << " "
                << sys->mesh->verts[n_i][0] << " "
                << sys->mesh->verts[n_i][1] << " "
                << sys->mesh->verts[n_i][2] << " "
                << "\n";

    /* Swap old coordinates back in temporarily. */
    newCrds = sys->mesh->verts;
    sys->mesh->verts = oldCrds;

    /* Update the displacement field for the remeshing. */
    if(sys->dispField)
      {
        double disp[3];
        for(n_i = 0; n_i<nDomainNodes; n_i++)
          {
            FeVariable_GetValueAtNode(sys->dispField, n_i, disp);
            for(unsigned dof_i = 0; dof_i<nDims; dof_i++)
              {
                disp[dof_i] += newCrds[n_i][dof_i] - oldCrds[n_i][dof_i];
              }
            FeVariable_SetValueAtNode(sys->dispField, n_i, disp);
          }
      }

    /* Swap back coordinates and free memory. */
    sys->mesh->verts = newCrds;
    FreeArray( oldCrds );
    
    /* Swap back old algorithms. */
    if( tmpAlgs ) {
      Mesh_SetAlgorithms( sys->mesh, oldAlgs );
    }

    /* Re-sync with new coordinates. */
    Mesh_Sync( sys->mesh );
    Mesh_DeformationUpdate( sys->mesh );

    /* Reset the coordinates of the pressure meshes. */

    if(sys->p_mesh!=NULL)
      InnerGenerator_SetCoordinates((InnerGenerator*)(sys->p_mesh->generator),
                                    (FeMesh*)(sys->p_mesh));
  }
}
