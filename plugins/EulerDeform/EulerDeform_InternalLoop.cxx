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

void EulerDeform_InternalLoop(EulerDeform_System* sys, Grid* grm,
                              double** oldCrds, unsigned* ijk,
                              unsigned curDim, const bool &top)
{
  unsigned nDims;
  XYZ newCrd, oldCrd;
  double* crds[4];
  double crds0[3], crds1[3], crds2[3], crds3[3];
  unsigned centerInd;
  Mesh* mesh;
  unsigned nLocalNodes;
  unsigned ind;
  double fudge_factor;
        
  /* fudge_factor nudges the coordinate inside the mesh a little
     to avoid failed interpolations */
  fudge_factor=1.0e-10;

  if(curDim<grm->nDims)
    {
      if(curDim==1)
        {
          if(top)
            {
              ijk[1]=grm->sizes[curDim] - 1;
            }
          else
            {
              ijk[1]=0;
            }
          EulerDeform_InternalLoop(sys,grm,oldCrds,ijk,curDim + 1,top);
        }
      else
        {
          for(ijk[curDim]=0; ijk[curDim]<grm->sizes[curDim]; ijk[curDim]++)
            {
              EulerDeform_InternalLoop(sys,grm,oldCrds,ijk,curDim + 1,top);
            }
        }
    }
  else
    {
      if(grm->nDims==2)
        {
          mesh = sys->mesh;
          nDims = Mesh_GetDimSize( mesh );

          crds[0] = crds0;
          crds[1] = crds1;
          nLocalNodes = Mesh_GetLocalSize( mesh, MT_VERTEX );

          /* Skip corners. */
          if(ijk[0]==0 || ijk[0]==grm->sizes[0] - 1)
            {
              return;
            }

          /* Get old and new coordinate. */
          centerInd = Grid_Project( grm, ijk );
          if(!Mesh_GlobalToDomain(mesh,MT_VERTEX,centerInd,&centerInd)
             || centerInd >= nLocalNodes)
            return;

          newCrd[0] = mesh->verts[centerInd][0];
          newCrd[1] = mesh->verts[centerInd][1];
          oldCrd[0] = oldCrds[centerInd][0];
          oldCrd[1] = oldCrds[centerInd][1];

          /* Are we left or right? */
          if(newCrd[0]<oldCrd[0])
            {
              ijk[0]--; ind = Grid_Project( grm, ijk ); ijk[0]++;
              insist( Mesh_GlobalToDomain( mesh, MT_VERTEX, ind, &ind ), != 0 );
              memcpy( crds[0], oldCrds[ind], nDims * sizeof(double) );
              memcpy( crds[1], oldCrd, nDims * sizeof(double) );
            }
          else if(newCrd[0]>oldCrd[0])
            {
              ijk[0]++; ind = Grid_Project( grm, ijk ); ijk[0]--;
              insist( Mesh_GlobalToDomain( mesh, MT_VERTEX, ind, &ind ), != 0 );
              memcpy( crds[1], oldCrds[ind], nDims * sizeof(double) );
              memcpy( crds[0], oldCrd, nDims * sizeof(double) );
            }

          if(newCrd[0]==oldCrd[0])
            {
              mesh->verts[centerInd][1]=oldCrd[1];
            }
          else
            {

              /* Interpolate. */
#ifndef NDEBUG
              Bool result=
#endif
                _EulerDeform_LineInterp(crds,newCrd[0],0,1,
                                        &mesh->verts[centerInd][1]);
#ifndef NDEBUG
              assert(result);
#endif
              if((mesh->verts[centerInd][1]>0) ^ top)
                {
                  mesh->verts[centerInd][1] *= 1+fudge_factor;
                }
              else
                {
                  mesh->verts[centerInd][1] *= 1-fudge_factor;
                }
            }
        }
      else if(grm->nDims==3)
        {
          mesh = sys->mesh;
          nDims = Mesh_GetDimSize( mesh );

          crds[0] = crds0; crds[1] = crds1; crds[2] = crds2; crds[3] = crds3;
          nLocalNodes = Mesh_GetLocalSize( mesh, MT_VERTEX );

          /* Skip corners. */
          if( (ijk[0] == 0 || ijk[0] == grm->sizes[0] - 1) && 
              (ijk[2] == 0 || ijk[2] == grm->sizes[2] - 1))
            {
              return;
            }

          /* Get old and new coordinate. */
          centerInd = Grid_Project( grm, ijk );
          if(!Mesh_GlobalToDomain(mesh,MT_VERTEX,centerInd,&centerInd)
             || centerInd >= nLocalNodes)
            return;

          newCrd[0] = mesh->verts[centerInd][0];
          newCrd[1] = mesh->verts[centerInd][1];
          newCrd[2] = mesh->verts[centerInd][2];
          oldCrd[0] = oldCrds[centerInd][0];
          oldCrd[1] = oldCrds[centerInd][1];
          oldCrd[2] = oldCrds[centerInd][2];

          /* Handle internal nodes. */
          if( ijk[0] > 0 && ijk[2] > 0 ) {
            ijk[0]--;
            ijk[2]--;
            ind=Grid_Project(grm,ijk);
            insist(Mesh_GlobalToDomain(mesh,MT_VERTEX,ind,&ind ), != 0);
            memcpy(crds[0],oldCrds[ind],nDims*sizeof(double));

            ijk[0]++;
            ind=Grid_Project(grm,ijk);
            insist(Mesh_GlobalToDomain(mesh,MT_VERTEX,ind,&ind), != 0);
            memcpy(crds[1],oldCrds[ind],nDims*sizeof(double));

            ijk[0]--;
            ijk[2]++;
            ind=Grid_Project(grm,ijk);
            insist(Mesh_GlobalToDomain(mesh,MT_VERTEX,ind,&ind), != 0);
            memcpy(crds[2],oldCrds[ind],nDims*sizeof(double));

            ijk[0]++;
            ind=Grid_Project(grm,ijk);
            insist(Mesh_GlobalToDomain(mesh,MT_VERTEX,ind,&ind), != 0);
            memcpy(crds[3],oldCrds[ind],nDims*sizeof(double));

            if(_EulerDeform_QuadYInterp(crds,newCrd,&mesh->verts[centerInd][1]))
              {
                if((mesh->verts[centerInd][1]>0) ^ top)
                  {
                    mesh->verts[centerInd][1] *= 1+fudge_factor;
                  }
                else
                  {
                    mesh->verts[centerInd][1] *= 1-fudge_factor;
                  }
                return;
              }
          }

          if(ijk[0]>0 && ijk[2]<grm->sizes[2] - 1 )
            {
              ijk[0]--;
              ind=Grid_Project(grm,ijk);
              insist( Mesh_GlobalToDomain( mesh, MT_VERTEX, ind, &ind ), != 0 );
              memcpy( crds[0], oldCrds[ind], nDims * sizeof(double) );

              ijk[0]++;
              ind=Grid_Project(grm,ijk);
              insist( Mesh_GlobalToDomain( mesh, MT_VERTEX, ind, &ind ), != 0 );
              memcpy( crds[1], oldCrds[ind], nDims * sizeof(double) );

              ijk[0]--;
              ijk[2]++;
              ind=Grid_Project(grm,ijk);
              insist( Mesh_GlobalToDomain( mesh, MT_VERTEX, ind, &ind ), != 0 );
              memcpy( crds[2], oldCrds[ind], nDims * sizeof(double) );

              ijk[0]++;
              ind=Grid_Project(grm,ijk);
              insist( Mesh_GlobalToDomain( mesh, MT_VERTEX, ind, &ind ), != 0 );
              memcpy( crds[3], oldCrds[ind], nDims * sizeof(double) );

              ijk[2]--;
              if(_EulerDeform_QuadYInterp(crds,newCrd,
                                          &mesh->verts[centerInd][1]))
                {
                  if((mesh->verts[centerInd][1]>0) ^ top)
                    {
                      mesh->verts[centerInd][1] *= 1+fudge_factor;
                    }
                  else
                    {
                      mesh->verts[centerInd][1] *= 1-fudge_factor;
                    }
                  return;
                }
            }

          if(ijk[0]<grm->sizes[0]-1 && ijk[2]>0)
            {
              ijk[2]--;
              ind=Grid_Project(grm,ijk);
              insist( Mesh_GlobalToDomain( mesh, MT_VERTEX, ind, &ind ), != 0 );
              memcpy( crds[0], oldCrds[ind], nDims * sizeof(double) );

              ijk[0]++;
              ind=Grid_Project(grm,ijk);
              insist( Mesh_GlobalToDomain( mesh, MT_VERTEX, ind, &ind ), != 0 );
              memcpy( crds[1], oldCrds[ind], nDims * sizeof(double) );

              ijk[0]--;
              ijk[2]++;
              ind=Grid_Project(grm,ijk);
              insist( Mesh_GlobalToDomain( mesh, MT_VERTEX, ind, &ind ), != 0 );
              memcpy( crds[2], oldCrds[ind], nDims * sizeof(double) );

              ijk[0]++;
              ind=Grid_Project(grm,ijk);
              insist( Mesh_GlobalToDomain( mesh, MT_VERTEX, ind, &ind ), != 0 );
              memcpy( crds[3], oldCrds[ind], nDims * sizeof(double) );

              ijk[0]--;
              if(_EulerDeform_QuadYInterp(crds,newCrd,
                                          &mesh->verts[centerInd][1]))
                {
                  if((mesh->verts[centerInd][1]>0) ^ top)
                    {
                      mesh->verts[centerInd][1] *= 1+fudge_factor;
                    }
                  else
                    {
                      mesh->verts[centerInd][1] *= 1-fudge_factor;
                    }
                  return;
                }
            }

          if(ijk[0]<grm->sizes[0]-1 && ijk[2]<grm->sizes[2]-1)
            {
              ind=Grid_Project(grm,ijk);
              insist( Mesh_GlobalToDomain( mesh, MT_VERTEX, ind, &ind ), != 0 );
              memcpy( crds[0], oldCrds[ind], nDims * sizeof(double) );

              ijk[0]++;
              ind=Grid_Project(grm,ijk);
              insist( Mesh_GlobalToDomain( mesh, MT_VERTEX, ind, &ind ), != 0 );
              memcpy( crds[1], oldCrds[ind], nDims * sizeof(double) );

              ijk[0]--;
              ijk[2]++;
              ind=Grid_Project(grm,ijk);
              insist( Mesh_GlobalToDomain( mesh, MT_VERTEX, ind, &ind ), != 0 );
              memcpy( crds[2], oldCrds[ind], nDims * sizeof(double) );

              ijk[0]++;
              ind=Grid_Project(grm,ijk);
              insist( Mesh_GlobalToDomain( mesh, MT_VERTEX, ind, &ind ), != 0 );
              memcpy( crds[3], oldCrds[ind], nDims * sizeof(double) );

              ijk[0]--;
              ijk[2]--;
              if(_EulerDeform_QuadYInterp(crds,newCrd,
                                          &mesh->verts[centerInd][1]))
                {
                  if((mesh->verts[centerInd][1]>0) ^ top)
                    {
                      mesh->verts[centerInd][1] *= 1+fudge_factor;
                    }
                  else
                    {
                      mesh->verts[centerInd][1] *= 1-fudge_factor;
                    }
                  return;
                }
            }
          abort();
        }
      else
        {
          abort();
        }
    }
}
