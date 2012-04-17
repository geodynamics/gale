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

/* Remesh the top corners in 2D or 3D */

void EulerDeform_Remesh_Corner(Mesh *mesh, const int &corner, const int &inside,
                               const double &side_coord, const int &boundary_dim,
                               const int &height_dim, const int &tangent_dim)
{
  IJK		ijk;
  Grid *grid;
  grid =
    *(Grid**)ExtensionManager_Get( mesh->info, mesh, 
                                   ExtensionManager_GetHandle( mesh->info,
                                                               "vertexGrid" ) );
  ijk[height_dim]=grid->sizes[height_dim]-1;
  if(Mesh_GetDimSize(mesh)==2)
    {
      unsigned n_corner, n_interior, n, n_in;
      ijk[boundary_dim]=corner;
      n_corner=RegularMeshUtils_Node_3DTo1D(mesh,ijk);
      ijk[boundary_dim]=inside;
      n_interior=RegularMeshUtils_Node_3DTo1D(mesh,ijk);
      if(Mesh_GlobalToDomain(mesh,MT_VERTEX,
                             n_corner,&n)
         && Mesh_GlobalToDomain(mesh,MT_VERTEX,
                                n_interior,&n_in))
        {
          double *crds[2];

          crds[0]=mesh->verts[n];
          crds[1]=mesh->verts[n_in];
          if(!_EulerDeform_LineInterp(crds,side_coord,boundary_dim,height_dim,
                                       &(mesh->verts[n][height_dim])))
            {
              printf("The side is moving in the wrong direction.\n");
              printf("%g %g %g %g %g\n",side_coord,crds[0][0],crds[0][1],
                     crds[1][0],crds[1][1]);
              abort();
            }
          mesh->verts[n][boundary_dim]=side_coord;
        }
    }
  else /* 3D */
    {
      for(ijk[tangent_dim]=0; ijk[tangent_dim]<grid->sizes[tangent_dim]; ++ijk[tangent_dim])
        {
          unsigned n_corner, n_interior, n, n_in;
          ijk[boundary_dim]=corner;
          n_corner=RegularMeshUtils_Node_3DTo1D(mesh,ijk);
          ijk[boundary_dim]=inside;
          n_interior=RegularMeshUtils_Node_3DTo1D(mesh,ijk);
          if(Mesh_GlobalToDomain(mesh,MT_VERTEX,n_corner,&n)
             && Mesh_GlobalToDomain(mesh,MT_VERTEX,n_interior,&n_in))
            {
              double *crds[2];
              crds[0]=mesh->verts[n];
              crds[1]=mesh->verts[n_in];

              if(!_EulerDeform_LineInterp(crds,side_coord,boundary_dim,
                                          height_dim,
                                          &(mesh->verts[n][height_dim])))
                {
                  printf("The side is moving in the wrong direction.\n");
                  printf("%g %g %g %g %g %g %g\n",side_coord,crds[0][0],
                         crds[0][1],crds[0][2],crds[1][0],crds[1][1],
                         crds[1][2]);
                  abort();
                }
              mesh->verts[n][boundary_dim]=side_coord;
            }
        }
    }
}
