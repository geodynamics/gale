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

void EulerDeform_FloatLeftTop(EulerDeform_System* sys, Grid *grid,
                              double** crds)
{
  int i, max_z;
  IJK ijk, inside_ijk;
  unsigned ind, nLocalNodes, inside;

  nLocalNodes = Mesh_GetLocalSize( sys->mesh, MT_VERTEX );
  ijk[0]=0;
  ijk[1]=grid->sizes[1]-1;

  if(Mesh_GetDimSize(sys->mesh)==2)
    {
      max_z=1;
    }
  else
    {
      max_z=grid->sizes[2];
    }

  for(i=0;i<max_z;++i)
    {
      ijk[2]=i;
      ind=Grid_Project(grid,ijk);
      if( !(!Mesh_GlobalToDomain( sys->mesh, MT_VERTEX, ind, &ind )
            || ind >= nLocalNodes ))
        {
          inside_ijk[0]=ijk[0]+1;
          inside_ijk[1]=ijk[1];
          inside_ijk[2]=ijk[2];
          inside=Grid_Project(grid,inside_ijk);
          if(!(!Mesh_GlobalToDomain( sys->mesh, MT_VERTEX, inside, &inside )
               || inside >= nLocalNodes ))
            crds[ind][1]=crds[inside][1];
        }
    }
}


