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

/* This creates a set which makes sure not to include the corners if
   we have not made both sides of the corner static. */

IndexSet* EulerDeform_CreateStaticSet(EulerDeform_System* sys)
{
  Grid*		grid;
  unsigned	nNodes;
  unsigned	n_i;
  IndexSet	*set;
  IJK			ijk;

  grid = *(Grid**)ExtensionManager_Get ( sys->mesh->info, sys->mesh, ExtensionManager_GetHandle( sys->mesh->info, (Name)"vertexGrid" )  );

  nNodes = Mesh_GetDomainSize( sys->mesh, MT_VERTEX );
  set = IndexSet_New( nNodes );

  for( n_i = 0; n_i < nNodes; n_i++ ) {
    Bool add;
    add=False;
    RegularMeshUtils_Node_1DTo3D
      ( sys->mesh, Mesh_DomainToGlobal( sys->mesh, MT_VERTEX, n_i ), ijk );

    /* 2D */
    if(sys->mesh->topo->nDims == 2)
      {
        /* Left side and Corner */
        if(ijk[0]==0)
          {
            if(ijk[1]==0)
              {
                if(sys->staticLeftBottom)
                  add=True;
              }
            else if(ijk[1]==grid->sizes[1]-1)
              {
                if(sys->staticLeftTop)
                  add=True;
              }
            else if(sys->staticLeft)
              add=True;
          }
        /* Right side and corner */
        else if(ijk[0]==grid->sizes[0]-1)
          {
            if(ijk[1]==0)
              {
                if(sys->staticRightBottom)
                  add=True;
              }
            else if(ijk[1]==grid->sizes[1]-1)
              {
                if(sys->staticRightTop)
                  add=True;
              }
            else if(sys->staticRight)
              add=True;
          }
        /* Top and Bottom */
        else if((ijk[1]==0 && sys->staticBottom)
                || (ijk[1]==grid->sizes[1]-1 && sys->staticTop))
          add=True;
      }
    /* 3D */
    else if(sys->mesh->topo->nDims == 3)
      {
        /* Left side */
        if(ijk[0]==0)
          {
            /* Left Bottom */
            if(ijk[1]==0)
              {
                if(ijk[2]==0)
                  {
                    if(sys->staticLeftBottomBack)
                      add=True;
                  }
                else if(ijk[2]==grid->sizes[2]-1)
                  {
                    if(sys->staticLeftBottomFront)
                      add=True;
                  }
                else if(sys->staticLeftBottom)
                  add=True;
              }
            /* Left Top */
            else if(ijk[1]==grid->sizes[1]-1)
              {
                if(ijk[2]==0)
                  {
                    if(sys->staticLeftTopBack)
                      add=True;
                  }
                else if(ijk[2]==grid->sizes[2]-1)
                  {
                    if(sys->staticLeftTopFront)
                      add=True;
                  }
                else if(sys->staticLeftTop)
                  add=True;
              }
            /* Left Back */
            else if(ijk[2]==0)
              {
                if(sys->staticLeftBack)
                  add=True;
              }
            /* Left Front */
            else if(ijk[2]==grid->sizes[2]-1)
              {
                if(sys->staticLeftFront)
                  add=True;
              }
            /* Left */
            else if(sys->staticLeft)
              add=True;
          }
        /* Right side */
        else if(ijk[0]==grid->sizes[0]-1)
          {
            /* Right Bottom */
            if(ijk[1]==0)
              {
                if(ijk[2]==0)
                  {
                    if(sys->staticRightBottomBack)
                      add=True;
                  }
                else if(ijk[2]==grid->sizes[2]-1)
                  {
                    if(sys->staticRightBottomFront)
                      add=True;
                  }
                else if(sys->staticRightBottom)
                  add=True;
              }
            /* Right Top */
            else if(ijk[1]==grid->sizes[1]-1)
              {
                if(ijk[2]==0)
                  {
                    if(sys->staticRightTopBack)
                      add=True;
                  }
                else if(ijk[2]==grid->sizes[2]-1)
                  {
                    if(sys->staticRightTopFront)
                      add=True;
                  }
                else if(sys->staticRightTop)
                  add=True;
              }
            /* Right Back */
            else if(ijk[2]==0)
              {
                if(sys->staticRightBack)
                  add=True;
              }
            /* Right Front */
            else if(ijk[2]==grid->sizes[2]-1)
              {
                if(sys->staticRightFront)
                  add=True;
              }
            /* Right */
            else if(sys->staticRight)
              add=True;
          }
        /* Bottom */
        else if(ijk[1]==0)
          {
            if(ijk[2]==0)
              {
                if(sys->staticBottomBack)
                  add=True;
              }
            else if(ijk[2]==grid->sizes[2]-1)
              {
                if(sys->staticBottomFront)
                  add=True;
              }
            else if(sys->staticBottom)
              add=True;
          }
        /* Top */
        else if(ijk[1]==grid->sizes[1]-1)
          {
            if(ijk[2]==0)
              {
                if(sys->staticTopBack)
                  add=True;
              }
            else if(ijk[2]==grid->sizes[2]-1)
              {
                if(sys->staticTopFront)
                  add=True;
              }
            else if(sys->staticTop)
              add=True;
          }
        /*  Front and Back */
        else if((ijk[2]==0 && sys->staticBack)
                || (ijk[2]==grid->sizes[2]-1 && sys->staticFront))
          add=True;
      }

    if( add )
      IndexSet_Add( set, n_i );
  }
  return set;
}
