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

void _EulerDeform_TriBarycenter(double** tri, const double* pnt, double* dst)
{
  double a = tri[0][0] - tri[2][0];
  double b = tri[1][0] - tri[2][0];
  double c = tri[2][0] - pnt[0];
  double d = tri[0][1] - tri[2][1];
  double e = tri[1][1] - tri[2][1];
  double f = tri[2][1] - pnt[1];
  double g = tri[0][2] - tri[2][2];
  double h = tri[1][2] - tri[2][2];
  double i = tri[2][2] - pnt[2];

  dst[0] = (b * (f + i) - c * (e + h)) / (a * (e + h) - b * (d + g));
  dst[1] = (a * (f + i) - c * (d + g)) / (b * (d + g) - a * (e + h));
  dst[2] = 1.0 - dst[0] - dst[1];
}


Bool _EulerDeform_QuadYInterp(double** crds, const double* pnt, double* val)
{
  double* modCrds[4];
  double modCrds0[2], modCrds1[2], modCrds2[2], modCrds3[2];
  unsigned* inds[2];
  unsigned inds0[3], inds1[3];
  unsigned inc[4];
  double modPnt[3];
  unsigned inside;
  double bc[3];

  modCrds[0] = modCrds0;
  modCrds[1] = modCrds1;
  modCrds[2] = modCrds2;
  modCrds[3] = modCrds3;
  modCrds[0][0] = crds[0][0]; modCrds[0][1] = crds[0][2];
  modCrds[1][0] = crds[1][0]; modCrds[1][1] = crds[1][2];
  modCrds[2][0] = crds[2][0]; modCrds[2][1] = crds[2][2];
  modCrds[3][0] = crds[3][0]; modCrds[3][1] = crds[3][2];
  modPnt[0] = pnt[0]; modPnt[1] = pnt[2];

  inds[0] = inds0;
  inds[1] = inds1;
  inds[0][0] = 0; inds[0][1] = 1; inds[0][2] = 2;
  inds[1][0] = 1; inds[1][1] = 3; inds[1][2] = 2;
  inc[0] = 0; inc[1] = 1; inc[2] = 2; inc[3] = 3;

  if( Simplex_Search2D( modCrds, inc, 2, inds, 
                        modPnt, bc, &inside ) )
    {
      *val = bc[0] * crds[inds[inside][0]][1]
        + bc[1] * crds[inds[inside][1]][1]
        + bc[2] * crds[inds[inside][2]][1];
      return True;
    }
  else
    return False;
}


