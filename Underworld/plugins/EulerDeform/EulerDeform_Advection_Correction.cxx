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

extern ExtensionInfo_Index EulerDeform_ContextHandle;

namespace {
  void restore_velocity(FeVariable* velocityField, double* oldVelocity)
  {
    FeMesh* mesh=velocityField->feMesh;
    int lNodeCount=FeMesh_GetNodeLocalSize(mesh);
    int dof=velocityField->fieldComponentCount;
    double oldV[3];
    for(int lNode_I = 0; lNode_I<lNodeCount; lNode_I++)
      {
        memcpy(oldV,&oldVelocity[lNode_I*dof],dof*sizeof(double));
        FeVariable_SetValueAtNode(velocityField,lNode_I,oldV);
      }
  }

  void add_correction(FeVariable* velocityField,
                      double* artVelocity)
  {
    FeMesh* mesh=velocityField->feMesh;
    int lNodeCount=FeMesh_GetNodeLocalSize(mesh);
    int dof=velocityField->fieldComponentCount;

    for(int lNode_I=0; lNode_I<lNodeCount; lNode_I++)
      {
        double v[3];
        FeVariable_GetValueAtNode(velocityField,lNode_I,v);            
        for(int i=0;i<dof;i++)
          {
            v[i]-=artVelocity[lNode_I*dof+i];
          }
        FeVariable_SetValueAtNode(velocityField,lNode_I,v);
      }
  }

  void store_current_velocity(FeVariable *velocityField, double *oldVelocity)
  {
    /* save the current values of the velocity field in the oldVelocity array */
    FeMesh *mesh = velocityField->feMesh;
    unsigned dof = velocityField->fieldComponentCount;
    unsigned numLocalNodes = FeMesh_GetNodeLocalSize( mesh );
    double vel[3];

    for(unsigned lNode_I=0; lNode_I<numLocalNodes; lNode_I++)
      {
        FeVariable_GetValueAtNode(velocityField,lNode_I,vel);
        memcpy(&oldVelocity[lNode_I*dof],vel,dof*sizeof(double));
      }
  }

  void compute_velocity_from_displacement(FeVariable *artDField,
                                          double *artVelocity,
                                          double dt)
  {
    /* save the purely artificial bit of the remeshing in the artVelocity */
    FeMesh* mesh=artDField->feMesh;
    Dof_Index dof=artDField->fieldComponentCount;
    double artV[3], artD[3];
    unsigned numLocalNodes=FeMesh_GetNodeLocalSize(mesh);
    Dof_Index dof_I, lNode_I;

    /* INITIAL CONDITION: artV = 0 */
    if(dt==0)
      {
        for(lNode_I=0; lNode_I<numLocalNodes; lNode_I++)
          {
            for(dof_I=0; dof_I<dof; dof_I++) 
              artV[dof_I]=0;

            memcpy( &artVelocity[lNode_I*dof] , artV, dof*sizeof(double) );
          }
        return;
      }
	
    /* artV = artD / dt. */
    for(lNode_I=0; lNode_I<numLocalNodes; lNode_I++)
      {
        FeVariable_GetValueAtNode(artDField,lNode_I,artD);
        for(dof_I=0; dof_I<dof; dof_I++) 
          artV[dof_I]= artD[dof_I]/dt;
        memcpy(&artVelocity[lNode_I*dof],artV,dof*sizeof(double));
      }
  }
}

void EulerDeform_Advection_Correction(void* sle, void* data)
{
  UnderworldContext* context=(UnderworldContext*)data;

  if(context->timeStep==context->restartTimestep)
    return;

  EulerDeform_Context* edCtx=(EulerDeform_Context*)
    ExtensionManager_Get(context->extensionMgr,context,
                         EulerDeform_ContextHandle);

  EulerDeform_System* sys;
  for(unsigned sys_i=0; sys_i<edCtx->nSystems; sys_i++)
    {
      sys=edCtx->systems + sys_i;
      if(sys->dispField)
        break;
    }

  FeVariable* velocityField=(FeVariable*)
    LiveComponentRegister_Get(context->CF->LCRegister,(Name)"VelocityField");
  double dt=context->dt;
  double *artVelocity(NULL), *oldVelocity(NULL);
  int lNodeCount=FeMesh_GetNodeLocalSize(velocityField->feMesh);

  /* store the current velocity in oldVelocity */
  oldVelocity=Memory_Alloc_Array(double, 
                                 lNodeCount * velocityField->fieldComponentCount, 
                                 "artificial nodal velocities");

  store_current_velocity(velocityField,oldVelocity);

  artVelocity=
    Memory_Alloc_Array(double,lNodeCount*velocityField->fieldComponentCount,
                       "artificial nodal velocities");
  compute_velocity_from_displacement(sys->dispField,artVelocity,dt);
                                     
  add_correction(velocityField,artVelocity);

  FeVariable_SyncShadowValues(velocityField);

  /* Solve Energy equation */
  sys->energySolverExecute(sle,context);

  /* Reverse correction and re-sync */
  restore_velocity(velocityField,oldVelocity);
  FeVariable_SyncShadowValues(velocityField);

  Memory_Free(artVelocity);
  Memory_Free(oldVelocity);
}
