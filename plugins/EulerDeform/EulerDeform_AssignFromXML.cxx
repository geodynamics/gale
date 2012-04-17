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

void EulerDeform_AssignFromXML(void* component,
                               Stg_ComponentFactory* cf,
                               void* data)
{
  Codelet* ed=(Codelet*)component;
  EulerDeform_Context* edCtx;

  assert(component);
  assert(cf);

  Journal_DPrintf(Underworld_Debug, "In: %s( void* )\n", __func__);

  UnderworldContext* uwCtx=(UnderworldContext*)
    Stg_ComponentFactory_ConstructByName(cf,(Name)"context",UnderworldContext,
                                         True,data);
  ed->context=(AbstractContext*)uwCtx;

  /* Create new context. */
  EulerDeform_ContextHandle=
    ExtensionManager_Add(uwCtx->extensionMgr,(Name)Underworld_EulerDeform_Type,
                         sizeof(EulerDeform_Context));
  edCtx=(EulerDeform_Context*)ExtensionManager_Get(uwCtx->extensionMgr,
                                                   uwCtx,
                                                   EulerDeform_ContextHandle);
  memset(edCtx,0,sizeof(EulerDeform_Context));
  edCtx->ctx=(AbstractContext*)uwCtx;

  /* Get the time integrator. */
  edCtx->timeIntegrator=
    Stg_ComponentFactory_ConstructByName(cf,(Name)"timeIntegrator",
                                         TimeIntegrator,True,data);
}


