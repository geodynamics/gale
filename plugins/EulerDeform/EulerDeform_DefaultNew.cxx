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

void* EulerDeform_DefaultNew(Name name)
{
  /* Variables set in this function */
  SizeT _sizeOfSelf = sizeof(Codelet);
  Type type = Underworld_EulerDeform_Type;
  Stg_Class_DeleteFunction* _delete = _Codelet_Delete;
  Stg_Class_PrintFunction* _print = _Codelet_Print;
  Stg_Class_CopyFunction* _copy = _Codelet_Copy;
  Stg_Component_DefaultConstructorFunction* _defaultConstructor =
    EulerDeform_DefaultNew;
  Stg_Component_ConstructFunction* _construct =
    EulerDeform_AssignFromXML;
  Stg_Component_BuildFunction* _build = EulerDeform_Build;
  Stg_Component_InitialiseFunction* _initialise = _Codelet_Initialise;
  Stg_Component_ExecuteFunction* _execute = _Codelet_Execute;
  Stg_Component_DestroyFunction* _destroy = EulerDeform_Destroy;

  /* Variables that are set to ZERO are variables that will be set
     either by the current _New function or another parent _New
     function further up the hierachy */
  AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;
  
  return _Codelet_New(  CODELET_PASSARGS   );
}


