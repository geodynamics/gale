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

void EulerDeform_Build(void* component, void* data)
{
  Codelet* ed=(Codelet*)component;
  UnderworldContext* uwCtx=(UnderworldContext*)ed->context;
  EulerDeform_Context* edCtx;
  Variable* crdVar;
  TimeIntegrand* crdAdvector;
  Stg_Component* tiData[2];
  unsigned sys_i;
  Dictionary_Entry_Value* edDict;
  Dictionary_Entry_Value* sysLst;

  assert(component);
  assert(uwCtx);

  edCtx=(EulerDeform_Context*)
    ExtensionManager_Get(uwCtx->extensionMgr,uwCtx,EulerDeform_ContextHandle);

  /* Get the dictionary. */
  edDict=Dictionary_Get(uwCtx->dictionary,"EulerDeform");
  if(!edDict)
    return;

  /* Read system list. */
  sysLst=Dictionary_Entry_Value_GetMember(edDict,"systems");
  if(sysLst)
    {
      unsigned sys_i;

      /* Allocate for systems. */
      edCtx->nSystems=Dictionary_Entry_Value_GetCount(sysLst);
      edCtx->systems=Memory_Alloc_Array(EulerDeform_System,edCtx->nSystems,
                                        "EulerDeform->systems");
      memset(edCtx->systems,0,sizeof(EulerDeform_System)*edCtx->nSystems);

      for(sys_i=0; sys_i<edCtx->nSystems; sys_i++)
        {
          EulerDeform_System* sys=edCtx->systems+sys_i;
          Dictionary* sysDict;
          char* meshName;
          char* remesherName;
          char* name;

          /* Get the dictionary for this system. */
          sysDict=Dictionary_Entry_Value_AsDictionary
            (Dictionary_Entry_Value_GetElement(sysLst,sys_i));
          assert(sysDict);

          /* Read contents. */
          meshName=Dictionary_GetString(sysDict,"mesh");
          sys->mesh=Stg_ComponentFactory_ConstructByName(uwCtx->CF,
                                                         (Name)meshName,
                                                         Mesh,True,data);
          char *p_MeshName=Dictionary_GetString(sysDict,"p-Mesh");
                                                
          if(strcmp(p_MeshName,""))
            sys->p_mesh=Stg_ComponentFactory_ConstructByName(uwCtx->CF,
                                                             (Name)p_MeshName,
                                                             Mesh,True,data);

          name=Dictionary_GetString(sysDict,"displacementField");
                                    
          if(strcmp(name, ""))
            sys->dispField=
              Stg_ComponentFactory_ConstructByName(uwCtx->CF,(Name)name,
                                                   FeVariable,True,data);
          else
            sys->dispField=NULL;

          AdvectionDiffusionSLE* energySLE=(AdvectionDiffusionSLE*)
            Stg_ComponentFactory_ConstructByName(uwCtx->CF,(Name)"EnergyEqn",
                                                 AdvectionDiffusionSLE,
                                                 False,data);
          Journal_Firewall(!(energySLE!=NULL && sys->dispField==NULL),
                           Journal_Register(Error_Type,
                                            (Name)Underworld_EulerDeform_Type),
                           "If you enable thermal evolution, you need to add"
                           "a displacement field to EulerDeform\n");

          if(sys->dispField!=NULL)
            {
              Journal_Firewall
                (energySLE!=NULL, 
                 Journal_Register(Error_Type,
                                  (Name)Underworld_EulerDeform_Type),
                 "The required energy SLE component has not been created "
                 "or placed on the context.\n");

              /* Replace the energy SLE's execute with this one. Save
                 the old value for use later. */
              sys->energySolverExecute = energySLE->_execute;
              energySLE->_execute = EulerDeform_Advection_Correction;
            }

          remesherName=Dictionary_GetString(sysDict,"remesher");
                                            
          if(strcmp(remesherName,""))
            sys->remesher=
              Stg_ComponentFactory_ConstructByName(uwCtx->CF,(Name)remesherName,
                                                   Remesher,True,data);

          sys->interval=Dictionary_GetInt_WithDefault(sysDict,"interval",-1);
                                          
          sys->wrapTop=Dictionary_GetBool_WithDefault(sysDict,"wrapTop",False);
          sys->wrapBottom=Dictionary_GetBool_WithDefault(sysDict,"wrapBottom",
                                                         False);
          sys->wrapLeft=Dictionary_GetBool_WithDefault(sysDict,"wrapLeft",
                                                       False);
          /* This line is currently not working, have to manually set
             the velocity field name.  This should be fixed once this
             plugin has been converted to a component. */

          /*sys->velField = Stg_ComponentFactory_ConstructByName( uwCtx->CF, (Name)velFieldName, FieldVariable, True, data  );*/
          sys->velField=
            Stg_ComponentFactory_ConstructByName(uwCtx->CF,
                                                 (Name)"VelocityField",
                                                 FieldVariable,True,data);

          sys->staticTop=Dictionary_GetBool_WithDefault(sysDict,"staticTop",
                                                        False);
          sys->staticBottom=Dictionary_GetBool_WithDefault(sysDict,
                                                           "staticBottom",
                                                           False);
          sys->staticLeft=Dictionary_GetBool_WithDefault(sysDict,"staticLeft",
                                                         False);
          sys->staticRight=Dictionary_GetBool_WithDefault(sysDict,"staticRight",
                                                          False);
          sys->staticFront=Dictionary_GetBool_WithDefault(sysDict,"staticFront",
                                                          False);
          sys->staticBack=Dictionary_GetBool_WithDefault(sysDict,"staticBack",
                                                         False);

          sys->staticLeftTop=
            Dictionary_GetBool_WithDefault(sysDict,"staticLeftTop",
                                           (sys->staticLeft && sys->staticTop)
                                           ? True : False);
          sys->staticRightTop=
            Dictionary_GetBool_WithDefault(sysDict,"staticRightTop",
                                           (sys->staticRight && sys->staticTop)
                                           ? True : False);
          sys->staticLeftTopFront=
            Dictionary_GetBool_WithDefault(sysDict,"staticLeftTopFront",
                                           (sys->staticLeft && sys->staticTop
                                            && sys->staticFront)
                                           ? True : False);
          sys->staticRightTopFront=
            Dictionary_GetBool_WithDefault(sysDict,"staticRightTopFront",
                                           (sys->staticRight && sys->staticTop
                                            && sys->staticFront)
                                           ? True : False);
          sys->staticLeftTopBack=
            Dictionary_GetBool_WithDefault(sysDict,"staticLeftTopBack",
                                           (sys->staticLeft && sys->staticTop
                                            && sys->staticBack)
                                           ? True : False);
          sys->staticRightTopBack=
            Dictionary_GetBool_WithDefault(sysDict,"staticRightTopBack",
                                           (sys->staticRight && sys->staticTop
                                            && sys->staticBack)
                                           ? True : False);

          sys->staticLeftBottom=
            Dictionary_GetBool_WithDefault(sysDict,"staticLeftBottom",
                                           (sys->staticLeft && sys->staticBottom)
                                           ? True : False);
          sys->staticRightBottom=
            Dictionary_GetBool_WithDefault(sysDict,"staticRightBottom",
                                           (sys->staticRight && sys->staticBottom)
                                           ? True : False);
          sys->staticLeftBottomFront=
            Dictionary_GetBool_WithDefault(sysDict,"staticLeftBottomFront",
                                           (sys->staticLeft && sys->staticBottom
                                            && sys->staticFront)
                                           ? True : False);
          sys->staticRightBottomFront=
            Dictionary_GetBool_WithDefault(sysDict,"staticRightBottomFront",
                                           (sys->staticRight && sys->staticBottom
                                            && sys->staticFront)
                                           ? True : False);
          sys->staticLeftBottomBack=
            Dictionary_GetBool_WithDefault(sysDict,"staticLeftBottomBack",
                                           (sys->staticLeft && sys->staticBottom
                                            && sys->staticBack)
                                           ? True : False);
          sys->staticRightBottomBack=
            Dictionary_GetBool_WithDefault(sysDict,"staticRightBottomBack",
                                           (sys->staticRight && sys->staticBottom
                                            && sys->staticBack)
                                           ? True : False);

          sys->staticLeftFront=
            Dictionary_GetBool_WithDefault(sysDict,"staticLeftFront",
                                           (sys->staticLeft && sys->staticFront)
                                           ? True : False);
          sys->staticRightFront=
            Dictionary_GetBool_WithDefault(sysDict,"staticRightFront",
                                           (sys->staticRight && sys->staticFront)
                                           ? True : False);
          sys->staticLeftBack=
            Dictionary_GetBool_WithDefault(sysDict,"staticLeftBack",
                                           (sys->staticLeft && sys->staticBack)
                                           ? True : False);
          sys->staticRightBack=
            Dictionary_GetBool_WithDefault(sysDict,"staticRightBack",
                                           (sys->staticRight && sys->staticBack)
                                           ? True : False);

          sys->staticTopFront=
            Dictionary_GetBool_WithDefault(sysDict,"staticTopFront",
                                           (sys->staticTop && sys->staticFront)
                                           ? True : False);
          sys->staticBottomFront=
            Dictionary_GetBool_WithDefault(sysDict,"staticBottomFront",
                                           (sys->staticBottom && sys->staticFront)
                                           ? True : False);
          sys->staticTopBack=
            Dictionary_GetBool_WithDefault(sysDict,"staticTopBack",
                                           (sys->staticTop && sys->staticBack)
                                           ? True : False);
          sys->staticBottomBack=
            Dictionary_GetBool_WithDefault(sysDict,"staticBottomBack",
                                           (sys->staticBottom && sys->staticBack)
                                           ? True : False);

          sys->floatLeftTop=
            Dictionary_GetBool_WithDefault(sysDict,"floatLeftTop",False);
          sys->floatRightTop=
            Dictionary_GetBool_WithDefault(sysDict,"floatRightTop",False);

          sys->staticSides = 
            (sys->staticLeft
             || sys->staticRight
             || sys->staticTop
             || sys->staticBottom
             || sys->staticFront
             || sys->staticBack
             || sys->staticLeftTop
             || sys->staticRightTop
             || sys->staticLeftTopFront
             || sys->staticRightTopFront
             || sys->staticLeftTopBack
             || sys->staticRightTopBack
             || sys->staticLeftBottom
             || sys->staticRightBottom
             || sys->staticLeftBottomFront
             || sys->staticRightBottomFront
             || sys->staticLeftBottomBack
             || sys->staticRightBottomBack
             || sys->staticLeftFront
             || sys->staticRightFront
             || sys->staticLeftBack
             || sys->staticRightBack
             || sys->staticTopFront
             || sys->staticBottomFront
             || sys->staticTopBack
             || sys->staticBottomBack)
            ? True : False;


          /* TODO: fix to allow equations */
          if(sys->staticRight && sys->wrapTop
             && !sys->staticRightTop)
            sys->static_right_coord =
              Dictionary_GetDouble( uwCtx->dictionary, "maxX");

          if(sys->staticLeft && sys->wrapTop
             && !sys->staticLeftTop)
            sys->static_left_coord =
              Dictionary_GetDouble( uwCtx->dictionary, "minX");

          if(sys->staticFront && sys->wrapTop
             && !sys->staticTopFront)
            sys->static_front_coord =
              Dictionary_GetDouble( uwCtx->dictionary, "minZ");

          if(sys->staticBack && sys->wrapTop
             && !sys->staticTopBack)
            sys->static_back_coord =
              Dictionary_GetDouble( uwCtx->dictionary, "maxZ");
        }
    }

  for(sys_i=0; sys_i<edCtx->nSystems; sys_i++)
    {
      EulerDeform_System* sys=edCtx->systems+sys_i;

      /* Create a time integrand for the mesh's coordinates. */
      crdVar=
        EulerDeform_RegisterLocalNodeCoordsAsVariables(sys,
                                                       uwCtx->variable_Register,
                                                       NULL);
      Stg_Component_Build(crdVar,data,False);

      tiData[0] = (Stg_Component*)sys->velField;
      tiData[1] = (Stg_Component*)&sys->mesh->verts;
      crdAdvector = TimeIntegrand_New("EulerDeform_Velocity",
                                      (DomainContext*)uwCtx,
                                      edCtx->timeIntegrator, crdVar, 2,
                                      tiData, True
                                      /* Presume we need to allow
                                         fallback on edges of
                                         stretching mesh -
                                         PatrickSunter, 7 June 2006 */ );
      crdAdvector->_calculateTimeDeriv=EulerDeform_TimeDeriv;

      /* Add to live component register... */
      LiveComponentRegister_Add(uwCtx->CF->LCRegister,
                                (Stg_Component*)crdAdvector);
      Stg_Component_Build(crdAdvector,data,False);
    }

  if(edCtx->nSystems>0)
    {
      /* Insert the sync step. */
      TimeIntegrator_PrependSetupEP(edCtx->timeIntegrator,
                                    "EulerDeform_IntegrationSetup",
                                    (void*)EulerDeform_IntegrationSetup,
                                    "EulerDeform",edCtx);
    }

  /* Insert the remesh step. Ideally this would look for the
     surface process plugin's time integrator finish routine and
     ensure we enter the remesh step after that one but before the
     particle updating routines. */
  TimeIntegrator_PrependFinishEP(edCtx->timeIntegrator,"EulerDeform_Execute",
                                 (void*)EulerDeform_Remesh,"EulerDeform",edCtx);
}

