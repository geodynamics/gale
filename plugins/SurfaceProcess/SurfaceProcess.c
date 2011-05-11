/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street, Melbourne, 3053, Australia.
** Copyright (C) 2008, 2010 California Institute of Technology
**
** Authors:
**	Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**	Patrick D. Sunter, Software Engineer, VPAC. (pds@vpac.org)
**	David May, PhD Student Monash University, VPAC. (david.may@sci.maths.monash.edu.au)
**	Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**	Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
**	Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
**      Walter Landry, CIG
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
#include <math.h>

#include "types.h"
#include "Context.h"
#include "SurfaceProcess.h"


const Type		Underworld_SurfaceProcess_Type = "SurfaceProcess";
ExtensionInfo_Index	Underworld_SurfaceProcess_ContextHandle;

void Underworld_SurfaceProcess_Execute( TimeIntegrand* crdAdvector,
                                        Underworld_SurfaceProcess_Context* spCtx)
{
  double				dt;
  
  /*
  ** SURFACE PROCESS CODE GOES HERE, SHOULD MODIFY THE VELOCITIES ONLY.
  */
  FeVariable *velocity;
  double K = spCtx->K;
  Grid* grid;
  Mesh* mesh;
  unsigned nDims;
  Node_LocalIndex n_i;
  unsigned nNodes;

  assert( spCtx );
  
  /* Extract information from contexts. */
  dt = spCtx->ctx->dt;

  mesh=spCtx->mesh;
  velocity=spCtx->v;
    
  nDims = Mesh_GetDimSize( mesh );
  grid =
    *(Grid**)ExtensionManager_Get( mesh->info, mesh, 
                                   ExtensionManager_GetHandle( mesh->info, 
                                                               "vertexGrid" ) );
    
  if( nDims != 2 )
    abort();

  nNodes = FeMesh_GetNodeLocalSize( mesh);
      
  for( n_i = 0; n_i < nNodes; n_i++ )
    {
      IJK ijk;
      RegularMeshUtils_Node_1DTo3D
        ( mesh, Mesh_DomainToGlobal( mesh, MT_VERTEX, n_i ), ijk );
                  
      if(ijk[1]==grid->sizes[1]-1
         && ijk[0]!=0 && ijk[0]!=grid->sizes[0]-1)
        {
          IJK ijk_minus, ijk_plus;
          Node_LocalIndex n_plus, n_minus;
          double y, y_plus, y_minus, new_y, average_y, dx, delta_v, v[3];

          Vec_Set2D(ijk_minus,ijk);
          ijk_minus[0]-=1;
          if(!Mesh_GlobalToDomain
             (mesh,MT_VERTEX,
              RegularMeshUtils_Node_3DTo1D(mesh,ijk_minus),&n_minus))
            {
              printf("Can not map to local domain %d %d %d\n",
                     ijk_minus[0],ijk_minus[1],n_i);
              abort();
            }
          Vec_Set2D(ijk_plus,ijk);
          ijk_plus[0]+=1;
          if(!Mesh_GlobalToDomain
             (mesh,MT_VERTEX,
              RegularMeshUtils_Node_3DTo1D(mesh,ijk_plus),&n_plus))
            {
              printf("Can not map to local domain %d %d %d\n",
                     ijk_plus[0],ijk_plus[1],n_i);
              abort();
            }

          y_plus=mesh->verts[n_plus][1];
          y_minus=mesh->verts[n_minus][1];
          y=mesh->verts[n_i][1];

          dx=mesh->verts[n_i][0]-mesh->verts[n_minus][0];

          /* If the diffusion overcorrects, then we set the
             new value to the average */
          new_y=y+dt*K*(y_plus+y_minus-2*y)/(dx*dx);
          average_y=(y_plus+y_minus)/2;
          if((average_y>y && average_y>new_y)
             || (average_y<y && average_y<new_y))
            {
              delta_v=(new_y-y)/dt;
            }
          else
            {
              delta_v=(average_y-y)/dt;
            }
          FeVariable_GetValueAtNode(velocity,n_i,v);
          v[1]+=delta_v;
          FeVariable_SetValueAtNode(velocity,n_i,v);
        }
    }
  /*
  ** END SURFACE PROCESS CODE.
  */
  FeVariable_SyncShadowValues(velocity);
}

Index Underworld_SurfaceProcess_Register( PluginsManager* pluginsMgr ) {
	return PluginsManager_Submit( pluginsMgr, 
				      Underworld_SurfaceProcess_Type, 
				      "0", 
				      _Underworld_SurfaceProcess_DefaultNew );
}

void* _Underworld_SurfaceProcess_DefaultNew( Name name ) {

	/* Variables set in this function */
	SizeT                                              _sizeOfSelf = sizeof(Codelet);
	Type                                                      type = Underworld_SurfaceProcess_Type;
	Stg_Class_DeleteFunction*                              _delete = _Codelet_Delete;
	Stg_Class_PrintFunction*                                _print = _Codelet_Print;
	Stg_Class_CopyFunction*                                  _copy = _Codelet_Copy;
	Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _Underworld_SurfaceProcess_DefaultNew;
	Stg_Component_ConstructFunction*                    _construct = _Underworld_SurfaceProcess_AssignFromXML;
	Stg_Component_BuildFunction*                            _build = _Underworld_SurfaceProcess_Build;
	Stg_Component_InitialiseFunction*                  _initialise = _Codelet_Initialise;
	Stg_Component_ExecuteFunction*                        _execute = _Codelet_Execute;
	Stg_Component_DestroyFunction*                        _destroy = _Underworld_SurfaceProcess_Destroy;

	/* Variables that are set to ZERO are variables that will be
           set either by the current _New function or another parent
           _New function further up the hierachy */
	AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

	return _Codelet_New(  CODELET_PASSARGS   );
}

void _Underworld_SurfaceProcess_AssignFromXML( void* component,
                                               Stg_ComponentFactory* cf,
                                               void* data ) {
  Codelet* sp = (Codelet*)component;
  UnderworldContext*			uwCtx;
  Underworld_SurfaceProcess_Context*	spCtx;
  Dictionary*			spDict;
  char*				meshName;
  char*                         velocityName;

  assert( component );
  assert( cf );

  Journal_DPrintf( Underworld_Debug, "In: %s( void* )\n", __func__ );

  /* Retrieve context. */
  uwCtx =
    (UnderworldContext*)Stg_ComponentFactory_ConstructByName(cf,"context",
                                                             UnderworldContext,
                                                             True, data );
  sp->context=(AbstractContext* )uwCtx;
  /* Create new context. */
  Underworld_SurfaceProcess_ContextHandle =
    ExtensionManager_Add( uwCtx->extensionMgr, 
                          Underworld_SurfaceProcess_Type, 
                          sizeof(Underworld_SurfaceProcess_Context) );

  spCtx = (Underworld_SurfaceProcess_Context*)ExtensionManager_Get( uwCtx->extensionMgr, uwCtx,
                                Underworld_SurfaceProcess_ContextHandle );
  memset( spCtx, 0, sizeof(Underworld_SurfaceProcess_Context) );
  spCtx->ctx = (AbstractContext*)uwCtx;

  /* Get the time integrator. */
  spCtx->timeIntegrator =
    Stg_ComponentFactory_ConstructByName( cf, (Name)"timeIntegrator",
                                          TimeIntegrator, True, data  );

  /* Get the dictionary. */
  spDict = Dictionary_GetDictionary( uwCtx->dictionary, "SurfaceProcess" );
  if( !spDict )
    return;

  /* Read in the variables. */
  meshName = Dictionary_GetString( spDict, "mesh" );
  assert( meshName && strcmp( meshName, "" ) );
  spCtx->mesh = Stg_ComponentFactory_ConstructByName( cf, meshName, Mesh, True,
                                                      NULL );
  velocityName = Dictionary_GetString( spDict, "VelocityField" );
  assert( velocityName && strcmp( velocityName, "" ) );
  spCtx->v = Stg_ComponentFactory_ConstructByName( cf, velocityName, FeVariable,
                                                   True, NULL );
  assert( spCtx->v);
  spCtx->K = Dictionary_GetDouble( spDict, "diffusionCoefficient" );
}

void _Underworld_SurfaceProcess_Build( void* codelet, void* data ) {
	Codelet* sp= (Codelet*)codelet;
	UnderworldContext* UnderworldCtx = (UnderworldContext*)sp->context;
	Underworld_SurfaceProcess_Context*	spCtx;

	assert( codelet );
	assert( UnderworldCtx );

	/* Get the context. */
	spCtx = (Underworld_SurfaceProcess_Context*)ExtensionManager_Get( UnderworldCtx->extensionMgr, UnderworldCtx, Underworld_SurfaceProcess_ContextHandle );

	if( !spCtx->mesh )
		return;

	/* Append to the list of time integratee finish routines.  It
           should come last, because EulerDeform will reset the
           values. */
	TimeIntegrator_AppendSetupEP( spCtx->timeIntegrator, 
					"Underworld_SurfaceProcess_Execute", 
                                      (void*)Underworld_SurfaceProcess_Execute, 
					"SurfaceProcess", 
					spCtx );
}

void _Underworld_SurfaceProcess_Destroy( void* codelet, void* data ) {
	UnderworldContext*	UnderworldCtx = (UnderworldContext*)data;

	assert( codelet );
	assert( UnderworldCtx );

	/* Clear the lot. */
	/* TODO */
}

