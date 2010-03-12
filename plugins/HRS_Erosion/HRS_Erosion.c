/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
** An implementation of the erosion rule in 
**
** Hilley, G. E., M. R. Strecker, and V. A. Ramos (2004), Growth and
** erosion of fold-and- thrust belts with an application to the
** Aconcagua fold-and-thrust belt, Argentina, J. Geophys. Res., 109,
** B01410, doi:10.1029/2002JB002282.
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
#include <float.h>

#include "types.h"
#include "Context.h"
#include "HRS_Erosion.h"


const Type		Underworld_HRS_Erosion_Type = "HRS_Erosion";
ExtensionInfo_Index	Underworld_HRS_Erosion_ContextHandle;

static double a_mean_old=0;
static double next_t_erosion=0;
static double t_erosion=0;
static int first_erosion_flag = 1;

void Underworld_HRS_Erosion_Execute( TimeIntegrand* crdAdvector,
                                        Underworld_HRS_Erosion_Context* spCtx)
{
  double				dt;
  FeVariable *velocity;
  double K = spCtx->K;
  Grid* grid;
  Mesh* mesh;
  unsigned nDims;
  Node_LocalIndex n_i, n_left;
  unsigned nNodes;
  MPI_Comm comm;
  double DT = spCtx->DT;
  double first_t_erosion = spCtx->first_t_erosion;
  double vT = spCtx->vT;
  IJK ijk_right;
  IJK ijk_left;
  Bool on_top, on_right, found;
  double x_right, y_right, x_left, y_left, temp, base_height, threshold;

  assert( spCtx );
  
  dt = spCtx->ctx->dt;
  mesh=spCtx->mesh;
  velocity=spCtx->v;
  t_erosion = t_erosion + dt;
  comm = Comm_GetMPIComm( Mesh_GetCommTopology( mesh, MT_VERTEX ) );
		
  if(next_t_erosion < first_t_erosion) {
    next_t_erosion = first_t_erosion;
  }
		
  if(t_erosion < next_t_erosion)
    return;

  next_t_erosion = next_t_erosion + DT;

  if(first_erosion_flag) {
    DT = t_erosion;
    first_erosion_flag = 0;
  }

  nDims = Mesh_GetDimSize( mesh );
  grid =
    *(Grid**)ExtensionManager_Get( mesh->info, mesh, 
                                   ExtensionManager_GetHandle( mesh->info, 
                                                               "vertexGrid" ) );
    
  if( nDims != 2 )
    abort();

  nNodes = FeMesh_GetNodeLocalSize( mesh);

  /* Get the coordinates for the top right corner */

  RegularMeshUtils_Node_1DTo3D
    ( mesh, Mesh_DomainToGlobal( mesh, MT_VERTEX, nNodes-1 ), ijk_right );
  on_top=(ijk_right[1]==grid->sizes[1]-1);
  on_right=(ijk_right[0]==grid->sizes[0]-1);
  if(on_top && on_right)
    {
      x_right=mesh->verts[nNodes-1][0];
      y_right=mesh->verts[nNodes-1][1];
    }
  else
    {
      x_right=DBL_MAX;
      y_right=DBL_MAX;
    }
  MPI_Allreduce( &x_right, &temp, 1, MPI_DOUBLE, MPI_MIN, comm );
  x_right=temp;
  MPI_Allreduce( &y_right, &temp, 1, MPI_DOUBLE, MPI_MIN, comm );
  y_right=temp;

  /* Get the coordinates for where the left side rises above the base
     height */
  base_height=DBL_MAX;

  if(on_top)
    {
      int i_right, i_left, i, j;
      i_right=ijk_right[0];
      j=ijk_right[1];
      
      RegularMeshUtils_Node_1DTo3D
        ( mesh, Mesh_DomainToGlobal( mesh, MT_VERTEX, 0 ), ijk_left );

      i_left=ijk_left[0];
      ijk_left[1]=j;

      if(!Mesh_GlobalToDomain
         (mesh,MT_VERTEX,
          RegularMeshUtils_Node_3DTo1D(mesh,ijk_left),&n_left)) 
        {
          printf("Can not map to local domain %d %d\n",
                 ijk_left[0],ijk_left[1]);
          abort();
        }
      base_height=mesh->verts[n_left][1];
      found=False;
      for(i=i_left+1; i<=i_right; ++i)
        {
          ijk_left[0]=i;
          if(!Mesh_GlobalToDomain
             (mesh,MT_VERTEX,
              RegularMeshUtils_Node_3DTo1D(mesh,ijk_left),&n_left)) 
            {
              printf("Can not map to local domain %d %d\n",
                     ijk_left[0],ijk_left[1]);
              abort();
            }
          if(mesh->verts[n_left][1]>base_height+threshold)
            {
              found=True;
              x_left=mesh->verts[n_left][0];
              y_left=mesh->verts[n_left][1];
              break;
            }
        }
      if(!found)
        {
          y_left=DBL_MAX;
          x_left=DBL_MAX;
        }
    }
  MPI_Allreduce( &x_left, &temp, 1, MPI_DOUBLE, MPI_MIN, comm );
  x_left=temp;
  MPI_Allreduce( &y_left, &temp, 1, MPI_DOUBLE, MPI_MIN, comm );
  y_left=temp;

  /* Leave if no one found the rise */

  if(x_left==DBL_MAX)
    {
      printf("Could not locate the edge of the wedge.  Not eroding.\n");
      return;
    }

  /* Start the erosion calculation */

  double W = x_right-x_left;
  double a_mean = atan((y_right-base_height)/W);
				
  /* Insert equation here and calculate new slope: */

  const double m=0.4; 
  const double n=1; 
  const double ka=4;
  const double h=1.4; 
  double S=tan(a_mean_old);

  /* alpha_calculated=(alpha1+atan((2.*vT./W1.^2 - 2*K*ka*W1^(h*m-1)*S^n/(h*m+1))*dt)) */
  double alpha_calculated=
    (a_mean_old+atan((2*vT/pow(W,2)
                      - (2*K*ka*pow(W,(h*m-1))*pow(S,n))/(h*m+1))*DT));

  /* At the end of this calculation, we need a variable called alpha_calculated (in radians): */
				
  double slope_calculated = tan(alpha_calculated);
  double ypred;
				
  printf("a_mean_old | a_mean | alpha_calculated | width | DT\n%11.9g %11.9g %11.9g %11.9g %11.9g\n",a_mean_old,a_mean,alpha_calculated, W, DT);
				
  if (a_mean>alpha_calculated)
    {
      a_mean_old=alpha_calculated;
      
      for( n_i = 0; n_i < nNodes; n_i++ ) 
        {
          IJK ijk;
          RegularMeshUtils_Node_1DTo3D
            ( mesh, Mesh_DomainToGlobal( mesh, MT_VERTEX, n_i ), ijk );
					
          if(ijk[1]==grid->sizes[1]-1
             && ijk[0]!=0 && ijk[0]!=grid->sizes[0]-1) 
            {
              double y_old,delta_v, v[3];
              double x; 
              
              x=mesh->verts[n_i][0];
              ypred = slope_calculated * (x-x_left) + y_left;
              y_old = mesh->verts[n_i][1];
              if(ypred < y_old && x>x_left) 
                {
                  delta_v = (ypred-y_old)/dt;
                }
              else 
                {
                  delta_v = 0;
                }
						
              FeVariable_GetValueAtNode(velocity,n_i,v);

              v[1]+=delta_v;
              FeVariable_SetValueAtNode(velocity,n_i,v);
            }
        }
    }
  else
    {	
      a_mean_old=a_mean;
    }
  FeVariable_SyncShadowValues(velocity);
}

Index Underworld_HRS_Erosion_Register( PluginsManager* pluginsMgr ) {
	return PluginsManager_Submit( pluginsMgr, 
				      Underworld_HRS_Erosion_Type, 
				      "0", 
				      _Underworld_HRS_Erosion_DefaultNew );
}

void* _Underworld_HRS_Erosion_DefaultNew( Name name ) {

	/* Variables set in this function */
	SizeT                                              _sizeOfSelf = sizeof(Codelet);
	Type                                                      type = Underworld_HRS_Erosion_Type;
	Stg_Class_DeleteFunction*                              _delete = _Codelet_Delete;
	Stg_Class_PrintFunction*                                _print = _Codelet_Print;
	Stg_Class_CopyFunction*                                  _copy = _Codelet_Copy;
	Stg_Component_DefaultConstructorFunction*  _defaultConstructor = _Underworld_HRS_Erosion_DefaultNew;
	Stg_Component_ConstructFunction*                    _construct = _Underworld_HRS_Erosion_AssignFromXML;
	Stg_Component_BuildFunction*                            _build = _Underworld_HRS_Erosion_Build;
	Stg_Component_InitialiseFunction*                  _initialise = _Codelet_Initialise;
	Stg_Component_ExecuteFunction*                        _execute = _Codelet_Execute;
	Stg_Component_DestroyFunction*                        _destroy = _Underworld_HRS_Erosion_Destroy;

	/* Variables that are set to ZERO are variables that will be
           set either by the current _New function or another parent
           _New function further up the hierachy */
	AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

	return _Codelet_New(  CODELET_PASSARGS   );
}

void _Underworld_HRS_Erosion_AssignFromXML( void* component,
                                               Stg_ComponentFactory* cf,
                                               void* data ) {
  Codelet* sp = (Codelet*)component;
  UnderworldContext*			uwCtx;
  Underworld_HRS_Erosion_Context*	spCtx;
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
  Underworld_HRS_Erosion_ContextHandle =
    ExtensionManager_Add( uwCtx->extensionMgr, 
                          Underworld_HRS_Erosion_Type, 
                          sizeof(Underworld_HRS_Erosion_Context) );

  spCtx = ExtensionManager_Get( uwCtx->extensionMgr, uwCtx,
                                Underworld_HRS_Erosion_ContextHandle );
  memset( spCtx, 0, sizeof(Underworld_HRS_Erosion_Context) );
  spCtx->ctx = (AbstractContext*)uwCtx;

  /* Get the time integrator. */
  spCtx->timeIntegrator =
    Stg_ComponentFactory_ConstructByName( cf, (Name)"timeIntegrator",
                                          TimeIntegrator, True, data  );

  /* Get the dictionary. */
  spDict = Dictionary_GetDictionary( uwCtx->dictionary, "HRS_Erosion" );
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
  spCtx->K = Dictionary_GetDouble( spDict, "K" );

  spCtx->DT = Dictionary_GetDouble(spDict, "dt_erosion");
  spCtx->first_t_erosion = Dictionary_GetDouble(spDict, "first_t_erosion");
  spCtx->vT = Dictionary_GetDouble(spDict, "vT");
}

void _Underworld_HRS_Erosion_Build( void* codelet, void* data ) {
	Codelet* sp= (Codelet*)codelet;
	UnderworldContext* UnderworldCtx = (UnderworldContext*)sp->context;
	Underworld_HRS_Erosion_Context*	spCtx;

	assert( codelet );
	assert( UnderworldCtx );

	/* Get the context. */
	spCtx = ExtensionManager_Get( UnderworldCtx->extensionMgr, UnderworldCtx, Underworld_HRS_Erosion_ContextHandle );

	if( !spCtx->mesh )
		return;

	/* Append to the list of time integratee finish routines.  It
           should come last, because EulerDeform will reset the
           values. */
	TimeIntegrator_AppendSetupEP( spCtx->timeIntegrator, 
					"Underworld_HRS_Erosion_Execute", 
					Underworld_HRS_Erosion_Execute, 
					"HRS_Erosion", 
					spCtx );
}

void _Underworld_HRS_Erosion_Destroy( void* codelet, void* data ) {
	UnderworldContext*	UnderworldCtx = (UnderworldContext*)data;

	assert( codelet );
	assert( UnderworldCtx );

	/* Clear the lot. */
	/* TODO */
}

