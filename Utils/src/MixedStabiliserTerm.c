/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003, Victorian Partnership for Advanced Computing (VPAC) Ltd,
** 110 Victoria Street, Melbourne, 3053, Australia.
**
** Authors:
**	Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**	Patrick D. Sunter, Software Engineer, VPAC. (pds@vpac.org)
**	Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**	Siew-Ching Tan, Software Engineer, VPAC. (siew@vpac.org)
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
** $Id: MixedStabiliserTerm.c 2192 2004-10-15 02:45:38Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <string.h>
#include <assert.h>
#include <mpi.h>

#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include <Underworld/Underworld.h>

#include "types.h"
#include "MixedStabiliserTerm.h"


const Type MixedStabiliserTerm_Type = "MixedStabiliserTerm";


void* _MixedStabiliserTerm_DefaultNew( Name name ) {
   return (void*)_MixedStabiliserTerm_New(
      sizeof(MixedStabiliserTerm), 
      MixedStabiliserTerm_Type,
      _MixedStabiliserTerm_Delete,
      _StiffnessMatrixTerm_Print,
      NULL,
      _MixedStabiliserTerm_DefaultNew,
      MixedStabiliserTerm_Construct,
      MixedStabiliserTerm_Build,
      MixedStabiliserTerm_Initialise,
      MixedStabiliserTerm_Execute,
      MixedStabiliserTerm_Destroy,
      MixedStabiliserTerm_AssembleElement,
      name );
}

/* Creation implementation / Virtual constructor */
MixedStabiliserTerm* _MixedStabiliserTerm_New( 
   SizeT                                               sizeOfSelf,  
   Type                                                type,
   Stg_Class_DeleteFunction*                           _delete,
   Stg_Class_PrintFunction*                            _print,
   Stg_Class_CopyFunction*                             _copy, 
   Stg_Component_DefaultConstructorFunction*           _defaultConstructor,
   Stg_Component_ConstructFunction*                    _construct,
   Stg_Component_BuildFunction*                        _build,
   Stg_Component_InitialiseFunction*                   _initialise,
   Stg_Component_ExecuteFunction*                      _execute,
   Stg_Component_DestroyFunction*                      _destroy,
   StiffnessMatrixTerm_AssembleElementFunction*        _assembleElement,
   Name                                                name )
{
   MixedStabiliserTerm* self;

   /* Allocate memory */
   assert( sizeOfSelf >= sizeof(MixedStabiliserTerm) );
   self = (MixedStabiliserTerm*) _StiffnessMatrixTerm_New( 
      sizeOfSelf, 
      type, 
      _delete, 
      _print, 
      _copy,
      _defaultConstructor,
      _construct,
      _build, 
      _initialise,
      _execute,
      _destroy,
      _assembleElement,
      name );

   _MixedStabiliserTerm_Init(self );

   return self;
}

void _MixedStabiliserTerm_Init( MixedStabiliserTerm* self ) {
   self->picSwarm = NULL;
   self->storeVisc = NULL;
   self->Ni = NULL;
   self->GNx = NULL;
   self->elStiffMat = NULL;
}

void _MixedStabiliserTerm_Delete( void* _self ) {
   MixedStabiliserTerm* self = (MixedStabiliserTerm*)_self;

   if( self->Ni )
      Memory_Free( self->Ni );
   if( self->GNx )
      Memory_Free( self->GNx );
   if( self->elStiffMat )
      Memory_Free( self->elStiffMat );

   _StiffnessMatrixTerm_Delete( self );
}

void MixedStabiliserTerm_Construct( void* _self, Stg_ComponentFactory* cf, void* data ) {
   MixedStabiliserTerm* self = (MixedStabiliserTerm*)_self;

   _MixedStabiliserTerm_Init( self );
   _StiffnessMatrixTerm_Construct( self, cf, data );

   self->picSwarm = Stg_ComponentFactory_ConstructByKey(
      cf, self->name, "picSwarm", IntegrationPointsSwarm,
      True, data );
   self->storeVisc =
     Stg_ComponentFactory_ConstructByKey(cf, self->name, "storeVisc",
                                         StoreVisc, True, data );
}

void MixedStabiliserTerm_Build( void* _self, void* data ) {
   MixedStabiliserTerm* self = (MixedStabiliserTerm*)_self;

   _StiffnessMatrixTerm_Build( self, data );
}

void MixedStabiliserTerm_Initialise( void* _self, void* data ) {
   MixedStabiliserTerm* self = (MixedStabiliserTerm*)_self;
   _StiffnessMatrixTerm_Initialise( self, data );
}

void MixedStabiliserTerm_Execute( void* _self, void* data ) {
   _StiffnessMatrixTerm_Execute( _self, data );
}

void MixedStabiliserTerm_Destroy( void* _self, void* data ) {
   _StiffnessMatrixTerm_Destroy( _self, data );
}

void MixedStabiliserTerm_AssembleElement( void* _self,
					  StiffnessMatrix* stiffMat,
					  int elementIndex,
					  SystemLinearEquations* _sle,
					  FiniteElementContext* ctx,
					  double** elStiffMat )
{
   MixedStabiliserTerm* self = (MixedStabiliserTerm*)_self;
   Stokes_SLE* sle = Stg_DCheckType( _sle, Stokes_SLE );
   StiffnessMatrix* gMatrix = sle->gStiffMat;
   FeVariable* colVar  = gMatrix->columnVariable;
   FeMesh* mesh = colVar->feMesh;
   int nDims = Mesh_GetDimSize( mesh );
   double *xi, weight, jacDet, weightJacDet, cellArea = 0.0;
   double *Ni, **GNx, viscFac;
   int nParticles, nElNodes, cellIndex;
   ElementType* elementType;
   IntegrationPointsSwarm* swarm;
   IntegrationPoint* integrationPoint;
   double geometric_factor;
   double** localElStiffMat;
   int ii, jj, kk;
   double sumVisc = 0.0, visc;

   /* Get the current element's type and cache the number
      of nodes. */
   elementType = FeMesh_GetElementType( mesh, elementIndex );
   nElNodes = elementType->nodeCount;

   /* If we don't already have arrays, allocate them now. */
   if( !self->Ni && !self->GNx && !self->elStiffMat ) {
      self->Ni = Memory_Alloc_Array( double, nElNodes, MixedStabiliserTerm_Type );
      self->GNx = Memory_Alloc_2DArray( double, nDims, nElNodes, MixedStabiliserTerm_Type );
      self->elStiffMat = Memory_Alloc_2DArray( double, nElNodes, nElNodes, MixedStabiliserTerm_Type );
   }

   /* Cache array pointers stored on the mesh. */
   Ni = self->Ni;
   GNx = self->GNx;
   localElStiffMat = self->elStiffMat;

   /* Zero the local element stiffness matrix. */
   for( ii = 0; ii < nElNodes; ii++ )
      memset( localElStiffMat[ii], 0, nElNodes * sizeof(double) );

   /* Assemble the mass matrix portion. */
   swarm = self->integrationSwarm;
   cellIndex = CellLayout_MapElementIdToCellId( swarm->cellLayout,
                                                elementIndex );
   nParticles = swarm->cellParticleCountTbl[cellIndex];
   for( ii = 0; ii < nParticles; ii++ ) {

      /* Cache information from the current integration point. */
      integrationPoint = (IntegrationPoint*)Swarm_ParticleInCellAt(
         swarm, cellIndex, ii );
      xi = integrationPoint->xi;
      weight = integrationPoint->weight;
      ElementType_EvaluateShapeFunctionsAt( elementType, xi, Ni );
      ElementType_ShapeFunctionsGlobalDerivs(
         elementType, mesh, elementIndex, xi, nDims, &jacDet, GNx );
      weightJacDet = weight * jacDet;

      /* Loop over element nodes. */
      for( jj = 0 ; jj < nElNodes; jj++ ) {
         for ( kk = 0 ; kk < nElNodes ; kk++ ) {
            localElStiffMat[jj][kk] += weightJacDet * Ni[jj] * Ni[kk];
         }
      }
   }

   /* Calculate the cell's area and viscosity. */
   swarm = self->picSwarm;
   cellIndex = CellLayout_MapElementIdToCellId( swarm->cellLayout, elementIndex );
   nParticles = swarm->cellParticleCountTbl[cellIndex];
   for( ii = 0; ii < nParticles; ii++ ) {
     StoreVisc_ParticleExt*            particleExt;
     MaterialPointsSwarm*    mSwarm;
     MaterialPoint*          materialparticle;

      /* Cache information from the current integration point. */
      integrationPoint = (IntegrationPoint*)Swarm_ParticleInCellAt(
         swarm, cellIndex, ii );
      ElementType_ShapeFunctionsGlobalDerivs(
         elementType, mesh, elementIndex, integrationPoint->xi,
         nDims, &jacDet, GNx );

      /* Add this particle's value to the area. */
      cellArea += integrationPoint->weight * jacDet;

      materialparticle =
        OneToOneMapper_GetMaterialPoint( swarm->mapper,
                                         integrationPoint, &mSwarm );
      particleExt=
        ExtensionManager_Get( mSwarm->particleExtensionMgr,
                              materialparticle,
                              self->storeVisc->particleExtHandle );

      visc=particleExt->effVisc;
      sumVisc += visc*integrationPoint->weight * jacDet;
   }

   /* Normalize the viscosity factor by dividing by cell area. */
   viscFac = cellArea / sumVisc;

   /* Adjust the calculated mass matrix by the 'special operator'. The
      nElNodes term comes from an averaging operator applied twice. */

   for( ii = 0; ii < nElNodes; ii++ )
      for( jj = 0; jj < nElNodes; jj++ )
        localElStiffMat[ii][jj] -= cellArea/(nElNodes*nElNodes);

   /* Apply the viscosity factor and negate the calculated value. */
   for( ii = 0; ii < nElNodes; ii++ ) {
      for( jj = 0; jj < nElNodes; jj++ )
         elStiffMat[ii][jj] += -localElStiffMat[ii][jj] * viscFac;
   }
}
