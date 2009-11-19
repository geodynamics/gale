/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
** Copyright (c) 2005, Monash Cluster Computing
** All rights reserved.
** Redistribution and use in source and binary forms, with or without modification,
** are permitted provided that the following conditions are met:
**
**       * Redistributions of source code must retain the above copyright notice,
**          this list of conditions and the following disclaimer.
**       * Redistributions in binary form must reproduce the above copyright
**       notice, this list of conditions and the following disclaimer in the
**       documentation and/or other materials provided with the distribution.
**       * Neither the name of the Monash University nor the names of its contributors
**       may be used to endorse or promote products derived from this software
**       without specific prior written permission.
**
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
** THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
** PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS
** BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
** CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
** SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
** HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
** LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
** OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
**
**
** Contact:
*%    Louis Moresi - Louis.Moresi@sci.monash.edu.au
*%
** Contributors:
*+    Robert Turnbull
*+    Vincent Lemiale
*+    Louis Moresi
*+    David May
*+    David Stegman
*+    Mirko Velic
*+    Patrick Sunter
*+    Julian Giordani
*+
** $Id: Byerlee.c 733 2008-05-14 04:55:40Z LukeHodkinson $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>

#include "types.h"
#include "RheologyClass.h"
#include "StrainWeakening.h"
#include "YieldRheology.h"
#include "VonMises.h"
#include "Byerlee.h"
#include "ConstitutiveMatrix.h"

#include <assert.h>

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type Byerlee_Type = "Byerlee";


/* Public Constructor */
Byerlee* Byerlee_New(
      Name                  name,
      AbstractContext*      context,
      StrainWeakening*      strainWeakening, 
      MaterialPointsSwarm*  materialPointsSwarm, 
      double                minVisc, 
      FeVariable*           strainRateField,
      SwarmVariable*        swarmStrainRate,
      double                cohesion,
      double                cohesionAfterSoftening,
      Bool                  strainRateSoftening,
      FeMesh*               mesh,
      double                depthCoefficient )
{
   Byerlee* self = (Byerlee*) _Byerlee_DefaultNew( name );

   _Rheology_Init( self, context );
   _YieldRheology_Init( self, strainWeakening, materialPointsSwarm, minVisc ); 
   _VonMises_Init( self, strainRateField, swarmStrainRate, cohesion, cohesionAfterSoftening, strainRateSoftening );
   _Byerlee_Init( self, mesh, depthCoefficient );
   self->isConstructed = True;
   return self;
}

/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
Byerlee* _Byerlee_New(
      SizeT                                              sizeOfSelf,
      Type                                               type,
      Stg_Class_DeleteFunction*                          _delete,
      Stg_Class_PrintFunction*                           _print,
      Stg_Class_CopyFunction*                            _copy,
      Stg_Component_DefaultConstructorFunction*          _defaultConstructor,
      Stg_Component_ConstructFunction*                   _construct,
      Stg_Component_BuildFunction*                       _build,
      Stg_Component_InitialiseFunction*                  _initialise,
      Stg_Component_ExecuteFunction*                     _execute,
      Stg_Component_DestroyFunction*                     _destroy,
      Rheology_ModifyConstitutiveMatrixFunction*         _modifyConstitutiveMatrix,
      YieldRheology_GetYieldCriterionFunction*           _getYieldCriterion,
      YieldRheology_GetYieldIndicatorFunction*           _getYieldIndicator,
      YieldRheology_HasYieldedFunction*                  _hasYielded,
      Name                                               name )
{
   Byerlee*             self;

   /* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
   assert( sizeOfSelf >= sizeof(Byerlee) );
   self = (Byerlee*) _VonMises_New(
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
         _modifyConstitutiveMatrix,
         _getYieldCriterion,
         _getYieldIndicator,
         _hasYielded,
         name );

   /* Function pointers for this class that are not on the parent class should be set here */

   return self;
}

void _Byerlee_Init( Byerlee* self, FeMesh* mesh, double depthCoefficient ) {
   self->mesh = mesh;
   self->depthCoefficient = depthCoefficient;
}

void* _Byerlee_DefaultNew( Name name ) {
   return (void*) _Byerlee_New(
      sizeof(Byerlee),
      Byerlee_Type,
      _YieldRheology_Delete,
      _YieldRheology_Print,
      _YieldRheology_Copy,
      _Byerlee_DefaultNew,
      _Byerlee_AssignFromXML,
      _YieldRheology_Build,
      _YieldRheology_Initialise,
      _YieldRheology_Execute,
      _Byerlee_Destroy,
      _YieldRheology_ModifyConstitutiveMatrix,
      _Byerlee_GetYieldCriterion,
      _VonMises_GetYieldIndicator,
      _VonMises_HasYielded,
      name );
}

void _Byerlee_AssignFromXML( void* rheology, Stg_ComponentFactory* cf, void* data ){
   Byerlee*            self           = (Byerlee*)rheology;
   FeMesh*         mesh;

   /* Construct Parent */
   _VonMises_AssignFromXML( self, cf, data );
   /* TODO: "KeyFallback' soon to be deprecated/updated */
   /*
   geometry = Stg_ComponentFactory_ConstructByNameWithKeyFallback(
                        cf, self->name, "geometry", "BlockGeometry", BlockGeometry, True, data );
   */
   mesh     = Stg_ComponentFactory_ConstructByNameWithKeyFallback(
                        cf, self->name, "mesh-linear", "FeMesh", FeMesh, True, data );
   /*geometry = Stg_ComponentFactory_ConstructByKey(
         cf, self->name, "BlockGeometry", BlockGeometry, True );
   mesh     = Stg_ComponentFactory_ConstructByKey(
         cf, self->name, "FiniteElement_Mesh", FiniteElement_Mesh, True );*/

   _Byerlee_Init(
         self,
         mesh,
         Stg_ComponentFactory_GetDouble( cf, self->name, "depthCoefficient", 0.0  ) );
}

void _Byerlee_Destroy( void* rheology, void* data ) {
   Byerlee*  self = (Byerlee*) rheology;
   
   Stg_Component_Destroy( self->mesh, data, False );
   
   _VonMises_Destroy( self, data );
   
}


double _Byerlee_GetYieldCriterion(
      void*                            rheology,
      ConstitutiveMatrix*              constitutiveMatrix,
      MaterialPointsSwarm*             materialPointsSwarm,
      Element_LocalIndex               lElement_I,
      MaterialPoint*                   materialPoint,
      Coord                            xi )
{
   Byerlee*                         self             = (Byerlee*) rheology;
   double                           depth = 0.0;
   double                           height;
   Coord           min, max;
   Coord                            coord;

   /* Calculate Depth */
   Mesh_GetGlobalCoordRange( self->mesh, min, max );
   height = max[ J_AXIS ];

   /* This rheology assumes particle is an integration points thats can be mapped to a particle
    * that has no meaningful coord. Best to re-calc the global from local */
   FeMesh_CoordLocalToGlobal( self->mesh, lElement_I, xi, coord );
   depth = height - coord[ J_AXIS ];

   return self->cohesion + self->depthCoefficient * depth;
}
