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
** $Id: DepthDependentViscosity.c 610 2007-10-11 08:09:29Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>

#include "types.h"
#include "RheologyClass.h"
#include "ConstitutiveMatrix.h"
#include "DepthDependentViscosity.h"

#include <assert.h>
#include <string.h>

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type DepthDependentViscosity_Type = "DepthDependentViscosity";

/* Public Constructor */
DepthDependentViscosity* DepthDependentViscosity_New(
      Name                  name,
      AbstractContext*      context,
      FeMesh*               feMesh,
      double                eta0,
      double                gamma,
      Axis                  variationAxis,
      double                referencePoint )
{
   DepthDependentViscosity* self = (DepthDependentViscosity*) _DepthDependentViscosity_DefaultNew( name );

   _Rheology_Init( self, (PICelleratorContext*)context );
   _DepthDependentViscosity_Init( self, feMesh, eta0, gamma, variationAxis, referencePoint );
   self->isConstructed = True;
   return self;
}

/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
DepthDependentViscosity* _DepthDependentViscosity_New(  DEPTHDEPENDENTVISCOSITY_DEFARGS  )
{
   DepthDependentViscosity*               self;

   /* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
   assert( _sizeOfSelf >= sizeof(DepthDependentViscosity) );
   self = (DepthDependentViscosity*) _Rheology_New(  RHEOLOGY_PASSARGS  );

   return self;
}

void _DepthDependentViscosity_Init( DepthDependentViscosity* self, FeMesh* feMesh, double eta0, double gamma, Axis variationAxis, double referencePoint ) {
   self->feMesh          = feMesh;
   self->eta0            = eta0;
   self->gamma           = gamma;
   self->variationAxis   = variationAxis;
   self->referencePoint  = referencePoint;
}

void* _DepthDependentViscosity_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                                     _sizeOfSelf = sizeof(DepthDependentViscosity);
	Type                                                             type = DepthDependentViscosity_Type;
	Stg_Class_DeleteFunction*                                     _delete = _Rheology_Delete;
	Stg_Class_PrintFunction*                                       _print = _Rheology_Print;
	Stg_Class_CopyFunction*                                         _copy = _Rheology_Copy;
	Stg_Component_DefaultConstructorFunction*         _defaultConstructor = _DepthDependentViscosity_DefaultNew;
	Stg_Component_ConstructFunction*                           _construct = _DepthDependentViscosity_AssignFromXML;
	Stg_Component_BuildFunction*                                   _build = _Rheology_Build;
	Stg_Component_InitialiseFunction*                         _initialise = _Rheology_Initialise;
	Stg_Component_ExecuteFunction*                               _execute = _Rheology_Execute;
	Stg_Component_DestroyFunction*                               _destroy = _DepthDependentViscosity_Destroy;
	Rheology_ModifyConstitutiveMatrixFunction*  _modifyConstitutiveMatrix = _DepthDependentViscosity_ModifyConstitutiveMatrix;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

   return (void*) _DepthDependentViscosity_New(  DEPTHDEPENDENTVISCOSITY_PASSARGS  );
}

void _DepthDependentViscosity_AssignFromXML( void* rheology, Stg_ComponentFactory* cf, void* data ){
   DepthDependentViscosity*  self                   = (DepthDependentViscosity*)rheology;
   FeMesh*          feMesh;
   Axis                      variationAxis          = 0;
   Name                      variationAxisName;
   Stream*                   errorStream            = Journal_Register( Error_Type, self->type );

   /* Construct Parent */
   _Rheology_AssignFromXML( self, cf, data );

   feMesh = Stg_ComponentFactory_ConstructByKey( cf, self->name, "Mesh", FeMesh, True, data ) ;

   variationAxisName = Stg_ComponentFactory_GetString( cf, self->name, "variationAxis", "Y" );

   Journal_Firewall(
         variationAxisName != NULL && strlen( variationAxisName ) >= 1,
         errorStream,
         "Error in func %s for %s '%s' - Bad 'variationAxis'\n",
         __func__, self->type, self->name );

   switch ( variationAxisName[0] ) {
      case 'X': case 'x': case 'I': case 'i': case '0':
         variationAxis = I_AXIS; break;
      case 'Y': case 'y': case 'J': case 'j': case '1':
         variationAxis = J_AXIS; break;
      case 'Z': case 'z': case 'K': case 'k': case '2':
         variationAxis = K_AXIS; break;
      default:
         Journal_Firewall(
            False,
            errorStream,
            "Error in func %s for %s '%s' - Bad 'variationAxis' in dictionary.\n",
            __func__, self->type, self->name );
   }

   _DepthDependentViscosity_Init(
         self,
         feMesh,
         Stg_ComponentFactory_GetDouble( cf, self->name, "eta0", 1.0 ),
         Stg_ComponentFactory_GetDouble( cf, self->name, "gamma", 0.0 ),
         variationAxis,
         Stg_ComponentFactory_GetDouble( cf, self->name, "referencePoint", 0.0 ) );
}

void _DepthDependentViscosity_Destroy( void* rheology, void* data ) {
	DepthDependentViscosity* self = (DepthDependentViscosity*) rheology;

	Stg_Component_Destroy( self->feMesh, data, False );
	/* Destroy parent */
	_Rheology_Destroy( self, data );

}


void _DepthDependentViscosity_ModifyConstitutiveMatrix(
      void*                                              rheology,
      ConstitutiveMatrix*                                constitutiveMatrix,
      MaterialPointsSwarm*                               swarm,
      Element_LocalIndex                                 lElement_I,
      MaterialPoint*                                     materialPoint,
      Coord                                              xi )
{
   DepthDependentViscosity*     self              = (DepthDependentViscosity*) rheology;
   double                            distance;
   double                            viscosity;
   Coord                             coord;

   /* This rheology assumes particle is an integration points thats can be mapped to a particle
    * that has no meaningful coord. Best to re-calc the global from local */

   /* Calculate distance */
   FeMesh_CoordLocalToGlobal( swarm->mesh, lElement_I, xi, coord );
   distance = coord[ self->variationAxis ];

   /* Calculate New Viscosity */
   viscosity = self->eta0 * exp( self->gamma * ( distance - self->referencePoint ) );
   ConstitutiveMatrix_SetIsotropicViscosity( constitutiveMatrix, viscosity );
}


