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
** $Id: Anisotropic.c 610 2007-10-11 08:09:29Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>

#include "types.h"
#include "RheologyClass.h"
#include "Anisotropic.h"
#include "ConstitutiveMatrix.h"
#include "Director.h"

#include <assert.h>

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type Anisotropic_Type = "Anisotropic";

/* Public Constructor */
Anisotropic* Anisotropic_New(
	Name					name,
	AbstractContext*	context,
	Director*			director,
	double				viscosityRatio )
{
   Anisotropic* self = (Anisotropic*) _Anisotropic_DefaultNew( name );

   _Rheology_Init( self, (PICelleratorContext*)context );
   _Anisotropic_Init( self, director, viscosityRatio ) ;

   self->isConstructed = True;
   return self;
}

/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
Anisotropic* _Anisotropic_New(  ANISOTROPIC_DEFARGS  )
{
   Anisotropic*               self;

   /* Call private constructor of parent - this will set virtual functions of parent and continue up the hierarchy tree. At the beginning of the tree it will allocate memory of the size of object and initialise all the memory to zero. */
   assert( _sizeOfSelf >= sizeof(Anisotropic) );
   self = (Anisotropic*) _Rheology_New(  RHEOLOGY_PASSARGS  );

   return self;
}

void _Anisotropic_Init( Anisotropic* self, Director* director, double viscosityRatio ) {
   self->director        = director;
   self->viscosityRatio  = viscosityRatio;
}

void* _Anisotropic_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                                     _sizeOfSelf = sizeof(Anisotropic);
	Type                                                             type = Anisotropic_Type;
	Stg_Class_DeleteFunction*                                     _delete = _Rheology_Delete;
	Stg_Class_PrintFunction*                                       _print = _Rheology_Print;
	Stg_Class_CopyFunction*                                         _copy = _Rheology_Copy;
	Stg_Component_DefaultConstructorFunction*         _defaultConstructor = _Anisotropic_DefaultNew;
	Stg_Component_ConstructFunction*                           _construct = _Anisotropic_AssignFromXML;
	Stg_Component_BuildFunction*                                   _build = _Anisotropic_Build;
	Stg_Component_InitialiseFunction*                         _initialise = _Anisotropic_Initialise;
	Stg_Component_ExecuteFunction*                               _execute = _Rheology_Execute;
	Stg_Component_DestroyFunction*                               _destroy = _Anisotropic_Destroy;
	Rheology_ModifyConstitutiveMatrixFunction*  _modifyConstitutiveMatrix = _Anisotropic_ModifyConstitutiveMatrix;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

   return (void*) _Anisotropic_New(  ANISOTROPIC_PASSARGS  );
}

void _Anisotropic_AssignFromXML( void* rheology, Stg_ComponentFactory* cf, void* data ){
   Anisotropic*     self = (Anisotropic*)rheology;
   Director*        director;

   /* Construct Parent */
   _Rheology_AssignFromXML( self, cf, data );

   director =  Stg_ComponentFactory_ConstructByKey(  cf,  self->name,  "Director", Director,  True, data ) ;

   _Anisotropic_Init(
         self,
         director,
         Stg_ComponentFactory_GetDouble( cf, self->name, "viscosityRatio", 1.0 ) );
}

void _Anisotropic_Destroy( void* _self, void* data ) {
   Anisotropic* self = (Anisotropic*) _self;

   Stg_Component_Destroy( self->director, data, False );

   /* Destroy Parent */
   _Rheology_Destroy( self, data );

}

void _Anisotropic_Build( void* _self, void* data ) {
   Anisotropic* self = (Anisotropic*) _self;

   Stg_Component_Build( self->director, data, False );

   /* Build Parent */
   _Rheology_Build( self, data );

}

void _Anisotropic_Initialise( void* _self, void* data ) {
   Anisotropic* self = (Anisotropic*) _self;

   Stg_Component_Initialise( self->director, data, False );

   /* Build Parent */
   _Rheology_Initialise( self, data );

}


void _Anisotropic_ModifyConstitutiveMatrix(
      void*                                              rheology,
      ConstitutiveMatrix*                                constitutiveMatrix,
      MaterialPointsSwarm*                               swarm,
      Element_LocalIndex                                 lElement_I,
      MaterialPoint*                                     materialPoint,
      Coord                                              xi )
{
   Anisotropic*                   self               = (Anisotropic*) rheology;
   double                          isotropicViscosity = ConstitutiveMatrix_GetIsotropicViscosity( constitutiveMatrix );
   double                          deltaViscosity;
   XYZ                             normal;

   deltaViscosity = isotropicViscosity * (1.0 - self->viscosityRatio);
   Director_GetNormal( self->director, materialPoint, normal );

   ConstitutiveMatrix_SetSecondViscosity( constitutiveMatrix, deltaViscosity, normal );
}


