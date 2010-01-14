/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
** Copyright (c) 2005, Monash Cluster Computing 
** All rights reserved.
** Redistribution and use in source and binary forms, with or without modification,
** are permitted provided that the following conditions are met:
**
** 		* Redistributions of source code must retain the above copyright notice, 
** 			this list of conditions and the following disclaimer.
** 		* Redistributions in binary form must reproduce the above copyright 
**			notice, this list of conditions and the following disclaimer in the 
**			documentation and/or other materials provided with the distribution.
** 		* Neither the name of the Monash University nor the names of its contributors 
**			may be used to endorse or promote products derived from this software 
**			without specific prior written permission.
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
*%		Louis Moresi - Louis.Moresi@sci.monash.edu.au
*%
** Contributors:
*+		Robert Turnbull
*+		Vincent Lemiale
*+		Louis Moresi
*+		David May
*+		David Stegman
*+		Mirko Velic
*+		Patrick Sunter
*+		Julian Giordani
*+
** $Id: ViscosityFieldRheology.c 736 2008-05-23 06:05:21Z RobertTurnbull $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>

#include "types.h"
#include "RheologyClass.h"
#include "ViscosityFieldRheology.h"
#include "ConstitutiveMatrix.h"

#include <assert.h>

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type ViscosityFieldRheology_Type = "ViscosityFieldRheology";

/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
ViscosityFieldRheology* _ViscosityFieldRheology_New(  VISCOSITYFIELDRHEOLOGY_DEFARGS  ) 
{
	ViscosityFieldRheology*					self;

	assert( _sizeOfSelf >= sizeof(ViscosityFieldRheology) );
	self = (ViscosityFieldRheology*) _Rheology_New(  RHEOLOGY_PASSARGS  );
	
	return self;
}

void _ViscosityFieldRheology_Init( ViscosityFieldRheology* self, FeVariable* viscosityField ) 
{
	self->viscosityField = viscosityField;
}

void* _ViscosityFieldRheology_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                                     _sizeOfSelf = sizeof(ViscosityFieldRheology);
	Type                                                             type = ViscosityFieldRheology_Type;
	Stg_Class_DeleteFunction*                                     _delete = _Rheology_Delete;
	Stg_Class_PrintFunction*                                       _print = _Rheology_Print;
	Stg_Class_CopyFunction*                                         _copy = _Rheology_Copy;
	Stg_Component_DefaultConstructorFunction*         _defaultConstructor = _ViscosityFieldRheology_DefaultNew;
	Stg_Component_ConstructFunction*                           _construct = _ViscosityFieldRheology_AssignFromXML;
	Stg_Component_BuildFunction*                                   _build = _ViscosityFieldRheology_Build;
	Stg_Component_InitialiseFunction*                         _initialise = _ViscosityFieldRheology_Initialise;
	Stg_Component_ExecuteFunction*                               _execute = _Rheology_Execute;
	Stg_Component_DestroyFunction*                               _destroy = _ViscosityFieldRheology_Destroy;
	Rheology_ModifyConstitutiveMatrixFunction*  _modifyConstitutiveMatrix = _ViscosityFieldRheology_ModifyConstitutiveMatrix;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

	return (void*) _ViscosityFieldRheology_New(  VISCOSITYFIELDRHEOLOGY_PASSARGS  );
}

void _ViscosityFieldRheology_AssignFromXML( void* rheology, Stg_ComponentFactory* cf, void* data ){
	ViscosityFieldRheology*     self = (ViscosityFieldRheology*)rheology;
   FeVariable*                 viscosityField = NULL;
	/* Construct Parent */
	_Rheology_AssignFromXML( self, cf, data );
   
	viscosityField = Stg_ComponentFactory_ConstructByName( cf, (Name)Stg_ComponentFactory_GetString( cf, self->name, (Dictionary_Entry_Key)"ViscosityField", "ViscosityField"   ), FeVariable, True, data );
	_ViscosityFieldRheology_Init( self, viscosityField );
}

void _ViscosityFieldRheology_Build( void* rheology, void* data ){
	ViscosityFieldRheology*     self = (ViscosityFieldRheology*)rheology;

	_Rheology_Build( self, data );

	Stg_Component_Build( self->viscosityField, data, False );

}

void _ViscosityFieldRheology_Initialise( void* rheology, void* data ){
	ViscosityFieldRheology*     self = (ViscosityFieldRheology*)rheology;

	_Rheology_Initialise( self, data );

	Stg_Component_Initialise( self->viscosityField, data, False );

}

void _ViscosityFieldRheology_Destroy( void* rheology, void* data ){
	ViscosityFieldRheology*     self = (ViscosityFieldRheology*)rheology;

	_Rheology_Destroy( self, data );

	Stg_Component_Destroy( self->viscosityField, data, False );

}

void _ViscosityFieldRheology_ModifyConstitutiveMatrix( 
		void*                                              rheology, 
		ConstitutiveMatrix*                                constitutiveMatrix,
		MaterialPointsSwarm*                               swarm,
		Element_LocalIndex                                 lElement_I,
		MaterialPoint*                                     materialPoint,
		Coord                                              xi )
{
	ViscosityFieldRheology*	                  self              = (ViscosityFieldRheology*) rheology;
	FeMesh*                           mesh              = ConstitutiveMatrix_GetMesh( constitutiveMatrix );
	double viscosity;

	FeVariable_InterpolateFromMeshLocalCoord( self->viscosityField, mesh, lElement_I, xi, &viscosity );
	ConstitutiveMatrix_SetIsotropicViscosity( constitutiveMatrix, viscosity );
}


