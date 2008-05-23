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
ViscosityFieldRheology* _ViscosityFieldRheology_New( 
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
		Name                                               name ) 
{
	ViscosityFieldRheology*					self;

	assert( sizeOfSelf >= sizeof(ViscosityFieldRheology) );
	self = (ViscosityFieldRheology*) _Rheology_New( 
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
			name );
	
	return self;
}

void _ViscosityFieldRheology_Init( ViscosityFieldRheology* self, Name viscosityFieldName ) 
{
	self->viscosityFieldName = viscosityFieldName;
}

void* _ViscosityFieldRheology_DefaultNew( Name name ) {
	return (void*) _ViscosityFieldRheology_New(
		sizeof(ViscosityFieldRheology),
		ViscosityFieldRheology_Type,
		_Rheology_Delete,
		_Rheology_Print,
		_Rheology_Copy,
		_ViscosityFieldRheology_DefaultNew,
		_ViscosityFieldRheology_Construct,
		_ViscosityFieldRheology_Build,
		_Rheology_Initialise,
		_Rheology_Execute,
		_Rheology_Destroy,
		_ViscosityFieldRheology_ModifyConstitutiveMatrix,
		name );
}

void _ViscosityFieldRheology_Construct( void* rheology, Stg_ComponentFactory* cf, void* data ){
	ViscosityFieldRheology*     self = (ViscosityFieldRheology*)rheology;

	/* Construct Parent */
	_Rheology_Construct( self, cf, data );
	
	_ViscosityFieldRheology_Init( self, 
		Stg_ComponentFactory_GetString( cf, self->name, "ViscosityField", "ViscosityField" ) );
}

void _ViscosityFieldRheology_Build( void* rheology, void* data ){
	ViscosityFieldRheology*     self = (ViscosityFieldRheology*)rheology;
	AbstractContext*            context     = Stg_CheckType( data, AbstractContext );
	Stg_ComponentFactory*       cf          = context->CF;

	_Rheology_Build( self, data );

	self->viscosityField = Stg_ComponentFactory_ConstructByName( cf, self->viscosityFieldName, FeVariable, True, data );

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
