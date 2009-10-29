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
** $Id: Arrhenius.c 618 2007-10-29 07:53:04Z RobertTurnbull $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>

#include "types.h"
#include "RheologyClass.h"
#include "Arrhenius.h"
#include "ConstitutiveMatrix.h"

#include <assert.h>

/* Textual name of this class - This is a global pointer which is used for times when you need to refer to class and not a particular instance of a class */
const Type Arrhenius_Type = "Arrhenius";

/* Private Constructor: This will accept all the virtual functions for this class as arguments. */
Arrhenius* _Arrhenius_New( 
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
	Arrhenius*					self;

	assert( sizeOfSelf >= sizeof(Arrhenius) );
	self = (Arrhenius*) _Rheology_New( 
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

void _Arrhenius_Init( Arrhenius* self, FeVariable* temperatureField, double eta0, double activationEnergy, double activationVolume, double referenceTemp ) 
{
	self->temperatureField = temperatureField;
	self->eta0             = eta0;
	self->activationEnergy = activationEnergy;
	self->activationVolume = activationVolume;
	self->referenceTemp    = referenceTemp;
}

void* _Arrhenius_DefaultNew( Name name ) {
	return (void*) _Arrhenius_New(
		sizeof(Arrhenius),
		Arrhenius_Type,
		_Rheology_Delete,
		_Rheology_Print,
		_Rheology_Copy,
		_Arrhenius_DefaultNew,
		_Arrhenius_AssignFromXML,
		_Rheology_Build,
		_Rheology_Initialise,
		_Rheology_Execute,
		_Rheology_Destroy,
		_Arrhenius_ModifyConstitutiveMatrix,
		name );
}

void _Arrhenius_AssignFromXML( void* rheology, Stg_ComponentFactory* cf, void* data ){
	Arrhenius*     self = (Arrhenius*)rheology;
	FeVariable*    temperatureField;

	/* Construct Parent */
	_Rheology_AssignFromXML( self, cf, data );
	
	/* TODO: 'KeyFallback' soon to be deprecated/updated */
	temperatureField = Stg_ComponentFactory_ConstructByNameWithKeyFallback(
                        cf, self->name, "TemperatureField", "TemperatureField", FeVariable, True, data );
	/*temperatureField = Stg_ComponentFactory_ConstructByKey( 
			cf, self->name, "TemperatureField", FeVariable, True );*/

	_Arrhenius_Init( 
			self, 
			temperatureField, 
			Stg_ComponentFactory_GetDouble( cf, self->name, "eta0", 1.0 ),
			Stg_ComponentFactory_GetDouble( cf, self->name, "activationEnergy", 0.0 ),
			Stg_ComponentFactory_GetDouble( cf, self->name, "activationVolume", 0.0 ),
			Stg_ComponentFactory_GetDouble( cf, self->name, "referenceTemperature", 1.0 ) );
}

void _Arrhenius_ModifyConstitutiveMatrix( 
		void*                                              rheology, 
		ConstitutiveMatrix*                                constitutiveMatrix,
		MaterialPointsSwarm*                               swarm,
		Element_LocalIndex                                 lElement_I,
		MaterialPoint*                                     materialPoint,
		Coord                                              xi )
{
	Arrhenius*	                  self              = (Arrhenius*) rheology;
	FeVariable*                       temperatureField  = self->temperatureField;  
	double                            eta0;
	double                            activationEnergy;
	double                            activationVolume;
	double                            referenceTemp;
	double                            temperature;
	double                            viscosity;
	double                            depth;
	double                            height;
	Coord				  min, max;	
	Coord                             coord;
	FeMesh*                           mesh              = ConstitutiveMatrix_GetMesh( constitutiveMatrix );

	eta0             = self->eta0;
	activationEnergy = self->activationEnergy;
	activationVolume = self->activationVolume;
	referenceTemp    = self->referenceTemp;

	/* Extract geometric extents. */
	Mesh_GetGlobalCoordRange( mesh, min, max );

	/* Calculate Parameters */
	FeVariable_InterpolateFromMeshLocalCoord( temperatureField, mesh, lElement_I, xi, &temperature );
	
	/* If activationVolume is 0 there is no need to calculate the depth of the particle see viscosity line below. */
	if( activationVolume > (0.0 + 1e-12 )  ) {
		/* Calculate Depth */
		height = max[ J_AXIS ];

		/* This rheology assumes particle is an integration points thats can be mapped to a particle
		 * that has no meaningful coord. Best to re-calc the global from local */
		FeMesh_CoordLocalToGlobal( mesh, lElement_I, xi, coord );
		depth = height - coord[ J_AXIS ];
		/* Calculate New Viscosity */
		viscosity = eta0 * exp(( activationEnergy + activationVolume * depth)/ (temperature + referenceTemp));
	}
	else {
		/* Calculate New Viscosity */
		viscosity = eta0 * exp(( activationEnergy )/ (temperature + referenceTemp));
	}
	ConstitutiveMatrix_SetIsotropicViscosity( constitutiveMatrix, viscosity );
}
