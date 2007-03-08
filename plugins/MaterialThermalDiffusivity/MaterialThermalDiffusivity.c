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
** $Id: MaterialThermalDiffusivity.c 358 2006-10-18 06:17:30Z SteveQuenette $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include <Underworld/Underworld.h>
#include <assert.h>

void Underworld_MaterialThermalDiffusivity_Setup( UnderworldContext* context ) {
	AdvectionDiffusionSLE*    energySLE           = context->energySLE;
	IntegrationPointsSwarm*   picIntegrationPoints      = context->picIntegrationPoints;
	ExtensionInfo_Index       particleExtHandle;
	SwarmVariable*            swarmVariable;
	ForceVector*              residual;
	AdvDiffResidualForceTerm* residualForceTerm;

	assert( energySLE );
	assert( picIntegrationPoints );

	/* Add Material Extension */
	particleExtHandle = ExtensionManager_Add( picIntegrationPoints->particleExtensionMgr, 
				CURR_MODULE_NAME, sizeof( double ) );

	swarmVariable = Swarm_NewScalarVariable(
			picIntegrationPoints,
			"thermalDiffusivity",
			(ArithPointer) ExtensionManager_Get( picIntegrationPoints->particleExtensionMgr, 0, particleExtHandle ),
			Variable_DataType_Double );

	/* Set Pointers */
	residual = energySLE->residual;
	residualForceTerm = Stg_CheckType( Stg_ObjectList_At( residual->forceTermList, 0 ), AdvDiffResidualForceTerm );
	residualForceTerm->diffusivityVariable = swarmVariable->variable;
	residualForceTerm->integrationSwarm    = (Swarm*) picIntegrationPoints;
}

void Underworld_MaterialThermalDiffusivity_Assign( UnderworldContext* context ) {
	Materials_Register*               materialRegister = context->materials_Register;
	Material*                         material;
	Material_Index                    material_I;
	Material_Index                    materialsCount = Materials_Register_GetCount( materialRegister );
	Dictionary*                       dictionary;
	double*                           materialThermalDiffusivity;
	Particle_Index                    lParticle_I;
	MaterialPointsSwarm*           	  materialPoints;               
	AdvectionDiffusionSLE*            energySLE           = context->energySLE;
	ForceVector*                      residual;
	AdvDiffResidualForceTerm*         residualForceTerm;
	Variable*                         variable;

	materialPoints = (MaterialPointsSwarm*)Stg_ComponentFactory_ConstructByName( context->CF, "materialPoints", MaterialPointsSwarm, True, 0 /* dummy */ );
	assert(materialPoints);
	residual = energySLE->residual;
	residualForceTerm = Stg_CheckType( Stg_ObjectList_At( residual->forceTermList, 0 ), AdvDiffResidualForceTerm );
	variable = residualForceTerm->diffusivityVariable;

	materialThermalDiffusivity = Memory_Alloc_Array( double, materialsCount, "materialThermalDiffusivity" );
	
	/* Loop over materials and get material properties from dictionary */
	for ( material_I = 0 ; material_I < materialsCount ; material_I++ ) {
		material    = Materials_Register_GetByIndex( materialRegister, material_I );
		dictionary  = material->dictionary;

		materialThermalDiffusivity[ material_I ] = 
			Dictionary_GetDouble_WithDefault( dictionary, "thermalDiffusivity", 1.0 );
	}

	/* Assign value to particle */
	for ( lParticle_I = 0 ; lParticle_I < materialPoints->particleLocalCount ; lParticle_I++ ) {
		Variable_SetValueDouble( 
			variable, 
			lParticle_I, 
			materialThermalDiffusivity[ MaterialPointsSwarm_GetMaterialIndexAt( materialPoints, lParticle_I ) ] );
	}

	Memory_Free( materialThermalDiffusivity );
}


const Type Underworld_MaterialThermalDiffusivity_Type = "Underworld_MaterialThermalDiffusivity";

typedef struct { 
	__Codelet
} Underworld_MaterialThermalDiffusivity;

void _Underworld_MaterialThermalDiffusivity_Construct( void* component, Stg_ComponentFactory* cf, void* data ) {
	UnderworldContext* context;

	context = (UnderworldContext*)Stg_ComponentFactory_ConstructByName( cf, "context", UnderworldContext, True, data ); 

	/* Add functions to entry points */
	ContextEP_Append( context, AbstractContext_EP_ConstructExtensions, Underworld_MaterialThermalDiffusivity_Setup );
	ContextEP_Append( context, AbstractContext_EP_Initialise,          Underworld_MaterialThermalDiffusivity_Assign );

}

void* _Underworld_MaterialThermalDiffusivity_DefaultNew( Name name ) {
	return Codelet_New(
			Underworld_MaterialThermalDiffusivity_Type, 
			_Underworld_MaterialThermalDiffusivity_DefaultNew,
			_Underworld_MaterialThermalDiffusivity_Construct,
			_Codelet_Build,
			_Codelet_Initialise,
			_Codelet_Execute,
			_Codelet_Destroy,
			name );
	}

Index Underworld_MaterialThermalDiffusivity_Register( PluginsManager* pluginsManager ) {
	return PluginsManager_Submit( 
			pluginsManager, 
			Underworld_MaterialThermalDiffusivity_Type, 
			"0",
			_Underworld_MaterialThermalDiffusivity_DefaultNew );
}
