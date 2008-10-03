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
** $Id:  $
** 
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include <Underworld/Underworld.h>
#include <assert.h>


const Type Underworld_DensityChange_Type = "Underworld_DensityChange";

typedef struct { 
	__Codelet
	BuoyancyForceTerm_MaterialExt* bftExt;
	IntegrationPointsSwarm*     swarm;
	Material*  material;
	double     height;
	double     newDensity;
} Underworld_DensityChange;

double Underworld_DensityChange_GetDiffusivityFromIntPoint( void* _residualForceTerm, void* particle ) {

	/* A function which should will only ever be called via a function pointer 
	 * hanging off the residualForceTerm data structure.
	 *
	 * Role: 
	 * Input: residualForceTerm and IntegrationPoint although both these are cast as void*
	 * Return: the diffusivity which is 'mapped' to that materialPoint
	 */
	AdvDiffResidualForceTerm* residualForceTerm   = (AdvDiffResidualForceTerm*)_residualForceTerm;
	IntegrationPointsSwarm*   integrationSwarm    = (IntegrationPointsSwarm*)residualForceTerm->integrationSwarm;
	IntegrationPoint*         integrationParticle = (IntegrationPoint*)particle;
	MaterialPointRef*         materialRef;

	materialRef = OneToOneMapper_GetMaterialRef( (OneToOneMapper*)integrationSwarm->mapper, integrationParticle );
	return Variable_GetValueDouble( residualForceTerm->diffusivityVariable, materialRef->particle_I );
}

void Underworld_lalala_Setup( UnderworldContext* context ) {
	/* 
	 * Role:
	 * 1) create a new variable on the material Swarm, called thermalDiffusivity.
	 * 2) modify the residualForceTerm used in the advection diffusion solver. The
	 * 	modifications enables the materialSwarm parameters for thermalDiffusivity to
	 * 	be used by the ResidualForceTerm instead of the default diffusivity
	 *
	 * Assumptions:
	 * 	A OneToOneMapper is used
	 *	AdvDiffResidualForceTerm setup and execution is compatible with it
	 */


	AdvectionDiffusionSLE*    energySLE           = context->energySLE;
        /* TODO: This assumes OneToOne mapping of intPoints to matPoints, should be fixed in future */
	OneToOneMapper*           mapper              = (OneToOneMapper*)context->picIntegrationPoints->mapper;

	MaterialPointsSwarm*      materialSwarm = mapper->materialSwarm;
	ExtensionInfo_Index       particleExtHandle;
	SwarmVariable*            swarmVariable;
	ForceVector*              residual;
	AdvDiffResidualForceTerm* residualForceTerm;

	assert( energySLE );
	assert( materialSwarm );

	 	

	/* Add Material Extension */
	particleExtHandle = ExtensionManager_Add( materialSwarm->particleExtensionMgr, 
				CURR_MODULE_NAME, sizeof( double ) );

	swarmVariable = Swarm_NewScalarVariable(
			materialSwarm,
			"thermalDiffusivity",
			(ArithPointer) ExtensionManager_Get( materialSwarm->particleExtensionMgr, 0, particleExtHandle ),
			Variable_DataType_Double );

	/* Set Pointers */
	residual = energySLE->residual;
	residualForceTerm = Stg_CheckType( Stg_ObjectList_At( residual->forceTermList, 0 ), AdvDiffResidualForceTerm );
	residualForceTerm->diffusivityVariable = swarmVariable->variable;
	residualForceTerm->integrationSwarm    = (Swarm*) context->picIntegrationPoints;
	/* Important that this function is defined here, but is used in 
	 * StgFEM/SLE/ProvidedSystems/AdvectionDiffusion/src/Residual.c
	 * see its _AdvDiffResidualForceTerm_AssembleElement for details
	 */
	residualForceTerm->_getDiffusivityFromIntPoint = Underworld_DensityChange_GetDiffusivityFromIntPoint;
}

void Underworld_DensityChange_Assign( UnderworldContext* context ) {
	/* 
	 * Role: 
	 * Assign diffusivity values given via the input dictionary to the materialPoints.
	 */ 
	Materials_Register*               materialRegister = context->materials_Register;
	Material*                         material;
	Material_Index                    material_I;
	Material_Index                    materialsCount = Materials_Register_GetCount( materialRegister );
	Dictionary*                       dictionary;
	double*                           materialThermalDiffusivity;
	Particle_Index                    lParticle_I;
	MaterialPointsSwarm*           	  materialSwarm;               
	AdvectionDiffusionSLE*            energySLE           = context->energySLE;
	ForceVector*                      residual;
	AdvDiffResidualForceTerm*         residualForceTerm;
	Variable*                         variable;

	materialSwarm = (MaterialPointsSwarm*)Stg_ComponentFactory_ConstructByName( context->CF, "materialSwarm", MaterialPointsSwarm, True, 0 /* dummy */ );
	assert(materialSwarm);
	residual = energySLE->residual;
	residualForceTerm = Stg_CheckType( Stg_ObjectList_At( residual->forceTermList, 0 ), AdvDiffResidualForceTerm );
	variable = residualForceTerm->diffusivityVariable;

	materialThermalDiffusivity = Memory_Alloc_Array( double, materialsCount, "materialThermalDiffusivity" );
	
	/* Loop over materials and get material properties from dictionary */
	for ( material_I = 0 ; material_I < materialsCount ; material_I++ ) {
		material    = Materials_Register_GetByIndex( materialRegister, material_I );
		/* Get the material's dictionary */
		dictionary  = material->dictionary;

		materialThermalDiffusivity[ material_I ] = 
			Dictionary_GetDouble_WithDefault( dictionary, "thermalDiffusivity", 1.0 );
	}

	/* Assign value to particle */
	for ( lParticle_I = 0 ; lParticle_I < materialSwarm->particleLocalCount ; lParticle_I++ ) {
		Variable_SetValueDouble( 
			variable, 
			lParticle_I, 
			materialThermalDiffusivity[ MaterialPointsSwarm_GetMaterialIndexAt( materialSwarm, lParticle_I ) ] );
	}

	Memory_Free( materialThermalDiffusivity );
}

void Underworld_DensityChange_Check( UnderworldContext* context ) {

	double volume;
	double centroid[3];
	static int densityChangeIsDone = 0;

	/* test if this has already happened, if so leave function */
	if ( densityChangeIsDone )
		return;

	/* Get self (the plugin) */
	Underworld_DensityChange* self = (Underworld_DensityChange*)LiveComponentRegister_Get(
					context->CF->LCRegister,
					Underworld_DensityChange_Type );

	/* get centroid coordinate */
	volume = Material_Volume( self->material, (IntegrationPointsSwarm*)self->swarm, centroid );

	/* test if centroid height (y-coord) is >= input height
	 * 	Note following code is only run once */
	if( centroid[1] >= self->height ) {
		self->bftExt->density = self->newDensity;
		densityChangeIsDone = 1;
		Journal_RPrintf( Journal_Register( Info_Type, "DensityChange" ), 
				"************** Density Change for material %s, new density is %g**************\n",
				self->material->name, self->bftExt->density );
	}

}
void Underworld_DensityChange_Setup( UnderworldContext* context ) {

	
	BuoyancyForceTerm*  bft = NULL;
	BuoyancyForceTerm_MaterialExt* materialExt = NULL;;
	Materials_Register*  materialRegister = context->materials_Register;
	Stream* stream = Journal_Register( Info_Type, "cows" );
	IntegrationPointsSwarm* swarm = NULL;
	Name   materialName = NULL;
	int materialIndex;
	double oldDensity;

	/* Get self (the plugin) */
	Underworld_DensityChange* self = (Underworld_DensityChange*)LiveComponentRegister_Get(
					context->CF->LCRegister,
					Underworld_DensityChange_Type );

	/* Initialise plugin data */
	self->bftExt = NULL;
	self->swarm = NULL;
	self->material = NULL;
	self->height = 0;
	self->newDensity = 0;

	/* Need buoyancy force term, to get density infomation from bft extension + the integration swarm */
	bft = Stg_ComponentFactory_ConstructByName( context->CF, "buoyancyForceTerm", BuoyancyForceTerm, True, NULL );

	/* Read in input from root xml dictionary */
	self->height = Dictionary_GetDouble_WithDefault( context->dictionary, "materialDensityChangeHeight", 0.5 );
	self->newDensity = Dictionary_GetDouble_WithDefault( context->dictionary, "materialDensityNewDensity", 0.3 );
	self->swarm = (IntegrationPointsSwarm*)bft->integrationSwarm;
	materialName = Dictionary_GetString( context->dictionary, "materialDensityToChange" );
	self->material = Materials_Register_GetByName( self->swarm->materials_Register, materialName );

	/* check if material index exists */
	if( self->material==NULL ) {
		printf("Error\nCounld find the material with index %d\n", materialIndex ); exit(0);
	}
	materialExt = ExtensionManager_Get( self->material->extensionMgr, self->material, bft->materialExtHandle );
	Journal_RPrintf( stream, "Will change %s's density at height %g from %g to %g\n", 
			self->material->name, self->height, materialExt->density, self->newDensity );

	self->bftExt = materialExt;
}



void _Underworld_DensityChange_Construct( void* component, Stg_ComponentFactory* cf, void* data ) {
	UnderworldContext* context;

	context = (UnderworldContext*)Stg_ComponentFactory_ConstructByName( cf, "context", UnderworldContext, True, data ); 

	/* Add functions to entry points */
	ContextEP_Append( context, AbstractContext_EP_Initialise,  Underworld_DensityChange_Setup );
	ContextEP_Append( context, AbstractContext_EP_FrequentOutput, 			Underworld_DensityChange_Check );
	
}

void* _Underworld_DensityChange_DefaultNew( Name name ) {
	return _Codelet_New(
			sizeof( Underworld_DensityChange ),
			Underworld_DensityChange_Type, 
			_Codelet_Delete, 
			_Codelet_Print, 
			_Codelet_Copy, 
			_Underworld_DensityChange_DefaultNew,
			_Underworld_DensityChange_Construct,
			_Codelet_Build,
			_Codelet_Initialise,
			_Codelet_Execute,
			_Codelet_Destroy,
			name );
	}

Index Underworld_DensityChangeAtDepth_Register( PluginsManager* pluginsManager ) {
	return PluginsManager_Submit( 
			pluginsManager, 
			Underworld_DensityChange_Type, 
			"0",
			_Underworld_DensityChange_DefaultNew );
}
