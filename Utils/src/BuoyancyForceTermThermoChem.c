/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
**
** Copyright (C), 2003-2006, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street,
**	Melbourne, 3053, Australia.
** Copyright (c) 2005-2006, Monash Cluster Computing, Building 28, Monash University Clayton Campus,
**	Victoria, 3800, Australia
**
** Primary Contributing Organisations:
**	Victorian Partnership for Advanced Computing Ltd, Computational Software Development - http://csd.vpac.org
**	Australian Computational Earth Systems Simulator - http://www.access.edu.au
**	Monash Cluster Computing - http://www.mcc.monash.edu.au
**
** Contributors:
**	Robert Turnbull, Research Assistant, Monash University. (robert.turnbull@sci.monash.edu.au)
**	Patrick D. Sunter, Software Engineer, VPAC. (patrick@vpac.org)
**	Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
**	Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
**	David May, PhD Student, Monash University (david.may@sci.monash.edu.au)
**	Vincent Lemiale, Postdoctoral Fellow, Monash University. (vincent.lemiale@sci.monash.edu.au)
**	Julian Giordani, Research Assistant, Monash University. (julian.giordani@sci.monash.edu.au)
**	Louis Moresi, Associate Professor, Monash University. (louis.moresi@sci.monash.edu.au)
**	Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
**	Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
**	David Stegman, Postdoctoral Fellow, Monash University. (david.stegman@sci.monash.edu.au)
**	Wendy Sharples, PhD Student, Monash University (wendy.sharples@sci.monash.edu.au)
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
** $Id: BuoyancyForceTermThermoChem.c 376 2006-10-18 06:58:41Z SteveQuenette $
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/


#include <mpi.h>
#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PopulationControl/PopulationControl.h>
#include <PICellerator/Weights/Weights.h>
#include <PICellerator/MaterialPoints/MaterialPoints.h>

#include "types.h"
#include "BuoyancyForceTermThermoChem.h"
#include "MaterialSwarmVariable.h"

#include <assert.h>
#include <string.h>

/* Textual name of this class */
const Type BuoyancyForceTermThermoChem_Type = "BuoyancyForceTermThermoChem";

BuoyancyForceTermThermoChem* BuoyancyForceTermThermoChem_New( 
	Name							name,
	FiniteElementContext*	context,
	ForceVector*				forceVector,
	Swarm*						integrationSwarm,
	FeVariable*					temperatureField,
	double						RaT,
	double						RaC,
	Bool							adjust,
	Materials_Register*		materials_Register )
{
	BuoyancyForceTermThermoChem* self = (BuoyancyForceTermThermoChem*) _BuoyancyForceTermThermoChem_DefaultNew( name );

	self->isConstructed = True;
	_ForceTerm_Init( self, context, forceVector, integrationSwarm, NULL );
	_BuoyancyForceTermThermoChem_Init( self, temperatureField, RaT, RaC, adjust, materials_Register );

	return self;
}

void* _BuoyancyForceTermThermoChem_DefaultNew( Name name ) {
	/* Variables set in this function */
	SizeT                                                 _sizeOfSelf = sizeof(BuoyancyForceTermThermoChem);
	Type                                                         type = BuoyancyForceTermThermoChem_Type;
	Stg_Class_DeleteFunction*                                 _delete = _BuoyancyForceTermThermoChem_Delete;
	Stg_Class_PrintFunction*                                   _print = _BuoyancyForceTermThermoChem_Print;
	Stg_Class_CopyFunction*                                     _copy = NULL;
	Stg_Component_DefaultConstructorFunction*     _defaultConstructor = _BuoyancyForceTermThermoChem_DefaultNew;
	Stg_Component_ConstructFunction*                       _construct = _BuoyancyForceTermThermoChem_AssignFromXML;
	Stg_Component_BuildFunction*                               _build = _BuoyancyForceTermThermoChem_Build;
	Stg_Component_InitialiseFunction*                     _initialise = _BuoyancyForceTermThermoChem_Initialise;
	Stg_Component_ExecuteFunction*                           _execute = _BuoyancyForceTermThermoChem_Execute;
	Stg_Component_DestroyFunction*                           _destroy = _BuoyancyForceTermThermoChem_Destroy;
	ForceTerm_AssembleElementFunction*               _assembleElement = _BuoyancyForceTermThermoChem_AssembleElement;
	BuoyancyForceTermThermoChem_CalcRaTFunction*             _calcRaT = _BuoyancyForceTermThermoChem_CalcRaT;
	BuoyancyForceTermThermoChem_CalcRaCFunction*             _calcRaC = _BuoyancyForceTermThermoChem_CalcRaC;

	/* Variables that are set to ZERO are variables that will be set either by the current _New function or another parent _New function further up the hierachy */
	AllocationType  nameAllocationType = NON_GLOBAL /* default value NON_GLOBAL */;

	return (void*)_BuoyancyForceTermThermoChem_New(  BUOYANCYFORCETERMTHERMOCHEM_PASSARGS  );
}

/* Creation implementation / Virtual constructor */
BuoyancyForceTermThermoChem* _BuoyancyForceTermThermoChem_New(  BUOYANCYFORCETERMTHERMOCHEM_DEFARGS  )
{
	BuoyancyForceTermThermoChem* self;
	
	/* Allocate memory */
	assert( _sizeOfSelf >= sizeof(BuoyancyForceTermThermoChem) );
	/* The following terms are parameters that have been passed into this function but are being set before being passed onto the parent */
	/* This means that any values of these parameters that are passed into this function are not passed onto the parent function
	   and so should be set to ZERO in any children of this class. */
	nameAllocationType = NON_GLOBAL;

	self = (BuoyancyForceTermThermoChem*) _ForceTerm_New(  FORCETERM_PASSARGS  );
	
	/* Virtual info */
	self->_calcRaT = _calcRaT;
	self->_calcRaC = _calcRaC;
	
	return self;
}

void _BuoyancyForceTermThermoChem_Init( 
	BuoyancyForceTermThermoChem*	self, 
	FeVariable*							temperatureField,
	double								RaT,
	double								RaC,
	Bool									adjust,
	Materials_Register*				materials_Register )
{
	self->temperatureField    = temperatureField;
	self->RaT                 = RaT;
	self->RaC                 = RaC;
	self->adjust              = adjust;
	self->materials_Register  = materials_Register;
}

void _BuoyancyForceTermThermoChem_Delete( void* forceTerm ) {
	BuoyancyForceTermThermoChem* self = (BuoyancyForceTermThermoChem*)forceTerm;

   Memory_Free( self->densitySwarmVariables );

	_ForceTerm_Delete( self );
}

void _BuoyancyForceTermThermoChem_Print( void* forceTerm, Stream* stream ) {
	BuoyancyForceTermThermoChem* self = (BuoyancyForceTermThermoChem*)forceTerm;
	
	_ForceTerm_Print( self, stream );

	/* General info */
	Journal_PrintPointer( stream, self->temperatureField );
	Journal_PrintDouble( stream, self->RaT );
	Journal_PrintDouble( stream, self->RaC );
}

void _BuoyancyForceTermThermoChem_AssignFromXML( void* forceTerm, Stg_ComponentFactory* cf, void* data ) {
	BuoyancyForceTermThermoChem*	self = (BuoyancyForceTermThermoChem*)forceTerm;
	FeVariable*							temperatureField;
	double								RaT;
	double								RaC;
	Bool									adjust;
	Materials_Register*				materials_Register;
	PICelleratorContext*				context;

	/* Construct Parent */
	_ForceTerm_AssignFromXML( self, cf, data );

	temperatureField = Stg_ComponentFactory_ConstructByKey( cf, self->name, (Dictionary_Entry_Key)"TemperatureField", FeVariable, False, data  ) ;
	RaT = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"RaT", 0.0  );
	RaC = Stg_ComponentFactory_GetDouble( cf, self->name, (Dictionary_Entry_Key)"RaC", 0.0  );
	adjust = Stg_ComponentFactory_GetBool( cf, self->name, (Dictionary_Entry_Key)"adjust", False );

	context = (PICelleratorContext* )self->context;
	assert( Stg_CheckType( context, PICelleratorContext ) );
	materials_Register = context->materials_Register; 
	assert( materials_Register );

	_BuoyancyForceTermThermoChem_Init( self, temperatureField, RaT, RaC, adjust, materials_Register );
}

void _BuoyancyForceTermThermoChem_Build( void* forceTerm, void* data ) {
	BuoyancyForceTermThermoChem*               self               = (BuoyancyForceTermThermoChem*)forceTerm;
	BuoyancyForceTermThermoChem_MaterialExt*   materialExt;
	Material_Index                   material_I;
	Material*                        material;
	Materials_Register*              materials_Register = self->materials_Register;
	IntegrationPointsSwarm*          swarm              = (IntegrationPointsSwarm*)self->integrationSwarm;
	MaterialPointsSwarm**            materialSwarms;
	Index                            materialSwarm_I;
	Name                             name;

	_ForceTerm_Build( self, data );

	if ( self->temperatureField )
		Stg_Component_Build( self->temperatureField, data, False );

	/* Sort out material extension stuff */
	self->materialExtHandle = Materials_Register_AddMaterialExtension( 
			self->materials_Register, 
			self->type, 
			sizeof(BuoyancyForceTermThermoChem_MaterialExt) );
	for ( material_I = 0 ; material_I < Materials_Register_GetCount( materials_Register ) ; material_I++) {
		material = Materials_Register_GetByIndex( materials_Register, material_I );
		materialExt = ExtensionManager_GetFunc( material->extensionMgr, material, self->materialExtHandle );

		materialExt->density = Dictionary_GetDouble_WithDefault( material->dictionary, (Dictionary_Entry_Key)"density", 0.0  );
	}
	
	/* Create Swarm Variables of each material swarm this ip swarm is mapped against */
	materialSwarms = IntegrationPointMapper_GetMaterialPointsSwarms( swarm->mapper, &(self->materialSwarmCount) );
	self->densitySwarmVariables = Memory_Alloc_Array( MaterialSwarmVariable*, self->materialSwarmCount, "DensityVariables" );
	
	for ( materialSwarm_I = 0; materialSwarm_I < self->materialSwarmCount; ++materialSwarm_I ) {
		name = Stg_Object_AppendSuffix( materialSwarms[materialSwarm_I], (Name)"Density"  );
		self->densitySwarmVariables[materialSwarm_I] = MaterialSwarmVariable_New( 
				name,
				(AbstractContext*) self->context,
				materialSwarms[materialSwarm_I], 
				1, 
				self->materials_Register, 
				self->materialExtHandle, 
				GetOffsetOfMember( *materialExt, density ) );
		Memory_Free( name );
		
		/* Build new Swarm Variables */
		Stg_Component_Build( self->densitySwarmVariables[materialSwarm_I], data, False );
	}


}

void _BuoyancyForceTermThermoChem_Initialise( void* forceTerm, void* data ) {
	BuoyancyForceTermThermoChem*             self             = (BuoyancyForceTermThermoChem*)forceTerm;
	Index                          i;

	_ForceTerm_Initialise( self, data );

	if ( self->temperatureField )
		Stg_Component_Initialise( self->temperatureField, data, False );
	
	for ( i = 0; i < self->materialSwarmCount; ++i ) {
		Stg_Component_Initialise( self->densitySwarmVariables[i], data, False );
	}
}

void _BuoyancyForceTermThermoChem_Execute( void* forceTerm, void* data ) {
	_ForceTerm_Execute( forceTerm, data );
}

void _BuoyancyForceTermThermoChem_Destroy( void* forceTerm, void* data ) {
	BuoyancyForceTermThermoChem* self = (BuoyancyForceTermThermoChem*)forceTerm;
	Index i;

	for ( i = 0; i < self->materialSwarmCount; ++i ) {
		Stg_Component_Destroy( self->densitySwarmVariables[i], data, False );
	}

	if ( self->temperatureField )
		Stg_Component_Destroy( self->temperatureField, data, False );

	_ForceTerm_Destroy( self, data );
}


void _BuoyancyForceTermThermoChem_AssembleElement( void* forceTerm, ForceVector* forceVector, Element_LocalIndex lElement_I, double* elForceVec ) {
	BuoyancyForceTermThermoChem*               self               = (BuoyancyForceTermThermoChem*) forceTerm;
	IntegrationPoint*                particle;
	BuoyancyForceTermThermoChem_MaterialExt*   materialExt;
	Particle_InCellIndex             cParticle_I;
	Particle_InCellIndex             cellParticleCount;
	Element_NodeIndex                elementNodeCount;
	Dimension_Index                  dim                = forceVector->dim;
	IntegrationPointsSwarm*          swarm              = (IntegrationPointsSwarm*)self->integrationSwarm;
	FeMesh*		                 mesh               = forceVector->feVariable->feMesh;
	Node_ElementLocalIndex           eNode_I;
	Cell_Index                       cell_I;
	ElementType*                     elementType;
	Dof_Index                        nodeDofCount;
	double                           RaT;
	double                           RaC;
	double                           detJac             = 0.0;
	double                           factor;
	double                           Ni[27];
	double                           force;
	double*                          xi;
	Material*                        material;
	FeVariable*                      temperatureField   = self->temperatureField;
	double                           temperature        = 0.0;

	double totalWeight = 0.0;
	double adjustFactor = 0.0;

	elementType       = FeMesh_GetElementType( mesh, lElement_I );
	elementNodeCount  = elementType->nodeCount;
	nodeDofCount      = dim;
	cell_I            = CellLayout_MapElementIdToCellId( swarm->cellLayout, lElement_I );
	cellParticleCount = swarm->cellParticleCountTbl[cell_I];

	/* adjust & adjustFactor -- 20060411 Alan
	 *
	 * The adjust decides whether an adjustFactor should be applied to the resulting factor.
	 * If on, the total weight of the particles in the cell are scaled to the cell local volume.
	 *
	 * This is designed to be used when integrating with swarms which do not cover the whole domain
	 * (ie - I used it to do dave.m's test of 1 swarm for blob, 1 swarm for background)
	 */ 
	if ( self->adjust ) {
		for( cParticle_I = 0 ; cParticle_I < cellParticleCount ; cParticle_I++ ) {
			particle = (IntegrationPoint*) Swarm_ParticleInCellAt( swarm, cell_I, cParticle_I );
			totalWeight += particle->weight;
		}
		adjustFactor = swarm->weights->cellLocalVolume / totalWeight;
	}
	else {
		adjustFactor = 1.0;
	}
			
	for( cParticle_I = 0 ; cParticle_I < cellParticleCount ; cParticle_I++ ) {
		particle = (IntegrationPoint*) Swarm_ParticleInCellAt( swarm, cell_I, cParticle_I );
		xi       = particle->xi;

		detJac = ElementType_JacobianDeterminant( elementType, mesh, lElement_I, xi, dim );
		ElementType_EvaluateShapeFunctionsAt( elementType, xi, Ni );

		/* Get parameters */
		if ( temperatureField )
			FeVariable_InterpolateFromMeshLocalCoord( temperatureField, mesh, lElement_I, xi, &temperature );

		material = IntegrationPointsSwarm_GetMaterialOn( (IntegrationPointsSwarm*) swarm, particle );
		materialExt = ExtensionManager_Get( material->extensionMgr, material, self->materialExtHandle );

		/* Calculate Force */
		RaT = BuoyancyForceTermThermoChem_CalcRaT( self, (Swarm*)swarm, lElement_I, particle );
		RaC = BuoyancyForceTermThermoChem_CalcRaC( self, (Swarm*)swarm, lElement_I, particle );

		force =  ( RaT * temperature) - ( materialExt->density * RaC );
		factor = detJac * particle->weight * adjustFactor * force;

		/* Apply force in verticle direction */
		for( eNode_I = 0 ; eNode_I < elementNodeCount; eNode_I++ ) { 		
			elForceVec[ eNode_I * nodeDofCount + J_AXIS ] += factor * Ni[ eNode_I ] ;
		}
	}
	
}

/* The default implementation is for the RaC to be constant. */
double _BuoyancyForceTermThermoChem_CalcRaT( void* forceTerm, Swarm* swarm, Element_DomainIndex dElement_I, void* particle ) {
	BuoyancyForceTermThermoChem*               self               = (BuoyancyForceTermThermoChem*) forceTerm;

	return self->RaT;
}

double _BuoyancyForceTermThermoChem_CalcRaC( void* forceTerm, Swarm* swarm, Element_DomainIndex dElement_I, void* particle ) {
	BuoyancyForceTermThermoChem*               self               = (BuoyancyForceTermThermoChem*) forceTerm;

	return self->RaC;
}


