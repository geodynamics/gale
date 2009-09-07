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
** $Id: BuoyancyForceTerm.c 518 2007-10-11 08:07:50Z SteveQuenette $
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
#include "BuoyancyForceTerm.h"
#include "MaterialSwarmVariable.h"

#include <assert.h>
#include <string.h>

/* Textual name of this class */
const Type BuoyancyForceTerm_Type = "BuoyancyForceTerm";

BuoyancyForceTerm* BuoyancyForceTerm_New( 
		Name                                                name,
		ForceVector*                                        forceVector,
		Swarm*                                              integrationSwarm,
		FeVariable*                                         temperatureField,
		double                                              gravity,
		Bool                                                adjust,
		Materials_Register*                                 materials_Register )
{
	BuoyancyForceTerm* self = (BuoyancyForceTerm*) _BuoyancyForceTerm_DefaultNew( name );

	BuoyancyForceTerm_InitAll( 
			self,
			forceVector,
			integrationSwarm,
			temperatureField,
			gravity,
			adjust,
			materials_Register );

	return self;
}

/* Creation implementation / Virtual constructor */
BuoyancyForceTerm* _BuoyancyForceTerm_New( 
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
		ForceTerm_AssembleElementFunction*                  _assembleElement,		
		BuoyancyForceTerm_CalcGravityFunction*              _calcGravity,
		Name                                                name )
{
	BuoyancyForceTerm* self;
	
	/* Allocate memory */
	assert( sizeOfSelf >= sizeof(BuoyancyForceTerm) );
	self = (BuoyancyForceTerm*) _ForceTerm_New( 
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
	
	/* Virtual info */
	self->_calcGravity = _calcGravity;
	
	return self;
}

void _BuoyancyForceTerm_Init( 
		BuoyancyForceTerm*                                  self, 
		FeVariable*                                         temperatureField,
		double                                              gravity,
		Bool                                                adjust,
		Materials_Register*                                 materials_Register )
{
	self->temperatureField    = temperatureField;
	self->gravity             = gravity;
	self->gHat		  = NULL;
	self->adjust              = adjust;
	self->materials_Register  = materials_Register;
}

void BuoyancyForceTerm_InitAll( 
		void*                                               forceTerm,
		ForceVector*                                        forceVector,
		Swarm*                                              integrationSwarm,
		FeVariable*                                         temperatureField,
		double                                              gravity,
		Bool                                                adjust,
		Materials_Register*                                 materials_Register )
{
	BuoyancyForceTerm* self = (BuoyancyForceTerm*) forceTerm;

	ForceTerm_InitAll( self, forceVector, integrationSwarm, NULL );
	_BuoyancyForceTerm_Init( self, temperatureField, gravity, adjust, materials_Register );
}

void _BuoyancyForceTerm_Delete( void* forceTerm ) {
	BuoyancyForceTerm* self = (BuoyancyForceTerm*)forceTerm;
	Index i;

	FreeArray( self->gHat );
	for ( i = 0; i < self->materialSwarmCount; ++i ) {
		Stg_Class_Delete( self->densitySwarmVariables[i] );
		Stg_Class_Delete( self->alphaSwarmVariables[i] );
	}
	Stg_Class_Delete( self->densitySwarmVariables );
	Stg_Class_Delete( self->alphaSwarmVariables );

	_ForceTerm_Delete( self );
}

void _BuoyancyForceTerm_Print( void* forceTerm, Stream* stream ) {
	BuoyancyForceTerm* self = (BuoyancyForceTerm*)forceTerm;
	
	_ForceTerm_Print( self, stream );

	/* General info */
	Journal_PrintPointer( stream, self->temperatureField );
	Journal_PrintDouble( stream, self->gravity );
}

void* _BuoyancyForceTerm_DefaultNew( Name name ) {
	return (void*)_BuoyancyForceTerm_New( 
		sizeof(BuoyancyForceTerm), 
		BuoyancyForceTerm_Type,
		_BuoyancyForceTerm_Delete,
		_BuoyancyForceTerm_Print,
		NULL,
		_BuoyancyForceTerm_DefaultNew,
		_BuoyancyForceTerm_Construct,
		_BuoyancyForceTerm_Build,
		_BuoyancyForceTerm_Initialise,
		_BuoyancyForceTerm_Execute,
		_BuoyancyForceTerm_Destroy,
		_BuoyancyForceTerm_AssembleElement,
		_BuoyancyForceTerm_CalcGravity,
		name );
}

void _BuoyancyForceTerm_Construct( void* forceTerm, Stg_ComponentFactory* cf, void* data ) {
	BuoyancyForceTerm*          self             = (BuoyancyForceTerm*)forceTerm;
	Dictionary*		dict;
	Dictionary_Entry_Value*	tmp;
	char*			rootKey;
	FeVariable*                 temperatureField;
	double                      gravity;
	Bool                        adjust;
	Materials_Register*         materials_Register;
	unsigned		    nDims;
	Dictionary_Entry_Value*	    direcList;
	double*			    direc;
	unsigned		    d_i;
	PICelleratorContext*	    context;

	/* Construct Parent */
	_ForceTerm_Construct( self, cf, data );

	dict = Dictionary_Entry_Value_AsDictionary( Dictionary_Get( cf->componentDict, self->name ) );
	temperatureField = Stg_ComponentFactory_ConstructByKey( cf, self->name, "TemperatureField", FeVariable, False, data ) ;
	gravity          = Stg_ComponentFactory_GetDouble( cf, self->name, "gravity", 0.0 );
	adjust           = Stg_ComponentFactory_GetBool( cf, self->name, "adjust", False );

	direcList = Dictionary_Get( dict, "gravityDirection" );
	if( direcList ) {
		nDims = Dictionary_Entry_Value_GetCount( direcList );
		direc = AllocArray( double, nDims );
		for( d_i = 0; d_i < nDims; d_i++ ) {
			tmp = Dictionary_Entry_Value_GetElement( direcList, d_i );
			rootKey = Dictionary_Entry_Value_AsString( tmp );
			if( !Stg_StringIsNumeric( rootKey ) )
				tmp = Dictionary_Get( cf->rootDict, rootKey );
			direc[d_i] = Dictionary_Entry_Value_AsDouble( tmp );
		}
		if( nDims == 2 )
			Vec_Norm2D( direc, direc );
		else
			Vec_Norm3D( direc, direc );
	}
	else
		direc = NULL;

	context = (PICelleratorContext*)self->context;
	assert( Stg_CheckType( context, PICelleratorContext ) );
	materials_Register = context->materials_Register;
	assert( materials_Register );

	_BuoyancyForceTerm_Init( self, temperatureField, gravity, adjust, materials_Register );
	self->gHat = direc;
}

void _BuoyancyForceTerm_Build( void* forceTerm, void* data ) {
	BuoyancyForceTerm*               self               = (BuoyancyForceTerm*)forceTerm;
	BuoyancyForceTerm_MaterialExt*   materialExt;
	Material_Index                   material_I;
	Material*                        material;
	Materials_Register*              materials_Register = self->materials_Register;
	IntegrationPointsSwarm*          swarm              = (IntegrationPointsSwarm*)self->integrationSwarm;
	MaterialPointsSwarm**            materialSwarms;
	Index                            materialSwarm_I;
	Name                             name;
	AbstractContext*                 context;
	Stg_ComponentFactory*            cf;

	cf = self->context->CF;

	_ForceTerm_Build( self, data );

	if ( self->temperatureField )
		Stg_Component_Build( self->temperatureField, data, False );

	/* Sort out material extension stuff */
	self->materialExtHandle = Materials_Register_AddMaterialExtension( 
			self->materials_Register, 
			self->type, 
			sizeof(BuoyancyForceTerm_MaterialExt) );
	for ( material_I = 0 ; material_I < Materials_Register_GetCount( materials_Register ) ; material_I++) {
		material = Materials_Register_GetByIndex( materials_Register, material_I );
		materialExt = ExtensionManager_GetFunc( material->extensionMgr, material, self->materialExtHandle );

		materialExt->density = Stg_ComponentFactory_GetDouble( cf, material->name, "density", 0.0 );
		materialExt->alpha   = Stg_ComponentFactory_GetDouble( cf, material->name, "alpha",   0.0 );
	}
	
	/* Create Swarm Variables of each material swarm this ip swarm is mapped against */
	materialSwarms = IntegrationPointMapper_GetMaterialPointsSwarms( swarm->mapper, &(self->materialSwarmCount) );
	self->densitySwarmVariables = Memory_Alloc_Array( MaterialSwarmVariable*, self->materialSwarmCount, "DensityVariables" );
	self->alphaSwarmVariables   = Memory_Alloc_Array( MaterialSwarmVariable*, self->materialSwarmCount, "AlphaVariables" );
	
	for ( materialSwarm_I = 0; materialSwarm_I < self->materialSwarmCount; ++materialSwarm_I ) {
		name = Stg_Object_AppendSuffix( materialSwarms[materialSwarm_I], "Density" );
		self->densitySwarmVariables[materialSwarm_I] = MaterialSwarmVariable_New( 
				name, 
				materialSwarms[materialSwarm_I], 
				1, 
				self->materials_Register, 
				self->materialExtHandle, 
				GetOffsetOfMember( *materialExt, density ) );
		Memory_Free( name );

		name = Stg_Object_AppendSuffix( materialSwarms[materialSwarm_I], "Alpha" );
		self->alphaSwarmVariables[materialSwarm_I] = MaterialSwarmVariable_New( 
				name, 
				materialSwarms[materialSwarm_I], 
				1, 
				self->materials_Register, 
				self->materialExtHandle, 
				GetOffsetOfMember( *materialExt, alpha ) );
		Memory_Free( name );
	
		/* Build new Swarm Variables */
		Stg_Component_Build( self->densitySwarmVariables[materialSwarm_I], data, False );
		Stg_Component_Build( self->alphaSwarmVariables[materialSwarm_I],   data, False );
	}
}

void _BuoyancyForceTerm_Initialise( void* forceTerm, void* data ) {
	BuoyancyForceTerm*             self             = (BuoyancyForceTerm*)forceTerm;
	Index                          i;

	_ForceTerm_Initialise( self, data );

	if ( self->temperatureField )
		Stg_Component_Initialise( self->temperatureField, data, False );
	
	for ( i = 0; i < self->materialSwarmCount; ++i ) {
		Stg_Component_Initialise( self->densitySwarmVariables[i], data, False );
		Stg_Component_Initialise( self->alphaSwarmVariables[i],   data, False );
	}
}

void _BuoyancyForceTerm_Execute( void* forceTerm, void* data ) {
	_ForceTerm_Execute( forceTerm, data );
}

void _BuoyancyForceTerm_Destroy( void* forceTerm, void* data ) {
	_ForceTerm_Destroy( forceTerm, data );
}


void _BuoyancyForceTerm_AssembleElement( void* forceTerm, ForceVector* forceVector, Element_LocalIndex lElement_I, double* elForceVec ) {
	BuoyancyForceTerm*               self               = (BuoyancyForceTerm*) forceTerm;
	IntegrationPoint*                particle;
	BuoyancyForceTerm_MaterialExt*   materialExt;
	Particle_InCellIndex             cParticle_I;
	Particle_InCellIndex             cellParticleCount;
	Element_NodeIndex                elementNodeCount;
	Dimension_Index                  dim                = forceVector->dim;
	IntegrationPointsSwarm*          swarm              = (IntegrationPointsSwarm*)self->integrationSwarm;
	FeMesh*              mesh               = forceVector->feVariable->feMesh;
	Node_ElementLocalIndex           eNode_I;
	Cell_Index                       cell_I;
	ElementType*                     elementType;
	Dof_Index                        nodeDofCount;
	double                           gravity;
	double                           detJac             = 0.0;
	double                           factor;
	double                           Ni[27];
	double                           force;
	double*                          xi;
	Material*                        material;
	FeVariable*                      temperatureField   = self->temperatureField;
	double                           temperature        = 0.0;
	double*				 gHat;
	unsigned			d_i;

	double totalWeight = 0.0;
	double adjustFactor = 0.0;

	elementType       = FeMesh_GetElementType( mesh, lElement_I );
	elementNodeCount  = elementType->nodeCount;
	nodeDofCount      = dim;
	cell_I            = CellLayout_MapElementIdToCellId( swarm->cellLayout, lElement_I );
	cellParticleCount = swarm->cellParticleCountTbl[cell_I];
	gHat		  = self->gHat;

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
		gravity = BuoyancyForceTerm_CalcGravity( self, (Swarm*)swarm, lElement_I, particle );
		force = materialExt->density * gravity * (1.0 - materialExt->alpha * temperature);
		factor = detJac * particle->weight * adjustFactor * force;

		/* Apply force in the correct direction */
		for( eNode_I = 0 ; eNode_I < elementNodeCount; eNode_I++ ) {
			if( gHat ) {
				for( d_i = 0; d_i < dim; d_i++ )
					elForceVec[ eNode_I * nodeDofCount + d_i ] += gHat[d_i] * factor * Ni[ eNode_I ] ;
			}
			else {
				elForceVec[ eNode_I * nodeDofCount + J_AXIS ] += -1.0 * factor * Ni[ eNode_I ] ;
			}
		}
	}
	
}

/* The default implementation is for the gravity to be constant. */
double _BuoyancyForceTerm_CalcGravity( void* forceTerm, Swarm* swarm, Element_DomainIndex dElement_I, void* particle ) {
	BuoyancyForceTerm*               self               = (BuoyancyForceTerm*) forceTerm;

	return self->gravity;
}
